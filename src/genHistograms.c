#include <assert.h>
#include <string.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <sys/time.h>
#include <mpi.h>
#include <omp.h>
#include <unistd.h>
//#include "metag.h"

typedef struct {
        size_t offset;
        size_t offset2;
	uint8_t readIdExtra;
        uint32_t readId;
        int fileNo;
        uint32_t kmerFreqCount[1048576+2];
}fileIndex;

char **fileNames=NULL;
size_t *fileSize=NULL, totalFastqSize=0;
uint32_t curReadId=0;
size_t fastqBufferSize=500000000;
int numPEFiles, numSEFiles, numInterleavedFiles;
uint64_t *kmerList=NULL;

struct stat buf;
uint32_t K=0;

uint64_t kmerMask=0;
int rotater=0;
uint64_t totalKmers=0;
void getKmerMask(uint32_t k)
{
//	kmerMask = (((1ul<<(k<<1))-1) >> ((k<<1)-20)) << ((k<<1)-20);
	kmerMask = (((1ul<<k)-1) >> (k-20)) << (32+k-20);
	rotater = 32+k-20;
}

extern int getCanonKmers (char *readStr, int readLen, uint64_t *rangeHist, uint32_t *kmerFreqCount);
extern void getFullStr (uint64_t kmer, FILE *fp);
extern void init(int kLen, int mLen);
FILE *fp1_copy=NULL, *fp2_copy;
size_t ofs_copy=0;
int eof=0;
void bufSeek (int bufNo, FILE *fp, long ofs, char **fqBuffer, size_t *fqInd, size_t *fqSize)
{
	// printf ("bufSeek bufNo %d, ofs %lu fp_copy %p\n", bufNo, ofs, fp);
	fqSize[bufNo]=fqInd[bufNo]=0;
	assert (fseek(fp, ofs, SEEK_SET) == 0);
	
	fqSize[bufNo]=fread(fqBuffer[bufNo], 1, fastqBufferSize, fp);
	if (bufNo == 0)
		ofs_copy=ofs;
	printf ("After Read bufNo %d size %lu\n", bufNo, fqSize[bufNo]);
}
char * getNextRead(int bufNo, int *len, int print, char **fqBuffer, size_t *fqInd, size_t *fqSize)
{
	*len=0;
	if (fqInd[bufNo] >= fqSize[bufNo])
		return NULL;
	//int ret=0;
	size_t curInd = fqInd[bufNo];
	if (print)
		printf ("getNextLine bufSize %lu, curInd %lu\n", fqSize[bufNo], curInd);
	while (fqBuffer[bufNo][curInd] != '\n') curInd++;
	curInd++;
	char *readPtr = &fqBuffer[bufNo][curInd];
	while (fqBuffer[bufNo][curInd] != '\n') {curInd++; (*len)++;}
	fqBuffer[bufNo][curInd]='\0';
	while (fqBuffer[bufNo][curInd] != '\n') curInd++; curInd++;
	curInd+=*len+1;
	fqInd[bufNo]=curInd;
//	sleep(20);
//	printf ("%d Return NULL size %lu, curIndd %lu\n", bufNo, fqSize[bufNo], curInd);
	return readPtr;
}

int main(int argc, char **argv)
{
	char *line1_copy=NULL, *line2_copy=NULL;
        int i, j, k, len;
	if  (argc < 5)
        {
                printf ("Usage: %s kmerSize numPEFiles numSEFiles "
						"numInterleavedFiles a.1.fastq a.2.fastq "
						"b.1.fastq b.2.fastq ...\n", argv[0]);
                exit(1);
        }
	int numTasks = 1, rank = 0, provided=0;
	MPI_Init_thread(&argc,&argv, MPI_THREAD_MULTIPLE, &provided);
        MPI_Comm_size(MPI_COMM_WORLD,&numTasks);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);

//	numTasks=64;
//	rank=14;
        char hostname[100];
        gethostname (hostname, 100);
	fprintf (stderr, "numTasks %d rank %d hostname %s, Thread support %d, requested %d\n", numTasks, rank, hostname, provided, MPI_THREAD_MULTIPLE);
	
	K=(uint32_t)atoi(argv[1]);
	init(K, 10);
	getKmerMask(K);

	printf ("kmer Mask %" PRIu64 ", Rotater %d\n", kmerMask, rotater);
        numPEFiles = atoi(argv[2])*2;
        numSEFiles = atoi(argv[3]);
        numInterleavedFiles = atoi(argv[4]);
        argv = &argv[4];
        argc-=4;
        printf ("Argc %d\n", argc);
        fileNames = (char **)malloc((argc-1) * sizeof(char *));
        fileSize = (size_t *)malloc((argc-1) * sizeof(size_t));
        assert (fileNames != NULL);
        assert (fileSize != NULL);
	for (i=1; i<argc; i++)
	{
		fileNames[i-1]=(char *)malloc(500);
                assert (fileNames[i-1] != NULL);
                sprintf (fileNames[i-1], "%s", argv[i]);
                int status = stat (fileNames[i-1], &buf);
		if (status != 0) printf ("NOT FOUND: %s\n", fileNames[i-1]);
                assert (status == 0);
                fileSize[i-1]=buf.st_size;
                printf ("file %s Size %lu\n", fileNames[i-1], buf.st_size);
		totalFastqSize += buf.st_size;
        }
	char idxFile[200];
	char outFile[200];
#ifdef MIN3
	sprintf (outFile, "%s.KmerIdx_3min%d", argv[1],K);
#else
	#ifdef ALLMINS
		sprintf (outFile, "%s.KmerIdx_allmins%d", argv[1],K);
	#else
		sprintf (outFile, "%s.KmerIdx%d", argv[1],K);
	#endif
#endif
//	FILE *fp_out = fopen(outFile, "w");
//	assert (fp_out != NULL);
	sprintf (idxFile, "%s.idx", argv[1]);
        int status = stat (idxFile, &buf);
        assert (status == 0);
        int numRecs = buf.st_size/sizeof(fileIndex)-1;
	printf ("NumPartitions %d\n", numRecs);
	uint32_t *readCounts = (uint32_t *)calloc ((numRecs+2), sizeof(uint32_t));
	assert (readCounts != NULL);

	MPI_Status mpistatus;
	MPI_File fh, fhout;
	MPI_Offset offset;
	MPI_File_open  (MPI_COMM_WORLD, idxFile, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_File_open  (MPI_COMM_WORLD, outFile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fhout);
	uint64_t totalKmers = 0, totalReads = 0;

	uint64_t **rangeHistCombined = NULL;
	int nThreads=0;
	#pragma omp parallel private (i,j,k) reduction(+:totalKmers, totalReads)
	{
		int numthreads = omp_get_num_threads();
		int tid = omp_get_thread_num();
		int chunksPerThread = numRecs/numthreads/numTasks;
		if (numRecs % numthreads != 0) chunksPerThread++;
		int readOffset = (rank * numthreads + tid) * chunksPerThread;
		if ((readOffset + chunksPerThread) > numRecs) 
			chunksPerThread = numRecs-readOffset;
		printf ("tid %d startOfs %d numChunks %d\n", rank*numthreads + tid, readOffset, chunksPerThread);
		int count=0;
//		printf ("Per Thread Memory %lu\n", sizeof(fileIndex)*(chunksPerThread+1));
		fileIndex *fI = (fileIndex *)malloc((chunksPerThread+1) * sizeof(fileIndex));
		assert (fI != NULL);
		MPI_File_read_at (fh, readOffset*sizeof(fileIndex), fI, (chunksPerThread+1) * sizeof(fileIndex), MPI_BYTE, &mpistatus);
		MPI_Get_count (&mpistatus, MPI_BYTE, &count);
/*		FILE  *fp = fopen (idxFile, "r");
		assert (fp != NULL);
		fseek (fp, readOffset*sizeof(fileIndex), SEEK_SET);
		count = fread (fI, sizeof(fileIndex), chunksPerThread+1, fp);
		fclose (fp);*/
		printf ("count %d chunksPerThread+1 %d  sizeof(fileIndex) %lu\n", count, chunksPerThread+1, sizeof(fileIndex));
		assert (count == ((chunksPerThread+1)  * sizeof(fileIndex) ));

		#pragma omp master 
		{
			nThreads = omp_get_num_threads();
			rangeHistCombined = (uint64_t **) malloc (numthreads *sizeof (uint64_t *));
			assert (rangeHistCombined != NULL);
			for (i=0; i<numthreads; i++)
			{
				rangeHistCombined[i] = (uint64_t *) calloc (1048576+2, sizeof(uint64_t));
				assert (rangeHistCombined[i] != NULL);
			}
		}
		#pragma omp barrier
		uint64_t *rangeHist=rangeHistCombined[tid];
		assert (rangeHist != NULL);
		
		for (i=0; i<chunksPerThread; i++)
		{
			bzero (fI[i].kmerFreqCount, (1048576+2)*sizeof(uint32_t));
			char *fqBuffer[2]={NULL,NULL};
			size_t fqInd[2]={0,0};
			size_t fqSize[2]={0,0};

			int startfile = fI[i].fileNo;
			int endfile = fI[i+1].fileNo;
			size_t bytesToRead = 0;
			for (j=startfile; j<=endfile; j++)
			{
				size_t startOfs = fI[i].offset;
				size_t endOfs = fI[i+1].offset;
				if (j != startfile) startOfs=0;
				if (j != endfile) endOfs = fileSize[j];
				bytesToRead += endOfs-startOfs;
			}
			printf ("tid %d (%d,%lu)-(%d,%lu) bytestoread %lu\n", /*rank*numthreads +*/ tid, fI[i].fileNo, fI[i].offset, fI[i+1].fileNo, fI[i+1].offset, bytesToRead);
			fqBuffer[0] = (char *)malloc(bytesToRead+10);
			assert (fqBuffer[0] != NULL);
			int fNo = fI[i].fileNo;
			FILE *fpFq = fopen (fileNames[fNo], "r");
			assert (fpFq != NULL);
			fseek (fpFq, fI[i].offset, SEEK_SET);
			size_t nread = 0;
			while (nread < bytesToRead)
			{
				nread += fread (&fqBuffer[0][nread], 1, bytesToRead-nread, fpFq);
				if (nread < bytesToRead)
				{
					fclose (fpFq);
					fNo++;
					fpFq = fopen (fileNames[fNo], "r");
					assert (fpFq != NULL);
				}
			}

			fclose (fpFq);
			if (fNo != fI[i+1].fileNo)
			{
//				printf ("\tTid %d (%d,%lu)-(%d,%lu) bytestoread %lu fNo %d\n", tid, fI[i].fileNo, fI[i].offset, fI[i+1].fileNo, fI[i+1].offset, bytesToRead, fNo);
			}
/*			if (tid == 37)
			{
				FILE *fout = fopen ("tid37", "w");
				assert (fout != NULL);
				fwrite (fqBuffer[0], 1, bytesToRead, fout);
				fclose (fout);
			}*/
			assert (fNo == fI[i+1].fileNo);
			fqSize[0] = bytesToRead;
			int readLen=0;
			char * readPtr = getNextRead(0, &readLen, 0, fqBuffer, fqInd, fqSize);
			while (readPtr != NULL)
			{
				totalKmers += (uint64_t)getCanonKmers(readPtr, readLen, rangeHist, fI[i].kmerFreqCount);
				totalReads++;
				readCounts [readOffset+ i+2]++;
//				if (tid == 37)
//					printf ("-%d,%lu:%s-\n", readLen, totalKmers, readPtr);
//				if (totalKmers > 500) break;
				readPtr = getNextRead(0, &readLen, 0, fqBuffer, fqInd, fqSize);
			}
			free (fqBuffer[0]);
		}
		printf ("THREAD %d TotalKmers %lu, TotalReads %lu\n", rank*numthreads + tid, totalKmers, totalReads);
		#pragma omp barrier
		
		#pragma omp master
		{
			MPI_Allreduce (MPI_IN_PLACE, readCounts, numRecs+2, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
			if (rank==0)
			{
				for (i=0; i<=numRecs+1; i++)
					printf ("%lu, ", readCounts[i]);
				printf ("\n");
			}
		}
		#pragma omp barrier
		uint64_t readcount=0;
		for (i=1; i<=numRecs; i++)
		{
			readcount += (uint64_t)readCounts[i+1];
			if ((i >= readOffset) && (i<= readOffset+chunksPerThread))
			{
				fI[i-readOffset].readId = (uint32_t)(readcount % UINT32_MAX);
				uint32_t remain = (uint32_t)(readcount / UINT32_MAX);
				assert (remain < 256);
				fI[i-readOffset].readIdExtra = (uint8_t)remain;
			}
		}
		#pragma omp barrier
		if ((rank < (numTasks-1)) || (tid != numthreads-1))
			MPI_File_write_at (fhout, readOffset*sizeof(fileIndex), fI, chunksPerThread * sizeof(fileIndex), MPI_BYTE, &mpistatus);
		else
		{
			printf ("write %lu bytes at %lu offset\n", chunksPerThread*sizeof(fileIndex), readOffset*sizeof(fileIndex));
			MPI_File_write_at (fhout, readOffset*sizeof(fileIndex), fI, (chunksPerThread+1) * sizeof(fileIndex), MPI_BYTE, &mpistatus);
		}
//		MPI_Get_count (&mpistatus, MPI_BYTE, &count);
//		assert (count == (int)((chunksPerThread+1) * sizeof(fileIndex)));
		#pragma omp master
		{
			for (i=1; i<numthreads; i++)
			{
				for (j=0; j<1048576+2; j++)
				{
					rangeHistCombined [0][j] += rangeHistCombined[i][j];
				}
			}
		}
		
//		fqBuffer[1] = (char *)malloc(fastqBufferSize+10);
//		assert ((fqBuffer[0] != NULL) && (fqBuffer[1] != NULL));
//		free (fqBuffer[0]);
//		free (fqBuffer[1]);
	}

	printf ("[%d]TOTALKMERS %lu, TOTALREADS %lu\n", rank, totalKmers, totalReads);

	if (rank == 0)
	{
		MPI_Reduce (MPI_IN_PLACE, &totalKmers, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce (MPI_IN_PLACE, rangeHistCombined[0], 1048576+2, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		for (i=1; i<=1048576; i++)
			rangeHistCombined[0][i] = rangeHistCombined[0][i-1] + rangeHistCombined[0][i+1];
		printf ("TOTALKMERS2 %lu TOTALKMERS3 %lu\n", rangeHistCombined[0][1048576], totalKmers);
		assert (rangeHistCombined[0][1048576]==totalKmers);

		char histFile[200];
#ifdef MIN3
		sprintf (histFile, "%s.kmerHist_3min%d", argv[1], K);
#else
	#ifdef ALLMINS
		sprintf (histFile, "%s.kmerHist_allmins%d", argv[1], K);
	#else
		sprintf (histFile, "%s.kmerHist%d", argv[1], K);
	#endif
#endif
		FILE *fp_hist = fopen(histFile, "w");
		assert (fp_hist != NULL);
		fwrite (rangeHistCombined[0], sizeof(uint64_t), 1048576+1, fp_hist);
		fclose (fp_hist);
	}
	else
	{
		MPI_Reduce (&totalKmers, &totalKmers, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce (rangeHistCombined[0], rangeHistCombined[0], 1048576+2, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	}
//	MPI_Barrier (MPI_COMM_WORLD);
/*        fileIndex *fI = (fileIndex *)malloc(numRecs * sizeof(fileIndex));
        assert (fI != NULL);
        fread(fI, sizeof(fileIndex), numRecs, fp);
        fclose (fp);

	//char line1[200],line2[200];
	size_t ofs1=0,ofs2=0;
	fp1_copy = fopen (fileNames[file1], "r");
//	printf ("fp1_copy %p\n", fp1_copy);
	assert ((fp != NULL) && (fp1_copy != NULL));
	if (file2 < numPEFiles)
	{
		fp2_copy = fopen (fileNames[file2], "r");
		assert (fp2_copy != NULL);
	}
        for (i=0; i<numRecs-1; i++)
	{
//                printf ("%d(%s)-%d %lu(ofs1 %lu),%lu FILE %p %p\n", i, fileNames[file1], fI[i].fileNo, fI[i].offset, ofs1, fI[i].offset2, fp1_copy, fp2_copy);
//		assert (ofs1 == fI[i].offset);
		ofs1 = fI[i].offset;
		ofs2 = fI[i].offset2;
		file1 = fI[i].fileNo;
		bufSeek (0, fp1_copy, ofs1);
		if (fp2_copy != NULL)
		{
			file2 = fI[i].fileNo+1;
//			assert (fseek (fp2, ofs2, SEEK_SET) == 0);
			bufSeek (1, fp2_copy, ofs2);
		}
		else
			file2 = -1;
		fI[i].readId = curReadId;
//		printf ("chunk[%d] %u\n", i, curReadId);
		line1_copy=getNextLine (0, fp1_copy, &len, 0);
//		line1[len]='\0';
//		assert (strcmp(line1_copy, line1)==0);
//		printf ("Actual *%s*\n", line1);
		if (fp2_copy != NULL)
		{
			line2_copy=getNextLine (1, fp2_copy, &len, 0);
		}
		printf ("[%d]ReadId %u totalKmersSoFar %lu\n", i, curReadId, totalKmers);
		while ((file1 < fI[i+1].fileNo) || ((file1 == fI[i+1].fileNo) && (ofs_copy < fI[i+1].offset)))
		{
			line1_copy=getNextLine (0, fp1_copy, &len, 0);
//			line1_copy[len]='\0';
//			printf ("*%s*\n", line1_copy);
			totalKmers += (uint64_t)getCanonKmers(line1_copy, len, rangeHist, fI[i].kmerFreqCount);
//			line1[len]='\0';
//			printf ("Actual *%s*\n", line1);
//			printf ("Compare %d\n", strcmp(line1_copy, line1));
//			assert (strcmp(line1_copy, line1)==0);
			
//			printf ("%s len %d", line1, strlen(line1));
			//genKmers
			if (fp2_copy != NULL)
			{
				line2_copy=getNextLine (1, fp2_copy, &len, 0);
//				line2_copy[len]='\0';
//				printf ("*%s*\n", line2_copy);
				totalKmers += (uint64_t)getCanonKmers(line2_copy, len, rangeHist, fI[i].kmerFreqCount);
			}
			//genKmers
			line1_copy=getNextLine (0, fp1_copy, &len, 0);
//			line1[len]='\0';
//			printf ("Actual *%s*\n", line1);
//			printf ("Compare %d\n", strcmp(line1_copy, line1));
//			assert (strcmp(line1_copy, line1)==0);
			line1_copy=getNextLine (0, fp1_copy, &len, 0);
//			line1[len]='\0';
//			printf ("Actual *%s*\n", line1);
//			printf ("Compare %d\n", strcmp(line1_copy, line1));
//			assert (strcmp(line1_copy, line1)==0);
			if (eof==0)
			{
				{
					line1_copy=getNextLine (0, fp1_copy, &len, 0);
				}
//				printf ("line1 %s\n", line1_copy);
//				if (s != NULL)	
//				{
//					line1[len]='\0';
//					if (curReadId>=55073606)
//					{
//						printf ("3:curReadId %u Read %s\n", curReadId, line1);
//						printf ("line1_copy %p\n", line1_copy);
//					}
//					printf ("Actual *%s*\n", line1);
//					printf ("Compare %d\n", strcmp(line1_copy, line1));
//					assert (strcmp(line1_copy, line1)==0);
//				}
			}
			if (fp2_copy != NULL)
			{
				line2_copy=getNextLine (1, fp2_copy, &len, 0);
				line2_copy=getNextLine (1, fp2_copy, &len, 0);
				if (eof==0)
				{
					line2_copy=getNextLine (1, fp2_copy, &len, 0);
				}
			}
//			if (curReadId>=55073606)
//				printf ("4:curReadId %u Read %s\n", curReadId, line1);
			if (file1 >= (numPEFiles + numSEFiles))
			{
				line1_copy=getNextLine (0, fp1_copy, &len, 0);
//				line1[len]='\0';
//				printf ("Actual *%s*\n", line1);
//				printf ("Compare %d\n", strcmp(line1_copy, line1));
//				assert (strcmp(line1_copy, line1)==0);
//				line1_copy[len]='\0';
//				printf ("*%s*\n", line1_copy);
				totalKmers += (uint64_t)getCanonKmers(line1_copy, len, rangeHist, fI[i].kmerFreqCount);
				//genKmers
				line1_copy=getNextLine (0, fp1_copy, &len, 0);
//				line1[len]='\0';
//				printf ("Actual *%s*\n", line1);
//				printf ("Compare %d\n", strcmp(line1_copy, line1));
//				assert (strcmp(line1_copy, line1)==0);
				line1_copy=getNextLine (0, fp1_copy, &len, 0);
//				line1[len]='\0';
//				printf ("Actual *%s*\n", line1);
//				printf ("Compare %d\n", strcmp(line1_copy, line1));
//				assert (strcmp(line1_copy, line1)==0);
//				if (!feof(fp))
				if (eof==0)
				{
					line1_copy=getNextLine (0, fp1_copy, &len,0);
//					line1[len]='\0';
//					printf ("Actual *%s*\n", line1);
//					printf ("Compare %d\n", strcmp(line1_copy, line1));
//					assert (strcmp(line1_copy, line1)==0);
				}
			}
			curReadId++;
			if (eof)
			{
				printf ("eof\n");
				printf ("fp1_copy %p file1 %d argc %d kmersSoFar %lu feof (fp1)?%d ind %lu size %lu\n", 
						((void *) fp1_copy), file1, argc, totalKmers, feof(fp1_copy), fqInd[0], fqSize[0]);
				fclose (fp1_copy);
				if (fp2_copy != NULL)
				{
					fclose(fp2_copy);
					file1+=2;
					if (file1 < (argc-1))
						fp1_copy = fopen (fileNames[file1], "r");
					else
						fp1_copy=NULL;
					file2+=2;
					if (file2 < numPEFiles)
					{
						fp2_copy = fopen (fileNames[file2], "r");
						if (fp2_copy != NULL)
						{
//							printf ("file2 %d numPEFiles %d CallBuffSeek 1\n", file2, numPEFiles);
							bufSeek (1, fp2_copy, 0);
							line2_copy=getNextLine (1, fp2_copy, &len, 0);
						}
					}
					else
						fp2_copy = NULL;
				}
				else
				{
					file1++;
					fp1_copy = fopen (fileNames[file1], "r");
				}
//				ofs1=0;
				if (fp1_copy != NULL)
				{
					bufSeek (0, fp1_copy, 0);
					line1_copy=getNextLine (0, fp1_copy, &len, 0);
				}
//				sleep(10);
			}
		}
	}
	printf ("[%d]ReadId %u totalKmersSoFar %lu\n", i, curReadId, totalKmers);
	fI[i].readId = curReadId;
	fwrite (fI, sizeof (fileIndex), numRecs, fp_out);
	for (i=1; i<=1048576; i++)
	{
		rangeHist[i]=rangeHist[i-1]+rangeHist[i+1];
//		if (i <= 0x1717c)
//			printf ("%d - %lu\n", i, rangeHist[i]);
	}
	assert (rangeHist[1048576]==totalKmers);
	FILE *fp_hist = fopen(histFile, "w");
	assert (fp_hist != NULL);
	fwrite (rangeHist, sizeof(uint64_t), 1048576+1, fp_hist);
	fclose (fp_hist);

	uint64_t newTotalKmers=0;
	for (i=0; i<numRecs; i++)
	{
		for (k=0; k<1048576; k++)
		{
			newTotalKmers += fI[i].kmerFreqCount[k];
		}
	}
	printf ("totalKmers %lu\n", totalKmers);
	assert (newTotalKmers == totalKmers);
	free (fI);*/
//	fclose (fp_out);
	for (i=0; i<(argc-1); i++)
	{
		free (fileNames[i]);
	}
	free(fileNames);
#if 0
	fclose(fp1_copy);
	fclose(fp2_copy);
#endif
	free(fileSize);
	for (i=0; i<nThreads; i++)
		free (rangeHistCombined[i]);
	free (rangeHistCombined);
	free (readCounts);
	MPI_File_close (&fh);
	MPI_File_close (&fhout);
	MPI_Finalize();

	return 0;
}

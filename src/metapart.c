#define _GNU_SOURCE
#include <sched.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/stat.h>
#include <stdint.h>
#include <inttypes.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <getopt.h>
#include <pthread.h>
#define RADIX_OPT 0
#define USE_MPI 1
#if USE_MPI
#include <mpi.h>
#endif
#define USE_OMP 1
//Comp-Comm overlap
#define NONBLOCKING_KE 0
//Kmers or Reads to be sent in 1 call
#define MERGE_COMM_OPT 0

#define ASYNC_SIZE 400000000
#if MERGE_COMM_OPT
	#define MERGE_COMM_SIZE 1000000000
#else
#endif
//Multi Scan optimization
#define SCAN_OPT  1
#define MERGE_LOG 1
#if USE_OMP
#include <omp.h>
#endif
// #include "metag.h"
#include "unionfind.h"

typedef struct {
        size_t offset;
        size_t offset2;
	uint8_t readIdExtra;
        uint32_t readId;
        int fileNo;
        uint32_t kmerFreqCount[1048576+2];
}fileIndex;

typedef struct {
	uint32_t **sendKmers;
	size_t *sendThreadOfs;
	size_t *sendThreadCount;
	size_t *threadGenerateCount;
	size_t *rcvThreadOfs;
	size_t *rcvThreadCount;
}asyncData;

char **fileNames=NULL;
size_t *fileSize=NULL, totalFastqSize=0;
size_t fastqBufferSize=500000000;
char *fastqBuffer1=NULL, *fastqBuffer2=NULL;

int numPEFiles, numSEFiles, numInterleavedFiles, K=0;
int numTasks=1, numThreads, rank, scans;
int totalFastqPartitions=0;
fileIndex *fI = NULL;
uint64_t *rangeHist;
uint64_t *scanKmerRanges=NULL, *processKmerRanges=NULL, *threadKmerRanges=NULL;
int *fastqRanges=NULL;
size_t *recvThreadCount=NULL;
size_t *recvThreadOfs=NULL;
uint64_t *threadKmerCounts=NULL;
uint64_t numKmers=0;
uint64_t numReads=0, numReads8Mul=0;
int hwThreads=0;
int numCommIters=1;
long *threadKmerGenTime=NULL;
FILE **openFp=NULL;//2*numThreads
int *openFileNo=NULL;
uint8_t *filteredReads = NULL;

uint32_t *kmerDegree=NULL;
uint8_t *parent=NULL, *otherRanksParent=NULL;
uint32_t *compSizes=NULL;
long kmerGenWorkTime=0, kmerGenTime=0, commTime=0, memcpyTime=0, sortTime1=0, sortTime2=0, filterTime=0, unionFindTime1=0, unionFindTime2=0, mergeTime=0, peUnionTime=0, outputTime=0;

//extern void getCanonKmers (char *readStr, uint32_t **kR, uint64_t *curSendOfs, int readLen, uint64_t startKmer, uint64_t endKmer, uint64_t kmerMask, int rotater, uint32_t curReadId, uint64_t *scanProcessKmerRange, int tid, int rank);
extern void getCanonKmers (char *readStr, uint16_t **kR, uint64_t *curSendOfs, int readLen, uint64_t startKmer, uint64_t endKmer, uint64_t kmerMask, int rotater, uint64_t curReadId, uint64_t *scanProcessKmerRange, int tid, int rank/*, uint64_t MAX_TUPLES*/);
extern uint64_t getMmerPrefix (uint64_t kH, uint64_t kL, int m);
extern void init(int kLen, int numP, int numThreads, int mmerSize);
extern void getKmerStr64 (uint64_t kmerl, uint64_t kmerh, FILE *fp, int ind);
size_t *sendSizes=NULL;
size_t *sendThreadSizes=NULL;
uint64_t kmerMask=0;
int rotater=0;

char oPrefix[100];

uint16_t *commonSendBuffer=NULL, *rcvKmers=NULL, *sortedKmers=NULL;
void getKmerMask(uint32_t k)
{
//      kmerMask = (((1ul<<(k<<1))-1) >> ((k<<1)-20)) << ((k<<1)-20);
      kmerMask = (((1ul<<k)-1) >> (k-20)) << (32+k-20);
      rotater = 32+k-20;
}

static int hardware_threads(void)
{
        char name[40];
        struct stat st;
        int cpus = -1;
        do {
                sprintf(name, "/sys/devices/system/cpu/cpu%d", ++cpus);
        } while (stat(name, &st) == 0);
        return cpus;
}

#if 0
static void cpu_bind(int cpu_id)
{
        int cpus = hardware_threads();
        size_t size = CPU_ALLOC_SIZE(cpus);

        cpu_set_t *cpu_set = CPU_ALLOC(cpus);
        assert(cpu_set != NULL);
        CPU_ZERO_S(size, cpu_set);
        CPU_SET_S(cpu_id, size, cpu_set);
        assert(pthread_setaffinity_np(pthread_self(),
               size, cpu_set) == 0);
        CPU_FREE(cpu_set);
}
#endif

static int parentCompare(const void *p1, const void *p2)//5 bytes per readid
{
/*        uint8_t *f1=(uint8_t *)p1;
        uint8_t *f2=(uint8_t *)p2;

	uint64_t r1 = (uint64_t)(*((uint32_t *)&f1[1]));
	uint64_t r1extra = (uint64_t)f1[0];
	r1 = r1extra*UINT32_MAX + r1;

	uint64_t r2 = (uint64_t)(*((uint32_t *)&f2[1]));
	uint64_t r2extra = (uint64_t)f2[0];
	r2 = r2extra * UINT32_MAX + r2;*/
	uint32_t r1 = *((uint32_t *)p1);
	uint32_t r2 = *((uint32_t *)p2);
        if (r1 < r2)
                return -1;
        else if (r1 > r2)
                return 1;
        else
                return 0;
}

void bufSeek (int bufNo, FILE *fp, long ofs, size_t *fqSize, size_t *fqInd, char **fqBuffer)
{
        fqSize[bufNo]=fqInd[bufNo]=0;
        assert (fseek(fp, ofs, SEEK_SET) == 0);

        fqSize[bufNo]=fread(fqBuffer[bufNo], 1, fastqBufferSize, fp);
//	printf ("[%d] fqSize[%d]=%lu\n", rank, bufNo, fqSize[bufNo]);
}

char * getNextLine(int bufNo, FILE *fp, int *len, int print, size_t *fqSize, size_t *fqInd, char **fqBuffer, int *eof, size_t *ofs, int f1)
{
        size_t curInd = fqInd[bufNo];
        if (bufNo == 0)
                *eof=0;
        while (fqBuffer[bufNo][curInd] != '\n')
        {
                if (curInd == fqSize[bufNo])
                {
                        memcpy (fqBuffer[bufNo], &fqBuffer[bufNo][fqInd[bufNo]], curInd-fqInd[bufNo]);
                        fqSize[bufNo] = curInd-fqInd[bufNo];
                        fqInd[bufNo] = 0;
                        curInd = fqSize[bufNo];
                        if (!feof(fp))
                        {
                                size_t nRead = fread(&fqBuffer[bufNo][fqSize[bufNo]], 1, fastqBufferSize-fqSize[bufNo], fp);
                                fqSize[bufNo] += nRead;
                        }
                        else
                        {
                                *eof=1;
                        }
                }
                if (*eof) break;
                curInd++;
        }
        if ((*eof==0) && (fqInd[bufNo] < fqSize[bufNo]))
        {
                *(fqBuffer[bufNo]+curInd) = '\0';
                *len=curInd-fqInd[bufNo];
                char *retPtr = fqBuffer[bufNo]+fqInd[bufNo];
                if (bufNo == 0)
                {
                        *ofs += (curInd+1-fqInd[bufNo]);
                        if (*ofs == fileSize[f1])
                        {
                                *eof=1;
                        }
                }
                fqInd[bufNo]=curInd+1;
                return retPtr;
        }
        return NULL;
}
//return type - 0 (paired end), 1 - single end , 2 - interleaved, -1 - EO chunk
int getNextFullRead(size_t *fqInd, size_t *fqSize, char **fqBuffer, char **line1, char **line2, int *len1, int *len2, int interleave, size_t interleaveStart)
{
        int ret=1;
        *line1 = NULL; *line2=NULL;
        *len1 = 0; *len2 = 0;
        if (fqInd[0] >= fqSize[0])
                return -1;
        *line1 = &(fqBuffer[0][fqInd[0]]);
	int oldFqInd=fqInd[0];

	if (fqInd[0] < fqSize[0])
        {
		while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
                fqInd[0]++;
//		fqBuffer[0][fqInd[0]-3]='/';
//		fqBuffer[0][fqInd[0]-2]='1';
	        while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n'))
	        {
        	        fqInd[0]++;
                	*len1 = *len1 +1;
	        }
		fqInd[0]++;
		while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
        	fqInd[0]++;
		fqInd[0]+=*len1;
        	fqInd[0]++;	
		*len1 = fqInd[0]-oldFqInd;
		if (interleave && (fqInd[0] >= interleaveStart))
                {
			oldFqInd=fqInd[0];
			*line2 = &(fqBuffer[0][fqInd[0]]);
                        while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
                        fqInd[0]++;
			while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n'))
                        {
                                fqInd[0]++;
                                *len2 = *len2 +1;
                        }
			fqInd[0]++;
			while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
                        fqInd[0]++;
			fqInd[0]+=*len2;
			fqInd[0]++;
			*len2 = fqInd[0]-oldFqInd;
			ret=2;
			if (fqInd[0] < fqSize[0])
                        {
			}
                }
        }
	if (ret==1)
        {
                if (fqInd[1] < fqSize[1])
                {
			oldFqInd=fqInd[1];
                        *line2 = &(fqBuffer[1][fqInd[1]]);
                        while ((fqInd[1] < fqSize[1]) && (fqBuffer[1][fqInd[1]] != '\n')) fqInd[1]++;
                        fqInd[1]++;
			fqBuffer[1][fqInd[1]-3]='/';
			fqBuffer[1][fqInd[1]-2]='2';
			while ((fqInd[1] < fqSize[1]) && (fqBuffer[1][fqInd[1]] != '\n'))
                        {
                                fqInd[1]++;
                                *len2 = *len2 +1;
                        }
			ret=0;
			fqInd[1]++;
			while ((fqInd[1] < fqSize[1]) && (fqBuffer[1][fqInd[1]] != '\n')) fqInd[1]++;
                        fqInd[1]++;
			fqInd[1]+=*len2;
			fqInd[1]++;
			*len2 = fqInd[1]-oldFqInd;
			if (fqInd[1] < fqSize[1])
                        {
			}
                }
        }
/*    if (ret != 2)
    {
        printf ("[%d] getNextFullRead interleave %d interleaveStart %lu\n", rank, interleave, interleaveStart);
    }
    assert (ret == 2);*/
        return ret;
}
//return type - 0 (paired end), 1 - single end , 2 - interleaved, -1 - EO chunk
int getNextRead(size_t *fqInd, size_t *fqSize, char **fqBuffer, char **line1, char **line2, int *len1, int *len2, int interleave, size_t interleaveStart)
{
        int ret=1;
        *line1 = NULL; *line2=NULL;
        *len1 = 0; *len2 = 0;
        if (fqInd[0] >= fqSize[0])
                return -1;
        *line1 = &(fqBuffer[0][fqInd[0]]);

        while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n'))
        {
                fqInd[0]++;
                *len1 = *len1 +1;
        }
	fqInd[0]++;
	while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
        fqInd[0]++;
	fqInd[0]+=*len1;
        fqInd[0]++;	
	if (fqInd[0] < fqSize[0])
        {
		while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
                fqInd[0]++;
		if (interleave && (fqInd[0] >= interleaveStart))
                {
			*line2 = &(fqBuffer[0][fqInd[0]]);
			while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n'))
                        {
                                fqInd[0]++;
                                *len2 = *len2 +1;
                        }
			fqInd[0]++;
			while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
                        fqInd[0]++;
			fqInd[0]+=*len2;
			fqInd[0]++;
			ret=2;
			if (fqInd[0] < fqSize[0])
                        {
                                while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
                                fqInd[0]++;
			}
                }
        }
	if (ret==1)
        {
                if (fqInd[1] < fqSize[1])
                {
                        *line2 = &(fqBuffer[1][fqInd[1]]);
			while ((fqInd[1] < fqSize[1]) && (fqBuffer[1][fqInd[1]] != '\n'))
                        {
                                fqInd[1]++;
                                *len2 = *len2 +1;
                        }
			ret=0;
			fqInd[1]++;
			while ((fqInd[1] < fqSize[1]) && (fqBuffer[1][fqInd[1]] != '\n')) fqInd[1]++;
                        fqInd[1]++;
			fqInd[1]+=*len2;
			fqInd[1]++;
			if (fqInd[1] < fqSize[1])
                        {
                                while ((fqInd[1] < fqSize[1]) && (fqBuffer[1][fqInd[1]] != '\n')) fqInd[1]++;
                                fqInd[1]++;
			}
                }
        }
/*    if (ret != 2)
    {
        printf ("[%d] getNextRead interleave %d interleaveStart %lu\n", rank, interleave, interleaveStart);
    }
    assert (ret == 2);*/
        return ret;
}
	
void loadChunk (int chunk, size_t *fqInd, size_t *fqSize, char **fqBuffer, int *interleave, size_t *interleaveStart, int tid)
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
	FILE *fp1=NULL, *fp2=NULL;
        fqInd[0] = fqInd[1]=0;
        fqSize[0] = fqSize[1]=0;
        int curFile=fI[chunk].fileNo;
        size_t curOfs1 = fI[chunk].offset;
        size_t curOfs2 = fI[chunk].offset2;
        int endFile=fI[chunk+1].fileNo;
        int startFile = curFile;
        size_t startOfs1 = curOfs1;
        size_t startOfs2 = curOfs2;

        size_t endOfs1 = fI[chunk+1].offset;
        size_t endOfs2 = fI[chunk+1].offset2;

        size_t bytesToRead1=0;
        size_t bufferBytesRemaining1 = fastqBufferSize;
        size_t bufferOfs1=0;

        size_t bytesToRead2=0;
        size_t bufferBytesRemaining2 = fastqBufferSize;
        size_t bufferOfs2=0;
        while ((curFile < endFile) || ((curFile == endFile) && (curOfs1 < endOfs1)))
        {
		if (curFile >= (numPEFiles + numSEFiles))//Interleaved
                {
                    if (*interleave != 1)
                    {
                        *interleave=1;
                        *interleaveStart = bufferOfs1;
//                        printf ("[%d] loadChunk interleaveStart=%lu\n", rank, bufferOfs1);
                    }
                }
/*		if (openFileNo[2*tid] == curFile)
		{
			fp1=openFp[2*tid];
		}
		else
		{*/
                	fp1 = fopen(fileNames[curFile], "r");
/*			openFileNo[2*tid] = curFile;
			if (openFp[2*tid] != NULL)
				fclose (openFp[2*tid]);
			openFp[2*tid] = fp1;
		}*/
		if (fp1==NULL) 
		{
			printf ("Cant open file %s\n", fileNames[curFile]);
			int pid = getpid();
			char line[2048];
			sprintf (line, "/proc/%d/status", pid);
			FILE *statFile = fopen(line, "r");
			assert (statFile != NULL);
			fgets (line, 2048, statFile);
			while (!feof (statFile))
			{
				if (strstr(line,"VmPeak") || strstr(line,"VmHWM"))
				{
					printf ("[%d] %s", rank, line);
				}
				fgets (line, 2048, statFile);
			}
			fclose (statFile);
printf ("[%d] tid %d fp1 NULL, fileNo %d\n", rank, tid, curFile);
			sleep(60);
		}
                assert (fp1 != NULL);
                assert (fseek (fp1, curOfs1, SEEK_SET) == 0);
                fp2 = NULL;
                if (curFile < numPEFiles)
                {
/*			if (openFileNo[2*tid + 1] == (curFile+1))
				fp2 = openFp[2*tid + 1];
			else
			{*/
	                        fp2 = fopen(fileNames[curFile+1], "r");
/*				openFileNo[2*tid + 1] = curFile+1;
				if (openFp[2*tid + 1] != NULL)
					fclose (openFp[2*tid + 1]);
				openFp[2*tid + 1] = fp2;
			}*/
                        assert (fp2 != NULL);
                        assert (fseek (fp2, curOfs2, SEEK_SET) == 0);
                }
                if (curFile < endFile)
                {
                        bytesToRead1 = fileSize[curFile]-curOfs1;
                        if (curFile < numPEFiles)
                                bytesToRead2 = fileSize[curFile+1]-curOfs2;
                }
                else
                {
                        bytesToRead1 = endOfs1-curOfs1;
                        if (curFile < numPEFiles)
                                bytesToRead2 = endOfs2-curOfs2;
                }
                assert (bufferBytesRemaining1 > bytesToRead1);
		size_t curReadOfs = ftell (fp1);
                size_t numBytesRead = fread (&(fqBuffer[0][bufferOfs1]), 1, bytesToRead1, fp1);
		if (numBytesRead != bytesToRead1)
		{
			printf ("ERROR: Chunk %d (File %d,Ofs %lu to File %d, Ofs %lu) curFile %d (%s) Ofs %lu bytesToRead1 %lu read %lu\n", chunk, fI[chunk].fileNo, fI[chunk].offset, fI[chunk+1].fileNo, fI[chunk+1].offset, curFile, fileNames[curFile], curReadOfs, bytesToRead1, numBytesRead);
		}
		assert (numBytesRead == bytesToRead1);
                bufferBytesRemaining1 -= bytesToRead1;
                bufferOfs1 += bytesToRead1;
		if (curFile < numPEFiles)
                {
                        assert (bufferBytesRemaining2 > bytesToRead2);
                        assert (fread (&(fqBuffer[1][bufferOfs2]), 1, bytesToRead2, fp2) == bytesToRead2);
                        bufferBytesRemaining2 -= bytesToRead2;
                        bufferOfs2 += bytesToRead2;
                }

                fclose (fp1);
                if (fp2 != NULL)
                        fclose (fp2);
                if (curFile < numPEFiles) curFile++;
                curFile++;
                curOfs1=curOfs2=0;
        }
        fqSize[0]=bufferOfs1;
        fqSize[1]=bufferOfs2;
	gettimeofday (&eTime, NULL);
      //  fastqTime +=(eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
}
uint64_t getReadAndKmerCount (int chunk, size_t *fqInd, size_t *fqSize, char **fqBuffer, int interleave, size_t interleaveStart, uint16_t **sendKmers, uint64_t *curSendOfs, uint64_t startKmer, uint64_t endKmer, uint64_t kmerMask, int rotater, uint64_t curReadId, size_t *lastKmerAddedByThread, uint64_t *scanProcessKmerRange, int tid)
{
        int len1=0, len2=0;
        char *read1=NULL, *read2=NULL;
        uint64_t numPEReads=0, numSEReads=0, numInterleavedReads=0;
        int ret = getNextRead (fqInd, fqSize, fqBuffer, &read1, &read2, &len1, &len2, interleave, interleaveStart);
#if SCAN_OPT
	if (parent==NULL)
	{
#endif
	        while (ret != -1)
        	{
	               	if (read1 != NULL)
			{
				{
					if (curReadId == 18173149)
						printf ("--read1 len1 %d curSendOfs %lu\n", len1, curSendOfs[0]);
    					getCanonKmers (read1, sendKmers, curSendOfs, len1, startKmer, endKmer, kmerMask, rotater, curReadId, scanProcessKmerRange, tid, rank);
				}
			}
	                if (read2 != NULL)
			{
				if (curReadId == 18173149)
					printf ("--read2 len2 %d curSendOfs %lu\n", len2, curSendOfs[0]);
				if (ret == 2)
		    			getCanonKmers (read2, sendKmers, curSendOfs, len2, startKmer, endKmer, kmerMask, rotater, curReadId+1, scanProcessKmerRange, tid, rank);
				else
		    			getCanonKmers (read2, sendKmers, curSendOfs, len2, startKmer, endKmer, kmerMask, rotater, curReadId, scanProcessKmerRange, tid, rank);
			}
                	if (ret == 0)
	                {
/*      	                read2[len2]='\0';
                	      read1[len1]='\0';
	                      if (numPEReads < 100)
        	                    printf ("*%s*%s*\n", read1, read2);*/
                	      numPEReads++;
	                }
        	        if (ret == 1) numSEReads++;
                	if (ret == 2) {numInterleavedReads++; curReadId++;}
			curReadId++;
        	        ret = getNextRead (fqInd, fqSize, fqBuffer, &read1, &read2, &len1, &len2, interleave, interleaveStart);
	        }
#if SCAN_OPT
	}
	else
	{
            int numRandAccesses=0;
	        while (ret != -1)
        	{
                	if (read1 != NULL)
				getCanonKmers (read1, sendKmers, curSendOfs, len1, startKmer, endKmer, kmerMask, rotater, findRoot(parent, curReadId, &numRandAccesses), scanProcessKmerRange, tid, rank);
	                if (read2 != NULL)
			{
				if (ret == 2)
					getCanonKmers (read2, sendKmers, curSendOfs, len2, startKmer, endKmer, kmerMask, rotater, findRoot(parent, curReadId+1, &numRandAccesses), scanProcessKmerRange, tid, rank);
				else
					getCanonKmers (read2, sendKmers, curSendOfs, len2, startKmer, endKmer, kmerMask, rotater, findRoot(parent, curReadId, &numRandAccesses), scanProcessKmerRange, tid, rank);
			}
                	if (ret == 0)
	                {
                	      numPEReads++;
	                }
        	        if (ret == 1) numSEReads++;
                	if (ret == 2) {numInterleavedReads++; curReadId++;}
			curReadId++;
        	        ret = getNextRead (fqInd, fqSize, fqBuffer, &read1, &read2, &len1, &len2, interleave, interleaveStart);
	        }
	}
#endif
//        printf ("[%d]Num Reads: PE %u SE %u IR %u\n", chunk, numPEReads, numSEReads, numInterleavedReads);
        return curReadId;

}

#if USE_MPI
void Alltoallv_p2p (uint16_t **sendKmers, size_t *sendCount, uint64_t *curSendOfs, size_t *rcvCount, uint64_t *rcvOfs)
{
	int maxSendSize = 999999999*2/numTasks;

//	int maxSendSize = 100000000;
	while (maxSendSize % 7 != 0) //send recv unit - uint16_t
		maxSendSize++;
	int commStep=0,i=0;
	MPI_Request reqs[1000];
//	MPI_Request reqs2[1000];
	MPI_Status stats[1000];
    struct timeval sTime1, eTime1;
	for (commStep=0; commStep < numTasks; commStep++)
	{
        gettimeofday (&sTime1, NULL);
		int fromProc = rank-commStep;
		if (fromProc < 0)
			fromProc += numTasks;
		int toProc = (rank+commStep) % numTasks;
		size_t totalSendSize = sendCount[toProc];
		size_t totalRcvSize = rcvCount[fromProc];
		int sendIters = totalSendSize/maxSendSize;
		if ((totalSendSize % maxSendSize) != 0)
			sendIters++;
		int recvIters = totalRcvSize/maxSendSize;
		if ((totalRcvSize % maxSendSize) != 0)
			recvIters++;
//		if (recvIters > 1000)
//			printf ("[%d] RecvSize %lu Iters %d\n", rank, totalRcvSize, recvIters);
		assert (recvIters <= 1000);
		size_t rcvSoFar=0;
		uint16_t *curRcvPtr=&rcvKmers[rcvOfs[fromProc]];
		uint16_t *curSendPtr=sendKmers[numThreads*toProc];
		if (commStep == 0)
		{
			assert ((fromProc == toProc) && (totalSendSize == totalRcvSize));
			memcpy (curRcvPtr, curSendPtr, totalSendSize*sizeof(uint16_t));
      		gettimeofday (&eTime1, NULL);
            long elapsed =(eTime1.tv_sec * 1000000 + eTime1.tv_usec) - (sTime1.tv_sec * 1000000 + sTime1.tv_usec);
//		    printf ("[%d] CommStep %d - TotalSendSize=%lu bytes, TotalRcvSize=%lu Bytes %lf sec\n", rank, commStep, totalSendSize*sizeof(uint16_t), totalRcvSize*sizeof(uint16_t), (double)elapsed/1000000);
			continue;
		}
		for (i=0; i<recvIters; i++)
		{
			int curRecv=0;
			size_t remToRcv = totalRcvSize-rcvSoFar;
			if (remToRcv > (size_t)maxSendSize)
				curRecv=maxSendSize;
			else
				curRecv=(int)remToRcv;
//            printf ("[%d] p2p curRcvPtr %p curRecv %d from %d\n", rank, curRcvPtr, curRecv, fromProc);
			MPI_Irecv(curRcvPtr, curRecv, MPI_UNSIGNED_SHORT, fromProc, i, MPI_COMM_WORLD, &reqs[i]);
			rcvSoFar += curRecv;
			curRcvPtr+= curRecv;
		}
		assert (rcvSoFar == totalRcvSize);
		size_t sendSoFar=0;
		for (i=0; i<sendIters; i++)
		{
			int curSend=0;
			size_t remToSend = totalSendSize-sendSoFar;
			if (remToSend > (size_t)maxSendSize)
				curSend = maxSendSize;
			else
				curSend = (int)remToSend;
//            printf ("[%d] p2p curSendPtr %p curSend %d toProc %d\n", rank, curSendPtr, curSend, toProc);
			MPI_Send (curSendPtr, curSend, MPI_UNSIGNED_SHORT, toProc, i, MPI_COMM_WORLD);
			curSendPtr += curSend;
			sendSoFar+= curSend;
		}
		assert (sendSoFar == totalSendSize);
		MPI_Waitall (recvIters, reqs, stats);
  		gettimeofday (&eTime1, NULL);
        long elapsed =(eTime1.tv_sec * 1000000 + eTime1.tv_usec) - (sTime1.tv_sec * 1000000 + sTime1.tv_usec);
//		printf ("[%d] CommStep %d - TotalSendSize=%lu bytes, TotalRcvSize=%lu Bytes %lf sec\n", rank, commStep, totalSendSize*sizeof(uint32_t), totalRcvSize*sizeof(uint32_t), (double)elapsed/1000000);
//		MPI_Waitall (vIters, reqs, stats);
	}
}

void Alltoallv (uint16_t **sendKmers, size_t *sendCount, uint64_t *curSendOfs, size_t *rcvCount, uint64_t *rcvOfs)
{
	struct timeval sTime, eTime;
	int maxSendSize = 999999999*2/numTasks;
    while (maxSendSize % 3 != 0)
            maxSendSize++;
	int curIter=0, toSendRcv=1,i=0, j=0, totalIters=0, curRcvProc;
	size_t *rcvSoFar=(size_t *)calloc(numTasks, sizeof(size_t));
	assert (rcvSoFar != NULL);
	int *rcvCurIter=(int *)malloc(numTasks * sizeof(int));
	assert (rcvCurIter != NULL);
	int *rcvCurIterOfs=(int *)malloc(numTasks * sizeof(int));
	assert (rcvCurIterOfs != NULL);
	uint16_t *rcvBuffer=NULL;
	size_t totalRcvd=0;
	size_t *iterRcvOfs = NULL;
	for (curRcvProc=0; curRcvProc<numTasks; curRcvProc++)
	{
		int numIters = sendCount[curRcvProc]/maxSendSize;
		if ((sendCount[curRcvProc] % maxSendSize) != 0)
			numIters++;
		MPI_Allreduce (MPI_IN_PLACE, &numIters, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		if (rank==0)
			printf ("RcvRank %d TotalCommIters %d\n", curRcvProc, numIters);
		if (rank == curRcvProc)
		{
			totalIters=numIters;
			iterRcvOfs = (size_t *)calloc(numIters*numTasks + 2, sizeof(size_t));
			assert (iterRcvOfs != NULL);
		}
		size_t sentSoFar=0;
		for (i=0; i<numIters; i++)
		{
			size_t remaining = sendCount[curRcvProc]-sentSoFar;
			int sendCountCurIter=0;
			if (remaining > maxSendSize)
				sendCountCurIter = maxSendSize;
			else
				sendCountCurIter = (int) remaining;
			size_t totRcvCopy=0;	
			if (rank == curRcvProc)
			{
				rcvBuffer = &rcvKmers[totalRcvd];
				totRcvCopy=totalRcvd;
				for (j=0; j<numTasks; j++)
				{
					remaining = rcvCount[j]-rcvSoFar[j];
					if (remaining > maxSendSize)
					{
						rcvCurIter[j]=maxSendSize;
					}
					else
						rcvCurIter[j]=(int)remaining;
					totalRcvd+=rcvCurIter[j];
					rcvSoFar[j]+=rcvCurIter[j];
				}
				rcvCurIterOfs[0]=0;
				for (j=1; j<numTasks; j++)
					rcvCurIterOfs[j]=rcvCurIterOfs[j-1]+rcvCurIter[j-1];
			}
			if (rank == curRcvProc)
			{
				printf ("[%d] Bef rcvKmers[%lu]=%llx\n", rank, totRcvCopy, *((uint64_t *)&rcvKmers[totRcvCopy]));
				for (j=0; j<numTasks; j++)
				{
					printf ("Rank %d: Iter %d: Rcv %d at %p from %d\n", rank, i, rcvCurIter[j], &rcvKmers[totRcvCopy+rcvCurIterOfs[j]], j); 
					iterRcvOfs[i*numTasks + j + 2] = rcvCurIter[j];
				}
			}
			printf ("Rank %d: Iter %d: Send %d at %p(%llx) to %d\n", rank, i, sendCountCurIter, &(sendKmers[curRcvProc*numThreads][sentSoFar]), *((uint64_t *)&commonSendBuffer[curSendOfs[curRcvProc]+sentSoFar]), curRcvProc);
			int status=MPI_Gatherv (&(sendKmers[curRcvProc*numThreads][sentSoFar]), sendCountCurIter, MPI_UNSIGNED, &rcvKmers[totRcvCopy], rcvCurIter, rcvCurIterOfs, MPI_UNSIGNED, curRcvProc, MPI_COMM_WORLD);
			if (rank == curRcvProc)
				printf ("[%d] Aft rcvKmers[%lu]=%llx\n", rank, totRcvCopy, *((uint64_t *)&rcvKmers[totRcvCopy]));
//			MPI_Barrier (MPI_COMM_WORLD);
			printf ("rank %d, status %d\n", rank, status);
			sentSoFar += sendCountCurIter;
		}
	}
//	printf ("TotalRecvd %lu\n", totalRcvd);
	for (j=1; j<=totalIters*numTasks; j++)
	{
		iterRcvOfs[j]=iterRcvOfs[j-1]+iterRcvOfs[j+1];
	}
	gettimeofday (&sTime, NULL);
	size_t curOutOfs=0;
	for (i=0; i<numTasks; i++)
	{
		for (j=0; j<totalIters; j++)
		{
			size_t sizeToCopy = iterRcvOfs[j*numTasks + i +1] - iterRcvOfs[j*numTasks + i];
			size_t from = iterRcvOfs[j*numTasks + i];
			memcpy (&sortedKmers[curOutOfs], &rcvKmers[from], sizeToCopy*sizeof(uint32_t));
			printf ("[%d] memcpy %lu from ofs %lu to ofs %lu\n", rank, sizeToCopy, from, curOutOfs);
			curOutOfs += sizeToCopy;
		}
	}
	gettimeofday (&eTime, NULL);
	memcpyTime +=(eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	uint16_t *tmp = rcvKmers;
	rcvKmers = sortedKmers;
	sortedKmers = tmp;
	numKmers = totalRcvd / 3;
	free (iterRcvOfs);
	free (rcvSoFar);
	free (rcvCurIter);
	free (rcvCurIterOfs);
}
#endif

int getNumRecvIters (size_t *rcvThreadCount, int maxSendSize)
{
	int numIters = 0;
	int i=0;
	for (i=0; i<numTasks*numThreads; i++)
	{
//		printf ("[%d] rcvThreadCount[%d]= %lu, numIters %d\n", rank, i, rcvThreadCount[i], numIters);
		if ((i/numThreads) == rank)
			continue;
		numIters += rcvThreadCount[i]/maxSendSize;
		if ((rcvThreadCount[i] % maxSendSize) != 0)
			numIters ++;
	}
	return numIters;
}

#if USE_MPI
void issueRecvs (size_t *rcvThreadCount, size_t *rcvThreadOfs, int maxSendSize, MPI_Request *reqs, int totalIters)
{
    printf ("[%d]Before Issue\n", rank);
	int numIters = 0;
	int rcvNo=0;
	int i=0,j=0;
//    MPI_Request *reqs = (MPI_Request *)malloc(totalIters *sizeof(MPI_Request));
//    assert (reqs != NULL);
	for (i=0; i<numTasks*numThreads; i++)
	{
		if ((i/numThreads) == rank)
			continue;
		numIters = rcvThreadCount[i]/maxSendSize;
		if ((rcvThreadCount[i] % maxSendSize) != 0)
			numIters ++;
		size_t rcvSoFar=0;
		int curRecv=0;
		uint16_t *curRcvPtr = &rcvKmers[rcvThreadOfs[i]*3];
        printf ("[%d] rcvThreadOfs[%d]=%lu\n", rank, i, rcvThreadOfs[i]);
		for (j=0; j<numIters; j++)
                {
                        int curRecv=0;
                        size_t remToRcv = rcvThreadCount[i]-rcvSoFar;
                        if (remToRcv > (size_t)maxSendSize)
                                curRecv=maxSendSize;
                        else
                                curRecv=(int)remToRcv;
                        printf ("[%d] Ptr %p Rcv %d(%p) from %d tag %d, rcvNo %d\n", rank, curRcvPtr, curRecv, curRcvPtr+curRecv*3, i/numThreads, i%numThreads+rank*numThreads+j*numTasks*numThreads, rcvNo);
                        MPI_Irecv(curRcvPtr, curRecv*3, MPI_UNSIGNED, i/numThreads, (i%numThreads)+(rank*numThreads)+(j*numTasks*numThreads), MPI_COMM_WORLD, &reqs[rcvNo]);
                        rcvSoFar += curRecv;
                        curRcvPtr+= curRecv*3;
			rcvNo++;
                }		
	}
//    *reqs1 = reqs;
	printf ("[%d]rcvNo %d\n", rank, rcvNo);
    printf ("[%d]After Issue\n", rank);
}
void *asyncSendRcv (void *arg)
{
	int i=0;
	asyncData *ad = (asyncData *)arg;
//	size_t *sendThreadOfs = ad->sendThreadOfs;
	size_t *sendThreadCount = ad->sendThreadCount;
	size_t *threadGenerateCount = ad->threadGenerateCount;
	uint32_t **sendKmers = ad->sendKmers;
	int nbytes=0;
	size_t *sentSoFar = (size_t *)calloc(numThreads*numTasks, sizeof(size_t));
	assert (sentSoFar != NULL);

	size_t *rcvThreadOfs = ad->rcvThreadOfs;
	size_t *rcvThreadCount = ad->rcvThreadCount;
	size_t *rcvSoFar = (size_t *)calloc(numThreads*numTasks, sizeof(size_t));
	assert (rcvSoFar != NULL);

/*	printf ("[%d] AsyncThread sendThreadOfs(Count):\n", rank);
	for (i=0; i<numThreads*numTasks; i++)
	{
		printf ("%lu(%lu)\t", sendThreadOfs[i], sendThreadCount[i]*3);
	}
	printf ("\n");
	printf ("[%d] AsyncThread rcvThreadOfs(Count):\n", rank);
	for (i=0; i<numThreads*numTasks; i++)
	{
		printf ("%lu(%lu)\t", rcvThreadOfs[i]*3, rcvThreadCount[i]*3);
	}
	printf ("\n");*/

	int toSend=0, toRecv=1, flag=0;
	MPI_Request request;
	MPI_Status status;
	while (toSend || toRecv)
	{
//		printf ("[%d]ToSend %d ToRecv %d\n", rank, toSend, toRecv);
		toSend=0;
		toRecv=0;
/*		for (i=0; i<numTasks*numThreads; i++)
		{
			if (sentSoFar[i] < sendThreadCount[i])
			{
				toSend=1;
				if (((sentSoFar[i]+100000000) < threadGenerateCount[i]) || (threadGenerateCount[i] == sendThreadCount[i]))
				{
					size_t kmersToSend=threadGenerateCount[i]-sentSoFar[i];
					if (kmersToSend > 100000000)
						kmersToSend=100000000;
					if (i/numThreads == rank)
					{
						size_t rcvOfs = (rcvThreadOfs[i]+rcvSoFar[i])*3;
						memcpy (&rcvKmers[rcvOfs], &(sendKmers[i][sentSoFar[i]*3]), kmersToSend*3*sizeof(uint32_t));
//						printf ("[%d] Cpy %lu bytes from i%d Ofs %lu to %lu\n", rank, kmersToSend*3*sizeof(uint32_t), i, sentSoFar[i]*3, rcvOfs);
						sentSoFar[i] += kmersToSend;
						rcvSoFar[i] += kmersToSend;
						
					}
					else
					{
						MPI_Isend (&(sendKmers[i][sentSoFar[i]*3]), (int)kmersToSend*3, MPI_UNSIGNED, i/numThreads, i%numThreads, MPI_COMM_WORLD, &request); 
						sentSoFar[i] += kmersToSend;
//						printf ("[%d] Send %d Kmers at i %d ofs %lu to %d(%d) SentSoFar %lu/%lu\n", rank, (int)kmersToSend, i, (sentSoFar[i]-kmersToSend)*3, i/numThreads, i%numThreads, sentSoFar[i], sendThreadCount[i]);
					}
				}	
			}
		}*/
		for (i=0; i<numTasks*numThreads; i++)
		{
			if (rcvSoFar[i] < rcvThreadCount[i])
				toRecv = 1;
		}
		if (toRecv)
		{
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
			if (flag)
			{
				MPI_Get_count(&status, MPI_UNSIGNED, &nbytes);
//				printf ("[%d] IRecv(tag %d) %lu kmers (at Ofs %lu/%lu) from %dP%dT\n", rank, status.MPI_TAG, nbytes/3, rcvSoFar[status.MPI_SOURCE*numThreads+status.MPI_TAG], rcvThreadCount[status.MPI_SOURCE*numThreads+status.MPI_TAG], status.MPI_SOURCE, status.MPI_TAG);
				size_t rcvOfs = (rcvThreadOfs[status.MPI_SOURCE*numThreads+status.MPI_TAG]+rcvSoFar[status.MPI_SOURCE*numThreads+status.MPI_TAG])*3;
			
				MPI_Recv (&rcvKmers[rcvOfs], nbytes, MPI_UNSIGNED, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
				rcvSoFar[status.MPI_SOURCE*numThreads+status.MPI_TAG] += nbytes/3;
//				printf ("[%d] Recv %d Kmers at %lu from %d(Tid %d) RcvSoFar %lu/%lu\n", rank, nbytes/3, rcvOfs, status.MPI_SOURCE, status.MPI_TAG, rcvSoFar[status.MPI_SOURCE*numThreads+status.MPI_TAG], rcvThreadCount[status.MPI_SOURCE*numThreads+status.MPI_TAG]);
			}
		}
	} 
	printf ("[%d] Async Comm Thread Out of While\n", rank);
	free (sentSoFar);
	free (rcvSoFar);
}
#endif

static int compSizeComp1 (const void *p1, const void *p2)
{
        uint32_t *k1 = (uint32_t *)p1;
        uint32_t *k2 = (uint32_t *)p2;
        if (k1[1] < k2[1])
                return -1;
        else if (k1[1] > k2[2])
                return 1;
        else
                return 0;
}
static int compSizeComp2 (const void *p1, const void *p2)
{
        uint32_t *k1 = (uint32_t *)p1;
        uint32_t *k2 = (uint32_t *)p2;
        if (k1[0] < k2[0])
                return -1;
        else if (k1[0] > k2[0])
                return 1;
        else
                return 0;
}
static int edgeComp (const void *p1, const void *p2)
{
        uint32_t *k1 = (uint32_t *)p1;
        uint32_t *k2 = (uint32_t *)p2;
        if ((k1[0] < k2[0]) || ((k1[0]==k2[0]) && (k1[1] < k2[1])))
                return -1;
        else if ((k1[0] > k2[0]) || ((k1[0]==k2[0]) && (k1[1] > k1[1])))
                return 1;
        else
                return 0;
}
static int edgeCountComp (const void *p1, const void *p2)
{
        uint32_t *k1 = (uint32_t *)p1;
        uint32_t *k2 = (uint32_t *)p2;
        if (k1[2] > k2[2])
                return -1;
        else if (k1[2] < k2[2])
                return 1;
        else
                return 0;
}
#define OUTBUF_SIZE 10000000
void writeOutput3(uint64_t *partitions, int numOutPartitions)
{
	printf ("writeOutput3: Rank %d\n", rank);
	struct timeval sTime, eTime;
        gettimeofday (&sTime, NULL);
	int numP=numOutPartitions;
	char fName[100];
	FILE **fpOut = (FILE **)malloc (numP * numThreads * sizeof(FILE *));
	assert (fpOut != NULL);
	int j=0;
	for (j=0;j<numP*numThreads; j++)
	{
		int pid = j/numThreads;
		int tid = j%numThreads;

//		sprintf (fName, "%s/Fk27_ToRank%d_From%dTid%d.paired.fastq", oPrefix, pid, rank, tid);
		sprintf (fName, "%s/Rank%d/FromTid%d_ToRank%d.paired.fastq", oPrefix, rank, tid, pid);
		fpOut[j] = fopen (fName, "w");
		assert (fpOut[j] != NULL);
	}
//	int tid=0;
//	for (tid=0; tid<numThreads; tid++)
	#pragma omp parallel num_threads(numThreads)
	{
/*		char *outP = (char *)malloc(OUTBUF_SIZE);
		char *outU = (char *)malloc(OUTBUF_SIZE);
		assert ((outP != NULL) && (outU != NULL));
		size_t outP_n=0, outU_n=0;*/
	
		int tid = omp_get_thread_num();
		uint32_t i=0;
	       	char *read1=NULL, *read2=NULL;
	    	char *line1=NULL,*line2=NULL;
	    	size_t fqSize[2]={0,0}, fqInd[2]={0,0}, interleaveStart=0;
	    	int len1=0, len2=0;
	    	char *fqBuffer[2];
	    	fqBuffer[0]=&fastqBuffer1[(size_t)tid*(fastqBufferSize+10)];//(char *)malloc(fastqBufferSize+10);
		if (fastqBuffer2)
		    	fqBuffer[1] = &fastqBuffer2[(size_t)tid*(fastqBufferSize+10)];//(char *)malloc(fastqBufferSize+10);
	    	int interleave=0;
	    	uint64_t curReadId=0, pairReadId=0;

    		for (i=fastqRanges[tid]; i<fastqRanges[tid+1]; i++)
    		{
    			curReadId = (uint64_t)fI[i].readId;
			curReadId += (uint64_t)fI[i].readIdExtra * UINT32_MAX;
    			loadChunk (i, fqInd, fqSize, fqBuffer, &interleave, &interleaveStart, tid);
		        int ret = getNextFullRead (fqInd, fqSize, fqBuffer, &read1, &read2, &len1, &len2, interleave, interleaveStart);
		        while (ret != -1)
	        	{
	        	        if (ret == 2)
				{
					uint64_t p_u = (uint64_t)(*((uint32_t *)&parent[curReadId*5 + 1])) + (uint64_t)parent[curReadId*5]*UINT32_MAX;
					if (curReadId%2 == 0)
						pairReadId = curReadId+1;
					else
						pairReadId = curReadId-1;
					
	/*				if (curReadId < numReads/2)
						pairReadId = curReadId+numReads/2;
					else
						pairReadId = curReadId-numReads/2;*/
					uint64_t p_v = (uint64_t)(*((uint32_t *)&parent[pairReadId*5 +1])) + (uint64_t)parent[pairReadId*5]*UINT32_MAX;
					int process_u = 0, process_v = 0;
					for (process_u = 0; process_u < numP; process_u++)
						if (p_u < partitions[process_u+1]) break;
					for (process_v = 0; process_v < numP; process_v++)
						if (p_v < partitions[process_v+1]) break;
					if ((p_u < partitions[process_u]) || (p_u >= partitions[process_u+1]) || (p_v < partitions[process_v]) || (p_v >= partitions[process_v+1]))
						printf ("p_u %lu p_v %lu process_u %d(%lu,%lu) process_v %d(%lu,%lu) \n", p_u, p_v, process_u, partitions[process_u], partitions[process_u+1], process_v, partitions[process_v], partitions[process_v+1]);
					assert ((p_u >= partitions[process_u]) && (p_u < partitions[process_u+1]));
					assert ((p_v >= partitions[process_v]) && (p_v < partitions[process_v+1]));
	
						fwrite (read1, 1, len1, fpOut[process_u*numThreads+tid]);
						fwrite (read2, 1, len2, fpOut[process_u*numThreads+tid]);
						if (p_u != p_v)
						{
							fwrite (read1, 1, len1, fpOut[process_v*numThreads+tid]);
							fwrite (read2, 1, len2, fpOut[process_v*numThreads+tid]);
						}
				}
				else if (ret == 1)
				{
					uint64_t p_u = (uint64_t)(*((uint32_t *)&parent[curReadId*5 + 1])) + (uint64_t)parent[curReadId*5]*UINT32_MAX;
					int process_u = 0;
					for (process_u = 0; process_u < numP; process_u++)
						if (p_u < partitions[process_u+1]) break;
					if ((p_u < partitions[process_u]) || (p_u >= partitions[process_u+1]))
						printf ("p_u %lu process_u %d(%lu,%lu) \n", p_u, process_u, partitions[process_u], partitions[process_u+1]);
					assert ((p_u >= partitions[process_u]) && (p_u < partitions[process_u+1]));
					fwrite (read1, 1, len1, fpOut[process_u*numThreads+tid]);
				}
				curReadId++;
				if (ret == 2) curReadId++;
	        	        ret = getNextFullRead (fqInd, fqSize, fqBuffer, &read1, &read2, &len1, &len2, interleave, interleaveStart);
				
		        }
		}
		printf ("curReadId %u\n", curReadId);
//		free (outP);
	//	free (outU);
	}

	for (j=0; j<numP*numThreads; j++)
	{
		fclose (fpOut[j]);
	}
	free (fpOut);
        gettimeofday (&eTime, NULL);
	outputTime += (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
}

/*void writeOutput(uint64_t maxCompId)
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
	int i=0, j=0;
//#define OUTBUF_SIZE 100000000
#if USE_OMP
	#pragma omp parallel private(i) num_threads (numThreads)
	{
		int tid = omp_get_thread_num();
		char fName[100];
         	sprintf (fName, "%s/Fk27_Rank%d_Tid%d.paired_LC.fastq", oPrefix, rank, tid);
		FILE *fp_paired_lc = fopen(fName, "w");
		assert (fp_paired_lc != NULL);
         	sprintf (fName, "%s/Fk27_Rank%d_Tid%d.paired_other.fastq", oPrefix, rank, tid);
		FILE *fp_paired_other = fopen(fName, "w");
		assert (fp_paired_other != NULL);
         	sprintf (fName, "%s/Fk27_Rank%d_Tid%d.unpaired_LC.fastq", oPrefix, rank, tid);
		FILE *fp_unpaired = fopen(fName, "w");
		assert (fp_unpaired != NULL);

		char *outP = (char *)malloc(OUTBUF_SIZE);
		char *outU = (char *)malloc(OUTBUF_SIZE);
		assert ((outP != NULL) && (outU != NULL));
		size_t outP_n=0, outU_n=0;

        	char *read1=NULL, *read2=NULL;
    #else 
    		int tid = 0;
    #endif
    		char *line1=NULL,*line2=NULL;
    		size_t fqSize[2]={0,0}, fqInd[2]={0,0}, interleaveStart=0;
    		int len1=0, len2=0;
    		char *fqBuffer[2];
    		fqBuffer[0]=&fastqBuffer1[(size_t)tid*(fastqBufferSize+10)];//(char *)malloc(fastqBufferSize+10);
    		fqBuffer[1] = &fastqBuffer2[(size_t)tid*(fastqBufferSize+10)];//(char *)malloc(fastqBufferSize+10);
    		int interleave=0;
    		uint64_t curReadId=0;
    		for (i=fastqRanges[tid]; i<fastqRanges[tid+1]; i++)
    		{
    			curReadId = fI[i].readId;
    			loadChunk (i, fqInd, fqSize, fqBuffer, &interleave, &interleaveStart, tid);
		        int ret = getNextFullRead (fqInd, fqSize, fqBuffer, &read1, &read2, &len1, &len2, interleave, interleaveStart);
		        while (ret != -1)
	        	{
				uint64_t curCompId = (uint64_t)(*((uint32_t *)&parent[curReadId*5 +1])) + (uint64_t)parent[curReadId*5]*UINT32_MAX;
				if (curCompId==maxCompId)
				{
		                	if (ret == 0)
			                {
						if ((outP_n+len1+len2) >= OUTBUF_SIZE)
						{
							{
								fwrite (outP, 1, outP_n, fp_paired_lc);
							}
							outP_n=0;
						}
					      memcpy (&outP[outP_n], read1, len1);
						outP_n+=len1;
					      memcpy (&outP[outP_n], read2, len2);
						outP_n+=len2;
			                }
		        	        if (ret == 1) 
					{
						if ((outU_n+len1) >= OUTBUF_SIZE)
						{
							{
								fwrite (outU, 1, outU_n, fp_unpaired);
							}
							outU_n=0;
						}
					      memcpy (&outU[outU_n], read1, len1);
						outU_n+=len1;
					}
		                	if (ret == 2) 
					{
						if ((outP_n+len1+len2) >= OUTBUF_SIZE)
						{
							{
								fwrite (outP, 1, outP_n, fp_paired_lc);
							}
							outP_n=0;
						}
					      memcpy (&outP[outP_n], read1, len1);
						outP_n+=len1;
					      memcpy (&outP[outP_n], read2, len2);
						outP_n+=len2;
					}
				}
				else
				{
//					if (ret != 1)
//					{
//						fwrite (read1, 1, len1, fp_paired_other);
//						fwrite (read2, 1, len2, fp_paired_other);
//					}
//					else
//					{
//						printf ("unpaired other\n");
//						assert (0);
//					}
				}
				if (ret == 2) curReadId++;
				curReadId++;
	        	        ret = getNextFullRead (fqInd, fqSize, fqBuffer, &read1, &read2, &len1, &len2, interleave, interleaveStart);
		        }
    		}
		{
			if (outP_n > 0)
				fwrite (outP, 1, outP_n, fp_paired_lc);
			if (outU_n > 0)
				fwrite (outU, 1, outU_n, fp_unpaired);
		}
		free (outP);
		free (outU);
    		fclose (fp_paired_lc);
    		fclose (fp_paired_other);
	    fclose (fp_unpaired);
#if USE_OMP
	}
#endif
	gettimeofday (&eTime, NULL);
//	kmerGenTime +=(eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
}*/
void generateKmersInScan(int scanNo)
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
	int i=0, j=0;
	size_t *sendThreadCount = (size_t *)malloc(numTasks*numThreads*sizeof(size_t));
	assert (sendThreadCount != NULL);
	recvThreadCount = (size_t *)malloc((numTasks*numThreads + 2)*sizeof(size_t));
	assert (recvThreadCount != NULL);
	recvThreadOfs = (size_t *)malloc((numTasks*numThreads + 2)*sizeof(size_t));
	assert (recvThreadOfs != NULL);
	size_t *sendCount = (size_t *)malloc(numTasks*sizeof(size_t));
	assert (sendCount != NULL);
	size_t *rcvCount = (size_t *)malloc(numTasks * sizeof(size_t));
	assert (rcvCount != NULL);
	uint16_t **sendKmers = (uint16_t **)malloc(numTasks*numThreads * sizeof(uint16_t *));
	assert (sendKmers != NULL);
	uint64_t *curSendOfs = (uint64_t *)calloc (numTasks*numThreads, sizeof(uint64_t));
	assert (curSendOfs != NULL);
//	uint64_t *sendThreadOfs = (uint64_t *)calloc (numTasks*numThreads, sizeof(uint64_t));
//	assert (sendThreadOfs != NULL);
//    printf ("[%d] sendThreadOfs %p\n", rank, sendThreadOfs);
	int *curSendIter = (int *)calloc (numTasks*numThreads, sizeof(int));
	assert (curSendIter != NULL);

	size_t totalSendSize=0;
		
	for (i=0; i<numTasks; i++)
	{
		sendCount[i] = sendSizes[scanNo*numTasks + i];
		totalSendSize += sendCount[i];
//		printf ("Scan %d ToProc %d TotalSendSize %lu\t", scanNo, i, sendCount[i]);
		for (j=0; j<numThreads; j++)
		{
			sendThreadCount[i*numThreads + j] = sendThreadSizes[scanNo*numTasks*numThreads + i*numThreads + j];
//			printf ("%lu + ", sendThreadCount[i*numThreads + j]);
		}
	}
#if USE_MPI
	MPI_Alltoall (sendThreadCount, numThreads, MPI_UNSIGNED_LONG, recvThreadCount, numThreads, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
	memcpy (&recvThreadOfs[2], recvThreadCount, numThreads*numTasks*sizeof(size_t));
#else
	memcpy (&recvThreadOfs[2], sendThreadCount, numThreads*sizeof(size_t));
	memcpy (recvThreadCount, sendThreadCount, numThreads*sizeof(size_t));
#endif
	
	recvThreadOfs[0]=0;
	for (i=1; i<=numTasks * numThreads; i++)
	{
		recvThreadOfs[i]=recvThreadOfs[i-1]+recvThreadOfs[i+1];
	}
	totalSendSize *= 7* sizeof(uint16_t);
//	printf ("totalSendSize %lu\n", totalSendSize);
	assert (commonSendBuffer != NULL);
	uint64_t curOfs=0;
	for (i=0; i<numTasks*numThreads; i++)
	{
		sendKmers[i] = &commonSendBuffer[curOfs *7];
//		printf ("[%d] curSendOfs[%d]=%lu(%p)\n", rank, i, curOfs, sendKmers[i]);
//		sendThreadOfs[i] = curOfs *3;
		curOfs += sendThreadCount[i];
	}

#if USE_MPI
	MPI_Alltoall (sendCount, 1, MPI_UNSIGNED_LONG, rcvCount, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
	rcvCount[0] = sendCount[0];
#endif

/*	if (rank==0)
		printf ("Scan Range [%lu,%lu)\n", scanKmerRanges[scanNo], scanKmerRanges[scanNo+1]);*/
	uint64_t *scanProcessKmerRange = &processKmerRanges[scanNo*numTasks]; 
	int curTask=0;
	size_t totalRcv=0;
	for (i=0; i<numTasks; i++)
	{
		totalRcv += rcvCount[i];
	}
	bzero (threadKmerGenTime, numThreads*sizeof(long));
	size_t *lastKmerAddedByThread = (size_t *)calloc(numThreads*numTasks , sizeof(size_t));
	assert (lastKmerAddedByThread != NULL);
	numKmers = totalRcv;
	assert (rcvKmers!=NULL);
//	printf ("[%d] rcvKmers %lu bytes\n", rank, totalRcv * 3 * sizeof(uint32_t));
	int numOut=0;

/*	asyncData ad;
	ad.sendKmers = sendKmers;
//	ad.sendThreadOfs = sendThreadOfs;
	ad.sendThreadCount = sendThreadCount;
	ad.rcvThreadOfs = recvThreadOfs;
	ad.rcvThreadCount = recvThreadCount;
	ad.threadGenerateCount = curSendOfs;*/
	int numRcvIters = getNumRecvIters(recvThreadCount, ASYNC_SIZE);
//	printf ("NumRcvIters %d\n", numRcvIters);
#if USE_MPI
	MPI_Request *reqs = (MPI_Request *)malloc(numRcvIters * sizeof(MPI_Request));
	assert (reqs != NULL);
	MPI_Status *stats = (MPI_Status *)malloc(numRcvIters * sizeof(MPI_Status));
	assert (stats != NULL);
//	pthread_t thr;
	#if NONBLOCKING_KE
	issueRecvs (recvThreadCount, recvThreadOfs, ASYNC_SIZE, reqs, numRcvIters);
    MPI_Barrier (MPI_COMM_WORLD);
	#endif
#endif
/*	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);*/
//	if (numTasks > 1)
//		pthread_create(&thr, &attr, asyncSendRcv, &ad);
#if USE_OMP
	#pragma omp parallel private(i, j) num_threads (numThreads)
	{
		int tid = omp_get_thread_num();
         
/*        if (tid != numThreads)
        {*/
//            if (tid == 0)
//               printf ("numThreads spawned %d\n", omp_get_num_threads());
    #else 
    		int tid = 0;
    #endif
    	    struct timeval sTime1, eTime1;
    		char *line1=NULL,*line2=NULL;
    		size_t fqSize[2]={0,0}, fqInd[2]={0,0}, interleaveStart=0;
    		int len1=0, len2=0;
    		char *fqBuffer[2];
    		fqBuffer[0]=&fastqBuffer1[(size_t)tid*(fastqBufferSize+10)];//(char *)malloc(fastqBufferSize+10);
    //		assert (fqBuffer[0] != NULL);
		if (fastqBuffer2 != NULL)
	    		fqBuffer[1] = &fastqBuffer2[(size_t)tid*(fastqBufferSize+10)];//(char *)malloc(fastqBufferSize+10);
		else
			fqBuffer[1] = NULL;
    //		assert (fqBuffer[1] != NULL);
//            if ((scanNo==0) && (rank ==0))
//                printf ("tid %d fqBuffer_0 = %p, fqBuffer_1 %p\n", tid, fqBuffer[0], fqBuffer[1]);
    		int interleave=0;
    //		openFp[2*tid] = openFp[2*tid +1] = NULL;
    //		openFileNo[2*tid] = openFileNo[2*tid + 1] = -1;
    		uint64_t curReadId=0;
            size_t cpySoFar=0;
    		for (i=fastqRanges[tid]; i<fastqRanges[tid+1]; i++)
    		{
    			curReadId = (uint64_t)fI[i].readIdExtra*UINT32_MAX+fI[i].readId;
    //			printf ("Tid %d proces fq chunk %d, curReadId %u\n", tid, i, curReadId);
    			loadChunk (i, fqInd, fqSize, fqBuffer, &interleave, &interleaveStart, tid);
		        while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
		        fqInd[0]++;
		//	printf ("Chunk %d Size (%lu,%lu), interleaved? %d(Ofs %lu) Start(%d-%lu, %lu) End (%d-%lu, %lu)\n", chunk, fqSize[0], fqSize[1], *interleave, *interleaveStart, startFile, startOfs1, startOfs2, endFile, endOfs1, endOfs2);
			assert (fqInd[0] < fqSize[0]);
		        if (fqInd[1] < fqSize[1])
		        {
		                while ((fqInd[1] < fqSize[1]) && (fqBuffer[1][fqInd[1]] != '\n')) fqInd[1]++;
                		fqInd[1]++;
		                assert (fqInd[1] < fqSize[1]);
		        }
    			gettimeofday (&sTime1, NULL);
    			curReadId = getReadAndKmerCount (i, fqInd, fqSize, fqBuffer, interleave, interleaveStart, sendKmers, curSendOfs, scanKmerRanges[scanNo], scanKmerRanges[scanNo+1], kmerMask, rotater, curReadId, &lastKmerAddedByThread[rank*numThreads + tid], scanProcessKmerRange, tid);

    			gettimeofday (&eTime1, NULL);
    			threadKmerGenTime[tid] +=(eTime1.tv_sec * 1000000 + eTime1.tv_usec) - (sTime1.tv_sec * 1000000 + sTime1.tv_usec);
#if NONBLOCKING_KE
			printf ("Not modified for 5-byte kmer tuple\n");
			assert (0);
			for (j=0; j<numTasks; j++)
			{
				int curSend=0;
				if (curSendOfs[j*numThreads + tid] > ((size_t)(curSendIter[j*numThreads + tid])*ASYNC_SIZE))
				{
					size_t kmersToSend = curSendOfs[j*numThreads + tid] - (size_t)(curSendIter[j*numThreads + tid])*ASYNC_SIZE;
					uint32_t *sendPtr = sendKmers[j*numThreads + tid];
					sendPtr += (size_t)(curSendIter[j*numThreads + tid])*ASYNC_SIZE*3;
					if (kmersToSend > ASYNC_SIZE)
					{
						curSend=ASYNC_SIZE;
					}
					else if (curSendOfs[j*numThreads + tid] == sendThreadCount[j*numThreads + tid])
					{
						curSend = (int)kmersToSend;
					}
					if (curSend > 0)
					{
						MPI_Request req;
						if (j != rank)
						{
//						    printf ("[%d] Tid %d SendPtr %p Send %d(%lu) to %d Tag %d\n", rank, tid, sendPtr, curSend, sendThreadCount[j*numThreads + tid], j, curSendIter[j*numThreads + tid]*numTasks*numThreads+j*numThreads+tid);
//                            uint32_t spVal = sendPtr[0];
                            int flag=0;
//                            if (curSendIter[j*numThreads + tid] > 0)
							MPI_Isend (sendPtr, curSend*3, MPI_UNSIGNED, j, curSendIter[j*numThreads+tid]*numTasks*numThreads+j*numThreads+tid, MPI_COMM_WORLD, &req);
//                            sendPtr[1]=spVal;
//                            sendNo++; 
						}
						else
						{
                            memcpy (&rcvKmers[(recvThreadOfs[rank*numThreads + tid] + cpySoFar)*3], sendPtr, curSend*3*sizeof(uint32_t));
                            cpySoFar+=curSend;
						}
						curSendIter[j*numThreads + tid]++;
					}
				}
			}
#endif
    		}
/*        }
        else
        {
            asyncSendRcv(&ad);
        }*/
#if USE_OMP
		for (j=0; j<numTasks; j++)
			printf ("Tid %d NumKmers for %d = %lu\n", tid, j, curSendOfs[j*numThreads + tid]);
	}
#endif
	long maxWorkTime=0;
	for (i=0; i<numThreads; i++)
	{
		if (threadKmerGenTime[i] > maxWorkTime)
			maxWorkTime = threadKmerGenTime[i];
	}
	kmerGenWorkTime += maxWorkTime;
	gettimeofday (&eTime, NULL);
	kmerGenTime +=(eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	gettimeofday (&sTime, NULL);
#if USE_MPI
#if NONBLOCKING_KE
//    printf ("[%d] Before Waitall\n", rank);
//    int waitIters = numRcvIters/100;
//    if ((numRcvIters % 100) != 0)
//        waitIters++;
    MPI_Waitall(numRcvIters, reqs, stats);
/*    for (i=0; i<(waitIters-1); i++)
    {
        MPI_Waitall(100, &reqs[100*i], &stats[100*i]);
    }
    if ((numRcvIters % 100) == 0)
    	MPI_Waitall(100, &reqs[100*i], &stats[100*i]);
    else
        MPI_Waitall(numRcvIters % 100, &reqs[100*i], &stats[100*i]);*/
//    for (i=0; i<numRcvIters; i++)
//        printf ("[%d] RcvStat[%d]=%d\n", rank, i, stats[i].MPI_ERROR);
//    printf ("[%d] After Waitall\n", rank);
//	if (numTasks > 1)
//		pthread_join (thr, 0); 
        MPI_Barrier (MPI_COMM_WORLD); 
#else
	if (numTasks > 1)
	{
		uint64_t *rcvOfs = (uint64_t *)malloc(numTasks * sizeof(uint64_t));
		assert (rcvOfs != NULL);
		rcvOfs[0]=0; curSendOfs[0]=0;
		sendCount[0]*=7; rcvCount[0]*=7;
		for (i=1; i<numTasks; i++)
		{
			rcvOfs[i]=rcvOfs[i-1]+rcvCount[i-1];
			curSendOfs[i]=curSendOfs[i-1]+sendCount[i-1];
			sendCount[i] *= 7;
			rcvCount[i] *= 7;
		}
	
		for (i=0; i<numTasks; i++)
		{
//			printf ("[%d] Scan %d Send %lu toProcess %d at %lu, Rcv %lu from process %d at %lu\n", rank, scanNo, sendCount[i], i, curSendOfs[i], rcvCount[i], i, rcvOfs[i]);
		}
		Alltoallv_p2p (sendKmers, sendCount, curSendOfs, rcvCount, rcvOfs);
		free (rcvOfs);
	}
#endif
	free (reqs);
	free (stats);
#endif
	gettimeofday (&eTime, NULL);
	commTime +=(eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);

	free (sendKmers);
	free (lastKmerAddedByThread);
	free (curSendOfs);
	free (rcvCount);
	free (sendCount);
	free (sendThreadCount);
//	free (sendThreadOfs);
	free (curSendIter);

}

#if 0
static int kmerCompare (const void *p1, const void *p2)
{
	uint64_t k1 = *(uint64_t *)p1;
	uint64_t k2 = *(uint64_t *)p2;
	uint32_t r1 = *((uint32_t *)p1+2);
	uint32_t r2 = *((uint32_t *)p2+2);
	if ((k1 < k2) || ((k1==k2) && (r1 < r2)))
		return -1;
	else if ((k1 > k2) || ((k1 == k2) && (r1 > r2)))
		return 1;
	else
		return 0;
}
#endif

uint64_t radixPass (int pass, uint16_t *in_list, uint16_t *out_list, uint64_t startOfs, uint64_t numKmersToProcess, int scanNo, int tid)
{
	int passOrder[14]={10,11,12,13,9,8,0,1,2,3,4,5,6,7};
	uint16_t radixMask=0xFFu;
	int radixRotate=0;
	uint32_t i=0,j=0;
	radixRotate = 8 * (passOrder[pass] & 1);
	radixMask = radixMask << radixRotate;
	int uintOfs = passOrder[pass]/2;

	uint64_t curKmerMask = 0;

	uint64_t *rangeOfsCopy = (uint64_t *)calloc(256+2 , sizeof(uint64_t));
	assert (rangeOfsCopy != NULL);
//	printf ("before for 0\n");
//	printf ("before for 0\n");
	uint64_t kmerInd=0;
	int curM=0;
	uint64_t beginInd=startOfs, beginInd2=startOfs;
	for (kmerInd = startOfs; kmerInd<(startOfs+numKmersToProcess); kmerInd++)
	{
		uint16_t bucket = (in_list[kmerInd*7 + uintOfs] & radixMask) >> radixRotate;
		rangeOfsCopy[bucket+2]++;
	}

#define LMAX 100
#define LSIZE 700
	uint16_t lRadixKmers[256][LSIZE];
	int lCount[256];
//	printf ("Bef for1\n");
	for (i=1; i<=256; i++)
	{
		rangeOfsCopy[i]=rangeOfsCopy[i-1]+rangeOfsCopy[i+1];
		lCount[i-1]=0;
//		if (tid==0)
//			printf ("%lu\t", rangeOfsCopy[i-1]);
	}
//	if (tid==0)
//		printf ("\n");
	uint64_t retVal = rangeOfsCopy[256];
	for (kmerInd = startOfs; kmerInd<(startOfs+numKmersToProcess); kmerInd++)
	{
		uint16_t bucket = (in_list[kmerInd*7 + uintOfs] & radixMask) >> radixRotate;
/*		if (*((uint64_t *)&in_list[kmerInd*3]) == 0x6e6000000000000lu)
		{
			printf ("tid %d kmerind %lu(startOfs %lu) Kmer %llx ReadId %u\n", tid, kmerInd, startOfs,  *((uint64_t *)&in_list[kmerInd*3]), in_list[kmerInd*3 +2]); 
		}*/
//		if (in_list[kmerInd*5 +4] != 0xFFFFFFFF)
		{
			if (lCount[bucket] == LMAX)
			{
				uint64_t curOfs = startOfs+rangeOfsCopy[bucket];
				for (i=0; i<LMAX; i++)
				{
					out_list[curOfs*7]=lRadixKmers[bucket][i*7];
					out_list[curOfs*7 + 1]=lRadixKmers[bucket][i*7 + 1];
					out_list[curOfs*7 + 2]=lRadixKmers[bucket][i*7 + 2];
					out_list[curOfs*7 + 3]=lRadixKmers[bucket][i*7 + 3];
					out_list[curOfs*7 + 4]=lRadixKmers[bucket][i*7 + 4];
					out_list[curOfs*7 + 5]=lRadixKmers[bucket][i*7 + 5];
					out_list[curOfs*7 + 6]=lRadixKmers[bucket][i*7 + 6];
					curOfs++;
				}
				rangeOfsCopy[bucket]+=LMAX;
				lCount[bucket]=0;
			}
			int curOfs = lCount[bucket];
			lRadixKmers[bucket][curOfs*7]=in_list[kmerInd*7];
			lRadixKmers[bucket][curOfs*7 + 1]=in_list[kmerInd*7 + 1];
			lRadixKmers[bucket][curOfs*7 + 2]=in_list[kmerInd*7 + 2];
			lRadixKmers[bucket][curOfs*7 + 3]=in_list[kmerInd*7 + 3];
			lRadixKmers[bucket][curOfs*7 + 4]=in_list[kmerInd*7 + 4];
			lRadixKmers[bucket][curOfs*7 + 5]=in_list[kmerInd*7 + 5];
			lRadixKmers[bucket][curOfs*7 + 6]=in_list[kmerInd*7 + 6];
			lCount[bucket]++;
/*			out_list[(startOfs+rangeOfsCopy[bucket])*3]=in_list[kmerInd*3];
			out_list[(startOfs+rangeOfsCopy[bucket])*3 +1]=in_list[kmerInd*3 +1];
			out_list[(startOfs+rangeOfsCopy[bucket])*3 +2]=in_list[kmerInd*3 +2];
			rangeOfsCopy[bucket]++;*/
		}
	}
	for (i=0; i<256; i++)
	{
		uint64_t curOfs = startOfs+rangeOfsCopy[i];
		for (j=0; j<lCount[i]; j++)
		{
			out_list[curOfs*7]=lRadixKmers[i][j*7];
			out_list[curOfs*7 +1]=lRadixKmers[i][j*7 +1];
			out_list[curOfs*7 +2]=lRadixKmers[i][j*7 +2];
			out_list[curOfs*7 +3]=lRadixKmers[i][j*7 +3];
			out_list[curOfs*7 +4]=lRadixKmers[i][j*7 +4];
			out_list[curOfs*7 +5]=lRadixKmers[i][j*7 +5];
			out_list[curOfs*7 +6]=lRadixKmers[i][j*7 +6];
			curOfs++;
		}
		rangeOfsCopy[i]+=lCount[i];
//		if (tid==0)
//			printf ("%lu\t", rangeOfsCopy[i]);
	}
	assert (rangeOfsCopy[255] == rangeOfsCopy[256]);
//	if (tid==0)
//		printf ("\n");
	free (rangeOfsCopy);
	return retVal;
}
uint64_t radixPass1 (int pass, uint32_t *in_list, uint32_t *out_list, uint64_t startOfs, uint64_t numKmersToProcess, int scanNo, int tid)
{
	int passOrder[3]={2,0,1};
	uint32_t radixMask=0xFFFFu;
	int radixRotate=0;
	uint32_t i=0,j=0;
	if ((pass & 1) == 1)
	{
		radixMask = radixMask << 16;
		radixRotate = 16;
	}
	int uintOfs = passOrder[pass/2];

	uint64_t curKmerMask = 0;

/*	if (rank==0)
		printf ("pass %d, start %lu, endOfs %lu, radixMask %" PRIu32 " radixRotate %d\n", pass, startOfs, startOfs+numKmersToProcess, radixMask, radixRotate);*/
	uint64_t *rangeOfs = (uint64_t *)calloc(65536+2 , sizeof(uint64_t));
	assert (rangeOfs != NULL);
	uint64_t *rangeOfsCopy = (uint64_t *)calloc(65536+2 , sizeof(uint64_t));
	assert (rangeOfsCopy != NULL);
	uint64_t maxCompCount=0;
	uint32_t *rangeLastSeenComp=NULL;
	uint64_t *rangeCompCount=NULL;
	uint64_t *suffixCompCount=NULL;
	uint64_t currentSuffix=0;
#if RADIX_OPT
	if ((scanNo >= 0) && (pass > 3)) 
//	if ((scanNo == 0) && (pass == 4)) 
	{
		curKmerMask = (1lu << (16*(pass-2) - 1)) | ((1lu << (16*(pass-2)-1)) - 1);
//		if (tid == 0)
//			printf ("scan %d tid %d pass %d curKmerMAsk %llx, startReadId %u endReadId %u, startOfs %lu, endOfs %lu, startKmer %llx, endKmer %llx\n", scanNo, tid, pass, curKmerMask, in_list[startOfs*3 + 2], in_list[(startOfs+numKmersToProcess-1)*3 + 2], startOfs, startOfs+numKmersToProcess, *((uint64_t *)&in_list[startOfs*3]), *((uint64_t *)&in_list[(startOfs+numKmersToProcess-1)*3]));

		rangeLastSeenComp = (uint32_t *)malloc(65536 * sizeof(uint32_t));
		assert (rangeLastSeenComp != NULL);
		rangeCompCount = (uint64_t *)calloc(65536, sizeof(uint64_t));
		assert (rangeCompCount != NULL);
		suffixCompCount = (uint64_t *)calloc(65536, sizeof(uint64_t));
		assert (suffixCompCount != NULL);
		uint32_t bucket = (in_list[startOfs*3 + uintOfs] & radixMask) >> radixRotate;
		currentSuffix = *((uint64_t *)&in_list[startOfs*3]);
		currentSuffix = currentSuffix & curKmerMask;
		rangeCompCount[bucket] = 1;
		rangeLastSeenComp[bucket]=in_list[startOfs*3 + 2];
	}
#endif
//	printf ("before for 0\n");
	uint64_t kmerInd=0;
	uint64_t kmerRemoveCount=0;
	uint64_t kmerCount=0;
	int curM=0;
	uint64_t beginInd=startOfs, beginInd2=startOfs;
	for (kmerInd = startOfs; kmerInd<(startOfs+numKmersToProcess); kmerInd++)
	{
		uint32_t bucket = (in_list[kmerInd*3 + uintOfs] & radixMask) >> radixRotate;
/*		if ((bucket==0) && (rank==0))
		{
			printf ("[%d] pass %d Kmer %llx ReadId %lx, uintOfs %d, in_list = %lx,%lx,%lx\n", rank, pass, *((uint64_t *)&in_list[kmerInd*3]),in_list[kmerInd*3+2],  uintOfs, in_list[kmerInd*3], in_list[kmerInd*3 + 1], in_list[kmerInd*3 + 2]);
		}*/
#if RADIX_OPT
		if ((scanNo >= 0) && (pass > 3))
//		if ((scanNo == 0) && (pass == 3)) 
		{	
			uint64_t newSuffix = *((uint64_t *)&in_list[kmerInd*3]);
			newSuffix = newSuffix & curKmerMask;
			if (newSuffix != currentSuffix)
			{
//				if (tid==0)
//					printf ("Old Suffix %llx, New Suffix %llx, beginInd %lu, kmerInd %lu, kmerRemoveCount %lu\n", currentSuffix, newSuffix, beginInd, kmerInd, kmerRemoveCount);
				while (beginInd < kmerInd)
				{
					i=(in_list[beginInd*3 + uintOfs] & radixMask) >> radixRotate;
					if (rangeCompCount[i]==1)
					{
						kmerRemoveCount++;//=rangeCompCount[bucket];
						in_list[beginInd*3 +2] = 0xFFFFFFFF;
						rangeOfsCopy[i+2]--;
//						if (tid == 0)
//							printf ("Kmer %llx bucket %x(curPassMask %llx) - Single COmp %u, Actual Kmer Comp %u CompCount %lu\n", *((uint64_t *)&in_list[beginInd*3]), i, curKmerMask, rangeLastSeenComp[i], in_list[beginInd*3 + 2], suffixCompCount[i]);
						if (maxCompCount < suffixCompCount[i])
							maxCompCount=suffixCompCount[i];
					}
					beginInd++;
				}
				while (beginInd2 < kmerInd)
				{
					i=(in_list[beginInd*3 + uintOfs] & radixMask) >> radixRotate;
					rangeCompCount[i]=0;
					suffixCompCount[i]=0;
					beginInd2++;
				}
				currentSuffix = newSuffix;
				kmerCount++;
/*				if ((tid==0) && (scanNo==1) && (pass==4) && (((kmerInd-startOfs)%1000000)==0))
				{
					printf ("%lu/%lu curSuffix %llx\n", kmerInd, startOfs+numKmersToProcess, currentSuffix);
				}*/
				rangeCompCount[bucket]=1;
				suffixCompCount[bucket]=1;
				rangeLastSeenComp[bucket]=in_list[kmerInd*3 + 2];
			}
			else /*if (rangeCompCount[bucket]<2) */
			{
				suffixCompCount[bucket]++;
				if (rangeLastSeenComp[bucket] != in_list[kmerInd*3 + 2])
				{
					rangeCompCount[bucket]++;
					rangeLastSeenComp[bucket]=in_list[kmerInd*3 + 2];
				}
			}
		}
#endif
		rangeOfs[bucket+2]++;
		rangeOfsCopy[bucket+2]++;
	}

//	int lCount[65536];
//	printf ("Bef for1\n");
	for (i=1; i<=65536; i++)
	{
/*		if (pass == 2)
		{
			if (rangeCompCount[i-1]==1)
				printf ("Kmer Ending with %x Count %lu - Single COmp %u\n", i-1, rangeOfs[i+1], rangeLastSeenComp[i-1]);
		}*/
		rangeOfs[i]=rangeOfs[i-1]+rangeOfs[i+1];
		rangeOfsCopy[i]=rangeOfsCopy[i-1]+rangeOfsCopy[i+1];
//		printf ("lCount i %d\n", i);
//		lCount[i-1]=0;
	}
	uint64_t retVal = rangeOfsCopy[65536];
#if RADIX_OPT
	if ((scanNo >= 0) && (pass > 3))
//	if ((scanNo == 0) && (pass == 3)) 
		printf ("Rank %d Tid %d Scan %d Pass %d Reduce %lu kmers, uniqueKmers %lu, maxCompCount %lu, numKmersToProcess %lu, lastRangeOfsCopy %lu\n", rank, tid, scanNo, pass, kmerRemoveCount, kmerCount, maxCompCount, numKmersToProcess, rangeOfsCopy[65536]);
#endif
//	if (rank==0)
//		printf ("rangeOfs[65536]=%lu numKmersToProcess %lu\n", rangeOfs[65536], numKmersToProcess);
	for (kmerInd = startOfs; kmerInd<(startOfs+numKmersToProcess); kmerInd++)
	{
		uint32_t bucket = (in_list[kmerInd*3 + uintOfs] & radixMask) >> radixRotate;
//		if ((kmerInd % 10000000) == 0)
//			printf ("kmerINd %lu rangeOfsCopy[0]=%lu\n", kmerInd, rangeOfsCopy[0]);
		if (in_list[kmerInd*3 +2] != 0xFFFFFFFF)
		{
/*			if (lCount[bucket] == 300)
			{
				uint64_t curOfs = startOfs+rangeOfsCopy[bucket];
//				assert (bucket < 65536);
				for (i=0; i<300; i++)
				{
					printf ("Access lkmer[%d][%u] \n", bucket, i*3);
//					out_list[curOfs*3]=lRadixKmers[tid][bucket][i*3];
//					out_list[curOfs*3 + 1]=0;//lRadixKmers[tid][bucket][i*3 + 1];
//					out_list[curOfs*3 + 2]=0;//lRadixKmers[tid][bucket][i*3 + 2];
					curOfs++;
				}
				rangeOfsCopy[bucket]+=300;
				if (curOfs > numKmers)
				{
					printf ("StartOfs %lu, bucket %d, rangeOfsCopy[%d]=%lu, curOfs %lu,  numKmers %lu\n", startOfs, bucket, bucket, rangeOfsCopy[bucket], curOfs, numKmers);
				}
				assert (curOfs <= numKmers);
				lCount[bucket]=0;
			}
			int curOfs = lCount[bucket];
			printf ("bucket %d, curOfs %d\n", bucket, curOfs);
			lRadixKmers[tid][bucket][curOfs*3]=in_list[kmerInd*3];
			lRadixKmers[tid][bucket][curOfs*3 + 1]=in_list[kmerInd*3 + 1];
			lRadixKmers[tid][bucket][curOfs*3 + 2]=in_list[kmerInd*3 + 2];
			lCount[bucket]++;*/
			out_list[(startOfs+rangeOfsCopy[bucket])*3]=in_list[kmerInd*3];
			out_list[(startOfs+rangeOfsCopy[bucket])*3 +1]=in_list[kmerInd*3 +1];
			out_list[(startOfs+rangeOfsCopy[bucket])*3 +2]=in_list[kmerInd*3 +2];
			rangeOfsCopy[bucket]++;
		}

		rangeOfs[bucket]++;
	}
/*	for (i=0; i<65536; i++)
	{
		uint64_t curOfs = startOfs+rangeOfsCopy[i];
		for (j=0; j<lCount[i]; j++)
		{
			out_list[curOfs*3]=lRadixKmers[tid][i][j*3];
			out_list[curOfs*3 +1]=lRadixKmers[tid][i][j*3 +1];
			out_list[curOfs*3 +2]=lRadixKmers[tid][i][j*3 +2];
			curOfs++;
		}
		rangeOfsCopy[i]+=lCount[i];
	}*/
/*	if (rank==0)
		printf ("out_list[%lu]=%" PRIu64 ", readId %u\n", startOfs, *((uint64_t *)&out_list[startOfs*3]), out_list[startOfs*3 + 2]);
	if (rank==0)
		printf ("out_list[%lu]=%" PRIu64 ", readId %u\n", startOfs+numKmersToProcess-1, *((uint64_t *)&out_list[(startOfs+numKmersToProcess-1)*3]), out_list[(startOfs+numKmersToProcess-1)*3 + 2]);*/
	free (rangeOfs);
	free (rangeOfsCopy);
#if RADIX_OPT
	
	if ((scanNo >= 0) && (pass > 3))
//	if ((scanNo == 0) && (pass == 3)) 
	{
		free (rangeLastSeenComp);
		free (rangeCompCount);
		free (suffixCompCount);
	}
#endif
	return retVal;
}

void calcRangeCounts (int tid, size_t *threadRangeCounts, uint64_t *localThreadKmerRanges)
{
	int startThread =  tid * numTasks;
	int endThread =  (tid+1) * numTasks;
	int fastqPartitionsPerThread = totalFastqPartitions/numTasks/numThreads;
	int startFastqPartition = startThread*fastqPartitionsPerThread;
	int endFastqPartition = endThread * fastqPartitionsPerThread;
	size_t startRcvBufInd = recvThreadOfs[startThread];
	size_t endRcvBufInd = recvThreadOfs[endThread];
//	printf ("Rank %d Tid %d startFilePart %d endFilePart %d startRcvBufInd %lu, endRcvBufInd %lu\n", rank, tid, startFastqPartition, endFastqPartition, startRcvBufInd, endRcvBufInd);
/*	if (tid == 0)
	{
		int i=0;
		for (i=0; i<numTasks * numThreads; i++)
			printf ("%lu\t", recvThreadCount[i+1]-recvThreadCount[i]);
		printf ("\n");
	}*/
	int f,t,p;
//	printf ("calcRangeCounts %d : startFastqPartition %d endFastqPartition %d\n", tid, startFastqPartition, endFastqPartition);
/*	for (f=startFastqPartition; f<endFastqPartition; f++)
	{
		for (t=0; t<numThreads; t++)
		{
			for (p=localThreadKmerRanges[t]; p<localThreadKmerRanges[t+1]; p++)
			{
				threadRangeCounts[tid*numThreads + t] += fI[f].kmerFreqCount[p];
			}
		}
	}*/
//	size_t *threadRangeCountscpy = (size_t *)calloc (numThreads*numThreads, sizeof(uint64_t));
	for (t=0; t<numThreads; t++)
	{
		for (p=localThreadKmerRanges[t]; p<localThreadKmerRanges[t+1]; p++)
		{
			threadRangeCounts[tid*numThreads + t] += threadKmerCounts[tid * 1048576 + p];
		}
//		printf ("Tid %d t %d threadRangeCountscpy %lu threadRangeCounts %lu\n", tid, t, threadRangeCountscpy[tid*numThreads + t], threadRangeCounts[tid*numThreads + t]);
//		assert (threadRangeCountscpy[tid*numThreads + t] == threadRangeCounts[tid*numThreads + t]);
	}
//	free (threadRangeCountscpy);
}

void sortKmerList (int scanNo, uint64_t *localRangeCounts, uint64_t *rangeOffsets)
{
	typedef struct 
	{
		uint16_t kmertuple[7];
//		uint32_t readId;
//		int range;
	}kmerTemp;
//	printf ("[%d] processKmerRanges=[%lu,%lu)\n", rank, processKmerRanges[scanNo*numTasks+rank], processKmerRanges[scanNo*numTasks+rank + 1]);
	size_t *threadRangeCounts =(size_t *)calloc (numThreads * numThreads ,sizeof(size_t));
	assert (threadRangeCounts != NULL);
	size_t *threadRangeOfs =(size_t *)calloc (numThreads * numThreads ,sizeof(size_t));
	assert (threadRangeOfs != NULL);

	struct timeval sTime, eTime;
	int *threadMask = (int *)calloc(1048576, sizeof(int));
	assert (threadMask!=NULL);
	int f=0, i;
	uint64_t i1=0;
	uint64_t processLimitUp = processKmerRanges[scanNo*numTasks+rank + 1];
	uint64_t processLimitDn = processKmerRanges[scanNo*numTasks+rank];
//	if (rank == 0)
#ifdef DEBUG
	{
		for (i1=0; i1<numKmers; i1++)
		{
			uint64_t curKmerL = *((uint64_t *)&rcvKmers[i1 * 7]);
//			uint64_t curKmerH = *((uint64_t *)&rcvKmers[i1 * 5 + 2]);
			uint64_t curKmerH = 0;
			uint64_t bucket=getMmerPrefix(curKmerH, curKmerL, 10);
			if ((bucket >= processLimitUp) || (bucket < processLimitDn))
			{
				printf ("[%d]i1 %lu curKmer %llx,%llx bucket %lu(%p) procLimUp %lu, procLimDn %lu", rank, i1, curKmerH, curKmerL, bucket, &rcvKmers[i1 * 7], processLimitUp, processLimitDn);
				printf ("*\n");
			}
			assert ((bucket < processLimitUp) && (bucket >= processLimitDn));
		}
	}
#endif
//	uint32_t prevReadId = rcvKmers[2];
/*	if (scanNo == 0)
	{
		for (i1=1; i1<numKmers; i1++)
		{
			assert (rcvKmers[i1 * 3 + 2] >= prevReadId);
			prevReadId = rcvKmers[i1 * 3 + 2];
		}
	}*/
//	printf ("[%d] Sort %lu kmers\n", rank, numKmers);
	
	assert (sortedKmers != NULL);
	uint64_t *fastqRangeCounts = (uint64_t *)calloc(1048576, sizeof(uint64_t));
	assert (fastqRangeCounts != NULL);
	uint64_t *localThreadKmerRanges = &threadKmerRanges[scanNo*numTasks*numThreads + rank*numThreads];
	uint64_t *localRangeOffsets = (uint64_t *)calloc(numThreads, sizeof(uint64_t));
	assert (localRangeOffsets != NULL);
/*	for (f=0; f<totalFastqPartitions; f++)
	{
		int curThread=0;
		uint64_t startRange = (localThreadKmerRanges[0]);// & kmerMask) >> rotater;
		uint64_t endRange = (localThreadKmerRanges[numThreads]);// & kmerMask) >> rotater;
//		printf ("[%d] startRange %llx endRange %llx\n", rank, startRange, endRange);
		for (i1=startRange; i1<endRange; i1++)
		{
//			printf ("chrThread %d i1 %llx\n", curThread, i1);
			if (i1 >= localThreadKmerRanges[curThread+1])
			{
				curThread++;
				localRangeOffsets[curThread]=localRangeOffsets[curThread-1]+localRangeCounts[curThread-1];
				rangeOffsets[curThread]=localRangeOffsets[curThread];
			}
			threadMask[i1]=curThread;
			localRangeCounts[curThread]+=fI[f].kmerFreqCount[i1];
		}
	}
	uint64_t *tmp = (uint64_t *)calloc(numThreads, sizeof(uint64_t));
	bzero (localRangeOffsets, numThreads*sizeof(uint64_t));*/
	for (f=0; f<numThreads; f++)
	{
		int curThread=0;
		uint64_t startRange = (localThreadKmerRanges[0]);// & kmerMask) >> rotater;
		uint64_t endRange = (localThreadKmerRanges[numThreads]);// & kmerMask) >> rotater;
//		printf ("[%d] startRange %llx endRange %llx\n", rank, startRange, endRange);
		for (i1=startRange; i1<endRange; i1++)
		{
//			printf ("chrThread %d i1 %llx\n", curThread, i1);
			if (i1 >= localThreadKmerRanges[curThread+1])
			{
				curThread++;
				localRangeOffsets[curThread]=localRangeOffsets[curThread-1]+localRangeCounts[curThread-1];
				rangeOffsets[curThread]=localRangeOffsets[curThread];
			}
//			tmp[curThread]+=threadKmerCounts[f*1048576+i1];
			localRangeCounts[curThread]+=threadKmerCounts[f*1048576+i1];
		}
	}

//	for (i=0; i<numThreads; i++)
//	{
//		assert (localRangeOffsets[i] == rangeOffsets[i]);
//	}
//	free (tmp);
	for (i=0; i<numThreads; i++)
	{
//		printf ("[%d]  threadKmerRanges[%d]=[%llx,%llx) Count %lu Ofs %lu\n", rank, i, localThreadKmerRanges[i], localThreadKmerRanges[i+1], localRangeCounts[i], localRangeOffsets[i]);
		calcRangeCounts (i, threadRangeCounts, localThreadKmerRanges);
	}
//	if (rank == 0)
	{
		uint64_t curCount=0;
		for (i=0; i<numThreads; i++)//column
		{
			for (f=0; f<numThreads; f++)//row
			{
				threadRangeOfs[f*numThreads + i]=curCount;
				curCount += threadRangeCounts[f*numThreads + i];
			}	
		}

/*		for (i=0; i<numThreads; i++)
		{
			for (f=0; f<numThreads; f++)
			{
				printf ("%lu,", threadRangeOfs[i*numThreads + f]);
			}	
			printf ("\n");
		}*/
	}

	// qsort (rcvKmers, numKmers, 3*sizeof(uint32_t), kmerCompare);
	// memcpy (sortedKmers, rcvKmers, numKmers*3*sizeof(uint32_t));
	gettimeofday (&sTime, NULL);
	printf ("Before range partition\n");

#if 1
#if USE_OMP
	#pragma omp parallel num_threads(numThreads)
	{
		int tid=omp_get_thread_num();
#else
	{
		int tid=0;
#endif
		int startThread =  tid * numTasks;
		int endThread = (tid+1) * numTasks;
		size_t startRcvBufInd = recvThreadOfs[startThread];
		size_t endRcvBufInd = recvThreadOfs[endThread];
		int f=0, i;
		uint64_t j1=0;
		uint64_t prevInd=0;
#define BUFSIZE 1000
		kmerTemp lKmers[256][BUFSIZE];
		int rangeWiseCounts[256];
		assert (numThreads <= 256); 
		for (i=0; i<numThreads; i++)
			rangeWiseCounts[i]=0;

		size_t *localPtr = &threadRangeOfs[tid*numThreads];
//		printf ("RangePartition tid %d StartRcvBufInd %lu end %lu, numKmersToProcess %lu, localPtr[0]=%lu\n", tid, startRcvBufInd, endRcvBufInd, endRcvBufInd-startRcvBufInd, localPtr[0]);
		for (j1=startRcvBufInd; j1 < endRcvBufInd; j1++)
		{
//			rcvKmers[j1*5 +3]=0;
//			if (tid==0)
//				printf ("tid %d j1 %lu\n", tid, j1);
			uint64_t maskedKmerL=*(uint64_t *)(&rcvKmers[j1*7]);
			uint64_t maskedKmerH=0;// *(uint64_t *)(&rcvKmers[j1*5 +2]);
//			maskedKmerH = maskedKmerH & 0x00FFFFFFFFFFFFFFul;
			uint64_t maskedKmer = getMmerPrefix(maskedKmerH, maskedKmerL, 10);
			int currentBin=0;
			while (maskedKmer >= localThreadKmerRanges[currentBin])
				currentBin++;
			currentBin--;
			if ((currentBin < 0) || (currentBin >= numThreads))
			{
				printf ("[%d] Tid %d numThreads %d j1 %lu curKmer %llx, maskedKmer %lu, localThreadKmerRanges[0]=%llx, currentBin %d kmerMask %llx, cond1 %d, cond2 %d\n",  rank, tid, numThreads, j1, *(uint64_t *)(&rcvKmers[j1*7]), maskedKmer, localThreadKmerRanges[0], currentBin, kmerMask, (currentBin < 0), (currentBin >= numThreads));
			}
			assert ((currentBin >= 0) && (currentBin < numThreads));
			uint16_t kmer[7];
			kmer[0] = rcvKmers[j1*7];
			kmer[1] = rcvKmers[j1*7 +1];
			kmer[2] = rcvKmers[j1*7 +2];
			kmer[3] = rcvKmers[j1*7 +3];
			kmer[4] = rcvKmers[j1*7 +4];
			kmer[5] = rcvKmers[j1*7 +5];
			kmer[6] = rcvKmers[j1*7 +6];
//			uint32_t readId = rcvKmers[j1*5 + 4];
			
			if (rangeWiseCounts[currentBin]==BUFSIZE)
			{
//			if (insertPos >= numKmers)
//			{
//				printf ("insertPos %lu currentBin %d\n", insertPos, currentBin);
//				for (i=0; i<numThreads; i++)
//				{
//					for (f=0; f<numThreads; f++)
//					{
//						printf ("%lu\t", threadRangeOfs[i*numThreads + f]);
//					}	
//					printf ("\n");
//				}
//			}
//			assert (insertPos < numKmers);
				uint64_t insertPos = localPtr[currentBin];
				for (i=0; i<BUFSIZE; i++)
				{
					sortedKmers[7*insertPos]=lKmers[currentBin][i].kmertuple[0];
					sortedKmers[7*insertPos + 1]=lKmers[currentBin][i].kmertuple[1];
					sortedKmers[7*insertPos + 2]=lKmers[currentBin][i].kmertuple[2];
					sortedKmers[7*insertPos + 3]=lKmers[currentBin][i].kmertuple[3];
					sortedKmers[7*insertPos + 4]=lKmers[currentBin][i].kmertuple[4];
					sortedKmers[7*insertPos + 5]=lKmers[currentBin][i].kmertuple[5];
					sortedKmers[7*insertPos + 6]=lKmers[currentBin][i].kmertuple[6];
					insertPos++;
				}
				localPtr[currentBin]+=BUFSIZE;
				lKmers[currentBin][0].kmertuple[0] = kmer[0];
				lKmers[currentBin][0].kmertuple[1] = kmer[1];
				lKmers[currentBin][0].kmertuple[2] = kmer[2];
				lKmers[currentBin][0].kmertuple[3] = kmer[3];
				lKmers[currentBin][0].kmertuple[4] = kmer[4];
				lKmers[currentBin][0].kmertuple[5] = kmer[5];
				lKmers[currentBin][0].kmertuple[6] = kmer[6];
//				lKmers[currentBin][0].range = currentBin;
				rangeWiseCounts[currentBin]=1;
			}
			else
			{
				int ofs=rangeWiseCounts[currentBin];
				lKmers[currentBin][ofs].kmertuple[0] = kmer[0];
				lKmers[currentBin][ofs].kmertuple[1] = kmer[1];
				lKmers[currentBin][ofs].kmertuple[2] = kmer[2];
				lKmers[currentBin][ofs].kmertuple[3] = kmer[3];
				lKmers[currentBin][ofs].kmertuple[4] = kmer[4];
				lKmers[currentBin][ofs].kmertuple[5] = kmer[5];
				lKmers[currentBin][ofs].kmertuple[6] = kmer[6];
//				lKmers[currentBin][ofs].range = currentBin;
				rangeWiseCounts[currentBin]++;
			}
		}
		for (i=0; i<numThreads; i++)
		{
			if (rangeWiseCounts[i] > 0)
			{
				uint64_t insertPos = localPtr[i];
				for (f=0; f<rangeWiseCounts[i]; f++)
				{
					sortedKmers[7*insertPos]=lKmers[i][f].kmertuple[0];
					sortedKmers[7*insertPos + 1]=lKmers[i][f].kmertuple[1];
					sortedKmers[7*insertPos + 2]=lKmers[i][f].kmertuple[2];
					sortedKmers[7*insertPos + 3]=lKmers[i][f].kmertuple[3];
					sortedKmers[7*insertPos + 4]=lKmers[i][f].kmertuple[4];
					sortedKmers[7*insertPos + 5]=lKmers[i][f].kmertuple[5];
					sortedKmers[7*insertPos + 6]=lKmers[i][f].kmertuple[6];
					insertPos++;
				}
				localPtr[i]=insertPos;
			}
		}
//		printf ("[%d] Tid %d After rangePartition localPtr[0]=%lu\n", rank, tid, localPtr[0]);
	} 
	gettimeofday (&eTime, NULL);
//		for (i=0; i<numThreads; i++)
//		{
//			for (f=0; f<numThreads; f++)
//			{
//				printf ("%lu,", threadRangeOfs[i*numThreads + f]);
//			}	
//			printf ("\n");
//		}
	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
    if (rank==0)
    	printf ("[%d] Range Partition Done (%lf sec)\n", rank, (double)elapsed/1000000);
	sortTime1 += elapsed;
//	if (rank == 0)
#ifdef DEBUG
	{
		for (i=0; i<numThreads; i++)
		{
			printf ("Verifying range partitions %d.... TotalKmers %lu p[%lu,%lu)\n", i, recvThreadOfs[numTasks*numThreads], localRangeOffsets[i], localRangeOffsets[i]+localRangeCounts[i]);
			for (i1=localRangeOffsets[i]; i1<(localRangeOffsets[i]+localRangeCounts[i]); i1++)
			{
				uint64_t curKmerL = *((uint64_t *)&sortedKmers[i1 * 7]);
				uint64_t curKmerH = 0;//*((uint64_t *)&sortedKmers[i1 * 5 +2]);
				uint64_t bucket=getMmerPrefix(curKmerH, curKmerL, 10);
				if ((bucket >= localThreadKmerRanges[i+1]) || (bucket < localThreadKmerRanges[i]))
				{
					printf ("[%d]i1 %lu curKmer %llx, %llx bucket %lu, expected[%llx,%llx)", rank, i1, curKmerH, curKmerL, bucket, localThreadKmerRanges[i], localThreadKmerRanges[i+1]);
					printf ("*\n");
				}
				assert ((bucket < localThreadKmerRanges[i+1]) && (bucket >= localThreadKmerRanges[i]));
			}
		}
	}
#endif
	int curThread=0;
//	for (curThread=0; curThread<numThreads; curThread++)
	gettimeofday (&sTime, NULL);
#if USE_OMP
	#pragma omp parallel num_threads(numThreads)
	{
		int curThread=omp_get_thread_num();
		assert (omp_get_num_threads() == numThreads);
#else
	{
		int curThread=0;
#endif
//		printf ("[%d] Tid %d SORT\n", rank, curThread);
//		cpu_bind ((rank*numThreads+curThread)%hwThreads);

//		qsort (&sortedKmers[rangeOffsets[curThread]*3], localRangeCounts[curThread], 3*sizeof(uint32_t), kmerCompare);
//		if (scanNo > 5)
//		{
//			localRangeCounts[curThread] = radixPass (0, sortedKmers, rcvKmers, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
//			localRangeCounts[curThread] = radixPass (1, rcvKmers, sortedKmers, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
//		}
//		localRangeCounts[curThread] = radixPass (2, sortedKmers, rcvKmers, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
//		localRangeCounts[curThread] = radixPass (3, rcvKmers, sortedKmers, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		localRangeCounts[curThread] = radixPass (5, sortedKmers, rcvKmers, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		if (curThread==0)
			printf ("5 Done\n");
		localRangeCounts[curThread] = radixPass (6, rcvKmers, sortedKmers, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		if (curThread==0)
			printf ("6 Done\n");
		localRangeCounts[curThread] = radixPass (7, sortedKmers, rcvKmers, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		if (curThread==0)
			printf ("7 Done\n");
		localRangeCounts[curThread] = radixPass (8, rcvKmers, sortedKmers, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		if (curThread==0)
			printf ("8 Done\n");
		localRangeCounts[curThread] = radixPass (9, sortedKmers, rcvKmers, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		if (curThread==0)
			printf ("9 Done\n");
		localRangeCounts[curThread] = radixPass (10, rcvKmers, sortedKmers, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		if (curThread==0)
			printf ("10 Done\n");
		localRangeCounts[curThread] = radixPass (11, sortedKmers, rcvKmers, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		if (curThread==0)
			printf ("11 Done\n");
		localRangeCounts[curThread] = radixPass (12, rcvKmers, sortedKmers, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		if (curThread==0)
			printf ("12 Done\n");
		localRangeCounts[curThread] = radixPass (13, sortedKmers, rcvKmers, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		if (curThread==0)
			printf ("13 Done\n");
//		break;
	}
	gettimeofday (&eTime, NULL);
	elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	memcpy (sortedKmers, rcvKmers, (rangeOffsets[numThreads-1]+localRangeCounts[numThreads-1])*7*sizeof(uint16_t));
    if (rank == 0)
    	printf ("[%d] Sort Done (%lf sec)\n", rank, (double)elapsed/1000000);
	sortTime2 += elapsed;
	uint64_t *lastAddr = (uint64_t *)(&(sortedKmers[(rangeOffsets[numThreads-1]+localRangeCounts[numThreads-1]-1)*7]));
//	if (rank == 0)
#ifdef DEBUG
	{
		for (i1=0; i1<numKmers; i1++)
		{
			uint64_t curKmer = *((uint64_t *)&sortedKmers[i1 * 3]);
			uint64_t bucket=(curKmer & kmerMask) >> rotater;
			if ((bucket >= processLimitUp) || (bucket < processLimitDn))
			{
				printf ("[%d]i1 %lu curKmer %llx bucket %lu", rank, i1, curKmer, bucket);
				printf ("*\n");
			}
			assert ((bucket < processLimitUp) && (bucket >= processLimitDn));
		}
	}
#endif
		printf ("[%d]LeastKmer %llx(%llx), "
			"MaxKmer %llx(%llx)\n", rank, ((uint64_t *)sortedKmers)[0], (((uint64_t *)sortedKmers)[0]&kmerMask) >> rotater, *lastAddr, (*lastAddr & kmerMask)>>rotater);

#ifdef DEBUG
	printf ("Verifying sort output...\n");
	int maxDegree=1;
	for (i=0; i<numThreads; i++)
	{
		uint64_t prevKmerL=*((uint64_t *)&sortedKmers[rangeOffsets[i]*7]);
		uint64_t prevKmerH=0;// *((uint64_t *)&sortedKmers[rangeOffsets[i]*5 +2]);
//		prevReadId = sortedKmers[rangeOffsets[i]*5 +4];
		int curDegree=1;
		for (i1=rangeOffsets[i]+1; i1<(rangeOffsets[i]+localRangeCounts[i]); i1++)
		{		
			uint64_t curKmerL = *((uint64_t *)&sortedKmers[i1 * 7]);
			uint64_t curKmerH = 0;// *(uint64_t *)(&(sortedKmers[i1*5 +2]));
//			uint32_t curReadId = sortedKmers[i1*5 +4];
/*			if (i1 < 1000) printf ("%lu\t%llx\n", i1, curKmerL);
			else*/
			{
				if ((curKmerH < prevKmerH) || ((curKmerH == prevKmerH) && (curKmerL < prevKmerL)))
					printf ("Sort Error cur %llx, prev %llx ofs %lu[%lu,%lu)\n", curKmerL, prevKmerL, i1, rangeOffsets[i], rangeOffsets[i]+localRangeCounts[i]);
				assert ((curKmerH > prevKmerH) || ((curKmerH == prevKmerH) && (curKmerL >= prevKmerL))); 
				if (curKmerL == prevKmerL) curDegree++;
				else
				{
					if (curDegree > maxDegree) maxDegree = curDegree;
					curDegree = 1;
				}
			}
			prevKmerL = curKmerL;
			prevKmerH = curKmerH;
//			prevReadId = curReadId;
		}
	}
	printf ("MaxDegree %d Done.\n", maxDegree);
#endif

//	for (i=0; i<numThreads; i++)
//		printf ("[%d]  threadKmerRanges[%d]=[%llx,%llx) Count %lu Ofs %lu\n", rank, i, localThreadKmerRanges[i], localThreadKmerRanges[i+1], localRangeCounts[i], localRangeOffsets[i]);
#endif
	printf ("Exit sort\n");
	free (localRangeOffsets);
	free (threadMask);
	free (threadRangeCounts);
	free (threadRangeOfs);
}


void markMischiefKmers (int scanNo, uint64_t *localRangeCounts, uint64_t *rangeOffsets )
{
	
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
//	for (i=0; i<numThreads; i++)
	#pragma omp parallel num_threads (numThreads) 
	{
		uint64_t i=0, j=0, k=0, l=0;
		i = omp_get_thread_num();
		uint64_t firstOfs = rangeOffsets[i];

		uint64_t curKmerL = *((uint64_t *)&(sortedKmers[firstOfs*7]));
		uint64_t curKmerH = 0;//*((uint64_t *)&(sortedKmers[2]));
		uint64_t curOfs=firstOfs;
		uint64_t numTotalKmers=0ul, numUniqueKmers = 0ul, numMischiefKmers=0;
		int beforeKmer[8]={0,0,0,0,0,0,0,0}, afterKmer[8]={0,0,0,0,0,0,0,0};
		int curKmerDegree=0;
//                printf ("Mark kmers %llu-%llu\n", rangeOffsets[i], rangeOffsets[i]+localRangeCounts[i]);
		for (j=rangeOffsets[i]; j<(rangeOffsets[i]+localRangeCounts[i]); j++)
		{
			uint64_t newKmerL = *((uint64_t *)&(sortedKmers[7*j]));
			uint64_t newKmerH = 0;//*((uint64_t *)&(sortedKmers[5*j + 2]));
			numTotalKmers ++;
			if (curKmerL != newKmerL)
			{
				numUniqueKmers++;
				int numInEdges=0, numOutEdges=0, maxInIndex = 0, maxOutIndex = 0, maxInCount = 0, maxOutCount = 0, sMaxInIndex=0, sMaxOutIndex = 0, sMaxInCount = 0, sMaxOutCount = 0;
				for (k=0; k<7; k++)
				{
					if (beforeKmer[k]) 
					{
						if (beforeKmer[k] > maxInCount)
						{
							sMaxInCount = maxInCount;
							sMaxInIndex = maxInIndex;
							maxInCount = beforeKmer[k];
							maxInIndex = k;
						}
						else if (beforeKmer[k] > sMaxInCount)
						{
							sMaxInCount = beforeKmer[k];
							sMaxInIndex = k;
						}
						numInEdges++;
					}
					if (afterKmer[k]) 
					{
						if (afterKmer[k] > maxOutCount)
						{
							sMaxOutCount = maxOutCount;
							sMaxOutIndex = maxOutIndex;
							maxOutCount = afterKmer[k];
							maxOutIndex = k;
						}
						else if (afterKmer[k] > sMaxOutCount)
						{
							sMaxOutCount = afterKmer[k];
							sMaxOutIndex = k;
						}
						numOutEdges++;
					}
//					beforeKmer[k] = afterKmer[k] = 0;
				}
				int disp=0;
/*				for (k=curOfs; k<j; k++)
				{
					if (sortedKmers[5*k+4] == 2616165)
					{
						disp=1;
						break;
					}
				}
				if (disp)
				{
					printf ("ReadId 2616166: kmer ");
					getKmerStr64 (*((uint64_t *)&sortedKmers[5*curOfs]), *((uint64_t *)&sortedKmers[5*curOfs+2]), NULL, 0);
					for (k=curOfs; k<j; k++)
					{
						printf ("\t%u\n", sortedKmers[5*k+4]);
					}
				}*/
				
				beforeKmer[7] = afterKmer[7] = 0;
				if (/*(curKmerDegree > 5) && */((numInEdges > 1) || (numOutEdges > 1)))
				{
					float maxInR = (float)maxInCount/curKmerDegree;
					float sMaxInR = (float)sMaxInCount/curKmerDegree;
					float maxOutR = (float)maxOutCount/curKmerDegree;
					float sMaxOutR = (float)sMaxOutCount/curKmerDegree;
					if ((sMaxInR > 0.4) && (sMaxOutR > 0.4))
					{
						numMischiefKmers++;
					}
					for (k=curOfs; k<j; k++)
					{
						uint8_t edge = (uint8_t)(sortedKmers[7*k+4] & 0xFF);
						assert (*((uint64_t *)&sortedKmers[7*k]) == curKmerL);
						uint8_t bIndex = (uint8_t)((sortedKmers[7*k + 4]) & 7);
						uint8_t aIndex = (uint8_t)((sortedKmers[7*k + 4] >> 3) & 0x7);
						if ((numInEdges > 1) && (numOutEdges > 1))
						{
							for (l=0; l<4; l++) sortedKmers[7*k + l] = 0;
							continue;
						}
/*						if ((sMaxInR > 0.4) && (sMaxOutR > 0.4))
						{
							printf ("Case0: Read %d numInEdges %d numOutEdges %d maxInCounts %d(%.2f),%d(%.2f) maxOutCounts %d(%.2f),%d(%.2f) kmerFreq %d Kmer %llx Prev %c Next %c - ", sortedKmers[5*k+4], numInEdges, numOutEdges, maxInCount, sMaxInCount, maxInR, sMaxInR, maxOutCount, sMaxOutCount, maxOutR, sMaxOutR, curKmerDegree, *((uint64_t *)&sortedKmers[5*k]), "ACTG    "[(edge>>3) & 7], "ACTG    "[edge&7]);
							for (l=0; l<7; l++)
								if (beforeKmer[l])
									printf ("%d, ",  beforeKmer[l]);
							printf ("\t");
							for (l=0; l<7; l++)
								if (afterKmer[l])
									printf ("%d, ",  afterKmer[l]);
	
							getKmerStr64 (*((uint64_t *)&sortedKmers[5*k]), *((uint64_t *)&sortedKmers[5*k+2]), NULL, 0);
						}*/
//						if (((sMaxInR > 0.2) || (sMaxOutR > 0.2))/* && (curKmerDegree > 20)*/)
						if (((numInEdges > 1) && (bIndex != maxInIndex)) || ((numOutEdges > 1) && (aIndex != maxOutIndex))) //filter out reads supporting lower kmer frequencies
//						if (((numInEdges > 1) && (beforeKmer[bIndex] == 1)) || ((numOutEdges > 1) && (afterKmer[aIndex] == 1)) /*|| ((numInEdges > 1) && (maxInCount < (curKmerDegree-10))) || ((numOutEdges > 1) && (maxOutCount < (curKmerDegree-10)))*/)
						{
							uint64_t creadid = (uint64_t)(*((uint32_t *)&sortedKmers[7*k + 5]));
							uint64_t crem = (uint64_t)(sortedKmers[7*k +4]>>8);
							creadid += crem*UINT32_MAX;
//							assert (sortedKmers[5*k+4] != 0xBADBADBADl);
							uint64_t badReadId = 0xBADBADBADl;
							uint32_t badReadIdRem = (uint32_t)(badReadId % UINT32_MAX);
							uint8_t badReadIdQuot = (uint8_t) (badReadId / UINT32_MAX);
							memcpy (&parent[creadid*5 +1], &badReadIdRem, sizeof(uint32_t));
							parent[creadid*5] = badReadIdQuot;
							filteredReads[creadid] = 1;

/*							if (i==0)
							{
								printf ("Case1: Read %u numInEdges %d numOutEdges %d maxInCounts %d(%.2f),%d(%.2f) maxOutCounts %d(%.2f),%d(%.2f) kmerFreq %d Kmer %llx Prev %c Next %c - ", *((uint32_t *)&sortedKmers[7*k+5]), numInEdges, numOutEdges, maxInCount, sMaxInCount, maxInR, sMaxInR, maxOutCount, sMaxOutCount, maxOutR, sMaxOutR, curKmerDegree, *((uint64_t *)&sortedKmers[7*k]), "ACTG    "[(edge>>3) & 7], "ACTG    "[edge&7]);
								for (l=0; l<7; l++)
									if (beforeKmer[l])
										printf ("%d, ",  beforeKmer[l]);
								printf ("\t");
								for (l=0; l<7; l++)
									if (afterKmer[l])
										printf ("%d, ",  afterKmer[l]);
	
								getKmerStr64 (*((uint64_t *)&sortedKmers[7*k]), *((uint64_t *)&sortedKmers[7*k+4]), NULL, 0);
							}*/
						}
						else
						{
							if ((sMaxInR > 0.2) || (sMaxOutR > 0.2)) //filter out reads supporting max frequency edge only if seconnd max frequency is << max frequency edge
							{
								uint64_t creadid = (uint64_t)(*((uint32_t *)&sortedKmers[7*k + 5]));
								uint64_t crem = (uint64_t)(sortedKmers[7*k +4]>>8);
								creadid += crem*UINT32_MAX;
//								assert (sortedKmers[5*k+4] != 0xBADBADBADl);
								uint64_t badReadId = 0xBADBADBADl;
								uint32_t badReadIdRem = (uint32_t)(badReadId % UINT32_MAX);
								uint8_t badReadIdQuot = (uint8_t) (badReadId / UINT32_MAX);
								memcpy (&parent[creadid*5 +1], &badReadIdRem, sizeof(uint32_t));
								parent[creadid*5] = badReadIdQuot;
								filteredReads[creadid] = 1;
//								parent[creadid] = 0xBADBADBADl;
/*								if (i==0)//tid
								{
									printf ("Case2: Read %u numInEdges %d numOutEdges %d maxInCounts %d(%.2f),%d(%.2f) maxOutCounts %d(%.2f),%d(%.2f) kmerFreq %d Kmer %llx Prev %c Next %c - ", *((uint32_t *)&sortedKmers[7*k+5]), numInEdges, numOutEdges, maxInCount, sMaxInCount, maxInR, sMaxInR, maxOutCount, sMaxOutCount, maxOutR, sMaxOutR, curKmerDegree, *((uint64_t *)&sortedKmers[7*k]), "ACTG    "[(edge>>3) & 7], "ACTG    "[edge&7]);
									for (l=0; l<7; l++)
										if (beforeKmer[l])
											printf ("%d, ",  beforeKmer[l]);
									printf ("\t");
									for (l=0; l<7; l++)
										if (afterKmer[l])
											printf ("%d, ",  afterKmer[l]);
			
									getKmerStr64 (*((uint64_t *)&sortedKmers[7*k]), *((uint64_t *)&sortedKmers[7*k+4]), NULL, 0);
								}*/
							}
						}
//						if (((numInEdges > 1) && (maxInCount < (curKmerDegree-10))) || ((numOutEdges > 1) && (maxOutCount < (curKmerDegree-10))))
//							parent[sortedKmers[5*k+4]] = 0xBADBADB;
//						printf ("parent[%u]=0xBADBADB\n", sortedKmers[5*k+4]);
//						if (k>curOfs) sortedKmers[5*k+4] = sortedKmers[5*curOfs+4];
					}
				}
				curKmerL = newKmerL;
				curOfs = j;
				curKmerDegree=0;
				for (l=0; l<7; l++)
					beforeKmer[l] = afterKmer[l] = 0;
			}
//			else
			{
				uint8_t bIndex = (uint8_t)((sortedKmers[7*j + 4]) & 7);
				uint8_t aIndex = (uint8_t)((sortedKmers[7*j + 4] >> 3) & 0x7);
				beforeKmer[bIndex]++; 
				afterKmer[aIndex]++;
				curKmerDegree++;
			}
		}
		printf ("[%d,%d]numTotalKmers %lu NumUniqueKmers %lu NumMischiefKmers %lu\n", rank, i, numTotalKmers, numUniqueKmers, numMischiefKmers);
	}
	gettimeofday (&eTime, NULL);
        long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	if (rank == 0)
        	printf ("[%d] Filter Done (%lf sec)\n", rank, (double)elapsed/1000000);
        filterTime = elapsed;
	uint64_t numBadReadsBefore = 0, numBadReadsBefore2 = 0, numBadReadsBefore3 = 0, i=0, numBadReadsAfter= 0;
	for (i=0; i<numReads; i++)
	{
		uint64_t compId = (uint64_t)(*((uint32_t *)&parent[i*5 +1])) + (uint64_t)parent[i*5]*UINT32_MAX;
		if (compId==0xBADBADBADl) 
		{
			numBadReadsBefore++;
		}
		if (filteredReads[i]) numBadReadsBefore2++;
	}
#if USE_MPI
//	assert (numReads < 1000000000l);
    uint64_t curReadStart=0;
    int numBytesPerTransfer = 1000000000;
    while (curReadStart < numReads)
    {
        uint64_t bytesToSend = numReads-curReadStart;
        if (bytesToSend > numBytesPerTransfer)
        {
            bytesToSend = numBytesPerTransfer;
        }
        if (rank == 0)
            printf ("markM allreduce %d reads from offset %lu\n", (int)bytesToSend, curReadStart);
	    MPI_Allreduce (MPI_IN_PLACE, &filteredReads[curReadStart], (int)bytesToSend, MPI_BYTE, MPI_BOR, MPI_COMM_WORLD);
        curReadStart += bytesToSend;
    }

#endif
	for (i=0; i<numReads; i++)
	{
		uint64_t compId = (uint64_t)(*((uint32_t *)&parent[i*5 +1])) + (uint64_t)parent[i*5]*UINT32_MAX;
		if (compId==0xBADBADBADl) 
			assert (filteredReads[i]);
		else if (filteredReads[i]) 
		{
			numBadReadsBefore3++;
			uint64_t badReadId = 0xBADBADBADl;
			uint32_t badReadIdRem = (uint32_t)(badReadId % UINT32_MAX);
			uint8_t badReadIdQuot = (uint8_t) (badReadId / UINT32_MAX);
			memcpy (&parent[i*5 +1], &badReadIdRem, sizeof(uint32_t));
			parent[i*5] = badReadIdQuot;
		}
	}
	printf ("[%d] Local numBadReadsBefore %lu (%lu) Extra %lu\n", rank, numBadReadsBefore, numBadReadsBefore2, numBadReadsBefore3);
}
void unionFindPar2 (int scanNo, uint64_t *localRangeCounts, uint64_t *rangeOffsets)
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
//	uint64_t i=0;
/*#if RADIX_OPT
	FILE *fpEdge=fopen("edges.opt", "w");
#else
	FILE *fpEdge=fopen("edges.noopt", "w"); 
#endif*/
//	int tid;
//	for (tid=0; tid<numThreads; tid++)
	#pragma omp parallel num_threads (numThreads)
	{
		uint64_t i=0;
		int numUnions=0;
		uint32_t maxDegree = 0;
        	int numRandAccess=0, maxRandAccess=0;
		int tid = omp_get_thread_num();
		uint64_t curKmerL = *(uint64_t *)&(sortedKmers[(numKmers-1)*7]);
		uint64_t curKmerH = 0;//*(uint64_t *)&(sortedKmers[(numKmers-1)*5 +2]);
		uint64_t minReadId = (uint64_t)(*((uint32_t *)&sortedKmers[5]));
		uint64_t minReadIdExtra = (uint64_t)(sortedKmers[4]>>8);
		minReadId = minReadIdExtra*UINT32_MAX + minReadId;

		uint64_t prevReadId = minReadId;
		uint32_t curVertexDegree=0;
		int processThisKmer=1;
		int disp=0;
		uint64_t startReadId = numReads*tid/numThreads;
		uint64_t endReadId = numReads*(tid+1)/numThreads;
        if (tid == (numThreads-1))
            endReadId = numReads;
		uint64_t insPos = rangeOffsets[tid];
        struct timeval sTime1, eTime1;
    	gettimeofday (&sTime1, NULL);
		for (i=rangeOffsets[tid]; i<(rangeOffsets[tid]+localRangeCounts[tid]); i++)
		{
//			if ((i % 10000000)==0)
//				printf 	("[%d] %lu/%lu\n", rank, i, numKmers);
			uint64_t newKmerL = *((uint64_t *)&(sortedKmers[7*i]));
			uint64_t newKmerH = 0;//*((uint64_t *)&(sortedKmers[5*i +2]));
//			assert (newKmerH == 0);
			uint64_t curReadId = (uint64_t)(*((uint32_t *)&sortedKmers[7*i + 5]));
			uint64_t curReadIdExtra = (uint64_t)(sortedKmers[7*i + 4]>>8);
			curReadId = curReadIdExtra*UINT32_MAX + curReadId;

			uint64_t curCompId = (uint64_t)(*((uint32_t *)&parent[curReadId*5 +1])) + (uint64_t)parent[curReadId*5]*UINT32_MAX;
/*			if ((curReadId==2892405) || (curReadId == 631217))
			{
				printf ("ReadId %lu cop %lu\n", curReadId, curCompId);
			}*/

			if (curCompId == 0xBADBADBADl) continue;
			if (newKmerL == 0) continue;
//			if ((curReadId < startReadId) || (curReadId >= endReadId)) continue;
			if ((newKmerH != curKmerH) || (newKmerL != curKmerL))
			{
				disp=0;
				curKmerH = newKmerH;
				curKmerL = newKmerL;
				minReadId = curReadId;
				if (minReadId==1663384)
				{
					disp=1;
//					printf ("tid %d kmer %llx, read %u\n", curKmer, minReadId);
				}
/*				if (kmerDegree != NULL)
				{*/
					prevReadId = curReadId;
					curVertexDegree=0;
					processThisKmer = 1;
					uint64_t i1=i;
					while ((i1 < numKmers) && (*((uint64_t *)&(sortedKmers[7*i1])) == newKmerL)/* && (*((uint64_t *)&(sortedKmers[5*i1 +2])) == newKmerH)*/)
					{
						uint64_t tmpReadId = (uint64_t)(*((uint32_t *)&sortedKmers[7*i1 + 5]));
						uint64_t tmpReadIdExtra = (uint64_t)(sortedKmers[7*i1 + 4]>>8);
						tmpReadId += tmpReadIdExtra*UINT32_MAX;
						if ((i1 == i) || (tmpReadId != prevReadId))
						{
							curVertexDegree++;
						}
						prevReadId = tmpReadId;
						i1++;
					}
/*					kmerDegree [curVertexDegree]++;*/
/*					if ((curVertexDegree >= 30) || (curVertexDegree < 10))
					{
						processThisKmer=0;
					}*/
					if ((processThisKmer) && (curVertexDegree > maxDegree))
					{
//						printf ("maxDegree %u\n", curVertexDegree);
						maxDegree = curVertexDegree;
					}
/*				}*/
			}	
			else
			{
				if (processThisKmer == 1) 
				{
					uint64_t p_u = minReadId;//findRoot (parent, minReadId, &numRandAccess);
					uint64_t p_v = curReadId;//findRoot (parent, curReadId, &numRandAccess);
                                        if (p_u != p_v)
                                        {
                                            int rAccess=0;
                    //                        if (tid == 0)
                    //                        {
                    //                            printf ("RandAcc%d u %u v %u\n", tid, p_u, p_v);
                    //                        }
                    				assert ((p_u != 0xBADBADBADl) && (p_v != 0xBADBADBADl));
						p_u = findRoot (parent, p_u, &numRandAccess);
						p_v = findRoot (parent, p_v, &numRandAccess);
                    				assert ((p_u != 0xBADBADBADl) && (p_v != 0xBADBADBADl));
/*                                            while (parent[p_u] != p_u)
                                            {
                                                if (parent[p_u] != parent[parent[p_u]])
                                                    parent[p_u] = parent[parent[p_u]];
                                                p_u = parent[p_u];
                                                rAccess = rAccess + 1;
                                            }
                                            numRandAccess = numRandAccess + rAccess;
                                            if (maxRandAccess < rAccess)
                                                maxRandAccess = rAccess;
                                            rAccess=0;
                    
                                            while (parent[p_v] != p_v)
                                            {
                                                if (parent[p_v] != parent[parent[p_v]])
                                                    parent[p_v] = parent[parent[p_v]];
                                                p_v = parent[p_v];
                                                rAccess = rAccess + 1;
                                            }
                                            numRandAccess = numRandAccess + rAccess;
                                            if (maxRandAccess < rAccess)
                                                maxRandAccess = rAccess;*/
                                        }

/*					if ((disp) && ((curReadId==14900807) || (curReadId==14987420) || (curReadId==12950003)))
					{
						printf ("tid %d  ofs %lu kmer %llx u %u v %u p_u %u p_v %u\n", tid, i, curKmer, minReadId, curReadId, p_u, p_v);
					}*/
//					uint64_t p_5908271 = findRoot (parent, 5908271, &numRandAccess);
//					uint64_t p_8329065 = findRoot (parent, 8329065, &numRandAccess);
					if (p_u != p_v)
					{
/*						if (scanNo==1)
							fprintf (fpEdge, "%u,%u\n", minReadId, curReadId);*/
						uint32_t mreadid = (uint32_t)(minReadId % UINT32_MAX);
						uint16_t mrem = (uint16_t)(minReadId / UINT32_MAX);
						memcpy (&sortedKmers[7*insPos + 6*numUnions], &mreadid, sizeof(uint32_t));
						sortedKmers[7*insPos + 6*numUnions + 2] = mrem;
						uint32_t creadid = (uint32_t)(curReadId % UINT32_MAX);
						uint16_t crem = (uint16_t)(curReadId / UINT32_MAX);
						memcpy (&sortedKmers[7*insPos + 6*numUnions + 3], &creadid, sizeof(uint32_t));
						sortedKmers[7*insPos + 6*numUnions + 5] = crem;
						numUnions++;
/*						if ((minReadId == 383552) || (curReadId == 383552))
						{
							printf ("%u -> %u Kmer %llx-%llx ", minReadId, curReadId, curKmerH, curKmerL);
							getKmerStr64 (curKmerL, curKmerH, NULL, 0);
						}*/
		
						if (p_u > p_v)
						{
//							printf ("Rank %d union: %lu (%lu)-->%lu (%lu)\n", rank, minReadId, p_v, curReadId, p_u);
//							if ((p_v == p_5908271) || (p_v == p_8329065))
//							{
//								printf ("%u (P %u) <-- %u (%u) p_5908271 %u p_8329065 %u kmer ", minReadId, p_u, curReadId, p_v, p_5908271, p_8329065);
//								getKmerStr64 (curKmerL, curKmerH, NULL, 0);
//							}
							uint32_t tmpRem = (uint32_t)(p_u % UINT32_MAX);
							uint8_t tmpQuot = (uint8_t) p_u/UINT32_MAX;
							memcpy (&parent[p_v*5 + 1], &tmpRem, sizeof(uint32_t));
							parent[p_v*5] = tmpQuot;
//				                        parent[p_v]=p_u;
							if (p_u == 0)
								printf ("%lu-->%lu\n", p_v, p_u);
//							unionOp (parent, compSizes, p_u, p_v);
						}
						else
						{
//							if ((p_u == p_5908271) || (p_u == p_8329065))
//							{
//								printf ("%u (P %u) --> %u (%u) p_5908271 %u p_8329065 %u kmer ", minReadId, p_u, curReadId, p_v, p_5908271, p_8329065);
//								getKmerStr64 (curKmerL, curKmerH, NULL, 0);
//							}
//							printf ("Rank %d union: %lu (%lu)-->%lu (%lu)\n", rank, minReadId, p_u, curReadId, p_v);
							uint32_t tmpRem = (uint32_t)(p_v % UINT32_MAX);
							uint8_t tmpQuot = (uint8_t) p_v/UINT32_MAX;
							memcpy (&parent[p_u*5 + 1], &tmpRem, sizeof(uint32_t));
							parent[p_u*5] = tmpQuot;
//                        			    	parent[p_u]=p_v;
							if (p_v == 0)
								printf ("%lu-->%lu\n", p_u, p_v);
//							unionOp (parent, compSizes, p_v, p_u);
						}
					}
				}
			}
			prevReadId = curReadId;
		}
    	gettimeofday (&eTime1, NULL);
	    long elapsed1 = (eTime1.tv_sec * 1000000 + eTime1.tv_usec) - (sTime1.tv_sec * 1000000 + sTime1.tv_usec);
		printf ("[%d]UnionFindKL Tid %d, Scan %d NumUnions %d Prcessed %lu numRandAccesses %d maxRandAccess %d Time %lf sec\n", rank, tid, scanNo, numUnions, localRangeCounts[tid], numRandAccess, maxRandAccess, (double)elapsed1/1000000);
		localRangeCounts[tid]=numUnions;
	}
//	fclose(fpEdge);
//	fclose (fpOut);
//	free (sortedKmers);
//	sortedKmers=NULL;
	gettimeofday (&eTime, NULL);
	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
//	printf ("[%d]UnionFindKL2- %lf sec\n", rank, (double) elapsed/1000000);
	unionFindTime1 += elapsed;
}
void unionFind (int scanNo, uint64_t *localRangeCounts, uint64_t *rangeOffsets)
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
	uint64_t curKmer = *(uint64_t *)&(sortedKmers[(numKmers-1)*3]);
	uint32_t curVertexDegree=0;
	int processThisKmer=1;//, tid=0;
	int disp=0;
/*#if RADIX_OPT
	FILE *fpEdge=fopen("edges.opt", "w");
#else
	FILE *fpEdge=fopen("edges.noopt", "w"); 
#endif*/
//	for (tid=0; tid<numThreads; tid++)
	#pragma omp parallel num_threads(numThreads)
	{
		int tid=omp_get_thread_num();
        int numRandAccess = 0;
//		printf ("Tid %d start %lu, end %lu\n", tid, rangeOffsets[tid], rangeOffsets[tid]+localRangeCounts[tid]);
		while (localRangeCounts[tid] > 0)
		{
			int numUnions=0;
			uint64_t i=0;
			for (i=0; i<localRangeCounts[tid]; i++)
			{
	//			if ((i % 10000000)==0)
	//				printf 	("[%d] %lu/%lu\n", rank, i, numKmers);
				uint64_t fromReadId = (uint64_t)(*((uint32_t *)&sortedKmers[7*rangeOffsets[tid] + 6*i]));
				uint64_t fromReadIdExtra = (uint64_t)sortedKmers[7*rangeOffsets[tid] + 6*i + 2];
				fromReadId += fromReadIdExtra*UINT32_MAX;
				uint64_t toReadId = (uint64_t)(*((uint32_t *)&sortedKmers[7*rangeOffsets[tid] + 6*i +3]));
				uint64_t toReadIdExtra = (uint64_t)sortedKmers[7*rangeOffsets[tid] + 6*i +5];
				toReadId += toReadIdExtra*UINT32_MAX;

				uint64_t parent_from = (uint64_t)(*((uint32_t *)&parent[fromReadId*5 +1])) + (uint64_t)parent[fromReadId*5]*UINT32_MAX;
				uint64_t parent_to = (uint64_t)(*((uint32_t *)&parent[toReadId*5 +1])) + (uint64_t)parent[toReadId*5]*UINT32_MAX;
				assert ((parent_from != 0xBADBADBADl) && (parent_to != 0xBADBADBADl));
				assert ((fromReadId != 0xBADBADBADl) && (toReadId != 0xBADBADBADl));
				uint64_t p_u = findRoot (parent, fromReadId, &numRandAccess);
				uint64_t p_v = findRoot (parent, toReadId, &numRandAccess);
				if (p_u != p_v)
				{
					numUnions++;
//					sortedKmers[5*rangeOffsets[tid] + 2*numUnions] = fromReadId;
//					sortedKmers[5*rangeOffsets[tid] + 2*numUnions + 1] = toReadId;
					memcpy (&sortedKmers[7*rangeOffsets[tid] + 6*numUnions], &sortedKmers[7*rangeOffsets[tid] + 6*i], 6*sizeof(uint16_t));
					if (p_u > p_v)
						unionOp (parent, compSizes, p_u, p_v);
					else
						unionOp (parent, compSizes, p_v, p_u);
				}
			}
			localRangeCounts[tid] = numUnions;
			printf ("[%d]UnionFind Tid %d, Scan %d NumUnions %d numRandAccesses %d\n", rank, tid, scanNo, numUnions, numRandAccess);
		}
	}
//	fclose(fpEdge);
//	fclose (fpOut);
//	free (sortedKmers);
//	sortedKmers=NULL;
	gettimeofday (&eTime, NULL);
	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
//	printf ("[%d]UnionFind- %lf sec\n", rank, (double) elapsed/1000000);
	unionFindTime1 += elapsed;
}
void unionFindPeReads (int upto)
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
	uint64_t i=0;
	int numRandAccess=0, numUnions=0;
//	for (i=0; i<numReads; i++)
//		if (parent[i]==0xBADBADB) parent[i] = i;

	for (i=0; i<numReads; i+=2)
	{
		uint64_t fromReadId = i;
		uint64_t toReadId = i+1;
		uint64_t parent_from = (uint64_t)(*((uint32_t *)&parent[fromReadId*5 +1])) + (uint64_t)parent[fromReadId*5]*UINT32_MAX;
		uint64_t parent_to = (uint64_t)(*((uint32_t *)&parent[toReadId*5 +1])) + (uint64_t)parent[toReadId*5]*UINT32_MAX;
		if ((parent_from == 0xBADBADBADl) || (parent_to == 0xBADBADBADl)) 
		{
			uint32_t from_rem = (uint32_t)(fromReadId % UINT32_MAX);
			uint8_t from_quot = (uint8_t) (fromReadId / UINT32_MAX);

			uint32_t to_rem = (uint32_t) (toReadId % UINT32_MAX);
			uint8_t to_quot = (uint8_t) (toReadId / UINT32_MAX);
			if (parent_from == 0xBADBADBADl)
			{
				memcpy (&parent[fromReadId*5 + 1], &from_rem, sizeof(uint32_t));
				parent[fromReadId*5] = from_quot;
			}
			if (parent_to == 0xBADBADBADl) 
			{
				memcpy (&parent[toReadId*5 + 1], &to_rem, sizeof(uint32_t));
				parent[toReadId*5] = to_quot;
			}

		}
			parent_from = findRoot (parent, fromReadId, &numRandAccess);
			parent_to = findRoot (parent, toReadId, &numRandAccess);
//		parent_from = (uint64_t)(*((uint32_t *)&parent[fromReadId*5 +1])) + (uint64_t)parent[fromReadId*5]*UINT32_MAX;
//		parent_to = (uint64_t)(*((uint32_t *)&parent[toReadId*5 +1])) + (uint64_t)parent[toReadId*5]*UINT32_MAX;
		assert ((parent_from != 0xBADBADBADl) && (parent_to != 0xBADBADBADl));

//		if ((parent_from == fromReadId) || (parent_to == toReadId))
		if ((compSizes[parent_from] <= upto) || (compSizes[parent_to] <= upto))
		{
			assert ((fromReadId != 0xBADBADBADl) && (toReadId != 0xBADBADBADl));
			uint64_t p_u = findRoot (parent, fromReadId, &numRandAccess);
			uint64_t p_v = findRoot (parent, toReadId, &numRandAccess);
			if (p_u != p_v)
			{
				numUnions++;
				if (p_u > p_v)
				{
					unionOp (parent, compSizes, p_u, p_v);
	//				printf ("Exp %u-->%u But Act %u\n", p_v, p_u, parent[p_v]);
				}
				else
				{
					unionOp (parent, compSizes, p_v, p_u);
	//				printf ("Exp %u-->%u But Act %u\n", p_u, p_v, parent[p_u]);
				}
			}
		}
	}
	gettimeofday (&eTime, NULL);
	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	peUnionTime = elapsed;
	printf ("[%d]UnionFindPeReads- %lf sec NumUnions %d\n", rank, (double) elapsed/1000000, numUnions);
}
/*void unionFindSmallComps (uint32_t upto)
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
	uint32_t i=0;
	int numRandAccess=0, numUnions=0;
//	for (i=0; i<numReads; i++)
//		if (parent[i]==0xBADBADB) parent[i] = i;

	for (i=0; i<numReads/2; i++)
	{
		uint32_t fromReadId = i;
		uint32_t toReadId = numReads/2 + i;
//		if ((compSizes[fromReadId] < upto) || (compSizes[toReadId] < upto))
		{
			uint32_t p_u = findRoot (parent, fromReadId, &numRandAccess);
			uint32_t p_v = findRoot (parent, toReadId, &numRandAccess);
			assert ((p_u != 0xBADBADBADl) && (p_v != 0xBADBADBADl));
			if (p_u != p_v)
			{
				if (p_u > p_v)
				{
					sortedKmers[3*numUnions] = p_v;
					sortedKmers[3*numUnions + 1] = p_u;
	//				unionOp (parent, compSizes, p_u, p_v);
	//				printf ("Exp %u-->%u But Act %u\n", p_v, p_u, parent[p_v]);
				}
				else
				{
					sortedKmers[3*numUnions] = p_u;
					sortedKmers[3*numUnions + 1] = p_v;
	//				unionOp (parent, compSizes, p_v, p_u);
	//				printf ("Exp %u-->%u But Act %u\n", p_u, p_v, parent[p_u]);
				}
				numUnions++;
			}
		}
	}
	qsort (sortedKmers, numUnions, 3*sizeof(uint32_t), edgeComp);
	uint32_t fromComp=sortedKmers[0], toComp=sortedKmers[1], count=0, insIndex=0;
	
	for (i=0; i<numUnions; i++)
	{
//			unionOp (parent, compSizes, sortedKmers[3*i+1], sortedKmers[3*i]);
		if ((sortedKmers[3*i] != fromComp) || (sortedKmers[3*i +1] != toComp))
		{
			if (count > 0)
			{
				assert ((3*insIndex+2) < (3*i));
				{
					sortedKmers[3*insIndex] = fromComp;
					sortedKmers[3*insIndex + 1] = toComp;
					sortedKmers[3*insIndex + 2] = count;
					insIndex++;
				}
//				printf ("Edge: %u, %u - %u\n", fromComp, toComp, count);
			}
			fromComp = sortedKmers[3*i];
			toComp = sortedKmers[3*i +1];
			count=0;
		}
		count++;
	}
	if (count > 0)
	{
		sortedKmers[3*insIndex] = fromComp;
		sortedKmers[3*insIndex + 1] = toComp;
		sortedKmers[3*insIndex + 2] = count;
		insIndex++;
//		printf ("Edge: %u, %u - %u\n", fromComp, toComp, count);
	}
	qsort (sortedKmers, insIndex, 3*sizeof(uint32_t), edgeCountComp);
	numUnions=0;
	for (i=0; i<insIndex; i++)
	{
//		printf ("Edge2: %u, %u - %u\n", sortedKmers[3*i], sortedKmers[3*i+1], sortedKmers[3*i+2]);
		if (sortedKmers[3*i+2] > upto)
		{
			uint64_t p_u = findRoot (parent, sortedKmers[3*i], &numRandAccess);
			uint64_t p_v = findRoot (parent, sortedKmers[3*i+1], &numRandAccess);
			if (p_u != p_v)
			{
				if (p_u > p_v)
					unionOp (parent, compSizes, p_u, p_v);
				else
					unionOp (parent, compSizes, p_v, p_u);
				numUnions++;
			}
		}
	}
	gettimeofday (&eTime, NULL);
	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	printf ("[%d]UnionFindSmallComps(%u)- %lf sec NumUnions %u, insIndex %u\n", rank, upto, (double) elapsed/1000000, numUnions, insIndex);
}*/
/*void mergeStep1 (uint64_t *otherRankParent)
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
	uint32_t i=0;
	uint32_t numUnions=0;
    int numRandAccesses=0;
	for (i=0; i<numReads; i++)
	{
        if (parent[i]==otherRankParent[i])
            continue;
		uint32_t p_u = findRoot (parent, i, &numRandAccesses);
		uint32_t p_v = findRoot (parent, otherRankParent[i], &numRandAccesses);
		if (p_u != p_v)
		{
            if (p_u > p_v)
    			unionOp (parent, compSizes, p_u, p_v);
            else
                unionOp (parent, compSizes, p_v, p_u);
			numUnions++;
		}
	}
	gettimeofday (&eTime, NULL);
	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
//	printf ("[%d]UnionFindMerge %lf sec, NumUnions %d\n", rank, (double) elapsed/1000000, numUnions);
//	unionFindTime1 += elapsed;
}*/
uint64_t mergeStepOpt (uint8_t *otherRankParent)
{
	uint64_t numEdges=0;
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
	uint8_t *globalEdgeBuffer = (uint8_t *)compSizes;
	uint64_t globalOffset = 0;
        #pragma omp parallel num_threads(numThreads)
        {
            uint32_t numUnions=1;
            int numRandAccesses = 0;
            int tid = omp_get_thread_num();
            uint64_t startOfs = (uint64_t)tid * numReads / numThreads;
            uint64_t endOfs = (uint64_t)(tid+1) * numReads / numThreads;
            uint64_t i=0;
#define EDGE_BUF_SIZE 20000
	    uint8_t tmpEdge[EDGE_BUF_SIZE*5*2];
	    int tmpCount=0;
    
            while (numUnions > 0)
            {
                numUnions = 0;
       	        for (i=startOfs; i<endOfs; i++)
                {
                    uint64_t thisParent = (uint64_t)(*((uint32_t *)&parent[i*5 + 1])) + (uint64_t)parent[i*5]*UINT32_MAX; 
		    if  ((thisParent==i) || (thisParent == 0xBADBADBADl))
			continue;
                    uint64_t otherParent1 = (uint64_t)(*((uint32_t *)&otherRankParent[i*5 + 1])) + (uint64_t)otherRankParent[i*5]*UINT32_MAX; ///from Rank0
                    if (otherParent1 == 0xBADBADBADl)
                    {
			uint32_t numread_rem = (uint32_t)(i % UINT32_MAX);
			uint8_t numread_quot = (uint8_t)(i / UINT32_MAX);
                        memcpy (&otherRankParent[i*5 +1], &numread_rem, sizeof(uint32_t));
			otherRankParent[i*5] = numread_quot;
                        otherParent1 = i;
                    }
                    uint64_t otherParent2 = (uint64_t)(*((uint32_t *)&otherRankParent[thisParent*5 + 1])) + (uint64_t)otherRankParent[thisParent*5]*UINT32_MAX; ///from Rank0
                    if (otherParent2 == 0xBADBADBADl)
                    {
			uint32_t numread_rem = (uint32_t)(thisParent % UINT32_MAX);
			uint8_t numread_quot = (uint8_t)(thisParent / UINT32_MAX);
                        memcpy (&otherRankParent[thisParent*5 +1], &numread_rem, sizeof(uint32_t));
			otherRankParent[thisParent*5] = numread_quot;
                        otherParent2 = thisParent;
                    }
//		    assert ((i != 0xBADBADBADl) && (otherParent != 0xBADBADBADl));
                    uint64_t p_u = findRoot (otherRankParent, otherParent1, &numRandAccesses);
                    uint64_t p_v = findRoot (otherRankParent, otherParent2, &numRandAccesses);
//		    assert ((p_u != 0xBADBADBADl) && (p_v != 0xBADBADBADl)) printf ("merge - %llx %llx\n", p_u, p_v);
                    if (p_u != p_v)
                    {
			if (tmpCount == EDGE_BUF_SIZE)
			{
				uint64_t insertOfs=0;
				#pragma omp atomic capture
				{
					insertOfs = globalOffset;
					globalOffset += EDGE_BUF_SIZE;
				}
				assert ((insertOfs+EDGE_BUF_SIZE) <= (numReads/2));
				memcpy (&globalEdgeBuffer[insertOfs*5*2], tmpEdge, EDGE_BUF_SIZE*5*2);
				tmpCount=0;
			}
			uint32_t tmpRem = (uint32_t)(p_u % UINT32_MAX);
			uint8_t tmpQuot = (uint8_t) p_u/UINT32_MAX;
			memcpy (&tmpEdge[tmpCount*10+1], &tmpRem, sizeof(uint32_t));
			tmpEdge[tmpCount*10] = tmpQuot;

			tmpRem = (uint32_t)(p_v % UINT32_MAX);
			tmpQuot = (uint8_t) p_v/UINT32_MAX;
			memcpy (&tmpEdge[tmpCount*10+5+1], &tmpRem, sizeof(uint32_t));
			tmpEdge[tmpCount*10+5] = tmpQuot;

			tmpCount++;

                        if (p_u > p_v)
                		unionOp (otherRankParent, NULL, p_u, p_v);
                        else
                            unionOp (otherRankParent, NULL, p_v, p_u);
                        numUnions++;
                    }
                    else
                    {
/*			uint32_t numread_rem = (uint32_t)(numReads % UINT32_MAX);
			uint8_t numread_quot = (uint8_t)(numReads / UINT32_MAX);
                        memcpy (&otherRankParent[i*5 +1], &numread_rem, sizeof(uint32_t));
			otherRankParent[i*5] = numread_quot;*/
//                        otherRankParent[i] = numReads;
                    }
                }
		if (tmpCount > 0)
		{
			uint64_t insertOfs=0;
			#pragma omp atomic capture
			{
				insertOfs = globalOffset;
				globalOffset += tmpCount;
			}
			assert ((insertOfs+tmpCount) <= (numReads/2));
			memcpy (&globalEdgeBuffer[insertOfs*5*2], tmpEdge, tmpCount*5*2);
			tmpCount=0;
		}
		printf ("%d\t%d\tmergeStepOpt: numUnions %d, globalOffset %lu\n", rank, tid, numUnions, globalOffset);
            }
        }
	gettimeofday (&eTime, NULL);
	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	printf ("[%d]UnionFindMergeOpt %lf sec globalOffset %lu\n", rank, (double) elapsed/1000000, globalOffset);
	return globalOffset;
//	unionFindTime1 += elapsed;
}

void mergeStep (uint8_t *otherRankParent)
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
        #pragma omp parallel num_threads(numThreads)
        {
            uint32_t numUnions=1;
            int numRandAccesses = 0;
            int tid = omp_get_thread_num();
            uint64_t startOfs = (uint64_t)tid * numReads / numThreads;
            uint64_t endOfs = (uint64_t)(tid+1) * numReads / numThreads;
	    uint64_t writeOfs = startOfs;
            uint64_t i=0;
	    int iter = 0;
    
//            while (numUnions > 0)
            {
                numUnions = 0;
       	        for (i=startOfs; i<endOfs; i++)
                {
                    uint64_t thisParent = (uint64_t)(*((uint32_t *)&parent[i*5 + 1])) + (uint64_t)parent[i*5]*UINT32_MAX; 
                    uint64_t otherParent = (uint64_t)(*((uint32_t *)&otherRankParent[i*5 + 1])) + (uint64_t)otherRankParent[i*5]*UINT32_MAX; 
                    if ((thisParent==otherParent) || (otherParent == 0xBADBADBADl) || (thisParent > numReads) || (otherParent > numReads))
                    {
/*			uint32_t numread_rem = (uint32_t)(numReads % UINT32_MAX);
			uint8_t numread_quot = (uint8_t)(numReads / UINT32_MAX);
                        memcpy (&otherRankParent[i*5 +1], &numread_rem, sizeof(uint32_t));
			otherRankParent[i*5] = numread_quot;*/
			if (((thisParent > numReads) && (thisParent != 0xBADBADBADl)) || ((otherParent > numReads) && (otherParent != 0xBADBADBADl)))
			{
				printf ("[%d]ERROR: mergeStep parent[%lu]=%lu, otherRankParent[%lu]=%lu numReads %lu\n", rank, i, thisParent, i, otherParent, numReads);
			}
                        continue;
                    }
		    if  (thisParent == 0xBADBADBADl)
		    {
//                        uint32_t otherParentRem = (uint32_t)(otherParent % UINT32_MAX);
//			uint8_t otherParentQuot = (uint8_t)(otherParent / UINT32_MAX);
//			memcpy (&parent[i*5 + 1], &otherParentRem, sizeof(uint32_t));
//			parent[i*5] = otherParentQuot;
//			continue;
                    }
//		    assert ((i != 0xBADBADBADl) && (otherParent != 0xBADBADBADl));
                    uint64_t p_u = findRoot (parent, thisParent, &numRandAccesses);
                    uint64_t p_v = findRoot (parent, otherParent, &numRandAccesses);
//		    assert ((p_u != 0xBADBADBADl) && (p_v != 0xBADBADBADl)) printf ("merge - %llx %llx\n", p_u, p_v);
                    if (p_u != p_v)
                    {
			if (writeOfs < (i-2))
			{
				uint32_t u_rem = (uint32_t)(p_u % UINT32_MAX);
				uint8_t u_quot = (uint8_t)(p_u / UINT32_MAX);
				uint32_t v_rem = (uint32_t)(p_v % UINT32_MAX);
				uint8_t v_quot = (uint8_t)(p_v / UINT32_MAX);
				memcpy (&otherRankParent[writeOfs*5 +1], &u_rem, sizeof(uint32_t));
				otherRankParent[writeOfs*5] = u_quot;
				memcpy (&otherRankParent[(writeOfs+1)*5 +1], &v_rem, sizeof(uint32_t));
				otherRankParent[(writeOfs+1)*5] = v_quot;
				writeOfs+=2;
			}
                        if (p_u > p_v)
			{
                			unionOp (parent, compSizes, p_u, p_v);
//                			unionOp (otherRankParent, compSizes, p_u, p_v);
			}
                        else
			{
                            unionOp (parent, compSizes, p_v, p_u);
//                            unionOp (otherRankParent, compSizes, p_v, p_u);
			}
                        numUnions++;
                    }
                    else
                    {
/*			uint32_t numread_rem = (uint32_t)(numReads % UINT32_MAX);
			uint8_t numread_quot = (uint8_t)(numReads / UINT32_MAX);
                        memcpy (&otherRankParent[i*5 +1], &numread_rem, sizeof(uint32_t));
			otherRankParent[i*5] = numread_quot;*/
//                        otherRankParent[i] = numReads;
                    }
                }
		if ((rank == 0) && (numUnions > 0))
			printf ("Rank %d Tid %d mergeStep[%d]: numUnions %d NewRange[%lu,%lu)\n", rank, tid, iter++, numUnions, startOfs, writeOfs);
            }
/*	    while (startOfs < writeOfs)
	    {
		endOfs = writeOfs;
		writeOfs = startOfs;
		numUnions = 0;
		for (i=startOfs; i<endOfs; i+=2)
		{
                    uint64_t thisParent = (uint64_t)(*((uint32_t *)&otherRankParent[i*5 + 1])) + (uint64_t)otherRankParent[i*5]*UINT32_MAX; 
                    uint64_t otherParent = (uint64_t)(*((uint32_t *)&otherRankParent[(i+1)*5 + 1])) + (uint64_t)otherRankParent[(i+1)*5]*UINT32_MAX; 
                    uint64_t p_u = findRoot (parent, thisParent, &numRandAccesses);
                    uint64_t p_v = findRoot (parent, otherParent, &numRandAccesses);
                    if (p_u != p_v)
                    {
			if (writeOfs < (i-2))
			{
				uint32_t u_rem = (uint32_t)(p_u % UINT32_MAX);
				uint8_t u_quot = (uint8_t)(p_u / UINT32_MAX);
				uint32_t v_rem = (uint32_t)(p_v % UINT32_MAX);
				uint8_t v_quot = (uint8_t)(p_v / UINT32_MAX);
				memcpy (&otherRankParent[writeOfs*5 +1], &u_rem, sizeof(uint32_t));
				otherRankParent[writeOfs*5] = u_quot;
				memcpy (&otherRankParent[(writeOfs+1)*5 +1], &v_rem, sizeof(uint32_t));
				otherRankParent[(writeOfs+1)*5] = v_quot;
				writeOfs+=2;
			}
                        if (p_u > p_v)
                			unionOp (parent, compSizes, p_u, p_v);
                        else
                            unionOp (parent, compSizes, p_v, p_u);
                        numUnions++;
                    }
		}
		if ((rank == 0) && (numUnions > 0))
			printf ("Rank %d Tid %d mergeStep[%d]: numUnions %d NewRange[%lu,%lu)\n", rank, tid, iter++, numUnions, startOfs, writeOfs);
            }*/
        }
	gettimeofday (&eTime, NULL);
	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	printf ("[%d]UnionFindMerge %lf sec\n", rank, (double) elapsed/1000000/*, numUnions*/);
//	unionFindTime1 += elapsed;
}

uint64_t printStats()
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
/*#if RADIX_OPT
	FILE *largeFp = fopen("large.opt", "w");
#else
	FILE *largeFp = fopen("large.noopt", "w");
#endif
	assert (largeFp != NULL);*/
//    bzero (compSizes, numReads*sizeof(uint32_t));
	uint64_t readNo = 0, countReads=0, uniqComps=0;
	uint64_t largestCompSize=0, largestCompId=0, badReads=0;
    int numRandAccesses=0;
//	#pragma omp parallel for schedule(static)//private(readNo) //num_threads (numThreads)
/*	for (readNo=0; readNo<numReads; readNo++)
	{
				if (parent[readNo] == 0) 
					printf ("p[%u]=0 compSizes[0]=%lu\n", readNo, compSizes[0]);
        }*/
	uint64_t numPeReadsSep = 0;
	for (readNo=0; readNo<numReads; readNo++)
	{
                uint64_t thisParent = (uint64_t)(*((uint32_t *)&parent[readNo*5 + 1])) + (uint64_t)parent[readNo*5]*UINT32_MAX; 
		if (thisParent != 0xBADBADBADl)
		{
			if ((readNo==14) || (readNo == 17))
				printf ("Bef p[%lu]=%lu\n", readNo, thisParent);
			assert (readNo != 0xBADBADBADl);
			uint64_t newParent=findRoot(parent, readNo, &numRandAccesses);
			uint32_t newParentRem = (uint32_t)(newParent % UINT32_MAX);
			uint8_t newParentQuot = (uint8_t)(newParent / UINT32_MAX);
			memcpy (&parent[readNo*5 + 1], &newParentRem, sizeof(uint32_t));
			parent[readNo*5] = newParentQuot;
			thisParent = newParent;
			if ((readNo==14) || (readNo == 17))
				printf ("Aft p[%lu]=%lu\n", readNo, thisParent);
		}
		if (thisParent == 0xBADBADBADl) badReads++;
		if (thisParent == readNo)
		{
			uniqComps++;
//			countReads += compSizes[readNo];
//			if (compSizes[readNo]==1)
//				fprintf (singleNew, "%u\n", readNo);
		}	
		else
		{
			if (thisParent != 0xBADBADBADl)
			{
				compSizes[readNo]=0;
				compSizes[thisParent]++;
				if ((compSizes[thisParent] > largestCompSize)/* && (parent[readNo] != 0xBADBADBu)*/)
				{
					largestCompSize = compSizes[thisParent] ;
					largestCompId = thisParent;
//					printf ("largestCompId %u largestCompSize %u\n", largestCompId, largestCompSize);
				}
			}
		}
	}
	for (readNo=0; readNo<numReads; readNo+=2)
	{
                uint64_t thisParent = (uint64_t)(*((uint32_t *)&parent[readNo*5 + 1])) + (uint64_t)parent[readNo*5]*UINT32_MAX; 
                uint64_t mateParent = (uint64_t)(*((uint32_t *)&parent[(readNo+1)*5 + 1])) + (uint64_t)parent[(readNo+1)*5]*UINT32_MAX; 
		if ((thisParent != mateParent) && ((thisParent != readNo) || (mateParent != readNo+1))) numPeReadsSep++;
	}
	printf ("[%d] BadReads %u largestCompId %lu largestCompSize %u numPeReadsInSepComps %u\n", rank, badReads, largestCompId, compSizes[largestCompId], numPeReadsSep);
/*	for (readNo=0; readNo<numReads; readNo++)
	{
		if (compSizes[parent[readNo]] > 250)
		{
			parent[readNo]=largestCompId;
		}
	}*/
//	fclose (largeFp);
//	printf ("[%d]Number of COmponents %u, countReads %u, numReads %u LargestCompSize %u, LargestCompId %u\n", rank, uniqComps, countReads, numReads, largestCompSize, largestCompId);
//	assert (countReads == numReads);
	gettimeofday (&eTime, NULL);
	unionFindTime2 += (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	return largestCompId;
}

void printStats2 ()
{
	struct timeval sTime, eTime;
	uint64_t readNo = 0, countReads=0, uniqComps=0;
	gettimeofday (&sTime, NULL);
	qsort (compSizes, numReads, sizeof(uint32_t), parentCompare);
	printf ("[%d]largestComp Size %u\n", rank, compSizes[numReads-1]);
	uint32_t curCompSize=0, curCompSizeCount=0, numComps=0;
	for (readNo=0; readNo<numReads; readNo++)
	{
		if (compSizes[readNo] != curCompSize)
		{
			if (curCompSize > 0)
			{
				printf ("CompSize %u NumComponents %u\n", curCompSize, curCompSizeCount);
				numComps+=curCompSizeCount;
			}
			curCompSize = compSizes[readNo];
			curCompSizeCount=0;
		}
		curCompSizeCount++;
	}
	if (curCompSize > 0)
	{
		printf ("CompSize %lu NumComponents %lu\n", curCompSize, curCompSizeCount);
		numComps+=curCompSizeCount;
	}
//	assert (numComps == uniqComps);
	if (kmerDegree != NULL)
	{
		FILE *fpKmer = fopen("KmerDegDistr", "w");
		assert (fpKmer != NULL);
		for (readNo=0; readNo<numReads; readNo++)
			if (kmerDegree[readNo]>0)
				fprintf (fpKmer, "%u\t%u\n", readNo, kmerDegree[readNo]);
		fclose (fpKmer);
	}
	gettimeofday (&eTime, NULL);
//	unionFindTime2 += (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
}
#if USE_MPI
void mergeCommOpt ()
{
    struct timeval sTime, eTime;
	uint64_t i;
	int tasks=numTasks-1;
    int commIters = (int)(numReads / ASYNC_SIZE);
	uint64_t j1=0;
    if ((numReads % ASYNC_SIZE) != 0) commIters++;
	    gettimeofday (&sTime, NULL);
	if (rank == 0)
		memcpy (otherRanksParent, parent, numReads*5);
            for (j1=0; j1<(commIters-1); j1++)
	    {
//		MPI_Bcast (&otherRanksParent[j1*ASYNC_SIZE*5], ASYNC_SIZE*5, MPI_BYTE, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
	    		gettimeofday (&eTime, NULL);
		    	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
			printf ("[%d] mergeCommOpt Bcast iter %d time %lf sec\n", rank, j1, (double)elapsed/1000000);
		}
	    }
//    	    MPI_Bcast (&otherRanksParent[j1*ASYNC_SIZE*5], (int)(numReads%ASYNC_SIZE)*5, MPI_BYTE, 0, MPI_COMM_WORLD);
	    if (rank == 0)
	    {
    		gettimeofday (&eTime, NULL);
	    	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
		printf ("[%d] mergeCommOpt Bcast iter %d time %lf sec\n", rank, j1, (double)elapsed/1000000);
	    }
	uint64_t globalOffset =0;
	uint64_t *edgeOfs = (uint64_t *) calloc (numTasks+2, sizeof(uint64_t));
	uint64_t *numEdges = (uint64_t *) calloc (numTasks, sizeof(uint64_t));
	assert ((edgeOfs != NULL) && (numEdges != NULL));
	if (rank > 0)
	{
		globalOffset = mergeStepOpt(otherRanksParent);
	}
	MPI_Gather (&globalOffset, 1, MPI_UNSIGNED_LONG, numEdges, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		for (i=1; i<=numTasks; i++)
			edgeOfs[i] = edgeOfs[i-1]+numEdges[i-1];
		for (i=0; i<=numTasks; i++)
		{
			printf ("EDGEOFS [%lu]=%lu\n", i, edgeOfs[i]);
		}
	}
	
	free (numEdges);
	free (edgeOfs);
//			mergeStep(otherRanksParent);
//			printStats();
//		MPI_Barrier(MPI_COMM_WORLD);
	    gettimeofday (&eTime, NULL);
	    long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
        printf ("[%d] MergeOpt Time %lf sec\n", rank,  (double)elapsed/1000000);
        mergeTime += elapsed;
}
void mergeComm ()
{
    struct timeval sTime, eTime, sTime1;
	uint64_t i;
	int numMergeIters=0;
	int tasks=numTasks-1;
	while (tasks > 0)
	{
		numMergeIters ++;
		tasks /= 2;
	}
	if (rank == 0)
		printf ("[%d] numMergeIters %d\n", rank, numMergeIters);
#if MERGE_COMM_OPT
	uint64_t numReads1 = numReads8Mul/8;
	int commIters = (int)(numReads1 / MERGE_COMM_SIZE);
	if ((numReads1 % MERGE_COMM_SIZE) != 0) commIters++;
	uint64_t j1=0;
	for (i=0; i<numMergeIters; i++)
	{
	    gettimeofday (&sTime, NULL);
		int commDistance = 1<< i;
		int processDivider = commDistance<<1;
		if (rank % processDivider == commDistance)
		{
			int toRank = rank - commDistance;
	            for (j1=0; j1<(commIters-1); j1++)
		    {
			printf ("MIter %d, CIter %lu, SendRank %d Ofs %lu Size %d to rank %d\n", i, j1, rank, j1*MERGE_COMM_SIZE*8, MERGE_COMM_SIZE, toRank);
    			MPI_Send (&parent[j1*MERGE_COMM_SIZE*8], MERGE_COMM_SIZE, MPI_UNSIGNED_LONG, toRank, j1, MPI_COMM_WORLD);
		    }
			printf ("MIter %d, CIter %lu, SendRank %d Ofs %lu Size %d to rank %d\n", i, j1, rank, j1*MERGE_COMM_SIZE*8, (int)(numReads1-j1*MERGE_COMM_SIZE), toRank);
    			MPI_Send (&parent[j1*MERGE_COMM_SIZE*8], (int)(numReads1-j1*MERGE_COMM_SIZE), MPI_UNSIGNED_LONG, toRank, j1, MPI_COMM_WORLD);
/*			if (toRank == 0)
			{*/
				for (j1=0; j1<numReads; j1++)
				{
					uint64_t thisParent = (uint64_t)(*((uint32_t *)&parent[j1*5 + 1])) + (uint64_t)parent[j1*5]*UINT32_MAX;
					if ((thisParent != j1) && (thisParent != 0xBADBADBADl))
						printf ("ERROR: At Rank %d Send To %d j1 %lu thisParent %lu (%llx) \n", rank, toRank, j1, thisParent, thisParent);
					assert ((thisParent == j1) || (thisParent == 0xBADBADBADl));
				}
		/*	}*/
		}
		else if (rank % processDivider == 0)
		{
			MPI_Status status;
//			assert (numReads <= 2147483647u);
		        for (j1=0; j1<(commIters-1); j1++)
		 	{
			    printf ("MIter %d, CIter %lu, RecvRank %d Ofs %lu Size %d from rank %d\n", i, j1, rank, j1*MERGE_COMM_SIZE*8, MERGE_COMM_SIZE, rank+commDistance);
			    MPI_Recv (&otherRanksParent[j1*MERGE_COMM_SIZE*8], MERGE_COMM_SIZE, MPI_UNSIGNED_LONG, rank+commDistance, j1, MPI_COMM_WORLD, &status);
			}
			printf ("MIter %d, CIter %lu, RecvRank %d Ofs %lu Size %d from rank %d\n", i, j1, rank, j1*MERGE_COMM_SIZE*8, (int)(numReads1-j1*MERGE_COMM_SIZE), rank+commDistance);
			MPI_Recv (&otherRanksParent[j1*MERGE_COMM_SIZE*8], (int)(numReads1-j1*MERGE_COMM_SIZE), MPI_UNSIGNED_LONG, rank+commDistance, j1, MPI_COMM_WORLD, &status);
/*			if (rank == 0)
			{*/
				for (j1=0; j1<numReads; j1++)
				{
					uint64_t thisParent = (uint64_t)(*((uint32_t *)&otherRanksParent[j1*5 + 1])) + (uint64_t)otherRanksParent[j1*5]*UINT32_MAX;
					if ((thisParent != j1) && (thisParent != 0xBADBADBADl))
						printf ("ERROR: At Rank %d Recv From %d j1 %lu thisParent %lu (%llx) \n", rank, rank+commDistance, j1, thisParent, thisParent);
					assert ((thisParent == j1) || (thisParent == 0xBADBADBADl));
//					printf ("MIter %d RECV[%lu]=%lx\t", i, j1, thisParent);
				}
//			}
//			mergeStep(otherRanksParent);
//			printStats();
		}
//		MPI_Barrier(MPI_COMM_WORLD);
	    gettimeofday (&eTime, NULL);
	    long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
        printf ("[%d] Merge Iter %d Time %lf sec\n", rank, i, (double)elapsed/1000000);
        mergeTime += elapsed;
	}
#else
    int commIters = (int)(numReads / ASYNC_SIZE);
	uint64_t j1=0;
    if ((numReads % ASYNC_SIZE) != 0) commIters++;
	for (i=0; i<numMergeIters; i++)
	{
	    gettimeofday (&sTime, NULL);
		int commDistance = 1<< i;
		int processDivider = commDistance<<1;
		if (rank % processDivider == commDistance)
		{
			int toRank = rank - commDistance;
			printf ("Iter %d, Rank %d to rank %d\n", i, rank, toRank);
	            for (j1=0; j1<(commIters-1); j1++)
    			MPI_Send (&parent[j1*ASYNC_SIZE*5], ASYNC_SIZE*5, MPI_BYTE, toRank, j1, MPI_COMM_WORLD);
    		MPI_Send (&parent[j1*ASYNC_SIZE*5], (int)(numReads%ASYNC_SIZE)*5, MPI_BYTE, toRank, j1, MPI_COMM_WORLD);
		}
		else if (rank % processDivider == 0)
		{
			MPI_Status status;
//			assert (numReads <= 2147483647u);
            for (j1=0; j1<(commIters-1); j1++)
			{
				
				gettimeofday (&sTime1, NULL);
			    MPI_Recv (&otherRanksParent[j1*ASYNC_SIZE*5], ASYNC_SIZE*5, MPI_BYTE, rank+commDistance, j1, MPI_COMM_WORLD, &status);
				gettimeofday (&eTime, NULL);
				if (rank == 0)
				{
					long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime1.tv_sec * 1000000 + sTime1.tv_usec);
					printf ("Rank 0 MIter %d CIter %d RECVTIME %lf sec %d bytes received from %d\n", i, j1, (double)elapsed/1000000, ASYNC_SIZE*5, rank+commDistance);
				}
			}
			gettimeofday (&sTime1, NULL);
			MPI_Recv (&otherRanksParent[j1*ASYNC_SIZE*5], (int)(numReads%ASYNC_SIZE)*5, MPI_BYTE, rank+commDistance, j1, MPI_COMM_WORLD, &status);
				gettimeofday (&eTime, NULL);
				if (rank == 0)
				{
					long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime1.tv_sec * 1000000 + sTime1.tv_usec);
					printf ("Rank 0 MIter %d CIter %d RECVTIME %lf sec %d bytes received from %d\n", i, j1, (double)elapsed/1000000, (int)(numReads%ASYNC_SIZE)*5, rank+commDistance);
				}
			gettimeofday (&sTime1, NULL);
			mergeStep(otherRanksParent);
				gettimeofday (&eTime, NULL);
				if (rank == 0)
				{
					long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime1.tv_sec * 1000000 + sTime1.tv_usec);
					printf ("Rank 0 MIter %d CIter %d MERGETIME %lf sec\n", i, j1, (double)elapsed/1000000);
				}
//			printStats();
		}
//		MPI_Barrier(MPI_COMM_WORLD);
	    gettimeofday (&eTime, NULL);
	    long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
        printf ("[%d] Merge Iter %d Time %lf sec\n", rank, i, (double)elapsed/1000000);
        mergeTime += elapsed;
	}
#endif

}
#endif
int main(int argc, char **argv)
{
	char hostname[256];

	int i=0, j=0, k=0;
	static struct option long_options[] = {
		{"out-prefix", required_argument, 0, 'o'},
		{0, 0, 0, 0}
	};
	int c=0;
	while (1)
	{
		int option_index = 0;
		c = getopt_long (argc, argv, "o:", long_options, &option_index);
		if (c==-1) 
		{
			break;
		}
		switch (c)
		{
			case 'o':
				strcpy (oPrefix, optarg);
				printf ("optarg %s\n", optarg);
				break;
		}
	}
    int providedThreadSupport;
#if USE_MPI
	MPI_Init_thread(&argc,&argv, MPI_THREAD_MULTIPLE, &providedThreadSupport);
	MPI_Comm_size(MPI_COMM_WORLD,&numTasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if (rank==0)
        printf ("[%d] providedThreadSupport %d\n", rank, providedThreadSupport);
#else
	numTasks = 1;
	rank = 0;
#endif
	assert (numTasks < 1024);
	gethostname(hostname, sizeof(hostname));
	printf("[%d]PID %d on %s ready for attach\n", rank, getpid(), hostname);
//	while (0 == i)
//           sleep(5);

    if (rank==0)
    	printf ("numTasks %d\n", numTasks);
	struct stat buf;
	hwThreads = hardware_threads();
    if (rank==0)
    	printf ("HwThreads %d argc %d optind %d\n", hwThreads, argc, optind);
	if (argc < 11)
	{
		printf ("Usage: %s -o <output-prefix> fastqBufferSize kmerSize numScans numThreads "
				"numPEFiles numSEFiles numInterleavedFiles a.fastq "
				"b.fastq ...\n", argv[0]);
		exit(1);
	}
	rangeHist = (uint64_t *)malloc((1048576 +1) * sizeof(uint64_t));
	assert (rangeHist != NULL);
	
	argv = &argv[optind-1];
	fastqBufferSize = (unsigned long)atol(argv[1]);
	K = atoi(argv[2]);
	getKmerMask (K);
//	if (rank == 0)
//		printf ("kmerMask %" PRIu64 ", rotater %d\n", kmerMask, rotater);
	scans = atoi (argv[3]);
	numThreads = atoi(argv[4]);
	threadKmerGenTime = (long *)malloc (numThreads * sizeof (long));
	assert (threadKmerGenTime != NULL);

	numPEFiles = atoi(argv[5])*2;
	numSEFiles = atoi(argv[6]);
	numInterleavedFiles = atoi(argv[7]);
	argv = &argv[7];
	argc-=(7+optind-1);
//	printf ("fastqBufferSize %lu K %d scans %d numPEFiles %d numSEFiles %d optind %d\n", fastqBufferSize, K, scans, numPEFiles, numSEFiles, optind);
	init(K, numTasks, numThreads, 10);
//	printf ("Argc %d\n", argc);
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
		assert (status == 0);
		fileSize[i-1]=buf.st_size;
//		if (rank==0)
//			printf ("file %s Size %lu\n", fileNames[i-1], buf.st_size);
	}
	fastqRanges = (int *)malloc ((numThreads+1) * sizeof(int));
	assert (fastqRanges != NULL);

	char idxFile[200];
#ifdef MIN3
	sprintf (idxFile, "%s.KmerIdx_3min%d", argv[1], K);
#else
	#ifdef ALLMINS
		sprintf (idxFile, "%s.KmerIdx_allmins%d", argv[1], K);
	#else
		sprintf (idxFile, "%s.KmerIdx%d", argv[1], K);
	#endif
#endif
    if (rank==0)
    	printf ("idxFile %s\n", idxFile);
	int status = stat (idxFile, &buf);
	assert (status == 0);
	totalFastqPartitions = buf.st_size/sizeof(fileIndex) -1;
	printf ("totalFastqPartitions %d numTasks %d numThreads %d, sizeof(fileIndex) %lu\n", totalFastqPartitions, numTasks, numThreads, sizeof(fileIndex));
	assert (totalFastqPartitions % (numTasks*numThreads) == 0);
	int fastqPartitionsPerRank = totalFastqPartitions/numTasks;
	int fastqPartitionsPerThread = totalFastqPartitions/numTasks/numThreads;
//	printf ("fastqPartitionsPerThread %d\n", fastqPartitionsPerThread);

	fI = (fileIndex *)malloc((fastqPartitionsPerRank + 1) * sizeof(fileIndex));
	assert (fI != NULL);
	FILE *fp_idx = fopen (idxFile, "r");
	assert (fp_idx != NULL);
	assert (fseek (fp_idx, rank*fastqPartitionsPerRank*sizeof(fileIndex), SEEK_SET) == 0);
	assert (fread(fI, sizeof(fileIndex), fastqPartitionsPerRank+1, fp_idx) == (fastqPartitionsPerRank+1));

	fileIndex lastIdx;
	assert (fseek (fp_idx, totalFastqPartitions*sizeof(fileIndex), SEEK_SET) == 0);
	assert (fread(&lastIdx, sizeof(fileIndex), 1, fp_idx) == 1);
	numReads = (uint64_t)lastIdx.readIdExtra*UINT32_MAX + lastIdx.readId;
	fclose (fp_idx);
	numReads8Mul = numReads*5;
	while (numReads8Mul%sizeof(uint64_t) != 0) numReads8Mul++;

	assert (numReads8Mul % 8 == 0);

	
    if (rank==0)
        printf ("[%d] NUMREADS %lu numReads8Mul %lu\n", rank, numReads, numReads8Mul);
/*	openFp = (FILE **)malloc(2 * numThreads * sizeof(FILE *));
	openFileNo = (int *)malloc (2 * numThreads * sizeof(int));
	assert (openFp != NULL);*/

	for (i=0; i<=numThreads; i++)
	{
		fastqRanges[i]=/*rank*numThreads*fastqPartitionsPerThread +*/ i*fastqPartitionsPerThread;
/*		if (i < numThreads)
		{
			openFp[2*i] = NULL; openFp[2*i + 1]=NULL;
			openFileNo[2*i] = -1; openFileNo[2*i + 1] = -1;
		}*/
/*		if (rank==0)
			printf ("[%d]fastqRanges[%d] = %d\n", rank, i, fastqRanges[i]);*/
	}
//	fastqRanges[numThreads]=totalFastqPartitions;
//	if (rank==0)
//		printf ("fastqRanges[%d] = %d\n", numThreads, fastqRanges[numThreads]);

	/*Range Hist*/
#ifdef MIN3
	sprintf (idxFile, "%s.kmerHist_3min%d", argv[1], K);
#else
	#ifdef ALLMINS
		sprintf (idxFile, "%s.kmerHist_allmins%d", argv[1], K);
	#else
		sprintf (idxFile, "%s.kmerHist%d", argv[1], K);
	#endif
#endif
	FILE *fp_rangeHist = fopen(idxFile, "r");
	assert (fp_rangeHist != NULL);
	assert (fread(rangeHist, sizeof(uint64_t), 1048576+1, fp_rangeHist) == (1048576+1));
	fclose (fp_rangeHist);
//	for (i=0; i<=94588; i++)
//		printf ("rangehist[%d]=%lu\n", i, rangeHist[i]);
	if (rank==0)
		printf ("TotalKmers %lu\n", rangeHist[1048576]);

	scanKmerRanges = (uint64_t *)calloc((scans+1) , sizeof(uint64_t));
	processKmerRanges = (uint64_t *)calloc(scans*numTasks +1 , sizeof(uint64_t));
	threadKmerRanges = (uint64_t *)calloc((scans*numTasks*numThreads + 1) , sizeof(uint64_t));
	assert ((scanKmerRanges != NULL) && (processKmerRanges != NULL) && (threadKmerRanges != NULL));
	uint64_t avgKmersPerScan = rangeHist[1048576]/scans;
	if (rangeHist[1048576] % scans != 0) avgKmersPerScan++;
//	printf ("avgKmersPerScan %lu\n", avgKmersPerScan);
	int curRange=1;
	scanKmerRanges[0]=0;
	
	for (i=0; i<1048576; i++)
	{
		if (rangeHist[i] >= avgKmersPerScan*curRange)
		{
			printf ("curRange %d avgKmersPerScan %d curRange*avgKmersPerScan %lu rangeHist[%llx]=%lu\n", curRange, avgKmersPerScan, avgKmersPerScan*curRange, i, rangeHist[i]);
			scanKmerRanges[curRange]=i;
			curRange++;
		}
	}
//	if (rangeHist[i] >= avgKmersPerScan*curRange)
	if (curRange == scans)
	{
//		printf ("scanKmerRanges[%d]=%u\n", curRange, i);
		scanKmerRanges[curRange]=i;
	}
	for (i=0; i<=scans; i++)
	{
//		printf ("curRange %d\n", curRange);
		printf ("[%d]  scanKmerRanges[%d]=%lu\n", rank, i, scanKmerRanges[i]);
	}
	uint64_t avgKmersPerProcessPerScan = avgKmersPerScan/numTasks;
	
	curRange=1;
	for (i=0; i<1048576; i++)
	{
		if (rangeHist[i] >= curRange*avgKmersPerProcessPerScan)
		{
			if ((curRange % numTasks) == 0)
				processKmerRanges[curRange]=scanKmerRanges[curRange/numTasks];
			else
				processKmerRanges[curRange]=i;
			curRange++;
		}
	}
//	if (rangeHist[i] >= curRange*avgKmersPerProcessPerScan)
	if (curRange == scans*numTasks)
	{
		processKmerRanges[curRange]=i;
	}
	for (i=0; i<=scans*numTasks; i++)
		printf ("[%d]  processKmerRanges[%d]=%lu\n", rank, i, processKmerRanges[i]);
//	printf ("avgKmersPerProcessPerScan %lu/%lu\n", avgKmersPerProcessPerScan, rangeHist[1048576]);
	uint64_t avgKmersPerThreadPerScan = avgKmersPerProcessPerScan/numThreads;
	curRange=1;
	for (i=0; i<1048576; i++)
	{
		if (rangeHist[i] >= curRange*avgKmersPerThreadPerScan)
		{
			if ((curRange % numThreads) == 0)
				threadKmerRanges[curRange]=processKmerRanges[curRange/numThreads];
			else
				threadKmerRanges[curRange]=i;
			curRange++;
		}
	}
//	if (rangeHist[i] >= curRange*avgKmersPerThreadPerScan)
	if (curRange == scans*numTasks*numThreads)
	{
		threadKmerRanges[curRange]=i;
	}
	for (i=0; i<=scans*numTasks*numThreads; i++)
		printf ("[%d]  threadKmerRanges[%d]=%lu\n", rank, i, threadKmerRanges[i]);

	sendSizes = (size_t *)calloc(scans*numTasks , sizeof(size_t));
	assert (sendSizes != NULL);
	sendThreadSizes = (size_t *)calloc(scans*numTasks*numThreads , sizeof(size_t));
	assert (sendThreadSizes != NULL);
	
//	printf ("[%d] fastqRanges[%d]=%d\n", rank, numThreads, fastqRanges[numThreads]);
	int tid=0;
	size_t *sendSizeInScan = (size_t *)calloc(scans, sizeof(size_t));
	assert (sendSizeInScan != NULL);
	for (tid=0; tid<numThreads; tid++)
	{
		for (i=fastqRanges[tid]; i<fastqRanges[tid+1]; i++)
		{
	//		printf ("[%d] Chunk %d\n", rank, i);
			k=0;
			for (j=0; j<scans*numTasks; j++)
			{
				int curScan=j/numTasks;
				int curTask = j%numTasks;
				while (k < processKmerRanges[curScan*numTasks + curTask + 1])
				{
					sendSizes[curScan*numTasks + curTask]+=fI[i].kmerFreqCount[k];
					sendThreadSizes[curScan*numTasks*numThreads + curTask*numThreads + tid]+=fI[i].kmerFreqCount[k];
					sendSizeInScan[curScan] += fI[i].kmerFreqCount[k];
					k++;
				}
			}
		}
	}
	size_t maxSendSize=0, maxRcvSize=0;
	for (i=0; i<scans; i++)
	{
		if (sendSizeInScan[i] > maxSendSize)
			maxSendSize = sendSizeInScan[i];
		printf ("[%d] send %lu kmers in Scan %d\n", rank, sendSizeInScan[i], i);
	}
	free (sendSizeInScan);
	threadKmerCounts = (uint64_t *)calloc (numThreads*1048576, sizeof(uint64_t));
	assert (threadKmerCounts  != NULL);
	for (i=0; i<fastqPartitionsPerRank; i++)
	{
		int curThread = ((rank * fastqPartitionsPerRank) + i) * numThreads / totalFastqPartitions;
//		printf ("i %d curThread %d\n", i, curThread);
		for (j=0; j<1048576; j++)
		{
			threadKmerCounts[curThread*1048576 + j] += fI[i].kmerFreqCount[j];
		}
	}
#if USE_MPI
	MPI_Allreduce (MPI_IN_PLACE, threadKmerCounts, numThreads*1048576, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

	size_t *rcvSizeInScan = (size_t *)calloc(scans, sizeof(size_t));
	assert (rcvSizeInScan != NULL);
/*	size_t *rcvSizeInScanCpy = (size_t *)calloc(scans, sizeof(size_t));
	assert (rcvSizeInScanCpy != NULL);
	for (i=0; i<totalFastqPartitions; i++)
	{
		for (j=0; j<scans; j++)
		{
			for (k=processKmerRanges[j*numTasks + rank]; k<processKmerRanges[j*numTasks + rank + 1]; k++)
			{
				rcvSizeInScan[j]+=fI[i].kmerFreqCount[k];
			}
//			if (rcvSizeInScan[j] > maxRcvSize)
//				maxRcvSize = rcvSizeInScan[j];
		}
	}*/
	for (i=0; i<numThreads; i++)
	{
		for (j=0; j<scans; j++)
		{
			for (k=processKmerRanges[j*numTasks + rank]; k<processKmerRanges[j*numTasks + rank + 1]; k++)
			{
				rcvSizeInScan[j]+=threadKmerCounts[i*1048576 + k];
			}
			if (rcvSizeInScan[j] > maxRcvSize)
				maxRcvSize = rcvSizeInScan[j];
		}
	}
	for (i=0; i<scans; i++)
	{
		printf ("[%d] Rcv %lu kmers in Scan %d\n", rank, rcvSizeInScan[i], i);
//		assert (rcvSizeInScan[i] == rcvSizeInScanCpy[j]);
	}
//	free (rcvSizeInScanCpy);
	free (rcvSizeInScan);
	if (maxRcvSize > maxSendSize)
		maxSendSize = maxRcvSize;
	commonSendBuffer = (uint16_t *)malloc (maxSendSize * 7 * sizeof(uint16_t));
	assert (commonSendBuffer != NULL);
	if (numTasks > 1)
	{
		rcvKmers = (uint16_t *)malloc (maxRcvSize * 7 * sizeof (uint16_t));
		assert (rcvKmers != NULL);
		sortedKmers = commonSendBuffer;
	}
	else
	{
		sortedKmers = (uint16_t *)malloc (maxRcvSize * 7 * sizeof (uint16_t));
		assert (sortedKmers != NULL);
		rcvKmers = commonSendBuffer;
	}
	printf ("[%d]maxSendSize %lu, maxRcvSize %lu commonSendBuffer %p, rcvKmers %p\n", rank, maxSendSize, maxRcvSize, commonSendBuffer, rcvKmers);
/*	for (i=0; i<scans*numTasks; i++)
		printf ("[%d]  sendSizes[%d]=%lu\n", rank, i, sendSizes[i]);
	for (i=0; i<scans*numTasks*numThreads; i++)
		printf ("[%d]  sendThreadSizes[%d]=%lu\n", rank, i, sendThreadSizes[i]);*/

/*	int partPerScan = totalPartitions/scans;
	int partPerRankPerScan = totalPartitions/numTasks/scans;
	int partPerThreadPerScan = partPerRankPerScan/numThreads;
	printf ("PartPerRankPerScan %d, TotalPartitions %d\n", partPerRankPerScan, totalPartitions);
	assert ((((partPerRankPerScan * scans * numTasks) == (totalPartitions-1)) && ((partPerThreadPerScan*numThreads)==partPerRankPerScan));*/
    fastqBuffer1 = (char *)malloc((size_t)numThreads*(fastqBufferSize+10));
    assert (fastqBuffer1 != NULL);
    if (numPEFiles > 0)
    {
        fastqBuffer2 = (char *)malloc((size_t)numThreads*(fastqBufferSize+10));
        assert (fastqBuffer2 != NULL);
    }
//    if (rank == 0)
//        printf ("fastqBuffer1=%p, fastqBuffer2=%p\n", fastqBuffer1, fastqBuffer2);
#if USE_MPI
	int retsize=0;
	MPI_Type_size(MPI_UNSIGNED_SHORT, &retsize);

	printf ("sizeof uint16_t %lu sizeof unsigned sort int %lu\n", sizeof (uint16_t), retsize);
	assert (sizeof (uint16_t) == retsize);
#endif
	for (i=0; i<scans; i++)
	{
		generateKmersInScan(i);
//		if (i == 1) break;
//		printf ("[%d] numReads %u\n", rank, numReads);
#if USE_MPI
//		if (numTasks > 1)	
//			MPI_Bcast (&numReads, 1, MPI_UNSIGNED, numTasks-1, MPI_COMM_WORLD);
#endif
//		printf ("[%d] numReads %u\n", rank, numReads);
		if (parent == NULL)
		{
			parent = (uint8_t *)calloc(numReads8Mul, sizeof(uint8_t));
			assert ((parent != NULL));
			filteredReads = (uint8_t *)calloc(numReads, sizeof(uint8_t));
			assert (filteredReads != NULL);
			compSizes = (uint32_t *)malloc(numReads*5 * sizeof(uint8_t));
			assert (compSizes != NULL);
			uint64_t readNo=0;
			for (readNo=0; readNo<numReads; readNo++)
			{
				uint32_t readNoRem = (uint32_t)(readNo % UINT32_MAX);
				uint8_t readNoQuot = (uint8_t)(readNo / UINT32_MAX);
				memcpy (&parent[readNo*5 + 1], &readNoRem, sizeof(uint32_t));
				parent[readNo*5] = readNoQuot;
				compSizes[readNo]=1;
			}
#ifdef KMERHIST
//			kmerDegree = (uint32_t *)calloc(numReads , sizeof(uint32_t));
//			assert (kmerDegree != NULL);
#endif
		}
		uint64_t *localRangeCounts = (uint64_t *)calloc((numThreads+1), sizeof(uint64_t));
		assert(localRangeCounts != NULL);
		uint64_t *rangeOffsets = (uint64_t *)calloc(numThreads, sizeof(uint64_t));
		assert (rangeOffsets != NULL);	
		sortKmerList (i, localRangeCounts, rangeOffsets);

		markMischiefKmers(i, localRangeCounts, rangeOffsets);
                unionFindPar2(i, localRangeCounts, rangeOffsets);

		unionFind(i, localRangeCounts, rangeOffsets);
		uint64_t maxCompId = printStats();
		if (numTasks > 1)
		{
//			uint32_t *tmp = rcvKmers;
//			rcvKmers = sortedKmers;
//			sortedKmers = tmp;
		}
		free (localRangeCounts);
		free (rangeOffsets);
	}
//    MPI_Barrier(MPI_COMM_WORLD);
	free (commonSendBuffer);
	free (filteredReads);

	if (numTasks > 1)
	{
//		printf ("[%d] before free rcvKmers %p\n", rank, rcvKmers);
		free (rcvKmers);
//		printf ("[%d] after free rcvKmers %p\n", rank, rcvKmers);
	}
	else
	{
//		free (sortedKmers);
	}
	free (sendSizes);
	free (sendThreadSizes);
	free (scanKmerRanges);
	free (processKmerRanges);
	free (threadKmerRanges);
	free (rangeHist);
	free (threadKmerCounts);

	otherRanksParent = (uint8_t *)calloc(numReads8Mul, sizeof(uint8_t));
	assert (otherRanksParent != NULL);
	//mergeComps
	uint64_t readNo=0;
//	FILE *fpComps = fopen ("comps40", "w");
//	assert (fpComps != NULL);
//	fclose (fpComps);
    long writeComm=0;

#if USE_MPI
//	mergeCommOpt();
	mergeComm();
#endif
	for (readNo=0; readNo<numReads; readNo++)
	{
//		uint64_t compId = (uint64_t)(*((uint32_t *)&parent[readNo*5+1])) + (uint64_t)parent[readNo*5]*UINT32_MAX;
//		fprintf (fpComps, "%lu\n", compId);
		compSizes[readNo]=1;
	}
	int numOutPartitions = 16;
	uint64_t maxCompId = 0;
        struct timeval sTime, eTime;
	uint64_t *workPartitions = (uint64_t *)calloc(numOutPartitions+1, sizeof(uint64_t));
	if (rank == 0)
	{
		maxCompId = printStats();
		printf ("MaxCompId %lu(%lx) Size %u\n", maxCompId, maxCompId, compSizes[maxCompId]);
		unionFindPeReads(1000);
//		unionFindSmallComps(5);
		for (readNo=0; readNo<numReads; readNo++)
		{
			compSizes[readNo]=1;
		}
		maxCompId = printStats();
		printf ("MaxCompId %lu(%lx) Size %lu\n", maxCompId, maxCompId, compSizes[maxCompId]);
//		memcpy (otherRanksParent, compSizes, numReads*sizeof(uint64_t));
//		qsort (otherRanksParent, numReads, sizeof(uint64_t), parentCompare);//parentCompare compares uint32_t

//********writing fastq********
		uint64_t numReadsPerProcess = numReads/numOutPartitions;
		workPartitions[numOutPartitions] = numReads;
		uint64_t curSize=0, curTotal=0, prevTotal=0, curInd=0;
		for (i=1; i<numOutPartitions; )
		{
			while ((curInd<numReads) && (curTotal < i*numReadsPerProcess))
			{
				uint64_t curParent = (uint64_t)(*((uint32_t *)&parent[curInd*5 + 1])) + (uint64_t)parent[curInd*5]*UINT32_MAX;
				if (curParent == curInd)
				{
//					printf ("curInd %u, curTotal %u nextTotal %u \n", curInd, curTotal, i*numReadsPerProcess);
					curTotal += compSizes[curParent];
				}
				curInd++;
			}
			printf ("%d: curInd %u, curTotal %u nextTotal %u \n", i, curInd, curTotal, i*numReadsPerProcess);
			workPartitions[i++]=curInd;
		}
		for (i=0; i<=numOutPartitions; i++)
			printf ("%lu\t", workPartitions[i]);
		printf ("\n");
		if (numTasks == 1)
			free (sortedKmers);
	}
#if USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	gettimeofday (&sTime, NULL);
#if USE_MPI
	if (numTasks > 1)
	{
/*#if MERGE_COMM_OPT
		uint64_t numReads1 = numReads8Mul/8;
		int commIters = (int)(numReads1 / MERGE_COMM_SIZE);
		if ((numReads1 % MERGE_COMM_SIZE) != 0) commIters++;
		uint64_t j1=0;
		for (j1=0; j1<(commIters-1); j1++)
			MPI_Bcast (&parent[j1*MERGE_COMM_SIZE*8], MERGE_COMM_SIZE, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast (&parent[j1*MERGE_COMM_SIZE*8], (int)(numReads1 -j1*MERGE_COMM_SIZE), MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast (workPartitions, numOutPartitions+1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	    	gettimeofday (&eTime, NULL);
	        writeComm  = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
#else*/
		int commIters = (int)(numReads / ASYNC_SIZE);
		if ((numReads % ASYNC_SIZE) != 0) commIters++;
		uint64_t j1=0;
		for (j1=0; j1<(commIters-1); j1++)
			MPI_Bcast (&parent[j1*ASYNC_SIZE*5], ASYNC_SIZE*5, MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast (&parent[j1*ASYNC_SIZE*5], (int)(numReads -j1*ASYNC_SIZE)*5, MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast (workPartitions, numOutPartitions+1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	    	gettimeofday (&eTime, NULL);
	        writeComm  = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
//#endif
	}
#endif
	writeOutput3 (workPartitions, numOutPartitions);
//	writeOutput (maxCompId);
	gettimeofday (&eTime, NULL);
	outputTime += (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	if (rank == 0)
	{
		printStats2();
	}

	printf ("[%d] FQRead %lf sec(MaxKmerGen %lf sec), CommTime %lf sec, SortTime %lf sec(%lf + %lf), FilterTime %lf sec, unionFind %lf sec (%lf + %lf), mergeComps %lf sec, peUnionTime %lf sec, WriteFq %lf sec(%lf sec)\n", rank, (double)(kmerGenTime-kmerGenWorkTime)/1000000, (double)kmerGenWorkTime/1000000, (double)commTime/1000000, (double)sortTime1/1000000 + (double)sortTime2/1000000, (double)sortTime1/1000000, (double)sortTime2/1000000, (double)filterTime/1000000, (double)unionFindTime1/1000000 + (double)unionFindTime2/1000000, (double)unionFindTime1/1000000, (double)unionFindTime2/1000000, (double)mergeTime/1000000, (double)peUnionTime/1000000, (double)outputTime/1000000 ,(double)writeComm/1000000);
	free (parent);
	free (otherRanksParent);
	free (compSizes);
	free (fastqBuffer1);
	if (fastqBuffer2 != NULL)
		free (fastqBuffer2);
	free (fastqRanges);
	free (fI);
	for (i=0; i<argc-1; i++) {
		free(fileNames[i]);
	}
	free (fileNames);
	free (fileSize);
	if (recvThreadCount != NULL)
	{
		free (recvThreadCount);
		free (recvThreadOfs);
	}
	
#ifdef KMERHIST
//	free (kmerDegree);
#endif
//	free (openFp);
	free (openFileNo);
	free (threadKmerGenTime);

    if (rank == 0)
    {
    	int pid = getpid();
	    char line[2048];
    	sprintf (line, "/proc/%d/status", pid);
	    FILE *statFile = fopen(line, "r");
    	assert (statFile != NULL);
	    fgets (line, 2048, statFile);
    	while (!feof (statFile))
	    {
		    if (strstr(line,"VmPeak") || strstr(line,"VmHWM"))
    		{
	    		printf ("[%d] %s", rank, line);
		    }
    		fgets (line, 2048, statFile);
	    }
	    fclose (statFile);
    }

#if USE_MPI
	MPI_Finalize();
#endif
	return 0;
}

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
/*#include <emmintrin.h>
#include <tmmintrin.h>
#include <immintrin.h>*/
#include <math.h>
#include <assert.h>
	
//__m128i mask, zeroReg, both_mask2, msb_mask2, lsb_mask2, n_mask;
//__m128i zero32,msb_mask32,lsb_mask32,minus1,mask32,n_mask32;
int numP=0, numThr=0;
uint32_t hiMask=0;
int kmerSize=27;
int mmerSize=10;
//char in_bases[110]="NTGAAAAGTGGGATTCATGATTACTCCTGCCCTGCAGACCCCACACAACTCTTGAGAGGACCAAATTTTAAAAATGCATGTGAAGTCATTGTGAAAACTGT";
//char *in_bases;
char out_bases[110];

static inline uint32_t quicklog2(const uint32_t x) {
  uint32_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
  );
  return y;
}

void getKmerStr64 (uint64_t kmerl, uint64_t kmerh, FILE *fp, int ind)
{
	int i=0;
	if (ind == kmerSize) return;
	char ch = (ind < 32)? ("ACGT"[kmerl & 3]) : ("ACGT"[kmerh & 3]);
	if (ind < 32)
		getKmerStr64 (kmerl >> 2, kmerh, fp, ind+1);
	else
		getKmerStr64 (kmerl, kmerh >> 2, fp, ind+1);
	printf ("%c", ch);
	if (ind==0)
		printf ("\n");
}
void getFullStr (uint64_t kmer, FILE *fp)
{
	int i=0;
	for (i=26; i>=0; i--)
	{
		if (fp!=NULL)
			fprintf (fp, "%c", "ACTG"[(((kmer & 1lu<<(32+i)) >> (31+i)) & 2) | ((kmer & 1lu<<i) >> i)]);
		else
			printf ("%c", "ACTG"[(((kmer & 1lu<<(32+i)) >> (31+i)) & 2) | ((kmer & 1lu<<i) >> i)]);
	}
/*	if (fp != NULL)
		fprintf (fp, "\n");
	else
		printf ("\n");*/
}
void init(int kLen, int numProcs, int numT, int mLen)
{
//	readLen=rLen;
	kmerSize=kLen;
	mmerSize=mLen;
	numP = numProcs;
	numThr = numT;
	assert (numP < 1024);
	hiMask = ~(((1u << (32-kmerSize))-1)<<(kmerSize-16)) & 0xFFFF;
//	printf ("Mask32 %x\n", (1<<kmerSize) -1);
}
typedef struct
{
	uint32_t hi_kmer;
	uint32_t lo_kmer;
//	uint64_t kmer;
	int process;
} kmerMap_t;
typedef struct
{
	uint32_t kmer[5];
	uint8_t adj;
	int process;
} longKmerMap_t;

uint64_t getMmerPrefix (uint64_t kH, uint64_t kL, int m)
{
	uint64_t mmerPrefix=0;
	if (kmerSize > 32)
	{
		if ((kmerSize-32) < m)
		{
			mmerPrefix = kH <<  ((m-(kmerSize-32))<<1);
			int rem = m-(kmerSize-32);
			mmerPrefix |= kL >> (64-(rem << 1));
		}
		else
		{
			mmerPrefix = kH >> ((kmerSize-32-m)<<1);
		}
	}
	else
	{
		mmerPrefix = kL >> ((kmerSize-m) << 1);
	}
	return mmerPrefix;
}

int binSrch (int low, int high, uint64_t *ranges, uint64_t key)
{
	int mid=low;
	while (low < high)
	{
		mid=(low+high)/2;
		if (ranges[mid] > key)
		{
			if ((mid>0) && (ranges[mid-1] <= key)) return mid-1;
			high = mid-1;
		}
		else if (ranges[mid] <= key)
		{
			if ((ranges[mid] == key) || (ranges[mid+1] > key)) return mid;
			low = mid +1;
		}
	}
	assert (0);
	return mid;
}
//return 1 if kR has no space for next read
void getCanonKmers (char *readStr, uint16_t **kR, uint64_t *curSendOfs, int readLen, uint64_t startKmer, uint64_t endKmer, uint64_t kmerMask, int rotater, uint64_t curReadId, uint64_t *scanProcessKmerRange, int tid, int rank/*, uint64_t MAX_TUPLES*/)
{
/*	if (curReadId > 18173140)
	{
		printf ("getCanonKmers for %u readLen %d kR[0]=%p\n", curReadId, readLen, kR[0]);
		readStr[readLen]='\0';
		printf ("%s\n", readStr);
	}*/
	int overflow = 0;
	if ((curReadId == 0) && (tid==0))
		printf ("startKmer %llx endKmer %llx\n", startKmer, endKmer);
	int charToIntMap[4]={0,1,3,2};
	int charToRcMap[4]={3,2,0,1};
	int i=0,j;
	if (readLen <= kmerSize)
		return;
	int lastKmerPos, kmersPerLane;
	lastKmerPos = readLen - kmerSize;
	kmersPerLane = lastKmerPos/4;
	uint64_t curKmerInt[2]={0,0};
	uint64_t curRcInt[2]={0,0};
	int lastNPos = -1;
	uint64_t mask[2];
	int rotateBy[2];
	
	uint8_t lOfs[1024];
        for (i=0; i<numP; i++) lOfs[i]=0;
        longKmerMap_t lBuf[300];
	assert (lastKmerPos < 300);
        int bufOfs=0;

	if (kmerSize > 32)
	{
		if (kmerSize < 64)
			mask[1]=(1lu<<((kmerSize-32)<<1))-1;
		else
			mask[1]=0xFFFFFFFFFFFFFFFFlu;
		mask[0]=0xFFFFFFFFFFFFFFFFlu;
	}
	else
	{
		mask[1]=0;
		mask[0]=(1lu<<(kmerSize << 1)) - 1;
	}
//	printf ("Mask1 %llx Mask0 %llx\n", mask[1], mask[0]);
	
/*	if (curReadId == 18173149)
	{
		printf ("Bef Stage 1\n");
	}*/
	char prevChar='N', nextChar=readStr[kmerSize];
	for (j=0; j<kmerSize; j++)
	{
		int kmerInd = j/32;
		if (readStr[j]=='N')
			lastNPos=j;
		if (kmerInd==1)
		{
			curKmerInt[kmerInd] = (curKmerInt[kmerInd] << 2) | (curKmerInt[0]>>62);
			curRcInt[1] = ((uint64_t)charToRcMap[(readStr[j] & 6)>>1] <<  ((j-32)<<1)) | curRcInt[1];
		}
		else
		{
			curRcInt[0] = ((uint64_t)charToRcMap[(readStr[j] & 6)>>1] <<  (j<<1)) | curRcInt[0];
		}
		curKmerInt[0] = (curKmerInt[0] << 2) | charToIntMap[(readStr[j] & 6)>>1];
//		printf ("%llx %llx - %c\n ", curRcInt[1], curRcInt[0], readStr[j]);
		
	}
	if (lastNPos == -1)
	{
		if ((curKmerInt[1] < curRcInt[1]) || ((curKmerInt[1] == curRcInt[1]) && (curKmerInt[0] < curRcInt[0])))
		{
//			printf ("%llx %llx ", curKmerInt[1], curKmerInt[0]);
//			getKmerStr64 (curKmerInt[0], curKmerInt[1], NULL, 0);
			uint64_t bucket = getMmerPrefix(curKmerInt[1], curKmerInt[0], mmerSize);
			if ((bucket >= startKmer) && (bucket < endKmer))
			{
//				int process=0;
//				while (bucket >= scanProcessKmerRange[process])
//					process++;
//				process--;
				int process = binSrch (0, numP, scanProcessKmerRange, bucket);
//				assert (binFound == process);
//				printf ("Bucket %lu BinFound %d Actual %d\n", bucket, binFound, process);
				lOfs[process]++;
				memcpy (&(lBuf[bufOfs].kmer[0]), &curKmerInt[0], sizeof(uint64_t));
				memcpy (&(lBuf[bufOfs].kmer[2]), &curKmerInt[1], sizeof(uint64_t));
				lBuf[bufOfs].adj = ((prevChar<<2) & 0xF8) | ((nextChar>>1) & 0x7);
				lBuf[bufOfs].process=process;
				bufOfs++;
			}
		}
		else
		{
//			printf ("%llx %llx ", curRcInt[1], curRcInt[0]);
//			getKmerStr64 (curRcInt[0], curRcInt[1], NULL, 0);
			uint64_t bucket = getMmerPrefix(curRcInt[1], curRcInt[0], mmerSize);
			if ((bucket >= startKmer) && (bucket < endKmer))
			{
//				int process=0;
//				while (bucket >= scanProcessKmerRange[process])
//					process++;
//				process--;
				int process = binSrch (0, numP, scanProcessKmerRange, bucket);
//				assert (binFound == process);
//				printf ("Bucket %lu BinFound %d Actual %d\n", bucket, binFound, process);
				lOfs[process]++;
				memcpy (&(lBuf[bufOfs].kmer[0]), &curRcInt[0], sizeof(uint64_t));
				memcpy (&(lBuf[bufOfs].kmer[2]), &curRcInt[1], sizeof(uint64_t));
				prevChar = "TGAC   N"[(prevChar>>1) & 7];
				nextChar = "TGAC   N"[(nextChar>>1) & 7];
				lBuf[bufOfs].adj = ((nextChar<<2) & 0xF8) | ((prevChar>>1) & 0x7);
				lBuf[bufOfs].process=process;
				bufOfs++;
			}
		}
	}
/*	if (curReadId == 18173149)
	{
		printf ("Aft Stage 1\n");
	}*/
	for (i=1; i <= lastKmerPos; i++)
	{
		prevChar = readStr[i-1];
		if (i==lastKmerPos)
			nextChar='N';
		else
			nextChar = readStr[i+kmerSize];
		if (kmerSize > 32)
		{
			curKmerInt[1]=((curKmerInt[1]<<2) | (curKmerInt[0]>>62)) & mask[1];
			curRcInt[0]=curRcInt[0] >> 2;
			curRcInt[0] |= ((uint64_t)(curRcInt[1] & 3) << 62);
			curRcInt[1] = curRcInt[1]>>2;
			int remKmerSize = kmerSize-32;
			curRcInt[1] |= (((uint64_t)charToRcMap[(readStr[i+kmerSize-1] & 6) >> 1]) << ((remKmerSize-1)<<1));
		}
		else
		{
			curRcInt[0]=curRcInt[0] >> 2;
			curRcInt[0] |= (((uint64_t)charToRcMap[(readStr[i+kmerSize-1] & 6) >> 1]) << ((kmerSize-1)<<1));
		}
		curKmerInt[0] = (curKmerInt[0] << 2) | charToIntMap[(readStr[i+kmerSize-1] & 6) >> 1];
		curKmerInt[1] = curKmerInt[1] & mask[1];
		curKmerInt[0] = curKmerInt[0] & mask[0];
		if (readStr[i+kmerSize-1] == 'N')
			lastNPos = i+kmerSize-1;
		if (i > lastNPos)
		{
			if ((curKmerInt[1] < curRcInt[1]) || ((curKmerInt[1] == curRcInt[1]) && (curKmerInt[0] < curRcInt[0])))
			{
//				printf (" %c - %llx\t%llx %llx\t", readStr[i+kmerSize-1], curKmerInt[1], curKmerInt[0], getMmerPrefix(curKmerInt[1], curKmerInt[0], 10));
//				getKmerStr64 (curKmerInt[0], curKmerInt[1], NULL, 0);
				uint64_t bucket = getMmerPrefix(curKmerInt[1], curKmerInt[0], mmerSize);
				if ((bucket >= startKmer) && (bucket < endKmer))
				{
//					int process=0;
//					while (bucket >= scanProcessKmerRange[process])
//						process++;
//					process--;
					int process = binSrch (0, numP, scanProcessKmerRange, bucket);
//					assert (binFound == process);
//					printf ("Bucket %lu BinFound %d Actual %d\n", bucket, binFound, process);
					lOfs[process]++;
					memcpy (&(lBuf[bufOfs].kmer[0]), &curKmerInt[0], sizeof(uint64_t));
					memcpy (&(lBuf[bufOfs].kmer[2]), &curKmerInt[1], sizeof(uint64_t));
					lBuf[bufOfs].adj = ((prevChar<<2) & 0xF8) | ((nextChar>>1) & 0x7);
					lBuf[bufOfs].process=process;
					bufOfs++;
				}
			}
			else
			{
//				printf ("RC - %llx\t%llx %llx\t", curRcInt[1], curRcInt[0], getMmerPrefix(curRcInt[1], curRcInt[0], 10));
//				getKmerStr64 (curRcInt[0], curRcInt[1], NULL, 0);
				uint64_t bucket = getMmerPrefix(curRcInt[1], curRcInt[0], mmerSize);
				if ((bucket >= startKmer) && (bucket < endKmer))
				{
//					int process=0;
//					while (bucket >= scanProcessKmerRange[process])
//						process++;
//					process--;
					int process = binSrch (0, numP, scanProcessKmerRange, bucket);
//					assert (binFound == process);
//					printf ("Bucket %lu BinFound %d Actual %d\n", bucket, binFound, process);
					lOfs[process]++;
					memcpy (&(lBuf[bufOfs].kmer[0]), &curRcInt[0], sizeof(uint64_t));
					memcpy (&(lBuf[bufOfs].kmer[2]), &curRcInt[1], sizeof(uint64_t));
					prevChar = "TGAC   N"[(prevChar>>1) & 7];
					nextChar = "TGAC   N"[(nextChar>>1) & 7];
//		printf ("prevChar %c, nextChar %c\n", prevChar, nextChar);
					lBuf[bufOfs].adj = ((nextChar<<2) & 0xF8) | ((prevChar>>1) & 0x7);
					lBuf[bufOfs].process=process;
					bufOfs++;
				}
			}
		}
	}
/*	if (curReadId == 18173149)
	{
		printf ("Aft Stage 3 bufOfs %d\n", bufOfs);
		for (i=0; i<bufOfs; i++)
		{
			printf ("%d: %lx %lx %lx %lx - %d\n", i, lBuf[i].kmer[3], lBuf[i].kmer[2], lBuf[i].kmer[1], lBuf[i].kmer[0], lBuf[i].process);
		}
	}*/
	uint16_t *ptrs[1024];
	for (i=0; i<numP; i++)
	{
		if (lOfs[i] > 0)
		{
			uint64_t curOfs = curSendOfs[i*numThr + tid];
			ptrs[i]=kR[i*numThr + tid] + curOfs*7;
//			if (curReadId > 18173140)
//				printf ("curReadId %u curSendOfs %lu kR[%d]=%p, ptr[%d]=%p lOfs[%d]=%d\n", curReadId, curOfs, i*numThr + tid, kR[i*numThr + tid], i, ptrs[i], i, lOfs[i]);
			curSendOfs[i*numThr + tid]+=lOfs[i];
/*			if ((curSendOfs[i*numThr + tid] + readLen) > MAX_TUPLES)
				overflow=1;*/
		}
	}
	if (curReadId == 18173149)
	{
		printf ("Aft Stage 4 ptrs[0] %p\n", ptrs[0]);
	}
	for (i=0; i<bufOfs; i++)
	{
		int p=lBuf[i].process;
		memcpy (&ptrs[p][0], &(lBuf[i].kmer[0]), sizeof(uint32_t));
		memcpy (&ptrs[p][2], &(lBuf[i].kmer[1]), sizeof(uint32_t));
//		ptrs[p][2]=lBuf[i].kmer[2];
//		ptrs[p][3]=lBuf[i].kmer[3];
		uint32_t readRem = (uint32_t)(curReadId%UINT32_MAX);
		ptrs[p][4] = ((uint16_t)(curReadId/UINT32_MAX)) << 8;
		ptrs[p][4] |= (uint16_t)lBuf[i].adj;
//		getKmerStr64 (*(uint64_t *)&ptrs[p][0], *(uint64_t *)&ptrs[p][2], NULL, 0);
		memcpy (&ptrs[p][5], &readRem, sizeof(uint32_t));
		ptrs[p]+=7;
	}
	if (curReadId == 18173149)
	{
		printf ("End\n");
	}
//	return overflow;
}

/*int main()
{
	int i,j;
//	char in[110]="NTGAAAAGTGGGATTCATGATTACTCCTGCCCTGCAGACCCCACACAACTCTTGAGAGGACCAAATTTTAAAAATG";
//	char in[110]="NTGAAAAGTGNGGANTCNATGATTACTCCTGCCCTGCAGACCCCACACAACTCTTGAGNAGGACCAAATTTTAAAA";
	char in[252]="TNTAGGCATTTCTTTACTCCAAAGTGCATTCAGAAAATCAGATGTTCAGCATTTTTTTCACGCGCTTGGTGTCAGACTTGCTGACCAGCGCATCCCCCGCCAGCTTGCGTACGCGTTTGCGGGGTTTTTTTGTCAGGATGTGCCGGGCATGACTGTGATGGCGCACGATTTTCCCTGTCCCCGTCACTTTAAATCTTTTTGCTGCACTCTTTGATGTCTTAAGCTTCACAGGTATTTCCCCCAACGGATTT";
	uint64_t scanProcessKmerRange[101]={0, 2819, 7943, 13152, 18452, 24599, 31566, 37095, 42320, 49242, 55653, 62593, 68879, 77335, 84572, 90888, 96336, 101944, 107528, 113794, 125156, 133236, 141017, 148315, 153846, 158300, 165084, 170473, 177364, 187139, 197406, 208480, 215624, 221301, 226775, 232938, 239074, 247505, 255933, 262811, 270267, 277934, 285430, 293327, 300504, 305058, 314702, 321114, 329704, 336857, 342398, 349537, 353963, 360446, 365579, 368611, 371934, 381634, 387686, 395517, 401962, 407918, 414414, 417632, 421288, 428391, 432531, 436887, 446822, 454620, 475534, 484761, 495129, 501813, 514402, 524289, 533061, 543956, 552891, 565281, 578891, 588359, 599383, 609176, 616034, 624267, 631634, 644254, 655871, 669869, 680367, 692656, 710052, 734775, 761474, 789098, 839624, 868385, 898644, 941056, 1048576};
	int numTasksOut = 1;
	int k=27;
	int readLen = 251;
	uint64_t n=readLen-k+1;
	init(k, numTasksOut, 1, 10);
	uint16_t *kR = (uint16_t *)malloc(numTasksOut* 300 *7* sizeof(uint16_t));
	assert (kR != NULL);
	uint16_t **kRs = (uint16_t **) malloc (numTasksOut * sizeof (uint16_t *));
	assert (kRs != NULL);

	for (i=0; i<numTasksOut; i++)
		kRs[i] = &kR[i*300];

	uint64_t *curSendOfs = (uint64_t *)calloc (numTasksOut, sizeof(uint64_t));
	if (numTasksOut == 1)
		scanProcessKmerRange[1]=1048576;
	uint64_t kmerMask = (((1ul<<k)-1) >> (k-20)) << (32+k-20);
	int rotater = 32+k-20;
	int lastKmerPos, kmersPerLane;
        lastKmerPos = readLen - k;
        kmersPerLane = lastKmerPos/4;
	printf ("lastKmerPos %d, kmersPerLane %d\n", lastKmerPos, kmersPerLane);

	printf ("KmerMask %llx rotateBy %d\n", kmerMask, rotater);
	printf ("Read %s\n", in);
	char *in_bases = &in[0];
	uint64_t largest=0;
	getCanonKmers (in_bases, kRs, curSendOfs, readLen, 0ul, 1048576ul, kmerMask, rotater, 0, &scanProcessKmerRange[0], 0, 0);
        printf ("numKmers %lu\n", curSendOfs[0]);
	uint64_t totalKmers = 0;
	for (j=0; j<numTasksOut; j++)
	{
		printf ("TO %d, Num Kmers %lu\n", j, curSendOfs[j]);
		totalKmers += curSendOfs[j];
		for (i=0; i<curSendOfs[j]; i++)
		{
			printf ("%s\n", in);
			for (k=0; k<i; k++) printf (" ");
			uint8_t adj = (uint8_t)(kRs[j][i*7+4]&0xFF);
			printf ("%c-", "ACTG    "[adj>>3 & 7]);
			getKmerStr64 (*(uint64_t *)&kRs[j][i*7], *(uint64_t *)&kRs[j][i*7], NULL, 0);
			printf ("-%c\n", "ACTG    "[adj & 7]);
		}
	}
	printf ("Total Kmers %lu\n", totalKmers);

	free(kR);
	free (kRs);
	free (curSendOfs);
	return 0;
}*/

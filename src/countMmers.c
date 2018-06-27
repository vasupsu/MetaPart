#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
	
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
void init(int kLen, int mLen)
{
//	readLen=rLen;
	kmerSize=kLen;
	mmerSize=mLen;
//	printf ("Mask32 %x\n", (1<<kmerSize) -1);
}
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

int getCanonKmers (char *readStr, int readLen, uint64_t *rangeHist, uint32_t *kmerFreqCount)
{
	int charToIntMap[4]={0,1,3,2};
	int charToRcMap[4]={3,2,0,1};
	if (readLen <= kmerSize)
		return 0;
	int i=0,j;
	int lastKmerPos, kmersPerLane;
	lastKmerPos = readLen - kmerSize;
	kmersPerLane = lastKmerPos/4;
	uint64_t curKmerInt[2]={0,0};
	uint64_t curRcInt[2]={0,0};
	int lastNPos = -1;
	uint64_t mask[2];
	int rotateBy[2];
	int kmerIndex=0;
	
//	uint8_t lOfs[1024];
//        for (i=0; i<numP; i++) lOfs[i]=0;
//	assert (lastKmerPos < 200);
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
//			if (bucket < 0x1717clu)
			{
				rangeHist[bucket+2]++;
				kmerFreqCount[bucket]++;
				kmerIndex++;
			}
		}
		else
		{
//			printf ("%llx %llx ", curRcInt[1], curRcInt[0]);
//			getKmerStr64 (curRcInt[0], curRcInt[1], NULL, 0);
			uint64_t bucket = getMmerPrefix(curRcInt[1], curRcInt[0], mmerSize);
//			if (bucket < 0x1717clu)
			{
				rangeHist[bucket+2]++;
				kmerFreqCount[bucket]++;
				kmerIndex++;
			}
		}
//		kmerIndex++;
	}
	for (i=1; i <= lastKmerPos; i++)
	{
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
//				if (bucket < 0x1717clu)
				{
					rangeHist[bucket+2]++;
					kmerFreqCount[bucket]++;
					kmerIndex++;
				}
			}
			else
			{
//				printf ("RC - %llx\t%llx %llx\t", curRcInt[1], curRcInt[0], getMmerPrefix(curRcInt[1], curRcInt[0], 10));
//				getKmerStr64 (curRcInt[0], curRcInt[1], NULL, 0);
				uint64_t bucket = getMmerPrefix(curRcInt[1], curRcInt[0], mmerSize);
//				if (bucket < 0x1717clu)
				{
					rangeHist[bucket+2]++;
					kmerFreqCount[bucket]++;
					kmerIndex++;
				}
			}
//			kmerIndex++;
		}
	}
	return kmerIndex;
}

/*int main()
{
	int i,j;
	char in[110]="NTGAAAAGTGGGATTCATGATTACTCCTGCCCTGCAGACCCCACACAACTCTTGAGAGGACCAAATTTTAAAAATG";
//	char in[110]="NTGAAAAGTGNGGANTCNATGATTACTCCTGCCCTGCAGACCCCACACAACTCTTGAGNAGGACCAAATTTTAAAA";
	int k=64;
	init(k, 10);
	uint64_t *rangeHist = (uint64_t *)calloc((1048576+2), sizeof(uint64_t));
	assert (rangeHist != NULL);
	uint32_t *kmerFreqCount = (uint32_t *)calloc(1048576, sizeof(uint32_t));
	assert (kmerFreqCount != NULL);
	printf ("Read %s\n", in);
	char *in_bases = &in[0];
	getCanonKmers (in_bases, 76, rangeHist, kmerFreqCount);
	for (i=0; i<1048576; i++)
	{
		if (kmerFreqCount[i]>0)
			printf ("%x - %u\n", i, kmerFreqCount[i]);
	}
	
	free(rangeHist);
	free(kmerFreqCount);
	return 0;
}*/

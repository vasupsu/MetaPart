#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <sys/stat.h>

typedef struct {
        size_t offset;
        size_t offset2;
	uint8_t readIdExtra;
        uint32_t readId;
        int fileNo;
        uint32_t kmerFreqCount[1048576+2];
}fileIndex;

struct stat buf;

int  main(int argc, char **argv)
{
	if (argc != 2) {
		fprintf(stderr, "%s fqIdxFile\n", argv[0]);
		exit(1);
	}
	int status = stat (argv[1], &buf);
        assert (status == 0);
	int numRecs = buf.st_size/sizeof(fileIndex);
	int i=0,j=0;
	fileIndex fI;
//	assert (fI != NULL);
	FILE *fp = fopen (argv[1], "r");
	fread(&fI, sizeof(fileIndex), 1, fp);
	uint64_t numKmers=0;
	for (i=0; i<numRecs; i++)
	{
		printf ("%d %d %lu,%lu %u-%u\n", i, fI.fileNo, fI.offset, fI.offset2, fI.readIdExtra, fI.readId);
	        uint64_t chunkKmers = 0;
/*        	for (j=0; j<1048576; j++)
	        {
        	    chunkKmers += fI[i].kmerFreqCount[j];
	        }
        	printf ("%lu\n", chunkKmers);
	        numKmers += chunkKmers;*/
		fread(&fI, sizeof(fileIndex), 1, fp);
	}
	fclose (fp);
	printf ("TotalKmers %lu\n", numKmers);
//	free (fI);
	return 0;
}

echo "Compiling code..."
gcc -o bin/genFastqIdx src/genFastqIdx.c
mpicc -O3 -fopenmp -o bin/genHistograms_metapart src/genHistograms.c src/countMmers.c
mpicc -O3 -fopenmp -o bin/genHistograms_metapartmin src/genHistograms.c src/countMmers_3minimizers.c -DMIN3
mpicc -O3 -fopenmp -o bin/metapart src/metapart.c src/enumerateKmers.c
mpicc -O3 -fopenmp -o bin/metapartmin src/metapart.c src/enumerateKmers_3minimizers.c -DMIN3
echo "Compile Done."

echo "Creating FASTQ index file"
bin/genFastqIdx 2 1 0 0 data/a1.fq data/a2.fq
echo "Created 2 chunks"

echo "Creating merHist and FASTQPart tables"
export OMP_NUM_THREADS=1
mpiexec -np 1 bin/genHistograms_metapart 27 1 0 0 data/a1.fq data/a2.fq
mpiexec -np 1 bin/genHistograms_metapartmin 27 1 0 0 data/a1.fq data/a2.fq

echo "Run MetaPrep using k=27"
mpiexec -np 1 bin/metapart -o ./output 30000 27 1 2 1 0 0 data/a1.fq data/a2.fq
#mpiexec -np 1 bin/metapartmin -o ./output 30000 27 1 2 1 0 0 data/a1.fq data/a2.fq
echo "Done. Output FASTQ files present in output directory"

echo "Cleaning index files"
rm data/a1.fq.*

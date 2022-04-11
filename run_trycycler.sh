#!/bin bash


# Install Trycycler
mamba create --name trycycler -c bioconda trycycler any2fasta flye raven-assembler shasta wtdbg2 necat -y



# Activate
conda activate trycycler


# Fastq files
long=/home/bioinfo/analyses/salmonella_nanopore_wildbirds_merged/filtered/32_filtered.fastq.gz
out=/home/bioinfo/analyses/salmonella_nanopore_wildbirds_merged/32
size=4850000


# Logs
# trycycler any2fasta flye raven-assembler shasta wtdbg2




# Trycycler subsample
# Max 8 subsamples
# Since we're using 3 assemblers and there's a limit of 26 total assemblies (8x3 = 24)
trycycler subsample \
    --reads "$long" \
    --out_dir "${out}"/read_subsets \
    --count 8 \
    --genome_size "$size" \
    --threads $(nproc)


# Assemble each subsample with different assemblers
function assemble_shasta(){
    sample=$(basename "$1" ".fastq")

    shasta \
        --input "$1" \
        --config Nanopore-Oct2021 \
        --assemblyDirectory "$2" \
        --Reads.minReadLength 3000 \
        --threads $(nproc)

    shasta \
        --command cleanupBinaryData \
        --assemblyDirectory "$2"

    mv "${2}"/Assembly.fasta "${out}"/assemblies/"${sample}"_shasta.fasta
}

function assemble_flye(){
    sample=$(basename "$1" ".fastq")

    # --nano-hq
    flye \
        --nano-raw "$1" \
        --out-dir "$2" \
        --genome-size "$size" \
        --min-overlap 1000 \
        --threads $(nproc)

    mv "${2}"/assembly.fasta "${out}"/assemblies/"${sample}"_flye.fasta
}

function assemble_raven(){
    sample=$(basename "$1" ".fastq")

    raven \
        --threads $(nproc) \
        "$1" \
        > "${2}"/assembly.fasta

    mv "${2}"/assembly.fasta "${out}"/assemblies/"${sample}"_raven.fasta
    rm raven.cereal
}

function assemble_wtdbg2(){
    sample=$(basename "$1" ".fastq")

    wtdbg2 \
        -x ont \
        -g "$size" \
        -t $(nproc) \
        -i "$1" \
        -fo "${2}"/"$sample"

    wtpoa-cns \
        -t $(nproc) \
        -i "${2}"/"${sample}".ctg.lay.gz \
        -fo "${2}"/"${sample}"_wtdbg2.fasta

    mv "${2}"/"${sample}"_wtdbg2.fasta "${out}"/assemblies/"${sample}"_wtdbg2.fasta
}

function assemble_necat(){
    sample=$(basename "$1" ".fastq")

    # Create ONT read list
    echo "$long" > "${2}"/necat_readlist.txt

    # Create config file
    echo -n "\
PROJECT="$sample"
ONT_READ_LIST="${2}"/necat_readlist.txt
GENOME_SIZE="$size"
THREADS=$(nproc)
MIN_READ_LENGTH=3000
PREP_OUTPUT_COVERAGE=40
OVLP_FAST_OPTIONS=-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000
OVLP_SENSITIVE_OPTIONS=-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000
CNS_FAST_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
CNS_SENSITIVE_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
TRIM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400
ASM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400
NUM_ITER=2
CNS_OUTPUT_COVERAGE=30
CLEANUP=1
USE_GRID=false
GRID_NODE=0
GRID_OPTIONS=
SMALL_MEMORY=0
FSA_OL_FILTER_OPTIONS=
FSA_ASSEMBLE_OPTIONS=
FSA_CTG_BRIDGE_OPTIONS=
POLISH_CONTIGS=true
" > "${2}"/"${sample}"_necat.cfg

    cd "$2"

    necat correct \
        "${2}"/"${sample}"_necat.cfg

    necat assemble \
        "${2}"/"${sample}"_necat.cfg

    necat bridge \
        "${2}"/"${sample}"_necat.cfg

    mv "${2}"/"${sample}"/6-bridge_contigs/polished_contigs.fasta \
        "${out}"/assemblies/"${sample}"_necat.fasta
}


# Run all assemblies
# Max of 26 assemblies
[ -d "${out}"/assemblies ] || mkdir -p "${out}"/assemblies
counter=0
for i in $(find "${out}"/read_subsets -name "*.fastq"); do
    sample=$(basename "$1" ".fastq")
    # counter=$(($counter+1))
    # [ -d "${out}"/assembly"${counter}" ] || mkdir -p "${out}"/assembly"${counter}"
    # assemble_wtdbg2 "$i" "${out}"/assembly"${counter}"
    counter=$(($counter+1))
    [ -d "${out}"/assembly"${counter}" ] || mkdir -p "${out}"/assembly"${counter}"
    assemble_necat "$i" "${out}"/assembly"${counter}"
    counter=$(($counter+1))
    [ -d "${out}"/assembly"${counter}" ] || mkdir -p "${out}"/assembly"${counter}"
    assemble_shasta "$i" "${out}"/assembly"${counter}"
    counter=$(($counter+1))
    [ -d "${out}"/assembly"${counter}" ] || mkdir -p "${out}"/assembly"${counter}"
    assemble_flye "$i" "${out}"/assembly"${counter}"
    # counter=$(($counter+1))
    # [ -d "${out}"/assembly"${counter}" ] || mkdir -p "${out}"/assembly"${counter}"
    # assemble_raven "$i" "${out}"/assembly"${counter}"

done


[ -d "${out}"/logs ] ||  mkdir -p "${out}"/logs


# Trycycler cluster
# Default distance is 0.01 (1%)
trycycler cluster \
    --threads $(nproc) \
    --distance 0.01 \
    --assemblies "${out}"/assemblies/*.fasta \
    --reads "$long" \
    --out_dir "${out}"/clusters \
    2>&1 | tee "${out}"/logs/cluster.log


# Manual inspection of cluster required
# Bad clusters and/or contigs withing the cluster have to be deleted
blastn -task megablast \
    -db nt \
    -query /home/bioinfo/anaylses/salmonella_nanopore_wildbirds_merged/13/clusters/cluster_004/1_contigs/A_contig_2.fasta \
    -num_threads $(nproc) \
    -evalue 1e-10 \
    -max_target_seqs 5 \
    -html \
    -out /home/bioinfo/anaylses/salmonella_nanopore_wildbirds_merged/13/clusters/cluster_004/1_contigs/A_contig_2.fasta.blastn.html



# Save the good clusters into an array
declare -a cluster_array=(001 002)


# Trycycler reconcile
for c in "${cluster_array[@]}"; do
    trycycler reconcile \
        --threads $(nproc) \
        --reads "$long" \
        --cluster_dir "${out}"/clusters/cluster_"$c" \
        2>&1 | tee "${out}"/logs/cluster_"$c".reconcile.log
done


# Trycycler MSA
for c in "${cluster_array[@]}"; do
    trycycler msa \
        --threads $(nproc) \
        --cluster_dir "${out}"/clusters/cluster_"$c" \
        2>&1 | tee "${out}"/logs/cluster_"$c".msa.log
done


# Trycycler partition
trycycler partition \
    --threads $(nproc) \
    --reads "$long" \
    --cluster_dirs "${out}"/clusters/cluster_* \
    2>&1 | tee "${out}"/logs/cluster_"$c".partition.log


# Trycycler consensus
for c in "${cluster_array[@]}"; do
    trycycler consensus \
        --threads $(nproc) \
        --cluster_dir "${out}"/clusters/cluster_"$c" \
        2>&1 | tee "${out}"/logs/cluster_"$c".consensus.log
done


# Combine all consensus
cat "${out}"/clusters/cluster_*/7_final_consensus.fasta \
    > "${out}"/clusters/consensus.fasta

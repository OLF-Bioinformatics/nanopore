#!/bin bash


# Install Trycycler
mamba create --name trycycler -c bioconda -c conda-forge -c defaults \
    trycycler any2fasta flye raven-assembler shasta wtdbg necat canu \
    minipolish smartdenovo \
    -y


# Activate
conda activate trycycler


# Fastq files
long=/home/bioinfo/analyses/salmonella_nanopore_wildbirds_merged/filtered/37_filtered.fastq.gz
out=/home/bioinfo/analyses/salmonella_nanopore_wildbirds_merged/37
genome_size=4850000

# Read
min_read_size=5000
min_read_overlap=3000

# Performance
cpu=$(nproc)
mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB



# Logs
# trycycler any2fasta flye raven-assembler shasta wtdbg2




# Trycycler subsample
# Max 26 assemblies (A-Z)
# Since we're using 4 assemblers, we can only use 6 subsamples (5x5 = 26)
trycycler subsample \
    --reads "$long" \
    --out_dir "${out}"/read_subsets \
    --count 8 \
    --genome_genome_size "$genome_size" \
    --threads "$cpu"


# Assemble each subsample with different assemblers
function assemble_shasta(){
    sample=$(basename "$1" ".fastq")

    shasta \
        --input "$1" \
        --config Nanopore-Oct2021 \
        --assemblyDirectory "$2" \
        --Reads.minReadLength "$min_read_size" \
        --threads "$cpu"

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
        --genome-genome_size "$genome_size" \
        --min-overlap "$min_read_overlap" \
        --threads "$cpu"

    mv "${2}"/assembly.fasta "${out}"/assemblies/"${sample}"_flye.fasta
}

function assemble_raven(){
    sample=$(basename "$1" ".fastq")

    raven \
        --graphical-fragment-assembly "${2}"/assembly.gfa \
        --polishing-rounds 5 \
        --threads "$cpu" \
        "$1" \
        > "${2}"/assembly.fasta

    mv "${2}"/assembly.fasta "${out}"/assemblies/"${sample}"_raven.fasta
    rm raven.cereal
}

function assemble_wtdbg2(){
    sample=$(basename "$1" ".fastq")

    wtdbg2 \
        -x ont \
        -g "$genome_size" \
        -t "$cpu" \
        -i "$1" \
        -fo "${2}"/"$sample"

    wtpoa-cns \
        -t "$cpu" \
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
GENOME_SIZE="$genome_size"
THREADS="$cpu"
MIN_READ_LENGTH="$min_read_size"
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

function assemble_canu(){
    sample=$(basename "$1" ".fastq")

    canu \
        -p "$sample" \
        -d "$2" \
        genomeSize="$genome_size" \
        maxInputCoverage=50 \
        minReadLength="$min_read_size" \
        minThreads="$cpu" \
        -nanopore "$1"

    mv "${2}"/"${sample}".contigs.fasta \
        "${out}"/assemblies/"${sample}"_canu.fasta

}

function assemble_miniasm(){
    sample=$(basename "$1" ".fastq")

    minimap2 \
        -t "$cpu" \
        -x ava-ont \
        "$1" \
        "$1" \
        > "${2}"/"${sample}"_overlaps.paf

    miniasm \
        -f "$1" \
        "${2}"/"${sample}"_overlaps.paf \
        > "${2}"/"${sample}"_assembly.gfa

    minipolish \
        -t "$cpu" \
        "$1" \
        "${2}"/"${sample}"_assembly.gfa \
        > "${2}"/"${sample}"_polished.gfa

    any2fasta \
        "${2}"/"${sample}"_polished.gfa \
        > "${2}"/"${sample}".fasta

    mv "${2}"/"${sample}".fasta \
        "${out}"/assemblies/"${sample}"_miniasm.fasta
}


function assemble_nextdenovo(){
    sample=$(basename "$1" ".fastq")

    # Create ONT read list
    echo "$long" > "${2}"/nextdenovo_readlist.txt

    # Prep config file
    echo -n "\
[General]
job_type = local
job_prefix = nextDenovo
task = all # 'all', 'correct', 'assemble'
rewrite = yes # yes/no
deltmp = yes
rerun = 3
parallel_jobs = $(($mem/$cpu))
input_type = raw
read_type = "ont"
input_fofn = "${2}"/nextdenovo_readlist.txt
workdir = "$2"

[correct_option]
read_cutoff = "$min_read_size"
genome_genome_size = "$genome_size"
pa_correction = 2
sort_options = -m 1g -t 2
minimap2_options_raw =  -t 8
correction_options = -p 15

[assemble_option]
minimap2_options_cns =  -t $(("$cpu"/$(($mem/$cpu))))
nextgraph_options = -a 1
" > "${2}"/"${sample}"_nextdenovo.cfg

python ~/prog/NextDenovo/nextDenovo "${2}"/"${sample}"_nextdenovo.cfg

mv "${2}"/03.ctg_graph/nd.asm.fasta \
    "${out}"/assemblies/"${sample}"_nextdenovo.fasta
    
}


function assemble_smartdenovo(){
    sample=$(basename "$1" ".fastq")

    smartdenovo.pl \
        -p "${2}"/"${sample}" \
        -c 1 \
        -t "$cpu" \
        -J "$min_read_size" \
        "$1" \
        > "${2}"/"${sample}".mak
    
    make -f "${2}"/"${sample}".mak


    mv "${2}"/"${sample}".dmo.cns \
        "${out}"/assemblies/"${sample}"_smartdenovo.fasta


}



# Run all assemblies
# Max of 26 assemblies
[ -d "${out}"/assemblies ] || mkdir -p "${out}"/assemblies
job_id=0
for i in $(find "${out}"/read_subsets -name "*.fastq"); do
    sample=$(basename "$1" ".fastq")
    
    # job_id=$(($job_id+1))
    # [ -d "${out}"/assembly"${job_id}" ] || mkdir -p "${out}"/assembly"${job_id}"
    # assemble_wtdbg2 "$i" "${out}"/assembly"${job_id}"
    
    # job_id=$(($job_id+1))
    # [ -d "${out}"/assembly"${job_id}" ] || mkdir -p "${out}"/assembly"${job_id}"
    # assemble_necat "$i" "${out}"/assembly"${job_id}"

    # job_id=$(($job_id+1))
    # [ -d "${out}"/assembly"${job_id}" ] || mkdir -p "${out}"/assembly"${job_id}"
    # assemble_miniasm "$i" "${out}"/assembly"${job_id}"

    # job_id=$(($job_id+1))
    # [ -d "${out}"/assembly"${job_id}" ] || mkdir -p "${out}"/assembly"${job_id}"
    # assemble_canu "$i" "${out}"/assembly"${job_id}"

    # job_id=$(($job_id+1))
    # [ -d "${out}"/assembly"${job_id}" ] || mkdir -p "${out}"/assembly"${job_id}"
    # assemble_nextdenovo "$i" "${out}"/assembly"${job_id}"

    # job_id=$(($job_id+1))
    # [ -d "${out}"/assembly"${job_id}" ] || mkdir -p "${out}"/assembly"${job_id}"
    # assemble_smartdenovo "$i" "${out}"/assembly"${job_id}"

    job_id=$(($job_id+1))
    [ -d "${out}"/assembly"${job_id}" ] || mkdir -p "${out}"/assembly"${job_id}"
    assemble_raven "$i" "${out}"/assembly"${job_id}"

    job_id=$(($job_id+1))
    [ -d "${out}"/assembly"${job_id}" ] || mkdir -p "${out}"/assembly"${job_id}"
    assemble_shasta "$i" "${out}"/assembly"${job_id}"
    
    job_id=$(($job_id+1))
    [ -d "${out}"/assembly"${job_id}" ] || mkdir -p "${out}"/assembly"${job_id}"
    assemble_flye "$i" "${out}"/assembly"${job_id}"
done


[ -d "${out}"/logs ] ||  mkdir -p "${out}"/logs


# Trycycler cluster
# Default distance is 0.01 (1%)
trycycler cluster \
    --threads "$cpu" \
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
    -num_threads "$cpu" \
    -evalue 1e-10 \
    -max_target_seqs 5 \
    -html \
    -out /home/bioinfo/anaylses/salmonella_nanopore_wildbirds_merged/13/clusters/cluster_004/1_contigs/A_contig_2.fasta.blastn.html



# Save the good clusters into an array
declare -a cluster_array=(001 002 003)


# Trycycler reconcile
for c in "${cluster_array[@]}"; do
    trycycler reconcile \
        --threads "$cpu" \
        --reads "$long" \
        --cluster_dir "${out}"/clusters/cluster_"$c" \
        2>&1 | tee "${out}"/logs/cluster_"$c".reconcile.log
done


# Trycycler MSA
for c in "${cluster_array[@]}"; do
    trycycler msa \
        --threads "$cpu" \
        --cluster_dir "${out}"/clusters/cluster_"$c" \
        2>&1 | tee "${out}"/logs/cluster_"$c".msa.log
done


# Trycycler partition
trycycler partition \
    --threads "$cpu" \
    --reads "$long" \
    --cluster_dirs "${out}"/clusters/cluster_* \
    2>&1 | tee "${out}"/logs/cluster_"$c".partition.log


# Trycycler consensus
for c in "${cluster_array[@]}"; do
    trycycler consensus \
        --threads "$cpu" \
        --cluster_dir "${out}"/clusters/cluster_"$c" \
        2>&1 | tee "${out}"/logs/cluster_"$c".consensus.log
done


# Combine all consensus
cat "${out}"/clusters/cluster_*/7_final_consensus.fasta \
    > "${out}"/clusters/consensus.fasta

#!/bin bash


# Activate conda environment
conda activate nanopore


###############
#
#   Installation, once
#
###############


# Install fastp
mamba install -c bioconda fastp 


# Install NextPolish
sudo apt-get install python-dev libbz2-dev liblzma-dev
mamba install -c bioconda seqtk -y
cd $HOME/prog
wget https://github.com/Nextomics/NextPolish/releases/download/v1.4.0/NextPolish.tgz
tar zxvf NextPolish.tgz && rm NextPolish.tgz
cd NextPolish && make -j


# Install polypolish
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh  # Install rust
cd $HOME/prog
git clone https://github.com/rrwick/Polypolish.git
cd Polypolish
cargo build --release

# Install ntEdit
mamba install -c bioocnda ntedit -y


###############
#
#   Data
#
###############


# Folders
export fastq=/media/36tb/data/hongsheng_listeria/illumina/Listeria_Amit
export trimmed_illumina=/home/bioinfo/anaylses/GTA_listeria/trimmed_illumina
export genome_folder=/home/bioinfo/anaylses/GTA_listeria/polished
export polished_illumina=/home/bioinfo/anaylses/GTA_listeria/polished_illumina

# Create output folder
[ -d "$trimmed_illumina" ] || mkdir -p "$trimmed_illumina"
[ -d "$polished_illumina" ] || mkdir -p "$polished_illumina"


###############
#
#   Resources
#
###############


# Performance
export cpu=$(nproc)
export maxProc=1


###############
#
#   Trim
#
###############


# Pre-process short reads
function trim()
{
    sample=$(basename "$1" | cut -d "_" -f 1)
    r1="$1"
    r2=$(echo "$r1" | sed 's/_R1/_R2/')

    fastp \
        --in1 "$r1" \
        --in2 "$r2" \
        --out1 "${trimmed_illumina}"/"${sample}"_trimmed_R1.fastq.gz \
        --out2 "${trimmed_illumina}"/"${sample}"_trimmed_R2.fastq.gz \
        --detect_adapter_for_pe \
        --cut_right \
        --cut_right_mean_quality 10 \
        --length_required 64 \
        --html "${trimmed_illumina}"/"${sample}".html \
        --thread $((cpu/maxProc))
}

export -f trim

find "$fastq" -name "*_R1*.fastq.gz" | \
parallel    --bar \
            --env trim \
            --env cpu \
            --env trimmed_illumina \
            --env maxProc \
            --jobs $maxProc \
            trim {}


###############
#
#   Polish
#
###############


function run_nextpolish()
{
    genome="$1"
    r1="$2"
    r2="$3"
    out="$4"

    sample=$(basename "$r1" | cut -d "_" -f 1)


    # Index genome for bwa
    bwa index "$genome"

    # Map reads
    bwa mem -t $((cpu/maxProc)) "$genome" "$r1" "$r2" | \
    samtools view -@ $((cpu/maxProc)) -F 0x4 -b - | \
    samtools fixmate -@ $((cpu/maxProc)) -m - - | \
    samtools sort -@ $((cpu/maxProc)) -m 4g - | \
    samtools markdup -@ $((cpu/maxProc)) -r - "${out}"/"${sample}".bam

    # Index bam file
    samtools index -@ $((cpu/maxProc)) "${out}"/"${sample}".bam
    
    # Index genome for NextPolish
    samtools faidx "$genome"

    # Polish
    python $HOME/prog/NextPolish/lib/nextpolish1.py \
        --genome "$genome" \
        --task 1 \
        --process $((cpu/maxProc)) \
        --bam_sgs "${out}"/"${sample}".bam \
        -ploidy 1 \
        > "${out}"/"${sample}".nextpolish.fasta

    # Make uppercase and one line per sequence
    seqtk seq "${out}"/"${sample}".nextpolish.fasta | \
        awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' \
        > "${out}"/"${sample}".fasta

    # Cleanup
    rm "${genome}".amb \
        "${genome}".ann \
        "${genome}".bwt \
        "${genome}".pac \
        "${genome}".sa \
        "${genome}".fai \
        "${out}"/"${sample}".bam* \
        "${out}"/"${sample}".nextpolish.fasta

    # Update genome for next step
    genome="${out}"/"${sample}".fasta
}


function run_ntedit()
{
    genome="$1"
    r1="$2"
    r2="$3"
    out="$4"

    sample=$(basename "$r1" | cut -d "_" -f 1)


    # Polish
    nthits \
        -t $((cpu/maxProc)) \
        -k 40 \
        -p "${out}"/"${sample}" \
        --outbloom \
        --solid \
        "$r1" "$r2"

    ntedit \
        -t $((cpu/maxProc)) \
        -f "$genome" \
        -b "${out}"/"${sample}" \
        -r "${out}"/"${sample}"_k40.bf \
        -m 1

    # Make uppercase and one line per sequence.
    seqtk seq "${out}"/"${sample}"_edited.fa | \
        awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' \
        > "${out}"/"${sample}".fasta

    # Cleanup
    rm "${out}"/"${sample}"_variants.vcf \
        "${out}"/"${sample}"_changes.tsv \
        "${out}"/"${sample}"_k40.bf \
        "${out}"/"${sample}"_edited.fa

    # Update genome for next step
    genome="${out}"/"${sample}".fasta
}


function run_polypolish()
{
    genome="$1"
    r1="$2"
    r2="$3"
    out="$4"

    sample=$(basename "$r1" | cut -d "_" -f 1)


    # Index genome
    bwa index "$genome"

    # Map reads
    bwa mem -t 64 -a "$genome" "$r1" > "${out}"/"${sample}"_R1.sam
    bwa mem -t 64 -a "$genome" "$r2" > "${out}"/"${sample}"_R2.sam

    # Polish
    polypolish \
        "$genome" \
        "${out}"/"${sample}"*.sam \
        > "${out}"/"${sample}".polypolish.fasta

    # Make uppercase and one line per sequence.
    seqtk seq "${out}"/"${sample}".polypolish.fasta | \
        awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' \
        > "${out}"/"${sample}".fasta

    # Cleanup
    rm "${genome}".amb \
        "${genome}".ann \
        "${genome}".bwt \
        "${genome}".pac \
        "${genome}".sa \
        "${out}"/"${sample}"*.sam \
        "${out}"/"${sample}".polypolish.fasta

    # Update genome for next step
    genome="${out}"/"${sample}".fasta
}


function rename_fasta_headers()
{
    genome="$1"

    cat "$genome" | \
    awk 'BEGIN {a=0} m=/^>/ {a++;print ">contig_"a} !m {print}' \
    > "${genome}".tmp

    mv "${genome}".tmp "$genome"

}


function polish()
{
    sample=$(basename "$1" | cut -d "_" -f 1)

    genome="$1"
    r1="${trimmed_illumina}"/"${sample}"_trimmed_R1.fastq.gz
    r2="${trimmed_illumina}"/"${sample}"_trimmed_R2.fastq.gz \
    out="$polished_illumina"


    ### NextPolish ###
    # 2 rounds
    for i in {1..2}; do
        run_nextpolish "$genome" "$r1" "$r2" "$out"
    done

    ### ntEdit ###
    run_ntedit "$genome" "$r1" "$r2" "$out"
    

    ### Polypolish ###
    # 3 rounds
    for i in {1..3}; do
        run_polypolish "$genome" "$r1" "$r2" "$out"
    done

    # Fix header
    rename_fasta_headers "$genome"

    # Get avereage depth of coverage
    short_read_coverage "$genome" "$r1" "$r2" "$out"
}


function short_read_coverage()
{
    genome="$1"
    r1="${trimmed_illumina}"/"${sample}"_trimmed_R1.fastq.gz
    r2="${trimmed_illumina}"/"${sample}"_trimmed_R2.fastq.gz \
    out="$polished_illumina"

    sample=$(basename "$r1" | cut -d "_" -f 1)

    bwa index "$genome"

    bwa mem -t $((cpu/maxProc)) "$genome" "$r1" "$r2" | \
    samtools view -@ $((cpu/maxProc)) -F 0x4 -b - | \
    samtools fixmate -@ $((cpu/maxProc)) -m - - | \
    samtools sort -@ $((cpu/maxProc)) -m 4g - | \
    samtools markdup -@ $((cpu/maxProc)) -r - "${out}"/"${sample}".bam

    cov=$(samtools depth "${out}"/"${sample}".bam | awk '{sum+=$3} END {print sum/NR}')
    printf "%s\t%.*f\n" "$sample" 0 "$cov" | tee -a "${out}"/short_reads_coverage.tsv

    # Cleanup
    rm "${genome}".amb \
        "${genome}".ann \
        "${genome}".bwt \
        "${genome}".pac \
        "${genome}".sa \
        "${out}"/"${sample}"*.bam
}




# Polish
for i in $(find "$genome_folder" -type f -name "*.fasta"); do
    polish "$i"
done

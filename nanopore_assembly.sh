#!/bin/bash



# Process and assemble Oxford nanopore sequencer data in fast5 format
# Include file format conversion, stats, qc and de novo assembly



#script version
version="0.2.0"


######################
#                    #
#    User Defined    #
#                    #
######################


# Analysis folder
export baseDir=""${HOME}"/analyses/YpD1_nanopore"

# Reads
export fast5="/media/6tb_raid10/data/noriko/20170720_2026_YpD1_20JUL2017_NG/fast5"

# Database to use for metagomic analysis of raw data (contamination)
db="/media/6tb_raid10/db/centrifuge/2017-10-12_bact_vir_h"
# db="/media/6tb_raid10/db/centrifuge/nt"

# Program location
export prog=""${HOME}"/prog"
export scripts=""${HOME}"/scripts"

# Maximum number of cores used per sample for parallel processing
# A highier value reduces the memory footprint.
export maxProc=8

# Assembly name
export prefix="YpD1"

# Estimated genome size in bp
export size=5000000

# Set smallest contig size for assemblies
export smallest_contig=1000


######################
#                    #
#     Resources      #
#                    #
######################


# Computer performance
export cpu=$(nproc) #total number of cores
export mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
memJava="-Xmx"$mem"g"


#######################
#                     #
#   Data Stucture     #
#                     #
#######################


# Folder structure
export logs=""${baseDir}"/logs"
export qc=""${baseDir}"/qc"
export fastq=""${baseDir}"/fastq"
basecalled=""${baseDir}"/basecalled"
export assemblies=""${baseDir}"/assemblies"
export polished=""${baseDir}"/polished"
aligned=""${baseDir}"/aligned"


# Create folders if do not exist
# "||" if test is false
# "&&" if test is true
[ -d "$baseDir" ] || mkdir -p "$baseDir"
[ -d "$logs" ] || mkdir -p "$logs"
[ -d "$qc" ] || mkdir -p "$qc"
[ -d "$fastq" ] || mkdir -p "$fastq"
[ -d "$assemblies" ] || mkdir -p "$assemblies"
[ -d "$polished" ] || mkdir -p "$polished"
[ -d "$aligned" ] || mkdir -p "$aligned"


################
#              #
#     Log      #
#              #
################


# Date
echo -e "$(date)\n" | tee "${logs}"/log.txt
echo -e "User: $(whoami)" | tee -a "${logs}"/log.txt
echo -e "Processors: "$cpu"" | tee -a "${logs}"/log.txt
echo -e "Memory: "$mem"G" | tee -a "${logs}"/log.txt

#script version
echo -e "\nnanopore_assembly.sh version "$version"\n" | tee -a "${logs}"/log.txt  # $0


####################
#                  #
#   Dependencies   #
#                  #
####################


# Check for dependencies and log versions

# java
if hash java 2>/dev/null; then 
    java -version 2>&1 1>/dev/null | grep "java version" | tr -d '"' | tee -a "${logs}"/log.txt
else
    echo >&2 "java was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# poretools
if hash nanoplot 2>/dev/null; then  # if installed
    poretools -v 2>&1 1>/dev/null | tee -a "${logs}"/log.txt
else
    echo >&2 "poretools was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# FastQC
if hash fastqc 2>/dev/null; then 
    fastqc -v | tee -a "${logs}"/log.txt
else
    echo >&2 "fastQC was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Centrifuge
if hash centrifuge 2>/dev/null; then 
    v=$(centrifuge --version | grep -F "version")
    echo "centrifuge $v" | tee -a "${logs}"/log.txt
else
    echo >&2 "centrifuge was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi


# minimap2
if hash minimap2 2>/dev/null; then
    v=$(minimap2 --version 2>&1)
    echo "minimap2 "$v"" | tee -a "${logs}"/log.txt
else
    echo >&2 "minimap2 was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Samtools
if hash samtools 2>/dev/null; then
    samtools --version | grep -F 'samtools' | tee -a "${logs}"/log.txt
else
    echo >&2 "samtools was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Nanopolish
if hash nanopolish 2>/dev/null; then
    version=$(nanopolish --version | head -n 1)
    echo "$version" | tee -a "${logs}"/log.txt  # not version option with this software, yet
else
    echo >&2 "nanopolish was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi


#check Centrifuge datase
if [ -s "${db}.1.cf" ]; then
    echo -e "\nCentrifuge database: $(basename "$db")" | tee -a "${logs}"/log.txt
else
    echo ""
    exit 1
fi

# Check if reads folder has fast5 files
if [ $(find "$fast5" -type f -name "*.fast5" | wc -l) -eq 0 ]; then
    echo "No .fast5 files are present in the provided data folder"
    exit 1
fi



###################
#                 #
#   Basecalling   #
#                 #
###################


#Use albacore


##########
#        #
#   QC   #
#        #
##########


# NanoPlot



################
#              #
#   Trimming   #
#              #
################


# Trim adapters
porechop \
    -i "${fastq}"/"${prefix}".fastq.gz \
    -o "${fastq}"/"${prefix}"_t.fastq.gz \
    --threads "$cpu" \
    | tee -a "${logs}"/"${prefix}"_porechop.log

# Discard reads < 1000bp  ##############NOT
# Diecard reads with min average quality below 12
bbduk.sh "$memJava" \
    threads="$cpu" \
    in="${fastq}"/"${prefix}"_t.fastq.gz \
    minavgquality=12 \
    minlength=1000 \
    out="${fastq}"/"${prefix}"_trimmed.fastq.gz \
    ziplevel=9 \
    2> >(tee -a "${logs}"/"${prefix}"_bbduk.txt)


####################
#                  #
#   FastQc - Raw   #
#                  #
####################


# Create folder to store report
[ -d "${qc}"/fastqc/raw ] || mkdir -p "${qc}"/fastqc/raw
[ -d "${qc}"/fastqc/trimmed ] || mkdir -p "${qc}"/fastqc/trimmed

# Run fastQC
fastqc \
    --o  "${qc}"/fastqc/raw \
    --noextract \
    --threads "$cpu" \
    "${fastq}"/"${prefix}".fastq.gz

fastqc \
    --o  "${qc}"/fastqc/trimmed \
    --noextract \
    --threads "$cpu" \
    "${fastq}"/"${prefix}"_trimmed.fastq.gz


######################
#                    #
#   Centrifuge raw   #
#                    #
######################


# Create folder
[ -d "${qc}"/centrifuge ] || mkdir -p "${qc}"/centrifuge

# Run centrifuge
centrifuge \
    -p "$cpu" \
    -t \
    --seed "$RANDOM" \
    -x "$db" \
    -U "${fastq}"/"${prefix}"_trimmed.fastq.gz \
    --report-file "${qc}"/centrifuge/"${prefix}"_report.tsv \
    > "${qc}"/centrifuge/"${prefix}".tsv

# Prepare result for display with Krona
cat "${qc}"/centrifuge/"${prefix}".tsv | \
    cut -f 1,3 | \
    ktImportTaxonomy /dev/stdin -o  "${qc}"/centrifuge/"${prefix}".html

# Visualize the resutls in Firefow browser
firefox file://"${qc}"/centrifuge/"${prefix}".html &


########################
#                      #
#   de novo assembly   #
#                      #
########################


# Assemble with miniasm and polish with racon
for i in $(find "$fastq" -type f -name "*.fastq.gz"); do
    sample=$(cut -d "_" -f 1 <<< $(basename "$i"))
    
    # canu \
    #     -p "$prefix" \
    #     -d "$assemblies" \
    #     genomeSize="$size" \
    #     -nanopore-raw "${fastq}"/"${prefix}"_trimmed.fastq.gz

    unicycler \
        -l "$i" \
        -o $"{assemblies}"/unicycler/"$sample" \
        -t "$cpu" \
        --verbosity 2 \
        --mode normal

    mv "${assemblies}"/unicycler/"${sample}"/assembly.fasta \
        "${assemblies}"/unicycler/"${sample}"/"${sample}".fasta

    circlator fixstart \
        --verbose \
        "${assemblies}"/unicycler/"${sample}"/"${sample}".fasta \
        "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered

    # Adjust max line width of fasta file
    perl "${scripts}"/formatFasta.pl \
        -i "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta \
        -o "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta.tmp \
        -w 80

    # Rename header
    # sed -i '/^>/s/$/_unicycler/' "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta.tmp
    
    mv "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta.tmp \
        "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta
done


#################
#               #
#   Polishing   #
#               #
#################


# Polish

function polish()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    [ -d "${polished}"/"$sample" ] || mkdir -p "${polished}"/"$sample"

    # Convert fastq to fasta
    zcat "$1" \
        | sed -n '1~4s/^@/>/p;2~4p' \
        > "${polished}"/"${sample}"/"${sample}"_reads.fasta

    # Index draft genome
    bwa index \
        "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta

    minimap2 \
        -ax map-ont \
        -t $((cpu/maxProc)) \
        "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta \
        "${polished}"/"${sample}"/"${sample}"_reads.fasta | \
    samtools sort -@ $((cpu/maxProc)) -o "${polished}"/"${sample}"/"${sample}".bam -

    # Index bam file
    samtools index -@ $((cpu/maxProc)) "${polished}"/"${sample}"/"${sample}".bam

    # Index reads for nanopolish
    nanopolish index \
        -d /media/6tb_raid10/data/salmonella_human_birds/nanopore/20180503_1623_Wild_Bird_8samples_2018-05-03/fast5 \
        -s "${basecalled}"/sequencing_summary.txt \
        "${polished}"/"${sample}"/"${sample}"_reads.fasta

    # Run nanopolish
    python "${prog}"/nanopolish/scripts/nanopolish_makerange.py \
        "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta \
        | parallel --results "${polished}"/"${sample}"/nanopolish.results -j $((cpu/maxProc)) \
            nanopolish variants \
                --consensus "${polished}"/"${sample}"/"${sample}"_nanopolished.{1}.fa \
                -w {1} \
                -r "${polished}"/"${sample}"/"${sample}"_reads.fasta \
                -b "${polished}"/"${sample}"/"${sample}".bam \
                -g "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta \
                -t 1 \
                -q dcm,dam \
                --fix-homopolymers

    # Merge individual segments
    python  "${prog}"/nanopolish/scripts/nanopolish_merge.py \
        "${polished}"/"${sample}"/"${sample}"_nanopolished.*.fa \
        > "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta

    # Adjust max line width of fasta file
    perl "${scripts}"/formatFasta.pl \
        -i "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta \
        -o "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta.tmp \
        -w 80

    # Rename header
    # sed -i '/^>/s/$/_nanopolish/' "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta.tmp
    
    mv "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta.tmp \
        "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta

    # Cleanup all but fasta read file and polished assembly
    find "${polished}"/"$sample" -type f ! -name "*_nanopolished.fasta" ! -name "*_reads.fasta" -exec rm -rf {} \;
    rm -rf "${polished}"/"${sample}"/nanopolish.results
}

export -f polish

find "$fastq" -type f -name "*.fastq.gz" | \
parallel    --bar \
            --env polish \
            --env qc \
            --env polished \
            --env assemblies \
            --env prog \
            --env cpu \
            --env maxProc \
            --jobs "$maxProc"  \
            "polish {}"


# Compare pre- and post-nanopolished genomes
function compare_assemblies()
{
    export sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    [ -d "${qc}"/mummer/"$sample" ] || mkdir -p "${qc}"/mummer/"$sample"

    # mummer -mum \
    #     -l 20 \
    #     -b \
    #     -n \
    #     "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta \
    #     "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta \
    #     > "${qc}"/mummer/"${sample}"/"${sample}".mum

    contigs=($(cat "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta | grep -E "^>" | cut -d " " -f 1 | tr -d ">"))

    # for i in "${contigs[@]}"; do  # for every fasta entry
    #     mummerplot \
    #         -p "${qc}"/mummer/"${sample}"/"${sample}"_"${i}" \
    #         -title 'Pre- versus post-nanopolish' \
    #         -r "$i" -q "$i" \
    #         --png \
    #         "${qc}"/mummer/"${sample}"/"${sample}".mum
    # done

    nucmer --mum \
        --delta "${qc}"/mummer/"${sample}"/"${sample}".mum \
        "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta \
        "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta

    mummerplot \
        -p "${qc}"/mummer/"${sample}"/"${sample}" \
        -title 'Pre- versus post-nanopolish' \
        -r 1 -q 1 \
        -layout \
        -large \
        --png \
        "${qc}"/mummer/"${sample}"/"${sample}".mum

}

export -f compare_assemblies

find "$polished" -type f -name "*_nanopolished.fasta" | \
parallel    --bar \
            --env compare_assemblies \
            --env polished \
            --env assemblies \
            --env qc \
            --env cpu \
            --env maxProc \
            --jobs "$maxProc" \
            "compare_assemblies {}"


### Coverage ###
function get_coverage()  # unsing unmerged reads only
{
    sample=$(basename "$1" | cut -d '_' -f 1)

    [ -d "${qc}"/coverage/"$sample" ] || mkdir -p "${qc}"/coverage/"$sample"

    minimap2 \
        -ax map-ont \
        -t $((cpu/maxProc)) \
        "$1" \
        "${polished}"/"${sample}"/"${sample}"_reads.fasta | \
    samtools view -@ $((cpu/maxProc)) -b -h -F 4 - | \
    samtools sort -@ $((cpu/maxProc)) - | \
    samtools rmdup - "${qc}"/coverage/"${sample}"/"${sample}".bam

    # Index bam file
    samtools index -@ $((cpu/maxProc)) "${qc}"/coverage/"${sample}"/"${sample}".bam

    #Average genome depth of coverage
    average_cov=$(samtools depth \
        "${qc}"/coverage/"${sample}"/"${sample}".bam  \
        | awk '{sum+=$3} END { print sum/NR}')

    printf "%s\t%.*f\n" "$sample" 0 "$average_cov" | tee -a "${qc}"/coverage/average_cov.tsv

    #Remove reads in fasta format
    rm "${polished}"/"${sample}"/"${sample}"_reads.fasta
}

export -f get_coverage

[ -d "${qc}"/coverage ] || mkdir -p "${qc}"/coverage
echo -e "Sample\tAverage_Cov" > "${qc}"/coverage/average_cov.tsv

find "$polished" -type f -name "*_nanopolished.fasta" | \
    parallel    --bar \
                --env get_coverage \
                --env cpu \
                --env maxProc \
                --env qc \
                --jobs "$maxProc"  \
                "get_coverage {}"

# Qualimap
function run_qualimap()
{
    sample=$(cut -d "." -f 1 <<< $(basename "$1"))
    
    [ -d "${qc}"/coverage/qualimap/"$sample" ] || mkdir -p "${qc}"/coverage/qualimap/"$sample"

    qualimap bamqc \
        -bam "$1" \
        --java-mem-size="${mem}"g \
        -nt $((cpu/maxProc)) \
        -outdir "${qc}"/coverage/qualimap/"$sample" \
        -outfile "${sample}" \
        -outformat HTML \
        -ip

    # Remove bam files
    rm -rf "${qc}"/coverage/"$sample"
}

export -f run_qualimap

find "${qc}"/coverage -type f -name "*.bam" |
parallel    --bar \
            --env run_qualimap \
            --env qc \
            --env mem \
            --env cpu \
            --env maxProc \
            --jobs "$maxProc"
            "run_qualimap {}"

# Blast genomes on nr

function blast()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "${1%.fasta}"))

    [ -d "${qc}"/blast/"$sample" ] || mkdir -p "${qc}"/blast/"$sample"

    blastn \
        -task megablast \
        -db nt \
        -query "$1" \
        -out "${qc}"/blast/"${sample}"/"${sample}".blastn.all.tsv \
        -evalue "1e-10" \
        -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames sskingdoms staxids' \
        -num_threads $((cpu/maxProc)) \
        -culling_limit 5

    # Add header
    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tsscinames\tsskingdoms\tstaxids" \
        > "${qc}"/blast/"${sample}"/"${sample}".blastn.all.tsv.tmp

    cat "${qc}"/blast/"${sample}"/"${sample}".blastn.all.tsv \
        >> "${qc}"/blast/"${sample}"/"${sample}".blastn.all.tsv.tmp

    mv "${qc}"/blast/"${sample}"/"${sample}".blastn.all.tsv.tmp \
        "${qc}"/blast/"${sample}"/"${sample}".blastn.all.tsv

    # Best hit only
    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tsscinames\tsskingdoms\tstaxids" \
        > "${qc}"/blast/"${sample}"/"${sample}".blastn.bestHit.tsv.tmp

    cat "${qc}"/blast/"${sample}"/"${sample}".blastn.all.tsv \
        | sed '1d' \
        | sort -t $'\t' -k1,1g -k13,13gr \
        | sort -t $'\t' -uk1,1g \
        >> "${qc}"/blast/"${sample}"/"${sample}".blastn.bestHit.tsv.tmp

    mv "${qc}"/blast/"${sample}"/"${sample}".blastn.bestHit.tsv.tmp \
        "${qc}"/blast/"${sample}"/"${sample}".blastn.bestHit.tsv
}

export -f blast

find "$polished" -type f -name "*_nanopolished.fasta" \
    | parallel  --bar \
                --env blast \
                --env cpu \
                --env maxProc \
                --env qc \
                --jobs "$maxProc" \
                "blast {}"

#clean blast index files
find "$ordered" -type f ! -name "*.fasta" -exec rm {} \;


### QUAST ###

declare -a genomes=()
for i in $(find "$polished" -type f -name "*_nanopolished.fasta"); do 
    genomes+=("$i")
done

source activate quast

quast.py \
    --output-dir "${qc}"/quast/all \
    --threads "$cpu" \
    --min-contig "$smallest_contig" \
    --est-ref-size "$size" \
    ${genomes[@]}


# Make quast report on individual assembly
function run_quast()
{
    sample=$(cut -d '_' -f 1 <<< $(basename "$1"))

    quast.py \
        -m "$smallest_contig" \
        -t $((cpu/maxProc)) \
        -o "${qc}"/quast/"$sample" \
        --min-contig "$smallest_contig" \
        --est-ref-size "$size" \
        $1  # don't put in quotes
}

#make function available to parallel
export -f run_quast  # -f is to export functions

#run paired-end merging on multiple samples in parallel
find "$polished" -type f -name "*_nanopolished.fasta" \
    | parallel  --bar \
                --env run_quast \
                --env maxProc \
                --env cpu \
                --env qc \
                --env smallest_contig \
                --jobs "$maxProc" \
                "run_quast {}"

source deactivate

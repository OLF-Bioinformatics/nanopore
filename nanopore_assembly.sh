#!/bin/bash



# Process and assemble Oxford nanopore sequencer data in fast5 format
# Include file format conversion, stats, qc and de novo assembly



#script version
version="0.1"


######################
#                    #
#    User Defined    #
#                    #
######################


# Analysis folder
export baseDir=""${HOME}"/analyses/lambda_nanopore"

# Reads
export fast5="/media/3tb_hdd/data/lambda_nanopore/fast5"

# Database to use for metagomic analysis of raw data (contamination)
db="/media/3tb_hdd/db/centrifuge/p_compressed+h+v"
# db="/media/3tb_hdd/db/centrifuge/nt"

# Program location
export prog=""${HOME}"/prog/nanopolish/scripts"

# Maximum number of cores used per sample for parallel processing
# A highier value reduces the memory footprint.
export maxProc=12

# Assembly name
prefix="lambda"


#######################
#                     #
#   Data Stucture     #
#                     #
#######################


# Folder structure
logs=""${baseDir}"/logs"
qc=""${baseDir}"/qc"
export fastq=""${baseDir}"/fastq"
assemblies=""${baseDir}"/assemblies"
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


######################
#                    #
#     Resources      #
#                    #
######################


# Computer performance
export cpu=$(nproc) #total number of cores
mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB


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
if hash poretools 2>/dev/null; then  # if installed
    poretools -v | tee -a "${logs}"/log.txt
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

# Canu
if hash canu 2>/dev/null; then
    canu --version | tee -a "${logs}"/log.txt
else
    echo >&2 "canu was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# BWA
if hash bwa 2>/dev/null; then
    v=$(bwa 2>&1 1>/dev/null | grep -F "Version" | cut -d " " -f 2)
    echo "bwa $v" | tee -a "${logs}"/log.txt
else
    echo >&2 "bwa was not found. Aborting." | tee -a "${logs}"/log.txt
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
    echo "nanopolish v0.6" | tee -a "${logs}"/log.txt  # not version option with this software, yet
else
    echo >&2 "nanopolish was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi


#check Centrifuge datase
if [ -s "$db" ]; then
    echo -e "\nCentrifuge database: $(basename "$db")" | tee -a "${logs}"/log.txt
else
    echo ""
    exit 1
fi

# Check if reads folder has fast5 files
if [ $(ls "$fast5" | wc -l) -eq 0 ]; then
    echo "No .fast5 files are present in the provided data folder"
    exit 1
fi


##########
#        #
#   QC   #
#        #
##########


# Collector’s curve of the yield
poretools yield_plot \
    --plot-type reads \
    --saveas "${qc}"/yield.png \
    --savedf "${qc}"/yield.tsv \
    "$fast5"

# Histogram of read sizes
poretools hist \
    --saveas "${qc}"/sizes.png \
    "$fast5"


#########################
#                       #
#   Converto to fastq   #
#                       #
#########################


function Fast5_to_fastq()
{
    name=$(basename "${1%.fast5}")
    poretools fastq "$1" > "${fastq}"/"${name}".fastq
}

# Make function available to parallel
export -f Fast5_to_fastq  # -f is to export functions

# Run in parallel
find "$fast5" -type f -name "*.fast5" \
    | parallel  --env Fast5_to_fastq \
                --env fastq \
                --jobs "$maxProc" \
                'Fast5_to_fastq {}'

# Merge all fastq into one file
cat "${fastq}"/*.fastq | pigz > "${fastq}"/all_reads.fastq.gz

# Delete unmerged fastq
for i in $(find "$fast5" -type f -name "*.fast5"); do
    name=$(basename "${i%.fast5}")
    rm -f "${fastq}"/"${name}".fastq
done


####################
#                  #
#   FastQc - Raw   #
#                  #
####################


# Create folder to store report
[ -d "${qc}"/fastqc/raw ] || mkdir -p "${qc}"/fastqc/raw

# Run fastQC
fastqc \
    --o  "${qc}"/fastqc/raw \
    --noextract \
    --threads "$cpu" \
    "${fastq}"/all_reads.fastq.gz


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
    -U "${fastq}"/all_reads.fastq.gz \
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


# canu
canu \
    -p "$prefix" \
    -d "$assemblies" \
    genomeSize=48.5k \
    -nanopore-raw "${fastq}"/all_reads.fastq.gz

# Polish assembly

#extrac reads in fasta
#convert fastq to fasta
zcat "${fastq}"/all_reads.fastq.gz \
    | sed -n '1~4s/^@/>/p;2~4p' \
    > "${polished}"/reads.fasta

# Index draft genome
bwa index "${assemblies}"/"${prefix}".contigs.fasta
# bwa index "${assemblies}"/"${prefix}".unassembled.fasta

# Align the basecalled reads to the draft sequence
bwa mem \
    -x ont2d \
    -t "$scpu" \
    "${assemblies}"/"${prefix}".contigs.fasta \
    "${polished}"/reads.fasta | \
    samtools sort -o "${polished}"/reads.sorted.bam -

#index bam file
samtools index "${polished}"/reads.sorted.bam

#run nanopolish
python "${prog}"/nanopolish_makerange.py "${assemblies}"/"${prefix}".unassembled.fasta \
    | parallel --results "${polished}"/nanopolish.results -P $((cpu/maxProc)) \
        nanopolish variants \
            --consensus "${polished}"/polished.{1}.fa \
            -w {1} \
            -r "${polished}"/reads.fasta \
            -b "${polished}"/reads.sorted.bam \
            -g "${assemblies}"/"${prefix}".unassembled.fasta \
            -t $((cpu/maxProc)) \
            --min-candidate-frequency 0.1

#merge individual segments
python "${prog}"/nanopolish_merge.py \
    "${polished}"/polished.*.fa \
    > "${polished}"/"${prefix}"_polished_genome.fasta

#Cleanup
rm -f "${polished}"/polished.*.fa

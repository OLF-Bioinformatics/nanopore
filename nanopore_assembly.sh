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


#Analysis folder
export baseDir=""${HOME}"/analyses/lambda_nanopore"

#reads
export fast5="/media/3tb_hdd/data/lambda_nanopore/fast5"

#Database to use for metagomic analysis of raw data (contamination)
db="/media/3tb_hdd/db/centrifuge/p_compressed+h+v"
# db="/media/3tb_hdd/db/centrifuge/nt"

#program location
export prog=""${HOME}"/prog/nanopolish/scripts"

#Maximum number of cores used per sample for parallel processing
#A highier value reduces the memory footprint.
export maxProc=12

#Assembly name
prefix="lambda"


#######################
#                     #
#   Data Stucture     #
#                     #
#######################


#Folder structure
logs=""${baseDir}"/logs"
qc=""${baseDir}"/qc"
export fastq=""${baseDir}"/fastq"
assemblies=""${baseDir}"/assemblies"
export polished=""${baseDir}"/polished"
aligned=""${baseDir}"/aligned"


#create folders if do not exist
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


#computer performance
export cpu=$(nproc) #total number of cores
mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB


################
#              #
#     Log      #
#              #
################


#Date
echo -e "$(date)\n" | tee "${logs}"/log.txt
echo -e "User: $(whoami)" | tee -a "${logs}"/log.txt
echo -e "Processors: "$cpu"" | tee -a "${logs}"/log.txt
echo -e "Memory: "$mem"G" | tee -a "${logs}"/log.txt

#script version
echo -e "\nnanopore_assembly.sh version "$version"\n" | tee -a "${logs}"/log.txt  # $0

#log software versions
#poretools
if hash poretools 2>/dev/null; then  # if installed
    poretools -v | tee -a "${logs}"/log.txt
else
    echo >&2 "poretools was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi


##########
#        #
#   QC   #
#        #
##########


#collector’s curve of the yield
poretools yield_plot \
    --plot-type reads \
    --saveas "${qc}"/yield.png \
    --savedf "${qc}"/yield.tsv \
    "$fast5"

#histogram of read sizes
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

#make function available to parallel
export -f Fast5_to_fastq  # -f is to export functions


#run trimming on multiple samples in parallel
find "$fast5" -type f -name "*.fast5" \
    | parallel  --env Fast5_to_fastq \
                --env fastq \
                --jobs "$maxProc" \
                'Fast5_to_fastq {}'

#merge all fastq into one file
cat "${fastq}"/*.fastq | pigz > "${fastq}"/all_reads.fastq.gz

#delete unmerged fastq
for i in $(find "$fast5" -type f -name "*.fast5"); do
    name=$(basename "${i%.fast5}")
    rm -f "${fastq}"/"${name}".fastq
done


####################
#                  #
#   FastQc - Raw   #
#                  #
####################


[ -d "${qc}"/fastqc/raw ] || mkdir -p "${qc}"/fastqc/raw

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


#create folder
[ -d "${qc}"/centrifuge ] || mkdir -p "${qc}"/centrifuge

#run centrifuge
centrifuge \
    -p "$cpu" \
    -t \
    --seed "$RANDOM" \
    -x "$db" \
    -U "${fastq}"/all_reads.fastq.gz \
    --report-file "${qc}"/centrifuge/"${prefix}"_report.tsv \
    > "${qc}"/centrifuge/"${prefix}".tsv


#Prepare result for display with Krona
cat "${qc}"/centrifuge/"${prefix}".tsv | \
    cut -f 1,3 | \
    ktImportTaxonomy /dev/stdin -o  "${qc}"/centrifuge/"${prefix}".html

#visualize the resutls in Firefow browser
firefox file://"${qc}"/centrifuge/"${prefix}".html &




########################
#                      #
#   de novo assembly   #
#                      #
########################


#canu
canu \
    -p "$prefix" \
    -d "$assemblies" \
    genomeSize=48.5k \
    -nanopore-raw "${fastq}"/all_reads.fastq.gz

#polish assembly

#extrac reads in fasta
#convert fastq to fasta
zcat "${fastq}"/all_reads.fastq.gz \
    | sed -n '1~4s/^@/>/p;2~4p' \
    > "${polished}"/reads.fasta

#index draft genome
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

#!/bin/bash



# Process and assemble Oxford nanopore sequencer data in fast5 format
# Include file format conversion, stats, qc and de novo assembly



#script version
version="0.1.1"


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

# Maximum number of cores used per sample for parallel processing
# A highier value reduces the memory footprint.
export maxProc=48

# Assembly name
export prefix="YpD1"

# Estimated genome size in bp
size=4600000


######################
#                    #
#     Resources      #
#                    #
######################


# Computer performance
export cpu=$(nproc) #total number of cores
mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
memJava="-Xmx"$mem"g"


#######################
#                     #
#   Data Stucture     #
#                     #
#######################


# Folder structure
logs=""${baseDir}"/logs"
qc=""${baseDir}"/qc"
export fastq=""${baseDir}"/fastq"
basecalled=""${baseDir}"/basecalled"
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
    "$fast5/pass"

# Histogram of read sizes
poretools hist \
    --saveas "${qc}"/sizes.png \
    "$fast5/pass"


########################
#                      #
#   Convert to fastq   #
#                      #
########################


function Fast5_to_fastq()
{
    name=$(basename "${1%.fast5}")
    poretools fastq "$1" > "${fastq}"/"${name}".fastq
    # poretools fastq "$1" | gzip > "${fastq}"/"${name}".fastq  # gzip iis OK to use beacause they are all very small files
    # poretools fastq "$1" | pigz >> "${fastq}"/"${prefix}".fromfast5.fastq.gz  # seems to produce a corrupted fastq file
}

# Make function available to parallel
export -f Fast5_to_fastq  # -f is to export functions

# Run in parallel
# TODO -> Add log
find "$fast5/pass" -type f -name "*.fast5" \
    | parallel  --bar \
                -k \
                --env Fast5_to_fastq \
                --env fastq \
                --env prefix \
                --jobs "$maxProc" \
                'Fast5_to_fastq {}'

# Merge all fastq into one file
# cat "${fastq}"/*.fastq | pigz > "${fastq}"/"${prefix}".fastq.gz # Does not work because too many files: bash: /bin/cat: Argument list too long
find "$fastq" -type f -name "*.fastq" -exec cat {} + \
    | pigz > "${fastq}"/"${prefix}".fastq.gz

# Delete unmerged fastq
# rm -f "${fastq}"/*.fastq #bash: /bin/rm: Argument list too long
find "$fastq" -type f -name "*.fastq" -exec rm {} +


# # https://github.com/rrwick/Fast5-to-Fastq
# # Aim for a 100x coverage (5MB x 100) -> not a good idea if metagenomic sample...
# # 2000bp is longer that most repeats in bacteria.
# ##--target_bases 500000000 \
# python3 "${prog}"/Fast5-to-Fastq/fast5_to_fastq.py \
#     --min_length 2000 \
#     "$fast5/pass" \
#     | pigz > "${fastq}"/"${prefix}".fastq.gz


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
bbduk.sh "$memJava" \
    threads="$cpu" \
    in="${fastq}"/"${prefix}"_t.fastq.gz \
    minavgquality=12 \
    minlength=1000 \
    out="${fastq}"/"${prefix}"_trimmed.fastq.gz \
    ziplevel=9 \
    2> >(tee -a "${logs}"/"${prefix}"_bbduk.txt)


###################
#                 #
#   Basecalling   #
#                 #
###################


# Using nanonetcall
# To check what the failed sequences look like

# Select the class of read to call
export read_class="pass"
# export read_class="fail"
# export read_class="skip"

# The calling function
function basecall()
{
    in_name=$(basename "$1")
    out_name="${in_name%.fast5}.fastq"

    nanonetcall \
        --fastq \
        "$1" \
        > "${fast5}"/../fastq/"${read_class}"/"$out_name"

    if [[ ! -s "${fast5}"/../fastq/"${read_class}"/"$out_name" ]]; then
        rm "${fast5}"/../fastq/"${read_class}"/"$out_name"
    fi
}

# Make function available to parallel
export -f basecall

[ -d "${fast5}"/../fastq/"${read_class}" ] || mkdir -p "${fast5}"/../fastq/"${read_class}"

# On Failed reads
find "${fast5}"/"${read_class}" -type f -name "*.fast5" |
    parallel    --bar \
                --env basecall \
                --env fast5 \
                'basecall {}'

cat "${fast5}"/../fastq/"${read_class}"/*.fastq \
    | pigz > "${fast5}"/../fastq/"${read_class}"/salmonella_lettuce_"${read_class}".fastq.gz

rm "${fast5}"/../fastq/"${read_class}"/*.fastq


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


# canu
# canu \
#     -p "$prefix" \
#     -d "$assemblies" \
#     genomeSize="$size" \
#     -nanopore-raw "${fastq}"/"${prefix}"_trimmed.fastq.gz

# Assemble
for i in $(find "$fastq" -type f -name "*.fastq.gz"); do
    sample=$(cut -d "_" -f 1 <<< $(basename "$i"))
    
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
done


#################
#               #
#   Polishing   #
#               #
#################


# Polish
for i in $(find "$fastq" -type f -name "*.fastq.gz"); do
    sample=$(cut -d "_" -f 1 <<< $(basename "$i"))

    [ -d "${polished}"/"$sample" ] || mkdir -p "${polished}"/"$sample"

    # Convert fastq to fasta
    zcat "$i" \
        | sed -n '1~4s/^@/>/p;2~4p' \
        > "${polished}"/"${sample}"/"${sample}".fasta

    # Index draft genome
    bwa index \
        "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta

    # Align the basecalled reads to the draft sequence
    # bwa mem \
    #     -x ont2d \
    #     -t "$cpu" \
    #     "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta \
    #     "${polished}"/"${sample}"/"${sample}".fasta | \
    #     samtools sort -@ "$cpu" -o "${polished}"/"${sample}"/"${sample}".bam -

    minimap2 \
        -ax map-ont \
        -t "$cpu" \
        "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta \
        "${polished}"/"${sample}"/"${sample}".fasta | \
    samtools sort -@ "$cpu" -o "${polished}"/"${sample}"/"${sample}".bam -

    # Index bam file
    samtools index -@ "$cpu" "${polished}"/"${sample}"/"${sample}".bam

    # Index reads for nanopolish
    nanopolish index \
        -d /media/6tb_raid10/data/salmonella_human_birds/nanopore/20180503_1623_Wild_Bird_8samples_2018-05-03/fast5 \
        -s "${basecalled}"/sequencing_summary.txt \
        "${polished}"/"${sample}"/"${sample}".fasta

    #run nanopolish
    python /home/bioinfo/prog/nanopolish/scripts/nanopolish_makerange.py \
        "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta \
        | parallel --results "${polished}"/"${sample}"/nanopolish.results -j "$cpu" \
            nanopolish variants \
                --consensus "${polished}"/"${sample}"/"${sample}"_nanopolished.{1}.fa \
                -w {1} \
                -r "${polished}"/"${sample}"/"${sample}".fasta \
                -b "${polished}"/"${sample}"/"${sample}".bam \
                -g "${assemblies}"/unicycler/"${sample}"/"${sample}"_ordered.fasta \
                -t 1 \
                --min-candidate-frequency 0.1 \
                -q dcm,dam \
                --fix-homopolymers

    #merge individual segments
    python  /home/bioinfo/prog/nanopolish/scripts/nanopolish_merge.py \
        "${polished}"/"${sample}"/"${sample}"_nanopolished.*.fa \
        > "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta

    #Cleanup
    find "${polished}"/"$sample" -type f ! -name "*_nanopolished.fasta" -exec rm -rf {} \;
    rm -rf "${polished}"/"${sample}"/nanopolish.results
done


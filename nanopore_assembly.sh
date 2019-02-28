#!/bin/bash



# Process and assemble Oxford nanopore sequencer data in fast5 format
# Include file format conversion, stats, qc and de novo assembly



#script version
version="0.3.0"


######################
#                    #
#    User Defined    #
#                    #
######################


# Analysis folder
export baseDir=""${HOME}"/analyses/burkholderia_nanopore_batch1"

# Reads
export fast5="/media/2TB_NVMe/burkholderia_fast5_1/20180823_1514_Bmallei/fast5"
export fastq_reads="/media/30tb_raid10/data/burkholderia/illumina/"

# Database to use for metagomic analysis of raw data (contamination)
db="/media/30tb_raid10/db/centrifuge/2017-10-12_bact_vir_h"
# db="/media/6tb_raid10/db/centrifuge/nt"

# Program location
export prog=""${HOME}"/prog"
export scripts=""${HOME}"/scripts"

# Maximum number of cores used per sample for parallel processing
# A highier value reduces the memory footprint.
export maxProc=4

# Estimated genome size in bp
export size=4850000

# Set assembler to use: unicycler, unicycler_hybrid, canu or flye
export assembler="unicycler"

# Set smallest contig size for assemblies
export smallest_contig=1000

#Annotation
export kingdom="Bacteria"
export genus="Salmonella"
export species="enterica"
export gram="neg"
export locus_tag="TOCHANGE"
export centre="OLF"


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
export trimmed=""${baseDir}"/trimmed"
export filtered=""${baseDir}"/filtered"
export basecalled=""${baseDir}"/basecalled"
export assemblies=""${baseDir}"/assemblies"
export polished=""${baseDir}"/polished"
export annotation=""${baseDir}"/annotation"
export phaster=""${baseDir}"/phaster"
export amr=""${baseDir}"/amr"


# Create folders if do not exist
# "||" if test is false
# "&&" if test is true
[ -d "$baseDir" ] || mkdir -p "$baseDir"
[ -d "$logs" ] || mkdir -p "$logs"
[ -d "$qc" ] || mkdir -p "$qc"
[ -d "$fastq" ] || mkdir -p "$fastq"
[ -d "$trimmed" ] || mkdir -p "$trimmed"
[ -d "$filtered" ] || mkdir -p "$filtered"
[ -d "$assemblies" ] || mkdir -p "$assemblies"
[ -d "$polished" ] || mkdir -p "$polished"
[ -d "$annotation" ] || mkdir -p "$annotation"
[ -d "$phaster" ] || mkdir -p "$phaster"
[ -d "$amr" ] || mkdir -p "$amr"


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
    echo "$version" | tee -a "${logs}"/log.txt
else
    echo >&2 "nanopolish was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi


#check Centrifuge database
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


######################
#                    #
#   Demultiplexing   #
#                    #
######################


# Since Guppy, we need to convert multi-reads fast5 to single-read fast5
# about 50min to convert 2M reads.
multi_to_single_fast5 \
    --input_path "$fast5" \
    --save_path "${fast5}"/singles \
    --threads "$cpu"

rm "${fast5}"/*.fast5

# Using Deepbinner
# $((cpu/2)) because we want to set it to number of physical cores
# https://www.tensorflow.org/guide/performance/overview
# if used ligation barcoding kit
deepbinner realtime \
    --in_dir "${fast5}"/singles \
    --out_dir "${fast5}"/demultiplexed \
    --native \
    --intra_op_parallelism_threads $((cpu/2)) \
    --omp_num_threads $((cpu/2))

# if used rapid barcoding kit
deepbinner realtime \
    --in_dir "${fast5}"/singles \
    --out_dir "${fast5}"/demultiplexed \
    --rapid \
    --intra_op_parallelism_threads $((cpu/2)) \
    --omp_num_threads $((cpu/2))

# Cleanup
rm -rf "${fast5}"/singles

#convert back to multi
# Increase Guppy speed?
for i in $(find "${fast5}"/demultiplexed -mindepth 1 -maxdepth 1 -type d); do
    barcode="$(basename "$i")"
    outdir="${fast5}"/demul_multi/"$barcode"

    [ -d "$outdir" ] || mkdir -p "$outdir"

    single_to_multi_fast5 \
        --input_path "$i" \
        --save_path "$outdir" \
        --threads "$cpu"
done

# cleanup
rm -rf "${fast5}"/demultiplexed


###################
#                 #
#   Basecalling   #
#                 #
###################


# Use Guppy
function call_bases_guppy()
{
    kit="$1"  # "SQK-LSK108"
    flow="$2"  # "FLO-MIN106"
    fast5_folder="$3"  # top folder - will look recursively for fast5 files because "-r" option
    fastq_folder="$4"

    guppy_bin_folder=$(dirname $(which guppy_basecaller))
    guppy_data_folder="${guppy_bin_folder%/bin}/data"

    #####   CPU mode   #####

    # Using Flip-Flop mode
    guppy_basecaller \
        --input_path "$fast5_folder" \
        --save_path "$fastq_folder" \
        --kit "$kit" \
        --flowcell "$flow" \
        --recursive \
        --records_per_fastq 0 \
        --disable_pings \
        --model_file ""${guppy_data_folder}"/template_r9.4.1_450bps_large_flipflop.jsn" \
        --calib_detect \
        --calib_reference "lambda_3.6kb.fasta" \
        --hp_correct 1 \
        --enable_trimming 1 \
        --trim_strategy 'dna' \
        --qscore_filtering \
        --num_callers 2 \
        --cpu_threads_per_caller 24

    #####   CPU mode   #####  

    # ( time sudo nvidia-docker run \
    #     -v "${DATA_DIR}":/home/workspace \
    #     docker.io/duceppemo/guppy_gpu:"${BIN_VERSION}" \
    #     /usr/bin/guppy_basecaller \
    #         --input_path /home/workspace/"$BC" \
    #         --save_path /home/workspace/basecalled/"$BC" \
    #         --kit "$kit" \
    #         --flowcell "$flow" \
    #         --recursive \
    #         --records_per_fastq 0 \
    #         --disable_pings \
    #         --model_file template_r9.4.1_450bps_large_flipflop.jsn \
    #         --calib_detect \
    #         --calib_reference lambda_3.6kb.fasta \
    #         --hp_correct 1 \
    #         --enable_trimming 1 \
    #         --trim_strategy 'dna' \
    #         --trim_threshold 2.5 \
    #         --trim_min_events 3 \
    #         --qscore_filtering \
    #         --min_qscore 7 \
    #         --gpu_runners_per_device 2 \
    #         --chunk_size 1000 \
    #         --chunks_per_runner 1000 \
    #         --device "cuda:0"
    # ) >""${LOG_DIR}"/"${BC}".call_variants.log" 2>&1
}

export -f call_bases_guppy

for i in $(find "${fast5}"/demul_multi -mindepth 1 -maxdepth 1 -type d); do
    echo $(date)
    barcode=$(basename "$i"); echo $barcode
    echo "Basecalling "$barcode" with Guppy in CPU mode..."
    call_bases_guppy "SQK-LSK108" "FLO-MIN106" "$i" "${basecalled}"/"$barcode"
    # call_bases_guppy "SQK-RBK004" "FLO-MIN106" "$i" "${basecalled}"/"$barcode"
    echo $(date)
done

# Compress the basecalled fastq
for i in $(find "$basecalled" -mindepth 2 -maxdepth 2 -type d); do  # pass and fail
    flag=$(basename "$i")
    barcode=$(basename $(dirname "$i"))
    # echo "${i}"/"${barcode}"_"${flag}".fastq.gz
    cat "${i}"/*.fastq | pigz > "${i}"/"${barcode}"_"${flag}".fastq.gz
    # # Remove non compressed files
    rm "${i}"/*.fastq
done

### Merge baecalling results ###
# list summary files
sumfiles=($(find "$basecalled" -type f -name "*sequencing_summary.txt"))

#write header
echo -e "filename\tread_id\trun_id\tchannel\tstart_time\tduration\tnum_events\tpasses_filtering\
\ttemplate_start\tnum_events_template\ttemplate_duration\tsequence_length_template\
\tmean_qscore_template\tstrand_score_template\tmedian_template\tmad_template\tcalibration_strand_genome\
\tcalibration_strand_genome_start\tcalibration_strand_genome_end\tcalibration_strand_strand_start\
\tcalibration_strand_strand_end\tcalibration_strand_num_insertions\tcalibration_strand_num_deletions\
\tcalibration_strand_num_aligned\tcalibration_strand_num_correct\tcalibration_strand_identity\
\tcalibration_strand_accuracy\tcalibration_strand_score" \
    > "${basecalled}"/sequencing_summary.txt

# Merge sequencing_summary.txt files, skipping header
# add the barcode information at the 20th column
# if not demultiplexed while basecalling with Albacore
for f in "${sumfiles[@]}"; do
    barcode=$(basename $(dirname "$f"))
    cat "$f" | sed -e '1d' \
        | awk -v bc="$barcode" -F $'\t' 'BEGIN {OFS = FS} {$20 = bc; print}' \
        >> "${basecalled}"/sequencing_summary.txt
done

# Reorganize using parent "pass" and "fail" folders
state=('pass' 'fail')
for b in $(find "$basecalled" -mindepth 1 -maxdepth 1 -type d); do
    barcode=$(basename "$b")
    for s in "${state[@]}"; do
        for c in $(find "${b}"/"${s}" -type f -name "*.fastq.gz"); do
            # Copy and rename file
            [ -d "${basecalled}"/"${s}"/"${barcode}" ] || mkdir -p "${basecalled}"/"${s}"/"${barcode}"
            mv "$c" "${basecalled}"/"${s}"/"${barcode}"
        done
    done
done


##########
#        #
#   QC   #
#        #
##########


# nanoQC on "sequencing_summary.txt" file
[ -d "${qc}"/nanoQC/raw/summary ] || mkdir -p "${qc}"/nanoQC/raw/summary
python3 /home/bioinfo/PycharmProjects/nanoQC/nanoQC_Guppy_v1.py \
    -s "${basecalled}"/sequencing_summary.txt \
    -o "${qc}"/nanoQC/raw/summary

# nanoQC on fastq files
[ -d "${qc}"/nanoQC/raw/fastq ] || mkdir -p "${qc}"/nanoQC/raw/fastq
python3 /home/bioinfo/PycharmProjects/nanoQC/nanoQC_Guppy_v1.py \
    -f "$basecalled" \
    -o "${qc}"/nanoQC/raw/fastq


################
#              #
#   Trimming   #
#              #
################


# Trim adapters
function trim()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    porechop \
        -i "$1" \
        -o "${trimmed}"/"${sample}"_trimmed.fastq.gz \
        --threads $((cpu/maxProc)) \
        | tee -a "${logs}"/trimming/"${sample}"_porechop.log

    # cleanup the log file by removing the progress lines
    cat "${logs}"/trimming/"${sample}"_porechop.log \
        | grep -vF '%' \
        > "${logs}"/trimming/"${sample}"_porechop.log.tmp

    mv "${logs}"/trimming/"${sample}"_porechop.log.tmp \
        "${logs}"/trimming/"${sample}"_porechop.log
}

export -f trim

[ -d "${logs}"/trimming ] || mkdir -p "${logs}"/trimming

find "${basecalled}/pass" -type f -name "*.fastq.gz" ! -name "*unclassified*" \
    | parallel  --bar \
                --env trim \
                --env trimmed \
                --env logs \
                --env cpu \
                --env maxProc \
                --jobs "$maxProc" \
                'trim {}'


#################
#               #
#   Filtering   #
#               #
#################


### Careful here. Filtered reads sometimes lead to a poorer quality assembly.
function filter()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    # https://github.com/rrwick/Filtlong
    # --min_mean_q 90 -> phred 10
    # --min_length 1000 \
    # --target_bases -> max 100X coverage
    filtlong \
        --target_bases $((size*100)) \
        --keep_percent 90 \
        "$1" \
        | pigz > "${filtered}"/"${sample}"_filtered.fastq.gz
        # 2> >(tee "${logs}"/filtering/"${sample}".txt) \
}

export -f filter

[ -d "${logs}"/filtering ] || mkdir -p "${logs}"/filtering

find "$trimmed" -type f -name "*_trimmed.fastq.gz" \
    | parallel  --bar \
                --env filter \
                --env trimmed \
                --env size \
                --env filtered \
                --env logs \
                --jobs "$maxProc" \
                'filter {}'

# Cleanup trimmed reads
# rm "${trimmed}"/*


######################
#                    #
#   Quality Control  #
#                    #
######################


### nanoQC ###

# nanoQC on fastq file
[ -d "${qc}"/nanoQC/trimmed ] || mkdir -p "${qc}"/nanoQC/trimmed
python3 /home/bioinfo/PycharmProjects/nanoQC/nanoQC_Guppy_v1.py \
    -f "$trimmed" \
    -o "${qc}"/nanoQC/trimmed

# nanoQC on fastq file
[ -d "${qc}"/nanoQC/filtered ] || mkdir -p "${qc}"/nanoQC/filtered
python3 /home/bioinfo/PycharmProjects/nanoQC/nanoQC_Guppy_v1.py \
    -f "$filtered" \
    -o "${qc}"/nanoQC/filtered

### FastQC###

# Create folder to store report
[ -d "${qc}"/fastqc/raw ] || mkdir -p "${qc}"/fastqc/raw
[ -d "${qc}"/fastqc/trimmed ] || mkdir -p "${qc}"/fastqc/trimmed
[ -d "${qc}"/fastqc/filtered ] || mkdir -p "${qc}"/fastqc/filtered

function run_fastqc()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    [ -d "${2}"/"$sample" ] || mkdir -p "${2}"/"$sample"

    fastqc \
        --o "${2}"/"$sample" \
        --noextract \
        --threads $((cpu/maxProc)) \
        "$1"
}

export -f run_fastqc

#raw
find "$basecalled" -type f -name "*.fastq.gz" \
    | parallel --bar run_fastqc {} "${qc}"/fastqc/raw

#trimmed
find "$trimmed" -type f -name "*.fastq.gz" \
    | parallel --bar run_fastqc {} "${qc}"/fastqc/trimmed

#filtered
find "$filtered" -type f -name "*.fastq.gz" \
    | parallel --bar run_fastqc {} "${qc}"/fastqc/filtered

#Merge all FastQC reports together
multiqc \
    -o "${qc}"/fastqc/raw \
    -n merged_reports.html \
    "${qc}"/fastqc/raw
multiqc \
    -o "${qc}"/fastqc/trimmed \
    -n merged_reports.html \
    "${qc}"/fastqc/trimmed
multiqc \
    -o "${qc}"/fastqc/filtered \
    -n merged_reports.html \
    "${qc}"/fastqc/filtered

### Centrifuge ###

[ -d "${qc}"/centrifuge/trimmed ] || mkdir -p "${qc}"/centrifuge/trimmed

function run_centrifuge()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    [ -d "${2}"/"${sample}" ] || mkdir -p "${2}"/"${sample}"

    #build the command
    centrifuge \
        -p "$cpu" \
        -t \
        --seed "$RANDOM" \
        -x "$db" \
        -U "$1" \
        --report-file "${2}"/"${sample}"/"${sample}"_report.tsv \
        > "${2}"/"${sample}"/"${sample}".tsv

    cat "${2}"/"${sample}"/"${sample}".tsv | \
        cut -f 1,3 | \
        ktImportTaxonomy /dev/stdin -o "${2}"/"${sample}"/"${sample}".html

    # Visualize the resutls in Firefow browser
    firefox file://"${2}"/"${sample}"/"${sample}".html &
}

for i in $(find "$trimmed" -type f -name "*fastq.gz"); do
    run_centrifuge "$i" "${qc}"/centrifuge/trimmed
done



########################
#                      #
#   de novo assembly   #
#                      #
########################


function assemble()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))
    ass="$2"
    
    if [[ "$ass" == "unicycler" ]]; then
        unicycler \
            -l "$1" \
            -o "${assemblies}"/"${ass}"/"$sample" \
            -t $((cpu/maxProc)) \
            --verbosity 2 \
            --mode normal

        mv "${assemblies}"/"${ass}"/"${sample}"/assembly.fasta \
            "${assemblies}"/"${ass}"/"${sample}"/"${sample}".fasta
    elif [[ "$ass" == "unicycler_hybrid" ]]; then
        unicycler \
            -1 unmerged_R1.fastq.gz \
            -2 unmerged R2.fastq.gz \
            -s merged.fastq.gz \
            --no_correct \
            -l "$1" \
            -o "${assemblies}"/"${ass}"/"$sample" \
            -t $((cpu/maxProc)) \
            --verbosity 2 \
            --mode normal

        mv "${assemblies}"/"${ass}"/"${sample}"/assembly.fasta \
            "${assemblies}"/"${ass}"/"${sample}"/"${sample}".fasta
    elif [[ "$ass" == "canu" ]]; then
        canu \
            -p "$sample" \
            -d "${assemblies}"/"${ass}"/"$sample" \
            genomeSize="$size" \
            maxTreads=$((cpu/maxProc)) \
            -nanopore-raw "$1"

        # Maybe some cleanup
    elif [[ "$ass" == "flye" ]]; then
        flye \
            --nano-raw "$1" \
            --out-dir "${assemblies}"/"${ass}"/"$sample" \
            --genome-size "$size" \
            --iterations 3 \
            --threads $((cpu/maxProc))
        # Maybe some cleanup
    else
        echo 'Please chose one of the following assembler: "unicycler", "unicycler_hybrid", "canu" or "flye"'
    fi
}

export -f assemble

find "${trimmed}" -type f -name "*.fastq.gz" \
    | parallel  --bar \
                --env cpu \
                --env maxProc \
                --env assemble \
                --env assembler \
                --env assemblies \
                --env size \
                --jobs "$maxProc" \
                'assemble {} unicycler'


#################
#               #
#   Polishing   #
#               #
#################


# Polish with Medaka
# https://github.com/nanoporetech/medaka
source activate medaka

function polish_medaka()
{
    sample=$(basename "$1" ".fasta")

    [ -d "${polished}"/"$sample" ] || mkdir -p "${polished}"/"$sample"

    # -m r94 -> if default algorithm
    medaka_consensus \
        -i "${basecalled}"/pass/"${sample}"/"${sample}"_pass.fastq.gz \
        -d "$1" \
        -o "${polished}"/"$sample" \
        -t $((cpu/maxProc)) \
        -m r941_flip

    # Rename sample
    mv "${polished}"/"${sample}"/consensus.fasta \
        "${polished}"/"${sample}"/"${sample}"_medaka.fasta
}

export -f polish_medaka

find "$assemblies" -type f -name "*.fasta" \
    | parallel  --bar \
                --env polish_medaka \
                --env polished \
                --env cpu \
                --env maxProc \
                'polish_medaka {}'

source deactivate medaka

#Cleanup
find "${polished}" -mindepth 2 -type f ! -name "*.fasta" -exec rm -rf {} \;


# # Polish with nanopolish
# function polish()
# {
#     sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

#     [ -d "${polished}"/"$sample" ] || mkdir -p "${polished}"/"$sample"

#     # Convert fastq to fasta
#     zcat "$1" \
#         | sed -n '1~4s/^@/>/p;2~4p' \
#         > "${polished}"/"${sample}"/"${sample}"_reads.fasta

#     # Index draft genome
#     bwa index \
#         "${assemblies}"/unicycler/"${sample}"/"${sample}".fasta

#     minimap2 \
#         -ax map-ont \
#         -t $((cpu/maxProc)) \
#         "${assemblies}"/unicycler/"${sample}"/"${sample}".fasta \
#         "${polished}"/"${sample}"/"${sample}"_reads.fasta | \
#     samtools sort -@ $((cpu/maxProc)) -o "${polished}"/"${sample}"/"${sample}".bam -

#     # Index bam file
#     samtools index -@ $((cpu/maxProc)) "${polished}"/"${sample}"/"${sample}".bam

#     # Index reads for nanopolish
#     nanopolish index \
#         -d "$fast5" \
#         -s "${basecalled}"/sequencing_summary.txt \
#         "${polished}"/"${sample}"/"${sample}"_reads.fasta

#     # Run nanopolish
#     python "${prog}"/nanopolish/scripts/nanopolish_makerange.py \
#         "${assemblies}"/unicycler/"${sample}"/"${sample}".fasta \
#         | parallel  --results "${polished}"/"${sample}"/nanopolish.results \
#                     --env polished \
#                     --env assemblies \
#                     -P $((cpu/maxProc)) \
#                     nanopolish variants \
#                         --consensus "${polished}"/"${sample}"/"${sample}"_nanopolished.{1}.fa \
#                         --window={1} \
#                         --methylation-aware=dcm,dam \
#                         --fix-homopolymers \
#                         --threads=1 \
#                         --min-candidate-frequency 0.1 \
#                         --reads="${polished}"/"${sample}"/"${sample}"_reads.fasta \
#                         --bam="${polished}"/"${sample}"/"${sample}".bam \
#                         --genome="${assemblies}"/unicycler/"${sample}"/"${sample}".fasta

#                         'parallel_polishing {} "$sample"'

#             # nanopolish variants \
#             #     --consensus "${polished}"/"${sample}"/"${sample}"_nanopolished.{1}.fa \
#             #     --window={1} \
#             #     --reads="${polished}"/"${sample}"/"${sample}"_reads.fasta \
#             #     --bam="${polished}"/"${sample}"/"${sample}".bam \
#             #     --genome="${assemblies}"/unicycler/"${sample}"/"${sample}".fasta \
#             #     --threads=1 \
#             #     --methylation-aware=dcm,dam \
#             #     --fix-homopolymers

#     # Merge individual segments
#     python  "${prog}"/nanopolish/scripts/nanopolish_merge.py \
#         "${polished}"/"${sample}"/"${sample}"_nanopolished.*.fa \
#         > "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta

#     # Adjust max line width of fasta file
#     perl "${scripts}"/formatFasta.pl \
#         -i "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta \
#         -o "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta.tmp \
#         -w 80

#     # Rename header
#     # sed -i '/^>/s/$/_nanopolish/' "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta.tmp
    
#     mv "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta.tmp \
#         "${polished}"/"${sample}"/"${sample}"_nanopolished.fasta

#     # Cleanup all but fasta read file and polished assembly
#     find "${polished}"/"$sample" -type f ! -name "*_nanopolished.fasta" ! -name "*_reads.fasta" -exec rm -rf {} \;
#     rm -rf "${polished}"/"${sample}"/nanopolish.results
# }

# export -f polish

# find "$filtered" -type f -name "*.fastq.gz" | \
# parallel    --bar \
#             --env polish \
#             --env qc \
#             --env polished \
#             --env assemblies \
#             --env basecalled \
#             --env prog \
#             --env cpu \
#             --env maxProc \
#             --jobs "$maxProc"  \
#             "polish {}"


# Compare pre- and post-nanopolished genomes
function compare_assemblies()
{
    export sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    [ -d "${qc}"/mummer/"$sample" ] || mkdir -p "${qc}"/mummer/"$sample"

    nucmer \
        --mum \
        --threads=$((cpu/maxProc)) \
        --delta "${qc}"/mummer/"${sample}"/"${sample}".delta \
        "${polished}"/"${sample}"/"${sample}"_medaka.fasta \
        "${assemblies}"/"${assembler}"/"${sample}"/"${sample}".fasta

    mummerplot \
        -p "${qc}"/mummer/"${sample}"/"${sample}" \
        -title 'Pre versus Post medaka' \
        --layout \
        --large \
        --color \
        --png \
        "${qc}"/mummer/"${sample}"/"${sample}".delta

    dnadiff \
        -p "${qc}"/mummer/"${sample}"/"${sample}" \
        -d "${qc}"/mummer/"${sample}"/"${sample}".delta
}

export -f compare_assemblies

find "$polished" -type f -name "*.fasta" | \
parallel    --bar \
            --env compare_assemblies \
            --env polished \
            --env assemblies \
            --env assembler \
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
        "${trimmed}"/"${sample}"_trimmed.fastq.gz | \
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

    #Remove reads in fasta format  -> for nanopolish
    # rm "${polished}"/"${sample}"/"${sample}"_reads.fasta
}

export -f get_coverage

[ -d "${qc}"/coverage ] || mkdir -p "${qc}"/coverage
echo -e "Sample\tAverage_Cov" > "${qc}"/coverage/average_cov.tsv

find "$polished" -type f -name "*.fasta" | \
    parallel    --bar \
                --env get_coverage \
                --env cpu \
                --env maxProc \
                --env qc \
                --jobs "$maxProc" \
                "get_coverage {}"

# Qualimap
function run_qualimap()
{
    sample=$(cut -d "." -f 1 <<< $(basename "$1"))
    
    [ -d "${qc}"/coverage/qualimap/"$sample" ] || mkdir -p "${qc}"/coverage/qualimap/"$sample"

    qualimap bamqc \
        --paint-chromosome-limits \
        -bam "$1" \
        --java-mem-size="${mem}"G \
        -nt $((cpu/maxProc)) \
        -outdir "${qc}"/coverage/qualimap/"$sample" \
        -outfile "${sample}" \
        -outformat HTML

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
            --jobs "$maxProc" \
            'run_qualimap {}'


# Blast genomes on nr
function blast()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "${1%.fasta}"))

    [ -d "${qc}"/blast/"$sample" ] || mkdir -p "${qc}"/blast/"$sample"

    blastn \
        -db nt \
        -query "$1" \
        -out "${qc}"/blast/"${sample}"/"${sample}".blastn.all.tsv \
        -evalue "1e-10" \
        -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames sskingdoms staxids' \
        -num_threads $((cpu/maxProc)) \
        -max_target_seqs 20

    # Add header
    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tsscinames\tsskingdoms\tstaxids" \
        > "${qc}"/blast/"${sample}"/"${sample}".blastn.all.tsv.tmp

    # Sort hits
    cat "${qc}"/blast/"${sample}"/"${sample}".blastn.all.tsv \
        | sed '1d' \
        | sort -t $'\t' -k1,1g -k13,13gr \
        >> "${qc}"/blast/"${sample}"/"${sample}".blastn.all.tsv.tmp

    # Replace original file
    mv "${qc}"/blast/"${sample}"/"${sample}".blastn.all.tsv.tmp \
        "${qc}"/blast/"${sample}"/"${sample}".blastn.all.tsv

    # Best hits only
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

# export -f blast

# find "$polished" -type f -name "*_nanopolished.fasta" \
#     | parallel  --bar \
#                 --env blast \
#                 --env cpu \
#                 --env maxProc \
#                 --env qc \
#                 --jobs "$maxProc" \
#                 "blast {}"

# It's faster to run the blast on long contigs with all cores one sample at the time
# than to parallel blast the samples with fewer cores, at least for nt
c=0
n=$(find "$polished" -type f -name "*.fasta" | wc -l)
for i in $(find "$polished" -type f -name "*.fasta"); do
    let c+=1
    echo -ne "Blasting assembly "${c}"/"${n}" \\r"
    blast "$i"
done


### QUAST ###

declare -a genomes=()
for i in $(find "$polished" -type f -name "*.fasta"); do 
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
find "$polished" -type f -name "*.fasta" \
    | parallel  --bar \
                --env run_quast \
                --env maxProc \
                --env cpu \
                --env qc \
                --env smallest_contig \
                --jobs "$maxProc" \
                "run_quast {}"

source deactivate


###################
#                 #
#   Annotation    #
#                 #
###################


function annotate()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    #Prokka
    prokka  --outdir "${annotation}"/"$sample" \
            --force \
            --prefix "$sample" \
            --kingdom "$kingdom" \
            --genus "$genus" \
            --species "$species" \
            --strain "$sample" \
            --gram "$gram" \
            --locustag "$locustag" \
            --compliant \
            --centre "$centre" \
            --cpus $((cpu/maxProc)) \
            --rfam \
            "$1"

    #extract hypothetical proteins
    cat "${annotation}"/"${sample}"/"${sample}".faa | \
        awk '{if(substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' | \
        grep --no-group-separator -A 1 -F "hypothetical protein" \
        > "${annotation}"/"${sample}"/"${sample}"_hypoth.faa

    echo -e ""$sample" hypothetical proteins (round1): $(cat "${annotation}"/"${sample}"/"${sample}".faa | grep -ic "hypothetical")" \
        | tee -a "${logs}"/log.txt
}

export -f annotate

find "$polished" -type f -name "*.fasta" | \
    parallel    --bar \
                --env annotate \
                --env annotation \
                --env kingdom \
                --env genus \
                --env species \
                --env gram \
                --env locustag \
                --env centre \
                --env cpu \
                --env maxProc \
                --env scripts \
                --jobs "$maxProc" \
                "annotate {}"


### Resfinder

function run_resfinder ()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    [ -d "${amr}"/"${sample}"/"$resfinder_db" ] || mkdir -p "${amr}"/"${sample}"/"$resfinder_db"

    perl "${prog}"/resfinder/resfinder.pl \
        -d "${prog}"/resfinder/resfinder_db/ \
        -a "$resfinder_db" \
        -i "$1" \
        -o "${amr}"/"$sample"/"$resfinder_db" \
        -k 90 \
        -l 60
}

export -f run_resfinder

for h in $(find "${prog}"/resfinder/resfinder_db -type f -name "*.fsa"); do
    export resfinder_db=$(sed 's/\.fsa//' <<< $(basename "$h"))

    find "$polished" -type f -name "*.fasta" | \
        parallel --env run_resfinder \
            --env resfinder_db \
            "run_resfinder {}"
done

#Check if any hit
find "${amr}" -type f -name "results_tab.txt" \
    -exec cat {} \; | sed -n '1d' | tee "${amr}"/resfinder_hits.txt


### Phaster

#trim assemblies
function phaster_trim()
{
    sample=$(cut -d '_' -f 1 <<< $(basename "$1"))

    # http://phaster.ca/instructions
    if [ $(cat "$1" | grep -Ec "^>") -gt 1 ]; then  # if more than one contig
        #remove contigs smaller than 2000 bp from assembly
        perl "${prog}"/phage_typing/removesmallscontigs.pl \
            2000 \
            "$1" \
            > "${phaster}"/assemblies/"${sample}"_trimmed2000.fasta
    elif [ $(cat "$1" | grep -Ec "^>") -eq 1 ]; then  # if only one contig
        #remove contigs smaller than 2000 bp from assembly
        perl "${prog}"/phage_typing/removesmallscontigs.pl \
            1500 \
            "$1" \
            > "${phaster}"/assemblies/"${sample}"_trimmed1500.fasta
    else
        echo "No assembly for "$sample""  # Should not get here!
        # exit 1
    fi
}

#make function available to parallel
export -f phaster_trim  # -f is to export functions

 # To store trimmed assemblies for phaster submission
[ -d "${phaster}"/assemblies ] || mkdir -p "${phaster}"/assemblies 

#run trimming on multiple assemblies in parallel
find "$polished" -type f -name "*.fasta" \
    | parallel  --bar \
                --env phaster_trim \
                --env prog \
                'phaster_trim {}'

# Get phaster results
python3 ~/scripts/checkPhasterServer.py --submit --check \
    -i "${phaster}"/assemblies \
    -o "$phaster"

# extract phaster zip file to access results and rename contigs



#######
# assemble metagenomic data

function assemble()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    [ -d "${assemblies}"/flye_meta/"$sample" ] || mkdir -p "${assemblies}"/flye_meta/"$sample"

    flye \
        --meta \
        --nano-raw "$1" \
        --out-dir "${assemblies}"/flye_meta/"$sample" \
        --genome-size "$size" \
        --iterations 3 \
        --threads $((cpu/maxProc))
}

export -f assemble

find "${trimmed}"/to_assemble -type f -name "*.fastq.gz" \
    | parallel  --bar \
                --env cpu \
                --env maxProc \
                --env assemble \
                --env assembler \
                --env assemblies \
                --env size \
                --jobs "$maxProc" \
                'assemble {} unicycler'

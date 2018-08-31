#!/bin/bash



# Process and assemble Oxford nanopore sequencer data in fast5 format
# Include file format conversion, stats, qc and de novo assembly



#script version
version="0.2.4"


######################
#                    #
#    User Defined    #
#                    #
######################


# Analysis folder
export baseDir=""${HOME}"/analyses/burkholderia_nanopore_batch1"

# Reads
export fast5="/media/2TB_NVMe/burkholderia_fast5_1/20180823_1514_Bmallei/fast5"

# Database to use for metagomic analysis of raw data (contamination)
db="/media/6tb_raid10/db/centrifuge/2017-10-12_bact_vir_h"
# db="/media/6tb_raid10/db/centrifuge/nt"

# Program location
export prog=""${HOME}"/prog"
export scripts=""${HOME}"/scripts"

# Maximum number of cores used per sample for parallel processing
# A highier value reduces the memory footprint.
export maxProc=6

# Estimated genome size in bp
export size=5800000

# Set smallest contig size for assemblies
export smallest_contig=1000

#Annotation
export kingdom="Bacteria"
export genus="Burkholderia"
export species="mallei"
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


# Use albacore

function call_bases_1D()
{
    # Will not do basecalling on ".fast5.tmp" files. They need to be renamed

    kit="$1"  # "SQK-LSK108"
    flow="$2"  # "FLO-MIN106"
    fast5_folder="$3"  # top folder - will look recursively for fast5 files because "-r" option
    fastq_folder="$4"
    barcoded="$5"

    if [[ "$barcoded" == "barcoded" ]]; then
        read_fast5_basecaller.py \
            -i "$fast5_folder" \
            --barcoding \
            -t $((cpu/bins)) \
            -s "$fastq_folder" \
            -k "$kit" \
            -f "$flow" \
            -r \
            -o "fastq" \
            -q 0 \
            --disable_pings
    else
        read_fast5_basecaller.py \
            -i "$fast5_folder" \
            -t $((cpu/bins)) \
            -s "$fastq_folder" \
            -k "$kit" \
            -f "$flow" \
            -r \
            -o "fastq" \
            -q 0 \
            --disable_pings

        #files of most interest are found in "${fastq_folder}"/workplace/pass
    fi
}

export -f call_bases_1D

# Use max 12 threads max per instance for Albacore
export bins=$((cpu/12))  # 4
[ "$bins" -eq 0 ] && bins=1  # in case computer has less than 12 cores
folder_count=$(find "$fast5" -mindepth 1 -type d | wc -l)
folder_per_bin=$(((folder_count/bins)+1))  # the "+1" make sure we never have more than 4 batches.
find "$fast5" -mindepth 1 -type d | parallel -n "$folder_per_bin" 'mkdir -p  {//}/batch{#} && mv {} {//}/batch{#}'

find "$fast5" -mindepth 1 -maxdepth 1 -type d | \
    parallel    --env bins \
                --env cpu \
                --jobs $((cpu/4)) \
                'call_bases_1D "SQK-LSK108" "FLO-MIN106" {} "${basecalled}"/{/} "barcoded"'

### Merge baecalling results ###

# list summary files
sumfiles=($(find "$basecalled" -type f -name "*sequencing_summary.txt"))

#write header
echo -e "filename\tread_id\trun_id\tchannel\tstart_time\tduration\tnum_events\tpasses_filtering\
\ttemplate_start\tnum_events_template\ttemplate_duration\tnum_called_template\tsequence_length_template\
\tmean_qscore_template\tstrand_score_template\tcalibration_strand_genome_template\tcalibration_strand_identity_template\
\tcalibration_strand_accuracy_template\tcalibration_strand_speed_bps_template\tbarcode_arrangement\tbarcode_score\
\tbarcode_full_arrangement\tfront_score\trear_score\tfront_begin_index\tfront_foundseq_length\trear_end_index\
\trear_foundseq_length\tkit\tvariant" \
    > "${basecalled}"/sequencing_summary.txt

#Merge sequencing_summary.txt files, skipping header
for f in "${sumfiles[@]}"; do
    cat "$f" | sed -e '1d' >> "${basecalled}"/sequencing_summary.txt
done

[ -d "${basecalled}"/pass ] || mkdir "${basecalled}"/pass
[ -d "${basecalled}"/fail ] || mkdir "${basecalled}"/fail
state=('pass' 'fail')

#Merge Pass files
# and compress them?
for b in $(find "$basecalled" -mindepth 1 -maxdepth 1 -type d -name "*batch*"); do
    for s in "${state[@]}"; do
        for c in $(find "${b}"/workspace/"${s}" -mindepth 1 -maxdepth 1 -type d); do
            barcode=$(basename "$c")
            # Create output folder
            [ -d "${basecalled}"/"${s}"/"$barcode" ] || mkdir -p "${basecalled}"/"${s}"/"$barcode"
            for f in $(find "$c" -type f -name "*.fastq"); do
                # concaterate files from batches and rename output file accorfing to barcode
                cat "$f" >> "${basecalled}"/"${s}"/"${barcode}"/"${barcode}"_"${s}".fastq
            done
        done
    done
done
            
# Remove batch folder
# find "$basecalled" -maxdepth 1 -type d -name "*batch*" -exec rm -rf {} \;

# compress all fastq files
find "$basecalled" -type f -name "*.fastq" -exec pigz {} \;


##########
#        #
#   QC   #
#        #
##########


# NanoPlot on "sequence=ing_summary.txt" file
[ -d "${qc}"/NanoPlot ] || mkdir -p "${qc}"/NanoPlot
NanoPlot \
    -t "$cpu" \
    -o "${qc}"/NanoPlot \
    -p nanoplot \
    --loglength \
    --summary "${basecalled}"/sequencing_summary.txt \
    --plot kde hex dot

# nanoQC on fastq file
[ -d "${qc}"/nanoQC ] || mkdir -p "${qc}"/nanoQC
python3 "${scripts}"/nanoQC_mmap.py \
    -f "$basecalled" \
    -o "${qc}"/nanoQC


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
        | tee -a "${logs}"/"${sample}"_porechop.log
}

export -f trim

find "$fastq" -type f -name "*.fastq.gz" \
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


function filter()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    # https://github.com/rrwick/Filtlong
    # --min_mean_q 90 -> phred 10
    # --target_bases -> max 100X coverage
    filtlong \
        --target_bases $((size*100)) \
        --keep_percent 90 \
        --min_length 1000 \
        --min_mean_q 90 \
        "$1" \
        2> >(tee "${logs}"/filtering/"${sample}".txt) \
        | pigz > "${filtered}"/"${sample}"_filtered.fastq.gz 
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


######################
#                    #
#   Quality Control  #
#                    #
######################


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

# Raw
find "$fastq"
# Run fastQC
fastqc \
    --o "${qc}"/fastqc/raw \
    --noextract \
    --threads "$cpu" \
    "${fastq}"/"${sample}".fastq.gz

fastqc \
    --o  "${qc}"/fastqc/trimmed \
    --noextract \
    --threads "$cpu" \
    "${fastq}"/"${sample}"_trimmed.fastq.gz


### Centrifuge ###

function run_centrifuge()
{
    r1="$1"
    r2=$(sed 's/_R1/_R2/' <<< "$r1")
    sample=$(cut -d "_" -f 1 <<< $(basename "$r1"))

    [ -d "${2}"/"${sample}" ] || mkdir -p "${2}"/"${sample}"

    #build the command
    centrifuge \
        -p "$cpu" \
        -t \
        --seed "$RANDOM" \
        -x "$centrifuge_db" \
        -1 "$r1" \
        -2 "$r2" \
        --report-file "${2}"/"${sample}"/"${sample}"_report.tsv \
        > "${2}"/"${sample}"/"${sample}".tsv

    cat "${2}"/"${sample}"/"${sample}".tsv | \
        cut -f 1,3 | \
        ktImportTaxonomy /dev/stdin -o "${2}"/"${sample}"/"${sample}".html
}

for i in $(find "$fastq" -type f -name "*_R1*fastq.gz"); do
    run_centrifuge "$i" "${qc}"/centrifuge/raw
done


# Create folder
[ -d "${qc}"/centrifuge ] || mkdir -p "${qc}"/centrifuge

# Run centrifuge
centrifuge \
    -p "$cpu" \
    -t \
    --seed "$RANDOM" \
    -x "$db" \
    -U "${fastq}"/"${sample}"_trimmed.fastq.gz \
    --report-file "${qc}"/centrifuge/"${sample}"_report.tsv \
    > "${qc}"/centrifuge/"${sample}".tsv

# Prepare result for display with Krona
cat "${qc}"/centrifuge/"${sample}".tsv | \
    cut -f 1,3 | \
    ktImportTaxonomy /dev/stdin -o  "${qc}"/centrifuge/"${sample}".html

# Visualize the resutls in Firefow browser
firefox file://"${qc}"/centrifuge/"${sample}".html &


########################
#                      #
#   de novo assembly   #
#                      #
########################


# Assemble with miniasm and polish with racon
for i in $(find "$fastq" -type f -name "*.fastq.gz"); do
    sample=$(cut -d "_" -f 1 <<< $(basename "$i"))
    
    # canu \
    #     -p "$sample" \
    #     -d "${assemblies}"/canu/"$sample" \
    #     genomeSize="$size" \
    #     -nanopore-raw "${fastq}"/"${sample}"_trimmed.fastq.gz

    # flye \
    #     --nano-raw "$i" \
    #     --out-dir "${assemblies}"/flye/"$sample" \
    #     --genome-size "$size" \
    #     --iterations 3 \
    #     --threads "$cpu"

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
            --env basecalled \
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
        --paint-chromosome-limits \
        --output-genome-coverage "${qc}"/coverage/qualimap/"$sample"_coverage \
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
n=$(find "$polished" -type f -name "*_nanopolished.fasta" | wc -l)
for i in $(find "$polished" -type f -name "*_nanopolished.fasta"); do
    let c+=1
    echo -ne "Blasting assembly "${c}"/"${n}" \\r"
    blast "$i"
done


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
            --sample "$sample" \
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

find "$polished" -type f -name "*_nanopolished.fasta" | \
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

    find "$polished" -type f -name "*_nanopolished.fasta" | \
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
find "$polished" -type f -name "*_nanopolished.fasta" \
    | parallel  --bar \
                --env phaster_trim \
                --env prog \
                'phaster_trim {}'


# function phasterSubmit ()
# {
#     sample=$(basename "$1" | cut -d '_' -f 1)

#     # {"job_id":"ZZ_7aed0446a6","status":"You're next!..."}
#     wget --post-file="$i" \
#         "http://phaster.ca/phaster_api?contigs=1" \
#         -O "${phaster}"/"${sample}".json \
#         -o "${phaster}"/"${sample}"_wget.log
# }

# # Submit to phaster sequencially
# c=0
# n=$(find "${phaster}"/assemblies -type f -name "*.fasta" | wc -l)
# for i in $(find "${phaster}"/assemblies -type f -name "*.fasta"); do
#     sample=$(cut -d "_" -f 1 <<< $(basename "$i"))

#     let c+=1
#     echo -ne "Submitting assembly of sample \""${sample}"\" to PHASTER server ("${c}"/"${n}") \\r"

#     phasterSubmit "$i"
#     sleep 10
# done

# Get phaster results
python3 ~/scripts/checkPhasterServer.py --submit --check \
    -i "${phaster}"/assemblies \
    -o "$phaster"

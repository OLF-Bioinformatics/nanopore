#!/bin/bash


##########################   Start of installation procedure - only needed once   ##########################


# Set conda channels order
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels r

# Create all conda envrironments
conda install mamba -n base -c conda-forge -y
mamba install -n base -c conda-forge git -y
mamba create -n nanopore -c bioconda filtlong porechop flye=2.9 pilon parallel bandage blast circlator -y
mamba create -n qualimap -c bioconda qualimap parallel -y
mamba create -n medaka -c conda-forge -c bioconda medaka parallel -y
# mamba create -n medaka -c conda-forge -c bioconda parallel -y && conda activate medaka && pip install medaka
# mamba create -n medaka -c conda-forge -c bioconda python=3.6 medaka parallel -y && conda activate medaka && pip install tensorflow==2.2.2
mamba create -n prokka -c bioconda prokka parallel -y  # the prokka file needs to be changed for it to work. Blast versions. Also, signalp is not included and needed for Gram+
mamba create --name rgi -c conda-forge -c bioconda -c defaults rgi parallel -y
mamba create -n quast -c bioconda quast parallel -y
mamba create -n qualimap -c bioconda qualimap parallel samtools -y
mamba create -n refseq_masher -c bioconda refseq_masher parallel -y && conda activate refseq_masher && pip install importlib-metadata && conda deactivate
mamba create -n resfinder -c bioconda blast kma git parallel numpy -y && conda activate resfinder && pip install tabulate biopython cgecore gitpython python-dateutil && conda deactivate
mamba create -n nanoQC -c bioconda numpy python-dateutil cython pandas matplotlib seaborn parallel git -y


# mamba create -n nanopore git filtlong porechop \
#     flye=2.9 pilon parallel bandage blast circlator qualimap rgi quast qualimap \
#     samtools refseq_masher numpy python-dateutil cython pandas matplotlib seaborn blast kma \
#     cython pandas matplotlib seaborn -y && activate nanopore && pip install medaka && conda deactivate

# Install PGAP

# conda activate nanopore
# mamba install -c conda-forge -c bioconda -c defaults rgi -y
# pip install tabulate biopython cgecore gitpython python-dateutil #tensorflow==2.2.2


# Install resfinder
# https://bitbucket.org/genomicepidemiology/resfinder/src/master/
[ -d $HOME/prog ] || mkdir -p $HOME/prog
cd $HOME/prog
git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git db_resfinder
cd resfinder/db_resfinder
conda activate resfinder
python3 INSTALL.py
conda deactivate

# Install nanoQC
[ -d $HOME/prog ] || mkdir -p $HOME/prog
cd $HOME/prog
git clone https://github.com/duceppemo/nanoQC.git
cd nanoQC
conda activate nanoQC
python setup.py build_ext --inplace
conda deactivate


##########################   End of installation procedure - only needed once   ##########################



#################################   Start of User defined variables   #################################


# Activate required conda environment
conda activate nanopore


# Inputs
export baseDir=/home/bioinfo/analyses/mbovis_nanopore/NP_20220419_AG
export data=/media/36tb/data/Mbovis/nanopore/NP_20220419_AG
bc_desc=/media/36tb/data/Mbovis/nanopore/NP_20220419_AG/bc.txt # Leave empty if not barcodes used
samp=Multiplexed  # Only used to rename sample if no barcodes used

# Guppy config file
nano_config="dna_r10.4_e8.1_sup.cfg"  # dna_r9.4.1_450bps_sup.cfg; dna_r10.3_450bps_sup.cfg; dna_r10.4_e8.1_hac.cfg; dna_r10.4_e8.1_sup.cfg
bc_kit="SQK-NBD112-24"  # "EXP-NBD104" (most common); "EXP-NBD103" (retired); "SQK-NBD112-24" (Q20+)

# Medaka model file
export medaka_model=r103_sup_g507  # r941_min_sup_g507; r941_min_hac_g507; r10_min_high_g340; r103_sup_g507

# Parallel execution
export maxProc=6  # Set according to number of samples  in run and CPU available

# Sample info
export size=4350000 # Pseudomonas syringae: 5630000, Mbovis: 4350000, Lmono: 2850000 # Salmonella: 4850000; Stagonospora: 37000000
export kingdom=bacteria
export genus=Mycobactrium  # Salmonella
export species=bovis  # enterica
export gram=pos  # neg
export locustag=XXX
export centre=OLF

# Assembler [flye, shasta]
export ass=flye  # shasta
export reads_min_len=1000  # only for shasta

# Annotation [prokka, pgap]
export annotator=prokka

#################################   End of User defined variables   #################################


# Performance
export cpu=$(nproc)
export prog=$HOME/prog
export scripts=$HOME/scripts
export mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
export memJava="-Xmx"$mem"g"

# Folder structure
export logs=""${baseDir}"/logs"
export qc=""${baseDir}"/qc"
export trimmed=""${baseDir}"/trimmed"
export filtered=""${baseDir}"/filtered"
export basecalled=""${baseDir}"/basecalled"
export assemblies=""${baseDir}"/assemblies"
export polished_medaka=""${baseDir}"/polished_medaka"
export annotation=""${baseDir}"/annotation"
export phaster=""${baseDir}"/phaster"
export amr=""${baseDir}"/amr"

# Create folders
[ -d "$baseDir" ] || mkdir -p "$baseDir"
[ -d "$logs" ] || mkdir -p "$logs"
[ -d "$qc" ] || mkdir -p "$qc"
[ -d "$trimmed" ] || mkdir -p "$trimmed"
[ -d "$filtered" ] || mkdir -p "$filtered"
[ -d "$basecalled" ] || mkdir -p "$basecalled"
[ -d "$assemblies" ] || mkdir -p "$assemblies"
[ -d "$polished_medaka" ] || mkdir -p "$polished_medaka"
[ -d "$annotation" ] || mkdir -p "$annotation"
[ -d "$amr" ] || mkdir -p "$amr"


##############
#
# Log
#
##############


# Log software versions
echo -e "$(date)\n" | tee "${logs}"/log.txt
echo -e "User: $(whoami)" | tee -a "${logs}"/log.txt
echo -e "Processors: "$cpu"" | tee -a "${logs}"/log.txt
echo -e "Memory: "$mem"G\n" | tee -a "${logs}"/log.txt

java -version 2>&1 1>/dev/null | grep "java version" | tr -d '"' | tee -a "${logs}"/log.txt
v=$(guppy_basecaller --version | grep "Version" | cut -d "+" -f 1 | rev | cut -d " " -f 1 | rev) && echo "Guppy v"${v}"" | tee -a "${logs}"/log.txt
v=$(porechop --version) && echo "Porechop v"${v}"" | tee -a "${logs}"/log.txt
filtlong --version | tee -a "${logs}"/log.txt
v=$(flye --version) && echo "flye v"${v}"" | tee -a "${logs}"/log.txt
v=$(circlator version) && echo "circlator v"${v}"" | tee -a "${logs}"/log.txt
conda activate medaka && medaka --version | tee -a "${logs}"/log.txt && conda deactivate
conda activate quast && quast --version | tee -a "${logs}"/log.txt && conda deactivate
conda activate prokka && prokka --version | tee -a "${logs}"/log.txt && conda deactivate
v=$(minimap2 --version 2>&1) && echo "minimap2 "$v"" | tee -a "${logs}"/log.txt
samtools --version | grep -F 'samtools' | tee -a "${logs}"/log.txt
conda activate qualimap && qualimap --version | grep -F "QualiMap" | tee -a "${logs}"/log.txt && conda deactivate
conda activate refseq_masher && refseq_masher --version  | tee -a "${logs}"/log.txt && conda deactivate
echo "Resfinder v4.1.5" | tee -a "${logs}"/log.txt
echo "nanoQC v0.2" | tee -a "${logs}"/log.txt
python "${prog}"/pgap.py --version | tee -a "${logs}"/log.txt


##############
#
# Basecalling
#
##############


# Basecall and rename samples
if [ -z "$bc_desc" ]; then  # if no barcodes
    # Guppy Super accuracy basecalling
    # EXP-NBD104
    # EXP-NBD114
    # EXP-PBC001
    # RTX A6000

    # configs
    # dna_r9.4.1_450bps_sup.cfg
    # dna_r10.3_450bps_sup.cfg

    # Requires Guppy v6+

    # --detect_primer
    # --trim_primers
    # --detect_mid_strand_adapter
    guppy_basecaller \
        --input_path $data \
        --save_path $basecalled \
        --config "$nano_config" \
        --recursive \
        --device 'cuda:0' \
        --num_callers 4 \
        --gpu_runners_per_device 2 \
        --chunk_size 1000 \
        --chunks_per_runner 512 \
        --disable_pings \
        --calib_detect \
        --compress_fastq \
        --detect_adapter \
        --trim_adapters \

    # remove calibration strands, if presents
    [ -d "${basecalled}"/calibration_strands ] && rm -rf "${basecalled}"/calibration_strands

    # Merge the basecalled and demultiplexed fastq
    for i in $(find "$basecalled" -mindepth 1 -maxdepth 1 -type d); do  # pass and fail
        flag=$(basename "$i")

        find "$i" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -name "fastq_runid_*" \
            -exec cat {} \; > "${i}"/"${samp}"_"${flag}".fastq.gz
        # # Remove non compressed files
        find "$i" -type f -name "*fastq_runid_*" -exec rm {} \;
    done
else  # if barcodes
    # Guppy Super accuracy basecalling
    # EXP-NBD103 (older)
    # EXP-NBD104
    # EXP-NBD114
    # EXP-PBC001
    # RTX A6000

    # configs
    # dna_r9.4.1_450bps_sup.cfg
    # dna_r10.3_450bps_sup.cfg

    # Requires Guppy v6+

    # --detect_primer
    # --trim_primers
    # --detect_mid_strand_adapter
    guppy_basecaller \
        --input_path $data \
        --save_path $basecalled \
        --config "$nano_config" \
        --recursive \
        --device 'cuda:0' \
        --num_callers 4 \
        --gpu_runners_per_device 2 \
        --chunk_size 1000 \
        --chunks_per_runner 512 \
        --disable_pings \
        --calib_detect \
        --compress_fastq \
        --detect_adapter \
        --trim_adapters \
        --detect_barcodes \
        --trim_barcodes \
        --barcode_kits "$bc_kit"

    # Only hac available with this workflow for old R10 flowcells
    # guppy_basecaller \
    #     --input_path $data \
    #     --save_path $basecalled \
    #     --flowcell "FLO-MIN110" \
    #     --kit "SQK-LSK109" \
    #     --recursive \
    #     --device 'cuda:0' \
    #     --num_callers 4 \
    #     --gpu_runners_per_device 4 \
    #     --chunk_size 3000 \
    #     --chunks_per_runner 1000 \
    #     --disable_pings \
    #     --calib_detect \
    #     --compress_fastq \
    #     --detect_adapter \
    #     --trim_adapters \
    #     --detect_barcodes \
    #     --trim_barcodes \
    #     --barcode_kits "EXP-NBD104"

    # remove calibration strands, if presents
    [ -d "${basecalled}"/calibration_strands ] && rm -rf "${basecalled}"/calibration_strands

    # Merge the basecalled and demultiplexed fastq
    for i in $(find "$basecalled" -mindepth 2 -maxdepth 2 -type d); do  # pass and fail
        barcode=$(basename "$i")
        flag=$(basename $(dirname "$i"))

        find "$i" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -name "fastq_runid_*" \
            -exec cat {} \; > "${i}"/"${barcode}"_"${flag}".fastq.gz
        # # Remove non compressed files
        find "$i" -type f -name "*fastq_runid_*" -exec rm {} \;
    done

    # Rename files
    declare -A myArray=()

    while IFS= read -r line || [[ -n "$line" ]]; do
        line="$(sed -e 's/ /\t/' -e 's/\r//g' -e 's/\n//g' <<< "$line")" #transform the space output field separator from read into tabs, remove carriage return
        key="$(cut -f 1 <<< "$line")"
        value="$(cut -f 2 <<< "$line")"
        myArray["${key}"]="${value}"
    done < "$bc_desc"

    find "${basecalled}"/ -type f -name "*_pass*" -o -name "*_fail*" | while read i; do
        # echo "$i"
        pathPart="$(dirname "$i")"
        # echo "$pathPart"
        oldName="$(basename "$i")"
        # echo "$oldName"

        #for each file, check if a part of the name matches on
        for j in "${!myArray[@]}"; do
            # echo "$j"
            if [ "$(echo "$oldName" | grep "$j")" ]; then
                newName="$(echo "$oldName" | sed "s/"$j"/"${myArray["$j"]}"/")"
                fullNewName=""${pathPart}"/"${newName}""

                if [ -e "$rename" ]; then
                    echo "Cannot rename "$oldName" to "$newName", file already exists. Skipping"
                    continue
                    # exit 1
                fi

                echo ""$i" -> "$fullNewName""
                mv "$i" "$fullNewName"
            fi
            if [ "$(echo "$pathPart" | grep "$j")" ]; then
                # rename folder too
                echo ""$pathPart" -> $(dirname "$pathPart")/"${myArray["$j"]}""
                mv $pathPart $(dirname "$pathPart")/"${myArray["$j"]}"
            fi
        done
    done

    # Remove samples starting by "barcode", i.e. not renamed because not to added barcode
    find "$basecalled" -type d -name "*barcode*" -exec rm -rf {} \;
fi



# QC on raw reads
[ -d "${qc}"/nanoQC/raw/fastq ] || mkdir -p "${qc}"/nanoQC/raw/fastq
conda activate nanoQC
python3 "${prog}"/nanoQC/nanoQC.py \
    -f "$basecalled" \
    -o "${qc}"/nanoQC/raw/fastq
conda deactivate


##############
#
# Read trimming
#
##############


# Trim sequencing adapters
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


##############
#
# Read filtering
#
##############


# Filter out low quality reads
function filter()
{
    #sample=$(cut -d "_" -f 1 <<< $(basename "$1"))
    sample=$(basename "$1" "_trimmed.fastq.gz")

    filtlong \
        --keep_percent 95 \
        --target_bases $(($size*100)) \
        "$1" \
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


# QC on filtered reads
[ -d "${qc}"/nanoQC/raw/fastq ] || mkdir -p "${qc}"/nanoQC/raw/fastq
conda activate nanoQC
python3 "${prog}"/nanoQC/nanoQC.py \
    -f "$filtered" \
    -o "${qc}"/nanoQC/filtered/fastq
conda deactivate


##############
#
# Assembly
#
##############


function assemble_flye()
{
    sample=$(basename "${1%_filtered.fastq.gz}")

    [ -d "${assemblies}"/"${ass}"/"$sample" ] || mkdir -p "${assemblies}"/"${ass}"/"$sample"
    [ -d "${qc}"/assembly_graphs/"${ass}" ] || mkdir -p "${qc}"/assembly_graphs/"${ass}"

    flye \
        --nano-hq "$1" \
        --out-dir "${assemblies}"/"${ass}"/"$sample" \
        --genome-size "$size" \
        --iterations 3 \
        --threads $((cpu/maxProc))

    # Rename assembly
    mv "${assemblies}"/"${ass}"/"${sample}"/assembly.fasta \
        "${assemblies}"/"${ass}"/"${sample}"/"${sample}".fasta

    mv "${assemblies}"/"${ass}"/"${sample}"/assembly_graph.gfa \
        "${assemblies}"/"${ass}"/"${sample}"/"${sample}".gfa

    Bandage image \
        "${assemblies}"/"${ass}"/"${sample}"/"${sample}".gfa \
        "${qc}"/assembly_graphs/"${ass}"/"${sample}".png

    # Maybe some cleanup
    find "${assemblies}"/"${ass}"/"$sample" -mindepth 1 -maxdepth 1 \
        ! -name "*.fasta" \
        ! -name "*.gfa" \
        ! -name "*.log" \
        ! -name "*.txt" \
        -exec rm -rf {} \;
}

export -f assemble_flye

function assemble_shasta()
{
    sample=$(basename "${1%_filtered.fastq.gz}")

    [ -d "${qc}"/assembly_graphs/"${ass}" ] || mkdir -p "${qc}"/assembly_graphs/"${ass}"
    [ -d "${assemblies}"/"$ass" ] || mkdir -p "${assemblies}"/"$ass"

    pigz -dk "$1"

    shasta \
        --input "${1%.gz}" \
        --config Nanopore-Oct2021 \
        --assemblyDirectory "${assemblies}"/"${ass}"/"$sample" \
        --Reads.minReadLength "$reads_min_len" \
        --threads $((cpu/maxProc))

    shasta \
        --command cleanupBinaryData \
        --assemblyDirectory "${assemblies}"/"${ass}"/"$sample"

    # Rename output files
    mv "${assemblies}"/"${ass}"/"${sample}"/Assembly.fasta \
        "${assemblies}"/"${ass}"/"${sample}"/"${sample}".fasta

    mv "${assemblies}"/"${ass}"/"${sample}"/Assembly.gfa \
        "${assemblies}"/"${ass}"/"${sample}"/"${sample}".gfa

    # Display assembly graphs
    Bandage image \
        "${assemblies}"/"${ass}"/"${sample}"/"${sample}".gfa \
        "${qc}"/assembly_graphs/"${ass}"/"${sample}".png

    # Delete uncompress fastq
    rm "${1%.gz}"
}

export -f assemble_shasta


find -L "${filtered}" -type f -name "*.fastq.gz" \
    | parallel  --bar \
                --env cpu \
                --env maxProc \
                --env assemble_${ass} \
                --env ass \
                --env assemblies \
                --env size \
                --env qc \
                --env reads_min_len \
                --env reads_min_overlap \
                --jobs "$maxProc" \
                "assemble_${ass} {}"


##############
#
# Polishing
#
##############


# Polish with Medaka
function polish_medaka()
{
    sample=$(basename "$1" ".fasta")

    [ -d "${polished_medaka}"/"${ass}"/"$sample" ] || mkdir -p "${polished_medaka}"/"${ass}"/"$sample"

    # It is crucially important to specify the correct model (-m) according to the basecaller used.
    # Allowed values can be found by running "medaka tools list\_models"
    # -i "${basecalled}"/pass/"${sample}"/"${sample}"_pass.fastq.gz \

    # -h  show this help text.
    # -i  fastx input basecalls (required).
    # -d  fasta input assembly (required).
    # -o  output folder (default: medaka).
    # -g  don't fill gaps in consensus with draft sequence.
    # -m  medaka model, (default: ${MODEL}).
    #     ${modeldata[0]}.
    #     Alternatively a .hdf file from 'medaka train'.
    # -f  Force overwrite of outputs (default will reuse existing outputs).
    # -t  number of threads with which to create features (default: 1).
    # -b  batchsize, controls memory use (default: ${BATCH_SIZE})."

    medaka_consensus \
        -i "${filtered}"/"${sample}"_filtered.fastq.gz \
        -d "$1" \
        -o "${polished_medaka}"/"${ass}"/"$sample" \
        -t $((cpu/maxProc)) \
        -m "$medaka_model" \
        2>&1 | tee "${polished_medaka}"/"${ass}"/"${sample}"/medaka.log

    # Rename sample
    mv "${polished_medaka}"/"${ass}"/"${sample}"/consensus.fasta \
        "${polished_medaka}"/"${ass}"/"${sample}"/"${sample}"_medaka.fasta

    # Remove temporary files
    find "${polished_medaka}"/"${ass}"/"${sample}" -type f \
        ! -name "*_medaka.fasta" \
        ! -name "*.gfa" \
        ! -name "*.log" \
        -exec rm {} \;

    # Reformat fasta
    python "${scripts}"/format_fasta.py \
        "${polished_medaka}"/"${ass}"/"${sample}"/"${sample}"_medaka.fasta \
        "${polished_medaka}"/"${ass}"/"${sample}"/"${sample}"_medaka.fasta1
    mv "${polished_medaka}"/"${ass}"/"${sample}"/"${sample}"_medaka.fasta1 \
        "${polished_medaka}"/"${ass}"/"${sample}"/"${sample}"_medaka.fasta
}

export -f polish_medaka

conda activate medaka

find "${assemblies}"/"$ass" -maxdepth 2 -type f -name "*.fasta" \
    | parallel  --bar \
                --env polish_medaka \
                --env polished_medaka \
                --env scripts \
                --env cpu \
                --env ass \
                --env maxProc \
                'polish_medaka {}'

conda deactivate




##############
#
# Polish with short reads if available
#  export maxProc=1
#
##############



##############
#
# Orienting
#
##############


# fixstart
function fix_start()
{
    sample=$(basename "$1" ".fasta" | cut -d "_" -f 1)

    circlator fixstart \
        --verbose \
        "$1" \
        "${baseDir}"/fixstart/"$ass"/"$sample"

    # check if dnaA was found
    # if not, dont keep the new "fixed" assebmly because it has be rerooted from the middle
}

export -f fix_start

[ -d "${baseDir}"/fixstart/"$ass" ] || mkdir -p "${baseDir}"/fixstart/"$ass"
cd "${baseDir}"/fixstart/"$ass"

find "${polished_medaka}"/"$ass" -type f -name "*.fasta" \
    | parallel  --bar \
                --env fix_start \
                --env maxProc \
                --env cpu \
                --env qc \
                --env ass \
                --jobs "$maxProc" \
                "fix_start {}"

# Cleanup
find "${baseDir}"/fixstart/"$ass" -type f ! -name "*.fasta" ! -name "*.log" -exec rm {} \;



# If no orientation done (say not a bacterial genome)
# [ -d "${baseDir}"/fixstart/"$ass" ] || mkdir -p "${baseDir}"/fixstart/"$ass"
# find "${polished_medaka}"/"$ass"  -type f -name "*.fasta" \
#     -exec sh -c 'ln -s {} "$baseDir"/fixstart/"$ass"/$(basename {} "_medaka.fasta").fasta' \;



##############
#
# Assembly QC
#
##############


# Quast
conda activate quast

# Create list of all assemblies
declare -a genomes=()
for i in $(find -L "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta"); do 
    genomes+=("$i")
done

# Run Quast on all samples to compare them
quast.py \
    --output-dir "${qc}"/quast/"${ass}"/all \
    --threads "$cpu" \
    --min-contig 1000 \
    --est-ref-size "$size" \
    ${genomes[@]}


# Run quast on individual assembly
function run_quast()
{
    sample=$(basename "$1" ".fasta")

    quast.py \
        -t $((cpu/maxProc)) \
        -o "${qc}"/quast/"${ass}"/"$sample" \
        --min-contig 1000 \
        --est-ref-size "$size" \
        $1  # don't put in quotes
}

# make function available to parallel
export -f run_quast  # -f is to export functions

find -L "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta" \
    | parallel  --bar \
                --env run_quast \
                --env maxProc \
                --env cpu \
                --env qc \
                --env ass \
                --jobs "$maxProc" \
                "run_quast {}"

conda deactivate


# Coverage

function get_coverage()  # unsing unmerged reads only
{
    # sample=$(basename "$1" ".fasta" | cut -d '_' -f 1)
    sample=$(basename "$1" ".fasta")
    fastq_file="${trimmed}"/"${sample}"_trimmed.fastq.gz

    [ -d "${qc}"/coverage/"${ass}"/"$sample" ] || mkdir -p "${qc}"/coverage/"${ass}"/"$sample"

    minimap2 \
        -ax map-ont \
        -t $((cpu/maxProc)) \
        "$1" \
        "$fastq_file" | \
    samtools view -@ $((cpu/maxProc)) -b -h -F 4 - | \
    samtools sort -@ $((cpu/maxProc)) - | \
    samtools rmdup - "${qc}"/coverage/"${ass}"/"${sample}"/"${sample}".bam

    # Index bam file
    samtools index -@ $((cpu/maxProc)) "${qc}"/coverage/"${ass}"/"${sample}"/"${sample}".bam

    #Average genome depth of coverage
    average_cov=$(samtools depth \
        "${qc}"/coverage/"${ass}"/"${sample}"/"${sample}".bam  \
        | awk '{sum+=$3} END { print sum/NR}')

    printf "%s\t%.*f\n" "$sample" 0 "$average_cov" | tee -a "${qc}"/coverage/"${ass}"/average_cov.tsv
}

export -f get_coverage

[ -d "${qc}"/coverage/"${ass}" ] || mkdir -p "${qc}"/coverage/"${ass}"
echo -e "Sample\tAverage_Cov" > "${qc}"/coverage/"${ass}"/average_cov.tsv

find -L "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta" | \
    parallel    --bar \
                --env get_coverage \
                --env cpu \
                --env maxProc \
                --env ass \
                --env qc \
                --jobs "$maxProc" \
                "get_coverage {}"


# Qualimap

conda activate qualimap

function run_qualimap()
{
    sample=$(basename "$1" ".bam")
    
    [ -d "${qc}"/coverage/"${ass}"/qualimap/"$sample" ] || mkdir -p "${qc}"/coverage/"${ass}"/qualimap/"$sample"

    qualimap bamqc \
        --paint-chromosome-limits \
        -bam "$1" \
        --java-mem-size="${mem}"G \
        -nt $((cpu/maxProc)) \
        -outdir "${qc}"/coverage/"${ass}"/qualimap/"$sample" \
        -outfile "${sample}" \
        -outformat HTML

    # Remove bam files
    rm -rf "${qc}"/coverage/"${ass}"/"$sample"
}

export -f run_qualimap

find "${qc}"/coverage/$ass -type f -name "*.bam" |
parallel    --bar \
            --env run_qualimap \
            --env qc \
            --env mem \
            --env cpu \
            --env ass \
            --env maxProc \
            --jobs "$maxProc" \
            'run_qualimap {}'

conda deactivate


# ID samples
conda activate refseq_masher

function id_sample()
{
    sample=$(basename "$1" ".fasta")

    #         --output "${qc}"/id/"${sample}".tsv \
    refseq_masher matches \
        --output-type "tab" \
        -n 1 \
        -m 1 \
        $1 \
    | sed -n 2p >> "${qc}"/id/refseq_masher_ids.tsv
}

export -f id_sample

echo -e "sample\ttop_taxonomy_name\tdistance\tpvalue\tmatching\tfull_taxonomy\ttaxonomic_species\ttaxonomic_genus\ttaxonomic_family\ttaxonomic_order\ttaxonomic_class\ttaxonomic_phylum\ttaxonomic_superkingdom\tsubspecies\tserovar\tplasmid\tbioproject\tbiosample\ttaxid\tassembly_accession\tmatch_id" \
    > "${qc}"/id/refseq_masher_ids.tsv


[ -d "${qc}"/id ] || mkdir -p "${qc}"/id

find -L "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta" \
    | parallel    --bar \
                --env id_sample \
                --env qc \
                "id_sample {}"

conda deactivate


##############
#
# Annotation
#
##############


conda activate prokka


if ("$annotator" == 'pgap'); then

    # Annotate with PGAP
    # v2021-11-29.build5742
    # Create input yaml file
    [ -d "${annotation}"/pgap ] || mkdir -p "${annotation}"/pgap

    function annotate_pgap()
    {
        sample=$(basename "$1" ".fasta")
        file_name=$(basename "$1")

        input_dir=$(dirname "$1")
        out_dir="${annotation}"/pgap/"$sample"

        input_yaml="${sample}"_inputs.yaml
        metadata_yaml="${sample}"_metadata.yaml

        cd "${annotation}"/pgap

        # 3 files (1 fasta and 2 yaml) need to be in same folder
        cp "$1" "${annotation}"/pgap

        # Create yaml file for input files
        echo -n "\
fasta: 
    class: File
    location: "$file_name"
submol:
    class: File
    location: "$metadata_yaml"
" > "$input_yaml"

    # Create yaml file for sample metadata
    echo -n "\
organism:
    genus_species: '"$genus" "$species"'
    strain: 'my_strain'
" > "$metadata_yaml"

    
        # Run PGAP
        python "${prog}"/pgap.py \
            --no-self-update \
            --report-usage-false \
            --cpus "$cpu" \
            --memory "${mem}"g \
            --output "$out_dir" \
            "$input_yaml"

        # Cleanup?
        # rm "${annotation}"/pgap/"$file_name"
    }

    export -f annotate_pgap

    for i in $(find -L "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta" ); do
        annotate_pgap "$i"
    done
else
    # Annotate assemblies with Prokka
    function annotate_prokka()
    {
        sample=$(basename "$1" ".fasta")

        prokka  \
            --outdir "${annotation}"/"${ass}"/"$sample" \
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
        cat "${annotation}"/"${ass}"/"${sample}"/"${sample}".faa | \
            awk '{if(substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' | \
            grep --no-group-separator -A 1 -F "hypothetical protein" \
            > "${annotation}"/"${ass}"/"${sample}"/"${sample}"_hypoth.faa

        echo -e ""$sample" hypothetical proteins (round1): $(cat "${annotation}"/"${ass}"/"${sample}"/"${sample}".faa | grep -ic "hypothetical")" \
            | tee -a "${logs}"/log.txt
    }

    export -f annotate_prokka

    find -L "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta" | \
        parallel    --bar \
                    --env annotate_prokka \
                    --env annotation \
                    --env kingdom \
                    --env genus \
                    --env species \
                    --env gram \
                    --env locustag \
                    --env centre \
                    --env cpu \
                    --env ass \
                    --env maxProc \
                    --env ass \
                    --env scripts \
                    --jobs "$maxProc" \
                    "annotate_prokka {}"

    conda deactivate
fi


##############
#
# AMR
#
##############


### Resfinder
conda activate resfinder

function run_resfinder ()
{
    # https://bitbucket.org/genomicepidemiology/resfinder/overview
    sample=$(basename "$1" ".fasta")

    [ -d "${amr}"/resfinder/"${ass}"/"$sample" ] || mkdir -p "${amr}"/resfinder/"${ass}"/"$sample"

    python3 "${prog}"/resfinder/run_resfinder.py \
        -ifa "$1" \
        -o "${amr}"/resfinder/"${ass}"/"$sample" \
        -db_res_kma "${prog}"/resfinder/resfinder_db/ \
        -t 0.9 \
        -l 0.6 \
        -acq
        # 1> >(tee "${amr}"/resfinder/"${sample}"/"${sample}"_resfinder.txt)

    rm -rf "${amr}"/resfinder/"${ass}"/"${sample}"/tmp
}

export -f run_resfinder

find -L "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta" | \
    parallel    --env run_resfinder \
                --env resfinder_db \
                --env ass \
                --env prog \
                "run_resfinder {}"

# Create merged report
echo -e 'Sample\tResistance gene\tIdentity\tAlignment Length/Gene Length\tCoverage\tPosition in reference\tContig\tPosition in contig\tPhenotype\tAccession no.' \
    > "${amr}"/resfinder/"${ass}"/resfinder_merged.tsv.tmp
    
for i in $(find "${amr}"/resfinder/"$ass" -name "*results_tab.txt"); do
    # sample name is folder name
    sample=$(basename $(dirname "$i"))
    sample="${sample%.fasta}"

    # Add a leading column with sample name
    cat "$i" \
        | sed -e '1d' \
        | awk -F $'\t' -v s="$sample" 'BEGIN {OFS = FS} {print s,$0}' \
        >> "${amr}"/resfinder/"${ass}"/resfinder_merged.tsv.tmp
done

# sort by sample name (column 1), then by Identity (column 3)
(cat "${amr}"/resfinder/"${ass}"/resfinder_merged.tsv.tmp | head -n 1;
    cat "${amr}"/resfinder/"${ass}"/resfinder_merged.tsv.tmp | sed -e '1d' | sort -t $'\t' -k1,1 -k3,3) \
    > "${amr}"/resfinder/"${ass}"/resfinder_merged.tsv

rm "${amr}"/resfinder/"${ass}"/resfinder_merged.tsv.tmp

conda deactivate


### RGI (CARD)
conda activate rgi

function run_rgi()
{
    sample=$(basename "$1" ".fasta")

    [ -d "${amr}"/rgi/"${ass}"/"$sample" ] || mkdir -p "${amr}"/rgi/"${ass}"/"$sample"

    rgi main \
        -i "$1" \
        -o  "${amr}"/rgi/"${ass}"/"${sample}"/"$sample" \
        -t 'contig' \
        -a 'BLAST' \
        -n $((cpu/maxProc)) \
        --clean \
        -d "chromosome"
}

export -f run_rgi

find -L "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta" \
    | parallel    --bar \
                --env run_rgi \
                --env resfinder_db \
                --env ass \
                "run_rgi {}"

# Create merged report
echo -e 'Sample\tORF_ID\tContig\tStart\tStop\tOrientation\tCut_Off\tPass_Bitscore\tBest_Hit_Bitscore\tBest_Hit_ARO\tBest_Identities\tARO\tModel_typeSNPs_in_Best_Hit_ARO\tOther_SNPs\tDrug\tClass\tResistance\tMechanism\tAMR\tGene\tFamily\tPredicted_DNA\tPredicted_Protein\tCARD_Protein_Sequence\tPercentage_Length_of_Reference_Sequence\tID\tModel_ID' \
    > "${amr}"/rgi/"${ass}"/rgi_merged.tsv.tmp
    
for i in $(find "${amr}"/rgi/"${ass}" -name "*.txt"); do
    # sample name is folder name
    sample=$(basename $(dirname "$i"))
    sample="${sample%.fasta}"

    # Add a leading column with sample name
    cat "$i" \
        | sed -e '1d' \
        | awk -F $'\t' -v s="$sample" 'BEGIN {OFS = FS} {print s,$0}' \
        >> "${amr}"/rgi/"${ass}"/rgi_merged.tsv.tmp
done

# sort by Sample name (column 1), then by Cutt_off (column 7)
(cat "${amr}"/rgi/"${ass}"/rgi_merged.tsv.tmp | head -n 1;
    cat "${amr}"/rgi/"${ass}"/rgi_merged.tsv.tmp | sed -e '1d' | sort -t $'\t' -k1,1 -k7,7) \
    > "${amr}"/rgi/"${ass}"/rgi_merged.tsv

rm "${amr}"/rgi/"${ass}"/rgi_merged.tsv.tmp

conda deactivate

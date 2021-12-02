conda install mamba -n base -c conda-forge
mamba install -c conda-forge git
mamba create -n nanopore -c bioconda filtlong porechop flye pilon parallel bandage blast circlator -y
mamba create -n qualimap -c bioconda qualimap parallel -y
mamba create -n medaka -c conda-forge -c bioconda python=3.6 medaka parallel -y && conda activate medaka && pip install tensorflow==2.2.2
mamba create -n prokka -c bioconda prokka parallel  # the prokka file needs to be changed for it to work. Blast versions. Also, signalp is not included and needed for Gram+
mamba create --name rgi --channel conda-forge --channel bioconda --channel defaults rgi parallel -y
mamba create -n quast -c bioconda quast parallel -y
mamba create -n qualimap -c bioconda qualimap parallel samtools -y
mamba create -n refseq_masher -c bioconda refseq_masher parallel -y &&  pip install importlib-metadata

# Install nanoQC
mamba create -n nanoqc -c bioconda numpy python-dateutil cython pandas matplotlib seaborn parallel git -y
conda activte nanoqc
[ -d $HOME/prog ] || mkdir -p $HOME/prog
cd $HOME/prog
git clone https://github.com/duceppemo/nanoQC.git
cd nanoQC
python setup.py build_ext --inplace

# Install Guppy
sudo apt-get update
sudo apt-get install wget lsb-release
export PLATFORM=$(lsb_release -cs)
wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
echo "deb http://mirror.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
sudo apt-get update
sudo apt install ont-guppy 

# Install resfinder
# https://bitbucket.org/genomicepidemiology/resfinder/src/master/
[ -d $HOME/prog ] || mkdir -p $HOME/prog
cd $HOME/prog
git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git
mamba create -n resfinder -c bioconda blast kma git parallel numpy
conda activate resfinder
pip3 install tabulate biopython cgecore gitpython python-dateutil
# cd resfinder/cge
# git clone https://bitbucket.org/genomicepidemiology/kma.git
# cd kma && make
# cd ..
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git db_resfinder
cd resfinder/db_resfinder
python3 INSTALL.py


conda activate nanopore

# User defined
export baseDir=/home/bioinfo/anaylses/23-listeria_nanopore
export data=/media/bioinfo/DATA1/data/23-listeria_nanopore/fast5
export size=2850000
export maxProc=6
export kingdom=bacteria
export genus=Listeria
export species=monocytogenes
export gram=pos
export locustag=XXX
export centre=OLF


export cpu=$(nproc)
export ass=flye
export prog=$HOME/prog
export scripts=$HOME/scripts
export mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
export memJava="-Xmx"$mem"g"

export logs=""${baseDir}"/logs"
export qc=""${baseDir}"/qc"
export fastq=""${baseDir}"/fastq"
export trimmed=""${baseDir}"/trimmed"
export filtered=""${baseDir}"/filtered"
export basecalled=""${baseDir}"/basecalled_sup"
export assemblies=""${baseDir}"/assemblies"
export polished=""${baseDir}"/polished"
export annotation=""${baseDir}"/annotation"
export phaster=""${baseDir}"/phaster"
export amr=""${baseDir}"/amr"


[ -d "$baseDir" ] || mkdir -p "$baseDir"
[ -d "$logs" ] || mkdir -p "$logs"
[ -d "$qc" ] || mkdir -p "$qc"
[ -d "$fastq" ] || mkdir -p "$fastq"
[ -d "$trimmed" ] || mkdir -p "$trimmed"
[ -d "$filtered" ] || mkdir -p "$filtered"
[ -d "$basecalled" ] || mkdir -p "$basecalled"
[ -d "$assemblies" ] || mkdir -p "$assemblies"
[ -d "$polished" ] || mkdir -p "$polished"
[ -d "$annotation" ] || mkdir -p "$annotation"
[ -d "$phaster" ] || mkdir -p "$phaster"
[ -d "$amr" ] || mkdir -p "$amr"


# Version logging
echo -e "$(date)\n" | tee "${logs}"/log.txt
echo -e "User: $(whoami)" | tee -a "${logs}"/log.txt
echo -e "Processors: "$cpu"" | tee -a "${logs}"/log.txt
echo -e "Memory: "$mem"G" | tee -a "${logs}"/log.txt

java -version 2>&1 1>/dev/null | grep "java version" | tr -d '"' | tee -a "${logs}"/log.txt
v=$(guppy_basecaller --version | rev | cut -d " " -f 1 | rev) && echo "Guppy v"${v}"" | tee -a "${logs}"/log.txt
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
# Resfinder v4.1.5
# nanoQC v0.2


# Guppy Super accuracy basecalling
# EXP-NBD114
guppy_basecaller \
    -i $data \
    -s $basecalled \
    -c dna_r9.4.1_450bps_sup.cfg \
    --recursive \
    -x 'cuda:0' \
    --num_callers 4 \
    --gpu_runners_per_device 2 \
    --chunk_size 1000 \
    --chunks_per_runner 128 \
    --disable_pings \
    --calib_detect \
    --compress_fastq \
    --barcode_kits "EXP-NBD104" \
    --trim_barcodes \
    --num_barcode_threads 48


# Merge the basecalled and demultiplexed fastq
for i in $(find "$basecalled" -mindepth 2 -maxdepth 2 -type d); do  # pass and fail
    barcode=$(basename "$i")
    flag=$(basename $(dirname "$i"))
    # echo "${i}"/"${barcode}"_"${flag}".fastq.gz
    # cat "${i}"/*.fastq.gz > "${i}"/"${barcode}"_"${flag}".fastq.gz  #bash: /bin/cat: Argument list too long
    find "$i" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -name "fastq_runid_*" \
        -exec cat {} \; > "${i}"/"${barcode}"_"${flag}".fastq.gz
    # # Remove non compressed files
    find "$i" -type f -name "*fastq_runid_*" -exec rm {} \;
done


# QC on raw reads
[ -d "${qc}"/nanoQC/raw/fastq ] || mkdir -p "${qc}"/nanoQC/raw/fastq
conda activate nanoqc
python3 "${prog}"/nanoQC/nanoQC.py \
    -f "$basecalled" \
    -o "${qc}"/nanoQC/raw/fastq
conda deactivate

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


# Filter out low quality reads
function filter()
{
    #sample=$(cut -d "_" -f 1 <<< $(basename "$1"))
    sample=$(basename "$1" "_trimmed.fastq.gz")

    filtlong \
        --keep_percent 95 \
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
conda activate nanoqc
python3 "${prog}"/nanoQC/nanoQC.py \
    -f "$filtered" \
    -o "${qc}"/nanoQC/filtered/fastq
conda deactivate


# Assemble with Flye
function assemble_flye()
{
    sample=$(basename "${1%_filtered.fastq.gz}")

    [ -d "${assemblies}"/"${ass}"/"$sample" ] || mkdir -p "${assemblies}"/"${ass}"/"$sample"
    [ -d "${qc}"/assembly_graphs/"${ass}" ] || mkdir -p "${qc}"/assembly_graphs/"${ass}"

    flye \
        --nano-raw "$1" \
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

find -L "${filtered}" -type f -name "*.fastq.gz" \
    | parallel  --bar \
                --env cpu \
                --env maxProc \
                --env assemble_flye \
                --env ass \
                --env assemblies \
                --env size \
                --env qc \
                --jobs "$maxProc" \
                "assemble_flye {}"


# Polish with Medaka
function polish_medaka()
{
    sample=$(basename "$1" ".fasta")

    [ -d "${polished}"/"${ass}"/"$sample" ] || mkdir -p "${polished}"/"${ass}"/"$sample"

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
        -o "${polished}"/"${ass}"/"$sample" \
        -t $((cpu/maxProc)) \
        -m r941_min_sup_g507 \
        2>&1 | tee "${polished}"/"${ass}"/"${sample}"/medaka.log

    # Rename sample
    mv "${polished}"/"${ass}"/"${sample}"/consensus.fasta \
        "${polished}"/"${ass}"/"${sample}"/"${sample}"_medaka.fasta

    # Remove temporary files
    find "${polished}"/"${ass}"/"${sample}" -type f \
        ! -name "*_medaka.fasta" \
        ! -name "*.gfa" \
        ! -name "*.log" \
        -exec rm {} \;

    # Reformat fasta
    python "${scripts}"/format_fasta.py \
        "${polished}"/"${ass}"/"${sample}"/"${sample}"_medaka.fasta \
        "${polished}"/"${ass}"/"${sample}"/"${sample}"_medaka.fasta1
    mv "${polished}"/"${ass}"/"${sample}"/"${sample}"_medaka.fasta1 \
        "${polished}"/"${ass}"/"${sample}"/"${sample}"_medaka.fasta
}

export -f polish_medaka

conda activate medaka

find "${assemblies}"/"$ass" -maxdepth 2 -type f -name "*.fasta" \
    | parallel  --bar \
                --env polish_medaka \
                --env polished \
                --env scripts \
                --env cpu \
                --env ass \
                --env maxProc \
                'polish_medaka {}'

conda deactivate


# fixstart
function fix_start()
{
    sample=$(basename "$1" "_medaka.fasta")

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

find "${polished}"/"$ass" -type f -name "*.fasta" \
    | parallel  --bar \
                --env fix_start \
                --env maxProc \
                --env cpu \
                --env qc \
                --env ass \
                --jobs "$maxProc" \
                "fix_start {}"

# cleanup
find "${baseDir}"/fixstart/"$ass" -type f ! -name "*.fasta" ! -name "*.log" -exec rm {} \;


# Quast
conda activate quast

# Create list of all assemblies
declare -a genomes=()
for i in $(find "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta"); do 
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

find "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta" \
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
conda activate qualimap

function get_coverage()  # unsing unmerged reads only
{
    # sample=$(basename "$1" | cut -d '_' -f 1)
    sample=$(basename "$1" ".fasta")
    fastq_file="${filtered}"/"${sample}"_filtered.fastq.gz

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

find "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta" | \
    parallel    --bar \
                --env get_coverage \
                --env cpu \
                --env maxProc \
                --env ass \
                --env qc \
                --jobs "$maxProc" \
                "get_coverage {}"

# Qualimap
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

find "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta" \
    | parallel    --bar \
                --env id_sample \
                --env qc \
                "id_sample {}"

conda deactivate


# Annotate assemblies wit Prokka
conda activate prokka

function annotate()
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

export -f annotate

find "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta" | \
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
                --env ass \
                --env maxProc \
                --env ass \
                --env scripts \
                --jobs "$maxProc" \
                "annotate {}"

conda deactivate



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

find "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta" | \
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

find "${baseDir}"/fixstart/"$ass" -type f -name "*.fasta" \
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

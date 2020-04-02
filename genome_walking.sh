

function map()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$2"))

    [ -d "${clustered}"/"$sample" ] || mkdir -p "${clustered}"/"$sample"

    minimap2 \
        -x map-ont \
        -a \
        -t 4 \
        "$1" "$2" \
    | samtools view -@ 4 -bhF 4 - \
    | samtools sort -@ 4 > "${clustered}"/"${sample}"/"${sample}".bam

    # Index bam
    samtools index "${clustered}"/"${sample}"/"${sample}".bam
}

export -f map

export clustered="${baseDir}"/clustered
[ -d "$clustered" ] || mkdir -p "$clustered"
find -L "$filtered" -type f -name "*.fastq.gz" \
| parallel --bar \
            --env clustered \
            --env map \
            --jobs 4 \
            'map /home/marco/analyses/CAN-GB-GW-CBDAS/E55107_Cannabis_sativa_CBDAS.fasta {}'


function gwalk ()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    # extract soft clipped reads on the left side of the target based on cigar
    # not in reverse complement (which are in the right side of the target)

    # $6 -> Cigar string (S=soft clipped)
    # $4 -> 1-based leftmost mapping POSition (relative to the read)
    # $2 == 16 -> reverse complement (POS 1 is still the first base on the left)
    if [ "$2" = "left" ]; then
        bamtools convert -format sam -in "$1" \
        | awk -F $'\t' 'BEGIN {OFS = FS} {
        if($6 ~ /^[0-9][0-9][0-9]S/) {
            print
        }
        else if($6 ~ /^[0-9][0-9][0-9][0-9]S/) {
            print
        }
        else if($6 ~ /^[0-9][0-9][0-9][0-9][0-9]S/) {
            print
        }
        else if ($0 ~ /^@/) {
            print
        }
        }' | \
        samtools view -@ 4 -bhF 4 - | \
        samtools sort -@ 4 > "${1%.bam}_left.bam"
        samtools index "${1%.bam}_left.bam"
    elif [[ "$2" == "right" ]]; then
        bamtools convert -format sam -in "$1" \
        | awk -F $'\t' 'BEGIN {OFS = FS} {
        if($6 ~ /[0-9][0-9][0-9]S$/) {
            print
        }
        else if($6 ~ /[0-9][0-9][0-9][0-9]S$/) {
            print
        }
        else if($6 ~ /[0-9][0-9][0-9][0-9][0-9]S$/) {
            print
        }
        else if ($0 ~ /^@/) {
            print
        }
        }' | \
        samtools view -@ 4 -bhF 4 - | \
        samtools sort -@ 4 > "${1%.bam}_right.bam"
        samtools index "${1%.bam}_right.bam"
    else
        echo "Please provide a direction for the walking (\"left\" or \"right\")"
    fi

    #visualize filtered alignment
    # tablet "${output_folder}"/"${name}"_"${direction}"_overhang.bam "$ref" "${ref%.fasta}.gff" &
}

export -f gwalk

find -L "$clustered" -name "*.bam" \
| parallel --bar \
            --env gwalk \
            --jobs 4 \
            'gwalk {} left'

find -L "$clustered" -name "*.bam" ! -name "*_left.bam" \
| parallel --bar \
            --env gwalk \
            --jobs 4 \
            'gwalk {} right'


function bam2fastq()
{
     samtools fastq -@ 4 "$1" \
     | pigz > "${1%.bam}.fastq.gz"
}

export -f bam2fastq

find -L "$clustered" -name "*_left.bam" -o -name "*_right.bam" \
| parallel --bar --env bam2fastq --jobs 4 'bam2fastq {}'


# isonclust
function cluster()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))
    side=$(basename "$1" ".fastq.gz" | rev | cut -d "_" -f 1 | rev)

    [ -d "${clustered}"/"${sample}"/isonclust/"${side}" ] || mkdir -p "${clustered}"/"${sample}"/isonclust/"${side}"

    # Uncompress
    pigz -dkc "$1" > "${clustered}"/"${sample}"/isonclust/"${side}"/"${sample}"_"${side}".fastq

    isONclust \
        --fastq "${clustered}"/"${sample}"/isonclust/"${side}"/"${sample}"_"${side}".fastq \
        --t 4 \
        --ont \
        --consensus \
        --outfolder "${clustered}"/"${sample}"/isonclust/"${side}"

    # rename
    mv "${clustered}"/"${sample}"/isonclust/"${side}"/consensus_references.fasta \
        "${clustered}"/"${sample}"/isonclust/"${side}"/"${sample}"_"${side}".fasta

    # Cleanup
    rm "${clustered}"/"${sample}"/isonclust/"${side}"/"${sample}"_"${side}".fastq
    find "${clustered}"/"${sample}"/isonclust/"${side}" -mindepth 1 ! -name "*.fasta"\
        -exec rm -rf {} \;
}

export -f cluster

conda activate isonclust

find -L "$clustered" -type f -name "*_left.fastq.gz" -o -name "*_right.fastq.gz" \
| parallel --bar --env cluster --env clustered --jobs 4 \
    'cluster {}'
conda deactivate




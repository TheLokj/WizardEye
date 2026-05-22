#!/bin/bash
# ==================================================================================================
# Max Planck Institute for Evolutionary Anthropology 
# Lesage Louison | louison_lesage[at]eva.mpg.de
# Last update: 12.03.2026
# ==================================================================================================
# This script is designed to generate a cross-mappability filter for X given genomes.
# Each fasta are splitted into k-mers and cross-mappability is computed against the target.
# This allows the identification of area of the target where other k-mers align.
# ==================================================================================================
# This script is the BWA-version of the pipeline. 
# This script is useful to compute stringency filter and to know exactly
# the depth of unique overlapping k-mer on a position. 
# It however consumes a lot ot disk storage and is slower than the GenMap one.
# It requires seqkit, parallel, bwa, samtools, bedGraphToBigWig and bedtools.
# ==================================================================================================

usage() {
    echo "Usage: $0 -i=<input_fasta_directory> -t=<input_target> -k=<kmer_length> [options]"
    echo
    echo "Required Arguments:"
    echo "  -i, --input_fasta_directory    Path to directory containing FASTA files to split into k-mers."
    echo "  -t, --input_target             Path to the target reference FASTA (must be located in the input_fasta_directory)."
    echo "  -k, --kmer_length              Length of k-mers to produce."
    echo
    echo "Optional Arguments:"
    echo "  -bn, --bwa_missing_prob_err_rate    Missing probability error rate, bwa aln -n (default: 0.01)."
    echo "  -bo, --bwa_max_gap_opens            Maximum number or fraction of gap opens, bwa aln -o (default: 2)."
    echo "  -bl, --bwa_seed_length              Seed length, bwa aln -l (default: 16500)."
    echo "  -w,  --offset_step                  Offset/sliding window step for k-mers (default: 1)."
    echo "  -cs, --chunk_size                   Number of k-mers in a chunk. Chunks are then splitted between the --n_threads."
    echo "  -s,  --cross_stringency             Minimum ratio (0.0 to 1.0) of possible overlapping k-mers which require to be unique to mask a base during cross-mappability (default: 0.99)"
    echo "  -o,  --output_prefix                Output directory/prefix (default: ./output/)."
    echo "  -j,  --n_threads                    Number of threads for parallelisation (default: 1)."
    echo "  -ka, --keep_all                     Keep all k-mers Irmation and save positions in a BEDPE format (0 or 1, default: 0)."
    echo 
    echo "Description:"
    echo "  This script generates a mappability filter excluding the target and the region where k-mers from other FASTA align."
    echo "  It also produces also a bedgraph file per FASTA to allows futher analyses with a genome browser." 
    echo "  It requires samtools, bedGraphToBigWig, bedtools and bwa."
    exit 0
}
set -euo pipefail
output_prefix="./output/"
n_threads=1
offset_step=1
chunk_size=2000000
cross_stringency=0.99         

bwa_missing_prob_err_rate=0.01
bwa_max_gap_opens=2
bwa_seed_length=16500
keep_all=0

if [[ "$*" == *"-h"* || "$*" == *"--help"* || $# -eq 0 ]]; then usage; fi

for arg in "$@"; do
    case $arg in
        -i=*|--input_fasta_directory=*)
            input_fasta_directory="${arg#*=}"
            ;;
        -t=*|--input_target=*)
            input_target="${arg#*=}"
            ;;
        -k=*|--kmer_length=*)
            kmer_length="${arg#*=}"
            ;;
        -w=*|--offset_step=*)
            offset_step="${arg#*=}"
            ;;
        -bn=*|--bwa_missing_prob_err_rate=*)
            bwa_missing_prob_err_rate="${arg#*=}"
            ;;
        -bo=*|--bwa_max_gap_opens=*)
            bwa_max_gap_opens="${arg#*=}"
            ;;
        -bl=*|--bwa_seed_length=*)
            bwa_seed_length="${arg#*=}"
            ;;
        -s=*|--cross_stringency=*|-rc=*)
            cross_stringency="${arg#*=}"
            ;;
        -o=*|--output_prefix=*)
            output_prefix="${arg#*=}"
            ;;
        -j=*|--n_threads=*)
            n_threads="${arg#*=}"
            ;;
        -ka=*|--keep_all=*)
            keep_all="${arg#*=}"
            ;;
        -cs=*|--chunk_size=*)
            chunk_size="${arg#*=}"
            ;;
        *)
            echo "Unknown option $arg"
            exit 1
            ;;
    esac
done

if [[ -z "$input_fasta_directory" || -z "$input_target" || -z "$kmer_length" ]]; then
    echo "Error: Missing required arguments."; usage
fi

for tool in bwa samtools bedtools seqkit parallel bedGraphToBigWig; do
    if ! command -v "$tool" &>/dev/null; then
        echo "Error: '$tool' not found in PATH."; exit 1
    fi
done

if ! [[ "$kmer_length" =~ ^[0-9]+$ ]] || [[ "$kmer_length" -lt 1 ]]; then
    echo "Error: kmer_length must be a positive integer."; exit 1
fi
if ! awk -v v="$cross_stringency" 'BEGIN { exit !(v >= 0.0 && v <= 1.0) }'; then
    echo "Error: cross_stringency must be between 0.0 and 1.0."; exit 1
fi
if ! [[ "$n_threads" =~ ^[0-9]+$ ]] || [[ "$n_threads" -lt 1 ]]; then
    echo "Error: n_threads must be a positive integer."; exit 1
fi

LOCAL_SCRATCH="${TMPDIR:-/tmp}/${USER}_BWA_cross_map_$$"
trap 'kill $(jobs -p) 2>/dev/null || true; rm -rf "${LOCAL_SCRATCH}"' EXIT TERM INT
echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Pipeline temporary location: $LOCAL_SCRATCH/"
mkdir -p "$LOCAL_SCRATCH"
mkdir -p "$output_prefix"

target_basename=$(basename "$input_target")
target_filename="${target_basename%.*}"
overlap_dir="${output_prefix}tracks"
mkdir -p "$overlap_dir"

# AWK scripts for parsing XA tags in BWA output
# Extracts all mapping positions (primary + XA repetitions).
cat << 'EOF' > "$LOCAL_SCRATCH/parse_xa_all.awk"
    BEGIN { OFS="\t" }
    function get_ref_len(cigar, def_k,   len, tmp_num, char, i) {
        len = 0; tmp_num = "";
        for (i = 1; i <= length(cigar); i++) {
            char = substr(cigar, i, 1);
            if (char ~ /[0-9]/) { tmp_num = tmp_num char; }
            else { if (char ~ /[MDN=X]/) len += (tmp_num + 0); tmp_num = ""; }
        }
        return len > 0 ? len : def_k;
    }
    /^@/ { next }
    $4 > 0 {
        c_len = get_ref_len($6, k);
        print $3, $4 - 1, $4 - 1 + c_len;
        for (i=12; i<=NF; i++) {
            if ($i ~ /^XA:Z:/) {
                split(substr($i, 6), hits, ";");
                for (j=1; j<=length(hits); j++) {
                    if (hits[j] != "") {
                        split(hits[j], fields, ",");
                        pos = int(fields[2]); if (pos < 0) pos = -pos;
                        print fields[1], pos - 1, pos - 1 + get_ref_len(fields[3], k);
                    }
                }
            }
        }
    }
EOF

# Extracts only primary mapping positions for reads with MAPQ > 0 and no XA tag (i.e., uniquely mapping reads).
cat << 'EOF' > "$LOCAL_SCRATCH/parse_xa_uniq.awk"
    BEGIN { OFS="\t" }
    function get_ref_len(cigar, def_k,   len, tmp_num, char, i) {
        len = 0; tmp_num = "";
        for (i = 1; i <= length(cigar); i++) {
            char = substr(cigar, i, 1);
            if (char ~ /[0-9]/) { tmp_num = tmp_num char; }
            else { if (char ~ /[MDN=X]/) len += (tmp_num + 0); tmp_num = ""; }
        }
        return len > 0 ? len : def_k;
    }
    /^@/ { next }
    $4 > 0 && $5 > 0 {
        is_unique = 1;
        for (i=12; i<=NF; i++) if ($i ~ /^XA:Z:/ && length($i) > 5) { is_unique = 0; break; }
        if (is_unique) print $3, $4 - 1, $4 - 1 + get_ref_len($6, k);
    }
EOF

if [[ ! -f "${input_target}.bwt" ]]; then
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] BWA target index not found. Indexing..."
    bwa index "$input_target"
else
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] BWA target index found. Skipping this step."
fi

if [[ ! -f "${output_prefix}target.genome.sizes" ]]; then
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Genome sizes not found. Computing..."
    samtools faidx "$input_target"
    cut -f1,2 "${input_target}.fai" | LC_ALL=C sort -k1,1 > "${output_prefix}target.genome.sizes"
fi

echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Processing input one by one for cross-mappability..."

shopt -s nullglob
for other_path in "$input_fasta_directory"/*.{fa,fasta}; do
    other_base=$(basename "$other_path")
    other_file="${other_base%.*}"

    if [[ "$other_file" == "$target_filename" ]]; then continue; fi

    if [[ ! -f "$overlap_dir/cov_uniq_${other_file}.bg" ]]; then

        echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Comparing $target_filename vs $other_file"
        tmp_local="$LOCAL_SCRATCH/splits_${other_file}"
        mkdir -p "$tmp_local"

        echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Generating ${kmer_length}-mers from $other_file..."
        if [ "$keep_all" = 1 ]; then
            seqkit sliding -W "$kmer_length" -s "$offset_step" -i "$other_path" | \
            seqkit split2 -s "$chunk_size" -O "$tmp_local" --extension .fasta
        else
            seqkit sliding -W "$kmer_length" -s "$offset_step" "$other_path" | \
            seqkit rmdup -s | \
            seqkit split2 -s "$chunk_size" -O "$tmp_local" --extension .fasta
        fi

        echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Aligning k-mers from $other_file on $target_filename..."

        # Use GNU parallel to run BWA alignments in parallel, then convert to BAM and merge into a single BAM file for counting.
        find "$tmp_local/" -name "*.fasta" | parallel -j "$n_threads" \
            "bwa aln -t 1 -R 2147483647 -n $bwa_missing_prob_err_rate -o $bwa_max_gap_opens \
                -l $bwa_seed_length '$input_target' {} > {.}.sai && \
            bwa samse -n 2147483647 '$input_target' {.}.sai {} | \
            samtools view -b -F 4 - > {.}.bam && \
            rm -f {} {.}.sai"

        # To avoid I/O related error (i.e. an access to a file non perfectly indexed yet by the system), wait a few seconds 
        sleep 10

        samtools cat "$tmp_local/"*.bam | samtools sort -@ "$n_threads" -o "$tmp_local/merged.bam"

        # If required, get the exact mapping positions in a BEDPE format
        # This is useful to analyze precisely where the k-mers from the other FASTA are from and where they map on the target.
        if [ "$keep_all" -eq 1 ]; then
            echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Generating a BEDPE to save exact mapping positions..."

            COMPRESS_TOOL=$(command -v bgzip || echo "gzip")
            
            samtools view "$tmp_local/merged.bam" | \
            awk -v k="$kmer_length" 'BEGIN {OFS="\t"} {
                split($1, a, ":");
                src_chr = a[1];
                split(a[2], b, "-");
                src_start = b[1] - 1; 
                src_end = b[2];

                tgt_chr = $3;
                tgt_start = $4 - 1;
                tgt_end = tgt_start + k; 

                strand_tgt = (and($2, 16) ? "-" : "+");

                print src_chr, src_start, src_end, tgt_chr, tgt_start, tgt_end, $1, $5, "+", strand_tgt;
            }' | $COMPRESS_TOOL > "$overlap_dir/mapping_${other_file}_to_${target_filename}.bedpe.gz"

            if command -v tabix &>/dev/null && [ "$COMPRESS_TOOL" == "bgzip" ]; then
                tabix -p bed "$overlap_dir/mapping_${other_file}_to_${target_filename}.bedpe.gz"
            fi
        fi

        echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Computing total and unique depths..."

        # For each position, count the number of k-mers that map (n_map) and the number of k-mers that uniquely map (n_uniq)
        # Mismatches and gaps are already accounted for in the CIGAR string, so we can directly use the alignment positions and lengths.

        # N_MAP: all hits (primary + XA repetitions)
        samtools view -F 4 "$tmp_local/merged.bam" | \
        awk -v k="$kmer_length" -f "$LOCAL_SCRATCH/parse_xa_all.awk" | \
        LC_ALL=C sort -T "$tmp_local" -k1,1 -k2,2n > "$tmp_local/n_map.bed"

        # N_UNIQ : hits without XA, MAPQ > 0
        samtools view -F 4 "$tmp_local/merged.bam" | \
        awk -v k="$kmer_length" -f "$LOCAL_SCRATCH/parse_xa_uniq.awk" | \
        LC_ALL=C sort -T "$tmp_local" -k1,1 -k2,2n > "$tmp_local/n_uniq.bed"

        # Compute coverage for both n_map and n_uniq, save as bedGraph to be able to change stringency threshold without reprocessing the BAM and to be able to deal with depth
        bedtools genomecov -i "$tmp_local/n_map.bed" -g "${output_prefix}target.genome.sizes" -bga | \
            LC_ALL=C sort -k1,1 -k2,2n > "$overlap_dir/cov_map_${other_file}.bg"
        bedtools genomecov -i "$tmp_local/n_uniq.bed" \
            -g "${output_prefix}target.genome.sizes" -bga | \
            LC_ALL=C sort -k1,1 -k2,2n > "$overlap_dir/cov_uniq_${other_file}.bg"
            
        bedGraphToBigWig "$overlap_dir/cov_map_${other_file}.bg" "${output_prefix}target.genome.sizes" "$overlap_dir/cov_map_${other_file}.bw"
        bedGraphToBigWig "$overlap_dir/cov_uniq_${other_file}.bg" "${output_prefix}target.genome.sizes" "$overlap_dir/cov_uniq_${other_file}.bw"

        rm -rf "$tmp_local"

        sleep 10
    fi

    # Compute the ratio unique_depth / max_theoretical_depth and apply the stringency filter
    if [[ ! -f "$overlap_dir/overlap_${other_file}_s${cross_stringency}_k${kmer_length}.bed" ]]; then
        echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Applying cross-mappability stringency threshold of $cross_stringency to $other_file..."
        
        awk -v r="$cross_stringency" -v k="$kmer_length" -v s="$offset_step" 'OFS="\t" {
            n_uniq = $4 + 0;
            max_kmers = k / s;
            if (max_kmers > 0 && (n_uniq / max_kmers) >= r) print $1, $2, $3;
        }' "$overlap_dir/cov_uniq_${other_file}.bg" | \
        bedtools merge -i stdin > "$overlap_dir/overlap_${other_file}_s${cross_stringency}_k${kmer_length}.bed"
    fi

done

overlap_files=("$overlap_dir"/overlap_*_s"${cross_stringency}_k${kmer_length}.bed")
shopt -u nullglob

echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Merging all overlaps into final mask..."

if [ ${#overlap_files[@]} -gt 0 ]; then
    cat "${overlap_files[@]}" | sort -k1,1 -k2,2n | \
        bedtools merge > "${output_prefix}all_overlaps_mask_${cross_stringency}_k${kmer_length}.bed"
else
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] No overlap found."
fi

rm -rf "$LOCAL_SCRATCH"
echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Cross-mappability done!"
#!/bin/bash
#SBATCH --job-name=finemap
#SBATCH --partition=GPU2
#SBATCH --cpus-per-task=8
#SBATCH --output=/gpfs/hpc/home/lijc/zhangsy/Sci.Bull.Revision/Re-analysis/Interaction/Sentinel.Variants/FINEMAP_RESULTS/log/finemap_%A_%a.out
#SBATCH --error=/gpfs/hpc/home/lijc/zhangsy/Sci.Bull.Revision/Re-analysis/Interaction/Sentinel.Variants/FINEMAP_RESULTS/log/finemap_%A_%a.err
#SBATCH --array=0-1139   # adjust based on number of lines in manifest

set -euo pipefail

# ==============================================================================
#  Environment
# ==============================================================================
source ~/anaconda3/etc/profile.d/conda.sh
conda activate FINEMAP

# ==============================================================================
#  Inputs
# ==============================================================================
MANIFEST="${1:-master_manifest.txt}"
N_SAMPLES="${2:-303810}"   # total number of samples; can override at submit

# Get master.csv line corresponding to current array task
MASTER="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$MANIFEST")"
if [[ -z "${MASTER:-}" ]]; then
  echo "No master file for index ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi
if [[ ! -f "$MASTER" ]]; then
  echo "Missing master.csv: $MASTER" >&2
  exit 3
fi

# ==============================================================================
#  LD files
# ==============================================================================
DIR="$(dirname "$MASTER")"   # .../LD_matrix/chrXX/CVYY/region
LD_GZ="$DIR/region.ld.gz"
LD_BGZ="$DIR/region.ld.bgz"
LD_TXT="$DIR/region.ld"

if   [[ -s "$LD_GZ"  ]]; then LD_FILE="$LD_GZ"
elif [[ -s "$LD_BGZ" ]]; then LD_FILE="$LD_BGZ"
elif [[ -s "$LD_TXT" ]]; then LD_FILE="$LD_TXT"
else
  echo "No LD file found in $DIR (expected region.ld[.gz|.bgz|.txt])" >&2
  exit 4
fi

# ==============================================================================
#  Summary statistics
# ==============================================================================
SUMSTATS_DIR="${DIR/LD_matrix/SUMSTATS}"
if   [[ -s "$SUMSTATS_DIR/sumstats.txt.gz" ]]; then
  SUMSTATS="$SUMSTATS_DIR/sumstats.txt.gz"
elif [[ -s "$SUMSTATS_DIR/sumstats.txt" ]]; then
  SUMSTATS="$SUMSTATS_DIR/sumstats.txt"
else
  echo "Missing sumstats in $SUMSTATS_DIR (looked for sumstats.txt[.gz])" >&2
  exit 5
fi

# ==============================================================================
#  Output directory
# ==============================================================================
CHR="$(basename "$(dirname "$(dirname "$DIR")")")"   # chrXX
CV="$(basename "$(dirname "$DIR")")"                 # CVYY
REGION="$(basename "$DIR")"

OUT_DIR="/gpfs/hpc/home/lijc/zhangsy/Sci.Bull.Revision/Re-analysis/Interaction/Sentinel.Variants/FINEMAP_RESULTS/${CHR}/${CV}/${REGION}"
mkdir -p "$OUT_DIR"

echo "[START] ${CHR}/${CV}/${REGION}  $(date '+%F %T')"

# ==============================================================================
#  Run FINEMAP
# ==============================================================================
python3.9 /gpfs/hpc/home/lijc/zhangsy/software/fine-mapping-inf/run_fine_mapping.py \
  --sumstats "$SUMSTATS" \
  --beta-col-name beta \
  --se-col-name se \
  --ld-file "$LD_FILE" \
  --n "$N_SAMPLES" \
  --save-npz \
  --save-tsv \
  --output-prefix "$OUT_DIR/Finemapping_${REGION}"

echo "[DONE ] ${CHR}/${CV}/${REGION}  $(date '+%F %T')"

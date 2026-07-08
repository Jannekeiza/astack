#!/bin/bash
# Usage: ./move_if_greater.sh <source_dir> <destination_dir>
# Loops over all .aq files in source_dir.
# If the first line of a .aq file is a number > 2, moves both the .aq file
# and its matching SNR_*.txt file to destination_dir.
# SNR file pattern: SNR_<basename>.txt  e.g. SNR_230618_203021_0.02-1.0Hz.txt
#                   for .aq file:             230618_203021_0.02-1.0Hz.aq

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <source_dir> <destination_dir>"
  exit 1
fi

SRC="$1"
DEST="$2"

if [[ ! -d "$SRC" ]]; then
  echo "Error: Source directory '$SRC' does not exist."
  exit 1
fi

shopt -s nullglob
TTR_FILES=("$SRC"/*.ttr)

if [[ ${#TTR_FILES[@]} -eq 0 ]]; then
  echo "No .ttr files found in '$SRC'."
  exit 0
fi

moved=0
skipped=0
errors=0

for TTR_FILE in "${TTR_FILES[@]}"; do
  BASENAME=$(basename "$TTR_FILE" .ttr)   # e.g. 230618_203021_0.02-1.0Hz

  # Read the first line and strip whitespace
  FIRST_LINE=$(head -n 1 "$TTR_FILE" | tr -d '[:space:]')

  # Validate it's a number
  if ! [[ "$FIRST_LINE" =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
    echo "  [SKIP]  '$TTR_FILE' — first line is not a valid number ('$FIRST_LINE')."
    ((errors++)) || true
    continue
  fi

  # Check if value > 2
  if awk "BEGIN { exit !($FIRST_LINE > 2) }"; then
    mkdir -p "$DEST"

    cp "$TTR_FILE" "$DEST/"
    echo "  [MOVED] '$(basename "$AQ_FILE")' (value: $FIRST_LINE)"

    ((moved++)) || true
  else
    echo "  [SKIP]  '$(basename "$AQ_FILE")' (value: $FIRST_LINE ≤ 2)"
    ((skipped++)) || true
  fi
done

echo ""
echo "Done. Copied: $moved file(s), Skipped: $skipped file(s), Errors: $errors file(s)."

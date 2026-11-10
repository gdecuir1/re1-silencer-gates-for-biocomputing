#!/usr/bin/env bash
# Download target genomes as .2bit from UCSC.
#
# .2bit rather than .fa.gz on purpose: about a quarter the size, and
# `twoBitToFa -bed regions.bed genome.2bit out.fa` extracts the lifted
# intervals directly, so the full FASTA is never needed on disk.
#
# Verifies every URL resolves BEFORE downloading anything, so a wrong
# assembly name costs five seconds instead of forty minutes.
#
# Usage:  bash scripts/fetch_genomes.sh [--check-only] [--yes] [assembly ...]
#         bash scripts/fetch_genomes.sh mm39
#         bash scripts/fetch_genomes.sh --yes mm39
#
# With no assembly names, every target species is included. --yes skips the
# confirmation prompt for an explicit caller such as src.fetch_mouse_scoring_inputs.

set -euo pipefail

BASE="https://hgdownload.soe.ucsc.edu/goldenPath"
DEST="data/raw/genomes"
PROVENANCE="data/README.md"

# assembly:common name
ASSEMBLIES=(
  "rheMac10:rhesus macaque"
  "mm39:mouse"
  "canFam6:dog"
  "bosTau9:cow"
  "monDom5:opossum"
)

CHECK_ONLY=false
ASSUME_YES=false
REQUESTED=()
for arg in "$@"; do
  case "$arg" in
    --check-only) CHECK_ONLY=true ;;
    --yes) ASSUME_YES=true ;;
    -*)
      echo "Unknown option: $arg" >&2
      exit 2
      ;;
    *) REQUESTED+=("$arg") ;;
  esac
done

if [[ ${#REQUESTED[@]} -gt 0 ]]; then
  FILTERED=()
  seen=" "
  for name in "${REQUESTED[@]}"; do
    found=""
    for entry in "${ASSEMBLIES[@]}"; do
      if [[ "${entry%%:*}" == "$name" ]]; then
        found="$entry"
        break
      fi
    done
    if [[ -z "$found" ]]; then
      echo "Unknown assembly: $name" >&2
      echo "Choose from: rheMac10 mm39 canFam6 bosTau9 monDom5" >&2
      exit 2
    fi
    if [[ "$seen" != *" $name "* ]]; then
      FILTERED+=("$found")
      seen+="$name "
    fi
  done
  ASSEMBLIES=("${FILTERED[@]}")
fi

mkdir -p "$DEST"

# ---------- phase 1: verify every URL resolves ----------

echo "Checking URLs..."
echo

printf "%-12s %-18s %12s  %s\n" "ASSEMBLY" "SPECIES" "SIZE" "STATUS"
printf "%-12s %-18s %12s  %s\n" "--------" "-------" "----" "------"

total_bytes=0
missing=0

for entry in "${ASSEMBLIES[@]}"; do
  db="${entry%%:*}"
  species="${entry#*:}"
  url="$BASE/$db/bigZips/$db.2bit"

  # -s silent, -I head, -L follow redirects, -f fail on 4xx/5xx
  headers=$(curl -sIL -f --max-time 30 "$url" 2>/dev/null) || {
    printf "%-12s %-18s %12s  %s\n" "$db" "$species" "-" "NOT FOUND"
    missing=$((missing + 1))
    continue
  }

  # last Content-Length wins after redirects
  bytes=$(printf '%s' "$headers" \
    | awk 'BEGIN{IGNORECASE=1} /^content-length:/ {gsub(/\r/,""); n=$2} END{print n+0}')

  if [[ "$bytes" -gt 0 ]]; then
    human=$(awk -v b="$bytes" 'BEGIN{printf "%.2f GB", b/1073741824}')
    total_bytes=$((total_bytes + bytes))
  else
    human="unknown"
  fi

  printf "%-12s %-18s %12s  %s\n" "$db" "$species" "$human" "ok"
done

echo
if [[ "$missing" -gt 0 ]]; then
  echo "$missing assembly URL(s) did not resolve."
  echo "Check the current name at https://hgdownload.soe.ucsc.edu/downloads.html"
  echo "and fix the ASSEMBLIES list at the top of this script."
  exit 1
fi

awk -v b="$total_bytes" 'BEGIN{printf "Total to download: %.2f GB\n", b/1073741824}'

# warn on free space (macOS and Linux differ; df -k is portable enough)
avail_kb=$(df -k . | awk 'NR==2 {print $4}')
awk -v a="$avail_kb" 'BEGIN{printf "Free disk here:    %.2f GB\n", a/1048576}'
echo

$CHECK_ONLY && { echo "Check only, stopping."; exit 0; }

if $ASSUME_YES; then
  echo "Downloading without a prompt (--yes)."
else
  read -r -p "Download these to $DEST? [y/N] " reply
  [[ "$reply" =~ ^[Yy]$ ]] || { echo "Cancelled."; exit 0; }
fi

# ---------- phase 2: download, resumable ----------

echo
for entry in "${ASSEMBLIES[@]}"; do
  db="${entry%%:*}"
  url="$BASE/$db/bigZips/$db.2bit"
  out="$DEST/$db.2bit"

  echo "=== $db ==="
  # -C - resumes a partial file; safe to re-run after an interrupted download
  curl -L -C - --retry 3 --retry-delay 5 -o "$out" "$url"

  # UCSC publishes md5sum.txt per bigZips dir
  if curl -sL -f --max-time 30 "$BASE/$db/bigZips/md5sum.txt" -o "$DEST/$db.md5.txt" 2>/dev/null; then
    expected=$(awk -v f="$db.2bit" '$2 == f {print $1}' "$DEST/$db.md5.txt" | head -1)
    if [[ -n "$expected" ]]; then
      if command -v md5sum >/dev/null 2>&1; then
        actual=$(md5sum "$out" | awk '{print $1}')
      else
        actual=$(md5 -q "$out")   # macOS
      fi
      if [[ "$actual" == "$expected" ]]; then
        echo "md5 ok"
      else
        echo "MD5 MISMATCH for $db.2bit"
        echo "  expected $expected"
        echo "  actual   $actual"
        echo "  Delete $out and re-run."
        exit 1
      fi
    else
      echo "md5 not listed for $db.2bit, skipping check"
    fi
    rm -f "$DEST/$db.md5.txt"
  else
    echo "md5sum.txt unavailable, skipping check"
  fi
  echo
done

# ---------- phase 3: provenance ----------

stamp=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
{
  echo ""
  echo "<!-- appended by scripts/fetch_genomes.sh -->"
  for entry in "${ASSEMBLIES[@]}"; do
    db="${entry%%:*}"
    species="${entry#*:}"
    out="$DEST/$db.2bit"
    [[ -f "$out" ]] || continue
    if command -v md5sum >/dev/null 2>&1; then
      sum=$(md5sum "$out" | awk '{print $1}')
    else
      sum=$(md5 -q "$out")
    fi
    size=$(awk -v b="$(wc -c < "$out")" 'BEGIN{printf "%.2f GB", b/1073741824}')
    echo "| $db.2bit | $species | $BASE/$db/bigZips/$db.2bit | $stamp | $size | $sum |"
  done
} >> "$PROVENANCE"

echo "Done. Files in $DEST, provenance appended to $PROVENANCE."
echo
echo "Extract lifted regions later with:"
echo "  twoBitToFa -bed=regions.bed $DEST/mm39.2bit out.fa"

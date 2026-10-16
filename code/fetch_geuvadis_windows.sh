#!/bin/bash
# =============================================================================
# fetch_geuvadis_windows.sh
# Disk-frugal E-GEUV-1 genotype fetch. The full per-chromosome VCFs are
# 3-4 GB EACH (annotv2), far too large to keep 22 of them, so per chromosome
# (1-7, the PolyFun-covered ones used by the gene list):
#   1. complete the download if needed (curl resume; skips when local size
#      matches the remote Content-Length),
#   2. extract the merged cis windows of the selected genes by STREAMING
#      (gzcat | awk on POS): the E-GEUV-1 VCFs are plain gzip, NOT BGZF, so
#      they cannot be tabix-indexed directly (that is also why the FTP has
#      no .tbi files),
#   3. bgzip + tabix the small subset VCF
#      (vcf_subsets/GEUVADIS.chr<c>.subset.vcf.gz),
#   4. DELETE the full VCF.
# Peak transient disk ~4 GB; final footprint = the subsets (tens of MB).
# Idempotent: a chromosome whose subset already exists is skipped.
#
# Usage (from code/):   bash fetch_geuvadis_windows.sh
# =============================================================================
set -u
cd "$(dirname "$0")/../geuvadis_data"
export PATH="/opt/homebrew/bin:$PATH"
BASE="http://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1/genotypes"
REG="geuvadis_regions.txt"
OUT="vcf_subsets"; mkdir -p "$OUT"
# Full VCFs live OUTSIDE the iCloud-synced Desktop tree: under disk pressure,
# iCloud evicts large files to dataless placeholders whose reads time out
# (bit us on chr3). /private/tmp is on the same volume and never synced.
FULL="/private/tmp/geuvadis_full"; mkdir -p "$FULL"
[ -s "$REG" ] || { echo "FETCH_FAIL: missing $REG (run make_geuvadis_gene_list.R first)"; exit 1; }

for c in 1 2 3 4 5 6 7; do
  sub="$OUT/GEUVADIS.chr${c}.subset.vcf.gz"
  if [ -s "$sub" ] && [ -s "$sub.tbi" ]; then echo "SKIP chr$c (subset exists)"; continue; fi
  regc="$OUT/regions_chr${c}.txt"
  awk -v c="$c" '$1==c' "$REG" > "$regc"
  if [ ! -s "$regc" ]; then echo "SKIP chr$c (no gene window)"; rm -f "$regc"; continue; fi

  fname="GEUVADIS.chr${c}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz"
  f="$FULL/$fname"

  # legacy copy in the synced folder: wait for iCloud re-materialization if
  # dataless, then move it out; give up after 30 min and re-download instead
  if [ -f "$fname" ] && [ ! -f "$f" ]; then
    apparent=$(stat -f%z "$fname")
    for i in $(seq 1 60); do
      ondisk=$(( $(stat -f%b "$fname") * 512 ))
      [ "$ondisk" -ge $(( apparent * 95 / 100 )) ] && break
      brctl download "$fname" 2>/dev/null || true
      sleep 30
    done
    ondisk=$(( $(stat -f%b "$fname") * 512 ))
    if [ "$ondisk" -ge $(( apparent * 95 / 100 )) ]; then
      mv "$fname" "$f"; echo ">>> chr$c: legacy copy recovered ($(du -h "$f" | cut -f1))"
    else
      echo ">>> chr$c: legacy copy still dataless after 30 min; discarding it"
      rm -f "$fname"
    fi
  fi

  remote=$(curl -sIL "$BASE/$fname" | awk 'tolower($1)=="content-length:"{gsub("\r","",$2); print $2}' | tail -1)
  local_size=0; [ -f "$f" ] && local_size=$(stat -f%z "$f")
  if [ -z "$remote" ] || [ "$local_size" != "$remote" ]; then
    echo ">>> chr$c: downloading ($local_size/${remote:-?} bytes present) $(date '+%T')"
    curl -sS --retry 5 -C - -o "$f" "$BASE/$fname" || { echo "FETCH_FAIL: download chr$c"; exit 1; }
  else
    echo ">>> chr$c: full VCF already complete"
  fi

  # stream-extract the windows (plain-gzip source; single-chromosome file, so
  # filtering on POS alone is sufficient)
  echo ">>> chr$c: streaming window extraction $(date '+%T')"
  gzcat "$f" | awk -v OFS='\t' '
      NR==FNR { lo[++n]=$2; hi[n]=$3; next }
      /^#/    { print; next }
      { for (i=1; i<=n; i++) if ($2>=lo[i] && $2<=hi[i]) { print; break } }
    ' "$regc" - | bgzip > "$sub" || { echo "FETCH_FAIL: extract chr$c"; exit 1; }
  if ! gzcat "$sub" 2>/dev/null | grep -v "^#" | head -1 | grep -q .; then
    echo "FETCH_FAIL: empty subset chr$c"; exit 1
  fi
  tabix -f -p vcf "$sub" || { echo "FETCH_FAIL: index subset chr$c"; exit 1; }
  echo "OK chr$c subset: $(du -h "$sub" | cut -f1), $(bcftools view -H "$sub" | wc -l | tr -d ' ') records"
  rm -f "$f" "$f.tbi" "$regc"
done

rm -f GEUVADIS.chr*.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz \
      GEUVADIS.chr*.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz.tbi 2>/dev/null || true
rm -f "$FULL"/GEUVADIS.chr*.vcf.gz "$FULL"/GEUVADIS.chr*.vcf.gz.tbi 2>/dev/null || true
rmdir "$FULL" 2>/dev/null || true
echo "WINDOWS_FETCH_DONE ($(du -sh "$OUT" | cut -f1) of subsets)"

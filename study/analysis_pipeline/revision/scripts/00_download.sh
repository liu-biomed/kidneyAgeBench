#!/bin/bash
# Download public inputs into ~/kidney_aging_data (outside iCloud).
# NCBI via ftp:// with resume (-C -) and retries, because the local https route is throttled
# and ftp transfers are dropped intermittently. Completed .gz files are verified with gzip -t.
D=~/kidney_aging_data; N=ftp://ftp.ncbi.nlm.nih.gov/geo; L=$D/urls.txt; : > $L
stub(){ x=$1; n=${x//[!0-9]/}; p=${x%%[0-9]*}; if [ ${#n} -le 3 ]; then echo "${p}nnn"; else echo "${p}${n:0:$((${#n}-3))}nnn"; fi; }
sm(){ echo "$N/series/$(stub $1)/$1/matrix/${1}_series_matrix.txt.gz $D/geo/${1}_series_matrix.txt.gz" >> $L; }
for g in GSE36895 GSE60807 GSE21785 GSE147451 GSE66272 GSE232249 GSE167093 GSE93069 GSE217427 GSE175759 GSE66494; do sm $g; done
for p in GPL22945 GPL24120; do echo "$N/series/GSE104nnn/GSE104954/matrix/GSE104954-${p}_series_matrix.txt.gz $D/geo/GSE104954-${p}_series_matrix.txt.gz" >> $L; done
for p in GPL570 GPL6480 GPL96 GPL571; do echo "$N/platforms/$(stub $p)/$p/annot/$p.annot.gz $D/gpl/$p.annot.gz" >> $L; done
echo "$N/series/GSE93nnn/GSE93069/suppl/GSE93069_count_matrix.txt.gz $D/suppl/GSE93069_count_matrix.txt.gz" >> $L
echo "$N/series/GSE217nnn/GSE217427/suppl/GSE217427_expression_matrix.TPM.gct.gz $D/suppl/GSE217427_expression_matrix.TPM.gct.gz" >> $L
echo "$N/series/GSE167nnn/GSE167093/suppl/GSE167093_Pairs_normalization_matrix.xlsx $D/suppl/GSE167093_Pairs_normalization_matrix.xlsx" >> $L
fetch(){ url=$1; out=$2
  for i in $(seq 1 60); do
    curl --noproxy "*" -sS -C - --max-time 900 -o "$out" "$url" 2>/dev/null; rc=$?
    if [ $rc -eq 0 ] || [ $rc -eq 33 ]; then
      case "$out" in *.gz) gzip -t "$out" 2>/dev/null && return 0 ;; *) return 0 ;; esac
    fi; sleep 3; done; echo "FAIL $url"; }
export -f fetch
cat $L | xargs -P 8 -n 2 bash -c 'fetch "$0" "$1"'
# Platform tables without sample data (family SOFT files embed every sample and are tens of GB)
for p in GPL10558 GPL25759; do f=$D/gpl/${p}_table.txt; [ -s $f ] || curl --noproxy "*" -sS --retry 5 --max-time 1800 -o $f "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$p&targ=self&form=text&view=data"; done
echo DONE; find $D -type f | xargs ls -la | awk '{print $5,$9}'

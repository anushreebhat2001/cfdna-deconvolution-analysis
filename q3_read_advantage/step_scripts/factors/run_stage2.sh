#!/bin/bash
NP="${1:-12}"
O=results/factors
R=data/ref/31CT_U250_Netanel_merged_regions_file.bed
RM=data/ref/blood_epithelial_merged_regions_file_max_markers.txt
ONE=step_scripts/common/run_one.sh
ATL=step_scripts/factors/run_one_atlas.sh
{
  for f in data/mix_breast/*__d*.pat; do
    o="$O/t9_breast_dil/$(basename "$f" .pat)"
    echo "$ONE $f $o/whole u250"; echo "$ONE $f $o/cut u250 --size 1"
  done
  awk -F'\t' '$2=="control"{print $1}' meta/samples.tsv | while read s; do
    o="$O/t8_no_megak/$s"; f="data/marker_pat/$s.pat"
    echo "$ATL $f $o/whole $O/atlas/max_minus_megak.txt $RM"; echo "$ATL $f $o/cut $O/atlas/max_minus_megak.txt $RM --size 1"
  done
  for arm in off on; do
    { awk -F'\t' 'NR>1{print $1}' meta/samples.tsv; echo HU017.04.PL4596; echo HU017and19.05.PL4434
      echo HU017and19.06.PL4518; echo HU025.13.PL4704; } | while read s; do
      o="$O/t7_$arm/$s"; f="data/marker_pat/$s.pat"
      echo "$ATL $f $o/whole $O/atlas/u250_plus_$arm.txt $R"; echo "$ATL $f $o/cut $O/atlas/u250_plus_$arm.txt $R --size 1"
    done
  done
  for c in clumped_1 clumped_3 scattered_1 scattered_3 symmetric_5 symmetric_10; do
    for f in $O/sim/$c/*.pat; do
      o="$O/t10_$c/$(basename "$f" .pat)"
      echo "$ONE $f $o/whole u250"; echo "$ONE $f $o/cut u250 --size 1"
    done
  done
} | if [ -n "${DRY:-}" ]; then cat; else xargs -P "$NP" -L 1 bash; fi
[ -n "${DRY:-}" ] && exit 0
for t in t9_breast_dil t8_no_megak t7_off t7_on t10_clumped_1 t10_clumped_3 t10_scattered_1 t10_scattered_3 t10_symmetric_5 t10_symmetric_10; do
  echo "$t: $(find $O/$t -name '*.ish_clipped.txt' | wc -l) fits, $(find $O/$t -name '*.failed' | wc -l) failed, $(find $O/$t -name '*.empty' | wc -l) empty"
done

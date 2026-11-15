#!/bin/bash -euo pipefail
# Run GSEA

gsea-cli GSEA \
    -res host_therapy.gct \
    -cls host_therapy.cls#Recurrent_U2_versus_Primary_U2 \
    -gmx combined_rat.gmt \
    -chip Rattus_norvegicus.anno.feature_metadata.chip -collapse true \
    -out . \
    --rpt_label host_therapy.combined_rat \
    -nperm 1000 -permute gene_set -scoring_scheme weighted -metric Diff_of_Classes -sort real -order descending -set_max 500 -set_min 15 -norm meandiv -rnd_type no_balance -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed 1234 -save_rnd_lists false -zip_report false

# Un-timestamp the outputs for path consistency
mv host_therapy.combined_rat.Gsea.*/* .
timestamp=$(cat *.rpt | grep producer_timestamp | awk '{print $2}')

for pattern in _${timestamp} .${timestamp}; do
    find . -name "*${pattern}*" | sed "s|^\./||" | while read -r f; do
        mv $f ${f//$pattern/}
    done
done
sed -i.bak "s/[_\.]$timestamp//g" *.rpt *.html && rm *.bak

# Prefix files that currently lack it
ls -p | grep -v / | grep -v "host_therapy.combined_rat." | while read -r f; do
    mv $f host_therapy.combined_rat.${f}
    sed -i.bak "s/$f/host_therapy.combined_rat.${f}/g" *.rpt *.html && rm *.bak
done

cat <<-END_VERSIONS > versions.yml
"NFCORE_DIFFERENTIALABUNDANCE:DIFFERENTIALABUNDANCE:GSEA_GSEA":
    gsea: 4.3.2
END_VERSIONS

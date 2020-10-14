
for fn in $(ls /share/fsmresfiles/pd_project/wes/annovar_vcf/*.vcf.gz); do
    fout="/share/fsmresfiles/pd_project/wes/annovar_qualityfilter_AND/PPMI_SI_${fn:54:4}_filter.hg19_multianno.vcf.gz"
    bcftools filter --include "FMT/GQ>=10 & FORMAT/DP>=10 & QUAL>=10" --output "${fout}" --output-type z "${fn}"
    echo "Finished for subject ${fn:54:4}"
done

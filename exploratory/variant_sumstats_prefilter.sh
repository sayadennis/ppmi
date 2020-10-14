
dpre="/share/fsmresfiles/pd_project/wes/annovar_vcf/*.vcf.gz"

for fn in $(ls $dpre); do
    echo "##${fn:46:12}"
    bcftools plugin counts ${fn}
done

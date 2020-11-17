
for fn in $(ls /data/annovar_vcf/*.vcf); do
    bgzip $fn
    tabix -p vcf ${fn}.gz
    echo "Finished for subject ${fn:54:4}"
done


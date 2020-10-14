
dqual="/share/fsmresfiles/pd_project/wes/annovar_qualityfilter_AND/*.gz"

for fn in $(ls $dqual); do
    echo "##${fn:60:12}"
    bcftools plugin counts ${fn}
done

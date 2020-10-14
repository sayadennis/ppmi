
dfunc="/share/fsmresfiles/pd_project/wes/filter_delMuts_ExAC_Sept2020"

for fn in $(ls $dfunc); do
    echo "##${fn:0:12}"
    bcftools plugin counts ${dfunc}/${fn}
done

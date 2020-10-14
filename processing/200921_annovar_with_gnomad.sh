# The purpose of this script is to loop through the filtered VCF files and run annovar with the newly downloaded gnomAD annotation.

din="/share/fsmresfiles/pd_project/wes/nonref_filtered"
dout="/share/fsmresfiles/pd_project/wes/annovar_vcf_gnomad"

for i in $(ls $din); do
    substr=${i:0:12}
    if substr in completed
    then
        out="${substr}_gnom"
        perl /home/srd6051/annovar/table_annovar.pl -thread 4 $din/$i /home/srd6051/annovar/humandb/ -buildver hg19 -out $dout/$out -remove -protocol refGene,knownGene,ensGene,dbnsfp30a,dbnsfp31a_interpro,1000g2015aug_all,exac03,gnomad211_genome,gnomad211_exome -operation g,g,g,f,f,f,f,f,f -nastring . -vcfinput
    fi
done

# The purpose of this script is to loop through the filtered VCF files and run annovar on all of them.

input_dir = "../data/wes/filtered_nonref" # fix so that directory is recognized in shell
output_dir = "../data/wes/annovar_vcf" # same as above
annovar_dir = "~" # the path to where annovar is installed

for i in $(ls $input_dir); do # is the $input_dir variable working? 
    substr=${i:0:12}
    out="${substr}_anv"
    perl $annovar_dir/annovar/table_annovar.pl -thread 6 $input_dir/$i $annovar_dir/annovar/humandb/ -buildver hg19 -out $output_dir/$out -remove -protocol refGene,knownGene,ensGene,snp138,dbnsfp30a,dbnsfp31a_interpro,1000g2015aug_all,exac03 -operation g,g,g,r,f,f,f,f -nastring . -vcfinput
done

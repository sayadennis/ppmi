# The purpose of this script is to iterate through all .txt files which are ANNOVAR outputs and
# eliminate rows that are not deleterious mutations

# define column names, define function, take vcf, filter using function, save it as an rds file in new directory

colnames <- c(
    "CHROM", "POS", "ID", "REF", "ALT", 
    "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene", "Func.knownGene", "Gene.knownGene", "GeneDetail.knownGene", "ExonicFunc.knownGene", "AAChange.knownGene", 
    "Func.ensGene", "Gene.ensGene", "GeneDetail.ensGene", "ExonicFunc.ensGene", "AAChange.ensGene", "snp138", "SIFT_score", "SIFT_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_pred",
    "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred", "LRT_score", "LRT_pred", "MutationTaster_score", "MutationTaster_pred", "MutationAssessor_score", "MutationAssessor_pred", "FATHMM_score", "FATHMM_pred",
    "PROVEAN_score", "PROVEAN_pred", "VEST3_score", "CADD_raw", "CADD_phred", "DANN_score", "fathmm-MKL_coding_score", "fathmm-MKL_coding_pred", "MetaSVM_score", "MetaSVM_pred", 
    "MetaLR_score", "MetaLR_pred", "integrated_fitCons_score", "integrated_confidence_value", "GERP++_RS", "phyloP7way_vertebrate", "phyloP20way_mammalian", "phastCons7way_vertebrate", "phastCons20way_mammalian", "SiPhy_29way_logOdds",
    "Interpro_domain", "1000g2015aug_all", "ExAC_ALL", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS",
    "GT", "QUAL", "DP", "chr", "pos", "name", "ref", "alt", "qual", "FILTER", "INFO", "FORMAT", "FORMAT_INFO"
)

delMuts <- function(va) {
    delExonicFuncMuts = c("frameshift insertion", "frameshift deletion", "stopgain", "stoploss")
    delFuncMuts = c("exonic;splicing", "splicing", "ncRNA_splicing")
    
    ftidx = (va$DP >= 10 & # add "not is.na" in front of this ??
             (va$ExonicFunc.refGene %in% delExonicFuncMuts |
              va$ExonicFunc.ensGene %in% delExonicFuncMuts |
              va$ExonicFunc.knownGene %in% delExonicFuncMuts |
              va$Func.refGene %in% delFuncMuts |
              va$Func.ensGene %in% delFuncMuts |
              va$Func.knownGene %in% delFuncMuts))
    return (ftidx)
}

for (file in list.files(path="../data/wes/annovar_vcf", pattern="PPMI_SI_[0-9][0-9][0-9][0-9]_anv.hg19_multianno.txt", full.names=FALSE)) {
    va <- read.table(file, skip=1, header=F, sep="\t", quote="")
    colnames(va) = colnames

    new_va <- va[delMuts(va),]

    saveRDS(new_va, file = paste("../data/wes/filtered_delMuts/", substr(file, 1, 12), "_delmuts.rds", sep=""))
}

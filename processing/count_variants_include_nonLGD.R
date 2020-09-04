# The purpose of this script is to generate a variant count matrix from the PPMI WES VCF files.
# The variant count matrix should include both deleterious (likely-gene-disruptive / LGD) variants and non-LGD variants.

library(VariantAnnotation)

"""
Tasks:
1. loop through VCFs in the input directory 
2. loop through lines of VCF
3. retrieve gene name, snp_ct and indel_ct (thus all_ct)
4. for (mx, ct) in [(allmx, all_ct), (snpmx, snp_ct), (indelmx, indel_ct)], add count to matrix 
5. save allmx, snpmx, indelmx as CSVs 

"""

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

allMuts <- function(va) {
    ftidx <- (va$DP >= 10)
    return (ftidx)
}

# write function to retrieve gene name, snp_ct, and indel_ct
getCounts <- function(vcfrow) {
    gene_name <- refGene.Gene(vcfrow)
    snp_ct <- 0
    indel_ct <- 0
    gt <- geno(vcfrow)$GT # gt = "0/0", "1/0", "2/2" etc.
    # separate genotypes?? 
}


for (file in list.files(path="../data/wes/annovar_vcf", pattern="PPMI_SI_[0-9][0-9][0-9][0-9]_anv.hg19_multianno.txt", full.names=FALSE)) {
    va <- read.table(file, skip=1, header=F, sep="\t", quote="")
    colnames(va) = colnames

    new_va <- va[allMuts(va),]

    for (i in range(length(new_va))) {
        gene_name <- info(new_va)[i,]$Gene.refGene

        snp_ct <- 0
        indel_ct <- 0

        gt <- geno(new_va)$GT[i]
        gtvals <- c(as.integer(gt[1]), as.integer(gt[3]))

        if (width(rowRanges(vcf)$REF[i]) == 1)
    }

    saveRDS(new_va, file = paste("../data/wes/filtered_delMuts/", substr(file, 1, 12), "_delmuts.rds", sep=""))
}

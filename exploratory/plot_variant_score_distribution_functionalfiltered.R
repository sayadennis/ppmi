
library("VariantAnnotation")

vcf <- readVcf("/share/fsmresfiles/pd_project/wes/filter_delMuts_ExAC_Sept2020/PPMI_SI_3000_fullfilter.hg19_multianno.gz", "hg19")

score_dict <- list(
    list("SIFT_score", info(vcf)$SIFT_score),
    list("Polyphen2_HDIV_score", info(vcf)$Polyphen2_HDIV_score),
    list("Polyphen2_HVAR_score", info(vcf)$Polyphen2_HVAR_score),
    list("LRT_score", info(vcf)$LRT_score),
    list("MutationTaster_score", info(vcf)$MutationTaster_score),
    list("MutationAssessor_score", info(vcf)$MutationAssessor_score),
    list("FATHMM_score", info(vcf)$FATHMM_score),
    list("PROVEAN_score", info(vcf)$PROVEAN_score),
    list("VEST3_score", info(vcf)$VEST3_score),
    list("CADD_phred", info(vcf)$CADD_phred),
    list("DANN_score", info(vcf)$DANN_score),
    list("MetaSVM_score", info(vcf)$MetaSVM_score),
    list("MetaLR_score", info(vcf)$MetaLR_score),
    list("integrated_fitCons_score", info(vcf)$integrated_fitCons_score)
)

for (score in score_dict) {
    sn <- score[[1]]
    s_ls <- score[[2]]
    s_vec <- c(as.numeric(s_ls))
    missrate <- sum(is.na(s_vec))*100 / length(s_vec)
    s_vec <- s_vec[!is.na(s_vec)]
    jpeg(paste0("/home/srd6051/pd_project/", sn, "_patid3000_funcfilter_hist", ".jpg"))
    hist(
        s_vec,
        main=paste0("Histogram of ", sn, " but with ", formatC(missrate, digits=3), "% missing"),
        xlab="Score values",
        ylab="Counts",
        breaks=20
    )
    dev.off()
}

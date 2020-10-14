
library("VariantAnnotation")

dn <- "/share/fsmresfiles/pd_project/wes/annovar_qualityfilter_AND"

exacFilter <- function(vcf) {
    as.vector(as.numeric(info(vcf)$ExAC_ALL) < 0.9)
}

delMuts <- function(vcf) {
    delExonicFuncMuts = c("frameshift insertion", "frameshift deletion", "stopgain", "stoploss")
    delFuncMuts = c("exonic;splicing", "splicing", "ncRNA_splicing")
    unlist(as.vector(
        info(vcf)$ExonicFunc.refGene %in% delExonicFuncMuts |
        info(vcf)$ExonicFunc.ensGene %in% delExonicFuncMuts |
        info(vcf)$ExonicFunc.knownGene %in% delExonicFuncMuts |
        info(vcf)$Func.refGene %in% delFuncMuts |
        info(vcf)$Func.ensGene %in% delFuncMuts |
        info(vcf)$Func.knownGene %in% delFuncMuts
        ))
}

filters <- FilterRules(list(delMuts=delMuts, exacFilter=exacFilter))

completed <- list()
for (fn in list.files(path="/share/fsmresfiles/pd_project/wes/filter_delMuts_ExAC_Sept2020", full.names=FALSE)) {
    completed <- append(completed, substr(fn, 9, 12))
}

for (fn in list.files(path=dn, pattern="PPMI_SI_[0-9][0-9][0-9][0-9]_filter.hg19_multianno.vcf.gz$", full.names=TRUE)) {
    patid <- substr(fn, 69, 72)
    if (patid %in% completed) {
        next
    } else {
        destfile <- paste0("/share/fsmresfiles/pd_project/wes/filter_delMuts_ExAC_Sept2020/PPMI_SI_", patid, "_fullfilter.hg19_multianno.gz")
        print(paste0("Filtering for patient ", patid))
        filterVcf(
            fn,
            "hg19",
            destfile,
            filters=filters,
            verbose=TRUE
        )
    }
}

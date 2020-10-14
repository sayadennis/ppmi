
library("VariantAnnotation")

din <- "/share/fsmresfiles/pd_project/wes/annovar_vcf"
dout <- "/share/fsmresfiles/pd_project/wes/annovar_qualityfilter_VariantAnnotation"

## Define filters 

qualFilter <- function(vcf) {
    as.vector((rowRanges(vcf)$QUAL >= 10) & (geno(vcf)$GQ >= 10) & (geno(vcf)$DP >= 10))
}

## Wrap filters together with FilterRules
filters <- FilterRules(list(qualFilter=qualFilter))

## make a list of patient IDs whose filtering is already completed to avoid refiltering them
completed <- list()
for (fn in list.files(path=dout, full.names=FALSE)) {
    completed <- append(completed, substr(fn, 9, 12))
}

## loop through files and filter them 
for (fn in list.files(path=din, pattern="PPMI_SI_[0-9][0-9][0-9][0-9]_anv.hg19_multianno.vcf.gz$", full.names=TRUE)) {
    patid <- substr(fn, 55, 58)
    if (patid %in% completed) {
        next
    } else {
        destfile <- paste0(dout, "/", "PPMI_SI_", patid, "_qualityfilter.hg19_multianno.vcf")
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

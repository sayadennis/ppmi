
library("VariantAnnotation")

vcf <- readVcf("/share/fsmresfiles/pd_project/wes/annovar_qualityfilter_AND/PPNI_SI_3000_filter.hg19_multianno.vcf.gz", "hg19")

ref_ex <- info(vcf)$ExonicFunc.refGene
ref <- info(vcf)$Func.refGene

known_ex <- info(vcf)$ExonicFunc.knownGene
known <- info(vcf)$Func.knownGene

ens_ex <- info(vcf)$ExonicFunc.ensGene
ens <- info(vcf)$Func.ensGene

ex_agreement <- 0
ex_priority <- 0
ex_conflict <- 0

agreement <- 0
priority <- 0
conflict <- 0

for (i in 1:length(ref_ex)) {
    # exonic function
    ls <- list(ref_ex[[i]], known_ex[[i]], ens_ex[[i]])
    if (length(unique(ls)) == 1) {
        ex_agreement <- ex_agreement + 1
    } else if (length(unique(ls)) == 2) {
        # if the two are "." or "unknown" --> priority
        if (sum(ls == "." | ls == "unknown") == 2) {
            ex_priority <- ex_priority + 1
        } else if (sum(ls == "." | ls == "unknown") == 1) {
            other_two <- ls[(ls != ".") & (ls != "unknown")] # a list of the two elements that are not "." or "unknown" 
            if (other_two[[1]] == other_two[[2]]) {
                ex_priority <- ex_priority + 1
            } else {
                ex_conflict <- ex_conflict + 1
            }
        } else if (sum(ls == "." | ls == "unknown") == 0) {
            conflict <- conflict + 1
        } else{
            paste("Problem values:", ls[[1]], ls[[2]], ls[[3]])
        }
    } else if (length(unique(ls)) == 3) {
        ex_conflict <- ex_conflict + 1
    } else {
        paste("Problem values: ", ls[[1]], ls[[2]], ls[[3]])
    }

    # general function
    ls <- list(ref[[i]], known[[i]], ens[[i]])
    if (length(unique(ls)) == 1) {
        agreement <- agreement + 1
    } else if (length(unique(ls)) == 2) {
        # if the two are "." or "unknown" --> priority
        if (sum(ls == "." | ls == "unknown") == 2) {
            priority <- priority + 1
        } else if (sum(ls == "." | ls == "unknown") == 1) {
            other_two <- ls[(ls != ".") & (ls != "unknown")] # a list of the two elements that are not "." or "unknown" 
            if (other_two[[1]] == other_two[[2]]) {
                priority <- priority + 1
            } else {
                conflict <- conflict + 1
            }
        } else if (sum(ls == "." | ls == "unknown") == 0) {
            conflict <- conflict + 1
        } else {
            paste("Problem values: ", ls[[1]], ls[[2]], ls[[3]])
        }
    } else if (length(unique(ls)) == 3) {
        conflict <- conflict + 1
    } else {
        paste("Problem values: ", ls[[1]], ls[[2]], ls[[3]])
    }
}

print("")
print("Final counts for exonic function")
print(paste("Conflict:", ex_conflict))
print(paste("Priority:", ex_priority))
print(paste("Agreement:", ex_agreement))
print("")
print("Final counts for general function")
print(paste("Conflict:", conflict))
print(paste("Priority:", priority))
print(paste("Agreement:", agreement))
print("")

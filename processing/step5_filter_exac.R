# The purpose of this script is to remove lines from the PPMI WES tables where the ExAC values are >0.9

read_dir = "../data/wes/filtered_delMuts"
write_dir = "../data/wes/filtered_exac"

for (file in list.files(path=read_dir, pattern="*.rds", full.names=TRUE)) {
    table <- readRDS(file)
    colnames(table)[66] <- "DS"
    exac <- table$ExAC_ALL # ExAC frequencies - but in factors.
    levels(exac)[levels(exac)=="."] <- "0" # changing missing values to zero for the purpose of creating logical vector in line 12
    new_exac <- as.numeric(levels(exac))[exac] # changing factors to numeric and saving to new vector
    below_thres <- new_exac<0.9 # logical vector of TRUE or FALSE (no NA since all of "." is changed to 0 in line 10)
    new_table <- table[below_thres,] # new table with only rows that have ExAC_ALL values smaller than 0.9 OR missing.
    write.table(new_table, file = paste(substr(file, 1, 28), "filtered_exac/", substr(file, 46, 57), "_exac.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
}

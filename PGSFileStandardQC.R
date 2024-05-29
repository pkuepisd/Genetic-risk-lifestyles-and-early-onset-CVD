options(stringsAsFactors=F)

library("optparse")
library("data.table")
reverse_allele <- function(allele) {
    a <- gsub("A","1",allele)
    a <- gsub("T","2",a)
    a <- gsub("C","3",a)
    a <- gsub("G","4",a)
    a <- gsub("1","T",a)
    a <- gsub("2","A",a)
    a <- gsub("3","G",a)
    a <- gsub("4","C",a)
    return(a)
}

#------------------------------------------------------------------------------
## set list of cmd line arguments
option_list <- list(
  make_option("--pgsf", type="character", default="",
    help="PGS Standard file"),
  make_option("--target", type="character", default="",
    help="info.table of target dataset"),
  make_option("--delete_ambiguous", type="logical", default=TRUE,
    help="QC2: Deleting ambiguous SNP (A/T)"),
  make_option("--delete_indel", type="logical", default=TRUE,
    help="QC3: Deleting indel variants"),
  make_option("--ck_info", type="logical", default=TRUE,
    help="QC4: Deleting low imputation quality variants"),
  make_option("--INFO", type="numeric", default=0,
    help="low imputation INFO threshold"),
  make_option("--ck_freq", type="logical", default=TRUE,
    help="QC5: Deleting variants with low minor allele frequency"),
  make_option("--FREQ", type="numeric", default=1,
    help="MAF threshold"),
  make_option("--ck_hwe", type="logical", default=TRUE,
    help="QC6: Deleting variants not satisfying the Hardy-Weinberg equilibrium"),
  make_option("--HWP", type="numeric", default=0,
    help="P values threshold of Hardy-Weinberg equilibrium test")
  )
## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
#------------------------------------------------------------------------------

PGS_FILE=opt$pgsf
INFO_FILE=opt$target
delete_ambiguous=opt$delete_ambiguous
delete_indel=opt$delete_indel
ck_info=opt$ck_info
INFO=opt$INFO
ck_freq=opt$ck_freq
FREQ=opt$FREQ
ck_hwe=opt$ck_hwe
HWP=opt$HWP

#===============================================================================
pgs <- fread(PGS_FILE, header=T)
N0 = nrow(pgs)

info <- fread(INFO_FILE, header=T)

#===============================================================================
info_sub <- subset(info, BP %in% pgs$BP)
df <- merge(pgs, info_sub, by=c("CHR", "BP"))
df <- as.data.frame(df)

ck1 <- ifelse( (df[,"EA"]==df[,"A1"] & df[,"NEA"]==df[,"A2"]) | 
               (df[,"EA"]==df[,"A2"] & df[,"NEA"]==df[,"A1"]), TRUE, FALSE)
df1 <- subset(df, ck1)

df2 <- subset(df, !ck1)
if (nrow(df2) > 0) {
    fEA <- reverse_allele(df2$EA)
    fNEA <- reverse_allele(df2$NEA)
    ck2 <- ifelse( (df2$A1==fEA & df2$A2==fNEA) | 
                (df2$A1==fNEA & df2$A2==fEA), TRUE, FALSE)
    df2_1 <- subset(df2, ck2)
    if (nrow(df2_1) > 0) {
        df2_1$EA <- reverse_allele(df2_1$EA)
        df2_1$NEA <- reverse_allele(df2_1$NEA)
    }
    df2_2 <- subset(df2, !ck2)
}

if (nrow(df2) > 0) {
    df <- rbind(df1, df2_1)
}
N1 = nrow(df)  
D1 = N0 - N1
#-------------------------------------------------------------------------------

if (delete_ambiguous) {
    special <- ifelse( (df[,"EA"]=="C" & df[,"NEA"]=="G") |
                   (df[,"EA"]=="G" & df[,"NEA"]=="C") |
                   (df[,"EA"]=="A" & df[,"NEA"]=="T") |
                   (df[,"EA"]=="T" & df[,"NEA"]=="A"), TRUE, FALSE)
    df <- subset(df, !special)
}
N2 = nrow(df)
D2 = N1 - N2
#-------------------------------------------------------------------------------

if (delete_indel) {
    nchar <- nchar(df[,"EA"]) + nchar(df[,"NEA"])
    indel <- ifelse(nchar>2, TRUE, FALSE)
    df <- subset(df, !indel)
}
N3 = nrow(df)
D3 = N2 - N3
#-------------------------------------------------------------------------------

if (ck_info) {
    lowinfo <- ifelse(df[,"info"] < INFO, TRUE, FALSE)
    df <- subset(df, !lowinfo)
}
N4 = nrow(df)
D4 = N3 - N4
#-------------------------------------------------------------------------------

if (ck_freq) {
    lowfreq <- ifelse(df[,"A1_freq"] < FREQ | df[,"A1_freq"] > 1-FREQ, TRUE, FALSE)
    df <- subset(df, !lowfreq)
}
N5 = nrow(df)
D5 = N4 - N5
#-------------------------------------------------------------------------------

if (ck_hwe) {
    viohwe <- ifelse(df[,"P_HW"] < HWP, TRUE, FALSE)
    df <- subset(df, !viohwe)
}
N6 = nrow(df)
D6 = N5 - N6

index = c(0:6)
qclab = c("Original", "Allele mismatch", "Ambiguous SNP", "Ins/Del variant", 
    paste0("Imputation info < ", INFO), 
    paste0("MAF < ", FREQ), 
    paste0("P_HWE < ", HWP) )
qcon = c(NA, TRUE, delete_ambiguous, delete_indel, 
    delete_indel, ck_freq, ck_hwe)
del = c(NA, D1, D2, D3, D4, D5, D6)
keep = c(N0, N1, N2, N3, N4, N5, N6)

qclog <- data.frame(index=index, qclab=qclab, qcon=qcon, del=del, keep=keep)
print(qclog)

write.csv(qclog, file=paste0(PGS_FILE, ".stdqc.log"), row.names=F, quote=F)

df$final_EA <- ifelse(df$BETA>=0, df[,"EA"], df[,"NEA"])
df$final_NEA <- ifelse(df$BETA>=0, df[,"NEA"], df[,"EA"])
df$final_BETA <- abs(df$BETA)
out <- subset(df, select=c("SNP", "CHR", "BP", "final_EA", "final_NEA", "final_BETA"))
colnames(out) <- c("SNP", "CHR", "BP", "EA", "NEA", "BETA") 

DIR = paste0(PGS_FILE, ".stdqc_chr")
try ( dir.create(DIR), silent=TRUE)
setwd(DIR)

for (i in 1:22) {
    outchr <- subset(out, CHR == i)
    if (nrow(outchr)>0) {
        fwrite(outchr, file=paste0("chr", i, ".txt"), sep=" ", quote=F, na="NA")
    }
}


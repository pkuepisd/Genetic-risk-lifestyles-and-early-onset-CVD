# QC of the GWAS data (R file)
options(stringsAsFactors=F)

## loading packages
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
message("### Self-defined function 'reverse_allele' is loaded.")
  
#------------------------------------------------------------------------------
## set list of cmd line arguments
option_list <- list(
  make_option("--ssf", type="character", default="",
    help="Standard GWAS summary statistics file"),
  make_option("--ssinfo", type="numeric", default=0,
    help="low imputation INFO threshold in base dataset"),
  make_option("--ssmaf", type="numeric", default=0,
    help="MAF threshold in base dataset"),
  make_option("--target", type="character", default="",
    help="info.table of target dataset"),
  make_option("--delete_ambiguous", type="logical", default=TRUE,
    help="Deleting ambiguous SNP (A/T)"),
  make_option("--delete_indel", type="logical", default=TRUE,
    help="Deleting indel variants"),
  make_option("--ck_info", type="logical", default=TRUE,
    help="Deleting low imputation quality variants in target dataset"),
  make_option("--INFO", type="numeric", default=0,
    help="low imputation INFO threshold in target dataset"),
  make_option("--ck_freq", type="logical", default=TRUE,
    help="Deleting variants with low minor allele frequency in target dataset"),
  make_option("--FREQ", type="numeric", default=1,
    help="MAF threshold in target dataset"),
  make_option("--ck_hwe", type="logical", default=TRUE,
    help="Deleting variants not satisfying the Hardy-Weinberg equilibrium in target dataset"),
  make_option("--HWP", type="numeric", default=0,
    help="P values threshold of Hardy-Weinberg equilibrium test in target dataset"),
  make_option("--qctype", type="numeric", default=0,
    help="QC type")
  )
## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
#------------------------------------------------------------------------------

SS_FILE = opt$ssf
SSINFO = opt$ssinfo
SSMAF = opt$ssmaf
INFO_FILE = opt$target
delete_ambiguous = opt$delete_ambiguous
delete_indel = opt$delete_indel
ck_info = opt$ck_info
INFO = opt$INFO
ck_freq = opt$ck_freq
FREQ = opt$FREQ
ck_hwe = opt$ck_hwe
HWP = opt$HWP
qctype = opt$qctype

#===============================================================================
message(paste0("### Reading base dataset: "), SS_FILE)
ssf <- fread(SS_FILE, header=T)
N0 = nrow(ssf)

message(paste0("### Reading info.table of target dataset: "), INFO_FILE)
info <- fread(INFO_FILE, header=T)

message("### QC ing ... ")
#===============================================================================
ssf$INFO <- ifelse(is.na(ssf$INFO), 1, ssf$INFO)   
sck1 <- ifelse( ssf$INFO < SSINFO, FALSE, TRUE )
ssf <- subset(ssf, sck1)
n1 = nrow(ssf)
d1 = N0 - n1

sck2 <- ifelse( ssf$MAF < SSMAF, FALSE, TRUE )
ssf <- subset(ssf, sck2)
n2 = nrow(ssf)
d2 = n1 - n2
#===============================================================================
info_sub <- subset(info, BP %in% ssf$BP)
df <- merge(ssf, info_sub, by=c("CHR", "BP"))
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
}

if (nrow(df2) > 0) {
    df <- rbind(df1, df2_1)
}
N1 = nrow(df) 
D1 = n2 - N1
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

index = c(0:8)
qclab = c("Original", 
    paste0("Base imputation info < ", SSINFO), 
    paste0("Base MAF < ", SSMAF),
    "Not found in target", 
    "Ambiguous SNP", 
    "Ins/Del variant", 
    paste0("Target imputation info < ", INFO), 
    paste0("Target MAF < ", FREQ), 
    paste0("Target P_HWE < ", HWP) )
qcon = c(NA, TRUE, TRUE, 
    TRUE, delete_ambiguous, delete_indel, 
    delete_indel, ck_freq, ck_hwe)
del = c(NA, d1, d2, D1, D2, D3, D4, D5, D6)
keep = c(N0, n1, n2, N1, N2, N3, N4, N5, N6)

qclog <- data.frame(index=index, qclab=qclab, qcon=qcon, del=del, keep=keep)
print(qclog)

message(paste0("### QC log file is exported as: ", SS_FILE, ".stdqc",qctype,".log"))
write.csv(qclog, file=paste0(SS_FILE, ".stdqc",qctype,".log"), row.names=F, quote=F)

df$final_EA <- ifelse(df$BETA>=0, df[,"EA"], df[,"NEA"])
df$final_NEA <- ifelse(df$BETA>=0, df[,"NEA"], df[,"EA"])
df$final_BETA <- abs(df$BETA)
df$final_EAF <- ifelse(df$BETA>=0, df[,"EAF"], 1-df[,"EAF"])

out <- subset(df, select=c("SNP", "CHR", "BP", "final_EA", "final_NEA", "final_EAF", "final_BETA", "SE", "P"))
colnames(out) <- c("SNP", "CHR", "POS", "A1", "A2", "A1FREQ", "EFF", "SE", "PVAL")

if (length(table(duplicated(out$SNP)))==2) {
    message("There are duplicates in SSF! Stop generated")
} else {
    fwrite(out, file=paste0(SS_FILE, ".stdqc",qctype), sep=" ", quote=F, na="NA")
    message(paste0("### Standard QCed summary statistics file is exported as: ", SS_FILE, ".stdqc",qctype))

    DIR = paste0(SS_FILE, ".stdqc",qctype,"_chr")
    try ( dir.create(DIR), silent=TRUE)
    setwd(DIR)

    message(paste0("### Standard QCed summary statistics files by chromosome are also exported to", SS_FILE, ".stdqc",qctype,"_chr/"))

    for (i in 1:22) {
        outchr <- subset(out, CHR == i)
        if (nrow(outchr)>0) {
            fwrite(outchr, file=paste0("chr", i, ".txt"), sep=" ", quote=F, na="NA")
        }
    }
}
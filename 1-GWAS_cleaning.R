# Clean the GWAS summary data
library("data.table")
library("plyr")
source("checkna.r")
source("reverse_allele.r") 
BASE = "~/PhD/Project1/GWASdata/IS_EAS/GCST90018644_buildGRCh37.tsv"
OUT = "~/PhD/Project1/GWASdata/IS_EAS/workdir/IS_EAS.org"

ba <- fread(BASE, header=TRUE)
ba$MAF <- ifelse(ba$effect_allele_frequency<0.5, ba$effect_allele_frequency, 1-ba$effect_allele_frequency)
ba$N <- 174686
ba$INFO <- NA
ba <- ba[, c("chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "MAF", "beta", "standard_error", "p_value", "N", "INFO")]
colnames(ba) <- c("CHR", "BP", "EA", "NEA", "EAF", "MAF", "BETA", "SE", "P", "N", "INFO")
ba <- arrange(ba, CHR, BP)
ba$EA <- toupper(ba$EA)
ba$NEA <- toupper(ba$NEA)

fwrite(ba, file=OUT, sep=" ", quote=F, na="NA")

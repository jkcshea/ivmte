# This file creates a smaller version of some of the data used in
# Angrist and Evans (1998, American Economic Review).
#
# The original data in Stata format is available here:
#
# http://sites.bu.edu/ivanf/files/2014/03/m_d_806.dta_.zip
#
# That file is cleaned using this Stata code:
#
# http://sites.bu.edu/ivanf/files/2014/03/code.zip
#
# We have run this and saved it in .csv format, which this script downloads.
# The file is around 80 MB, so not included with the package.
#
# This script creates smaller version by removing unneeded columns.
rm(list = ls())

temp = "AE.csv"
download.file(
    "http://ivmte.data.s3-website.us-east-2.amazonaws.com/AE-1980-AllWomen.csv",
    temp
)
AE <- read.csv(temp)
unlink(temp) # delete the file

AE <- subset(AE, agefstm >= 20)
AE <- subset(AE, select=c(workedm, hourswm, morekids, samesex, YOBM,
                          blackm, hispm, othracem))
colnames(AE) <- c("worked", "hours", "morekids", "samesex", "yob",
                  "black", "hisp", "other")
rownames(AE) <- NULL
library("devtools")
use_data(AE, overwrite = TRUE)

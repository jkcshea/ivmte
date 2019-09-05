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
# This script creates smaller version by removing unneeded columns.
rm(list = ls())

#### WAITING FOR WEB SPACE FROM U CHICAGO

#download.file("http://sites.bu.edu/ivanf/files/2014/03/m_d_806.dta_.zip",temp)
#aedata <- read_dta(unz(temp, "m_d_806.dta"))
#unlink(temp)

AE = read.csv("1980-AllWomen.csv")
AE <- subset(AE, agefstm >= 20)
AE <- subset(AE, select=c(workedm, hourswm, morekids, samesex, YOBM))
colnames(AE) <- c("worked", "hours", "morekids", "samesex", "yob")
library("devtools")
use_data(AE, overwrite = TRUE)

# --------------------------------------------------------------------------------------------- #
# --------------------------------------- ALDH2 overview -------------------------------------- #
# programmer: Zhengting He
# date: May 12th, 2021
# --------------------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------------------- #
# --------------------- overlap-correct-analysis without managing overlap --------------------- #
# -- aim: meta-analysis against all included primary studies by genetics models and outcomes -- #
# --------------------------------------------------------------------------------------------- #

library(metafor)
library(dplyr)
library(ggplot2)

setwd("F:/Box Sync/Archives2020LLY/Zhengting/Duke Kunshan University Intern (zh133@duke.edu)/6 ALDH2 SR/data synthesis/overlap-correct-analysis")

# --------------------------------------------------------------------------------------------- #
# -------------------------- using data provided by original authors -------------------------- #
# --------------------------------------------------------------------------------------------- #

data <- read.csv("data/data extraction table (primary study)_original author's data.csv", encoding = "UTF-8")
colnames(data)[1] <- "Group"

# convert effect size and confidential intervial to log value, generate SE_lnOR
# formula is generated according to https://stats.stackexchange.com/questions/156597/how-to-calculate-se-of-an-odds-ratio
# CI = exp(lnOR +- 1.96 * sqrt(1/a + 1/b + 1/c + 1/d))
# ln(CI_lower) = ln(OR) - 1.96 * SE_lnOR
# ln(CI_upper) = ln(OR) + 1.96 * SE_lnOR
# SE_lnOR = 1/(2*1.96) * (ln(CI_upper) - ln(CI_lower))
data$lnOR <- log(data$OR)
data$SE <- 1/3.92 * (log(data$CI2) - log(data$CI1))
data$lnOR <- round(data$lnOR, 2) # round lnOR to 2-digit floating number
data$SE <- round(data$SE, 2) # round SE to 2-digit floating number

# separate ma_data by genetic model, generate ma_data_"model"
model <- c("allelic", "dominant", "heterozygous", "homozygous", "recessive")
data_model <- list()
for(m in model) {
  print(m)
  data_model[[m]] <- filter(data, model == m) # data_model[[m]] = data_"model"
}
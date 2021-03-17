# --------------------------------------------------------------------------------------------- #
# --------------------------------------- ALDH2 overview -------------------------------------- #
# programmer: Zhengting He
# date: Feburary 4th, 2021
# code is based on original code "ALDH2 overview code v2.R" by Zhengting He
# --------------------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------------------- #
# ------------------------ meta-meta-analysis without managing overlap ------------------------ #
# --- aim: meta-analysis against all included meta-analyses by genetics models and outcomes --- #
# --------------------------------------------------------------------------------------------- #

library(metafor)
library(dplyr)
library(ggplot2)

setwd("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/data")

# --------------------------------------------------------------------------------------------- #
# -------------------------- using data provided by original authors -------------------------- #
# --------------------------------------------------------------------------------------------- #

ma_data <- read.csv("data extraction table (meta-analysis)_original authors' data.csv")

# uniform all observations using "ALDH2*2" as exposure and "ALDH2*1" as comparison
ma_data$CI1_copy <- ma_data$CI1
for(i in 1:nrow(ma_data)) {
  if(ma_data$Exposure[i] == "ALDH2*1" & ma_data$Comparison[i] == "ALDH2*2") {
    ma_data$OR[i] <- 1/ma_data$OR[i]
    ma_data$CI1[i] <- 1/ma_data$CI2[i]
    ma_data$CI2[i] <- 1/ma_data$CI1_copy[i]
    ma_data$Exposure[i] <- "ALDH2*2"
    ma_data$Comparison[i] <- "ALDH2*1"
  }
}
ma_data <- subset(ma_data, select = - CI1_copy)

write.csv(ma_data, file = "F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/data extraction table (meta-analysis)_original authors' data_uniform exposure and comparison.csv")

# subset ma_data to measurement of effect size as "Odds Ratio"
# Only studies which measurement of effect size is "Odds Ratio" are eligable for synthesizing
ma_data <- filter(ma_data, ES_measurement == "OR")

# convert effect size and confidential intervial to log value, generate SE_lnOR
# formula is generated according to https://stats.stackexchange.com/questions/156597/how-to-calculate-se-of-an-odds-ratio
# CI = exp(lnOR +- 1.96 * sqrt(1/a + 1/b + 1/c + 1/d))
# ln(CI_lower) = ln(OR) - 1.96 * SE_lnOR
# ln(CI_upper) = ln(OR) + 1.96 * SE_lnOR
# SE_lnOR = 1/(2*1.96) * (ln(CI_upper) - ln(CI_lower))
ma_data$lnOR <- log(ma_data$OR)
ma_data$SE <- 1/3.92 * (log(ma_data$CI2) - log(ma_data$CI1))
ma_data$lnOR <- round(ma_data$lnOR, 2) # round lnOR to 2-digit floating number
ma_data$SE <- round(ma_data$SE, 2) # round SE to 2-digit floating number

# separate ma_data by genetic model, generate ma_data_"model"
model <- c("allelic", "dominant", "heterozygous", "homozygous", "recessive")
ma_data_model <- list()
for(m in model) {
  print(m)
  ma_data_model[[m]] <- filter(ma_data, model == m) # ma_data_model[[m]] = ma_data_"model"
}

ma_data$Group <- as.numeric(ma_data$Group)

# --------------------------------------------------------------------------------------------- #
# ----------------- Part 1: present the results of all included meta-analyses ----------------- #
# --------------------------------------------------------------------------------------------- #

# fit random-effect or fix-effect model
ma.model <- list()
for(m in model) {
  ma <- rma(data = ma_data_model[[m]],
            yi = lnOR,
            sei = SE,
            method = "DL",
            slab = paste(ma_data_model[[m]]$Author, ma_data_model[[m]]$Year, sep = ", "),
            digits = 2,
            level = 95)
  if(ma$I2 >= 50) {
    ma.model[[m]] <- ma       # if I2 >= 50, fit random-effect model (DL method)
  } else {                    # if I2 < 50, fit fix-effect model (DL method)
    ma.model[[m]] <- rma(data = ma_data_model[[m]],
                         yi = lnOR,
                         sei = SE,
                         method = "FE",
                         slab = paste(ma_data_model[[m]]$Author, ma_data_model[[m]]$Year, sep = ", "),
                         digits = 2,
                         level = 95)
    
  } # ma.model[[m]] = ma."model"
}

# ----------------------------------- Part 1.1: forest plot ----------------------------------- #
# refer to http://www.metafor-project.org/doku.php/plots:forest_plot_with_subgroups

# allelic model
pdf("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/forestplot/ma.allelic.forestplot.pdf", width = 7.23, height = 6.39)
par(mar=c(2, 1, 1, 1)) # decrease margins so the full space is used
# set up forest plot, the 'rows' argument is used to specify in which rows the outcomes will be plotted
forest(x = ma.model[["allelic"]], addfit=FALSE, showweights=FALSE,
       xlim=c(- 8, 4),
       ilab = cbind(ma_data_model[["allelic"]]$Outcome),
       ilab.xpos = c(-3), cex=0.55,
       xlab = "log Odds Ratio", mlab = "", psize = 1,
       order = order(ma_data_model[["allelic"]]$Group), ylim=c(1, 67), rows = c(63:55, 52:15, 12:8, 5:4, 1)) # ylim_forest_upper = ylim_subgroup_upper + 4
# add headings
par(font=2) # switch to bold font
text(c(-3), 66, c("Outcome of Interest"), cex=0.65) # Y-axis is ylim_forest_upper - 1
text(c(-7.05), 66, c("Author(s) and Year"), cex=0.65) # Y-axis is ylim_forest_upper - 1
text(c(2.82), 66, c("log Odds Ratio[95% CI]"), cex=0.65) # Y-axis is ylim_forest_upper - 1
# add text for subgroups
text(c(-6.335), 64.2, c("Cardio-cerebral Vascular Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-7.56), 53.2, c("Cancer"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-6.59), 13.2, c("Alocholic Digestive Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-6.12), 6.2, c("Nervous System Degenerative Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-7.12), 2.2, c("Diabetes Mellitus"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
dev.off()

# dominant model
pdf("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/forestplot/ma.dominant.forestplot.pdf", width = 7.23, height = 7.39)
par(mar=c(2, 1, 1, 1)) # decrease margins so the full space is used
# set up forest plot, the 'rows' argument is used to specify in which rows the outcomes will be plotted
forest(x = ma.model[["dominant"]], addfit=FALSE, showweights=FALSE,
       xlim=c(- 8, 4),
       ilab = cbind(ma_data_model[["dominant"]]$Outcome),
       ilab.xpos = c(-3), cex=0.55,
       xlab = "log Odds Ratio", mlab = "", psize = 1,
       order = order(ma_data_model[["dominant"]]$Group), ylim=c(1, 79), rows = c(73:58, 55:17, 14:11, 8:5, 2:1)) # ylim_forest_upper = ylim_subgroup_upper + 4
# add headings
par(font=2) # switch to bold font
text(c(-3), 78, c("Outcome of Interest"), cex=0.65) # Y-axis is ylim_forest_upper - 1
text(c(-7.05), 78, c("Author(s) and Year"), cex=0.65) # Y-axis is ylim_forest_upper - 1
text(c(2.82), 78, c("log Odds Ratio[95% CI]"), cex=0.65) # Y-axis is ylim_forest_upper - 1
# add text for subgroups
text(c(-6.335), 74.2, c("Cardio-cerebral Vascular Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-7.56), 56.2, c("Cancer"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-6.59), 15.2, c("Alocholic Digestive Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-6.12), 9.2, c("Nervous System Degenerative Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-7.12), 3.2, c("Diabetes Mellitus"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
dev.off()

# heterozygous model
pdf("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/forestplot/ma.heterozygous.forestplot.pdf", width = 7.23, height = 6.39)
par(mar=c(2, 1, 1, 1)) # decrease margins so the full space is used
# set up forest plot, the 'rows' argument is used to specify in which rows the outcomes will be plotted
forest(x = ma.model[["heterozygous"]], addfit=FALSE, showweights=FALSE,
       xlim=c(- 8, 4),
       ilab = cbind(ma_data_model[["heterozygous"]]$Outcome),
       ilab.xpos = c(-3), cex=0.55,
       xlab = "log Odds Ratio", mlab = "", psize = 1,
       order = order(ma_data_model[["heterozygous"]]$Group), ylim=c(1, 69), rows = c(65:57, 54:12, 9:8, 5:4, 1)) # ylim_forest_upper = ylim_subgroup_upper + 4
# add headings
par(font=2) # switch to bold font
text(c(-3), 68, c("Outcome of Interest"), cex=0.65) # Y-axis is ylim_forest_upper - 1
text(c(-7.05), 68, c("Author(s) and Year"), cex=0.65) # Y-axis is ylim_forest_upper - 1
text(c(2.82), 68, c("log Odds Ratio[95% CI]"), cex=0.65) # Y-axis is ylim_forest_upper - 1
# add text for subgroups
text(c(-6.335), 66.2, c("Cardio-cerebral Vascular Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-7.56), 55.2, c("Cancer"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-6.59), 10.2, c("Alocholic Digestive Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-6.12), 6.2, c("Nervous System Degenerative Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-7.12), 2.2, c("Diabetes Mellitus"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
dev.off()

# homozygous model
pdf("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/forestplot/ma.homozygous.forestplot.pdf", width = 7.23, height = 6.39)
par(mar=c(2, 1, 1, 1)) # decrease margins so the full space is used
# set up forest plot, the 'rows' argument is used to specify in which rows the outcomes will be plotted
forest(x = ma.model[["homozygous"]], addfit=FALSE, showweights=FALSE,
       xlim=c(- 8, 4),
       ilab = cbind(ma_data_model[["homozygous"]]$Outcome),
       ilab.xpos = c(-4), cex=0.55,
       xlab = "log Odds Ratio", mlab = "", psize = 1,
       order = order(ma_data_model[["homozygous"]]$Group), ylim=c(1, 69), rows = c(65:57, 54:11, 8, 5:4, 1)) # ylim_forest_upper = ylim_subgroup_upper + 4
# add headings
par(font=2) # switch to bold font
text(c(-4), 68, c("Outcome of Interest"), cex=0.65) # Y-axis is ylim_forest_upper - 1
text(c(-7.05), 68, c("Author(s) and Year"), cex=0.65) # Y-axis is ylim_forest_upper - 1
text(c(2.82), 68, c("log Odds Ratio[95% CI]"), cex=0.65) # Y-axis is ylim_forest_upper - 1
# add text for subgroups
text(c(-6.335), 66.2, c("Cardio-cerebral Vascular Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-7.56), 55.2, c("Cancer"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-6.59), 9.2, c("Alocholic Digestive Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-6.12), 6.2, c("Nervous System Degenerative Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-7.12), 2.2, c("Diabetes Mellitus"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
dev.off()

# recessive model
pdf("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/forestplot/ma.recessive.forestplot.pdf", width = 7.23, height = 6.39)
par(mar=c(2, 1, 1, 1)) # decrease margins so the full space is used
# set up forest plot, the 'rows' argument is used to specify in which rows the outcomes will be plotted
forest(x = ma.model[["recessive"]], addfit=FALSE, showweights=FALSE,
       xlim=c(- 8, 4),
       ilab = cbind(ma_data_model[["recessive"]]$Outcome),
       ilab.xpos = c(-4.4), cex=0.55,
       xlab = "log Odds Ratio", mlab = "", psize = 1,
       order = order(ma_data_model[["recessive"]]$Group), ylim=c(1, 64), rows = c(60:52, 49:11, 8, 5:4, 1)) # ylim_forest_upper = ylim_subgroup_upper + 4
# add headings
par(font=2) # switch to bold font
text(c(-4.4), 63, c("Outcome of Interest"), cex=0.65) # Y-axis is ylim_forest_upper - 1
text(c(-7.05), 63, c("Author(s) and Year"), cex=0.65) # Y-axis is ylim_forest_upper - 1
text(c(2.82), 63, c("log Odds Ratio[95% CI]"), cex=0.65) # Y-axis is ylim_forest_upper - 1
# add text for subgroups
text(c(-6.335), 61.2, c("Cardio-cerebral Vascular Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-7.56), 50.2, c("Cancer"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-6.59), 9.2, c("Alocholic Digestive Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-6.12), 6.2, c("Nervous System Degenerative Diseases"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
text(c(-7.12), 2.2, c("Diabetes Mellitus"), cex=0.65) # Y-axis is ylim_subgroup_upper + 1.2
dev.off()

# --------------------------------- Part 1.2: stop light plot --------------------------------- #

for(m in model) {
  ma_data_model[[m]]$lnCI1 <- ma_data_model[[m]]$lnOR - 1.96 * ma_data_model[[m]]$SE
  ma_data_model[[m]]$lnCI2 <- ma_data_model[[m]]$lnOR + 1.96 * ma_data_model[[m]]$SE
  ma_data_model[[m]]$lnCI1 <- round(ma_data_model[[m]]$lnCI1, 2)
  ma_data_model[[m]]$lnCI2 <- round(ma_data_model[[m]]$lnCI2, 2)
  
  ma_data_model[[m]]$Authoryear <- paste(ma_data_model[[m]]$Author, ma_data_model[[m]]$Year, sep = " ")
  ma_data_model[[m]]$plotlabel <- paste(format(ma_data_model[[m]]$lnOR, digits = 2), " [", format(ma_data_model[[m]]$lnCI1, digits = 2), ", ", format(ma_data_model[[m]]$lnCI2, digits = 2), "]", sep = "")
  
  for(i in 1:nrow(ma_data_model[[m]])) {
    if(ma_data_model[[m]]$lnOR[i]>=0 & ma_data_model[[m]]$lnCI1[i]>=0 & ma_data_model[[m]]$lnCI2[i]>=0) {
      ma_data_model[[m]]$plotlabel[i] <- paste("   ", ma_data_model[[m]]$plotlabel[i], sep = "")
    } else if(ma_data_model[[m]]$lnOR[i]>=0 & ma_data_model[[m]]$lnCI1[i]>=0 & ma_data_model[[m]]$lnCI2[i]<0) {
      ma_data_model[[m]]$plotlabel[i] <- paste("  ", ma_data_model[[m]]$plotlabel[i], sep = "")
    } else if(ma_data_model[[m]]$lnOR[i]>=0 & ma_data_model[[m]]$lnCI1[i]<0 & ma_data_model[[m]]$lnCI2[i]>=0) {
      ma_data_model[[m]]$plotlabel[i] <- paste("  ", ma_data_model[[m]]$plotlabel[i], sep = "")
    } else if(ma_data_model[[m]]$lnOR[i]>=0 & ma_data_model[[m]]$lnCI1[i]<0 & ma_data_model[[m]]$lnCI2[i]<0) {
      ma_data_model[[m]]$plotlabel[i] <- paste(" ", ma_data_model[[m]]$plotlabel[i], sep = "")
    }
  }
  
  ma_data_model[[m]]$no <- c(1:nrow(ma_data_model[[m]]))
  
  for(i in 1:nrow(ma_data_model[[m]])) {
    if(ma_data_model[[m]]$lnCI1[i] > 0) {
      ma_data_model[[m]]$group[i] <- "preventive"
    } else if(ma_data_model[[m]]$lnCI2[i] >= 0) {
      ma_data_model[[m]]$group[i] <- "no effect"
    } else if(ma_data_model[[m]]$lnCI2[i] < 0) {
      ma_data_model[[m]]$group[i] <- "promotive"
    }
  }
  
  pdf(file = paste("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/stoplight plot/trafficlight_", m, ".pdf", sep = ""), width = 20.013, height = 11.377)
  print(
    ggplot(ma_data_model[[m]], aes(x = lnOR, y = no, fill = group)) +
      geom_vline(xintercept = 0) +
      geom_crossbar(aes(xmin = lnCI1, xmax = lnCI2), alpha = 0.3, width = 0.95) +
      scale_fill_manual(values = c("#e7fb05", "#fb0516", "#05fb1c")) +
      geom_text(aes(label = Outcome, x = lnCI2 + 0.27), size = 4) +
      xlab("Log Odds Ratio") +
      scale_y_continuous("Author and year", breaks = c(1:nrow(ma_data_model[[m]])),
                         labels = ma_data_model[[m]]$Authoryear,
                         sec.axis = sec_axis(trans = ~., name = "Log Odds Ratio and Confidential Interval", breaks = c(1:nrow(ma_data_model[[m]])),
                                             labels = ma_data_model[[m]]$plotlabel)) +
      theme(axis.text = element_text(size = 14,color="black"),axis.title = element_text(size = 14,color="black")) +
      guides(fill=FALSE) +
      theme_classic()
  )
  dev.off()
}

# --------------------------------------------------------------------------------------------- #
# -- Part 2: synthesize the ES and CI of included meta-analyses under each outcome and model -- #
# --------------------------------------------------------------------------------------------- #

# ------------------ Part 2.1: generate label to subset data based on outcome ----------------- #

for(i in 1:nrow(ma_data)) {
    if(ma_data$Group[i] == 1.1) {
        ma_data$label[i] <- "htn"
    } else if(ma_data$Group[i] == 1.2) {
        if(ma_data$Author[i] == "Li Y" | ma_data$Author[i] == "Zhang L L") {
            if(ma_data$Outcome[i] == "coronary artery disease") {
                ma_data$label[i] <- "cad + mi"
            } else {
                ma_data$label[i] <- "mi"
            }
        } else if(ma_data$Author[i] == "Gu J" | ma_data$Author[i] == "Wang Q") {
            if(ma_data$Outcome[i] == "CHD & MI") {
                ma_data$label[i] <- "cad + mi"
            } else if(ma_data$Outcome[i] == "coronary artery disease" | ma_data$Outcome[i] == "coronary heart disease") {
                ma_data$label[i] <- "cad"
            } else {
              ma_data$label[i] <- "mi"
            }
        }
    } else if(ma_data$Group[i] == 1.3) {
        ma_data$label[i] <- "stroke"
    } else if(ma_data$Group[i] == 2.0) {
        if(ma_data$Outcome[i] == "cancer") {
            ma_data$label[i] <- "cancer"
        } else if(ma_data$Outcome[i] == "esophageal cancer") {
            ma_data$label[i] <- "esophageal"
        } else if(ma_data$Outcome[i] == "gastric cancer") {
            ma_data$label[i] <- "gastric"
        } else if(ma_data$Outcome[i] == "hepatocellular carcinoma") {
            ma_data$label[i] <- "hepatocellular"
        } else if(ma_data$Outcome[i] == "pancreatic cancer") {
            ma_data$label[i] <- "pancreatic"
        } else if(ma_data$Outcome[i] == "colorectal cancer") {
            ma_data$label[i] <- "colorectal"
        } else if(ma_data$Outcome[i] == "head and neck cancer") {
            ma_data$label[i] <- "head"
        } else if(ma_data$Outcome[i] == "breast cancer") {
            ma_data$label[i] <- "breast"
        } else {
            ma_data$label[i] <- "exclude for synthesis"
        }
    } else if(ma_data$Group[i] == 2.1) {
        ma_data$label[i] <- "esophageal"
    } else if(ma_data$Group[i] == 2.2) {
        ma_data$label[i] <- "gastric"
    } else if(ma_data$Group[i] == 2.3) {
        ma_data$label[i] <- "hepatocellular"
    } else if(ma_data$Group[i] == 2.4) {
        ma_data$label[i] <- "pancreatic"
    } else if(ma_data$Group[i] == 2.5) {
        ma_data$label[i] <- "colorectal"
    } else if(ma_data$Group[i] == 2.6) {
        ma_data$label[i] <- "head"
    } else if(ma_data$Group[i] == 2.7) {
        ma_data$label[i] <- "breast"
    } else if(ma_data$Group[i] == 3.0) {
        if(ma_data$Outcome[i] == "alcohol-induced diseases") {
            ma_data$label[i] <- "alcoholic"
        } else {
            ma_data$label[i] <- "cirrhosis"
        }
    } else if(ma_data$Group[i] == 3.1) {
        ma_data$label[i] <- "cirrhosis"
    } else if(ma_data$Group[i] == 3.2) {
        ma_data$label[i] <- "pancreatitis"
    } else if(ma_data$Group[i] == 4.1) {
        ma_data$label[i] <- "AD"
    } else if(ma_data$Group[i] == 4.2) {
        ma_data$label[i] <- "PD"
    } else if(ma_data$Group[i] == 5.1) {
        ma_data$label[i] <- "DM"
    } else if(ma_data$Group[i] == 5.2) {
        ma_data$label[i] <- "DR"
    }
}

ma_data <- filter(ma_data, label != "exclude for synthesis")

for(m in model) {
  print(m)
  ma_data_model[[m]] <- filter(ma_data, model == m) # ma_data_model[[m]] = ma_data_"model"
}

outcome <- c("htn", "cad + mi", "cad", "mi", "stroke", "cancer", "esophageal", "gastric", "hepatocellular", "pancreatic",
           "colorectal", "head", "breast", "alcoholic", "cirrhosis", "pancreatitis", "AD", "PD", "DM", "DR")

# ----------------------- Part 2.2: fit random-effect or fix-effect model --------------------- #

ma.outcome.model <- list()
ma.outcome.model.result <- data.frame(outcome = c(NA), model = c(NA), k = c(NA), ES = c(NA),
                                      z_val = c(NA), p_val = c(NA), I2 = c(NA), tau2 = c(NA), Q = c(NA))

for(o in outcome) {
  for(m in model) {
    if(nrow(subset(ma_data_model[[m]], label == o)) >1) {
      group <- paste(m, o, sep = ",")
      print(group)
      ma <- rma(data = subset(ma_data_model[[m]], label == o),
                yi = lnOR,
                sei = SE,
                method = "DL",
                slab = paste(subset(ma_data_model[[m]], label == o)$Author, subset(ma_data_model[[m]], label == o)$Year, sep = ", "),
                digits = 2,
                level = 95)
      if(ma$I2 >= 50) {
        ma.outcome.model[[group]] <- ma       # if I2 >= 50, fit random-effect model (DL method)
      } else {                                # if I2 < 50, fit fix-effect model (DL method)
        ma.outcome.model[[group]] <- rma(data = subset(ma_data_model[[m]], label == o),
                                         yi = lnOR,
                                         sei = SE,
                                         method = "FE",
                                         slab = paste(subset(ma_data_model[[m]], label == o)$Author, subset(ma_data_model[[m]], label == o)$Year, sep = ", "),
                                         digits = 2,
                                         level = 95)
      }
      result <- c(o, m, round(ma.outcome.model[[group]]$k, 2), 
                  paste(round(ma.outcome.model[[group]]$beta, 2), " (", 
                  round(ma.outcome.model[[group]]$ci.lb, 2), ", ", round(ma.outcome.model[[group]]$ci.ub, 2), ")", sep = ""), 
                  round(ma.outcome.model[[group]]$zval, 2), round(ma.outcome.model[[group]]$pval, 2),
                  round(ma.outcome.model[[group]]$I2, 2), round(ma.outcome.model[[group]]$tau2, 2), 
                  round(ma.outcome.model[[group]]$QE, 2))
      ma.outcome.model.result <- rbind(ma.outcome.model.result, result)
    }
  }
}

ma.outcome.model.result <- ma.outcome.model.result[2:nrow(ma.outcome.model.result),]

write.csv(ma.outcome.model.result, file = "F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/meta-meta ES and CI.csv")

# ------------------------------------- Part 2.3: forestplot ---------------------------------- #
# For details of exploring the format of the forestplot, please refer to
# "ALDH2 overview code v2.R" and "ALDH2_meta-meta-analysis_code_v1.R" by Zhengting He

# -------------------------------------- fit forestplot --------------------------------------- #

model_outcome <- c("allelic,htn",
         "dominant,htn",
         "heterozygous,htn",
         "homozygous,htn",
         "recessive,htn",
         "allelic,cad + mi",
         "dominant,cad + mi",
         "heterozygous,cad + mi",
         "homozygous,cad + mi",
         "recessive,cad + mi",
         "dominant,cad",
         "dominant,mi",
         "allelic,stroke",
         "dominant,stroke",
         "heterozygous,stroke",
         "homozygous,stroke",
         "recessive,stroke",
         "allelic,cancer",
         "dominant,cancer",
         "heterozygous,cancer",
         "homozygous,cancer",
         "recessive,cancer",
         "allelic,esophageal",
         "dominant,esophageal",
         "heterozygous,esophageal",
         "homozygous,esophageal",
         "recessive,esophageal",
         "allelic,gastric",
         "dominant,gastric",
         "heterozygous,gastric",
         "homozygous,gastric",
         "recessive,gastric",
         "allelic,hepatocellular",
         "dominant,hepatocellular",
         "heterozygous,hepatocellular",
         "homozygous,hepatocellular",
         "recessive,hepatocellular",
         "allelic,pancreatic",
         "dominant,pancreatic",
         "heterozygous,pancreatic",
         "homozygous,pancreatic",
         "recessive,pancreatic",
         "allelic,colorectal",
         "dominant,colorectal",
         "heterozygous,colorectal",
         "homozygous,colorectal",
         "recessive,colorectal",
         "allelic,head",
         "dominant,head",
         "heterozygous,head",
         "homozygous,head",
         "recessive,head",
         "allelic,breast",
         "dominant,breast",
         "heterozygous,breast",
         "homozygous,breast",
         "recessive,breast",
         "allelic,cirrhosis",
         "dominant,cirrhosis",
         "dominant,AD")

for(m_o in model_outcome) {

  m_o_split <- strsplit(model_outcome, ",")
  m <- m_o_split[[1]][1]
  o <- m_o_split[[1]][2]

  data <- subset(ma_data_model[[m]], label == o)
  min_CI_low <- min(data$lnOR - 1.96 * data$SE)
  max_CI_upp <- max(data$lnOR + 1.96 * data$SE)
  min_CI_low <- min(min_CI_low, ma.outcome.model[[m_o]]$ci.lb)
  max_CI_upp <- max(max_CI_upp, ma.outcome.model[[m_o]]$ci.ub)

  left_new <- min_CI_low - 0.6
  right_new <- max_CI_upp + 0.6

  author_new <- (-0.1704 + 0.25) / (0.4 + 0.25) * (right_new - left_new) + left_new
  or_new <- right_new - (0.4 - 0.3025) / (0.4 + 0.25) * (right_new - left_new)
  legend_new <- (-0.2498 + 0.25) / (0.4 + 0.25) * (right_new - left_new) + left_new
  
  if(ma.outcome.model[[m_o]]$method == "FE") {
    method_legend <- "FE Model for All Studies"
  } else if(ma.outcome.model[[m_o]]$method == "DL") {
    method_legend <- "RE Model for All Studies"
  }

  pdf(paste("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/forestplot by outcome_new/", m_o, ".pdf", sep = ""), width = 7.23, height = 4.39)
  par(mar=c(7, 4, 7, 2)) # increase margins to condense the plot
  print(
  # set up forest plot, the 'rows' argument is used to specify in which rows the outcomes will be plotted
  forest(x = ma.outcome.model[[m_o]], showweights = FALSE, addcred = FALSE, annotate = TRUE, addfit = TRUE,
         xlim = c(left_new, right_new), #at= c(-0.2, 0, 0.2),
         xlab = "log Odds Ratio", mlab = "", psize = 0.8,
         ylim = c(-2, ma.outcome.model[[m_o]]$k + 3)))
  # add headings
  par(font=2) # switch to bold font
  text(c(author_new), ma.outcome.model[[m_o]]$k + 1.9, c("Author(s) and Year"), cex=1)
  text(c(or_new), ma.outcome.model[[m_o]]$k + 1.9, c("log Odds Ratio[95% CI]"), cex=1)
  text(legend_new, - 1, pos= 4, cex= 0.75, method_legend)
  text(legend_new, - 1.8, pos= 4, cex= 0.75, bquote(paste( "(Q = ", .(formatC(ma.outcome.model[[m_o]]$QE, digits= 2, format= "f")),
                                                        ", df = ", .(ma.outcome.model[[m_o]]$k - ma.outcome.model[[m_o]]$p), ", p = ", .(formatC(ma.outcome.model[[m_o]]$QEp, digits= 2, format= "f")),
                                                        "; ", I^ 2, " = ", .(formatC(ma.outcome.model[[m_o]]$I2, digits= 1, format= "f")), "%)")))
  dev.off()
}

# - Following: adjust left_new and right_new for those plots which components are overlapped -- #

# ------------------------ adjust those not good in format: 1st round ------------------------- #

model_outcome <- c("dominant,htn",
                   "heterozygous,cad + mi",
                   "homozygous,cad + mi",
                   "recessive,cad + mi",
                   "dominant,mi",
                   "allelic,esophageal",
                   "dominant,esophageal",
                   "heterozygous,esophageal",
                   "homozygous,esophageal",
                   "recessive,esophageal",
                   "heterozygous,gastric",
                   "homozygous,gastric",
                   "recessive,gastric",
                   "homozygous,hepatocellular",
                   "recessive,hepatocellular",
                   "dominant,pancreatic",
                   "heterozygous,pancreatic",
                   "homozygous,pancreatic",
                   "recessive,pancreatic",
                   "homozygous,colorectal",
                   "recessive,colorectal",
                   "dominant,head",
                   "heterozygous,head",
                   "homozygous,head",
                   "recessive,head",
                   "homozygous,breast",
                   "recessive,breast",
                   "allelic,cirrhosis",
                   "dominant,cirrhosis",
                   "dominant,AD")

for(m_o in model_outcome) {

  m_o_split <- strsplit(model_outcome, ",")
  m <- m_o_split[[1]][1]
  o <- m_o_split[[1]][2]

  data <- subset(ma_data_model[[m]], label == o)
  min_CI_low <- min(data$lnOR - 1.96 * data$SE)
  max_CI_upp <- max(data$lnOR + 1.96 * data$SE)
  min_CI_low <- min(min_CI_low, ma.outcome.model[[m_o]]$ci.lb)
  max_CI_upp <- max(max_CI_upp, ma.outcome.model[[m_o]]$ci.ub)

  left_new <- min_CI_low - 0.6
  right_new <- max_CI_upp + 1

  author_new <- (-0.1704 + 0.25) / (0.4 + 0.25) * (right_new - left_new) + left_new
  or_new <- right_new - (0.4 - 0.3025) / (0.4 + 0.25) * (right_new - left_new)
  legend_new <- (-0.2498 + 0.25) / (0.4 + 0.25) * (right_new - left_new) + left_new

  if(ma.outcome.model[[m_o]]$method == "FE") {
    method_legend <- "FE Model for All Studies"
  } else if(ma.outcome.model[[m_o]]$method == "DL") {
    method_legend <- "RE Model for All Studies"
  }

  pdf(paste("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/forestplot by outcome_new/", m_o, ".pdf", sep = ""), width = 7.23, height = 4.39)
  par(mar=c(7, 4, 7, 2)) # increase margins to condense the plot
  print(
    # set up forest plot, the 'rows' argument is used to specify in which rows the outcomes will be plotted
    forest(x = ma.outcome.model[[m_o]], showweights = FALSE, addcred = FALSE, annotate = TRUE, addfit = TRUE,
           xlim = c(left_new, right_new), #at= c(-0.2, 0, 0.2),
           xlab = "log Odds Ratio", mlab = "", psize = 0.8,
           ylim = c(-2, ma.outcome.model[[m_o]]$k + 3)))
  # add headings
  par(font=2) # switch to bold font
  text(c(author_new), ma.outcome.model[[m_o]]$k + 1.9, c("Author(s) and Year"), cex=1)
  text(c(or_new), ma.outcome.model[[m_o]]$k + 1.9, c("log Odds Ratio[95% CI]"), cex=1)
  text(legend_new, - 1, pos= 4, cex= 0.75, method_legend)
  text(legend_new, - 1.8, pos= 4, cex= 0.75, bquote(paste( "(Q = ", .(formatC(ma.outcome.model[[m_o]]$QE, digits= 2, format= "f")),
                                                           ", df = ", .(ma.outcome.model[[m_o]]$k - ma.outcome.model[[m_o]]$p), ", p = ", .(formatC(ma.outcome.model[[m_o]]$QEp, digits= 2, format= "f")),
                                                           "; ", I^ 2, " = ", .(formatC(ma.outcome.model[[m_o]]$I2, digits= 1, format= "f")), "%)")))
  dev.off()
}

# ------------------------ adjust those not good in format: 2nd round ------------------------- #

model_outcome <- c("homozygous,cad + mi",
                   "dominant,mi",
                   "dominant,esophageal",
                   "heterozygous,esophageal",
                   "homozygous,esophageal",
                   "recessive,esophageal",
                   "heterozygous,gastric",
                   "recessive,gastric",
                   "homozygous,hepatocellular",
                   "recessive,hepatocellular",
                   "homozygous,pancreatic",
                   "homozygous,colorectal",
                   "recessive,colorectal",
                   "heterozygous,head",
                   "homozygous,breast",
                   "allelic,cirrhosis",
                   "dominant,cirrhosis",
                   "dominant,AD")

for(m_o in model_outcome) {

  m_o_split <- strsplit(model_outcome, ",")
  m <- m_o_split[[1]][1]
  o <- m_o_split[[1]][2]

  data <- subset(ma_data_model[[m]], label == o)
  min_CI_low <- min(data$lnOR - 1.96 * data$SE)
  max_CI_upp <- max(data$lnOR + 1.96 * data$SE)
  min_CI_low <- min(min_CI_low, ma.outcome.model[[m_o]]$ci.lb)
  max_CI_upp <- max(max_CI_upp, ma.outcome.model[[m_o]]$ci.ub)

  left_new <- min_CI_low - 1.3
  right_new <- max_CI_upp + 1.1

  author_new <- (-0.1704 + 0.25) / (0.4 + 0.25) * (right_new - left_new) + left_new
  or_new <- right_new - (0.4 - 0.3025) / (0.4 + 0.25) * (right_new - left_new)
  legend_new <- (-0.2498 + 0.25) / (0.4 + 0.25) * (right_new - left_new) + left_new

  if(ma.outcome.model[[m_o]]$method == "FE") {
    method_legend <- "FE Model for All Studies"
  } else if(ma.outcome.model[[m_o]]$method == "DL") {
    method_legend <- "RE Model for All Studies"
  }

  pdf(paste("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/forestplot by outcome_new/", m_o, ".pdf", sep = ""), width = 7.23, height = 4.39)
  par(mar=c(7, 4, 7, 2)) # increase margins to condense the plot
  print(
    # set up forest plot, the 'rows' argument is used to specify in which rows the outcomes will be plotted
    forest(x = ma.outcome.model[[m_o]], showweights = FALSE, addcred = FALSE, annotate = TRUE, addfit = TRUE,
           xlim = c(left_new, right_new), #at= c(-0.2, 0, 0.2),
           xlab = "log Odds Ratio", mlab = "", psize = 0.8,
           ylim = c(-2, ma.outcome.model[[m_o]]$k + 3)))
  # add headings
  par(font=2) # switch to bold font
  text(c(author_new), ma.outcome.model[[m_o]]$k + 1.9, c("Author(s) and Year"), cex=1)
  text(c(or_new), ma.outcome.model[[m_o]]$k + 1.9, c("log Odds Ratio[95% CI]"), cex=1)
  text(legend_new, - 1, pos= 4, cex= 0.75, method_legend)
  text(legend_new, - 1.8, pos= 4, cex= 0.75, bquote(paste( "(Q = ", .(formatC(ma.outcome.model[[m_o]]$QE, digits= 2, format= "f")),
                                                           ", df = ", .(ma.outcome.model[[m_o]]$k - ma.outcome.model[[m_o]]$p), ", p = ", .(formatC(ma.outcome.model[[m_o]]$QEp, digits= 2, format= "f")),
                                                           "; ", I^ 2, " = ", .(formatC(ma.outcome.model[[m_o]]$I2, digits= 1, format= "f")), "%)")))
  dev.off()
}

# ------------------------ adjust those not good in format: 3rd round ------------------------- #

model_outcome <- c("dominant,mi",
                   "homozygous,esophageal",
                   "recessive,esophageal",
                   "heterozygous,gastric",
                   "recessive,gastric",
                   "homozygous,hepatocellular",
                   "recessive,hepatocellular",
                   "homozygous,pancreatic",
                   "homozygous,colorectal",
                   "recessive,colorectal",
                   "heterozygous,head",
                   "homozygous,breast",
                   "dominant,cirrhosis",
                   "dominant,AD")

for(m_o in model_outcome) {
  
  m_o_split <- strsplit(model_outcome, ",")
  m <- m_o_split[[1]][1]
  o <- m_o_split[[1]][2]
  
  data <- subset(ma_data_model[[m]], label == o)
  min_CI_low <- min(data$lnOR - 1.96 * data$SE)
  max_CI_upp <- max(data$lnOR + 1.96 * data$SE)
  min_CI_low <- min(min_CI_low, ma.outcome.model[[m_o]]$ci.lb)
  max_CI_upp <- max(max_CI_upp, ma.outcome.model[[m_o]]$ci.ub)
  
  left_new <- min_CI_low - 1.2
  right_new <- max_CI_upp + 0.8
  
  author_new <- (-0.1704 + 0.25) / (0.4 + 0.25) * (right_new - left_new) + left_new
  or_new <- right_new - (0.4 - 0.3025) / (0.4 + 0.25) * (right_new - left_new)
  legend_new <- (-0.2498 + 0.25) / (0.4 + 0.25) * (right_new - left_new) + left_new
  
  if(ma.outcome.model[[m_o]]$method == "FE") {
    method_legend <- "FE Model for All Studies"
  } else if(ma.outcome.model[[m_o]]$method == "DL") {
    method_legend <- "RE Model for All Studies"
  }
  
  pdf(paste("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/forestplot by outcome_new/", m_o, ".pdf", sep = ""), width = 7.23, height = 4.39)
  par(mar=c(7, 4, 7, 2)) # increase margins to condense the plot
  print(
    # set up forest plot, the 'rows' argument is used to specify in which rows the outcomes will be plotted
    forest(x = ma.outcome.model[[m_o]], showweights = FALSE, addcred = FALSE, annotate = TRUE, addfit = TRUE,
           xlim = c(left_new, right_new), #at= c(-0.2, 0, 0.2),
           xlab = "log Odds Ratio", mlab = "", psize = 0.8,
           ylim = c(-2, ma.outcome.model[[m_o]]$k + 3)))
  # add headings
  par(font=2) # switch to bold font
  text(c(author_new), ma.outcome.model[[m_o]]$k + 1.9, c("Author(s) and Year"), cex=1)
  text(c(or_new), ma.outcome.model[[m_o]]$k + 1.9, c("log Odds Ratio[95% CI]"), cex=1)
  text(legend_new, - 1, pos= 4, cex= 0.75, method_legend)
  text(legend_new, - 1.8, pos= 4, cex= 0.75, bquote(paste( "(Q = ", .(formatC(ma.outcome.model[[m_o]]$QE, digits= 2, format= "f")),
                                                           ", df = ", .(ma.outcome.model[[m_o]]$k - ma.outcome.model[[m_o]]$p), ", p = ", .(formatC(ma.outcome.model[[m_o]]$QEp, digits= 2, format= "f")),
                                                           "; ", I^ 2, " = ", .(formatC(ma.outcome.model[[m_o]]$I2, digits= 1, format= "f")), "%)")))
  dev.off()
}

# ------------------------ adjust those not good in format: 4th round ------------------------- #

model_outcome <- c("homozygous,esophageal",
                   "recessive,esophageal",
                   "heterozygous,gastric",
                   "recessive,gastric",
                   "homozygous,hepatocellular",
                   "recessive,hepatocellular",
                   "homozygous,colorectal",
                   "recessive,colorectal",
                   "homozygous,breast",
                   "dominant,cirrhosis")

for(m_o in model_outcome) {
  
  m_o_split <- strsplit(model_outcome, ",")
  m <- m_o_split[[1]][1]
  o <- m_o_split[[1]][2]
  
  data <- subset(ma_data_model[[m]], label == o)
  min_CI_low <- min(data$lnOR - 1.96 * data$SE)
  max_CI_upp <- max(data$lnOR + 1.96 * data$SE)
  min_CI_low <- min(min_CI_low, ma.outcome.model[[m_o]]$ci.lb)
  max_CI_upp <- max(max_CI_upp, ma.outcome.model[[m_o]]$ci.ub)
  
  left_new <- min_CI_low - 1.25
  right_new <- max_CI_upp + 0.8
  
  author_new <- (-0.1704 + 0.25) / (0.4 + 0.25) * (right_new - left_new) + left_new
  or_new <- right_new - (0.4 - 0.3025) / (0.4 + 0.25) * (right_new - left_new)
  legend_new <- (-0.2498 + 0.25) / (0.4 + 0.25) * (right_new - left_new) + left_new
  
  if(ma.outcome.model[[m_o]]$method == "FE") {
    method_legend <- "FE Model for All Studies"
  } else if(ma.outcome.model[[m_o]]$method == "DL") {
    method_legend <- "RE Model for All Studies"
  }
  
  pdf(paste("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/forestplot by outcome_new/", m_o, ".pdf", sep = ""), width = 7.23, height = 4.39)
  par(mar=c(7, 4, 7, 2)) # increase margins to condense the plot
  print(
    # set up forest plot, the 'rows' argument is used to specify in which rows the outcomes will be plotted
    forest(x = ma.outcome.model[[m_o]], showweights = FALSE, addcred = FALSE, annotate = TRUE, addfit = TRUE,
           xlim = c(left_new, right_new), #at= c(-0.2, 0, 0.2),
           xlab = "log Odds Ratio", mlab = "", psize = 0.8,
           ylim = c(-2, ma.outcome.model[[m_o]]$k + 3)))
  # add headings
  par(font=2) # switch to bold font
  text(c(author_new), ma.outcome.model[[m_o]]$k + 1.9, c("Author(s) and Year"), cex=1)
  text(c(or_new), ma.outcome.model[[m_o]]$k + 1.9, c("log Odds Ratio[95% CI]"), cex=1)
  text(legend_new, - 1, pos= 4, cex= 0.75, method_legend)
  text(legend_new, - 1.8, pos= 4, cex= 0.75, bquote(paste( "(Q = ", .(formatC(ma.outcome.model[[m_o]]$QE, digits= 2, format= "f")),
                                                           ", df = ", .(ma.outcome.model[[m_o]]$k - ma.outcome.model[[m_o]]$p), ", p = ", .(formatC(ma.outcome.model[[m_o]]$QEp, digits= 2, format= "f")),
                                                           "; ", I^ 2, " = ", .(formatC(ma.outcome.model[[m_o]]$I2, digits= 1, format= "f")), "%)")))
  dev.off()
}

# ------------------------ adjust those not good in format: 5th round ------------------------- #

model_outcome <- c("homozygous,esophageal",
                   "recessive,gastric")

for(m_o in model_outcome) {
  
  m_o_split <- strsplit(model_outcome, ",")
  m <- m_o_split[[1]][1]
  o <- m_o_split[[1]][2]
  
  data <- subset(ma_data_model[[m]], label == o)
  min_CI_low <- min(data$lnOR - 1.96 * data$SE)
  max_CI_upp <- max(data$lnOR + 1.96 * data$SE)
  min_CI_low <- min(min_CI_low, ma.outcome.model[[m_o]]$ci.lb)
  max_CI_upp <- max(max_CI_upp, ma.outcome.model[[m_o]]$ci.ub)
  
  left_new <- min_CI_low - 1.25
  right_new <- max_CI_upp + 1.5
  
  author_new <- (-0.1704 + 0.25) / (0.4 + 0.25) * (right_new - left_new) + left_new
  or_new <- right_new - (0.4 - 0.3025) / (0.4 + 0.25) * (right_new - left_new)
  legend_new <- (-0.2498 + 0.25) / (0.4 + 0.25) * (right_new - left_new) + left_new
  
  if(ma.outcome.model[[m_o]]$method == "FE") {
    method_legend <- "FE Model for All Studies"
  } else if(ma.outcome.model[[m_o]]$method == "DL") {
    method_legend <- "RE Model for All Studies"
  }
  
  pdf(paste("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/forestplot by outcome_new/", m_o, ".pdf", sep = ""), width = 7.23, height = 4.39)
  par(mar=c(7, 4, 7, 2)) # increase margins to condense the plot
  print(
    # set up forest plot, the 'rows' argument is used to specify in which rows the outcomes will be plotted
    forest(x = ma.outcome.model[[m_o]], showweights = FALSE, addcred = FALSE, annotate = TRUE, addfit = TRUE,
           xlim = c(left_new, right_new), #at= c(-0.2, 0, 0.2),
           xlab = "log Odds Ratio", mlab = "", psize = 0.8,
           ylim = c(-2, ma.outcome.model[[m_o]]$k + 3)))
  # add headings
  par(font=2) # switch to bold font
  text(c(author_new), ma.outcome.model[[m_o]]$k + 1.9, c("Author(s) and Year"), cex=1)
  text(c(or_new), ma.outcome.model[[m_o]]$k + 1.9, c("log Odds Ratio[95% CI]"), cex=1)
  text(legend_new, - 1, pos= 4, cex= 0.75, method_legend)
  text(legend_new, - 1.8, pos= 4, cex= 0.75, bquote(paste( "(Q = ", .(formatC(ma.outcome.model[[m_o]]$QE, digits= 2, format= "f")),
                                                           ", df = ", .(ma.outcome.model[[m_o]]$k - ma.outcome.model[[m_o]]$p), ", p = ", .(formatC(ma.outcome.model[[m_o]]$QEp, digits= 2, format= "f")),
                                                           "; ", I^ 2, " = ", .(formatC(ma.outcome.model[[m_o]]$I2, digits= 1, format= "f")), "%)")))
  dev.off()
}

# --------------------------------------------------------------------------------------------- #
#  Part 3: sensitivity analysis for studies which ES and CI are not consistent with re-calculate results  #
# --------------------------------------------------------------------------------------------- #

ma_data_sen <- read.csv("data extraction table (meta-analysis)_sensitivity data.csv")

# uniform all observations using "ALDH2*2" as exposure and "ALDH2*1" as comparison
ma_data_sen$CI1_copy <- ma_data_sen$CI1
for(i in 1:nrow(ma_data_sen)) {
  if(ma_data_sen$Exposure[i] == "ALDH2*1" & ma_data_sen$Comparison[i] == "ALDH2*2") {
    ma_data_sen$OR[i] <- 1/ma_data_sen$OR[i]
    ma_data_sen$CI1[i] <- 1/ma_data_sen$CI2[i]
    ma_data_sen$CI2[i] <- 1/ma_data_sen$CI1_copy[i]
    ma_data_sen$Exposure[i] <- "ALDH2*2"
    ma_data_sen$Comparison[i] <- "ALDH2*1"
  }
}
ma_data_sen <- subset(ma_data_sen, select = - CI1_copy)

write.csv(ma_data_sen, file = "F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/data extraction table (meta-analysis)_sensitivity data_uniform exposure and comparison.csv")

# subset ma_data to measurement of effect size as "Odds Ratio"
# Only studies which measurement of effect size is "Odds Ratio" are eligable for synthesizing
ma_data_sen <- filter(ma_data_sen, ES_measurement == "OR")

# convert effect size and confidential intervial to log value, generate SE_lnOR
# formula is generated according to https://stats.stackexchange.com/questions/156597/how-to-calculate-se-of-an-odds-ratio
# CI = exp(lnOR +- 1.96 * sqrt(1/a + 1/b + 1/c + 1/d))
# ln(CI_lower) = ln(OR) - 1.96 * SE_lnOR
# ln(CI_upper) = ln(OR) + 1.96 * SE_lnOR
# SE_lnOR = 1/(2*1.96) * (ln(CI_upper) - ln(CI_lower))
ma_data_sen$lnOR <- log(ma_data_sen$OR)
ma_data_sen$SE <- 1/3.92 * (log(ma_data_sen$CI2) - log(ma_data_sen$CI1))
ma_data_sen$lnOR <- round(ma_data_sen$lnOR, 2) # round lnOR to 2-digit floating number
ma_data_sen$SE <- round(ma_data_sen$SE, 2) # round SE to 2-digit floating number

# separate ma_data by genetic model, generate ma_data_"model"
model <- c("allelic", "dominant", "heterozygous", "homozygous", "recessive")
ma_data_sen_model <- list()
for(m in model) {
  print(m)
  ma_data_sen_model[[m]] <- filter(ma_data_sen, model == m) # ma_data_model[[m]] = ma_data_"model"
}

ma_data_sen$Group <- as.numeric(ma_data_sen$Group)

# ------------------ Part 3.1: generate label to subset data based on outcome ----------------- #

for(i in 1:nrow(ma_data_sen)) {
  if(ma_data_sen$Group[i] == 1.1) {
    ma_data_sen$label[i] <- "htn"
  } else if(ma_data_sen$Group[i] == 1.2) {
    if(ma_data_sen$Author[i] == "Li Y" | ma_data_sen$Author[i] == "Zhang L L") {
      if(ma_data_sen$Outcome[i] == "coronary artery disease") {
        ma_data_sen$label[i] <- "cad + mi"
      } else {
        ma_data_sen$label[i] <- "mi"
      }
    } else if(ma_data_sen$Author[i] == "Gu J" | ma_data_sen$Author[i] == "Wang Q") {
      if(ma_data_sen$Outcome[i] == "CHD & MI") {
        ma_data_sen$label[i] <- "cad + mi"
      } else if(ma_data_sen$Outcome[i] == "coronary artery disease" | ma_data_sen$Outcome[i] == "coronary heart disease") {
        ma_data_sen$label[i] <- "cad"
      } else {
        ma_data_sen$label[i] <- "mi"
      }
    }
  } else if(ma_data_sen$Group[i] == 1.3) {
    ma_data_sen$label[i] <- "stroke"
  } else if(ma_data_sen$Group[i] == 2.0) {
    if(ma_data_sen$Outcome[i] == "cancer") {
      ma_data_sen$label[i] <- "cancer"
    } else if(ma_data_sen$Outcome[i] == "esophageal cancer") {
      ma_data_sen$label[i] <- "esophageal"
    } else if(ma_data_sen$Outcome[i] == "gastric cancer") {
      ma_data_sen$label[i] <- "gastric"
    } else if(ma_data_sen$Outcome[i] == "hepatocellular carcinoma") {
      ma_data_sen$label[i] <- "hepatocellular"
    } else if(ma_data_sen$Outcome[i] == "pancreatic cancer") {
      ma_data_sen$label[i] <- "pancreatic"
    } else if(ma_data_sen$Outcome[i] == "colorectal cancer") {
      ma_data_sen$label[i] <- "colorectal"
    } else if(ma_data_sen$Outcome[i] == "head and neck cancer") {
      ma_data_sen$label[i] <- "head"
    } else if(ma_data_sen$Outcome[i] == "breast cancer") {
      ma_data_sen$label[i] <- "breast"
    } else {
      ma_data_sen$label[i] <- "exclude for synthesis"
    }
  } else if(ma_data_sen$Group[i] == 2.1) {
    ma_data_sen$label[i] <- "esophageal"
  } else if(ma_data_sen$Group[i] == 2.2) {
    ma_data_sen$label[i] <- "gastric"
  } else if(ma_data_sen$Group[i] == 2.3) {
    ma_data_sen$label[i] <- "hepatocellular"
  } else if(ma_data_sen$Group[i] == 2.4) {
    ma_data_sen$label[i] <- "pancreatic"
  } else if(ma_data_sen$Group[i] == 2.5) {
    ma_data_sen$label[i] <- "colorectal"
  } else if(ma_data_sen$Group[i] == 2.6) {
    ma_data_sen$label[i] <- "head"
  } else if(ma_data_sen$Group[i] == 2.7) {
    ma_data_sen$label[i] <- "breast"
  } else if(ma_data_sen$Group[i] == 3.0) {
    if(ma_data_sen$Outcome[i] == "alcohol-induced diseases") {
      ma_data_sen$label[i] <- "alcoholic"
    } else {
      ma_data_sen$label[i] <- "cirrhosis"
    }
  } else if(ma_data_sen$Group[i] == 3.1) {
    ma_data_sen$label[i] <- "cirrhosis"
  } else if(ma_data_sen$Group[i] == 3.2) {
    ma_data_sen$label[i] <- "pancreatitis"
  } else if(ma_data_sen$Group[i] == 4.1) {
    ma_data_sen$label[i] <- "AD"
  } else if(ma_data_sen$Group[i] == 4.2) {
    ma_data_sen$label[i] <- "PD"
  } else if(ma_data_sen$Group[i] == 5.1) {
    ma_data_sen$label[i] <- "DM"
  } else if(ma_data_sen$Group[i] == 5.2) {
    ma_data_sen$label[i] <- "DR"
  }
}

ma_data_sen <- filter(ma_data_sen, label != "exclude for synthesis")

for(m in model) {
  print(m)
  ma_data_sen_model[[m]] <- filter(ma_data_sen, model == m) # ma_data_model[[m]] = ma_data_"model"
}

outcome <- c("htn", "cad + mi", "cad", "mi", "stroke", "cancer", "esophageal", "gastric", "hepatocellular", "pancreatic",
             "colorectal", "head", "breast", "alcoholic", "cirrhosis", "pancreatitis", "AD", "PD", "DM", "DR")

# ----------------------- Part 3.2: fit random-effect or fix-effect model --------------------- #

ma.outcome.model.sen <- list()
ma.outcome.model.result.sen <- data.frame(outcome = c(NA), model = c(NA), k = c(NA), ES = c(NA),
                                      z_val = c(NA), p_val = c(NA), I2 = c(NA), tau2 = c(NA), Q = c(NA))

for(o in outcome) {
  for(m in model) {
    if(isTRUE(table(subset(ma_data_sen_model[[m]], label == o)$Notes) == nrow(subset(ma_data_sen_model[[m]], label == o))) == FALSE) {
      if(nrow(subset(ma_data_sen_model[[m]], label == o)) >1) {
        group <- paste(m, o, sep = ",")
        print(group)
        ma <- rma(data = subset(ma_data_sen_model[[m]], label == o),
                  yi = lnOR,
                  sei = SE,
                  method = "DL",
                  slab = paste(subset(ma_data_sen_model[[m]], label == o)$Author, subset(ma_data_sen_model[[m]], label == o)$Year, sep = ", "),
                  digits = 2,
                  level = 95)
        if(ma$I2 >= 50) {
          ma.outcome.model.sen[[group]] <- ma       # if I2 >= 50, fit random-effect model (DL method)
        } else {                                # if I2 < 50, fit fix-effect model (DL method)
          ma.outcome.model.sen[[group]] <- rma(data = subset(ma_data_sen_model[[m]], label == o),
                                           yi = lnOR,
                                           sei = SE,
                                           method = "FE",
                                           slab = paste(subset(ma_data_sen_model[[m]], label == o)$Author, subset(ma_data_sen_model[[m]], label == o)$Year, sep = ", "),
                                           digits = 2,
                                           level = 95)
        }
        result <- c(o, m, round(ma.outcome.model.sen[[group]]$k, 2), 
                    paste(round(ma.outcome.model.sen[[group]]$beta, 2), " (", 
                          round(ma.outcome.model.sen[[group]]$ci.lb, 2), ", ", round(ma.outcome.model.sen[[group]]$ci.ub, 2), ")", sep = ""), 
                    round(ma.outcome.model.sen[[group]]$zval, 2), round(ma.outcome.model.sen[[group]]$pval, 2),
                    round(ma.outcome.model.sen[[group]]$I2, 2), round(ma.outcome.model.sen[[group]]$tau2, 2), 
                    round(ma.outcome.model.sen[[group]]$QE, 2))
        ma.outcome.model.result.sen <- rbind(ma.outcome.model.result.sen, result)
      }
    }
  }
}

ma.outcome.model.result.sen <- ma.outcome.model.result.sen[2:nrow(ma.outcome.model.result),]
ma.outcome.model.result.sen <- filter(ma.outcome.model.result.sen, is.na(ES) == FALSE)

write.csv(ma.outcome.model.result.sen, file = "F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/meta-meta ES and CI_sensitivity analysis.csv")

# ------------- Part 3.3: compare the sensitivity results with the original results ----------- #

for(m_o in c("allelic,cad + mi", "dominant,cad + mi", "heterozygous,cad + mi", "homozygous,cad + mi", "recessive,cad + mi",
             "allelic,cancer", "heterozygous,esophageal", "homozygous,esophageal", "recessive,gastric",
             "allelic,hepatocellular", "dominant,hepatocellular", "heterozygous,hepatocellular", "homozygous,hepatocellular",
             "recessive,hepatocellular", "dominant,pancreatic", "heterozygous,pancreatic", "allelic,colorectal", "dominant,colorectal")) {
  if(ma.outcome.model[[m_o]]$ci.ub < 0) {
    if(ma.outcome.model.sen[[m_o]]$ci.ub < 0) {
      print(paste(m_o, "not changed", sep = ", "))
    } else {
      print(paste(m_o, "changed", sep = ", "))
    }
  } else if(ma.outcome.model[[m_o]]$ci.ub >= 0 & ma.outcome.model[[m_o]]$ci.lb <= 0) {
    if(ma.outcome.model.sen[[m_o]]$ci.ub < 0) {
      print(paste(m_o, "changed", sep = ", "))
    } else {
      if(ma.outcome.model.sen[[m_o]]$ci.lb <= 0) {
        print(paste(m_o, "not changed", sep = ", "))
      } else {
        print(paste(m_o, "changed", sep = ", "))
      }
    }
  } else if(ma.outcome.model[[m_o]]$ci.lb > 0) {
    if(ma.outcome.model.sen[[m_o]]$ci.ub < 0) {
      print(paste(m_o, "changed", sep = ", "))
    } else {
      if(ma.outcome.model.sen[[m_o]]$ci.lb <= 0) {
        print(paste(m_o, "changed", sep = ", "))
      } else if(ma.outcome.model.sen[[m_o]]$ci.lb > 0) {
        print(paste(m_o, "not changed", sep = ", "))
      }
    }
  }
} # "heterozygous,pancreatic, changed"

# ----- Part 3.4: forestplot of the MA^2 which results changed during sensitivity analysis ---- #

# heterozygous,pancreatic
pdf("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/meta-meta-analysis/output/forestplot by outcome_new/heterozygous,pancreatic_sensitivity.pdf", width = 7.23, height = 4.39)
par(mar=c(7, 4, 7, 2)) # increase margins to condense the plot
# set up forest plot, the 'rows' argument is used to specify in which rows the outcomes will be plotted
forest(x = ma.outcome.model.sen[["heterozygous,pancreatic"]], showweights = FALSE, addcred = FALSE, annotate = TRUE, addfit = TRUE,
       xlim = c(- 0.5, 1), at= c(-0.1, 0, 0.2, 0.4, 0.6),
       xlab = "log Odds Ratio", mlab = "", psize = 0.8,
       ylim = c(-2, 6))
# add headings
par(font=2) # switch to bold font
text(c(-0.3163), 4.9, c("Author(s) and Year"), cex=1)
text(c(0.775), 4.9, c("log Odds Ratio[95% CI]"), cex=1)
text(-0.4995, - 1, pos= 4, cex= 0.75, bquote("FE Model for All Studies"))
text(-0.4995, - 1.8, pos= 4, cex= 0.75, bquote(paste( "(Q = ", .(formatC(ma.outcome.model.sen[["heterozygous,pancreatic"]]$QE, digits= 2, format= "f")), 
                                                      ", df = ", .(ma.outcome.model.sen[["heterozygous,pancreatic"]]$k - ma.outcome.model.sen[["heterozygous,pancreatic"]]$p), ", p = ", .(formatC(ma.outcome.model.sen[["heterozygous,pancreatic"]]$QEp, digits= 2, format= "f")), 
                                                      "; ", I^ 2, " = ", .(formatC(ma.outcome.model.sen[["heterozygous,pancreatic"]]$I2, digits= 1, format= "f")), "%)")))
dev.off()

# --------------------------------------------------------------------------------------------- #
# ----------------------------- part 4: citation matrix, CA, CCA ------------------------------ #
# ------------------------------ aim: present degree of overlap ------------------------------- #
# --------------------------------------------------------------------------------------------- #

# report overlap
# overlap is calculated using CA and CCA according to D. Pieper et al., 2014
# CA (Covered area) = N / (r * c)
# CCA (Corrected CA) = (N - r) / (r * c - r)
# interpret CCA: 0-5, slight overlap; 6-10, moderate overlap; 11-15, high overlap; >15, very high overlap
# setwd("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data synthesis/citation matrix/data")
# filenames = list.files(pattern="*.csv") # a vector of all file names of citation matrix csv document named as "outcome"
# overlap <- data.frame(row.names = c("CA", "CCA", "degree")) # empty data frame for CA, CCA
# for(i in 1:length(filenames)) { # for each outcome, import citation matrix, calculate CA, CCA, interpret overlap degree
#   dat <- read.csv(filenames[i])
#   apply(dat, 2, as.numeric)
#   CA <- sum(rowSums(dat, na.rm = TRUE))/(nrow(dat)*ncol(dat))
#   CCA <- (sum(rowSums(dat, na.rm = TRUE))- nrow(dat))/(nrow(dat)*ncol(dat)-nrow(dat))
#   if(CCA <= 0.05) {
#     degree <- "slight"
#   } else if(CCA <= 0.1) {
#     degree <- "moderate"
#   } else if(CCA <= 0.15) {
#     degree <- "high"
#   } else if(CCA > 0.15) {
#     degree <- "very high"
#   }
#   overlap[,i] <- c(CA, CCA, degree)
#   # CA <- c(CA, sum(rowSums(read.csv((i)), na.rm = TRUE))/(nrow(read.csv((i)))*ncol(read.csv((i)))))
# }
# colnames(overlap) <- make.names(gsub("*.csv$", "",list.files(pattern="*.csv"))) # change colnames as "outcome"
# write.csv(overlap, file = "F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/primary analysis/final decision/result/part 1.2 overlap/overlap.csv")
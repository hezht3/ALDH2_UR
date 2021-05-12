<<<<<<< HEAD
require(dplyr)

############################### method 1: use data extraction table of primary studies ###############################

setwd("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data preparation/original data")

data <- read.csv("data extraction table (primary study)_combine.csv")

data$MA <- paste(data$MA..First.author, ", ", data$MA..Year)
data$PR <- paste(data$PR..First.author, ", ", data$PR..Year)
data$MA_PR <- paste(data$MA, ", ", data$PR)

# exclude those studies that are not included in synthesis
data$note <- NA

for(i in 1:nrow(data)) {
  if(data$MA..First.author[i] == "Zuo W" | data$MA..First.author[i] == "Cai Q") {
    if(data$Outcome[i] != "esophageal cancer" & data$Outcome[i] != "gastric cancer" & data$Outcome[i] != "hepatocellular carcinoma" &
       data$Outcome[i] != "pancreatic cancer" & data$Outcome[i] != "colorectal cancer" & data$Outcome[i] != "head and neck cancer") {
      data$note[i] <- "exclude for synthesis"
    }
  }
  
  if(data$Measurement.of.effect.size[i] == "RR") {
    data$note[i] <- "exclude for synthesis"
  }
}

data <- filter(data, is.na(note) == TRUE)

# generate label for subset
data$label <- data$Outcome

for(i in 1:nrow(data)) {
  if(isTRUE(data$Group[i] == 1.2)) {
    if(data$MA..First.author[i] == "Li Y" | data$MA..First.author[i] == "Zhang L L") {
      if(data$Outcome[i] == "coronary artery disease") {
        data$label[i] <- "cad + mi"
      } else {
        data$label[i] <- "mi"
      }
    } else if(data$MA..First.author[i] == "Gu J" | data$MA..First.author[i] == "Wang Q") {
      if(data$Outcome[i] == "CHD & MI") {
        data$label[i] <- "cad + mi"
      } else if(data$Outcome[i] == "coronary artery disease" | data$Outcome[i] == "coronary heart disease") {
        data$label[i] <- "cad"
      } else {
        data$label[i] <- "mi"
      }
    }
  }
}

data$label <- recode(data$label, "essential hypertension" = "htn", "hypertension" = "htn",
                                 "esophageal cancer" = "esophageal", "gastric cancer" = "gastric",
                                 "hepatocellular carcinoma" = "hepatocellular", "pancreatic cancer" = "pancreatic",
                                 "colorectal cancer" = "colorectal", "head and neck cancer" = "head",
                                 "breast cancer" = "breast", "late-onset Alzheimer Disease" = "AD",
                                 "Alzheimer Disease" = "AD")

label <- as.data.frame(table(data$label))
label <- as.vector(label$Var1)

# subset and calculate CA, CCA
# r <- number of primary studies
# c <- number of meta-analyses
# N <- number of ticked boxes in the citation matrix
# CA <- N / (r * c)
# CCA <- (N - r) / (r * c - r)
# interpret CCA: 0-5, slight overlap; 6-10, moderate overlap; 11-15, high overlap; >15, very high overlap

data_outcome <- list()
overlap <- data.frame(row.names = c("label", "CA", "CCA", "degree"))

for(lab in label) {
  data_outcome[[lab]] <- filter(data, label == lab)
  
  c <- nrow(table(data_outcome[[lab]]$MA))
  r <- nrow(table(data_outcome[[lab]]$PR))
  N <- nrow(table(data_outcome[[lab]]$MA_PR))
  
  CA <- round(N / (r * c), 2)
  CCA <- round((N - r) / (r * c - r), 2)
  
  if(CCA <= 0.05) {
    degree <- "slight"
  } else if(CCA <= 0.1) {
    degree <- "moderate"
  } else if(CCA <= 0.15) {
    degree <- "high"
  } else if(CCA > 0.15) {
    degree <- "very high"
  }
  
  overlap <- bind_cols(overlap, c(lab, CA, CCA, degree))
}

colnames(overlap) <- overlap[1, ]

########################################### method 2: use citation matrixes ##########################################

setwd("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/primary analysis/citation matrix/citation matrix 20210221")
filenames = list.files(pattern="*.csv")
overlap2 <- data.frame(row.names = c("label", "CA", "CCA", "degree"))
data_outcome2 <- list()

for(file in filenames) {
  label2 <- strsplit(file, ".csv")[[1]]
  dat <- read.csv(file)
  data_outcome2[[label2]] <- dat
  
  dat <- dat[, 2:ncol(dat)]
  for(i in 1:ncol(dat)) {
    dat[, i] <- recode(dat[, i], "x" = "1")
    dat[, i] <- as.numeric(dat[, i])
  }
  
  CA <- sum(rowSums(dat, na.rm = TRUE))/(nrow(dat)*ncol(dat))
  CCA <- (sum(rowSums(dat, na.rm = TRUE))- nrow(dat))/(nrow(dat)*ncol(dat)-nrow(dat))
  CA <- round(CA, 2)
  CCA <- round(CCA, 2)
  if(CCA <= 0.05) {
    degree <- "slight"
  } else if(CCA <= 0.1) {
    degree <- "moderate"
  } else if(CCA <= 0.15) {
    degree <- "high"
  } else if(CCA > 0.15) {
    degree <- "very high"
  }
  
  overlap2 <- bind_cols(overlap2, c(label2, CA, CCA, degree))
}

=======
require(dplyr)

############################### method 1: use data extraction table of primary studies ###############################

setwd("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/data preparation/original data")

data <- read.csv("data extraction table (primary study)_combine.csv")

data$MA <- paste(data$MA..First.author, ", ", data$MA..Year)
data$PR <- paste(data$PR..First.author, ", ", data$PR..Year)
data$MA_PR <- paste(data$MA, ", ", data$PR)

# exclude those studies that are not included in synthesis
data$note <- NA

for(i in 1:nrow(data)) {
  if(data$MA..First.author[i] == "Zuo W" | data$MA..First.author[i] == "Cai Q") {
    if(data$Outcome[i] != "esophageal cancer" & data$Outcome[i] != "gastric cancer" & data$Outcome[i] != "hepatocellular carcinoma" &
       data$Outcome[i] != "pancreatic cancer" & data$Outcome[i] != "colorectal cancer" & data$Outcome[i] != "head and neck cancer") {
      data$note[i] <- "exclude for synthesis"
    }
  }
  
  if(data$Measurement.of.effect.size[i] == "RR") {
    data$note[i] <- "exclude for synthesis"
  }
}

data <- filter(data, is.na(note) == TRUE)

# generate label for subset
data$label <- data$Outcome

for(i in 1:nrow(data)) {
  if(isTRUE(data$Group[i] == 1.2)) {
    if(data$MA..First.author[i] == "Li Y" | data$MA..First.author[i] == "Zhang L L") {
      if(data$Outcome[i] == "coronary artery disease") {
        data$label[i] <- "cad + mi"
      } else {
        data$label[i] <- "mi"
      }
    } else if(data$MA..First.author[i] == "Gu J" | data$MA..First.author[i] == "Wang Q") {
      if(data$Outcome[i] == "CHD & MI") {
        data$label[i] <- "cad + mi"
      } else if(data$Outcome[i] == "coronary artery disease" | data$Outcome[i] == "coronary heart disease") {
        data$label[i] <- "cad"
      } else {
        data$label[i] <- "mi"
      }
    }
  }
}

data$label <- recode(data$label, "essential hypertension" = "htn", "hypertension" = "htn",
                                 "esophageal cancer" = "esophageal", "gastric cancer" = "gastric",
                                 "hepatocellular carcinoma" = "hepatocellular", "pancreatic cancer" = "pancreatic",
                                 "colorectal cancer" = "colorectal", "head and neck cancer" = "head",
                                 "breast cancer" = "breast", "late-onset Alzheimer Disease" = "AD",
                                 "Alzheimer Disease" = "AD")

label <- as.data.frame(table(data$label))
label <- as.vector(label$Var1)

# subset and calculate CA, CCA
# r <- number of primary studies
# c <- number of meta-analyses
# N <- number of ticked boxes in the citation matrix
# CA <- N / (r * c)
# CCA <- (N - r) / (r * c - r)
# interpret CCA: 0-5, slight overlap; 6-10, moderate overlap; 11-15, high overlap; >15, very high overlap

data_outcome <- list()
overlap <- data.frame(row.names = c("label", "CA", "CCA", "degree"))

for(lab in label) {
  data_outcome[[lab]] <- filter(data, label == lab)
  
  c <- nrow(table(data_outcome[[lab]]$MA))
  r <- nrow(table(data_outcome[[lab]]$PR))
  N <- nrow(table(data_outcome[[lab]]$MA_PR))
  
  CA <- round(N / (r * c), 2)
  CCA <- round((N - r) / (r * c - r), 2)
  
  if(CCA <= 0.05) {
    degree <- "slight"
  } else if(CCA <= 0.1) {
    degree <- "moderate"
  } else if(CCA <= 0.15) {
    degree <- "high"
  } else if(CCA > 0.15) {
    degree <- "very high"
  }
  
  overlap <- bind_cols(overlap, c(lab, CA, CCA, degree))
}

colnames(overlap) <- overlap[1, ]

########################################### method 2: use citation matrixes ##########################################

setwd("F:/Box Sync/Duke Kunshan University Intern/6 ALDH2 SR/primary analysis/citation matrix/citation matrix 20210221")
filenames = list.files(pattern="*.csv")
overlap2 <- data.frame(row.names = c("label", "CA", "CCA", "degree"))
data_outcome2 <- list()

for(file in filenames) {
  label2 <- strsplit(file, ".csv")[[1]]
  dat <- read.csv(file)
  data_outcome2[[label2]] <- dat
  
  dat <- dat[, 2:ncol(dat)]
  for(i in 1:ncol(dat)) {
    dat[, i] <- recode(dat[, i], "x" = "1")
    dat[, i] <- as.numeric(dat[, i])
  }
  
  CA <- sum(rowSums(dat, na.rm = TRUE))/(nrow(dat)*ncol(dat))
  CCA <- (sum(rowSums(dat, na.rm = TRUE))- nrow(dat))/(nrow(dat)*ncol(dat)-nrow(dat))
  CA <- round(CA, 2)
  CCA <- round(CCA, 2)
  if(CCA <= 0.05) {
    degree <- "slight"
  } else if(CCA <= 0.1) {
    degree <- "moderate"
  } else if(CCA <= 0.15) {
    degree <- "high"
  } else if(CCA > 0.15) {
    degree <- "very high"
  }
  
  overlap2 <- bind_cols(overlap2, c(label2, CA, CCA, degree))
}

>>>>>>> 82e8cc74bb8a630e988d065e339c566795dfd7b2
colnames(overlap2) <- overlap2[1, ]
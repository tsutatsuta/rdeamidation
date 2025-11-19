# =====R code=========================
# rDeamidation: deamidation, R version
# ==============================
# 2021-08-26: Tsutaya T: Modified from Python script.
# 2022-06-20: Tsutaya T: Temporarily finalized the script.
# 2025-05-17: Tsutaya T: Fixed an error of NaN processing
#  after the application of PercentDeamidation function.
# 2025-11-13: Tsutaya T: Finalized for publication.
# 2025-11-18: Tsutaya T: Added an example for multiple inputs.
# ==============================
# This is a command-line script to calculate Asparagine and
#  Glutamine deamidation rates of peptides. Original script
#  is written in Python by David Lyon and available via Github.
# https://github.com/dblyon/deamidation
#
# There are settings at the end of this script that users need to
#  modify according to their own environment. When using this script,
#  please scroll down and change the settings in the section
#  starting with "USER INPUT".
#
# Please modify scripts in the USER INPUT section for your
#  dataset, and then run commands in the EXECUTE section,
#  after inputting scripts in the FUNCTIONS section.
#  Packages "stringr" and "tidyverse" are needed,
#  and please download them using the following command.
# ----------
# install.packages(c("stringr", "tidyverse"))
# ----------
#
# When using this R script,
#  please indicate the name of this script "rDeamidation" and
#  cite the following original publication.
#  Mackie et al. 2018. Palaeoproteomic profiling of
#   conservation layers on a 14th Century Italian wall painting.
#   Angew Chem Int 57:7369-7374.
#  https://doi.org/10.1002/anie.201713020
#  https://www.ncbi.nlm.nih.gov/pubmed/29603563
#  If applicable, please also cite the subject version of
#   this script, uploaded to Zenodo with a persistent identifier.
# ==============================
# PURPOSE
#  Calculate N and Q deamidation rates of peptides.
#
# INPUTS
#  "evidence.txt" file from MaxQuant.
#
# OUTPUTS
#  Text files showing the deamidation rates and related results.
# ==============================
# PACKAGES
require("stringr")
require("tidyverse")

# ==============================
# FUNCTIONS
MQcolnames <- function(df,
  colname.abundance = "Intensity",
  colname.charge = "Charge",
  colname.deamNQ = "Deamidation__NQ_",
  colname.modseq = "Modified_sequence",
  colname.proteins = "Leading_razor_protein",
  colname.rawfile = "Raw_file",
  colname.sequence = "Sequence"){
# Change colnames of a data frame.
#
# args:
#  df: Data frame.
#  colname.*: Charactor scalar showing the original column names.
#
# return:
#  Data frame with the changed column names.
  cn.df <- colnames(df)
  cn.df[which(cn.df == colname.abundance)] <- "Intensity"
  cn.df[which(cn.df == colname.charge)] <- "Charge"
  cn.df[which(cn.df == colname.deamNQ)] <- "Deamidation__NQ_"
  cn.df[which(cn.df == colname.modseq)] <- "Modified_sequence"
  cn.df[which(cn.df == colname.proteins)] <- "Leading_razor_protein"
  cn.df[which(cn.df == colname.rawfile)] <- "Raw_file"
  cn.df[which(cn.df == colname.sequence)] <- "Sequence"

  colnames(df) <- cn.df

  return(df)
}

Count_N_Q <- function(aaseq){
# Count N and Q residues of a given aa sequence.
#
# args:
#  aaseq: Chracter vector of 1-letter aa abbreviation sequence.
#
# return:
#  Vector showing c(N, Q) count.
  N.count <- str_count(aaseq, "N")
  Q.count <- str_count(aaseq, "Q")

  return(c(num_N = N.count, num_Q = Q.count))
}

CalcNum_N_Q <- function(df){
# Count N and Q residues in a given data frame.
#
# args:
#  df: Data frame containing aa $Sequence column.
#
# return:
#  Data frame containing N and Q counts.
#
# dependence:
#  Count_N_Q
  aaseq <- df$Sequence
  aaseq2count.N.Q <- mapply(Count_N_Q, aaseq = aaseq)

  return(aaseq2count.N.Q)
}

Count_N2D_Q2E <- function(aaseq,
  char2find = c("\\(de\\)", "\\(Deamidation \\(NQ\\)\\)")){
# Count N2D and Q2E residues of a given modified aa sequence.
#
# args:
#  aaseq: Chracter vector of modified aa abbreviation sequence.
#  char2find: Character vector showing the deamidation in aaseq.
#   Default is "(de)" and "(Deamidation (NQ))".
#   These are used in older and newer versions of MaxQuant. 
#
# return:
#  Vector showing c(N2D, Q2E) counts.
#
# dependence:
#  Package "stringr"
  # Replacing "(de)" and "(Deamidation (NQ))" to "*"
  for(i in char2find){
    aaseq <- gsub(i, "*", aaseq)
  }
  # Counting deamidated NQ
  aaseq.split <- strsplit(aaseq, "\\*")
  aa.deamidated <- str_sub(aaseq.split[[1]], start = -1, end = -1)
  N2D.count <- sum(aa.deamidated == "N")
  Q2E.count <- sum(aa.deamidated == "Q")

  return(c(num_N2D = N2D.count, num_Q2E = Q2E.count))
}

CalcRatio_N2D_Q2E <- function(df,
  char2find = c("(de)", "\\(Deamidation \\(NQ\\)\\)")){
# Calculate N2D/N and Q2E/Q ratios in a given data frame.
#
# args:
#  df: Data frame containing aa $Modified_sequence column.
#   $Deamidation__NQ_: Count of deamidation sites.
#   $num_N, $num_Q: Count of N and Q residues.
#  char2find: Character vector showing the deamidation in aaseq.
#   Default is "(de)" and "(Deamidation (NQ))".
#
# return:
#  Data frame containing N2D and Q2E counts and ratios (weighted %).
#
# dependence:
#  Count_N2D_Q2E
  # Counting deamidation sites
  aaseq.mod <- df$Modified_sequence
  aaseq2count.N2D.Q2E <- mapply(
    function(a){Count_N2D_Q2E(aaseq = a, char2find = char2find)},
    a = aaseq.mod)

  # Asserting the results
  deamidation.count1 <- df$Deamidation__NQ_
  deamidation.count2 <- apply(aaseq2count.N2D.Q2E, 2, sum)
  if(!all(deamidation.count1 == deamidation.count2,
    na.rm = TRUE)){
    print("Something wrong in counting deamidation sites")
  }

  # Calculating ratio
  ratio_N2D <- aaseq2count.N2D.Q2E[1, ] * 1.0 / df$num_N
  ratio_Q2E <- aaseq2count.N2D.Q2E[2, ] * 1.0 / df$num_Q

  result <- rbind(aaseq2count.N2D.Q2E, ratio_N2D, ratio_Q2E)
  return(result)
}

PercentDeamidation <- function(group){
# Calculate abundance-weighted perdentage of deamidated N and Q
#  (%N2D and %Q2E, respectively)
#
# args:
#  group: Data frame divided by Raw_file, Sequence, and Charge.
#   $Intensity: Ion intensity of MS measurement.
#   $ratio_N2D, $ratio_Q2E: Ratio of deamidated N or deamidated Q residues.
#
# return:
#  Data frame containing %N2D and %Q2E.
#   If ratio_N2D or ratio_Q2E is NaN, the returned value is 0, not NaN.
#    So some 0 should be replaced with NaN if the ratio is NaN.
  # Calculating abundance-weighted perdentage of N2D and Q2E.
  intensity <- sum(group$Intensity, na.rm = TRUE)
  perc.deamN.sum <- 100.0 * 
    sum(group$ratio_N2D * group$Intensity, na.rm = TRUE) / intensity
  perc.deamQ.sum <- 100.0 * 
    sum(group$ratio_Q2E * group$Intensity, na.rm = TRUE) / intensity

  perc.deamN <- rep(perc.deamN.sum, length(group[,1]))
  perc.deamQ <- rep(perc.deamQ.sum, length(group[,1]))

  return(data.frame(id = group$id,
    perc_DeamN = perc.deamN, perc_DeamQ = perc.deamQ))
}

CountNaN <- function(group){
# Count NA and NaN
#
# args:
#  group: Data frame. Column 1 is omitted from the counting.
#
# return:
#  Data frame showing the count of NA or NaN.
  group <- group[ , -1]
  result <- apply(group, 2,
    function(x){sum(!is.na(x))})

  return(data.frame(result))
}

NumPeptide_N_Q_DeamPerc <- function(df,
  fn.out = "Number_of_Peptides_per_RawFile.txt"){
# Calculate and write the number of peptides that were used for
#  calculation of %N2D and %Q2E for each raw file.
#
# args:
#  df: Data frame containing $Raw_file, $perc_DeamN, and $perc_DeamQ.
#  fn.out: Charactor showing the complete path and file name.
#
# depend:
#  Package "tidyverse"
#
# export:
#  Number of peptides that were used for the calculation.
  # Reshaping the data frame.
  df2 <- data.frame(Raw_file = df$Raw_file,
    perc_DeamN = df$perc_DeamN, perc_DeamQ = df$perc_DeamQ)
  df2.groupby <- group_by(df2, Raw_file)
  dfm <- df2.groupby %>% do(CountNaN(.))
  N_Q <- rep(c("N", "Q"), length(unique(dfm$Raw_file)))

  result <- data.frame(Raw_file = dfm$Raw_file,
    N_Q, num_peptides = dfm$result)

  # Writing the result.
  write.table(result, file = fn.out,
    quote = FALSE, sep = "\t", row.names = FALSE)
}

MeanDeamidation <- function(group, mean.or.median = "mean"){
# Calculate mean %N2D and %Q2E for a given grouped data frame.
#
# args:
#  group: Data frame containing $perc_DeamN and $perc_DeamQ.
#  mean.or.median: Character vector c("mean", "median")
#   showing the averaging method.
#
# return:
#  Data frame showing mean %N2D and %Q2E.
  if(mean.or.median == "mean"){
    perc.deamN <- mean(group$perc_DeamN, na.rm = TRUE)
    perc.deamQ <- mean(group$perc_DeamQ, na.rm = TRUE)
  }else if(mean.or.median == "median"){
    perc.deamN <- median(group$perc_DeamN, na.rm = TRUE)
    perc.deamQ <- median(group$perc_DeamQ, na.rm = TRUE)
  }else{
    print("Method doesn't exist")
  }

  return(data.frame(
    perc_DeamN = perc.deamN, perc_DeamQ = perc.deamQ))
}

AvgPSM2Pep <- function(df, mean.or.median = "mean"){
# Calculate mean %N2D and %Q2E of peptides
#  among the different charge states.
#
# args:
#  df: Data frame containing $Raw_file, $Leading_razor_protein,
#    $Sequence, $Charge, $perc_DeamN, and $perc_DeamQ columns.
#  mean.or.median: Character vector c("mean", "median")
#   showing the averaging method. Passed to MeanDeamidation().
#
# depend:
#  MeanDeamidation, package "tidyverse"
#
# return:
#  Data frame showing mean %N2D and %Q2E for each peptides.
  # Reshaping data frame.
  dfx <- data.frame(Raw_file = df$Raw_file,
    Leading_razor_protein = df$Leading_razor_protein,
    Sequence = df$Sequence, Charge= df$Charge,
    perc_DeamN = df$perc_DeamN, perc_DeamQ = df$perc_DeamQ)
  
  # Deleting redundant rows.
  dfx <- dfx %>% distinct(Raw_file, Leading_razor_protein,
    Sequence, Charge, perc_DeamN, perc_DeamQ)
  
  # Calculating mean %N2D and %Q2E per RawFile, protein,
  #  and Sequence by applying MeanDeamidation().
  dfx.groupby <- group_by(dfx,
    Raw_file, Leading_razor_protein, Sequence)
  result <- dfx.groupby %>%
    do(MeanDeamidation(., mean.or.median))

  return(result)
}
  
AvgPep2Prot <- function(df, mean.or.median = "mean"){
# Calculate mean %N2D and %Q2E of proteins.
#
# args:
#  df: Data frame containing $Raw_file, $Leading_razor_protein,
#   $perc_DeamN, and $perc_DeamQ columns.
#  mean.or.median: Character scalar c("mean", "median")
#   showing the averaging method. Passed to MeanDeamidation().
#
# depend:
#  MeanDeamidation, package "tidyverse"
#
# return:
#  Data frame showing mean %N2D and %Q2E for each proteins.
  # Reshaping the data frame.
  dfx.groupby <- group_by(df,
    Raw_file, Leading_razor_protein)

  # Calculating mean %N2D and %Q2E per RawFile and protein
  #  by applying MeanDeamidation(). 
  result <- dfx.groupby %>%
    do(MeanDeamidation(., mean.or.median))

  return(result)
}

MeanOrMedian <- function(x, mean.or.median = "mean"){
# Calculate mean or median.
#
# args:
#  x: Vector of the values.
#  mean.or.median: Character vector c("mean", "median")
#   showing the averaging method.
#
# return:
#  Averaged value.
  if(mean.or.median == "mean"){
    result <- mean(x, na.rm = TRUE)
  }else if(mean.or.median == "median"){
    result <- median(x, na.rm = TRUE)
  }else{
    print("Method doesn't exist")
  }

  return(result)
}

BootStrap <- function(group, num.bootstraps = 1000,
  mean.or.median = "mean",
  average.method = "Peptides_per_groupby"){
# Apply bootstraping to calculate error range of %N2D and %Q2E.
#
# args:
#  group: Data frame containing $Raw_file, $perc_DeamN,
#   and $perc_DeamQ columns.
#  num.bootstrap: Numeric vector setting the number of bootstraps.
#  mean.or.median: Character vector c("mean", "median")
#   showing the averaging method.
#  average.method: Character scalar c("Peptides_per_groupby",
#   "Peptides_2_Proteins_per_groupby") to select
#   average Peptides per RawFile/Folderaverage OR
#   Peptides to Proteins and average Proteins per RawFile/Folder.
#
# depend: 
#  MeanDeamidation, package "tidyverse"
#
# return:
#  Data frame showing mean %N2D and %Q2E for each proteins.
  # Setting variables.
  group.len <- length(group$Raw_file)
  group.index <- 1:group.len
  group.N <- numeric(num.bootstraps)
  group.Q <- numeric(num.bootstraps)

  # Bootstrapping.
  for(i in 1:num.bootstraps){
    random.index <- sample(group.index,
      size = group.len, replace = TRUE)
    dft <- group[random.index, ]
    if(average.method == "Peptides_2_Proteins_per_groupby"){
      dft2 <- data.frame(
        Leading_razor_protein = dft$Leading_razor_protein,
        perc_DeamN = dft$perc_DeamN,
        perc_DeamQ = dft$perc_DeamQ)
      dft.ave <- group_by(dft2, Leading_razor_protein) %>%
        do(MeanDeamidation(., mean.or.median))
      group.N[i] <- MeanOrMedian(dft$perc_DeamN, mean.or.median)
      group.Q[i] <- MeanOrMedian(dft$perc_DeamQ, mean.or.median)
    }else if(average.method == "Peptides_per_groupby"){
      group.N[i] <- MeanOrMedian(dft$perc_DeamN, mean.or.median)
      group.Q[i] <- MeanOrMedian(dft$perc_DeamQ, mean.or.median)
    }else{
      print("Method doesn't exist")
    }
  }

  dfm <- data.frame(perc_DeamN = group.N, perc_DeamQ = group.Q)
  return(dfm)
}

BootStrapPepPerPro <- function(group, num.bootstraps = 1000,
  mean.or.median = "mean"){
# Apply bootstraping to calculate error range of %N2D and %Q2E
#  for individual proteins.
#
# args:
#  group: Data frame containing $Raw_file,
#   $Leading_razor_protein, $perc_DeamN, and $perc_DeamQ columns.
#  num.bootstrap: Numeric vector setting the number of bootstraps.
#  mean.or.median: Character vector c("mean", "median")
#   showing the averaging method.
#
# depend: 
#  MeanDeamidation
#
# return:
#  Data frame showing mean %N2D and %Q2E for each proteins.
  # Setting variables.
  group.len <- length(group$Raw_file)
  group.index <- 1:group.len
  dfm <- NULL

  # Bootstrapping.
  for(i in 1:num.bootstraps){
    random.index <- sample(group.index,
      size = group.len, replace = TRUE)
    dft <- group[random.index, ]
    dft2 <- data.frame(
      Leading_razor_protein = dft$Leading_razor_protein,
      perc_DeamN = dft$perc_DeamN,
      perc_DeamQ = dft$perc_DeamQ)
    dft.ave <- group_by(dft2, Leading_razor_protein) %>%
      do(MeanDeamidation(., mean.or.median))
    dfm <- rbind(dfm, dft.ave)
  }

  return(dfm)
}

CalcMeanCI <- function(df, ci = 95){
# Calculate mean and CI for a bootstrapped result.
#
# args:
#  df: Data frame containing bootstraped values for %N2D and %Q2E.
#  ci: Percentage for the confidence interval.
#
# return:
#  Vector containing mean, SD, lower and higher CI limits.
  ci.low = (100 - ci) / 2.0
  ci.up = ci + ci.low

  df.mean <- apply(df, 2, mean, na.rm = TRUE)
  df.sd <- apply(df, 2, sd, na.rm = TRUE)
  df.percentile.low <- apply(df, 2, quantile,
    probs = ci.low / 100, na.rm = TRUE)
  df.percentile.up <- apply(df, 2, quantile,
    probs = ci.up / 100, na.rm = TRUE)
  
  result <- cbind(df.mean, df.sd,
    df.percentile.low, df.percentile.up)
  return(result)
}

CalcMeanCIPepPerPro <- function(df, ci = 95){
# Calculate mean and CI for a bootstrapped result of
#  individual proteins.
#
# args:
#  df: Data frame containing $Leading_razor_protein column and
#   bootstraped values for %N2D and %Q2E.
#  ci: Percentage for the confidence interval.
#
# depend: 
#  CalcMeanCI
#
# return:
#  Vector containing mean, SD, lower and higher CI limits.
  ci.low = (100 - ci) / 2.0
  ci.up = ci + ci.low

  protein.unique <- unique(df$Leading_razor_protein)
  protein.unique <- protein.unique[protein.unique != ""]
  protein.len <- length(protein.unique)
  ss <- NULL

  # Calculating summary statistics for each protein.
  for(i in seq(protein.unique)){
    protein.temp <- subset(df,
      df$Leading_razor_protein == protein.unique[i])
    ss.temp <- CalcMeanCI(protein.temp[ , -1])
    ss.temp2 <- data.frame(
      Leading_razor_protein = rep(protein.unique[i], 2),
      N_Q = rownames(ss.temp), ss.temp)
    ss <- rbind(ss, ss.temp2)
    rownames(ss) <- NULL
  }

  return(ss)
}

BootStrapSummarizeDeamRawfile <- function(df,
  sampling = "peptides", num.bootstraps = 1000, ci = 95,
  mean.or.median = "mean",
  average.method = "Peptides_per_groupby"){
# Wrapper function for bootstrapping.
#
# args:
#  df: Data frame containing $Raw_file,
#   $perc_DeamN, and $perc_DeamQ columns.
#  sampling: Character sclar c("peptides", "proteins")
#   showing the target of sampling.
#  num.bootstraps: Numeric scalar setting the number of bootstraps.
#   Passed to BootStrap().
#  ci: Percentage for the confidence interval.
#  mean.or.median: Character vector c("mean", "median")
#   showing the averaging method.
#  average.method: Character scalar c("Peptides_per_groupby",
#   "Peptides_2_Proteins_per_groupby") to select
#   average Peptides per RawFile/Folder OR
#   average Peptides to Proteins and average Proteins per RawFile/Folder.
#
# depend: 
#  BootStrap, CalcMeanCI
#
# return:
#  List containing bootstrapped data and summary statistics.
  # Bootstrapping and calculating the summary statistics.
  rawfile.unique <- unique(df$Raw_file)
  rawfile.unique <- rawfile.unique[rawfile.unique != ""]
  df.bs.li <- list()
  ss.bs.li <- list()
  for(i in seq(rawfile.unique)){
    group <- subset(df, df$Raw_file == rawfile.unique[i])
    df.bs.li[[i]] <- BootStrap(group, num.bootstraps,
      mean.or.median, average.method)
    ss.bs.li[[i]] <- CalcMeanCI(df.bs.li[[i]], ci = ci)
  }
  
  # Reshaping the results (list -> data.frame).
  df.bs <- NULL
  ss.bs <- NULL
  for(i in seq(rawfile.unique)){
    # For bootstrapped values.
    df.bs.temp <- data.frame(
      Raw_file = rep(rawfile.unique[i], 2 * num.bootstraps),
      N_Q = c(rep("N", num.bootstraps), rep("Q", num.bootstraps)),
      percDeam = c(df.bs.li[[i]][ , 1], df.bs.li[[i]][ , 2]))
    df.bs <- rbind(df.bs, df.bs.temp)

    # For summary statistics.
    ss.bs.temp <- data.frame(
      rep(rawfile.unique[i], 2),
      c("N", "Q"),
      ss.bs.li[[i]])
    colnames(ss.bs.temp) <- c("Raw_file", "N_Q",
      "mean", "std", "CI_low", "CI_up")
    rownames(ss.bs.temp) <- NULL
    ss.bs <- rbind(ss.bs, ss.bs.temp)
  }
  
  return(list(df = df.bs, ss = ss.bs))
}

BootStrapSummarizeDeamRawfilePepPerPro <- function(df,
  sampling = "peptides", num.bootstraps = 1000, ci = 95,
  mean.or.median = "mean"){
# Wrapper function for bootstrapping of individual proteins.
#
# args:
#  df: Data frame containing $Raw_file,
#   $Leading_razor_protein, $perc_DeamN, and $perc_DeamQ columns.
#  sampling: Character sclar c("peptides", "proteins")
#   showing the target of sampling.
#  num.bootstraps: Numeric scalar setting the number of bootstraps.
#   Passed to BootStrapPepPerPro().
#  ci: Percentage for the confidence interval.
#  mean.or.median: Character vector c("mean", "median")
#   showing the averaging method.
#
# depend: 
#  BootStrapPepPerPro, CalcMeanCIPepPerPro
#
# return:
#  List containing bootstrapped data and summary statistics.
  # Bootstrapping and calculating the summary statistics.
  rawfile.unique <- unique(df$Raw_file)
  rawfile.unique <- rawfile.unique[rawfile.unique != ""]
  df.bs.li <- list()
  ss.bs.li <- list()
  for(i in seq(rawfile.unique)){
    group <- subset(df, df$Raw_file == rawfile.unique[i])
    df.bs.li[[i]] <- BootStrapPepPerPro(group, num.bootstraps,
      mean.or.median)
    ss.bs.li[[i]] <- CalcMeanCIPepPerPro(df.bs.li[[i]], ci = ci)
  }
  
  # Reshaping the results (list -> data.frame).
  df.bs <- NULL
  ss.bs <- NULL
  for(i in seq(rawfile.unique)){
    # For bootstrapped values.
    df.bs.len <- length(df.bs.li[[i]]$Leading_razor_protein)
    df.bs.temp <- data.frame(
      Raw_file = rep(rawfile.unique[i], df.bs.len * 2),
      Leading_razor_protein = rep(
        df.bs.li[[i]]$Leading_razor_protein, 2),
      N_Q = c(rep("N", df.bs.len), rep("Q", df.bs.len)),
      percDeam = c(df.bs.li[[i]]$perc_DeamN,
        df.bs.li[[i]]$perc_DeamQ))
    df.bs <- rbind(df.bs, df.bs.temp)

    # For summary statistics.
    ss.bs.len <- length(ss.bs.li[[i]]$Leading_razor_protein)
    ss.bs.temp <- data.frame(
      rep(rawfile.unique[i], ss.bs.len),
      ss.bs.li[[i]])
    ss.bs.temp$N_Q <- rep(c("N", "Q"), ss.bs.len / 2)
    colnames(ss.bs.temp) <- c("Raw_file",
      "Leading_razor_protein", "N_Q",
      "mean", "std", "CI_low", "CI_up")
    rownames(ss.bs.temp) <- NULL
    ss.bs <- rbind(ss.bs, ss.bs.temp)
  }
  
  return(list(df = df.bs, ss = ss.bs))
}

Deamidation <- function(fn.evidence,
  output.dir = "~",
  output.suffix = "",
  char2find = c("\\(de\\)", "\\(Deamidation \\(NQ\\)\\)"),
  colname.abundance = "Intensity",
  colname.charge = "Charge",
  colname.deamNQ = "Deamidation__NQ_",
  colname.modseq = "Modified_sequence",
  colname.proteins = "Leading_razor_protein",
  colname.rawfile = "Raw_file",
  colname.sequence = "Sequence",
  ci = 95, num.bootstraps = 1000,
  sampling = "peptides",
  average.method = "Peptides_per_groupby",
  protein.bootstraps = FALSE){
# Executive function to calculate deamidation ratios.
#
# args:
#  fn.evidence: File path to "evidence.txt".
#  output.dir: Directory path for the output files.
#  output.suffix: Suffix for the output files.
#  colname.*: Charactor scalar showing the original column names.
#   Note that " " and other special characters in the original
#   column names in "evidence.txt" will be comverted to "_" in
#   Run().
#  ci: Percentage for the confidence interval.
#  sampling: Character sclar c("peptides", "proteins")
#   showing the target of sampling.
#  num.bootstraps: Numeric scalar setting the number of bootstraps.
#  average.method: Character scalar c("Peptides_per_groupby",
#   "Peptides_2_Proteins_per_groupby") to select
#   average Peptides per RawFile/Folder OR
#   average Peptides to Proteins and average Proteins per RawFile/Folder.
#  protein.bootstraps: Boolean. If TRUE, an additional output of
#   bootstraps of peptides per protein per rawfile will be given.
#   Usually there are not enough data for this to be meaningful
#   and thus this can lead to misleading results, this is NOT
#   recommended as a default analysis.
#
# depend:
#  MQcolnames, CalcNum_N_Q, CalcRatio_N2D_Q2E, PercentDeamidation,
#  NumPeptide_N_Q_DeamPerc, AvgPSM2Pep, AvgPep2Prot,
#  BootStrapSummarizeDeamRawfile,
#  package "stringr" and "tidyverse"
#
# export:
#  Number of peptides that were used for the calculation.
#  Protein deamidation rates.
#  Bootstrapped values of deamidation rates.
#  Summary statistics of the result of bootstrapping.
  # Reading and renaming the data frame.
  print("Reading file")
  df <- read.table(fn.evidence, sep = "\t",
    header = TRUE, na.strings = NA)
  df.colnames <- colnames(df)
  df.colnames2 <- gsub("\\.", "_", df.colnames)
  colnames(df) <- df.colnames2
  df <- MQcolnames(df)

  # Caluculating NQ number and ratio
  print("Calculating deamidation ratio")
  num.N.Q <- CalcNum_N_Q(df)
  df <- data.frame(df, t(num.N.Q))
  ratio.N2D.Q2E <- CalcRatio_N2D_Q2E(df, char2find)
  df <- data.frame(df, t(ratio.N2D.Q2E))

  # Calculating deamidation per RawFile, Sequence and Charge
  df.groupby <- group_by(df, Raw_file, Sequence, Charge)
  list.percdeam <- df %>%
    group_by(Raw_file, Sequence, Charge) %>%
    do(PercentDeamidation(.)) %>%
    ungroup()
  df.percdeam <- left_join(df, list.percdeam, by = c("Raw_file", "Sequence", "Charge"))

  # Replacing 0 to NaN, if ratio_N2D or ratio_Q2E is NaN.
  df.percdeam$perc_DeamN[is.na(df.percdeam$ratio_N2D)] <- NaN
  df.percdeam$perc_DeamQ[is.na(df.percdeam$ratio_Q2E)] <- NaN

  # Writing the number of N/Q peptides with deamidation percent
  NumPeptide_N_Q_DeamPerc(df.percdeam,
    fn.out = paste(
      output.dir, "/Number_of_Peptides_per_RawFile", output.suffix, ".txt",
      sep = ""))
  
  # Protein level deamidation --> without CI
  df.pep <- AvgPSM2Pep(df.percdeam)
  df.pro <- AvgPep2Prot(df.pep)

  # Writing protein deamidation percent
  fn.out <- paste(
    output.dir, "/Protein_deamidation", output.suffix, ".txt",
    sep = "")
  write.table(df.pro, file = fn.out,
    quote = FALSE, sep = "\t", row.names = FALSE)

  # Bootstrapping
  print("Bootstrapping")
  if(sampling == "peptides"){
    df.ss.bs <- BootStrapSummarizeDeamRawfile(df.pep,
      sampling = sampling, num.bootstraps = num.bootstraps,
      ci = ci, mean.or.median = "mean",
      average.method = average.method)
  }else if(sampling == "proteins"){
    df.ss.bs <- BootStrapSummarizeDeamRawfile(df.pro,
      sampling = sampling, num.bootstraps = num.bootstraps,
      ci = ci, mean.or.median = "mean",
      average.method = average.method)
  }

  # Writing the boot strapped values and deamidation
  print("Writing results")
  fn.out <- paste(
    output.dir, "/Bootstrapped_values", output.suffix, ".txt",
    sep = "")
  write.table(df.ss.bs[[1]], file = fn.out,
    quote = FALSE, sep = "\t", row.names = FALSE)
  fn.out <- paste(
    output.dir, "/Deamidation", output.suffix, ".txt",
    sep = "")
  write.table(df.ss.bs[[2]], file = fn.out,
    quote = FALSE, sep = "\t", row.names = FALSE)

  if(protein.bootstraps){
    print("Calculating protein level deamidation per RawFile using bootstraps")
    df.ss.bs.pro <- BootStrapSummarizeDeamRawfilePepPerPro(df.pep,
      sampling = sampling, num.bootstraps = num.bootstraps,
      ci = ci, mean.or.median = "mean")
    print("Writing protein bootstrap results")
    fn.out <- paste(
      output.dir, "/Bootstrapped_values_proteins", output.suffix, ".txt",
      sep = "")
    write.table(df.ss.bs.pro[[1]], file = fn.out,
      quote = FALSE, sep = "\t", row.names = FALSE)
    fn.out <- paste(
      output.dir, "/Deamidation_proteins", output.suffix, ".txt",
      sep = "")
    write.table(df.ss.bs.pro[[2]], file = fn.out,
      quote = FALSE, sep = "\t", row.names = FALSE)
  }
}

# ==============================
# USER INPUT
# Please modify objects for your dataset and objectives.
## MaxQuant evidence.txt file (absolute path)
fn.evidence <- "~/rdeamidation/evidence_PF01E.txt"

## Output directory for results (absolute path)
output.dir <- "~/rdeamidation/"

## Suffix for output files
output.suffix <- "_PF01E"

## Character for deamidation in Modified_sequence.
#   (de) for older versions of MaxQuant
#   (Deamidation (NQ)) for recent versions of MaxuQuant
char2find = c("\\(de\\)", "\\(Deamidation \\(NQ\\)\\)")

## Corresponding column names of the R data frame
#   generated from evidence.txt file.
colname.abundance <- "Intensity"
colname.charge <- "Charge"
colname.deamNQ <- "Deamidation__NQ_"
colname.modseq <- "Modified_sequence"
colname.proteins <- "Leading_razor_protein"
colname.rawfile <- "Raw_file"
colname.sequence <- "Sequence"

## Confidence Interval (default = 95)
ci <- 95

## Number of bootstrap iterations (default = 1000)
num.bootstraps <- 1000

## Sample 'peptides' or 'proteins' for
#   Confidence Interval calculation (default='peptides')
sampling <- "peptides"

## Bootstrap Peptides per Protein per RawFile. Usually there
#   are not enough data for this to be meaningful and thus,
#   this can lead to misleading results. This is NOT
#   recommended as a default analysis. (default = FALSE)
protein.bootstraps <- FALSE

# ==============================
# EXECUTE
Deamidation(fn.evidence = fn.evidence,
  output.dir = output.dir,
  output.suffix = output.suffix,
  colname.abundance = colname.abundance,
  colname.charge = colname.charge,
  colname.deamNQ = colname.deamNQ,
  colname.modseq = colname.modseq,
  colname.proteins = colname.proteins,
  colname.rawfile = colname.rawfile,
  colname.sequence = colname.sequence,
  ci = ci,
  sampling = sampling,
  num.bootstraps = num.bootstraps,
  protein.bootstraps = protein.bootstraps)

# ==============================
# FOR MULTIPLE INPUTS
## Set variables to execute multiple evidence.txt files.
fns.evidence <- fn.evidence <- c(
  "~/rdeamidation/evidence_PF01E.txt",
  "~/rdeamidation/evidence_PF02E.txt")
outputs.suffix <- c("_PF01E", "_PF02E")

## Perform multiple calculations.
for(i in seq_along(fns.evidence)){
  Deamidation(fn.evidence = fns.evidence[i],
    output.dir = output.dir,
    output.suffix = outputs.suffix[i])
}


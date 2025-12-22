#===============================================================================
# NHANES Data Preparation for Glycohemoglobin Prediction
#===============================================================================
# This script prepares a pre-pandemic NHANES dataset for predictive modeling
# of glycohemoglobin and diabetes-related outcomes.
#
# Main steps:
# 1. Load and merge NHANES tables based on a predefined list of variables.
# 2. Recode key covariates (e.g., smoking status, ancestry).
# 3. Derive clinically relevant variables:
#    - non-HDL cholesterol
#    - prediabetes indicator (HbA1c ≥ 5.7%)
#    - estimated glomerular filtration rate (eGFR, CKD-EPI 2021 equation)
# 4. Apply inclusion and exclusion criteria consistent with the study design:
#    - fasting participants only
#    - adults (≥ 18 years)
#    - exclusion of participants on antidiabetic medication
# 5. Select variables used for prediction modeling and remove observations
#    with missing values.
# 6. Export the final, complete-case dataset for analysis
#
# Output:
# - A cleaned and fully observed NHANES dataset saved as
#   'nhanes_prepandemic_complete.csv'
#
# Authors: Lea Kronziel, Silke Szymczak
#===============================================================================



#------------------------------------------------------------------------------
# Load required packages
#------------------------------------------------------------------------------
library(knitr)
library(rmdformats)
library(rio)
library(gtsummary)
library(R.utils)   
library(dplyr)
library(kableExtra)
library(ggpubr)
library(this.path)
library(nhanesA)
library(kidney.epi)


#------------------------------------------------------------------------------
# Define input files
#------------------------------------------------------------------------------
project.dir = this.dir()
input.dir = file.path(project.dir, "external_data")

info.var.file <- file.path(
  input.dir,
  "nhanes_used_variables.csv"
)

info.drug.file <- file.path(
  input.dir,
  "RXQ_DRUG.xpt"
)


#------------------------------------------------------------------------------
# Load variable metadata
#------------------------------------------------------------------------------
info.var <- import(info.var.file)


#------------------------------------------------------------------------------
# Load and merge NHANES tables
#------------------------------------------------------------------------------
# Each NHANES table is loaded and merged by participant ID (SEQN)
tables <- unique(info.var$tablename)
dat.all <- NULL

for (tab in tables) {
  message("Loading table ", tab, " ...")
  
  # Select variables used from the current table
  cols.use <- c(
    "SEQN",
    subset(info.var, tablename == tab, "varname", drop = TRUE)
  )
  
  temp <- nhanes(paste0("P_", tab))
  
  if (is.null(dat.all)) {
    dat.all <- temp[, cols.use]
  } else {
    dat.all <- merge(
      dat.all,
      temp[, cols.use],
      all.x = TRUE,
      all.y = TRUE
    )
  }
}


#------------------------------------------------------------------------------
# Recode covariates
#------------------------------------------------------------------------------

# Smoking status
smoking <- rep(NA, nrow(dat.all))
smoking[dat.all$SMQ020 == "No"] <- "never"
smoking[
  dat.all$SMQ020 == "Yes" &
    dat.all$SMQ040 %in% c("Every day", "Some days")
] <- "current"
smoking[
  dat.all$SMQ020 == "Yes" &
    dat.all$SMQ040 == "Not at all"
] <- "former"

dat.all$smoking <- smoking


# Ancestry (restrict to major groups)
ancestry <- as.character(dat.all$RIDRETH3)
ancestry[!(ancestry %in% c("Non-Hispanic White", "Non-Hispanic Black"))] <- "Other"
dat.all$RIDRETH3 <- ancestry


#------------------------------------------------------------------------------
# Derive clinical variables
#------------------------------------------------------------------------------

# Non-HDL cholesterol
dat.all$non_hdl <- dat.all$LBXTC - dat.all$LBDHDD

# Prediabetes indicator (HbA1c ≥ 5.7%)
dat.all$prediabetes <- as.numeric(dat.all$LBXGH >= 5.7)

# Estimated glomerular filtration rate (eGFR)
# CKD-EPI Creatinine Equation (2021)
dat.all$egfr <- egfr.ckdepi.cr.2021(
  creatinine = dat.all$LBXSCR,
  age = dat.all$RIDAGEYR,
  sex = dat.all$RIAGENDR,
  creatinine_units = "mg/dL"
)


#------------------------------------------------------------------------------
# Apply inclusion and exclusion criteria
#------------------------------------------------------------------------------

# Restrict to fasting participants
dat.use <- subset(dat.all, WTSAFPRP > 0)

# Exclude participants younger than 18 years
ind.rm.age <- which(dat.use$RIDAGEYR < 18)
ids.rm.age <- unique(dat.use$SEQN[ind.rm.age])

# Exclude participants on antidiabetic medication
info.drugs <- import(info.drug.file)

id.antidiab <- unique(
  subset(
    info.drugs,
    RXDICN1B == "ANTIDIABETIC AGENTS",
    "RXDDRGID",
    drop = TRUE
  )
)

ind.rm.diab.med <- which(dat.use$RXDDRGID %in% id.antidiab)
ids.rm.diab.med <- unique(dat.use$SEQN[ind.rm.diab.med])

# Combine exclusion indices
ind.rm <- unique(c(
  ind.rm.age,
  ind.rm.diab.med
))

dat.use <- dat.use[-ind.rm, ]


#------------------------------------------------------------------------------
# Reduce dataset to variables used for prediction modeling
#------------------------------------------------------------------------------

# Rename variables using descriptive labels
ind.rename <- which(!is.na(info.var$label))

dat.pred.part <- dat.use[, info.var$varname[ind.rename]]
colnames(dat.pred.part) <- info.var$label[ind.rename]

dat.pred <- unique(
  cbind(
    dat.use[, c("SEQN", "smoking", "egfr", "non_hdl", "prediabetes")],
    dat.pred.part
  )
)


#------------------------------------------------------------------------------
# Create complete-case dataset
#------------------------------------------------------------------------------
dat.pred.complete <- na.omit(dat.pred)


#------------------------------------------------------------------------------
# Export prepared dataset
#------------------------------------------------------------------------------
export(
  dat.pred.complete,
  file = file.path(
    data.dir,
    "nhanes_prepandemic_complete.csv"
  )
)

paste0(nrow(dat.pred.complete), " participants left")

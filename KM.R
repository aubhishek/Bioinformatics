# Install and load TCGAbiolinks library
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

'''# Define the query for NSCLC (LUAD and LUSC are the TCGA codes for the two subtypes of NSCLC)
query <- GDCquery(
  project = c("TCGA-LUAD", "TCGA-LUSC"),
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  workflow.type = "Biospecimen"
)


# Download the data
GDCdownload(query)

# Load the data into R
clinical_data <- GDCprepare(query)

# Now you can proceed to filter and prepare the data for analysis'''

# Query and download clinical data for TCGA-LUAD
luad_data <- GDCquery_clinic(project = "TCGA-LUAD")
colnames(luad_data)

# Query and download clinical data for TCGA-LUSC
lusc_data <- GDCquery_clinic(project = "TCGA-LUSC")
colnames(lusc_data)

# Load necessary libraries
library(survival)
install.packages("survminer")
library(survminer)

# Select only the columns of interest
luad_data <- luad_data[, c("days_to_last_follow_up", "vital_status", "tumor_grade")]
lusc_data <- lusc_data[, c("days_to_last_follow_up", "vital_status", "tumor_grade")]

# Combine the data from both projects
combined_data <- rbind(luad_data, lusc_data)

# Convert vital_status to a binary censoring indicator
combined_data$censor <- ifelse(combined_data$vital_status == "Dead", 1, 0)


# Ensure the relevant columns are in the correct format
combined_data$days_to_last_follow_up <- as.numeric(combined_data$days_to_last_follow_up)  # Replace 'time' with the actual column name for survival time
#combined_data$vital_status <- as.numeric(combined_data$vital_status)  # Replace 'status' with the actual column name for censoring status
combined_data$tumor_grade <- as.factor(combined_data$tumor_grade)  # Replace 'grade' with the actual column name for cancer grade

head(combined_data)
write.csv(combined_data, file = "combined_data.csv", row.names = FALSE)


# Create a Surv object
surv_obj <- Surv(time = combined_data$days_to_last_follow_up, event = combined_data$censor)

# Fit the Kaplan-Meier survival curve for each grade
km_fit <- survfit(surv_obj ~ combined_data$tumor_grade)

# Plot the Kaplan-Meier survival curve
ggsurvplot(km_fit, data = combined_data,
           xlab = "Time (days)",
           ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curve by Cancer Grade",
           legend.title = "Grade",
           legend.labs = levels(combined_data$tumor_grade))
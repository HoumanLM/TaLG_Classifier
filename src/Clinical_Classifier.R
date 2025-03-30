# Libraries to load
library(survival)      
library(survminer)     
library(glmnet)        
library(dplyr)         
library(survcomp)      
library(caret)         
library(limma)         
library(sva)           
library(plsRcox)



# Load data
uromol_data <- readRDS("../data/UROMOL_TaLG.teachingcohort.rds")
knowles_data <- readRDS("../data/knowles_matched_TaLG_final.rds")


## Training cohort (uromol)

# Filter out samples with missing Recurrence data
uromol_data <- uromol_data[!is.na(uromol_data$Recurrence), ]

# Set row names using UROMOL.ID
rownames(uromol_data) <- uromol_data$UROMOL.ID

# Convert some variables to factors
uromol_data <- uromol_data %>%
    mutate(
        Sex = as.factor(Sex),
        Smoking = as.factor(Smoking),
        Tumor.stage = as.factor(Tumor.stage),
        Tumor.grade = as.factor(Tumor.grade),
        Concomitant.CIS = as.factor(Concomitant.CIS),
        Incident.tumor = as.factor(Incident.tumor),
        EAU.risk = as.factor(EAU.risk),
        BCG = as.factor(BCG),
        UROMOL2021.classification = as.factor(UROMOL2021.classification)
    )

# Create the survival object
uromol_data$Recurrence <- as.numeric(uromol_data$Recurrence)
surv_obj <- with(uromol_data, Surv(RFS_time, Recurrence))

# Extract common clinical features
clinical_vars <- uromol_data %>%
    select(Age, Sex, Tumor.stage, Tumor.grade, Concomitant.CIS, BCG, UROMOL2021.classification)

# Remove columns with only one value
single_level_vars <- sapply(clinical_vars, function(x) is.factor(x) && nlevels(x) < 2)
if(any(single_level_vars)) {
    clinical_vars <- clinical_vars[, !single_level_vars, drop = FALSE]
}

# Create the matrix
clinical_mat <- model.matrix(~ . - 1, data = clinical_vars)

# Extract gene expression data 
gene_expr <- uromol_data$exprs
if (is.matrix(gene_expr) && is.null(rownames(gene_expr))) {
    rownames(gene_expr) <- rownames(uromol_data)
}

# Perform normalization
gene_expr_norm <- normalizeBetweenArrays(gene_expr, method = "quantile")



## Training cohort (uromol)
# Variables to factors
knowles_data <- knowles_data %>%
    mutate(
        Sex = as.factor(Sex),
        Tumor.stage = as.factor(Tumor.stage),
        Tumor.grade = as.factor(Tumor.grade),
        Concomitant.CIS = as.factor(Concomitant.CIS),
        BCG = as.factor(BCG),
        UROMOL2021.classification = as.factor(UROMOL2021.classification)
    )

# Filter out samples with missing Recurrence data
knowles_data <- knowles_data %>% filter(!is.na(Recurrence))
knowles_data$Recurrence <- as.numeric(knowles_data$Recurrence)

# Set row names using UROMOL.ID
rownames(knowles_data) <- knowles_data$UROMOL.ID

# Extract clinical features.
clinical_vars_knowles <- knowles_data %>%
    select(Age, Sex, Tumor.stage, Tumor.grade, Concomitant.CIS, BCG, UROMOL2021.classification)
rownames(clinical_vars_knowles) <- rownames(knowles_data)

# Remove columns with only one value
single_level_vars <- sapply(clinical_vars_knowles, function(x) is.factor(x) && nlevels(x) < 2)
if(any(single_level_vars)){
    clinical_vars_knowles <- clinical_vars_knowles[, !single_level_vars, drop = FALSE]
}

# Create the matrix
clinical_mat_knowles <- model.matrix(~ . - 1, data = clinical_vars_knowles)

# Extract gene expression data from Knowles and quantile normalize it
gene_expr_knowles <- knowles_data$exprs

if(is.list(knowles_data$exprs)){
    gene_expr_knowles <- do.call(rbind, knowles_data$exprs)
} else if(is.matrix(knowles_data$exprs)){
    gene_expr_knowles <- knowles_data$exprs
}
rownames(gene_expr_knowles) <- rownames(knowles_data)

# Perform normalization
gene_expr_knowles_norm <- limma::normalizeBetweenArrays(gene_expr_knowles, method = "quantile")



## Batch Effect Correction
# Identify common genes between the two datasets
common_genes <- intersect(colnames(gene_expr_norm), colnames(gene_expr_knowles_norm))

# Matrices to the common genes
uromol_expr_common <- gene_expr_norm[, common_genes, drop = FALSE]
knowles_expr_common <- gene_expr_knowles_norm[, common_genes, drop = FALSE]

# Combine the two datasets by rows 
combined_expr <- rbind(uromol_expr_common, knowles_expr_common)

# Transpose
combined_expr_t <- t(combined_expr)

batch <- c(rep("UROMOL", nrow(uromol_expr_common)),
           rep("Knowles", nrow(knowles_expr_common)))
# batch effect correction
combat_expr_t <- ComBat(dat = combined_expr_t, batch = batch, par.prior = TRUE, prior.plots = FALSE)

# Transpose back 
combined_expr_corrected <- t(combat_expr_t)

# Split back into UROMOL and Knowles parts
gene_expr_norm_corrected_uromol <- combined_expr_corrected[rownames(uromol_expr_common), , drop = FALSE]
gene_expr_norm_corrected_knowles <- combined_expr_corrected[rownames(knowles_expr_common), , drop = FALSE]



## Univariate Filtering and PLS for UROMOL
# Use only the aligned UROMOL samples 
common_ids <- intersect(rownames(clinical_mat), rownames(gene_expr_norm_corrected_uromol))

# Expression data to aligned samples
expr_for_filtering <- gene_expr_norm_corrected_uromol[common_ids, , drop = FALSE]

# Perform univariate Cox regression
gene_pvals <- apply(expr_for_filtering, 2, function(gene) {
    res <- tryCatch({
        fit <- suppressWarnings(coxph(Surv(uromol_data[common_ids, "RFS_time"],
                                           uromol_data[common_ids, "Recurrence"]) ~ gene))
        pval <- summary(fit)$coefficients[,"Pr(>|z|)"][1]
        return(pval)
    }, error = function(e) {
        return(NA)
    })
    return(res)
})
gene_pvals <- gene_pvals[!is.na(gene_pvals)]

# Select genes with p-value < 0.01 
selected_genes <- names(gene_pvals)[gene_pvals < 0.01]

# Expression matrix to the selected genes
expr_filtered <- expr_for_filtering[, selected_genes, drop = FALSE]

# Perform Partial Least Squares (PLS) regression 
n_lv <- 15  
pls_model <- plsRcox(
    X = expr_filtered,
    time = uromol_data[common_ids, "RFS_time"],
    cens = uromol_data[common_ids, "Recurrence"],
    nt = n_lv
)

pls_scores <- pls_model$tt[, 1:n_lv, drop = FALSE]

# Combine clinical features 
clinical_mat_aligned <- clinical_mat[common_ids, , drop = FALSE]
X_train_new <- cbind(clinical_mat_aligned, pls_scores)

uromol_data_aligned <- uromol_data[common_ids, , drop = FALSE]
surv_obj_aligned <- with(uromol_data_aligned, Surv(RFS_time, Recurrence))




### Model: Cox Regression
set.seed(123)  
cvfit <- cv.glmnet(
    x = X_train_new,
    y = surv_obj_aligned,
    family = "cox",
    alpha = 1  
)
best_lambda <- cvfit$lambda.min
final_model <- glmnet(
    x = X_train_new,
    y = surv_obj_aligned,
    family = "cox",
    alpha = 1,
    lambda = best_lambda
)

## Internal Validation
risk_scores <- as.vector(predict(final_model, newx = X_train_new, type = "link"))

# Risk score 
uromol_data_aligned$risk_score <- risk_scores

# Quantile was considered to have more sensitivity
median_score <- quantile(uromol_data_aligned$risk_score, probs = 0.22) 
uromol_data_aligned$risk_group <- ifelse(risk_scores > median_score, "High", "Low")

### Performance matrices
# Kaplan-Meier curve
km_fit <- survfit(Surv(RFS_time, Recurrence) ~ risk_group, data = uromol_data_aligned)
ggsurvplot(
    km_fit,
    data = uromol_data_aligned,
    pval = TRUE,
    risk.table = TRUE,
    title = "Kaplan-Meier Survival Curves (UROMOL Cohort)",
    xlab = "Time (months)",
    ylab = "Recurrence-Free Survival Probability"
)

# c-index
c_index_train <- concordance.index(
    x = risk_scores,
    surv.time = uromol_data_aligned$RFS_time,
    surv.event = uromol_data_aligned$Recurrence,
    method = "noether"
)
cat("UROMOL Cohort c-index:", c_index_train$c.index, "\n")

# Confusion matrix
conf_matrix_train <- confusionMatrix(
    factor(uromol_data_aligned$risk_group, levels = c("Low", "High")),
    factor(ifelse(uromol_data_aligned$Recurrence == 1, "High", "Low"), levels = c("Low", "High"))
)
print(conf_matrix_train)




### Prediction on Knowles 
# Align clinical and gene expression data
common_ids_knowles <- intersect(rownames(clinical_mat_knowles), rownames(gene_expr_norm_corrected_knowles))
clinical_mat_knowles_aligned <- clinical_mat_knowles[common_ids_knowles, , drop = FALSE]

# Select the same genes as used in training
common_genes_selected <- intersect(selected_genes, colnames(gene_expr_norm_corrected_knowles))
common_genes_selected <- selected_genes[selected_genes %in% common_genes_selected]
expr_knowles_filtered <- gene_expr_norm_corrected_knowles[common_ids_knowles, common_genes_selected, drop = FALSE]

# PLS model
pls_proj <- predict(pls_model, newdata = expr_knowles_filtered, type = "scores")

# Check
if(nrow(pls_proj) == length(common_ids_knowles)){
    rownames(pls_proj) <- common_ids_knowles
} else {
    warning("Mismatch between number of predicted scores and common_ids_knowles; using rownames from expr_knowles_filtered.")
    rownames(pls_proj) <- rownames(expr_knowles_filtered)
}
pls_scores_knowles <- pls_proj[, 1:n_lv, drop = FALSE]

# Combine clinical features with scores
X_knowles_new <- cbind(clinical_mat_knowles_aligned, pls_scores_knowles)
knowles_data_aligned <- knowles_data[common_ids_knowles, , drop = FALSE]

# Check same columns as X_train_new
missing_cols <- setdiff(colnames(X_train_new), colnames(X_knowles_new))
if(length(missing_cols) > 0){
    add_mat <- matrix(0, nrow = nrow(X_knowles_new), ncol = length(missing_cols),
                      dimnames = list(rownames(X_knowles_new), missing_cols))
    X_knowles_new <- cbind(X_knowles_new, add_mat)
}
X_knowles_new <- X_knowles_new[, colnames(X_train_new), drop = FALSE]

# Predict risk scores 
risk_scores_knowles <- as.vector(predict(final_model, newx = X_knowles_new, type = "link"))
knowles_data_aligned$risk_group <- ifelse(risk_scores_knowles > median_score, "High", "Low")
knowles_data_aligned$risk_score <- risk_scores_knowles

## Performance matrices
km_fit_knowles <- survfit(Surv(RFS_time, Recurrence) ~ risk_group, data = knowles_data_aligned)
ggsurvplot(
    km_fit_knowles,
    data = knowles_data_aligned,
    pval = TRUE,
    risk.table = TRUE,
    title = "Kaplan-Meier Survival Curves (Knowles Cohort)",
    xlab = "Time (months)",
    ylab = "Recurrence-Free Survival Probability"
)

# c-index
df_cindex <- data.frame(
    score = risk_scores_knowles,
    time = knowles_data_aligned$RFS_time,
    event = knowles_data_aligned$Recurrence
)
# Remove rows with any NAs
df_cindex <- df_cindex[complete.cases(df_cindex), ]
# Compute c-index
c_index_knowles <- concordance.index(
    x = df_cindex$score,
    surv.time = df_cindex$time,
    surv.event = df_cindex$event,
    method = "noether"
)
cat("Knowles Cohort c-index:", c_index_knowles$c.index, "\n")

# Confusion matrix
conf_matrix_knowles <- confusionMatrix(
    factor(knowles_data_aligned$risk_group, levels = c("Low", "High")),
    factor(ifelse(knowles_data_aligned$Recurrence == 1, "High", "Low"), levels = c("Low", "High"))
)
print(conf_matrix_knowles)



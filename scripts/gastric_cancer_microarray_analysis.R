# All the necessary packages for this analysis
library(GEOquery)
library(limma)
library(corrplot)
library(MASS)

# ----------------------- DATA ACQUISITION & FORMATTING ------------------------

# Load series and platform data from GEO
cancer_dataset<-getGEO('GSE33335', GSEMatrix = TRUE, getGPL = FALSE) 
if (length(cancer_dataset) > 1) {
  idx<-grep('GPL5175', attr(gset, 'names')) 
} else {
  idx<-1
}
cancer_data<-cancer_dataset[[idx]]
exprs_data<-exprs(cancer_data)

# Let's assign class labels to the samples
num_patients<-25
class_labels<-c(rep(paste0('PT', 1:num_patients, '_Nor'), each = 1), rep(paste0('PT', 1:num_patients, '_Can'), each = 1)); colnames(exprs_data)<-class_labels
exprs_data<-as.data.frame(exprs_data)

# Now print the data table for examination
dim(exprs_data)
exprs_data

# ----------------------------- OUTLIER TESTING --------------------------------

# Reset plot window if re-running script 
par(mfrow = c(1, 1))

# Correlation matrix (heatmap)
gastric_corr_plot<-cor(exprs_data, method='pearson', use = 'pairwise.complete.obs')
corrplot(gastric_corr_plot, method = 'color', title = 'Correlation Matrix of Gene Expression in\nNormal vs. Tumor Gastric Tissues', tl.cex = 0.7, tl.col = 'black', mar = c(0,0,3,0))

# Hierarchical Clustering Dendrogram
exprs_data_tpose<-t(exprs_data)
exprs_data_dist<-dist(exprs_data_tpose, method = 'euclidean')
exprs_data_clust<-hclust(exprs_data_dist, method = 'single')
plot(exprs_data_clust, labels=names(exprs_data_tpose), main = 'Cluster Tree of 
     Gene Expression in\nNormal vs. Tumor Gastric Tissues', xlab = 'PatientNumber_DiseaseState', cex = 0.75)

# CV vs. mean plot
exprs_data_means<-apply(exprs_data, 2, mean)
exprs_data_sd<-apply(exprs_data, 2, sd)
exprs_data_cv<-(exprs_data_sd / exprs_data_means)
xlim_range<-c(3.995, 4.4)
ylim_range<-c(0.23, 0.32)
plot(exprs_data_means, exprs_data_cv, main = 'Mean vs. CV Plot of Gene Expression in\nNormal vs. Tumor Gastric Tissues', xlab = 'Mean', ylab = 'Coefficient of Variation (CV)', col = 'red', xlim = xlim_range, ylim = ylim_range)
text(exprs_data_means, exprs_data_cv, labels = colnames(exprs_data), pos = 1, cex = 0.7, col = 'blue')

# Average correlation matrix
exprs_data_avg<-apply(gastric_corr_plot, 1, mean)
par(mar = c(4, 4, 3, 2) + 0.1)
plot(c(1, length(exprs_data_avg)), range(exprs_data_avg), type = 'n', ylab = 'Average', xlab = 'PatientNumber_DiseaseState', main = 'Average Correlation of Gene Expression in\nNormal vs. Tumor Gastric Cells', xaxt = 'n')
points(exprs_data_avg, bg = 'black', col = 1, pch = 21, cex = 1.25)
text(1:length(exprs_data_avg), exprs_data_avg, labels=colnames(exprs_data), pos = 1, offset = 0.5, cex = 0.7, col = 'blue')
axis(2)
abline(v = seq(0.5,62.5,1), col = 'gray')

# I conclude that PT16_Nor is an outlier; removing from dataset
exprs_data_trimmed<-exprs_data[, !colnames(exprs_data) %in% 'PT16_Nor']

# ---------------------- EXPRESSION ANALYSIS & FILTERING -----------------------

# Re-define groups for the trimmed dataset
normal_samples<-exprs_data_trimmed[, 1:24] # Now has 24 samples
cancer_samples<-exprs_data_trimmed[, 25:49]

# Also keep track of indices
normal_indices<-1:24
cancer_indices<-25:49

# Calculate group means and execute the t-test
normal_means<-apply(normal_samples, 1, mean, na.rm = TRUE)
cancer_means<-apply(cancer_samples, 1, mean, na.rm = TRUE)

# Calculate log2 fold changes, report names/numbers of probes meeting threshold
log2_fold_changes<-cancer_means - normal_means
linear_fold_changes<-2^log2_fold_changes
probes_meeting_fold<-abs(linear_fold_changes) >= 2
altered_genes<-names(which(probes_meeting_fold))
num_altered_genes<-length(altered_genes)

# ---------------------- P-VALUE CALCULATION & FILTERING -----------------------

# Student's t-test function
t_test_all_genes <- function(x, s1, s2) {
  x1<-x[s1]
  x2<-x[s2]
  x1<-as.numeric(x1)
  x2<-as.numeric(x2)
  t.out<-t.test(x1, x2, alternative = 'two.sided', var.equal = TRUE)
  out<-as.numeric(t.out$p.value)
  return(out)
}

# Run the t-test function
computed_p_vals<-apply(exprs_data_trimmed, 1, t_test_all_genes, s1 = normal_indices, s2 = cancer_indices)

# Do a Bonferroni correction to minimize false positives
bonferroni_alpha<-(0.05 / length(computed_p_vals))
significant_genes_bonferroni<-sum(computed_p_vals <= bonferroni_alpha)

# Now get the genes that meet both the Bonferroni alpha and fold-change cutoff
probes_meeting_alpha<-computed_p_vals <= bonferroni_alpha
probes_meeting_fold<-abs(linear_fold_changes) >= 2
probes_meeting_alpha_and_fold<-(probes_meeting_alpha & probes_meeting_fold)
probes_meeting_both<-sum(probes_meeting_alpha_and_fold)

# Summary stats on the filtering process
print('Starting number of probes:'); length(computed_p_vals)
print('Number of probes meeting |log2FC| >= 2:'); num_altered_genes
print('Number of probes meeting Bonferroni alpha'); significant_genes_bonferroni
print('Number of probes meeting both FC and P-Val criteria:'); probes_meeting_both

# Histogram of Bonferroni p-values, -10log(p)
par(mfrow = c(1, 2))
bonferroni_pvals<-computed_p_vals[probes_meeting_alpha]
hist(-log10(bonferroni_pvals), col = 'lightblue', xlab = '-log10(P-Values)', main = 'Distribution of P-Values for Genes Meeting\nBonferroni-Corrected Alpha Thresholds', cex.main = 0.9, ylim = c(0, 900))
text(x = 12, y = 620, labels = "N = 2,914\np <= 2.271591E-06", pos = 4, cex = 1.0, col = 'darkred')

# Histogram of Bonferroni p-values also meeting fold-change watermark, -10log(p)
bonferroni_fold_change_pvals<-computed_p_vals[probes_meeting_alpha_and_fold]
hist(-log10(bonferroni_fold_change_pvals), col = 'lightblue', xlab = '-log10(P-Values)', main = 'Distribution of Genes Meeting Both\nSignificance and Fold-Change Thresholds', cex.main = 0.9, ylim = c(0, 900))
text(x = 12, y = 620, labels = 'N = 268\np <= 2.271591E-06\n|log2FC| >= 2', pos = 4, cex = 1.0, col = 'darkred')

# ------------------------- DIMENSIONALITY REDUCTION ---------------------------

# Subset data by probes that meet both criteria
candidate_cancer_genes<-exprs_data_trimmed[probes_meeting_alpha_and_fold, ]
numeric_cancer_data<-candidate_cancer_genes
numeric_cancer_data_tpose<-t(numeric_cancer_data)

# Add in the disease state variable
disease_state_vector<-c(rep('Normal', 24), rep('Cancer', 25))
disease_state<-data.frame(Sample = rownames(numeric_cancer_data_tpose), Disease_state = disease_state_vector)

# Conduct PCA test
pca_result<-prcomp(numeric_cancer_data_tpose, center = TRUE, scale. = TRUE)
pca_scores<-pca_result$x

# Frame PC data and merge it with disease state data
pca_data<-data.frame(Sample = rownames(pca_scores), PC1 = pca_scores[, 1], PC2 = pca_scores[, 2])
combined_data<-merge(disease_state, pca_data, by='Sample')

# Now plot the PCA results
par(mfrow = c(1, 1))  # Reset margins
disease_state_colors<-ifelse(combined_data$Disease_state == 'Cancer', 'firebrick2', 'dodgerblue2')
plot(combined_data$PC1, combined_data$PC2, col = disease_state_colors, xlab = 'Component 1', ylab = 'Component 2', main = 'PCA Plot of Gastric Cancer Expression Data by Disease State', pch = 19)
legend(-20, -7, legend = c('Normal Samples', 'Cancer Samples'), col = c('dodgerblue2', 'firebrick2'), pch = 19)

# PCA summary statistics
summary(pca_result)

# Scree plot to visualize PC variance
std_devs<-pca_result$sdev
variances<-std_devs^2
plot(variances, type = 'b', main = 'Scree Plot of PCA', xlab = 'Principal Component', ylab = 'Variance', pch = 19)

# ---------------------------- CLASSIFICATION ----------------------------------

# Extract disease state from sample names
sample_names<-rownames(numeric_cancer_data_tpose)
disease_state<-gsub('.*_(.*)', '\\1', sample_names)
disease_state<-as.factor(disease_state)

# Combine the numeric and class data
lda_data<-data.frame(Disease_state=disease_state, numeric_cancer_data_tpose)

# Set seed for consistent results
set.seed(42)

# Separate the features and class labels for training and test sets
sample_indices<-sample(seq_len(nrow(lda_data)), size = 0.7*nrow(lda_data))
training_set<-lda_data[sample_indices,]  # 70% for training
test_set<-lda_data[-sample_indices,]  # 30% for testing

# Rename 'Class' to 'Disease_state' for consistency
names(training_set)[names(training_set) == "Class"]<-'Disease_state'
names(test_set)[names(test_set) == 'Class']<-'Disease_state'

# Assign samples to the training or test set for prediction and testing
train_classes<-training_set$Disease_state
train_features<-training_set[, -which(names(training_set) == 'Disease_state')]
test_classes<-test_set$Disease_state
test_features<-test_set[, -which(names(test_set) == 'Disease_state')]

# Run the LDA function and predict the test set
lda_model<-lda(Disease_state ~ ., data = training_set)
test_set_prediction<-predict(lda_model, newdata = test_features)

# Output the prediction results
prediction_output<-table(Predicted = test_set_prediction$class, Actual = test_classes)
prediction_output

# Plot the discriminant function scores
plot_data<-data.frame(Score = test_set_prediction$x[, 1], Class = test_classes)
plot(plot_data$Score, col = as.factor(plot_data$Class), pch = 19, xlab = 'Discriminant Function 1', ylab = '', main = 'Discriminant Function Plot For\nNormal vs. Cancer Gastric Tissue Samples')
legend(2, -2, legend = c('Normal', 'Cancer'), col = 1:length(levels(plot_data$Class)), pch = 19)

# ---------------------- DISCRIMINANT GENE EVALUATION --------------------------

# Extract coefficients from the LDA model, convert to a data frame
lda_coefficients<-lda_model$scaling
lda_coefficients_df<-as.data.frame(lda_coefficients)
lda_coefficients_df$Gene<-rownames(lda_coefficients_df)

# Compute absolute values of coefficients for ranking
lda_coefficients_df$Absolute_Coeff<-abs(lda_coefficients_df$LD1)

# Get the top 5 positive and negative discriminant genes
top_genes<-lda_coefficients_df[order(lda_coefficients_df$Absolute_Coeff, decreasing = TRUE), ]
top_5_positive<-head(top_genes[top_genes$LD1 > 0, ], 5)
top_5_negative<-head(top_genes[top_genes$LD1 < 0, ], 5)

# Print top/bottom discriminant genes
top_5_positive; top_5_negative

# Load the GPL platform data
gpl<-getGEO('GPL5175', destdir = '.')
gpl_table<-Table(gpl)

# Remove leading X symbols from probe IDs to map them properly
top5_probe_ids_clean<-gsub('^X', '', top5_probe_ids)
bottom5_probe_ids_clean<-gsub('^X', '', bottom5_probe_ids)

# Match the probe IDs to the GPL annotation
top5_gene_info_raw<-gpl_table[gpl_table$ID %in% top5_probe_ids_clean, c('ID', 'gene_assignment')]
bottom5_gene_info_raw<-gpl_table[gpl_table$ID %in% bottom5_probe_ids_clean, c('ID', 'gene_assignment')]

top5_gene_info_raw; bottom5_gene_info_raw

# Function to extract first gene symbol and title from assignment string
parse_gene_assignment<-function(entry) {
  if (entry == '---' || is.na(entry)) {
    return(c(Symbol = NA, Title = NA))
  }
  parts<-strsplit(entry, ' /// ')[[1]]
  first<-strsplit(parts[1], ' // ')[[1]]
  gene_symbol<-ifelse(length(first) >= 2, first[2], NA)
  gene_title<-ifelse(length(first) >= 3, first[3], NA)
  return(c(Symbol = gene_symbol, Title = gene_title))
}

# Apply to top 5 and bottom 5
top5_clean<-t(apply(as.matrix(top5_gene_info_raw$gene_assignment), 1, parse_gene_assignment))
bottom5_clean<-t(apply(as.matrix(bottom5_gene_info_raw$gene_assignment), 1, parse_gene_assignment))

# Convert to dfs and add probe IDs
top5_clean_df<-data.frame(ID = top5_gene_info_raw$ID, top5_clean, row.names = NULL)
bottom5_clean_df<-data.frame(ID = bottom5_gene_info_raw$ID, bottom5_clean, row.names = NULL)

# Finally, print and evaluate!
print('Top 5 Discriminant genes:')
top5_clean_df
print('Bottom 5 Discriminant genes:')
bottom5_clean_df

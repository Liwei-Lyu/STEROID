# Load required libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Alpha Diversity Analysis ----
feature_table[is.na(feature_table)] <- 0
alpha_diversity <- data.frame(
  Richness = specnumber(feature_table),
  Shannon = diversity(feature_table, index = "shannon"),
  Simpson = diversity(feature_table, index = "simpson")
)
alpha_diversity$SampleID <- rownames(alpha_diversity)

merged_alpha <- merge(alpha_diversity, metadata, by = "SampleID")
long_alpha <- pivot_longer(merged_alpha, cols = c(Richness, Shannon, Simpson), 
                           names_to = "Metric", values_to = "Value")

ggplot(long_alpha, aes(x = Timepoint, y = Value, color = Group)) +
  geom_boxplot() +
  geom_line(aes(group = Subject), alpha = 0.5) +
  facet_grid(Metric ~ Group, scales = "free_y") +
  theme_bw()

# 2. Beta Diversity Analysis ----
bray_dist <- vegdist(feature_table, method = "bray")
pcoa_result <- cmdscale(bray_dist, k = 2, eig = TRUE)
pcoa_df <- data.frame(SampleID = rownames(feature_table),
                      PC1 = pcoa_result$points[,1],
                      PC2 = pcoa_result$points[,2])

pcoa_df <- merge(pcoa_df, metadata, by = "SampleID")
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group, shape = Timepoint)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.8) +
  theme_bw() +
  labs(title = "PCoA of Bray-Curtis Distances")


# 3. Differential Analysis ----
library(Maaslin2)

control <- meta %>% filter(Group == 'CTL') %>% rownames()
PO <- meta %>% filter(Group == 'PO') %>% rownames()
IM <- meta %>% filter(Group == 'IM') %>% rownames()

bac_mpa_control <- Maaslin2(input_data = mpa_species[control, ],
                            input_metadata = meta[control, ],
                            analysis_method = 'LM',
                            min_prevalence = 0.2,
                            min_abundance = 0.1,
                            normalization  = "NONE",
                            transform = 'LOG',
                            output = "control_output",
                            fixed_effects  =  c("Days", 'cell_count'),
                            random_effects = 'Subject',
                            max_significance = 0.1,
                            standardize = T,
                            correction = "BH",
                            plot_scatter = T)

bac_mpa_PO <- Maaslin2(input_data = mpa_species[PO, ],
                       input_metadata = meta[PO, ],
                       analysis_method = 'LM',
                       min_prevalence = 0.2,
                       min_abundance = 0.1,
                       normalization  = "NONE",
                       transform = 'LOG',
                       output = "treatment1_output",
                       fixed_effects  =  c("Days", 'cell_count'),
                       random_effects = 'Subject',
                       max_significance = 0.1,
                       standardize = T,
                       correction = "BH",
                       plot_scatter = T)

bac_mpa_IM <- Maaslin2(input_data = mpa_species[IM, ],
                       input_metadata = meta[IM, ],
                       analysis_method = 'LM',
                       min_prevalence = 0.2,
                       min_abundance = 0.1,
                       normalization  = "NONE",
                       transform = 'LOG',
                       output = "treatment2_output",
                       fixed_effects  =  c("Days", 'cell_count'),
                       random_effects = 'Subject',
                       max_significance = 0.1,
                       standardize = T,
                       correction = "BH",
                       plot_scatter = T)

control_result <- bac_mpa_control$results %>% filter(metadata == 'Days') %>% filter(qval < 0.1)
rownames(control_result) <- control_result$feature
PO_result <- bac_mpa_PO$results %>% filter(metadata == 'Days') %>% filter(qval < 0.1)
rownames(PO_result) <- PO_result$feature
IM_result <- bac_mpa_IM$results %>% filter(metadata == 'Days') %>% filter(qval < 0.1)
rownames(IM_result) <- IM_result$feature

PO_result <- PO_result[!row.names(PO_result) %in% intersect(control_result$feature, PO_result$feature), ]
IM_result <- IM_result[!row.names(IM_result) %in% intersect(control_result$feature, IM_result$feature), ]
# association----
# Load required libraries
library(lme4)
library(broom.mixed)
library(dplyr)

# Subset to significant species from prior analysis
significant_species <- PO_sig_species
feature_table_sub <- feature_table_scaled[, significant_species]

# Prepare data - ensure matching rows
df_combined <- cbind(feature_table_sub, metadata_scaled)
df_combined <- df_combined[complete.cases(df_combined), ]

# Define model components
predictors <- host_metabolic_markers
covariates <- "fecal_cell_count"
random_effect <- "(1|SubjectID)"

# Initialize results storage
all_results <- list()

# Perform LMM for each species-predictor pair
for(species in significant_species) {
  for(predictor in predictors) {
    
    # Construct model formula
    formula <- as.formula(paste0("`", species, "` ~ `", predictor, "` + ", 
                                 covariates, " + ", random_effect))
    
    # Fit linear mixed model
    model <- lmer(formula, data = df_combined)
    
    # Extract results
    results <- tidy(model) %>%
      filter(term == predictor) %>%
      mutate(Species = species,
             Predictor = predictor,
             Estimate = estimate,
             P_value = p.value)
    
    all_results[[paste(species, predictor)]] <- results
  }
}

# Combine results and adjust p-values
results_df <- bind_rows(all_results) %>%
  mutate(Q_value = p.adjust(P_value, method = "BH"))

# Filter significant results
significant_results <- results_df %>% filter(Q_value < 0.1)
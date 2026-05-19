#Author name: Nazmus Sadat   |Matriculation number: 218203123   |Email: sada8005@uni-landau.de

#All calculations, statistics and data visualizations were conducted with R (4.3.3, R Core Team,2024).

# Libraries ####
#The required package for data analysis needs to be loaded
library(plyr)
library(dplyr)
library(stats)
require(reshape2)
library(ggplot2)
library(vegan)
library(ggpubr)
library(lme4)
library(tidyverse)
library(lsmeans)
library(MuMIn)
library(car)
library(pairwiseAdonis)
library(cowplot)
library(Rmisc)

# Print package versions
packageVersion("plyr") #‘1.8.8’
packageVersion("dplyr") #‘1.1.4’
packageVersion("stats") #‘4.3.3’
packageVersion("reshape2") #‘1.4.4’
packageVersion("ggplot2") #‘3.4.2’
packageVersion("vegan")  #‘2.6.4’
packageVersion("ggpubr") #‘0.6.0’
packageVersion("lme4") #‘1.1.35.1’
packageVersion("tidyverse") #‘2.0.0’
packageVersion("lsmeans") #‘2.30.0’
packageVersion("MuMIn") #‘1.47.5’
packageVersion("car") #‘3.1.2’
packageVersion("pairwiseAdonis") #‘0.4.1’
packageVersion("cowplot") #‘1.1.1’
packageVersion("Rmisc") #‘1.5.1’

#Title: A mesocosm experiment assessing the effects of signal crayfish and artificial light at night on emergent aquatic insects

#Load Data ####
# Set Directory 
#This line needs to be replaced with a path on your computer
setwd ("D:/Documents/Msc/Thesis/Data")

# Read the "Taxa_abundance" sheet
taxa_abundance <- read.csv("Taxa_abundance.csv", header = T, sep=";", dec=",")

# Reshape the abundance data
abundance_reshaped <- taxa_abundance %>%
  pivot_longer(cols = -c(1:5), names_to = "Taxa", values_to = "Abundance")

head(abundance_reshaped)
#################
# explanation of data set:
# Treatment:The term "treatment" refers to the group of stressors that includes: control (CO), artificial light at night (ALAN), crayfish (CR), and a combination of both (ALAN + CR).
# Week: The term "week" describes the samples taken during the stressors' exposure. The experiment ran for six weeks, which was also the length of the stressors.
# Before the stressors started, we conducted one sampling (Week_0) during this time. After that, we conducted samples twice a week for the six weeks that the experiment was in progress.
# Week 0 (pre-sampling), week 1 (short-term response), week 4, and week 6 (long-term response) times will be taken into consideration for the thesis.
# Day: Sampling was performed twice a week, on day 1 and day 2.
# Flume: The RSM has 16 flumes; we employed 4 flumes for each treatment (replication), and this is considered a random factor.
# Trap: There were 3 sample locations in each flume: up, mid, and down, and they represented different mesohabitats.
# abundance will be standardize to one square meter
# taxa identification is most commonly on family-level

#############

# Calculate the percentage of Total_Abundance for each Taxa
Taxa_percentage <- abundance_reshaped %>%
  group_by(Taxa) %>%
  summarize(Total_Abundance = sum(Abundance))

Taxa_percentage <- Taxa_percentage %>%
  mutate(PercentageOfTotal = (Total_Abundance / sum(Total_Abundance)) * 100)

# Print the updated data frame
print(Taxa_percentage)

#############

##Data standardization ####
#The emergence trap was 1 square meter, but the sampling area was not, because of the low width of the stream. 
#Thus calculate the elevation data per square meter to be able to compare the data with other researchers. 
# 0.7 m2 for the trap upstream, 0.75 m2 for midstream and 0.8 m2 for downstream need to standardize to a 1 m^2 area.
## For abundance
abundance <- abundance_reshaped %>%
  group_by(Treatment, Week, Flume, Trap, Day) %>%
  summarise(across(-Taxa, sum)) %>%
  ungroup()


# Standardize trap areas to m2
# (0.7 m2 for the trap upstream, 0.75 m2 for midstream and 0.8 m2 for downstream) need to standardize to a 1 m^2 area.
abundance$Abundance[abundance$Trap == "up"] <- abundance$Abundance[abundance$Trap == "up"] / 0.7
abundance$Abundance[abundance$Trap == "mid"] <- abundance$Abundance[abundance$Trap == "mid"] / 0.75
abundance$Abundance[abundance$Trap == "down"] <- abundance$Abundance[abundance$Trap == "down"] / 0.8

# Sum replicates by week (After standardizing trap size, sum the Day 1 and Day 2 counts for each week)
abundance <- abundance %>%
  group_by(Treatment, Week, Flume, Trap) %>%
  summarise(Abundance = sum(Abundance))

# average replicates by week (average emergence rate of the 3 traps counts for each week)
abundance <- abundance %>%
  group_by(Treatment, Week, Flume) %>%
  summarise(Abundance = sum(Abundance)/3)


# Calculate emergence rate (Week 0 was only one sampling; the duration of the collection was 3 days. Week 0 was divided by 3 days, and weeks 1, 4, and 6 were divided by 7 days to get the emergence rate per m^2 per day.)
abundance <- abundance %>%
  group_by(Treatment, Week, Flume) %>%
  summarise(Rate = ifelse(Week == "week_0", Abundance / 3, Abundance / 7))

# abundance (ind/m2/day)
head(abundance)

#############
##checking normality## 
## For abundance
t <- abundance$Rate
hist(t)

# Apply logarithmic transformation to Abundance 
log10_abundance <- log10(abundance$Rate)
hist(log10_abundance)

##############
#H1- EAIs abundance and biomass decrease with respect to the control in the presence of signal crayfish, ALAN and similar to the both stressors combined.
#H2-  abundance decrease is more pronounced after the long-term (week 6) than the short-term (week 1) in the presence of signal crayfish, followed by ALAN and then the combined stressors.
# Create  plot using the abundance values
# This plot can show how each stressor group's abundance varies over time.
#Emergence rate calculations and plots
Rate_abundance <- group.CI(Rate~ Treatment+Week, data=abundance, ci = 0.95)
Rate_abundance
# Organize the data with the specified order for each variable
treatment_order <- c("CO",  "CR", "ALAN", "ALAN+CR")
week_order <- c("Week_0","Week_1", "Week_4", "Week_6")

Rate_abundance <- Rate_abundance %>%
  mutate(Treatment = factor(Treatment, levels = treatment_order),
         Week = factor(Week, levels = week_order))

Plot_abundance <- ggplot(Rate_abundance, aes(x = Treatment)) + 
  geom_pointrange(aes(y = Rate.mean, ymin = Rate.lower, ymax = Rate.upper, shape= Treatment)) + 
  facet_wrap(~Week, nrow = 1) +
  ylab("abundance (ind/m²/d)") + 
  xlab("Treatment") + 
  ylim(-20, 200) + 
  ggtitle("a) Total abundance") + 
  theme_classic() + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 12),
    legend.position = "right",
    strip.text.x = element_text(face = "bold"),
    text = element_text(size = 14),
    legend.title = element_blank()
  )
Plot_abundance
###############
#The count data becomes a continuous variable.
#As a continuous variable, the emergence rate can be compared across different treatments, weeks, flumes, and traps because the units of measurement are now the same for all data points.
#Given these characteristics, the usage of a LMM with Gaussian distribution involves a continuous variable.
# Organize the data with the specified order for each variable
abundance <- abundance %>%
  mutate(Treatment = factor(Treatment, levels = treatment_order),
         Week = factor(Week, levels = week_order))

# Fit the generalized linear mixed-effects model with Gaussian distribution
model_abundance <- lmer(log(Rate + 1) ~ Treatment + Week + Treatment:Week  + (1|Flume), 
                        data = abundance)

# Some model diagnostics are checked first.
hist(log ((abundance$Rate) + 1))
qqnorm(resid(model_abundance))
qqline(resid(model_abundance)) 
plot(model_abundance, type = "p") # residual vs fitted plot
shapiro.test(resid(model_abundance)) # check normality

#These values indicate how well the model fits the abundance data
r.squaredGLMM(model_abundance) # r-squared value

#Investigate overall significance of model
anova(model_abundance)
#calculate p value from f value
# Degrees of freedom
df_treatment_abundance <- 3   # Degrees of freedom for Treatment
df_week_abundance <- 3        # Degrees of freedom for Week
df_interaction_abundance <- 9 # Degrees of freedom for Treatment:Week

# Given F-values
f_treatment_abundance <- 2.7355
f_week_abundance <- 43.0786
f_interaction_abundance <- 1.6930

# Calculate the p-values
p_treatment_abundance <- 1 - pf(f_treatment_abundance, df_treatment_abundance, df_interaction_abundance)
p_week_abundance <- 1 - pf(f_week_abundance, df_week_abundance, df_interaction_abundance)
p_interaction_abundance <- 1 - pf(f_interaction_abundance, df_interaction_abundance, df_week_abundance * df_treatment_abundance)

# Print the p-values
print(p_treatment_abundance)
print(p_week_abundance)
print(p_interaction_abundance)

# Detailed look at parameters
# Let’s further investigate the differences using least-square means (lsmeans)
#To test hypotheses, h1 and h2 least-square means (lsmeans) were computed for the interaction between Week and Treatment. 
#Pairwise comparisons were performed, and a detailed summary with inferential statistics, considering multiple testing adjustments, was obtained.
result_A <- lsmeans(model_abundance, ~ Week | Treatment )
result_dunn_abundance <- contrast(result_A, 'pairwise')
result_abundance <- summary(result_dunn_abundance, infer = TRUE, adjust = 'mvt')

# Display the data frame for both Treatment and  Week
pairs(result_A, simple = "each")

################################################################################################
# Read the "Biomass" sheet
taxa_biomass <- read.csv("Taxa_biomass.csv", header = T, sep=";", dec=",")
# Reshape the biomass data
biomass_reshaped <- taxa_biomass %>%
  pivot_longer(cols = -c(1:6), names_to = "Taxa", values_to = "Biomass")
head(biomass_reshaped)

# explanation of data set:
#Biomass and abundance files are the same; only the biomass file does not have pre-stressor samples (week_0).

# Calculate the percentage of Total_Biomass for each Taxa(mg)
Biomass_percentage <- biomass_reshaped %>%
  group_by(Taxa) %>%
  summarize(Total_Biomass = sum(Biomass))

Biomass_percentage <- Biomass_percentage %>%
  mutate(PercentageOfTotal = (Total_Biomass / sum(Total_Biomass)) * 100)

# Print the updated data frame(mg)
print(Biomass_percentage)
################
##Data standardization ####
#The emergence trap was 1 square meter, but the sampling area was not, because of the low width of the stream. 
#Thus calculate the elevation data per square meter to be able to compare the data with other researchers. 
# 0.7 m2 for the trap upstream, 0.75 m2 for midstream and 0.8 m2 for downstream need to standardize to a 1 m^2 area.

## For biomass
biomass <- biomass_reshaped %>%
  group_by(Treatment, Week, Flume, Trap,Day) %>%
  summarise_at(vars(-Unit,-Taxa), sum) %>%
  ungroup()

# Standardize trap areas 
#(0.7 m2 for the trap upstream, 0.75 m2 for midstream and 0.8 m2 for downstream) need to standardize to a 1 m^2 area
biomass$Biomass[biomass$Trap == "up"] <- biomass$Biomass[biomass$Trap == "up"] /0.7 
biomass$Biomass[biomass$Trap == "mid"] <- biomass$Biomass[biomass$Trap == "mid"] /0.75
biomass$Biomass[biomass$Trap == "down"] <- biomass$Biomass[biomass$Trap == "down"] /0.8

# Sum replicates by week( After standardizing trap size, sum the Day 1 and Day 2 counts for each week) 
biomass <- biomass %>%
  group_by(Treatment, Week, Flume, Trap) %>%
  summarise(Biomass = sum(Biomass))

# average replicates by week (average emergence rate of the 3 traps counts for each week)
biomass <- biomass %>%
  group_by(Treatment, Week, Flume) %>%
  summarise(Biomass = sum(Biomass)/3)

# Calculate emergence rate (Then divide the total weekly biomass by 7 days to get the emergence rate per m^2 per day)
biomass <- biomass %>%
  group_by(Treatment, Week, Flume) %>%
  summarise(Rate = Biomass/7)
#biomass (mg/m2/day)
head(biomass)
##############
# Create  plot using the biomass values
# This plot can show how each stressor group's biomass varies over time.
#Emergence rate calculations and plots
Rate_biomass <- group.CI(Rate~ Treatment+Week, data=biomass, ci = 0.95)
Rate_biomass
# Organize the data with the specified order for each variable
treatment_order <- c("CO",  "CR", "ALAN", "ALAN+CR")
week_order <- c("Week_1", "Week_4", "Week_6")

Rate_biomass <- Rate_biomass %>%
  mutate(Treatment = factor(Treatment, levels = treatment_order),
         Week = factor(Week, levels = week_order))

Plot_biomass <- ggplot(Rate_biomass, aes(x = Treatment)) + 
  geom_pointrange(aes(y = Rate.mean, ymin = Rate.lower, ymax = Rate.upper, shape= Treatment)) + 
  facet_wrap(~Week, nrow = 1) +
  ylab("biomass (mg/ m²/d)") + 
  xlab("Treatment") + 
  ylim(-5,50) + 
  ggtitle("b) Total biomass") + 
  theme_classic() + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 12),
    legend.position = "right",
    strip.text.x = element_text(face = "bold"),
    text = element_text(size = 14),
    legend.title = element_blank()
  )
Plot_biomass

##############
##checking normality## 
## For biomass
t <- biomass$Rate
hist(t)

# Apply logarithmic transformation to Biomass 
log10_biomass <- log10(biomass$Rate)
hist(log10_biomass)
###################
#The count data becomes a continuous variable.
#As a continuous variable, the emergence rate can be compared across different treatments, weeks, flumes, and traps because the units of measurement are now the same for all data points.
#Given these characteristics, the usage of a LMM with Gaussian distribution involves a continuous variable.
# Organize the data with the specified order for each variable
biomass <- biomass %>%
  mutate(Treatment = factor(Treatment, levels = treatment_order),
         Week = factor(Week, levels = week_order))
# Fit the generalized linear mixed-effects model with Gaussian distribution
model_biomass <- lmer(log(Rate + 1) ~ Treatment + Week + Treatment:Week  + (1|Flume), 
                      data = biomass)

# Some model diagnostics are checked first.
hist(log ((biomass$Rate) + 1))
qqnorm(resid(model_biomass))
qqline(resid(model_biomass)) 
plot(model_biomass, type = "p") # residual vs fitted plot
shapiro.test(resid(model_biomass)) # check normality

#These values indicate how well the model fits the biomass data
r.squaredGLMM(model_biomass) # r-squared value

#Investigate overall significance of model
anova(model_biomass)

#calculate p value from f value
# Degrees of freedom
df_treatment_biomass <- 3    # Degrees of freedom for Treatment
df_week_biomass <- 2         # Degrees of freedom for Week
df_interaction_biomass <- 6  # Degrees of freedom for Treatment:Week

# Given F-values
f_treatment_biomass <- 0.9130
f_week_biomass <- 24.3939
f_interaction_biomass <- 1.7134

# Calculate the p-values
p_treatment_biomass <- 1 - pf(f_treatment_biomass, df_treatment_biomass, df_interaction_biomass)
p_week_biomass <- 1 - pf(f_week_biomass, df_week_biomass, df_interaction_biomass)
p_interaction_biomass <- 1 - pf(f_interaction_biomass, df_interaction_biomass, df_week_biomass * df_treatment_biomass)

# Print the p-values
print(p_treatment_biomass)
print(p_week_biomass)
print(p_interaction_biomass)

# Detailed look at parameters
# Let’s further investigate the differences using least-square means (lsmeans)
#To test hypotheses, h1 and h2 least-square means (lsmeans) were computed for the interaction between Week and Treatment. 
#Pairwise comparisons were performed, and a detailed summary with inferential statistics, considering multiple testing adjustments, was obtained.
result_B <- lsmeans(model_biomass, ~ Treatment | Week  )
result_dunn_biomass <- contrast(result_B, 'pairwise')
result_biomass <- summary(result_dunn_biomass, infer = TRUE, adjust = 'mvt')

# Display the data frame for both Treatment and  Week
pairs(result_B, simple = "each")

#################################################################################################
#### Community analysis #####
#H3: emergence insect composition deviated most from control after the long-term (week 6) than the short-term (week 1) in the presence of signal crayfish, followed by ALAN and then combined (ALAN+CR) stressor. 

#The emergence trap was 1 square meter, but the sampling area was not, because of the low width of the stream. 
#Thus calculate the elevation data per square meter to be able to compare the data with other researchers. 
# 0.7 m2 for the trap upstream, 0.75 m2 for midstream and 0.8 m2 for downstream need to standardize to a 1 m^2 area.
## For abundance

# Columns to standardize (column indices 6 to 28)
cols_to_standardize <- 6:28
#This line defines the indices of the columns (6 to 28) that need to be standardized. 
#These columns likely contain numerical data that needs to be adjusted based on the values in the "Trap" column.
# Standardizing function
standardize_abundance <- function(x) {
  ifelse(taxa_abundance$Trap == "up", x / 0.7,
         ifelse(taxa_abundance$Trap == "mid", x / 0.75,
                ifelse(taxa_abundance$Trap == "down", x / 0.8, x)))
}

# Apply the standardization function to selected columns
taxa_abundance[, cols_to_standardize] <- lapply(taxa_abundance[, cols_to_standardize], standardize_abundance)

# Sum replicates by week (After standardizing trap size, sum the Day 1 and Day 2 counts for each week)
# make subset for NMDS
NMDS_abundance <- taxa_abundance %>%
  group_by(Treatment, Week, Flume, Trap) %>%
  summarise_at(vars(-Day), sum) %>%
  ungroup()

community_abundance <- NMDS_abundance 
community_abundance$Treatment = as.factor(NMDS_abundance$Treatment)
community_abundance$Week = NMDS_abundance$Week
community_abundance$Flume = NMDS_abundance$Flume
community_abundance$Trap = NMDS_abundance$Trap

# average replicates by week (average emergence rate of the 3 traps counts for each week)
community_abundance = community_abundance %>%
  group_by(Treatment, Week, Flume) %>%
  get_summary_stats(type="mean") %>% 
  pivot_wider(names_from = variable, values_from = mean)
names(community_abundance)
structure(community_abundance)

# Calculate emergence rate (Week 0 was only one sampling; the duration of the collection was 3 days. Week 0 was divided by 3 days, and weeks 1, 4, and 6 were divided by 7 days to get the emergence rate per m^2 per day)
community_abundance <- community_abundance %>%
  mutate(
    across(
      5:27,
      ~ if_else(Week == "week_0", ./3, ./7)
    )
  )

# Calculate distance matrix using Bray-Curtis distance
#The Bray-Curtis dissimilarity is a popular metric for quantifying compositional dissimilarity between ecological communities.
#Bray-Curtis dissimilarities with permutational multivariate analysis of variance (PERMANOVA) in 9999 permutations were used to further examine changes in community composition of abundance, and biomass.
abundance_distance_matrix <- vegdist(log10(community_abundance[, -(1:4)]+1), method = 'bray', permutations = 999)
#Log transformation compresses the large values and expands the smaller ones, giving more weight to the contributions of less abundant species in the final dissimilarity measure.
#Log transformation can help normalize the skewed abundance data, making the results of PERMANOVA more reliable.
abundance_dist <- as.dist(abundance_distance_matrix)
# Set a random seed
set.seed(12345) #Using the same random seed will ensure that you get consistent results each time run the analysis
### PERMANOVA analysis for different factors
#PERMANOVA analysis was performed to investigate the overall impact of Treatment, Week along with their interactions, on community composition. 
permanova_effects_abundance <- adonis2(abundance_dist ~ Treatment + Week + Treatment:Week, data = community_abundance,permutations = 999)
permanova_effects_abundance

#To test hypotheses,h3 pairwise.adonis were computed for the interaction between Week and Treatment. 
#Pairwise comparisons were performed, and a detailed summary with inferential statistics, considering multiple testing adjustments, was obtained.
# Pairwise comparisons for Treatment:Week interaction
# Set a random seed
set.seed(123) #Using the same random seed will ensure that you get consistent results each time run the analysis
# Run pairwise adonis
pairwise.adonis(abundance_dist, interaction(community_abundance$Treatment, community_abundance$Week))

# Compute NMDS
set.seed(123) #Using the same random seed will ensure that you get consistent results each time run the analysis
nmds1 <- metaMDS(abundance_dist)
nmds1$stress
#stress value of around 0.09. The stress values are stated as a “goodness-of-fit” measure, indicating a good ordination for values below 0.1 (Clarke, 1993).
# Species scores needed to plot the NMDS.
scores_dist <- as.data.frame(scores(nmds1))
scores_dist$Treatment = community_abundance$Treatment
scores_dist$Week <- community_abundance$Week 
scores_dist$Flume <- community_abundance$Flume
# Organize the data with the specified order for each variable
treatment_order <- c("CO",  "CR", "ALAN", "ALAN+CR")
week_order <- c("Week_0","Week_1", "Week_4", "Week_6")

scores_dist <- scores_dist %>%
  mutate(Treatment = factor(Treatment, levels = treatment_order),
         Week = factor(Week, levels = week_order))

plot_nmds <- scores_dist %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  stat_stars(color = "grey65") +
  stat_conf_ellipse(geom = "polygon", aes(color = Treatment, fill = Treatment), alpha = 0.65) +
  geom_point(aes(color = Treatment, fill = Treatment, shape = Treatment), size = 1.5) +
  ggpubr::stat_mean(aes(color = Treatment, fill = Treatment, shape = Treatment), size = 2, position = "identity") +
  scale_color_manual(values = c("blue", "red", "pink", "violet"),
                     labels = c("CO", "CR", "ALAN", "ALAN+CR"),
                     name = "Treatment") +
  scale_fill_manual(values = c("blue", "red", "pink", "violet"),
                    labels = c("CO", "CR", "ALAN", "ALAN+CR"),
                    name = "Treatment") +
  scale_shape_manual(values = c(16, 17, 15, 18),
                     labels = c("CO", "CR", "ALAN", "ALAN+CR"),
                     name = "Treatment") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.3, vjust = 3,
           label = paste("stress:", nmds1$stress %>% round(4)), size = 4) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        legend.position = "right",
        plot.margin = margin(t = 0.2, b = 0.2, l = 0.2, r = 0.2, unit = "cm"),
        panel.grid.major = element_line(colour = "grey90", size = 0.1),
        title = element_text(size = 10)) +
  guides(color = guide_legend(title = "Treatment"),
         fill = guide_legend(title = "Treatment"),
         shape = guide_legend(title = "Treatment")) +
  facet_grid(. ~ Week) +
  xlim(-1.00, 0.5) +
  ylim(-0.30, 0.30)

# Print the NMDS plot with facets
print(plot_nmds)

# Facet by Week and Treatment
plot_nmds + facet_grid(Treatment ~ Week )
##################################################################################################
########preparing Figure for thesis

Plot_a <- ggplot(Rate_abundance, aes(x = Treatment)) + 
  geom_pointrange(aes(y = Rate.mean, ymin = Rate.lower, ymax = Rate.upper, shape= Treatment)) + 
  facet_wrap(~Week, nrow = 1) +
  ylab("abundance (ind/m²/d)") + 
  xlab("Treatment") + 
  ylim(-20, 200) + 
  ggtitle("a) Total abundance") + 
  theme_classic() + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 12),
    legend.position = "none",
    strip.text.x = element_text(face = "bold"),
    text = element_text(size = 14),
    legend.title = element_blank()
  )
Plot_b <- ggplot(Rate_biomass, aes(x = Treatment)) + 
  geom_pointrange(aes(y = Rate.mean, ymin = Rate.lower, ymax = Rate.upper, shape= Treatment)) + 
  facet_wrap(~Week, nrow = 1) +
  ylab("biomass (mg/ m²/d)") + 
  xlab("Treatment") + 
  ylim(-5,50) + 
  ggtitle("b) Total biomass") + 
  theme_classic() + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 12),
    legend.position = "right",
    strip.text.x = element_text(face = "bold"),
    text = element_text(size = 14),
    legend.title = element_blank()
  )

ggdraw() + draw_plot(Plot_a, x = 0, y = 0, width = .45, height = 1) +  draw_plot(Plot_b, x = .50, y = 0, width = .50, height = 1)
###################
plot_nmds1 <- function(scores_dist, title, xlim_range, ylim_range) {
  ggplot(scores_dist, aes(x = NMDS1, y = NMDS2, fill = Treatment, shape = Treatment, color = Treatment)) +
    stat_stars(color = "grey65") +
    stat_conf_ellipse(geom = "polygon", alpha = 0.65) +
    geom_point(size = 1.5) +
    ggpubr::stat_mean(size = 2, position = "identity") +
    scale_color_manual(
      values = c("blue", "red", "pink", "violet"),
      labels = c("CO", "CR", "ALAN", "ALAN+CR"),
      name = "Treatment"
    ) +
    scale_fill_manual(
      values = c("blue", "red", "pink", "violet"),
      labels = c("CO", "CR", "ALAN", "ALAN+CR"),
      name = "Treatment"
    ) +
    scale_shape_manual(
      values = c(16, 17, 15, 3),
      labels = c("CO", "CR", "ALAN", "ALAN+CR"),
      name = "Treatment"
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      hjust = -0.3,
      vjust = 3,
      label = paste("stress:", nmds1$stress %>% round(4)),
      size = 4
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12, colour = "black"),
      axis.title.y = element_text(size = 12, colour = "black"),
      axis.title.x = element_text(size = 12, colour = "black"),
      legend.position = "bottom",
      plot.margin = margin(t = 0.2, b = 0.5, l = 0.2, r = 0.2, unit = "cm"),
      panel.grid.major = element_line(colour = "grey90", size = 0.1),
      title = element_text(size = 10)
    ) +
    guides(
      fill = guide_legend(title = "Treatment", nrow = 1),
      shape = guide_legend(title = "Treatment", nrow = 1),
      color = guide_legend(title = "Treatment", nrow = 1)
    ) +
    xlim(xlim_range) +
    ylim(ylim_range) +
    labs(title = title)
}

# Generate plots
p1_data <- scores_dist[scores_dist$Week == "Week_0", ]
p1 <- plot_nmds1(p1_data, "Week_0", c(-1.00, 0.10), c(-0.30, 0.10))

p2_data <- scores_dist[scores_dist$Week == "Week_1", ]
p2 <- plot_nmds1(p2_data, "Week_1", c(0.01, 0.35), c(-0.25, 0.01))

p3_data <- scores_dist[scores_dist$Week == "Week_4", ]
p3 <- plot_nmds1(p3_data, "Week_4", c(-0.70, 0.40), c(0.05, 0.40))

p4_data <- scores_dist[scores_dist$Week == "Week_6", ]
p4 <- plot_nmds1(p4_data, "Week_6", c(-0.25, 0.50), c(-0.05, 0.10))

# Extract the legend from p4
legend <- get_legend(p1)

# Remove the legend from individual plots
a1 <- p1 + theme(legend.position = "none")
b2 <- p2 + theme(legend.position = "none")
c3 <- p3 + theme(legend.position = "none")
d4 <- p4 + theme(legend.position = "none")

# Arrange the plots without the legend
plot_grid <- plot_grid(a1, b2, c3, d4, ncol = 2)

# Combine the plot grid and the legend
legend_plot <- plot_grid(
  plotlist = list(plot_grid, ggdraw() + draw_grob(legend)),
  ncol = 1,
  rel_heights = c(0.9, 0.1)
)

# Display the combined plot with the legend
legend_plot
##############
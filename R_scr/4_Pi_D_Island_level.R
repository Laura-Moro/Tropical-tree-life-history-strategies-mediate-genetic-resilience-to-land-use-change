library(emmeans)
# =============================================================================
#COMPARE Pi AND D across species and Islands 
# =============================================================================

#load the data
I_stats <- read.csv("Data/PopGen_stats_islands.csv")

# =============================================================================
# Pi: Does Island explain variation in Pi, beyond species?
# =============================================================================

# model1: pi ~ island + species
mod1pi <- lm(pi_scaled ~ Island + sp_code, data = I_stats)

# model2: pi ~ species only
mod2pi <- lm(pi_scaled ~ sp_code, data = I_stats)

# Compare: does adding Island significantly improve the fit?
anova(mod2pi, mod1pi)
### ===> Island factor is significant!!

# =============================================================================
# Pairwise comparisons across islands 
# =============================================================================
# Fit the full model, with the islands as separate levels
mod1pi_full <- lm(pi_scaled ~ Island + sp_code, data = I_stats)
summary(mod1pi_full)

# All pairwise comparisons between island levels, adjusted for multiple
emm_pop_pi <- emmeans(mod1pi_full, ~ Island)
pairs(emm_pop_pi , adjust = "none")

#Result of pairwise comparisons 
# contrast estimate    SE df t.ratio p.value
# DR - PR  -0.000178 0.000220 22  -0.808  0.4277 
# DR - VI   0.000578 0.000241 22   2.397  0.0254 --> DR is significatly different form VI !!
# PR - VI   0.000756 0.000220 22   3.434  0.0024 --> PR is significatly different form VI !!


# =============================================================================
#  Does Island explain variation in D, beyond species?
# =============================================================================

# model1: D ~ island + species
mod1D <- lm(D ~ Island + sp_code, data = I_stats)

# model2: D ~ species only
mod2D <- lm(D ~ sp_code, data = I_stats)

# Compare: does adding Island significantly improve the fit?
anova(mod2D, mod1D)
### ===> Island factor is significant (as before)

# =============================================================================
# Pairwise comparisons across islands 
# =============================================================================
# Fit the full model on the 'pi' data, with all islands as separate levels
mod1D_full <- lm(D ~ Island + sp_code, data = I_stats)
summary(mod1pi_full)

# All pairwise comparisons between island levels, adjusted for multiple
# testing (Tukey by default; use adjust = "bonferroni" or "fdr" if preferred)
emm_pop_D <- emmeans(mod1D_full, ~ Island)
pairs(emm_pop_D, adjust = "none")

#Result of pairwise comparisons 
# contrast estimate    SE df t.ratio p.value
# DR - PR     0.512 0.164 22   3.127  0.0049 --> DR is significatly different form PR!!
# DR - VI     0.331 0.179 22   1.844  0.0786 --> DR is marginally different form PR!!
# PR - VI    -0.181 0.164 22  -1.106  0.2805

# =============================================================================
# Box plots comparing Pi and D across islands 
# =============================================================================

#cose colours
col <- c ("#009E73", "#E69F00", "#CC79A7")

# Boxplot Pi
par(mar = c(3, 5, 1, 1))
boxplot(pi_scaled ~ Island,
        data = I_stats,
        col = col,
        xlab = "",       
        ylab = expression(pi), 
        cex.lab = 2.5, 
        cex.axis =1.5,
        main = "")

#Boxplot D
par(mar = c(3, 5, 1, 1))
boxplot(D ~ Island,
        data = I_stats,
        col = col,
        xlab = "",       
        ylab = "Tajima's D",
        cex.lab = 2, 
        cex.axis =1.5,
        main = "")

library(corrplot)        

#### IBD versus traits ####
#load files
ibd<-read.csv("Data/IBD_acorss_islands.csv")
traits<-read.csv("Data/Original_trait_values.csv")
traits$thk.log.z[is.na(traits$thk.log.z)] <- 0

#rank IBD slope
ibd$rank<-rank(ibd$Slope)

#merge trait data to to IBD
ibd<-merge(ibd,traits,by='sp_code')

# ---------------------------------------------------------------------------
# Subset the data to keep only sites/populations with more than 1 or more
# than 2 islands, respectively. 
# ---------------------------------------------------------------------------
pops<-read.csv("Data/PopGen_stats_islands.csv")
n_islands <- aggregate(Island ~ sp_code,
          data = pops,
          FUN = function(x) length(unique(x)))
ibd$N_Islands <- n_islands$Island[match(ibd$sp_code,n_islands$sp_code)]
ibd$N_Islands[is.na(ibd$N_Islands)] <- 1

ibd_2 <- ibd[ibd$N_Islands > 1, ]
ibd_3 <- ibd[ibd$N_Islands > 2, ]

# ---------------------------------------------------------------------------
# Spearman correlations between Slope_maf and columns 24:30, using ibd_2
# (N_Islands > 1). For each column i, store the correlation estimate (rho)
# and p-value in vectors indexed by k.
# ---------------------------------------------------------------------------
slope_p.value_2 <- NULL
slope_estimate_2 <- NULL
k <- 1
  
for (i in c(17:22)) { #"ss.log.z","thk.log.z","sla.log.z","maxht.z","wd.z","PC1"
  cor <- cor.test(ibd_2$Slope, ibd_2[, i], method = 'spearman')
  slope_estimate_2[k] <- cor$estimate   # Spearman's rho
  slope_p.value_2[k] <- cor$p.value     # associated p-value
  k <- k + 1
}

# ---------------------------------------------------------------------------
# Same as above, but using ibd_3 (N_Islands > 2).
# ---------------------------------------------------------------------------
slope_p.value_3 <- NULL
slope_estimate_3 <- NULL
k <- 1
for (i in c(17:22)) {
  cor <- cor.test(ibd_3$Slope, ibd_3[, i], method = 'spearman')
  slope_estimate_3[k] <- cor$estimate
  slope_p.value_3[k] <- cor$p.value
  k <- k + 1
}

# ---------------------------------------------------------------------------
# Same procedure, but correlating Fst_maf_mean (instead of Slope_maf)
# against columns 24:30, using ibd_2 (N_Islands > 1).
# ---------------------------------------------------------------------------
Fst_p.value_2 <- NULL
Fst_estimate_2 <- NULL
k <- 1
for (i in c(17:22)) {
  cor <- cor.test(ibd_2$mean_fst, ibd_2[, i], method = 'spearman')
  Fst_estimate_2[k] <- cor$estimate
  Fst_p.value_2[k] <- cor$p.value
  k <- k + 1
}

# ---------------------------------------------------------------------------
# Same as above, but using ibd_3 (N_Islands > 2).
# ---------------------------------------------------------------------------
Fst_p.value_3 <- NULL
Fst_estimate_3 <- NULL
k <- 1
for (i in c(17:22)) {
  cor <- cor.test(ibd_3$mean_fst, ibd_3[, i], method = 'spearman')
  Fst_estimate_3[k] <- cor$estimate
  Fst_p.value_3[k] <- cor$p.value
  k <- k + 1
}

# ---------------------------------------------------------------------------
# Combine all four sets of correlation estimates (rho) into one matrix.
# Rows = the variables in columns 24:30 of ibd; columns = each
# combination of response variable (Slope_maf / Fst_maf_mean) and island
# threshold (>1 / >2).
# ---------------------------------------------------------------------------
cor_estimate <- as.matrix(cbind(slope_estimate_2, slope_estimate_3,
                                Fst_estimate_2, Fst_estimate_3))
rownames(cor_estimate) <- colnames(ibd)[c(17:22)]
colnames(cor_estimate) <- c('IBD_2_islands', 'IBD_3_islands',
                            'Fst_2_islands', 'Fst_3_islands')

# ---------------------------------------------------------------------------
# Same layout as cor_estimate, but storing the corresponding p-values
# instead of the correlation estimates.
# ---------------------------------------------------------------------------
cor_pvalues <- as.matrix(cbind(slope_p.value_2, slope_p.value_3,
                               Fst_p.value_2, Fst_p.value_3))
rownames(cor_pvalues) <- colnames(ibd)[c(17:22)]
colnames(cor_pvalues) <- c('IBD_2_islands', 'IBD_3_islands',
                           'Fst_2_islands', 'Fst_3_islands')

# ---------------------------------------------------------------------------
# Plot the correlation matrix as a corrplot, using p-values to mark
# significance.
#   - insig = 'label_sig'  -> non-significant cells are left blank/hidden,
#                              significant ones are flagged with sig symbols
#   - method = "circle"    -> circle size/color encodes correlation strength
#   - pch.cex / pch.col    -> size and color of the significance symbols
#   - addgrid.col          -> color of the grid lines between cells
#   - type = 'full'        -> display the full matrix (not just upper/lower
#                              triangle)
#   - diag = TRUE          -> include the diagonal
#   - sig.level            -> significance thresholds (0.001, 0.01, 0.05)
#                              used to determine number of sig. symbols
#   - tl.col/tl.srt/tl.cex -> text label color, rotation angle, and size
#                              for row/column labels
# ---------------------------------------------------------------------------
corrplot(cor_estimate,
         p.mat = cor_pvalues,
         insig = 'label_sig',
         method = "circle",
         pch.cex = .7,
         pch.col = "gray0",
         addgrid.col = "gray",
         type = 'full', diag = T, sig.level = c(0.001, 0.01, 0.05),
         tl.col = 'black', tl.srt = 45, tl.cex = 0.7
)

# ---------------------------------------------------------------------------
#Plot IBD rank versus traits species occuring in more then 1 islands
# ---------------------------------------------------------------------------

# Build a 10-step color gradient from dark green to light green
cols <- colorRampPalette(c("darkgreen", "lightgreen"))(10)

# Map PC1 values to colors
# Bin PC1 into 10 equal-width intervals and get, for each row, which bin
# (1-10) it falls into; this index is used to pick a color from 'cols'
col_index <- cut(ibd_2$PC1, breaks = 10, labels = FALSE)

# Plot IBD rank versus traits (PC1) with colored points
# Each point's color reflects its PC1 value (dark green = low, light green = high)
plot(ibd_2$PC1, ibd_2$rankMaf,
     col = cols[col_index],
     pch = 16,   # solid filled circle
     cex = 3     # point size
)

# Add regression line
# Fit and overlay a linear regression of rankMaf on PC1
abline(lm(ibd_2$rankMaf ~ ibd_2$PC1), lwd = 3)


# ---------------------------------------------------------------------------
# Plot IBD rank versus traits species occuring in more then 1 islands
# ---------------------------------------------------------------------------

# Rebuild the same 10-step color gradient (recomputed here, though it's
# identical to 'cols' above since the palette/breaks count don't change)
cols <- colorRampPalette(c("darkgreen", "lightgreen"))(10)

# Map PC1 values (for ibd_3 this time) to color bins
col_index <- cut(ibd_3$PC1, breaks = 10, labels = FALSE)

# Plot IBD rank versus traits (PC1) with colored points, for ibd_3
plot(ibd_3$PC1, ibd_3$rankMaf,
     col = cols[col_index],
     pch = 16,
     cex = 3
)

# Add regression line
# NOTE: bug — this regresses ibd_3$rankMaf on ibd_2$PC1 (mismatched
# datasets). It should almost certainly be:
#   abline(lm(ibd_3$rankMaf ~ ibd_3$PC1), lwd = 3)
abline(lm(ibd_3$rankMaf ~ ibd_2$PC1), lwd = 3)


ibd_2 <- ibd[ibd$N_Islands > 2, ]

cols <- colorRampPalette(c("darkgreen", "lightgreen"))(10)

# Map PC1 values to colors
col_index <- cut(ibd_2$PC1, breaks = 10, labels = FALSE)

# Plot IBD rank versus traits (PC1) with colored points
plot(ibd_2$PC1, ibd_2$rankMaf,
     col = cols[col_index],
     pch = 16, 
     cex=3
)

ibd_2 <- ibd[ibd$N_Islands > 2, ]

cols <- colorRampPalette(c("darkgreen", "lightgreen"))(10)

# Map PC1 values to colors
col_index <- cut(ibd_2$PC1, breaks = 10, labels = FALSE)

# Plot IBD rank versus traits (PC1) with colored points
plot(ibd_2$PC1, ibd_2$rankMaf,
     col = cols[col_index],
     pch = 16, 
     cex=3
)

# Add regression line
abline(lm(ibd_2$rankMaf ~ ibd_2$PC1), lwd = 3)


#species occuring in more then 2 islands 
ibd_3 <- ibd[ibd$N_Islands > 2, ]

cols <- colorRampPalette(c("darkgreen", "lightgreen"))(10)

# Map PC1 values to colors
col_index <- cut(ibd_3$PC1, breaks = 10, labels = FALSE)

# Plot IBD rank versus traits (PC1) with colored points
plot(ibd_3$PC1, ibd_3$rankMaf,
     col = cols[col_index],
     pch = 16, 
     cex=3
)

# Add regression line
abline(lm(ibd_3$rankMaf ~ ibd_3$PC1), lwd = 3)


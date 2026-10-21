library(reshape2)
library(geosphere)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)


#geographic distance 
samp <- read.table("Data/sampling_locations.csv", header = TRUE, sep = ",")

#geogrpahic distance matrix
dist_matrix <- distm(samp[, c("longitude","latitude")], fun = distHaversine)/1000
rownames(dist_matrix) <- samp$pop
colnames(dist_matrix) <- samp$pop
dist_matrix <- round(dist_matrix, 1)

#loading FST for all species - filtered form maf 0.5 for PR DR VI 
fst <- read.table("Data/allspecies_fst_hud.txt", header = TRUE, sep = "\t")

#make one data frame of all the fst files
species_list <- c("AMYELE", "CASARB", "CASSYL", "DACEXC", "GYMLUC", "KRUFER",
                  "MICGUY", "POIFLO", "SCHMOR", "TABHET", "BURSIM", "CASGUI",
                  "CECSCH", "FAROCC", "HIRRUG", "MANBID", "MYRDEF", "PREMON",
                  "SLOBER")

mantel_results <- data.frame(
  sp_code      = character(),
  Mantel_R     = numeric(),   
  P_value      = numeric(), 
  Intercept    = numeric(),
  Slope        = numeric(),
  stringsAsFactors = FALSE
)

residuals <- data.frame(sp_code = character(),
                        pop1 = character(),
                        pop2 = character(),
                        pair = character(),
                        residuals = numeric(),
                        fst = numeric(),
                        dist = numeric(),
                        stringsAsFactors = FALSE)


  pdf("IBD_fst_maf.pdf", width = 15, height = 15)
  
  par(mfrow = c(5, 4),        # <-- change depending on how many species
      oma = c(3, 3, 5, 3),    # outer margins for global title
      mar = c(4, 4, 2, 1),
      mgp = c(1.8, 0.6, 0))
  
  for(j in species_list){
    
    species_name=j
  
    fst_all_sp <- fst %>% 
    filter(species == j)
    
    df <- fst_all_sp[,1:3] 
    if (nrow(df) == 0) {
      next
    }
  
    # Create square matrix
    fst_matrix <- acast(df, pop1 ~ pop2, value.var = "fst")
    names_cols <- colnames(fst_matrix)
    names_rows <- rownames(fst_matrix)
  
    # Reorder matrix so rownames == colnames and same order
      missing_rows <- setdiff(names_cols, names_rows)
      missing_cols <-setdiff(names_rows, names_cols)
  
    if (length(missing_rows) > 0) {
    new_rows <- matrix(
      NA,
      nrow = length(missing_rows),
      ncol = ncol(fst_matrix),
      dimnames = list(missing_rows, names_cols)
    )
    fst_matrix <- rbind(fst_matrix, new_rows)
    }

    if (length(missing_cols) > 0) {
      new_cols <- matrix(
      NA,
      nrow = nrow(fst_matrix),
      ncol = length(missing_cols),
      dimnames = list(rownames(fst_matrix), missing_cols)
       )
      fst_matrix <- cbind(fst_matrix, new_cols)
      }

    # 3. Reorder so rownames == colnames and same order
    all_names <- sort(unique(c(rownames(fst_matrix), colnames(fst_matrix))))
    fst_matrix <- fst_matrix[all_names, all_names]

    diag(fst_matrix) <- 0

    # Make full symmetric matrix
    # 1) Mirror upper triangle into lower triangle
    fst_matrix_full <- fst_matrix
    
    fst_matrix_full <- pmax(fst_matrix, t(fst_matrix), na.rm = TRUE)
    
    selected_rows <- rownames(fst_matrix_full)
    selected_cols <- colnames(fst_matrix_full)

    dist_sp <- dist_matrix[match(selected_rows, rownames(dist_matrix)),
                       match(selected_cols, colnames(dist_matrix))]

    mantel_result <- mantel(
      as.dist(dist_sp),
      as.dist(fst_matrix_full),
      method = "pearson",
      permutations = 1000
    )

    
    #make one plot for each species
    idx <- upper.tri(fst_matrix_full)
    
    fst_vals  <- fst_matrix_full[idx]
    dist_vals <- dist_sp[idx]
    
    fst_vals <- fst_vals / (1 - fst_vals)
    dist_vals <- log(dist_vals)
  
    
    # Extract pop1–pop2 labels in the same order as fst_vals
    pop1_vec <- rownames(fst_matrix_full)[row(fst_matrix_full)[idx]]
    pop2_vec <- colnames(fst_matrix_full)[col(fst_matrix_full)[idx]]
    
    pair_vec <- paste(pop1_vec, pop2_vec, sep = "-")
    
    plot(dist_vals, fst_vals,
         xlab = "Geographic distance Km",
         ylab = "FST",
         pch = 19,
         col="gray40", 
         ylim = c(0, 0.6),
         main = j)
    
    
    if (!is.na(mantel_result$signif)) {
      model <- lm(fst_vals ~ dist_vals)
      line_col <- ifelse(mantel_result$signif < 0.05, "red", "black")
      abline(model, col = line_col, lwd = 2)
    }
     

    model <- lm(fst_vals ~ dist_vals)
    
    
    intercept <- unname(coef(model)["(Intercept)"])
    slope     <- unname(coef( model)["dist_vals"])
    
    # Add results to your results data frame
    
    mantel_results <- rbind(
      mantel_results,
      data.frame(
        sp_code = species_name,
        Mantel_R = mantel_result$statistic,
        P_value = mantel_result$signif,
        Intercept = intercept,
        Slope = slope,
        stringsAsFactors = FALSE
      )
    )
    
    residuals <- rbind(
      residuals,
      data.frame(
        sp_code   = j,
        pop1      = pop1_vec,
        pop2      = pop2_vec,
        pair      = pair_vec,
        residuals = resid(model),
        fst = fst_vals,
        dist=dist_vals,
        stringsAsFactors = FALSE
      )
    )
    
}   
    
mtext("Isolation by Distance maf", 
      outer = TRUE, cex = 2, font = 2, line = 1)

##add mean Fst acprss the PR FR and VI to the mantel_results
m_fst <- fst %>%
  group_by(species) %>%
  summarise(mean_fst = mean(fst, na.rm = TRUE))

mantel_results$mean_fst <- m_fst$mean_fst[match(mantel_results$sp_code,m_fst$species)]

#save the tow dataframes
write.csv(mantel_results,"data/IBD_acorss_islands.csv", row.names = F)
write.csv(residuals,"data/residuals_IBD_acorss_islands")


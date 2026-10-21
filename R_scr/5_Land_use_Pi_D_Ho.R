require(corrplot)
require(lmtest)
library(dplyr)

Pop_stats <- read.csv("Data/PopGen_stats_pops.csv")
Pop_Habitat <- read.csv("Data/Habitat_51_00_pop.csv")
I_stats <-read.csv("Data/PopGen_stats_islands.csv")
traits<-read.csv("Data/Species_traits_habitat.csv")


##########################################################################
######################### D and Pi versus traits PR only ##################
###########################################################################

#filter the PR stats
PR_I_stats<-I_stats[which(I_stats$Island=='PR'),]
traits<-traits[,-2]
PR_I_stats<-merge(PR_I_stats,traits,by='sp_code')

#calculate delta and ratio
PR_I_stats$delta <- PR_I_stats$SH_00 - PR_I_stats$SH_51
PR_I_stats$ratio <- PR_I_stats$SH_00 / PR_I_stats$SH_51

#select the columns
cols <- PR_I_stats[,c("delta","ratio","PC1","maxht.z","sla.log.z","ss.log.z","thk.log.z","wd.z")]

#runf ro PI 
coeff_Pi <- NULL
for (colname in colnames(cols)) {
  coeff <- summary(lm(scale(PR_I_stats$pi_scaled) ~ scale(PR_I_stats[, colname])))$coefficients[2, ]
  coeff_Pi <- rbind.data.frame(coeff_Pi, coeff)
}

rownames(coeff_Pi) <- paste(colnames(cols), "_pi", sep = "")
colnames(coeff_Pi) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
coeff_Pi <- as.data.frame(coeff_Pi)

#run for D
coeff_D <- NULL
for (colname in colnames(cols)) {
  coeff <- summary(lm(scale(PR_I_stats$D) ~ scale(PR_I_stats[, colname])))$coefficients[2, ]
  coeff_D <- rbind.data.frame(coeff_D, coeff)
}

rownames(coeff_D) <- paste(colnames(cols), "_D", sep = "")
colnames(coeff_D) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
coeff_D <- as.data.frame(coeff_D)

#cobine more 
coeffout<-rbind.data.frame(coeff_Pi,coeff_D)

#reorder the raws
base_order <- sub("_pi$", "", rownames(coeffout)[1:8])
new_order <- as.vector(rbind(paste0(base_order, "_pi"),
                             paste0(base_order, "_D")))
coeffout <- coeffout[new_order, ]

n_predictors <- nrow(coeffout) / 2   # 2 rows (pi, D) per predictor
predictor_names <- sub("_pi$", "", rownames(coeffout)[seq(1, nrow(coeffout), by = 2)])

# ---------------------------------------------------------------------------
# 'test' = y-position for each row of coeffout
offset <- 0.15
test <- rep(n_predictors:1, each = 2) + rep(c(offset, -offset), n_predictors)
# n_predictors:1 (rather than 1:n_predictors) so the first predictor in
# coeffout appears at the TOP of the plot, matching typical forest-plot layout

col_D  <- 'dodgerblue2'
col_pi <- 'darkorange2'

# ---------------------------------------------------------------------------
# Set up empty plot
# ---------------------------------------------------------------------------
alpha95 <- 0.05
z95 <- qnorm(1 - alpha95 / 2)
xr <- range(c(coeffout$Estimate - z95 * coeffout$`Std. Error`,
              coeffout$Estimate + z95 * coeffout$`Std. Error`))

par(mar = c(5, 8, 2, 6))  # extra left margin for predictor labels

plot(NA, xlim = c(-1, 1.5),, ylim = c(0.5, n_predictors + 0.5),
     yaxt = 'n', ylab = '', xlab = 'standardized coefficeint', cex.axis = 1.5, cex.lab = 1.5 )

axis(2, at = n_predictors:1, labels = predictor_names, las = 1, cex.axis = 1.5)
abline(v = 0, lty = 2, col = 'grey60')

# ---------------------------------------------------------------------------
# Draw CI segments (this is your original loop, with 'test' now defined)
# ---------------------------------------------------------------------------
for (ci in 1:3) {
  
  # Calculate confidence interval for the estimate
  alpha <- c(0.05, 0.1, 0.2)[ci]   # significance level (1 - confidence level)
  z_value <- qnorm(1 - alpha / 2)  # critical value for normal distribution
  
  lower_bound <- coeffout$Estimate - z_value * coeffout$`Std. Error`
  upper_bound <- coeffout$Estimate + z_value * coeffout$`Std. Error`
  
  segments(lower_bound, test, upper_bound, test,
           col = c(col_D, col_pi),      # recycled: odd rows = D, even rows = pi
           lwd = c(1, 3, 5)[ci])
}

# ---------------------------------------------------------------------------
# Point estimates on top
# ---------------------------------------------------------------------------
points(coeffout$Estimate, test, pch = 16, cex = 1.4,
       col = rep(c(col_D, col_pi), n_predictors))

# ---------------------------------------------------------------------------
# Legends
# ---------------------------------------------------------------------------
legend("bottomright", legend = c("95% CI", "90% CI", "80% CI"),
       lwd = c(1, 3, 5), bty = "n", cex = )

legend("bottom", inset = c(0, 0.98), xpd = TRUE, horiz = TRUE,
       legend = c("Tajima's D", expression(pi)),
       col = c(col_D, col_pi), pch = 16, bty = 'n', cex = 1.3)
# ---------------------------------------------------------------------------
# Plot PC1 ~ Ration 00/51
# ---------------------------------------------------------------------------
par(mar = c(4, 4.5, 0.5, 0.5))
plot(PR_I_stats$ratio, PR_I_stats$PC1,
     pch = 19,
     xlab = expression(Ratio["00/51"]),
     ylab = "PC1", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     cex=1.5, 
     col = "gray40")
abline(lm(PR_I_stats$PC1~PR_I_stats$ratio), lwd = 2.5)

summary(lm(PR_I_stats$PC1~PR_I_stats$ratio))

# ---------------------------------------------------------------------------
#correlation plot (figure S5)
# ---------------------------------------------------------------------------
cor_PR_I_stats <- PR_I_stats[,c("delta","ratio","PC1","maxht.z","sla.log.z","ss.log.z","thk.log.z","wd.z")]

matrice_spearman <- cor(cor_PR_I_stats[,c(1:8)], method = "spearman", use = "complete.obs")
res <- cor.mtest(cor_PR_I_stats[,c(1:8)], method = "spearman",conf.level = 0.95)
#signes <- ifelse(matrice_spearman > 0, "+", "-")
corrplot(matrice_spearman,
         p.mat = res$p,
         insig = 'label_sig',  # Cache les non-significatifs
         method = "circle",
         # Fond coloré pour le signe/force             # Utilise nos symboles + et -
         pch.cex = .7,            # Taille des symboles
         pch.col = "gray0",        # Couleur des symboles
         addgrid.col = "gray",
         type = 'lower',diag = F,sig.level = c(0.001,0.01,0.05)
         ,tl.col = 'black',tl.srt = 45,tl.cex = 0.7
)

############################# ############################# ############################# 
############################# Obs het~ ratio ############################# ###################
############################# ############################# #############################
#calcualte habitat change by population
Pop_Habitat$delta<-c((Pop_Habitat$habitat_00))-c((Pop_Habitat$habitat_51))
Pop_Habitat$ratio<-c(Pop_Habitat$habitat_00+1)/c(Pop_Habitat$habitat_51+1)

#filter and keep only the sampling locations in Puerto Rico 
Pop_stats_PR <-  Pop_stats[which(Pop_stats$island=='PR'),]

#add habitat change to the dataframe
Pop_stats_PR <- Pop_stats_PR %>%
left_join(
  Pop_Habitat %>% select(sp_code, pop_ID , delta , ratio),
  by = c("sp_code", "pop_ID")
)

coeff_Ho<-NULL
keep<-names(which(c(unlist(by(Pop_stats_PR$obs_het,Pop_stats_PR$sp_code,function(x) length(x))))>4))
for(i in keep){
  dataf2<-Pop_stats_PR[(Pop_stats_PR$sp_code==i),]
  coeff<-summary(lm(scale(dataf2$obs_het)~scale(dataf2$ratio)))$coefficients[2,]
  coeff_Ho<-rbind.data.frame(coeff_Ho,coeff)
}

rownames(coeff_Ho)<-keep
colnames(coeff_Ho)<-names(coeff)
coeff_Ho<-as.data.frame(coeff_Ho)

coeff_Ho$sp_code <- rownames(coeff_Ho)
rownames(coeff_Ho) <- NULL
coeff_Ho<-coeff_Ho[, c("sp_code" ,"Estimate", "Std. Error", "t value",  "Pr(>|t|)")]

coeff_Ho <- coeff_Ho %>%
  left_join(
    traits %>% select(sp_code,wd.z,ss.log.z,thk.log.z, sla.log.z,maxht.z,PC1),
    by = c("sp_code")
  )

coeff_Ho$sp_name <- Pop_stats_PR$sp_name[match(coeff_Ho$sp_code,Pop_stats_PR$sp_code)]

######plot
par(mar = c(4, 15, 0.5, 0.5))

# Order alphabetically by species name
ord <- order(coeff_Ho$sp_name)
coeff_Ho_ord <- coeff_Ho[ord, ]
sp_names_ord <- (coeff_Ho$sp_name[ord])

#oder by coefficient
ord <- order(coeff_Ho$Estimate)
coeff_Ho_ord <- coeff_Ho[ord, ]
sp_names_ord <- coeff_Ho$sp_name[ord]


## ---- Build italic Genus species labels ---------------------------------

labs <- as.expression(
  lapply(sp_names_ord, function(x) {
    parts <- strsplit(x, " ")[[1]]
    bquote(italic(.(parts[1]) ~ .(parts[2])))
  })
)

## ---- Main Plot ----------------------------------------------------------

plot(
  (coeff_Ho_ord$Estimate),
  1:nrow(coeff_Ho_ord),
  cex.lab = 1.5,
  ylab = "",
  xlim = c(-2, 2),
  pch = 19,
  col = "gray50",
  cex = 2,
  yaxt = "n",
  cex.axis = 1.5,
  xlab = expression(italic(H)[o] ~ "~ Ratio")
)

abline(v = 0, lty = 4, col = "gray65")


## ---- Add Confidence Intervals ------------------------------------------

for (ci in 1:3) {
  
  alpha <- c(0.05, 0.10, 0.20)[ci]          # significance levels
  z_value <- qnorm(1 - alpha / 2)
  
  lower_bound <- coeff_Ho_ord$Estimate - z_value * coeff_Ho_ord$`Std. Error`
  upper_bound <- coeff_Ho_ord$Estimate + z_value * coeff_Ho_ord$`Std. Error`
  
  segments(
    lower_bound, 1:nrow(coeff_Ho_ord),
    upper_bound, 1:nrow(coeff_Ho_ord),
    col = "gray65",
    lwd = c(2, 4, 5)[ci]
  )
}

legend("bottomright",
       legend = c("95% CI", "90% CI", "80% CI"),
       lwd = c(1, 3, 5),
       bty = "n",
       cex = 1)


## ---- Add Y-axis with alphabetical italic labels ------------------------

axis(
  2,
  at = 1:length(labs),
  labels = labs,
  las = 1,
  cex.axis = 1.3
)

#####---------made 

plot(coeff_Ho_ord$Estimate,coeff_Ho_ord$PC1,
     pch = 19,
     xlab = expression(Ratio["00/51"]~"~"~Ho),
     ylab = "PC1", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     cex=1.5, 
     col = "gray40")
abline(lm(coeff_Ho_ord$PC1 ~ coeff_Ho_ord$Estimate), lwd = 2.5)
summary(lm(coeff_Ho_ord$PC1 ~ coeff_Ho_ord$Estimate))
##new results
# coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)         0.6226     0.4774   1.304    0.225
# coeff_Ho$Estimate   0.2534     1.1067   0.229    0.824

#old results 
# Estimate Std. Error t value Pr(>|t|)  
# (Intercept)            0.03457    0.35781   0.097   0.9252  
# coeff_Ho_ord$Estimate -2.83629    0.98556  -2.878   0.0182 *
# 
# cols <- colorRampPalette(c("darkgreen", "lightgreen"))(10)


# Map PC1 values to colors
col_index <- cut(coeff_Ho_ord$PC1, breaks = 10, labels = FALSE)

par(mar = c(5, 5, 1, 1))
plot(coeff_Ho_ord$Estimate,coeff_Ho_ord$PC1,
     col = cols[col_index],          # change color of points
     pch = 19,
     cex=2,# filled circles
     ylab = "",   # custom x‑axis label
     xlab = "", 
     xaxt = "n",
     cex.axis = 1,
     cex.lab =1
)
abline(lm(coeff_Ho_ord$Estimate ~ coeff_Ho_ord$PC1 ), lwd = 3)
abline(h = 0, lwd = 2, lty = 2, col='gray')

######correlation plot figure S6 --> I have to finish this! 
colnames(dataf2)
coeff3b<-NULL
for(j in c(18:21)){
  
  coeff3<-NULL
  keep<-names(which(c(unlist(by(Pop_stats_PR$obs_het,Pop_stats_PR$sp_code,function(x) length(x))))>4))
  for(i in keep){
    #dataf2<-dataf[which(dataf$sps==i),]
    coeff<-summary(lm(scale(dataf2$obs_het)~scale(dataf2[,j])))$coefficients[2,]
    coeff3<-rbind.data.frame(coeff3,coeff)
  }
  coeff3b <- cbind(coeff3b,coeff3[,1])
}


head(coeff3b)
rownames(coeff3b)<-keep
colnames(coeff3b)<-colnames(dataf2)[18:21]
coeff3b<-as.data.frame(coeff3b)

specs<-read.table('Fst_traits_across_island.txt',header = T)
head(specs)
rownames(specs)<-specs$Species

out<-merge(coeff3b,specs,by = 0 )
head(out)

require(corrplot)
matrice_spearman <- cor(out[,c(2:5,17,20,22,23,24,30)], method = "pearson", use = "complete.obs")
msub<-matrice_spearman[c(5:nrow(matrice_spearman)),c(1:4)]
msub<-msub[c(6,1:5),]
res <- cor.mtest(out[,c(2:5,17,20,22,23,24,30)], method = "pearson",conf.level = 0.95)
rsub<-res$p[c(5:nrow(matrice_spearman)),c(1:4)]
rsub<-rsub[c(6,1:5),]
#signes <- ifelse(matrice_spearman > 0, "+", "-")
corrplot(msub,
         p.mat = rsub,
         insig = 'label_sig',  # Cache les non-significatifs
         method = "circle",
         # Fond coloré pour le signe/force             # Utilise nos symboles + et -
         pch.cex = .7,            # Taille des symboles
         pch.col = "gray0",        # Couleur des symboles
         addgrid.col = "gray",
         diag = T,sig.level = c(0.001,0.01,0.05)
         ,tl.col = 'black',tl.srt = 45,tl.cex = 0.7
)

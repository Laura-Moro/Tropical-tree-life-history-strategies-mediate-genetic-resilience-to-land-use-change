######################################################
### Run PR-ONLY MPLE models (pairwise IBD/IBE) ###
######################################################

##########################################
### Load libraries
library(lme4) # For fitting MPLE models
library(MuMIn) # For R2 of mixed models
library(cluster) # For getting distance matrices
library(sf) # For vector spatial data
library(terra) # For raster spatial data
library(raster) # For raster spatial data

##########################################
### Prep FST and other data
dat <- readRDS("Data/LG/fst_PR.rds")

### Read and add forest proportion data
habprop <- read.csv("Data/Sp-specific-pairwise_forest_proporotion.csv")

dat <- cbind(dat, habprop[match(paste(dat$species, dat$pair), 
                                           paste(habprop$sp, habprop$pair)),
                                     c("geo", "cosDisth51", "resDisth51", 
                                       "cosDisth00", "resDisth00", 
                                       "cosDisth51t", "resDisth51t",
                                       "cosDisth00t", "resDisth00t", 
                                       "cosDistf51", "resDistf51", 
                                       "cosDistf00", "resDistf00",
                                       "cosDistsdm", "resDistsdm")])

### Get ENV site data
env <- read.csv("Data/LG/smapling_rda.csv", row.names=1)
env <- env[order(env$Location),]
env$GEO <- as.factor(env$GEO)

### Read site data
sites <- read.csv("Data/LG/eco_samp.csv", row.names = 1)
sites <- sites[order(sites$Location),]
pts_df <- sites[,1:3]
names(pts_df)[1] <- c("name")
pts_sf <- st_as_sf(pts_df, coords = c("x", "y"), crs = 4326)
pts_sf <- st_transform(pts_sf, 9001)
geo_distmat <- st_distance(pts_sf)
df_geo <- data.frame(site1=rep(env$Location, each=nrow(geo_distmat)),
                     site2=rep(env$Location, times=nrow(geo_distmat)),
                     envdist=as.vector(geo_distmat))
df_geo <- df_geo[df_geo$site1!=df_geo$site2,]
df_geo <- df_geo[!duplicated(df_geo$envdist),]
df_geo$pair <- paste(df_geo$site1, df_geo$site2, sep=":")
dat$geo_dist <- df_geo$envdist[match(dat$pair, df_geo$pair)]

### Environmental distance as Gower distance
env_distmat <- as.matrix(daisy(env[,c("Pm","PSI","TM","GEO")], metric="gower"))
clim_distmat <- as.matrix(daisy(env[,c("Pm","PSI","TM")], metric="euclidean"))

soilcols <- c("Alluvional", "limestone", "vulcanic", "Ultramafic")
for(i in seq_along(soilcols)){
  env[,soilcols[i]] <- as.factor(env[,soilcols[i]])
}
soil_distmat <- as.matrix(daisy(env[,soilcols], metric="gower"))

### Merge Env distance matrix with data.frame
df_env <- data.frame(site1=rep(env$Location, each=nrow(env_distmat)),
                     site2=rep(env$Location, times=nrow(env_distmat)),
                     envdist=as.vector(env_distmat))
df_env <- df_env[df_env$site1!=df_env$site2,]
df_env <- df_env[!duplicated(df_env$envdist),]
df_env$pair <- paste(df_env$site1, df_env$site2, sep=":")
dat$env_dist <- df_env$envdist[match(dat$pair, df_env$pair)]

df_clim <- data.frame(site1=rep(env$Location, each=nrow(clim_distmat)),
                     site2=rep(env$Location, times=nrow(clim_distmat)),
                     envdist=as.vector(clim_distmat))
df_clim <- df_clim[df_clim$site1!=df_clim$site2,]
df_clim$pair <- paste(df_clim$site1, df_clim$site2, sep=":")
dat$clim_dist <- df_clim$envdist[match(dat$pair, df_clim$pair)]

df_soil <- data.frame(site1=rep(env$Location, each=nrow(soil_distmat)),
                      site2=rep(env$Location, times=nrow(soil_distmat)),
                      envdist=as.vector(soil_distmat))
df_soil <- df_soil[df_soil$site1!=df_soil$site2,]
df_soil$pair <- paste(df_soil$site1, df_soil$site2, sep=":")
dat$soil_dist <- df_soil$envdist[match(dat$pair, df_soil$pair)]


### Scale variables
dat$geo_dist.z <- as.numeric(scale(dat$geo_dist))
dat$geo_dist_log.z <- as.numeric(scale(log(dat$geo_dist)))
dat$env_dist.z <- as.numeric(scale(dat$env_dist))
dat$clim_dist.z <- as.numeric(scale(dat$clim_dist))
dat$soil_dist.z <- as.numeric(scale(dat$soil_dist))

z <- apply(habprop[match(paste(dat$species, dat$pair), 
                         paste(habprop$sp, habprop$pair)),
                   c("cosDisth51", "resDisth51", 
                     "cosDisth00", "resDisth00", 
                     "cosDisth51t", "resDisth51t",
                     "cosDisth00t", "resDisth00t", 
                     "cosDistf51", "resDistf51", 
                     "cosDistf00", "resDistf00",
                     "cosDistsdm", "resDistsdm")], 2, scale)
colnames(z) <- paste0(colnames(z), ".z")
dat <- cbind(dat, z)

### Linearlize FST
# dat$fst_lin <- dat$fst_maf0.5 / (1-dat$fst_maf0.5)
dat$fst_lin <- dat$fst_sing / (1-dat$fst_sing)

### ADD TRAIT DATA
trait_dat <- read.csv("Data/LG/Species_level.csv", row.names = 1)
trait_dat$ss.log <- log(trait_dat$SM)
trait_dat$wd <- trait_dat$WD
trait_dat$sla.log <- log(trait_dat$SLA)
trait_dat$thk.log <- log(trait_dat$LTH)
trait_dat$maxht <- trait_dat$MH

foctraits <- trait_dat[trait_dat$Sp_code %in% dat$sp, 
                       c("wd","thk.log","sla.log","maxht","ss.log", "PC1", "PC2")]
foctraits.z <- apply(foctraits, 2, scale)
rownames(foctraits.z) <- rownames(foctraits) <- trait_dat$Sp_code[trait_dat$Sp_code %in% dat$sp]
foctraits.z[is.na(foctraits.z)] <- 0
(foctraits.z <- as.data.frame(foctraits.z))

### PRINCIPLE COMPONENTS ANALYSIS OF TRAITS
dat$pc1 <- foctraits.z$PC1[match(dat$species, rownames(foctraits.z))]
dat$pc2 <- foctraits.z$PC2[match(dat$species, rownames(foctraits.z))]
dat$ss <- foctraits.z$ss.log[match(dat$species, rownames(foctraits.z))]
dat$wd <- foctraits.z$wd[match(dat$species, rownames(foctraits.z))]
dat$sla <- foctraits.z$sla.log[match(dat$species, rownames(foctraits.z))]
dat$maxht <- foctraits.z$maxht[match(dat$species, rownames(foctraits.z))]

### FIT MODEL
# Select trait
dat$trait <- dat$pc1

###########################
### Community-level models

# Remove species with fewer than Nremove pairwise comparisons
table(dat$species)

Nremove <- 1 # 1 for keeping all species
dat2 <- dat[dat$species %in% names(table(dat$species))[table(dat$species) >= Nremove],]

table(dat2$species)


### First, we test if different resistance surfaces explain pairwise genetic differences among populations better than geographic distance along. We test resistance (connectivity) based on cost surfaces generated by suitable climate/soil conditions, forest cover in 1951, and a combined metric of habitat suitability + forest cover.

### GEO DIST ONLY MODEL
model_A <- lmerTest::lmer(fst_lin ~ geo_dist_log.z +
                  (1 | species) +
                  (1 | pop1) + (1 | pop2),
                REML = FALSE, data = dat2)

### HABITAT SUITABILITY ONLY
model_Bc <- lmerTest::lmer(fst_lin ~ geo_dist_log.z + cosDistsdm.z +
                            (1 | species) +
                            (1 | pop1) + (1 | pop2),
                          REML = FALSE, data = dat2)

### FOREST COVER ONLY
model_Cc <- lmerTest::lmer(fst_lin ~ geo_dist_log.z + cosDistf51.z +
                            (1 | species) +
                            (1 | pop1) + (1 | pop2),
                          REML = FALSE, data = dat2)


### HABITAT SUITABILITY + FOREST COVER
model_Dc <- lmerTest::lmer(fst_lin ~ geo_dist_log.z + cosDisth51.z +
                             (1 | species) +
                             (1 | pop1) + (1 | pop2),
                           REML = FALSE, data = dat2)

summary(model_A) # ns
summary(model_Bc) # sig
summary(model_Cc) # sig
summary(model_Dc) # sig

x <- anova(model_A, 
      model_Bc,
      model_Cc,
      model_Dc)

x$deltaAIC <- x$AIC-min(x$AIC)
x

round(MuMIn::r.squaredGLMM(model_A), 2)
round(MuMIn::r.squaredGLMM(model_Dc), 2)

### Based on AIC, including Least Cost Distance (LCD, based on a combination of forest cover in 1951 and habitat suitability) explains pairwise population genetic distance better than geographic distance alone. We do get only slightly worse results if LCD is based on only forest cover and again slightly worse if it is only based on habitat suitability. Notably, the geographic distance term does have a statistically significant effect in the geo-only model and, as a fixed effect, it explains only 0.8% of the variation (26% is explained by random effects). When we include the LCD, that term is significantly positively related to genetic distance, and the fixed (random) terms of the model explain 2.5% (31%) of the variation in pairwise genetic differences.

### Next, we see if adding pairwise climate and/or soil differences at each sampling location further improves the model fit. These models test for evidence of e.g. local adaptation / genetic differentiation is associated with pairwise differences between sites as opposed to the connectivity of sites. We build on the best selected model from above and test inclusion of climate-only, soil-only, or a distance metric including both climate and soil differences.

model_Dc_clim <- lmerTest::lmer(fst_lin ~ geo_dist_log.z + cosDisth51.z + clim_dist.z +
                             (1 | species) +
                             (1 | pop1) + (1 | pop2),
                           REML = FALSE, data = dat2)

model_Dc_soil <- lmerTest::lmer(fst_lin ~ geo_dist_log.z + cosDisth51.z + soil_dist.z  +
                                (1 | species) +
                                  (1 | pop1) + (1 | pop2),
                                REML = FALSE, data = dat2)

model_Dc_env <- lmerTest::lmer(fst_lin ~ geo_dist_log.z + cosDisth51.z + env_dist.z  +
                                (1 | species) +
                                  (1 | pop1) + (1 | pop2),
                                REML = FALSE, data = dat2)


xx <- anova(model_Dc, 
            model_Dc_clim, 
            model_Dc_soil, 
            model_Dc_env)

xx$deltaAIC <- xx$AIC - min(xx$AIC)
xx

summary(model_Dc_env)

round(MuMIn::r.squaredGLMM(model_Dc_env), 2)

anova(model_Dc_env, 
      model_Dc, 
      model_Dc_env)


### Including pairwise environmental differences (considering both climate and soil) further improves the model. Pairs of sites with larger environmental differences tend to have larger genetic differences. Explanatory power of the fixed terms raises to 21% and the full model to 56%.

## 

###########################
### Species-specific models

N_min_pairs <- 7

(focspp <- names(table(dat$species))[table(dat$species) >= N_min_pairs])

species_results <- list()

for(sp in seq_along(focspp)){
  
  focdat <- dat[dat$species == focspp[sp],]
  
  message(paste(focspp[sp], "N =", nrow(focdat)))
  
  # Prep data
  focdat$geo_dist_log.z <- as.vector(scale(log(focdat$geo_dist)))
  focdat$env_dist.z <- as.vector(scale(focdat$env_dist))
  cosDisth51.z <- as.vector(scale(focdat$cosDisth51))
  cosDistf51.z <- as.vector(scale(focdat$cosDistf51))
  cosDistsdm.z <- as.vector(scale(focdat$cosDistsdm))
  
  # Fit models
  model_A <- lmerTest::lmer(fst_lin ~ geo_dist_log.z +
                              (1 | pop1) + (1 | pop2),
                            REML = FALSE, data = focdat)
  
  model_B <- lmerTest::lmer(fst_lin ~ geo_dist_log.z + cosDistsdm.z +
                              (1 | pop1) + (1 | pop2),
                            REML = FALSE, data = focdat)
  
  model_C <- lmerTest::lmer(fst_lin ~ geo_dist_log.z + cosDistf51.z +
                              (1 | pop1) + (1 | pop2),
                            REML = FALSE, data = focdat)
  
  model_D <- lmerTest::lmer(fst_lin ~ geo_dist_log.z + cosDisth51.z +
                              (1 | pop1) + (1 | pop2),
                            REML = FALSE, data = focdat)
  
  
  model_D2 <- lmerTest::lmer(fst_lin ~ geo_dist_log.z + cosDisth51.z + env_dist.z +
                              (1 | pop1) + (1 | pop2),
                            REML = FALSE, data = focdat)
  
  species_results[[sp]] <- list(anova=anova(model_A, model_B, model_C, model_D), 
                                model_A=model_A, 
                                model_B=model_B, 
                                model_C=model_C,
                                model_D=model_D,
                                model_D2=model_D2)
  
  names(species_results)[sp] <- focspp[sp]
  
}

### Look at summary of individual models
lapply(species_results, function(x) x$anova)


modelD_res <- as.data.frame(do.call(rbind, lapply(species_results, function(x) fixef(x$model_D))))

plot(modelD_res$geo_dist_log.z,
         dat$pc1[match(rownames(modelD_res), dat$species)])
cor.test(modelD_res$geo_dist_log.z,
         dat$pc1[match(rownames(modelD_res), dat$species)])

plot(modelD_res$cosDisth51.z,
     dat$pc1[match(rownames(modelD_res), dat$species)])
cor.test(modelD_res$cosDisth51.z,
         dat$pc1[match(rownames(modelD_res), dat$species)])

modelD_res <- as.data.frame(do.call(rbind, lapply(species_results, function(x) fixef(x$model_D2))))

plot(modelD_res$env_dist.z,
     dat$pc1[match(rownames(modelD_res), dat$species)])
cor.test(modelD_res$env_dist.z,
     dat$pc1[match(rownames(modelD_res), dat$species)])


# lapply(species_results, function(x) summary(x$model_B))
lapply(species_results, function(x) summary(x$model_C))

### Look at LRTs for each species
lapply(species_results, function(x) anova(x$model_A, x$model_B))
lapply(species_results, function(x) anova(x$model_B, x$model_C))

rbind(sapply(species_results, function(x) round(anova(x$model_A, 
                                                      x$model_B)$`Pr(>Chisq)`[2], 3)),
      sapply(species_results, function(x) round(anova(x$model_B, 
                                                      x$model_C)$`Pr(>Chisq)`[2], 3)))

round(t(sapply(species_results, function(x) r.squaredGLMM(x$model_A))), 2)
round(t(sapply(species_results, function(x) r.squaredGLMM(x$model_B))), 2)
round(t(sapply(species_results, function(x) r.squaredGLMM(x$model_C))), 2)



###########################
### Trait-based models

# Read trait data

trait_dat <- read.csv("Data/Traits/Traits_fix.csv", row.names = 1)

dat2$ss.z <- as.numeric(scale(log(trait_dat$seed_mass[match(dat2$species, 
                                                            trait_dat$Code)])))
dat2$ss.z[is.na(dat2$ss.z)] <- 0

dat2$wd.z <- as.numeric(scale(trait_dat$WD[match(dat2$species, trait_dat$Code)]))
dat2$wd.z[is.na(dat2$wd.z)] <- 0


model_A <- lmer(fst_lin ~ geo_dist_log.z +
                  (1 | species) +
                  (1 | pop1) + (1 | pop2),
                REML = FALSE, data = dat2)

model_Ab <- lmer(fst_lin ~ geo_dist_log.z * ss.z +
                  (1 | species) +
                  (1 | pop1) + (1 | pop2),
                REML = FALSE, data = dat2)

# model_B <- lmer(fst_lin ~ geo_dist_log.z + env_dist.z +
#                   (1 | species) +
#                   (1 | pop1) + (1 | pop2),
#                 REML = FALSE, data = dat2)
# 
# model_C <- lmer(fst_lin ~ geo_dist_log.z + env_dist.z + forest_prop51.z +
#                   (1 | species) +
#                   (1 | pop1) + (1 | pop2),
#                 REML = FALSE, data = dat2)

anova(model_A, model_Ab)


r.squaredGLMM(model_A)
r.squaredGLMM(model_B)
r.squaredGLMM(model_C)


anova(model_A, model_B)
anova(model_B, model_C)







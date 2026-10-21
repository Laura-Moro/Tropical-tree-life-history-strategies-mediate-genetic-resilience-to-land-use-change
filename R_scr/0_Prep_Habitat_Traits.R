#prepare data 
#1# extract trait data and habitat data for each species
#2# caluclate habitata ammount at the sampling location 

library(terra)
library(dplyr)
library(tidyr)

#load data
Sp_H_T <- read.csv("Data/FIGSHARE/6b-output-20260830.csv") # habitata ammount data form figshare
samp_locations <- read.csv("Data/Sampling_location_Eco.csv") #all smapling location
pops <- read.csv("Data/PopGen_stats_pops.csv") #sampling loacations for each species 
traits <-  read.csv("Data/Original_trait_values.csv")

#match the habitat data to the trait data 
sp <- unique(pops$sp_code)
Sp_H_T_filt <- Sp_H_T[Sp_H_T$sp %in% sp, ]
habitat <- Sp_H_T_filt[,c("sp","species","fcover_51_raw","fcover_00_raw")]

df <- traits

df$SH_51 <- habitat$fcover_51_raw[match(df$sp_code,habitat$sp)]
df$SH_00 <- habitat$fcover_00_raw[match(df$sp_code,habitat$sp)]

write.csv(df,"Data/Species_traits_habitat.csv", row.names = FALSE)

# #### !!!! rerun the pca with the original trit values
# df_PCA <- df[,c("sp_code","ss.log.z","thk.log.z","sla.log.z","maxht.z","wd.z")]
# df_PCA[is.na(df_PCA)] <- 0
# pca <- prcomp(df_PCA[, -1], center = T, scale. = T)
# pca_scores <- cbind(
#   sp = df_PCA$sp,
#   as.data.frame(pca$x)
# )

  # ## !!!!!! This is is we re run the PCA from the published data!!!!
# #1# extract trait data and habitat data for each species
# df <- Sp_H_T_filt[,c("sp","species","fcover_51_raw","fcover_00_raw",
#   "ss.log.z","thk.log.z","sla.log.z","maxht.z","wd.z","pca1")]
# 
# df_PCA <- df[,c("sp","ss.log.z","thk.log.z","sla.log.z","maxht.z","wd.z")]
# df_PCA[is.na(df_PCA)] <- 0
# pca <- prcomp(df_PCA[, -1], center = F, scale. = FALSE)
# pca_scores <- cbind(
#   sp = df_PCA$sp,
#   as.data.frame(pca$x)
# )
# 
# df$PC1 <- pca_scores$PC1[match(df$sp,pca_scores$sp)]
# 
# df <- df %>%
#   rename(
#     sp_name = species,
#     sp_code = sp,
#     SH_51 =fcover_51_raw,
#     SH_00=fcover_00_raw
#   )
# 
# write.csv(df,"Data/Species_traits_habitat.csv", row.names = FALSE)


#2# caluclate habitata ammount at the sampling location 

# Sampling points
samp_locations <- read.csv("Data/Sampling_location_Eco.csv")
pops <- read.csv("Data/PopGen_stats_pops.csv")

pops_PR <- pops %>%
  filter(island == "PR") %>%
  select(sp_name, sp_code, pop_ID) %>%
  arrange(sp_code)

pops_PR$lat <- samp_locations$lat[
  match(pops_PR$pop_ID, samp_locations$Location)
]

pops_PR$lon <- samp_locations$lon[
  match(pops_PR$pop_ID, samp_locations$Location)
]

sp_keep <- unique(pops_PR$sp_code)

sp_keep <- intersect(sp_keep, names(habitat_raw_51))

habitat_raw_51 <- rast("Data/FIGSHARE/Pred_f51_raw-near-20260830.tif")
habitat_raw_00 <- rast("Data/FIGSHARE/Pred_f00_raw-near-20260830.tif")

# subset raster
habitat_raw_51 <- habitat_raw_51[[sp_keep]]
habitat_raw_00 <- habitat_raw_00[[sp_keep]]

results <- data.frame()

for (i in 1:nlyr(habitat_raw_51)) {
  
  sp <- names(habitat_raw_51[[i]])
  
  pts_df <- pops_PR %>%
    filter(sp_code == sp)
  
  if (nrow(pts_df) == 0) next
  
  pts <- vect(
    pts_df[, c("lon", "lat")],
    geom = c("lon", "lat"),
    crs = "EPSG:4326"
  )
  
  pts_proj <- project(pts, crs(habitat_raw_51[[i]]))
  
  buf <- buffer(pts_proj, width = 1000)
  
  vals_51 <- terra::extract(
    habitat_raw_51[[i]],
    buf,
    fun = sum,
    na.rm = TRUE
  )
  
  vals_00 <- terra::extract(
    habitat_raw_00[[i]],
    buf,
    fun = sum,
    na.rm = TRUE
  )
  
  tmp <- data.frame(
    sp_code = sp,
    pop_ID = pts_df$pop_ID,
    habitat_51 = vals_51[, 2],
    habitat_00 = vals_00[, 2]
  )
  
  results <- bind_rows(results, tmp)
}

results

write.csv(results,"Data/Habitat_51_00_pop.csv", row.names = FALSE)




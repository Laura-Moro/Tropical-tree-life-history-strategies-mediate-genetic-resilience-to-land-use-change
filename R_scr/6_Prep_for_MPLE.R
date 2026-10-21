###################################################
### Compute proportion of forest between points ###
###################################################

library(terra)
library(sf)


################################################
### Load forest cover rasters
f <- rast("Data/Derived/forest_change_resampled_modal-20251022.tif")
f51 <- f[[1]]
f00 <- f[[2]]


################################################
### Load SDM output habitat for genetic species

# Get genetic species names
filenames <- list.files("Data/LG/fst_sing")
genetic_sp <- sapply(filenames[grep("_fst_hud.txt", filenames)], function(x) substring(x, 1, 6))
combined_genetic_sp <- paste(genetic_sp, collapse = "|")

h51 <- rast("Data/Derived/Pred_f51_raw-near-20260830.tif")
h00 <- rast("Data/Derived/Pred_f00_raw-near-20260830.tif")
h51t <- rast("Data/Derived/Pred_f51-near-20260830.tif")
h00t <- rast("Data/Derived/Pred_f00-near-20260830.tif")

h51 <- h51[[names(h51) %in% genetic_sp]]
h00 <- h00[[names(h00) %in% genetic_sp]]
h51t <- h51t[[names(h51t) %in% genetic_sp]]
h00t <- h00t[[names(h00t) %in% genetic_sp]]

sdm <- rast("Data/Derived/raw_stack-20260830.tif")
sdm <- sdm[[names(sdm) %in% names(h51)]]


################################################
### Get coordiantes of sample locations and convert for downstream use
sites <- read.csv("Data/LG/eco_samp.csv", row.names = 1)

# Slight adjustment of a point because it falls off the shore with aggregated raster...
sites$y[sites$Location=="LP"] <- 17.958

sites <- sites[order(sites$Location),]
pts_df <- sites[,1:3]
(names(pts_df)[1] <- c("name"))
pts_sf <- st_as_sf(pts_df, coords = c("x", "y"), crs = 4326)
pts_sf <- st_transform(pts_sf, st_crs(f51))


################################################
### Extract coordinates and names, get pair index
coords <- st_coordinates(pts_sf)
pt_names <- pts_sf$name
n <- nrow(coords)
pair_index <- combn(n, 2)
n_pairs <- ncol(pair_index)


################################################
### Create pairwise lines
line_list <- vector("list", n_pairs)
pt1_names <- character(n_pairs)
pt2_names <- character(n_pairs)

for (i in seq_len(n_pairs)) {
  idx1 <- pair_index[1, i]
  idx2 <- pair_index[2, i]
  p1 <- coords[idx1, ]
  p2 <- coords[idx2, ]
  line_list[[i]] <- st_linestring(rbind(p1, p2))
  pt1_names[i] <- pt_names[idx1]
  pt2_names[i] <- pt_names[idx2]
}

lines <- st_sf(
  pair_id = seq_len(n_pairs),
  pt1 = pt1_names,
  pt2 = pt2_names,
  geometry = st_sfc(line_list, crs = st_crs(pts_sf))
)

lines$total_length <- st_length(lines)


# convert sf lines to SpatVector
lines_vect <- vect(lines)


################################################
### Generate resistance surface

library(gdistance)

sp_dist <- list()

for(sp in seq_along(names(h51))){
  
  message(paste(sp, names(h51)[sp]))
  
  focrast51 <- raster(h51[[sp]] + 1)
  focrast00 <- raster(h00[[sp]] + 1)
  focrast51t <- raster(h51t[[sp]] + 1)
  focrast00t <- raster(h00t[[sp]] + 1)
  focrastf51 <- raster(f51 + 1)
  focrastf00 <- raster(f00 + 1)
  focrastsdm <- raster(sdm + 1)

  values(focrast51)[is.na(values(focrast51))] <- 0
  values(focrast00)[is.na(values(focrast00))] <- 0
  values(focrast51t)[is.na(values(focrast51t))] <- 0
  values(focrast00t)[is.na(values(focrast00t))] <- 0
  values(focrastf51)[is.na(values(focrastf51))] <- 0
  values(focrastf00)[is.na(values(focrastf00))] <- 0
  values(focrastsdm)[is.na(values(focrastsdm))] <- 0
  
  trh51 <- transition(focrast51, mean, directions = 8)
  trh00 <- transition(focrast00, mean, directions = 8)
  trh51t <- transition(focrast51t, mean, directions = 8)
  trh00t <- transition(focrast00t, mean, directions = 8)
  tr51 <- transition(focrastf51, mean, directions = 8)
  tr00 <- transition(focrastf00, mean, directions = 8)
  trsdm <- transition(focrastsdm, mean, directions = 8)
  
  trCtrh51 <- geoCorrection(trh51, "c", scl = TRUE)
  trRtrh51 <- geoCorrection(trh51, "r", scl = TRUE)
  
  trCtrh00 <- geoCorrection(trh00, "c", scl = TRUE)
  trRtrh00 <- geoCorrection(trh00, "r", scl = TRUE)
  
  trCtrh51t <- geoCorrection(trh51t, "c", scl = TRUE)
  trRtrh51t <- geoCorrection(trh51t, "r", scl = TRUE)
  
  trCtrh00t <- geoCorrection(trh00t, "c", scl = TRUE)
  trRtrh00t <- geoCorrection(trh00t, "r", scl = TRUE)
  
  trCtrf51 <- geoCorrection(tr51, "c", scl = TRUE)
  trRtrf51 <- geoCorrection(tr51, "r", scl = TRUE)
  
  trCtrf00 <- geoCorrection(tr00, "c", scl = TRUE)
  trRtrf00 <- geoCorrection(tr00, "r", scl = TRUE)
  
  trCsdm <- geoCorrection(trsdm, "c", scl = TRUE)
  trRsdm <- geoCorrection(trsdm, "r", scl = TRUE)
  
# ---
  cosDisth51 <- costDistance(trCtrh51, st_coordinates(pts_sf))
  resDisth51 <- commuteDistance(trRtrh51, st_coordinates(pts_sf))

  cosDisth00 <- costDistance(trCtrh00, st_coordinates(pts_sf))
  resDisth00 <- commuteDistance(trRtrh00, st_coordinates(pts_sf))
  
  cosDisth51t <- costDistance(trCtrh51t, st_coordinates(pts_sf))
  resDisth51t <- commuteDistance(trRtrh51t, st_coordinates(pts_sf))
  
  cosDisth00t <- costDistance(trCtrh00t, st_coordinates(pts_sf))
  resDisth00t <- commuteDistance(trRtrh00t, st_coordinates(pts_sf))
  
  cosDistf51 <- costDistance(trCtrf51, st_coordinates(pts_sf))
  resDistf51 <- commuteDistance(trRtrf51, st_coordinates(pts_sf))
  
  cosDistf00 <- costDistance(trCtrf00, st_coordinates(pts_sf))
  resDistf00 <- commuteDistance(trRtrf00, st_coordinates(pts_sf))
  
  cosDistsdm <- costDistance(trCsdm, st_coordinates(pts_sf))
  resDistsdm <- commuteDistance(trRsdm, st_coordinates(pts_sf))
  
  sp_dist[[sp]] <- data.frame(sp = names(h51)[sp],
                              pt1 = lines$pt1,
                              pt2 = lines$pt2,
                              geo = as.vector(st_distance(pts_sf)[lower.tri(st_distance(pts_sf))]),
                              cosDisth51 = as.vector(cosDisth51),
                              resDisth51 = as.vector(resDisth51),
                              cosDisth00 = as.vector(cosDisth00),
                              resDisth00 = as.vector(resDisth00),
                              cosDisth51t = as.vector(cosDisth51t),
                              resDisth51t = as.vector(resDisth51t),
                              cosDisth00t = as.vector(cosDisth00t),
                              resDisth00t = as.vector(resDisth00t),
                              cosDistf51 = as.vector(cosDistf51),
                              resDistf51 = as.vector(resDistf51),
                              cosDistf00 = as.vector(cosDistf00),
                              resDistf00 = as.vector(resDistf00),
                              cosDistsdm = as.vector(cosDistsdm),
                              resDistsdm = as.vector(resDistsdm))
}
    

spdf <- do.call(rbind, sp_dist)

head(spdf)

corrplot::corrplot(cor(spdf[,-c(1:3)]), type='lower', addCoef.col=T, diag=F)

spdf$pair <- paste(spdf$pt1, spdf$pt2, sep=":")

write.csv(spdf, "Data/Derived/Sp-specific-pairwise_forest_proporotion-20261005.csv", row.names=F)


#data paper plots 
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(patchwork)
library(ape)
library(tidyverse)
library(terra)
library(viridis)
library(sf)
library(ggrepel)
library(ggtree)
library(treeio)
library(ggplot2)


files <- c("AMYELE", "CASARB", "CASSYL", "DACEXC", "GYMLUC", "KRUFER", 
           "POIFLO", "SCHMOR", "TABHET", "BURSIM", "CASGUI", "FAROCC",
           "HIRRUG", "MANBID", "MYRDEF", "PREMON", "SLOBER", "MICGUY") #,"CECSCH")

files <- ("CECSCH")

taxonomy <- read.csv("/Users/laumo791/Documents/PR/Data_paper/Data/species_info/Taxonomy.csv")
tax <- taxonomy[,c("Sp_code", "Genus", "Species_name")]
tax$species <- paste(tax$Genus, tax$Species_name)

samp <- read.csv("/Users/laumo791/Documents/PR/Data_paper/Data/species_info/Sampling_locations.csv")
# Keep PR rows unchanged
PR_rows <- samp[samp$island == "PR", ]

# Collapse DR and VI into one row each
DR_VI_means <- aggregate(
  cbind(Median_Latitude, Median_Longitude) ~ island,
  data = samp[samp$island %in% c("DR", "VI"), ],
  FUN = mean
)

# Set Location equal to island name
DR_VI_means$Location <- DR_VI_means$island

# Reorder columns to match original data frame
DR_VI_means <- DR_VI_means[, c(
  "Location",
  "Median_Latitude",
  "Median_Longitude",
  "island"
)]

# Combine everything
samp_final <- rbind(PR_rows, DR_VI_means)
samp_final <- samp_final[order(samp_final$Median_Longitude), ]
pop_order <- samp_final$Location

for (f in seq_along(files)) {
  
  setwd(paste0("/Users/laumo791/upp_new/nobackup/private/Laura/vcf_pascal/species/", files[f] ,"/admixture"))
  
  palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  sp <- dir(pattern = "\\cv.txt$")
  
  cv <- read.table(sp, header = T)
  
  name <-sub("_cv.txt","",sp)
  
  #find the best k 
  best_K<- cv$K[which.min(cv$CV_error)]
  if(best_K==1){
    K_plus_1 <- best_K + 1
    K_minus_1 <- best_K + 2}else{
      K_plus_1 <- best_K + 1  
      K_minus_1 <- best_K - 1}
  
  #list the Q files 
  list <- dir(pattern = "\\.Q$")
  
  #create the name 
  file_best_K <- paste0( name,".stitched.Sing.", best_K ,".Q")
  #import the file
  Q_best_K <- read.table(file_best_K, header = F)
  
  #add the names of popomp and the samples list 
  popmap <- dir(pattern = "pop_filtered.*\\.txt$")
  
  #import pop files the files
  pop <-read.table(popmap,header = F)
  
  #add the popmap info to the Q file 
  Q_K_pop<- cbind(pop, Q_best_K)
  
  Q_K_pop$V2 <- ifelse(
    str_sub(Q_K_pop$V2, 1, 2) %in% c("DR", "VI"),
    str_sub(Q_K_pop$V2, 1, 2),
    Q_K_pop$V2
  )
  
  #raname columns best K 
  num_cols <- ncol(Q_K_pop)
  new_colnames <- c("BMK_id", "pop") # First two are fixed
  if (num_cols > 2) { 
    cluster_names <- paste0("Cluster", 1:(num_cols - 2)) # Name the remaining columns
    new_colnames <- c(new_colnames, cluster_names)
  }
  colnames(Q_K_pop) <- new_colnames
  
  # change data to long format 
  Q_K_pop_long <- tidyr::pivot_longer(
    Q_K_pop,
    cols = starts_with("Cluster"),
    names_to = "Cluster",
    values_to = "Ancestry"
  )
  
  present_pops <- intersect(pop_order, unique(Q_K_pop_long$pop))
  
  Q_K_pop_long$pop <- factor(
    Q_K_pop_long$pop,
    levels = present_pops
  )
  
  
  #assign colours
  cluster_colors <-palette[c(1,3,7)] # Orange
  
  #plot
  K <- ggplot(Q_K_pop_long, aes(x = factor(BMK_id), y = Ancestry, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    facet_grid(~ pop, scales = "free_x", space = "free_x") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          panel.spacing = unit(0.1, "lines"),
          strip.text = element_text(size = 10)) +
    labs(x = "Individuals", y = "Ancestry Proportion") +
    scale_fill_manual(
      name = "Ancestry component", 
      values = cluster_colors, 
      labels = c("component 1", "component 2", "component 3"))

  
 ##Tree
  
  # Set working directory
  setwd(paste0("/Users/laumo791/upp_new/nobackup/private/Laura/vcf_pascal/species/", files[f], "/tree"))
  
  # Load tree and population data
  pop_I <- pop %>%
    mutate(island = str_sub(V2, 1, 2))
  
  pop_I <- pop_I %>%
    mutate(island = case_when(island %in% c("DR", "VI") ~ island,TRUE ~ "PR"))
  #recode.subs10
  tree_file <- paste0(name, ".stitched.Sing.subs10.recode.min4.phy.treefile")
  
  tree <- read.tree(tree_file)
  
  # Filter and reorder population data to match tree tips
  pop_I <- pop_I[pop_I$V1 %in% tree$tip.label, ]
  pop_I <- pop_I[order(match(pop_I$V1, tree$tip.label)), ]
  
  # Merge tree tip labels with population metadata
  tree_data <- tibble(label = tree$tip.label) %>%
    left_join(pop_I, by = c("label" = "V1"))
  
  # Assign colors and shapes
  custom_colors <- palette[c(4,2,8)]
  
  custom_colors <- c(
    PR = palette[2],
    DR = palette[4],
    VI = palette[8]
  )
  
  tree_data$island <- as.factor(tree_data$island)
  tree_data$island_shape <- as.numeric(as.factor(tree_data$island)) + 16  # shapes starting from 17
  
  p <- ggtree(tree, layout = "ape") %<+% tree_data +
    geom_tippoint(aes(color = island), size = 2) +  # Optional: keeps colored dots  # Labels from tree_data$V2
    scale_color_manual(values = custom_colors, name = "Island") +
    ggtitle("") +
    theme(
      legend.position = "right",
      legend.justification = "center",
      plot.title = element_text(hjust = 0.5, size = 16)
    )
  
  #SDMS
  raster_path <- paste0("/Users/laumo791/Documents/PR/C1/PR_Habitat_Change/Results/SDM_predictions/", name, ".tif")
  shp <- st_read("/Users/laumo791/Documents/PR/C1/PR_Habitat_Change/Data/PR_shapes/outline/PR_outline_Project.shp")
  df <- read.csv("/Users/laumo791/Documents/PR/Genetics/Data/collection/sampling_PR_DR_VI.csv", sep = ",", header  = T)
  r <- rast(raster_path)
  
  samp <- setdiff(unique(pop$V2), c("DR", "VI"))
  sampling<-as.data.frame(samp)
  sampling$lat <- df$Median_Latitude[match(sampling$samp, df$Location)]
  sampling$long <- df$Median_Longitude[match(sampling$samp, df$Location)]
  
  
  layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE))
  r_df <- as.data.frame(r, xy = TRUE)
  colnames(r_df) <- c("x", "y", "value")  # Rename for clarity
  
  M <- ggplot() +
    geom_raster(data = r_df, aes(x = x, y = y, fill = value)) +
    geom_sf(data = shp, fill = NA, color = "#424242", size = 0.5) +  # shapefile outline
    scale_fill_viridis(name = "Environmental suitability", na.value = "transparent") +
    geom_point(data = sampling, aes(x = long, y = lat), color = "brown4", size = 2.5) +
    coord_sf() +
    labs(title = "") +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.4, size = 18, face = "bold"),
      legend.position = "right",
      legend.justification = "center",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key.width = unit(0.4, "cm"),
      legend.key.height = unit(0.4, "cm")
    )
  
  setwd("/Users/laumo791/Documents/PR/Paper_4/multipannel")
  
  
  sp_name <- tax$species[match(name, tax$Sp_code)]
  title_expr <- as.expression(bquote(italic(.(sp_name))))
  
  #combined_plot <- ((M | p) / K) +
  #  plot_annotation(
  #    title = title_expr,
   #   theme = theme(
   #     plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
   #   )
   # )
  
  
  combined_plot <- ((M | p) / K) +
    plot_annotation(
      title = title_expr,
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 16, face = "bold")
      ))

  
  ggsave(
    filename = paste0(name , "_multipanel.png"),
    plot = combined_plot,
    width = 10, height = 5, dpi = 300
  )
}











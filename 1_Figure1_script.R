#________________________________________________________________________________________________
#
# PCAs, karyotype classification and figure 1
#
# Vinagre-Izquierdo et al. (2025) Chromosomal inversion associated with dietary 
# niche differences in common quails sharing wintering grounds. Molecular Ecology.
#
# by Celia Vinagre-Izquierdo (Feb 2025)
#
# R version 4.4.2 (2024-10-31 ucrt) -- "Pile of Leaves"
# RStudio 2024.12.0+467 "Kousa Dogwood" Release (cf37a3e5488c937207f992226d255be71f5e3f41, 2024-12-11) for windows
# Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) RStudio/2024.12.0+467 Chrome/126.0.6478.234 Electron/31.7.6 Safari/537.36, Quarto 1.5.57
#________________________________________________________________________________________________


# PACKAGES ----

  library(readr)
  
  library(ggplot2)
  
  library(tidyverse)
  
  library(jpeg) 
  
  library(grid)
  
  library(cowplot)
  
  library(magick)
  
  library(png)
  
  library(ggpubr)
  
  library(dplyr)
  
  library(svglite)



# DATA ----
  project.path <- getwd()

  ## PCA all genome ----
  pca <- read_table(file.path(project.path, "1_Genomic_Analyses", "wintering_maf003_miss075.eigenvec"), col_names = FALSE)
  eigenval <- scan(file.path(project.path,"1_Genomic_Analyses", "wintering_maf003_miss075.eigenval"))

  # Sort out the PCA data
  pca <- pca[,-1]  # Remove nuisance column
  names(pca)[1] <- "ind"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  
  # Read population data
  pop.data <- read.table(file.path(project.path, "1_Genomic_Analyses", "snp_wintering.txt"), sep = "\t", header = TRUE)
  combined_data_check <- merge(pop.data, pca, by = "ind")
  
  # Calculate percentage variance explained
  pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
  
  # Plot percentage variance explained
  ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()
  
  # Calculate cumulative sum of percentage variance explained
  cumsum(pve$pve)
  
  # Define groups of individuals
  jap <- c("453", "1124PS", "91PS", "93PS")
  hybrids <- c("639", "1094PS", "49CE", "86PS", "87PS", "88PS")
  selected_individuals <- c("41PS", "63PS", "80PS")
  
  # Subset data for groups
  jap <- subset(combined_data_check, ind %in% jap)
  jap$ind <- "Domestic Japanese quails"
  hybrids <- subset(combined_data_check, ind %in% hybrids)
  hybrids$ind <- "Admixed quails"
  gen_1 <- subset(combined_data_check, ind %in% selected_individuals)
  
  # Define custom colors and labels
  custom_colors <- c('80PS' = '#c2a545', '41PS' = '#677e39', '63PS' = '#b27a91', "Domestic Japanese quails" = "#c5793e", "Admixed quails" = "#2baad3")
  custom_labels <- c("80PS" = "AA", "63PS" = "AB", "41PS" = "BB", "Domestic Japanese quails", "Admixed quails")
  
  # Plot PCA
  ggplot(combined_data_check, aes(PC1, PC2)) +
    geom_point(size = 3, alpha = 0.75) +
    geom_point(data = jap, aes(col = ind), size = 3) +
    geom_point(data = hybrids, aes(col = ind), size = 3) +
    geom_point(data = gen_1, aes(col = ind), size = 5) +
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
    labs(col = "Previously genotyped quails") +  
    theme_bw(base_size = 20) +
    annotate("text", x = 0.1, y = 0.3, label = "Domestic Japanese quails", color = '#c5793e', size = 5, hjust = 0) +
    annotate("text", x = 0.25, y = -0.3, label = "Admixed quails", color = '#2baad3', size = 5, hjust = 0) +
    scale_color_manual(values = custom_colors, labels = custom_labels) +
    theme(legend.position = "none", legend.title = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.text = element_text(size = 60))
  
  

  ## PCA inversion chr 1 ----
  pca2 <- read_table(file.path(project.path, "1_Genomic_Analyses", "inv1_samples_control.eigenvec"), col_names = FALSE)
  eigenval2 <- scan(file.path(project.path, "1_Genomic_Analyses", "inv1_samples_control.eigenval"))
  
  # Sort out the PCA data
  pca2 <- pca2[,-1]  # Remove nuisance column
  names(pca2)[1] <- "ind"
  names(pca2)[2:ncol(pca2)] <- paste0("PC", 1:(ncol(pca2)-1))
  
  # Calculate percentage variance explained
  pve2 <- data.frame(PC = 1:20, pve2 = eigenval2/sum(eigenval2)*100)
  
  # Combine PCA and population data
  wint_inv1 <- merge(pop.data, pca2, by = "ind")
  
  # Define custom colors and labels for inversion
  custom_colors_inv <- c('80PS' = '#c2a545', '41PS' = '#677e39', '63PS' = '#b27a91')
  custom_labels_inv <- c("80PS" = "AA", "63PS" = "AB", "41PS" = "BB")
  
  # Plot PCA inversion
  ggplot(wint_inv1, aes(PC1, PC2, col = ind)) +
    geom_point(size = 5) +
    labs(col = "Previously genotyped quails") +
    xlab(paste0("PC1 (", signif(pve2$pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve2$pve[2], 3), "%)")) +
    theme_bw(base_size = 20) +
    scale_color_manual(values = custom_colors_inv, labels = custom_labels_inv) +
    theme(legend.position = "bottom") 

  

# Classify genotypes ----

  ## AA quails ----
  filtered_aa <- subset(pca2, PC1 < -0.1)
  filtered_aa$genotype <- "AA"
  
  ggplot(pca2, aes(PC1, PC2)) + 
    geom_point(size = 5) +
    geom_point(data = filtered_aa, aes(col = genotype), size = 5, col = "#c2a545", show.legend = FALSE) + 
    coord_equal() + 
    theme_light() 

  ## AB quails ----
  filtered_ab <- subset(pca2, PC1 < 0.15 & PC1 > -0.1)
  filtered_ab$genotype <- "AB"
  
  ggplot(pca2, aes(PC1, PC2)) + 
    geom_point(size = 5) +
    geom_point(data = filtered_ab, aes(col = genotype), size = 5, col = "#b27a91", show.legend = FALSE) + 
    coord_equal() + 
    theme_light()

  ## BB quails ----
  filtered_bb <- subset(pca2, PC1 > 0.15)
  filtered_bb$genotype <- "BB"

  ggplot(pca2, aes(PC1, PC2)) + 
    geom_point(size = 5) +
    geom_point(data = filtered_bb, aes(col = genotype), size = 5, col = "#677e39", show.legend = FALSE) + 
    coord_equal() + 
    theme_light()

  
  
# Plots for figure 1 ----
  
  ## Map and picture ----
  img <- readJPEG(file.path(project.path,"1_Genomic_Analyses", "4096-2731-max_bonita.JPG"))
  img_grob <- rasterGrob(img, interpolate = TRUE)
  
  image_plot <- ggdraw() + draw_image(file.path(project.path,"1_Genomic_Analyses", "4096-2731-max_bonita.JPG"), scale = 1)
  
  img1 <- readPNG(file.path(project.path,"1_Genomic_Analyses", "mapa_invCaracoles.png"))
  img_grob1 <- rasterGrob(img, interpolate = TRUE)
  
  mapa <- ggdraw() + draw_image(file.path(project.path,"1_Genomic_Analyses", "mapa_invCaracoles.png"), scale = 1)
  
  mapaimagen <- ggarrange(mapa, image_plot, labels = c("A", "B"), common.legend = TRUE, legend = "none", heights = c(3, 0.25), font.label = list(size = 40), widths = c(1, 1), align = "h")

  
  ## PCA all genome ----
  plot1 <- ggplot(combined_data_check, aes(PC1, PC2)) +
    geom_point(size = 3, alpha = 0.75) +
    geom_point(data = jap, aes(col = ind), size = 3) +
    geom_point(data = hybrids, aes(col = ind), size = 3) +
    geom_point(data = gen_1, aes(col = ind), size = 5) +
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
    labs(col = "Previously genotyped quails") +  
    theme_bw(base_size = 40) +
    annotate("text", x = 0.1, y = 0.3, label = "Domestic Japanese quails", color = '#c5793e', size = 10, hjust = 0) +
    annotate("text", x = 0.18, y = -0.3, label = "Admixed quails", color = '#2baad3', size = 10, hjust = 0) +
    scale_color_manual(values = custom_colors, labels = custom_labels) +
    theme(legend.position = "none", legend.title = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.text = element_text(size = 40))
  plot1

  
  ## PCA inversion chromosome 1 ----
  wint_inv1_genotype <- rbind(filtered_aa, filtered_ab, filtered_bb)
  wint_inv1_genotype_select <- wint_inv1_genotype %>% select(ind, PC1, PC2, PC3, genotype)
  
  custom_colors_gen <- c('AA' = '#c2a545', 'BB' = '#677e39', 'AB' = '#b27a91')
  selected_individuals <- c("41PS", "63PS", "80PS")
  gen <- subset(wint_inv1_genotype_select, ind %in% selected_individuals)
  
  plot2 <- ggplot(wint_inv1_genotype_select, aes(PC1, PC2)) +
    stat_ellipse(data = wint_inv1_genotype_select, aes(x = PC1, y = PC2, color = genotype, fill = genotype), 
                 geom = "polygon", alpha = 0.25, type = "norm", level = 0.97) +
    geom_point(size = 3, alpha = 0.75) +
    geom_point(data = gen, aes(col = ind), size = 5) +
    xlab(paste0("PC1 (", signif(pve2$pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve2$pve[2], 3), "%)")) +
    theme_bw(base_size = 40) + 
    scale_fill_manual(values = custom_colors_gen) +
    scale_color_manual(values = custom_colors_inv) +
    annotate("text", x = -0.26, y = 0.8, label = "AA", color = '#c2a545', size = 20, hjust = 0) +
    annotate("text", x = 0.01, y = 0.8, label = "AB", color = '#b27a91', size = 20, hjust = 0) +
    annotate("text", x = 0.22, y = 0.8, label = "BB", color = '#677e39', size = 20, hjust = 0) +
    theme(legend.position = "none", legend.title = element_blank(), 
          panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    ylim(-1, 1)
  plot2

  pca_plot <- ggarrange(plot1, plot2, labels = c("C", "D"), common.legend = TRUE, legend = "none", font.label = list(size = 40))


  ## All together in a file called figure1.svg ----
  res <- 144
  svglite("figure1.svg", width = 3500/res, height = 2000/res)
  ggarrange(mapaimagen, pca_plot, common.legend = F, legend = "none", heights = c(1, 2), nrow = 2, font.label = list(size = 60), align = "h")
  dev.off()

  
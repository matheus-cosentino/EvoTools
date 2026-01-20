#!/usr/bin/env Rscript

# -----------------
# 1. SETUP
# -----------------
required_packages <- c("argparse", "ggtree", "ggplot2", "ape", "dplyr",
                       "treeio", "tidytree", "scales", "phytools", "viridis")

# Define um mirror padrão para downloads automáticos
cran_mirror <- "https://cloud.r-project.org"

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing missing package:", pkg, "\n"))
    
    if (pkg %in% c("ggtree", "treeio", "tidytree")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = cran_mirror)
      }
      BiocManager::install(pkg, ask = FALSE)
    } else {
      # AQUI ESTAVA O ERRO: Adicionado 'repos = cran_mirror'
      install.packages(pkg, repos = cran_mirror)
    }
  }
}

library(argparse)
suppressPackageStartupMessages({
  for (pkg in required_packages) {
    library(pkg, character.only = TRUE)
  }
})

# -----------------
# 2. ARGUMENT PARSING
# -----------------
parser <- ArgumentParser(description = "Plot a phylogenetic tree with depth-sorted clade coloring.")

# Inputs
parser$add_argument("-t", "--tree", type = "character", required = TRUE, help = "Path to the input tree file")
parser$add_argument("-m", "--metadata", type = "character", required = TRUE, help = "Path to the input metadata TSV file")

# --- HIGHLIGHTING & CLADES ---
parser$add_argument("--tips_file", type="character", default=NULL, help="File to prune tree (subtree).")
parser$add_argument("--highlight_file", type="character", default=NULL, help="File with IDs to show tip points.")

parser$add_argument("--clades_file", type="character", default=NULL, 
                    help="Manual TSV (tip, label, color) for clade bars.")

parser$add_argument("--auto_clade_label", type="character", default=NULL, 
                    help="Column name in metadata to AUTOMATICALLY generate colored clade bars.")

parser$add_argument("--sort_by_depth", action="store_true", default=FALSE,
                    help="If set, clade colors are assigned based on distance from root (Basal -> Distal).")

# --- VISUAL OPTIONS ---
parser$add_argument("--hide_tip_labels", action="store_true", default=FALSE, help="Hide text labels.")
parser$add_argument("-l", "--layout", type="character", default="rectangular", choices=c("rectangular", "circular", "fan"), help="Layout type")
parser$add_argument("--flip", action="store_true", default=FALSE, help="Invert tree order.")

# --- COLUMNS ---
parser$add_argument("--id_column", type = "character", default = "Acc", help = "Metadata ID column")
parser$add_argument("--species_column", type = "character", default = "Species", help = "Metadata Label column")
parser$add_argument("--color_column", type = "character", default = "Genus", help = "Metadata Color column (for tip points)")

# --- SUPPORT VALUES ---
parser$add_argument("--support_mode", type="character", default="single", 
                    choices=c("single", "dual", "none"),
                    help="how to parse support: 'single', 'dual', or 'none'.")
parser$add_argument("--support_cutoff", type="double", default=75, help="Threshold for support values.")

# --- OUTPUT ---
parser$add_argument("-o", "--output", type = "character", default = "Final_Tree.svg", help = "Output file")
parser$add_argument("-W", "--width", type = "double", default = 10, help = "Width")
parser$add_argument("-H", "--height", type = "double", default = 12, help = "Height")

args <- parser$parse_args()

# -----------------
# 3. DATA LOADING
# -----------------
cat("Loading tree file:", args$tree, "\n")
tree <- read.tree(args$tree)

# Pruning
if (!is.null(args$tips_file)) {
  cat("Pruning tree...\n")
  tips_to_keep <- readLines(args$tips_file)
  tips_in_tree <- intersect(tips_to_keep, tree$tip.label)
  if (length(tips_in_tree) == 0) stop("No matching tips found for pruning.")
  tree <- ape::keep.tip(tree, tips_in_tree)
}

cat("Rooting tree (midpoint)...\n")
tree <- phytools::midpoint.root(tree)

cat("Loading metadata:", args$metadata, "\n")
metadata <- read.csv(args$metadata, sep = "\t", header = TRUE, comment.char="")
metadata[[args$id_column]] <- as.character(metadata[[args$id_column]])
metadata <- metadata[metadata[[args$id_column]] %in% tree$tip.label, ]
if (nrow(metadata) == 0) stop("No matching IDs found in metadata.")

# Highlight logic
if (!is.null(args$highlight_file)) {
  highlight_ids <- readLines(args$highlight_file)
  metadata$should_show_point <- metadata[[args$id_column]] %in% highlight_ids
} else {
  metadata$should_show_point <- TRUE
}

# Color Map (Tip Points)
all_color_levels <- sort(unique(na.omit(metadata[[args$color_column]])))
color_palette <- scales::viridis_pal()(length(all_color_levels))
color_map <- setNames(color_palette, all_color_levels)

metadata$FinalLabel <- metadata[[args$species_column]]

# -----------------
# 4. PLOTTING SETUP
# -----------------
cat("Plotting...\n")
p <- ggtree(tree, layout = args$layout) %<+% metadata

# 1. Tip Labels
if (!args$hide_tip_labels) {
  p <- p + geom_tiplab(aes(label = FinalLabel), size=5, offset=0.05, align=TRUE, linesize=0.2, hjust=0)
}

# 2. Tip Points
p <- p + geom_tippoint(aes(color = !!sym(args$color_column), subset = should_show_point), size = 8) +
  scale_color_manual(name = args$color_column, values = color_map, na.value="grey50")

# 3. Support Values
cutoff <- args$support_cutoff
if (args$support_mode == "single") {
  p <- p + geom_point2(aes(subset = !isTip & !is.na(as.numeric(label)) & as.numeric(label) >= cutoff),
                       shape=23, size=3, fill="black", color="black") +
           geom_point2(aes(subset = !isTip & !is.na(as.numeric(label)) & as.numeric(label) < cutoff),
                       shape=23, size=3, fill="white", color="black")
} else if (args$support_mode == "dual") {
  p <- p + geom_point2(aes(subset = !isTip & 
                             as.numeric(sub("/.*", "", label)) >= cutoff & 
                             as.numeric(sub(".*/", "", label)) >= cutoff),
                       shape=23, size=3, fill="black", color="black") +
           geom_point2(aes(subset = !isTip & 
                             ((as.numeric(sub("/.*", "", label)) >= cutoff & as.numeric(sub(".*/", "", label)) < cutoff) |
                              (as.numeric(sub("/.*", "", label)) < cutoff & as.numeric(sub(".*/", "", label)) >= cutoff))),
                       shape=23, size=3, fill="gray", color="black") +
           geom_point2(aes(subset = !isTip & 
                             as.numeric(sub("/.*", "", label)) < cutoff & 
                             as.numeric(sub(".*/", "", label)) < cutoff),
                       shape=23, size=3, fill="white", color="black")
}

# 4. Clade Labels Calculation
clade_df <- NULL

# Strategy A: Manual File (Overrides sorting logic usually, as colors are hardcoded)
if (!is.null(args$clades_file)) {
  cat("Adding Manual Clade Labels from file...\n")
  clade_df <- read.table(args$clades_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  if(!"color" %in% colnames(clade_df)) clade_df$color <- "black"

# Strategy B: Automatic with Optional Depth Sorting
} else if (!is.null(args$auto_clade_label)) {
  cat("Adding Automatic Clade Labels from column:", args$auto_clade_label, "\n")
  
  target_col <- args$auto_clade_label
  valid_meta <- metadata[!is.na(metadata[[target_col]]), ]
  unique_groups <- unique(valid_meta[[target_col]])
  
  # --- SORTING LOGIC ---
  if (args$sort_by_depth) {
    cat("  -> Calculating distances from root to sort colors (Basal -> Distal)...\n")
    
    # Calculate node depths from root
    node_depths <- ape::node.depth.edgelength(tree)
    
    group_depths <- numeric(length(unique_groups))
    names(group_depths) <- unique_groups
    
    for (grp in unique_groups) {
      tips <- valid_meta[[args$id_column]][valid_meta[[target_col]] == grp]
      valid_tips <- intersect(tips, tree$tip.label)
      if (length(valid_tips) > 1) {
        mrca <- ape::getMRCA(tree, valid_tips)
        group_depths[grp] <- node_depths[mrca]
      } else if (length(valid_tips) == 1) {
        # If single tip, depth is tip depth
        tip_idx <- which(tree$tip.label == valid_tips)
        group_depths[grp] <- node_depths[tip_idx]
      } else {
        group_depths[grp] <- 0
      }
    }
    
    # Sort groups by depth
    sorted_groups <- names(sort(group_depths))
    
    # Assign colors using Viridis (or Turbo) along this sorted list
    # Turbo is great for distinct spectral colors
    clade_colors <- viridis::turbo(length(sorted_groups)) 
    names(clade_colors) <- sorted_groups
    
  } else {
    # Default: Alphabetical/Random Hue
    clade_colors <- scales::hue_pal()(length(unique_groups))
    names(clade_colors) <- unique_groups
  }
  
  # Map colors back to dataframe
  color_lookup_keys <- as.character(valid_meta[[target_col]])
  
  clade_df <- data.frame(
    tip = valid_meta[[args$id_column]],
    label = valid_meta[[target_col]],
    color = clade_colors[color_lookup_keys],
    stringsAsFactors = FALSE
  )
}

# 5. Clade Plotting (with Local Environment Fix)
if (!is.null(clade_df)) {
  clade_df$label <- as.character(clade_df$label)
  
  for (lbl in unique(clade_df$label)) {
    local({
      current_lbl <- lbl
      tips_in_group <- clade_df$tip[clade_df$label == current_lbl]
      valid_tips <- intersect(tips_in_group, tree$tip.label)
      
      if (length(valid_tips) > 1) {
        mrca_node <- ape::getMRCA(tree, valid_tips)
        current_color <- unique(clade_df$color[clade_df$label == current_lbl])[1]
        
        p <<- p + geom_cladelabel(
          node = mrca_node,
          label = current_lbl,
          color = current_color,
          offset = 0.2,      
          barsize = 2,       
          fontsize = 6,      
          align = TRUE       
        )
      }
    })
  }
}

p <- p + geom_treescale(fontsize=3.5, linesize=0.7) +
  theme(plot.margin = margin(20, 100, 20, 20), 
        legend.position = "bottom", 
        legend.title = element_text(size=24), 
        legend.text = element_text(size=22))

if (args$flip) p <- p + scale_y_reverse()

cat("Saving to:", args$output, "\n")
ggsave(p, filename=args$output, width=args$width, height=args$height)

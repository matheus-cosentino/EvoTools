#!/usr/bin/env Rscript

# -----------------
# 1. SETUP
# -----------------
# Ensure required packages are installed and loaded
required_packages <- c("argparse", "ggtree", "ggplot2", "ape", "dplyr",
                       "treeio", "tidytree", "scales", "phytools")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("ggtree", "treeio", "tidytree")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
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
# Define command-line arguments
parser <- ArgumentParser(description = "Plot a phylogenetic tree with various layouts and optional subtree extraction.")
parser$add_argument("-t", "--tree", type = "character", required = TRUE, help = "Path to the input tree file")
parser$add_argument("-m", "--metadata", type = "character", required = TRUE, help = "Path to the input metadata TSV file")
parser$add_argument("--tips_file", type="character", default=NULL,
                    help="Optional path to a text file with one tip label per line to create a subtree.")
parser$add_argument("-l", "--layout", type="character", default="rectangular",
                    choices=c("rectangular", "circular", "fan", "slanted", "unrooted"),
                    help="Tree layout type [default: %(default)s]")
parser$add_argument("--id_column", type = "character", default = "Acc", help = "Column in metadata that matches tree tips")
parser$add_argument("--species_column", type = "character", default = "Species", help = "Column in metadata with species names to use as labels")
parser$add_argument("--color_column", type = "character", default = "Genus", help = "Column in metadata to color tips by")
parser$add_argument("-o", "--output", type = "character", default = "Final_Tree.svg", help = "Path for the output plot file")
parser$add_argument("-W", "--width", type = "double", default = 10, help = "Width of the output plot")
parser$add_argument("-H", "--height", type = "double", default = 12, help = "Height of the output plot")
parser$add_argument("--flip", action="store_true", default=FALSE, help="Rotates the tree (inverts tip order).")
args <- parser$parse_args()


# -----------------
# 3. DATA LOADING & PREPARATION
# -----------------
cat("Loading tree file:", args$tree, "\n")
tree <- read.tree(args$tree)

# --- SUBTREE CREATION (OPTIONAL) ---
# Check if a tips file was provided to create a subtree
if (!is.null(args$tips_file)) {
  cat("Tips file provided. Pruning tree to create a subtree...\n")
  
  tips_to_keep <- readLines(args$tips_file)
  tips_in_tree <- intersect(tips_to_keep, tree$tip.label)
  
  if (length(tips_in_tree) == 0) {
    stop("CRITICAL ERROR: None of the tips in your --tips_file were found in the main tree.")
  }
  
  cat("Found", length(tips_in_tree), "matching tips. Pruning tree...\n")
  tree <- ape::keep.tip(tree, tips_in_tree)
}

cat("Midpoint rooting the tree...\n")
tree <- phytools::midpoint.root(tree)

cat("Loading metadata file:", args$metadata, "\n")
metadata <- read.csv(args$metadata, sep = "\t", header = TRUE, comment.char="")
metadata[[args$id_column]] <- as.character(metadata[[args$id_column]])

# --- CONSISTENT COLOR MAP CREATION ---
cat("Creating a consistent color map for '", args$color_column, "'...\n")
# Sort levels to ensure the color assignment is stable every time the script is run
all_color_levels <- sort(unique(na.omit(metadata[[args$color_column]])))
color_palette <- scales::viridis_pal()(length(all_color_levels))
color_map <- setNames(color_palette, all_color_levels)

cat("Setting tip labels from the specified column...\n")
metadata$FinalLabel <- metadata[[args$species_column]]
tip_labels <- tree$tip.label

# --- DEBUGGING & VALIDATION ---
cat("\n--- DEBUGGING ---\n")
cat("First 5 tip labels from the final tree object:\n")
print(head(tip_labels, 5))
cat("\nFirst 5 IDs from METADATA file (in column '", args$id_column, "'):\n")
print(head(metadata[[args$id_column]], 5))

# Filter the metadata to only include tips present in the final tree
metadata <- metadata[metadata[[args$id_column]] %in% tip_labels, ]

cat("\nNumber of matching IDs found between final tree and metadata:", nrow(metadata), "\n\n")
if (nrow(metadata) == 0) {
  stop("CRITICAL ERROR: No matching IDs found. Halting script. Check that tree tip labels exactly match the IDs in your metadata file.")
}


# -----------------
# 4. PLOTTING
# -----------------
cat("Generating plot with '", args$layout, "' layout...\n")

p <- ggtree(tree, layout = args$layout) %<+%
  metadata

p <- p +
  geom_tiplab(aes(label = FinalLabel),
              size = 5,
              align = TRUE,
              linesize = 0.2,
              offset = 0.05,
		husjt = 0)  +
  xlim(0, 1.05) +
  geom_tippoint(aes(color = !!sym(args$color_column)),
                size = 8) +
  
  # Apply the manual color map for consistent colors
  scale_color_manual(name = args$color_column, values = color_map, na.value="grey50") +

  # Plot conditional node support symbols
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label)) > 75 & as.numeric(sub(".*/", "", label)) > 75 & !isTip),
              colour = "black", size=3, fill = "black", shape = 23) +
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label)) < 75 & as.numeric(sub(".*/", "", label)) > 75 & !isTip),
              colour = "black", size=3, fill = "gray", shape = 23) +
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label)) > 75 & as.numeric(sub(".*/", "", label)) < 75 & !isTip),
              colour = "black", size=3, fill = "white", shape = 23) +
  
  geom_treescale(fontsize = 3.5, linesize = 0.7) +
  
  # Customize theme elements like legend
  theme(
    plot.margin = margin(20, 20, 20, 20),
    legend.position = "bottom",
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    legend.key.size = unit(1, "cm")
  ) 

if (args$flip) {
  cat("Flipping tree order...\n")
  p <- p + scale_y_reverse()
}

# -----------------
# 5. SAVING THE PLOT
# -----------------
cat("Saving final plot to:", args$output, "\n")
ggsave(plot = p, filename = args$output, width = args$width, height = args$height)

cat("Done!\n")

# --- 1. Install and Load Packages ---
# install.packages(c("rotl", "ape", "phytools")) # Run once if not installed
library(rotl)
library(ape)       # For tree manipulation (read.tree, write.tree, drop.tip)
library(phytools)  # For read.newick, and other helpful functions (optional)

# --- 2. Define Your Species List ---
anemonefish_species_list <- c(
  "Amphiprion clarkii",
  "Amphiprion frenatus",
  "Amphiprion ocellaris",
  "Amphiprion perideraion",
  "Amphiprion polymnus",
  "Amphiprion sandaracinos"
)

# --- 3. Resolve Species Names to OTL Taxon IDs (OTT IDs) ---
cat("Resolving species names to Open Tree of Life IDs...\n")
# Context helps narrow down the search. Pomacentridae is the family.
# "Ray-finned fishes" or "Animals" are broader good contexts.
resolved_taxa <- tnrs_match_names(names = anemonefish_species_list,
                                  context_name = "Animals") # Or "Ray-finned fishes"

# Check resolution quality
print(resolved_taxa)

# Extract valid OTT IDs
ott_ids <- resolved_taxa$ott_id[!is.na(resolved_taxa$ott_id)]

if (length(ott_ids) < length(anemonefish_species_list)) {
  warning("Not all species were resolved to OTT IDs. Problematic species:")
  print(anemonefish_species_list[!anemonefish_species_list %in% resolved_taxa$search_string[resolved_taxa$ott_id %in% ott_ids]])
  # You might need to manually check names on opentreeoflife.org or try different context_name
}
if (length(ott_ids) < 2) {
  stop("Less than 2 species resolved. Cannot form a tree.")
}

cat("Successfully resolved OTT IDs:", paste(ott_ids, collapse=", "), "\n")

# --- 4. Fetch the Induced Subtree from OTL ---
# This gets a tree containing only your specified taxa (and any necessary internal nodes)
cat("Fetching phylogeny from Open Tree of Life...\n")
# label_format = "name" gives tip labels as "Genus_species"
# label_format = "id" gives tip labels as "ottID"
# label_format = "name_and_id" gives "Genus_species_ottID"

# OTL will try to return a rooted tree with branch lengths if available from its synthesis.
# If the synthetic tree is very large, this might take a moment.
anemonefish_tree_otl <- tryCatch({
  tol_induced_subtree(ott_ids = ott_ids, label_format = "name")
}, error = function(e) {
  cat("Error fetching tree from OTL:", conditionMessage(e), "\n")
  return(NULL)
})

if (is.null(anemonefish_tree_otl)) {
  stop("Failed to retrieve tree from Open Tree of Life.")
}

cat("Tree retrieved successfully.\n")
plot(anemonefish_tree_otl, cex = 0.8, main = "Anemonefish Phylogeny (from OTL)")
axisPhylo() # Adds a scale bar if branch lengths are present

# --- 5. Inspect and Potentially Prune the Tree ---
# The retrieved tree might include more taxa if OTL needed them to connect your species.
# Or, it might perfectly match your list if they form a monophyletic group.
print("Tip labels in the OTL tree:")
print(anemonefish_tree_otl$tip.label)

# Prepare your list for matching (OTL uses underscores)
anemonefish_species_formatted_for_tips <- gsub(" ", "_", anemonefish_species_list)

# Identify tips in the OTL tree that are NOT in your desired list
tips_to_remove <- setdiff(anemonefish_tree_otl$tip.label, anemonefish_species_formatted_for_tips)

if (length(tips_to_remove) > 0) {
  cat("Pruning extra tips from the OTL tree:", paste(tips_to_remove, collapse=", "), "\n")
  anemonefish_tree_pruned <- drop.tip(anemonefish_tree_otl, tips_to_remove)
} else {
  cat("OTL tree already matches the target species or no extra tips to prune.\n")
  anemonefish_tree_pruned <- anemonefish_tree_otl
}

# Ensure we still have a valid tree
if (Ntip(anemonefish_tree_pruned) < 2) {
  stop("Pruning resulted in a tree with less than 2 tips. Check species name matching.")
}

plot(anemonefish_tree_pruned, cex = 0.8, main = "Pruned Anemonefish Phylogeny")
axisPhylo()

# This is your final 'phylo' object
final_anemonefish_tree <- anemonefish_tree_pruned

# --- 6. Save the Tree in Newick Format (Standard) ---
newick_filename <- "anemonefish_phylogeny_OTL.nwk"
write.tree(final_anemonefish_tree, file = newick_filename)
cat("Tree saved in Newick format to:", newick_filename, "\n")

# --- 7. Create a Basic NEXUS File (Similar to your example) ---
nexus_filename <- "anemonefish_phylogeny_OTL.nex"

# Get the Newick string for the tree with branch lengths
newick_string_for_nexus <- write.tree(final_anemonefish_tree)

# Construct the NEXUS file content
nexus_content <- c(
  "#NEXUS",
  "",
  "[This phylogeny for Anemonefish was obtained from the Open Tree of Life synthesis]",
  "[Date retrieved:", Sys.Date(), "]",
  "",
  "BEGIN TAXA;",
  paste0("    TITLE Anemonefish_Taxa;"),
  paste0("    DIMENSIONS NTAX=", Ntip(final_anemonefish_tree), ";"),
  "    TAXLABELS",
  # Ensure labels are just the names, without quotes if write.tree added them for NEXUS
  paste0("        ", final_anemonefish_tree$tip.label, collapse = "\n"),
  "    ;",
  "END;",
  "",
  "BEGIN TREES;",
  paste0("    TITLE Anemonefish_Tree_from_OTL;"),
  paste0("    LINK TAXA = Anemonefish_Taxa;"), # Good practice
  # [&R] indicates a rooted tree, which OTL typically provides for induced subtrees
  # If OTL provided an unrooted tree, you might omit [&R] or root it first.
  # The `tol_induced_subtree` function from rotl should return a rooted tree.
  paste0("    TREE Amphiprion_OTL_Tree = [&R] ", newick_string_for_nexus, ";"),
  "END;"
)

# Write the NEXUS content to a file
writeLines(nexus_content, con = nexus_filename)
cat("Tree saved in basic NEXUS format to:", nexus_filename, "\n")

cat("\n--- Script Finished ---\n")
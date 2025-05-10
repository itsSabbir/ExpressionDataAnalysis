# BCB420 Assignment 1 - Sabbir Hossain
# Enhanced R Script Version

# --- PREAMBLE ---
# This script performs an analysis of RNA-sequencing data (GSE193417)
# related to Corticotropin-Releasing Hormone (CRH) expressing cells
# in the subgenual anterior cingulate cortex (sgACC) of individuals
# with Major Depressive Disorder (MDD) and control subjects.

# --- BACKGROUND ---
# (Paraphrased from the paper and assignment description)
# A prior transcriptome meta-analysis found decreased CRH mRNA in MDD patients,
# suggesting impairment in cortical CRH-expressing (CRH+) cells.
# This study investigates CRH+ cells in human sgACC, their characteristics,
# and their relationship to MDD.
# sgACC samples from human volunteers (controls and MDD patients, n=6/group)
# were analyzed. CRH+ cells were identified using FISH for CRH and markers
# for excitatory (SLC17A7), inhibitory (GAD1), and interneuron subpopulations
# (PVALB, SST, VIP).
# RNA-sequencing was performed on sgACC CRH+ interneurons.
# Key findings from the paper:
# - GABAergic cells: ~80% of CRH+ cells. Glutamatergic: ~17.5%.
# - CRH+ GABAergic interneurons co-expressed VIP (52%), SST (7%), PVALB (7%).
# - MDD patients: Lower CRH mRNA in GABAergic interneurons, no change in cell density.
# - CRH+ interneuron transcriptome in MDD: Suggests decreased excitability,
#   less GABA release/reuptake.
# - Molecular alterations are not due to altered glucocorticoid feedback but
#   likely downstream of a common neurotrophic function modulator.
#
# Dataset: GSE193417 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193417)

# --- Global Options ---
# options(timeout = 300) # Increase timeout for downloads if needed

# --- 1. INQUIRE AND DOWNLOAD REQUIRED FILES FROM GEO ---

# --- 1.1. Install and Load Packages ---
# Check and install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# List of required packages
required_packages <- c(
    "GEOmetadb" = "BiocManager",   # For querying GEO metadata database
    "limma" = "BiocManager",       # For differential expression analysis (used for MA plots)
    "org.Hs.eg.db" = "BiocManager",# For human gene annotations
    "edgeR" = "BiocManager",       # For differential expression of digital gene expression data
    "tidyverse" = "CRAN",      # For data manipulation and visualization (includes ggplot2, dplyr)
    "ggplot2" = "CRAN",        # For plotting
    "AnnotationDbi" = "BiocManager", # For annotation database interface
    "readxl" = "CRAN",         # For reading Excel files (not strictly used here but good to have)
    "RColorBrewer" = "CRAN",   # For color palettes
    "dplyr" = "CRAN",          # For data manipulation
    "Biobase" = "BiocManager",     # Core Bioconductor data structures
    "biomaRt" = "BiocManager",     # For accessing Ensembl BioMart
    "GEOquery" = "BiocManager",    # For downloading and parsing GEO data
    "RSQLite" = "CRAN"         # For GEOmetadb (if used extensively, not primary here)
)

# Install missing packages
message("Checking and installing required packages...")
for (pkg_name in names(required_packages)) {
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
        message("Installing ", pkg_name, "...")
        if (required_packages[[pkg_name]] == "BiocManager") {
            BiocManager::install(pkg_name, ask = FALSE, update = FALSE)
        } else {
            install.packages(pkg_name)
        }
    }
}

# Load libraries
message("Loading required libraries...")
suppressPackageStartupMessages({
    library(Biobase)
    library(GEOquery)
    library(GEOmetadb) # Though not directly used for main workflow, good for deeper GEO exploration
    library(limma)
    library(org.Hs.eg.db)
    library(edgeR)
    library(tidyverse)
    library(ggplot2)
    library(AnnotationDbi)
    library(biomaRt)
    library(RColorBrewer)
    library(dplyr)
    library(RSQLite)
})

message("All packages checked and loaded.")

# --- 1.2. Download GEO Metadata ---
# Define GEO accession ID
geo_accession_id <- "GSE193417"
rds_file_path <- paste0(geo_accession_id, ".rds")

if (!file.exists(rds_file_path)) {
    message("Downloading GEO data for ", geo_accession_id, "...")
    # Get the GEO data; GSEMatrix=TRUE gets the expression data if available in series matrix
    # getGPL=FALSE to avoid downloading platform information separately if not needed immediately
    GSE193417_eset <- getGEO(geo_accession_id, GSEMatrix = TRUE, getGPL = FALSE)
    
    # GEO objects can be lists if multiple platforms are involved.
    # We need to select the correct ExpressionSet.
    # The original Rmd implies there might be multiple platforms (GPL16791 was mentioned for filtering).
    # However, GSE193417 only uses GPL15520. Let's keep it simple.
    if (length(GSE193417_eset) > 1) {
        # If multiple platforms, you might need to select based on platform ID
        # For GSE193417, it seems to be a single ExpressionSet in the list
        idx <- 1 # Or use grep for a specific GPL if needed, e.g., grep("GPL15520", names(GSE193417_eset))
        GSE193417 <- GSE193417_eset[[idx]]
    } else if (is.list(GSE193417_eset) && length(GSE193417_eset) == 1) {
        GSE193417 <- GSE193417_eset[[1]]
    } else if (is(GSE193417_eset, "ExpressionSet")) { # If getGEO directly returns an ExpressionSet
        GSE193417 <- GSE193417_eset
    } else {
        stop("Downloaded GEO object is not in the expected format (list of ExpressionSet or ExpressionSet).")
    }
    
    message("Saving downloaded GEO data to ", rds_file_path, "...")
    saveRDS(GSE193417, rds_file_path)
} else {
    message("Loading GEO data from local file: ", rds_file_path, "...")
    GSE193417 <- readRDS(rds_file_path)
}

# --- 1.3. Initial Dataset Verification ---
message("Verifying downloaded GEO data...")
if (!exists("GSE193417") || !is(GSE193417, "ExpressionSet")) {
    stop("PANIC: GSE193417 was not loaded or is not a valid ExpressionSet.
         Clear RDS file and run code again from the start. Can't continue.")
} else {
    message("GSE193417 ExpressionSet loaded successfully.")
}

# Print basic information about the ExpressionSet
print(GSE193417)

# Check number of samples
num_samples <- ncol(exprs(GSE193417)) # exprs() extracts the expression matrix
message("Number of samples in the series: ", num_samples)
if (num_samples != 12) {
    warning("Expected 12 samples, but found ", num_samples, ". Please verify.")
}

# Display sample names (GEO accession numbers for samples)
message("Sample names (GSM IDs):")
print(Biobase::sampleNames(GSE193417))

# --- 1.4. Download Supplementary Files ---
# Supplementary files often contain raw counts or processed data not in the series matrix.
supp_file_dir <- geo_accession_id
if (!dir.exists(supp_file_dir)) {
    message("Downloading supplementary files for ", geo_accession_id, "...")
    # This downloads files into a subdirectory named after the GEO ID
    gse_supp_files_info <- getGEOSuppFiles(geo_accession_id, fetch_files = TRUE, baseDir = getwd())
    # rownames(gse_supp_files_info) gives the full paths to downloaded files
    fnames <- rownames(gse_supp_files_info)
} else {
    message("Supplementary files directory '", supp_file_dir, "' already exists. Listing files...")
    # If directory exists, list files. Need to be careful if not all files were downloaded previously.
    # The original Rmd uses gsefiles$fname which is not standard output of getGEOSuppFiles
    # Let's list files in the directory
    existing_files <- list.files(path = supp_file_dir, full.names = TRUE)
    if (length(existing_files) > 0) {
      fnames <- existing_files
    } else {
      message("Directory ", supp_file_dir, " exists but is empty. Attempting re-download...")
      gse_supp_files_info <- getGEOSuppFiles(geo_accession_id, fetch_files = TRUE, baseDir = getwd())
      fnames <- rownames(gse_supp_files_info)
    }
}

message("Supplementary file paths:")
print(fnames)

# The assignment expects two specific supplementary files:
# 1. Raw count matrix (e.g., GSE193417_Raw_count_matrix_CRH.csv.gz)
# 2. Sample metadata (e.g., GSE193417_Sample_metadata.csv.gz)
# We need to identify the raw count matrix file.
# The first file listed is typically the one used for counts.
raw_counts_file_path <- fnames[grep("Raw_count_matrix", fnames, ignore.case = TRUE)]
if (length(raw_counts_file_path) == 0) {
    stop("Raw count matrix file not found among supplementary files.")
} else if (length(raw_counts_file_path) > 1) {
    warning("Multiple raw count matrix files found. Using the first one: ", raw_counts_file_path[1])
    raw_counts_file_path <- raw_counts_file_path[1]
} else {
     message("Identified raw counts file: ", raw_counts_file_path)
}

sample_metadata_file_path <- fnames[grep("Sample_metadata", fnames, ignore.case = TRUE)]
if (length(sample_metadata_file_path) == 0) {
    warning("Sample metadata file not found among supplementary files. Will rely on pData from ExpressionSet.")
    sample_metadata_file_path <- NULL # Set to NULL if not found
} else if (length(sample_metadata_file_path) > 1) {
    warning("Multiple sample metadata files found. Using the first one: ", sample_metadata_file_path[1])
    sample_metadata_file_path <- sample_metadata_file_path[1]
} else {
    message("Identified sample metadata file: ", sample_metadata_file_path)
}


# --- 2. ASSESS, EXPLORE, AND CLEAN SAMPLE FILE(S) AND METADATA, THEN MAP TO HUGO SYMBOLS ---

# --- 2.1. Load Platform and Experiment Metadata (Alternative to ExpressionSet pData) ---
message("Fetching detailed metadata using getGEO (GSEMatrix = FALSE)...")
gse_full_metadata <- getGEO(geo_accession_id, GSEMatrix = FALSE)

# Experiment metadata
experiment_metadata <- Meta(gse_full_metadata)
message("\n--- Experiment Metadata (from Meta(gse_full_metadata)) ---")
# Print some key fields
cat("Title:", experiment_metadata$title, "\n")
cat("Submission Date:", experiment_metadata$submission_date, "\n")
cat("Last Update Date:", experiment_metadata$last_update_date, "\n")
cat("Organism:", experiment_metadata$organism_ch1, "\n") # Or organism if available directly
cat("Overall Design:", experiment_metadata$overall_design, "\n")
# cat("Summary:", experiment_metadata$summary, "\n") # Often very long

# Platform metadata
# GPLList extracts platform information
platform_gpl_id <- names(GPLList(gse_full_metadata))[1] # Assuming one platform
platform_metadata <- Meta(GPLList(gse_full_metadata)[[platform_gpl_id]])
message("\n--- Platform Metadata (", platform_gpl_id, ") ---")
cat("Platform Title:", platform_metadata$title, "\n")
cat("Platform Organism:", platform_metadata$organism, "\n")
cat("Platform Technology:", platform_metadata$technology, "\n")
cat("Number of Samples on this Platform for this GSE:", length(platform_metadata$sample_id), "\n")


# --- 2.2. Load and Inspect Raw Count Data ---
message("Loading raw count data from: ", raw_counts_file_path)
# The file is a CSV. check.names=FALSE preserves original column names.
mddInMe_raw_counts <- read.csv(raw_counts_file_path, header = TRUE, check.names = FALSE, row.names = 1)
# The Rmd renames the first column to "EntrezID" if it's not row names.
# Here, row.names=1 means the first column is used as row names.
# Let's ensure row names are Ensembl IDs and create an EntrezID column if needed,
# or use row names directly. The original Rmd script sets `colnames(mddInMe)[1] <- "EntrezID"`.
# This implies the first column was data, not row names. Let's re-read as per Rmd.
mddInMe_raw_counts <- read.csv(raw_counts_file_path, header = TRUE, check.names = FALSE)
colnames(mddInMe_raw_counts)[1] <- "EntrezID" # Assumes first column contains Ensembl IDs (used as Entrez)

message("\n--- Initial Raw Count Data (mddInMe_raw_counts) ---")
cat("Dimensions (rows, cols):", dim(mddInMe_raw_counts), "\n")
message("First 5 rows and 5 columns:")
print(head(mddInMe_raw_counts[, 1:min(5, ncol(mddInMe_raw_counts))], 5))
message("Column names:")
print(colnames(mddInMe_raw_counts))
cat("Number of unique gene IDs (EntrezID column):", length(unique(mddInMe_raw_counts$EntrezID)), "\n")
if (length(unique(mddInMe_raw_counts$EntrezID)) != nrow(mddInMe_raw_counts)) {
    warning("Gene IDs in 'EntrezID' column are not all unique.")
}

# Check for NA values
message("Checking for NA values in raw counts...")
if(any(is.na(mddInMe_raw_counts))) {
    warning("NA values found in raw count data. The Rmd used na.omit for display but not reassignment.")
    # mddInMe_raw_counts <- na.omit(mddInMe_raw_counts) # Uncomment to remove rows with NAs
    # cat("Dimensions after na.omit:", dim(mddInMe_raw_counts), "\n")
} else {
    message("No NA values found in raw counts.")
}
# The original Rmd just showed na.omit output. We'll proceed with mddInMe_raw_counts as is for now.


# --- 3. DATA CLEANING AND PREPROCESSING --- (Original Rmd had section "4. Clean as required")

# --- 3.1. Filter Low-Count Genes ---
# Rationale: Genes with very low counts across most samples provide little information
# and can interfere with statistical analyses (e.g., variance estimation).
# A common filtering strategy: keep genes with > X counts in at least Y samples.
# The paper mentions "deleting low-expressing genes with fewer than ten reads
# and not found in more than two-thirds of the samples."
# The Rmd uses: sum(rowSums(mddInMe > 1) >= 6)
# This means: keep genes that have >1 count in at least 6 samples (half of the 12 samples).
# This is a reasonable heuristic.

# Separate gene IDs from count data
gene_ids_mddInMe <- mddInMe_raw_counts$EntrezID
count_data_mddInMe <- mddInMe_raw_counts[, -1] # All columns except EntrezID

# Apply filtering
# A gene must have more than 1 read count in at least 6 samples.
keep_genes_idx <- rowSums(count_data_mddInMe > 1) >= 6
filt_mddInMe_counts <- count_data_mddInMe[keep_genes_idx, ]
filt_mddInMe_gene_ids <- gene_ids_mddInMe[keep_genes_idx]

# Combine filtered gene IDs and counts
filt_mddInMe <- data.frame(EntrezID = filt_mddInMe_gene_ids, filt_mddInMe_counts, check.names = FALSE)

message("\n--- Gene Filtering ---")
cat("Initial number of genes:", nrow(mddInMe_raw_counts), "\n")
cat("Number of genes after filtering (count > 1 in >= 6 samples):", nrow(filt_mddInMe), "\n")
cat("Number of genes removed:", nrow(mddInMe_raw_counts) - nrow(filt_mddInMe), "\n")
cat("Dimensions of filtered data (filt_mddInMe):", dim(filt_mddInMe), "\n")

# --- 3.2. Create Sample Metadata ---
# The Rmd had a complex way then hardcoded. Let's create it systematically.
# Sample names are column names of the count matrix (excluding EntrezID)
sample_column_names <- colnames(filt_mddInMe)[-1]

samples_meta <- data.frame(FullSampleName = sample_column_names, stringsAsFactors = FALSE)
# Extract short names like "Hu1001"
samples_meta$ShortName <- gsub("\\.bam$", "", gsub("CRH-", "", samples_meta$FullSampleName))

# Define groups based on paper/Rmd description.
# Controls (Group 1 / "CON" / "CTL")
controls_short_names <- c("Hu1031", "Hu1047", "Hu1086", "Hu615", "Hu789", "Hu852")
# MDD (Group 2 / "MDD")
mdd_short_names <- c("Hu1001", "Hu513", "Hu600", "Hu809", "Hu863", "Hu943")

samples_meta$Group <- "Unknown" # Default
samples_meta$Group[samples_meta$ShortName %in% controls_short_names] <- "Control"
samples_meta$Group[samples_meta$ShortName %in% mdd_short_names] <- "MDD"
samples_meta$Group <- factor(samples_meta$Group, levels = c("Control", "MDD"))

# Check if all samples were assigned a group
if (any(samples_meta$Group == "Unknown")) {
    warning("Some samples could not be assigned to a group: ",
            paste(samples_meta$FullSampleName[samples_meta$Group == "Unknown"], collapse=", "))
}

message("\n--- Sample Metadata (samples_meta) ---")
print(samples_meta)
table(samples_meta$Group)


# --- 3.3. Map Ensembl IDs to HUGO Gene Symbols ---
# Using biomaRt
message("\n--- Mapping Ensembl IDs to HUGO Symbols using biomaRt ---")
# Check if Ensembl mart object exists or create it
if (!exists('ensembl_mart_hsapiens')) {
    message("Connecting to Ensembl BioMart (hsapiens_gene_ensembl)...")
    ensembl_mart_hsapiens <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
} else {
    message("Using existing Ensembl BioMart connection.")
}

# Get gene symbols for the Ensembl IDs in our filtered dataset
# The 'EntrezID' column in filt_mddInMe actually contains Ensembl IDs.
message("Fetching HGNC symbols for Ensembl IDs in 'filt_mddInMe$EntrezID'...")
gene_mappings <- getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'),
    filters = 'ensembl_gene_id',
    values = filt_mddInMe$EntrezID, # Use Ensembl IDs from the filtered data
    mart = ensembl_mart_hsapiens
)

message("BioMart query returned ", nrow(gene_mappings), " mappings.")
cat("Number of unique Ensembl IDs with mappings:", length(unique(gene_mappings$ensembl_gene_id)), "\n")
cat("Number of unique HGNC symbols mapped:", length(unique(gene_mappings$hgnc_symbol[!is.na(gene_mappings$hgnc_symbol) & gene_mappings$hgnc_symbol != ""])), "\n")

# --- 3.4. Merge Mappings and Handle Duplicates/Missing Symbols ---
# Join the gene symbols with the filtered count data
# Renaming for clarity before join:
colnames(gene_mappings)[colnames(gene_mappings) == 'ensembl_gene_id'] <- 'EntrezID' # Match column name

# Inner join to keep only genes that have both counts and mappings
FinalGeneFilter_with_ids <- dplyr::inner_join(gene_mappings, filt_mddInMe, by = "EntrezID")
message("Dimensions after joining with gene mappings: ", paste(dim(FinalGeneFilter_with_ids), collapse=" x "))

# Data cleaning steps from Rmd:
# 1. Remove specific hgnc_symbols (STRA6LP, LINC00856, POLR2J3, TBCE)
symbols_to_remove <- c("STRA6LP", "LINC00856", "POLR2J3", "TBCE")
FinalGeneFilter_cleaned <- FinalGeneFilter_with_ids %>%
    filter(!(hgnc_symbol %in% symbols_to_remove))
message("Dimensions after removing specific symbols: ", paste(dim(FinalGeneFilter_cleaned), collapse=" x "))

# 2. Remove rows with empty hgnc_symbol
FinalGeneFilter_cleaned <- FinalGeneFilter_cleaned %>%
    filter(hgnc_symbol != "" & !is.na(hgnc_symbol))
message("Dimensions after removing empty/NA HGNC symbols: ", paste(dim(FinalGeneFilter_cleaned), collapse=" x "))

# 3. Handle duplicated HGNC symbols:
# If an HGNC symbol maps to multiple Ensembl IDs, or vice-versa, this can be an issue.
# The Rmd's logic implies keeping one mapping. A common approach is to prioritize, e.g., by highest average expression.
# For simplicity, let's check for duplicates and then decide.
# Duplicated Ensembl IDs (should not happen if EntrezID was unique and getBM returns one mapping per Ensembl)
# duplicated_ensembl <- FinalGeneFilter_cleaned$EntrezID[duplicated(FinalGeneFilter_cleaned$EntrezID)]
# if(length(duplicated_ensembl) > 0) {
#     warning(length(duplicated_ensembl), " duplicated Ensembl IDs found after filtering. This may need attention.")
# }

# Duplicated HGNC symbols:
hgnc_counts <- table(FinalGeneFilter_cleaned$hgnc_symbol)
duplicated_hgnc <- names(hgnc_counts[hgnc_counts > 1])
if (length(duplicated_hgnc) > 0) {
    message(length(duplicated_hgnc), " HGNC symbols are duplicated (map to multiple Ensembl IDs or vice-versa).")
    # Example: take the first occurrence or one with max row sum of counts
    # FinalGeneFilter_cleaned <- FinalGeneFilter_cleaned %>%
    #    group_by(hgnc_symbol) %>%
    #    filter(row_number() == 1) %>% # or some other criterion like filter(sum(c_across(where(is.numeric))) == max(sum(c_across(where(is.numeric)))))
    #    ungroup()
    # For now, keeping them as per Rmd (Rmd implicitly kept them if they passed earlier filters)
    # The Rmd's `smmryGeneCts[which(smmryGeneCts > 1)]` just showed them.
}

# 4. Final filtering step from Rmd: `keep = rowSums(FinalGeneFilter[<numeric_cols>] >1) >= 6`
# This step was already applied to `mddInMe` to create `filt_mddInMe`.
# Applying it again to `FinalGeneFilter_cleaned` might be redundant or for stringency.
# Let's identify numeric count columns in FinalGeneFilter_cleaned
count_column_indices_final <- which(sapply(FinalGeneFilter_cleaned, is.numeric))
# Ensure these are actual sample count columns, not e.g. numeric IDs
# Assuming sample columns are the only numeric ones after EntrezID, hgnc_symbol, description
sample_data_cols_final <- colnames(FinalGeneFilter_cleaned)[count_column_indices_final]
# Exclude EntrezID if it was numeric and accidentally included. Here EntrezID is char.

# The numeric columns should be the sample data columns.
# In FinalGeneFilter_cleaned: EntrezID, hgnc_symbol, description, SAMPLENAME1, SAMPLENAME2 ...
# Let's assume sample names start after 'description'.
first_sample_col_index <- which(colnames(FinalGeneFilter_cleaned) == sample_column_names[1])
numeric_sample_data_final <- FinalGeneFilter_cleaned[, first_sample_col_index:ncol(FinalGeneFilter_cleaned)]

# Re-apply the rowSums filter.
# This is effectively `rowSums(FinalGeneFilter_cleaned[, 4:15] > 1) >= 6` if structure is Ens,Hgnc,Desc,Samples...
keep_final_idx <- rowSums(numeric_sample_data_final > 1) >= 6
FinalGeneFilter_final <- FinalGeneFilter_cleaned[keep_final_idx, ]

message("Dimensions after final rowSums filter: ", paste(dim(FinalGeneFilter_final), collapse=" x "))
cat("Final number of unique Ensembl IDs:", length(unique(FinalGeneFilter_final$EntrezID)), "\n")
cat("Final number of unique HGNC symbols:", length(unique(FinalGeneFilter_final$hgnc_symbol)), "\n")

# --- 3.5. Prepare Data for edgeR and Visualization ---
# The final dataset for downstream analysis will use `FinalGeneFilter_final`.
# Extract count matrix (numeric part) and gene annotation.
counts_for_analysis <- FinalGeneFilter_final[, first_sample_col_index:ncol(FinalGeneFilter_final)]
rownames(counts_for_analysis) <- FinalGeneFilter_final$hgnc_symbol # Use HGNC symbols as row names for DGEList

# Check for duplicated HGNC symbols as rownames
if(any(duplicated(rownames(counts_for_analysis)))){
    warning("Duplicated HGNC symbols found as rownames. This will cause issues with DGEList. Aggregating or renaming needed.")
    # Simplistic fix: make them unique
    rownames(counts_for_analysis) <- make.unique(rownames(counts_for_analysis))
    message("Made rownames unique using make.unique().")
}

# Gene annotations for DGEList
gene_annotations_for_analysis <- FinalGeneFilter_final[, c("EntrezID", "hgnc_symbol", "description")]
colnames(gene_annotations_for_analysis)[colnames(gene_annotations_for_analysis) == "EntrezID"] <- "EnsemblID"


# Calculate CPM (Counts Per Million) for visualizations
# Using the final filtered and cleaned count matrix `counts_for_analysis`
cpm_values <- edgeR::cpm(counts_for_analysis)
log2_cpm_values <- log2(cpm_values + 0.1) # Add offset to avoid log(0)

message("\n--- CPM Values ---")
cat("Dimensions of CPM matrix:", dim(cpm_values), "\n")
message("Summary of log2(CPM+0.1) values (first 5 samples):")
print(summary(log2_cpm_values[,1:min(5, ncol(log2_cpm_values))]))


# --- 4. DATA VISUALIZATION ---

# --- 4.1. MDS Plot (Multidimensional Scaling) ---
# MDS plot visualizes sample relationships based on their expression profiles.
# Uses log2 CPM values.
message("\nGenerating MDS Plot...")
# Ensure colors are distinct for groups
group_levels <- levels(samples_meta$Group)
num_groups <- length(group_levels)
plot_colors <- RColorBrewer::brewer.pal(n = max(3, num_groups), name = "Set1")[1:num_groups]

# Prepare labels for MDS plot (e.g., ShortName - Group)
mds_labels <- paste(samples_meta$ShortName, samples_meta$Group, sep=" - ")

# Generate plot (open in a new window if run interactively)
# if(interactive()) dev.new()
plotMDS(log2_cpm_values, # Using log2 transformed CPMs
        labels = mds_labels,
        col = plot_colors[samples_meta$Group],
        main = "MDS Plot of Normalized CRH+ Neuron RNA-Seq Samples",
        xlab = "Leading logFC Dimension 1",
        ylab = "Leading logFC Dimension 2",
        cex=0.8)
legend("topright",
       legend = group_levels,
       fill = plot_colors,
       title = "Group",
       cex = 0.8)
message("MDS plot generated.")
# Interpretation from Rmd:
# "Shows how control (NR) and MDD (R) groups interact.
# Large variability within and between groups, but some clustering.
# Indicates genes responsible for MDD response. More samples would be beneficial.
# Controls denoted with '-1' (here 'Control'), affected with '-2' (here 'MDD')."


# --- 4.2. Box Plots ---
# To visualize expression distribution across samples.
# Non-normalized data (using mddInMe_raw_counts, original scale with EntrezID)
message("\nGenerating Box Plot for non-normalized (raw) data...")
# Calculate CPM for the initial raw data (before filtering but after EntrezID handling)
# Need to ensure mddInMe_raw_counts is numeric only for CPM
raw_counts_numeric_part <- mddInMe_raw_counts[, -which(colnames(mddInMe_raw_counts) == "EntrezID")]
cpm_raw <- edgeR::cpm(raw_counts_numeric_part)
log2_cpm_raw <- log2(cpm_raw + 0.1) # Add offset

# if(interactive()) dev.new()
boxplot(log2_cpm_raw,
        xlab = "Samples", ylab = "log2(CPM + 0.1)",
        main = "Box Plot of Raw (Non-Normalized) RNA-Seq Data",
        las = 2, # Rotate x-axis labels
        cex.axis = 0.7,
        col = plot_colors[samples_meta$Group]) # Color by group
legend("topright", legend = group_levels, fill = plot_colors, title = "Group", cex = 0.8)
message("Box plot for non-normalized data generated.")

# Normalized data (using log2_cpm_values from FinalGeneFilter_final)
message("Generating Box Plot for normalized data...")
# if(interactive()) dev.new()
boxplot(log2_cpm_values,
        xlab = "Samples", ylab = "log2(CPM + 0.1)",
        main = "Box Plot of Normalized RNA-Seq Data",
        las = 2,
        cex.axis = 0.7,
        col = plot_colors[samples_meta$Group])
legend("topright", legend = group_levels, fill = plot_colors, title = "Group", cex = 0.8)
message("Box plot for normalized data generated.")
# Interpretation from Rmd:
# "Normalized values show a more aligned distribution of medians,
# indicating successful normalization in reducing systematic variations."


# --- 4.3. Density Plots ---
# To visualize expression density distributions for each sample.
message("\nGenerating Density Plot for non-normalized (raw) data...")
# if(interactive()) dev.new()
plot(density(log2_cpm_raw[,1]),
     col = plot_colors[samples_meta$Group[1]], # Color of the first sample's line
     lwd = 2,
     xlab = "log2(CPM + 0.1)", ylab = "Density",
     main = "Density Plot of Raw (Non-Normalized) RNA-Seq Data",
     ylim = c(0, max(sapply(1:ncol(log2_cpm_raw), function(i) max(density(log2_cpm_raw[,i])$y))))) # Dynamic ylim

for (i in 2:ncol(log2_cpm_raw)) {
    lines(density(log2_cpm_raw[,i]), col = plot_colors[samples_meta$Group[i]], lwd = 2)
}
legend("topright", legend = group_levels, fill = plot_colors, title = "Group (color by sample group)", cex = 0.8)
message("Density plot for non-normalized data generated.")


message("Generating Density Plot for normalized data...")
# if(interactive()) dev.new()
plot(density(log2_cpm_values[,1]),
     col = plot_colors[samples_meta$Group[1]],
     lwd = 2,
     xlab = "log2(CPM + 0.1)", ylab = "Density",
     main = "Density Plot of Normalized RNA-Seq Data",
     ylim = c(0, max(sapply(1:ncol(log2_cpm_values), function(i) max(density(log2_cpm_values[,i])$y)))))

for (i in 2:ncol(log2_cpm_values)) {
    lines(density(log2_cpm_values[,i]), col = plot_colors[samples_meta$Group[i]], lwd = 2)
}
legend("topright", legend = group_levels, fill = plot_colors, title = "Group (color by sample group)", cex = 0.8)
message("Density plot for normalized data generated.")
# Interpretation from Rmd:
# "Density curves show how normalization brings sample distributions closer,
#  reducing technical variability. Peaks and shapes become more similar."

# --- 4.4. MA Plots (Ratio-Intensity Plots) ---
# To visualize expression differences between two samples (or groups).
# The Rmd compares first vs last sample of mddInMe (raw) and FinalGeneFilter (filtered).
# This corresponds to CRH-Hu1001.bam vs CRH-Hu943.bam
sample1_name <- "CRH-Hu1001.bam"
sample2_name <- "CRH-Hu943.bam"

# MA Plot for raw data (mddInMe_raw_counts)
message("\nGenerating MA Plot for raw data (", sample1_name, " vs ", sample2_name, ")...")
idx_s1_raw <- which(colnames(raw_counts_numeric_part) == sample1_name)
idx_s2_raw <- which(colnames(raw_counts_numeric_part) == sample2_name)

if (length(idx_s1_raw) == 1 && length(idx_s2_raw) == 1) {
    # if(interactive()) dev.new()
    # Using raw counts directly for plotMA (it does log transform internally)
    # Add pseudo-count of 0.5 to raw counts before log for stability
    plotMA(data.frame(S1 = raw_counts_numeric_part[, idx_s1_raw] + 0.5,
                      S2 = raw_counts_numeric_part[, idx_s2_raw] + 0.5),
           main = paste("MA Plot (Raw Counts):", sample1_name, "(y-axis) vs", sample2_name, "(x-axis avg)"),
           xlab = "Average log-expression (A)",
           ylab = "Expression log-ratio (M)")
    abline(h = 0, col = "red", lty = 2)
    message("MA plot for raw data generated.")
} else {
    warning("Could not find specified samples for MA plot in raw data.")
}


# MA Plot for normalized data (using `counts_for_analysis`, which has HGNC rownames)
message("Generating MA Plot for normalized data (", sample1_name, " vs ", sample2_name, ")...")
idx_s1_norm <- which(colnames(counts_for_analysis) == sample1_name)
idx_s2_norm <- which(colnames(counts_for_analysis) == sample2_name)

if (length(idx_s1_norm) == 1 && length(idx_s2_norm) == 1) {
    # if(interactive()) dev.new()
    # counts_for_analysis are filtered counts, not yet CPM or log transformed for this plotMA
    # plotMA can take a DGEList or a matrix of counts.
    # Add pseudo-count of 0.5 to these filtered counts
    plotMA(data.frame(S1 = counts_for_analysis[, idx_s1_norm] + 0.5,
                      S2 = counts_for_analysis[, idx_s2_norm] + 0.5),
           main = paste("MA Plot (Filtered Counts):", sample1_name, "(y-axis) vs", sample2_name, "(x-axis avg)"),
           xlab = "Average log-expression (A)",
           ylab = "Expression log-ratio (M)")
    abline(h = 0, col = "blue", lty = 2)
    message("MA plot for normalized (filtered) data generated.")
} else {
    warning("Could not find specified samples for MA plot in normalized data.")
}
# Interpretation from Rmd:
# "MA plots show log-fold change (M) vs average expression (A).
# Ideally, most points cluster around M=0. Normalization helps achieve this."

# --- 5. INTERPRETATION (from Rmd, as comments) ---
# This section summarizes findings and answers specific questions from the assignment.

# **What are the control and test conditions of the dataset?**
# Controls (Group 'Control' in this script, originally '-1'):
#   CRH-Hu1031, CRH-Hu1047, CRH-Hu1086, CRH-Hu615, CRH-Hu789, CRH-Hu852
# Test Conditions (MDD Group 'MDD', originally '-2'):
#   CRH-Hu1001, CRH-Hu513, CRH-Hu600, CRH-Hu809, CRH-Hu863, CRH-Hu943
# Meaning: Controls are individuals without MDD, expected to have normal CRH+ cell function.
# Test (MDD) individuals are those with Major Depressive Disorder, hypothesized to have
# altered CRH+ cell function/expression.
# The study labeled sgACC cells using FISH. MDD was thought to affect CRH due to
# glucocorticoid feedback, but this paper suggests neurotrophic function degradation.

# **Why is the dataset of interest to you?** (User's personal motivation)
# "As someone who struggles with mental health issues, I am always interested in these topics.
# It can only make us more aware and better for the long term once we are able to understand
# how to take care of our mind, and not just our physical bodies too.
# Broken bones almost always heal, but broken souls do not. Or so I've heard."

# **Were there expression values that were not unique for specific genes? How did you handle these?**
# (This refers to situations where one Ensembl ID might map to multiple gene symbols, or one symbol to multiple Ensembl IDs,
# or if the same gene ID appears multiple times in the input.)
# - Original gene IDs (Ensembl in 'EntrezID' column) in `mddInMe_raw_counts` were checked for uniqueness.
#   If not unique, this script would have warned.
# - During biomaRt mapping, multiple Ensembl IDs can map to the same HGNC symbol, or one Ensembl ID
#   might have multiple (aliased) HGNC symbols.
# - This script:
#   - Filters out specific problematic symbols (STRA6LP, etc.).
#   - Removes rows with empty or NA HGNC symbols.
#   - If duplicated HGNC symbols remain (e.g., from different Ensembl IDs mapping to the same HGNC),
#     this script currently keeps them if they pass filters (as Rmd did), but warns.
#     If they become rownames for DGEList, `make.unique()` is applied.
#   - The paper's count was 15,472 genes. Our `FinalGeneFilter_final` has `r nrow(FinalGeneFilter_final)` genes.

# **Were there expression values that could not be mapped to current HUGO symbols?**
# Yes, some Ensembl IDs might not have a corresponding HGNC symbol in BioMart,
# or the symbol might be empty/NA. These are typically filtered out.
# `gene_mappings` from biomaRt shows `r sum(gene_mappings$hgnc_symbol == "" | is.na(gene_mappings$hgnc_symbol))`
# Ensembl IDs that initially had empty/NA HGNC symbols. These were removed.
# Original number of Ensembl IDs in `filt_mddInMe`: `r length(unique(filt_mddInMe$EntrezID))`.
# Number of Ensembl IDs in `FinalGeneFilter_final` (mapped and kept): `r length(unique(FinalGeneFilter_final$EnsemblID))`.

# **How many outliers were removed?**
# The term "outliers" here might refer to genes removed during filtering, or sample outliers.
# Gene filtering:
#   Initial genes: `r nrow(mddInMe_raw_counts)`
#   Final genes after all filtering: `r nrow(FinalGeneFilter_final)`
#   Total genes removed: `r nrow(mddInMe_raw_counts) - nrow(FinalGeneFilter_final)`
#   Percentage of genes removed: `r round((1 - nrow(FinalGeneFilter_final) / nrow(mddInMe_raw_counts)) * 100, 2)`%
# Sample outliers were not explicitly removed in this script (MDS plot helps identify them).

# **How did you handle replicates?**
# This dataset does not appear to have technical replicates for the same biological sample.
# It has biological replicates (6 control individuals, 6 MDD individuals).
# These biological replicates are kept and used to define groups for comparison.

# **What is the final coverage of your dataset?**
# The final number of genes included in the analysis (rows in `counts_for_analysis`) is:
# `r nrow(counts_for_analysis)` unique (or made unique) HGNC symbols.
# The paper states "This study looked at 15,472 genes in total."
# Our final count is `r nrow(counts_for_analysis)`, which is `r ifelse(nrow(counts_for_analysis) > 15472, "higher", "lower")`
# than the paper's count. Differences can arise from versions of annotation databases,
# specific filtering criteria, and handling of ambiguous mappings.

# --- Final Processed Data Overview ---
message("\n--- Final Processed Data (counts_for_analysis) ---")
cat("Dimensions (genes, samples):", dim(counts_for_analysis), "\n")
message("First 5 rows and 3 columns of the final count matrix for analysis:")
print(head(counts_for_analysis[,1:min(3, ncol(counts_for_analysis))], 5))

message("Gene annotations for the final matrix (gene_annotations_for_analysis):")
print(head(gene_annotations_for_analysis))

# --- 6. CITATIONS (from Rmd) ---
# 1. R programming for data science: Peng, R. D. (2020). R Programming for Data Science. Leanpub.
#    (Specific link was to a case study: https://bookdown.org/rdpeng/rprogdatascience/data-analysis-case-study-changes-in-fine-particle-air-pollution-in-the-u-s-.html)
# 2. Lecture modules: (BCB420 course specific) https://q.utoronto.ca/courses/248455/modules
# 3. Oh, Hyunjung & Newton, Dwight & Lewis, David & Sibille, Etienne. (2022).
#    Lower Levels of GABAergic Function Markers in Corticotropin-Releasing Hormone-Expressing
#    Neurons in the sgACC of Human Subjects With Depression. Frontiers in Psychiatry. 13.
#    DOI: 10.3389/fpsyt.2022.827972.

message("\n--- Software Package Citations ---")
message("To cite GEOmetadb in publications, use:")
print(citation("GEOmetadb"))
message("\nTo cite GEOquery in publications, use:")
print(citation("GEOquery"))
message("\nTo cite biomaRt in publications, use:")
print(citation("biomaRt"))
message("\nTo cite edgeR in publications, use:")
print(citation("edgeR"))
message("\nTo cite limma in publications, use:")
print(citation("limma"))
message("\nTo cite dplyr in publications, use:")
print(citation("dplyr"))
message("\nTo cite ggplot2 in publications, use:")
print(citation("ggplot2"))

message("\n--- End of Analysis Script ---")
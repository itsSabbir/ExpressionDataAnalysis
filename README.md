# Expression Data Analysis: Human Gene Expression Normalization

Welcome to the "Expression Data Analysis" repository, where we delve into the normalization and analysis of human gene expression datasets. This project is part of my Bioinformatics course assignment and focuses on selecting, cleaning, normalizing, and analyzing a gene expression dataset to ensure accurate subsequent analyses. The final output is a well-documented R Notebook, complemented by a detailed README to guide you through the process and findings.

## Objective

The objective of this project is to enhance our understanding of selecting and analyzing gene expression datasets using tools like GEO2R, pivotal in bioinformatics. The specific focus for this analysis was on dataset GSE193417, which involves studying gene expression related to depressive disorders in human subjects.

## Time Management

- **Estimated Time**: 5 hours
- **Time Used**: 15 hours
- **Project Duration**: Started on 2022/02/11 and completed on 2022/03/01

## Repository Overview

This repository contains all the resources and outputs related to the expression dataset analysis. The analysis pipeline is implemented in an R Notebook, which is thoroughly documented to reflect each step of the process.

### Key Features

- Selection of a high-quality human gene expression dataset
- Data cleaning and mapping to HUGO gene symbols
- Application of normalization techniques to ensure data quality
- Detailed analysis and documentation of the dataset's control and test conditions

## Dataset Selection

The chosen dataset, [GSE193417](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193417), corresponds to a study titled "Lower Levels of GABAergic Function Markers in Corticotropin-Releasing Hormone-Expressing Neurons in the sgACC of Human Subjects With Depression." This dataset was selected based on the following criteria:

- Comprehensive gene coverage
- Recent data collection (within the last ten years)
- Sufficient biological replicates under interesting experimental conditions

## Data Cleaning and Mapping

The data was initially downloaded using the `GEOquery` package to maintain consistency and reproducibility. Key steps included:

- Assessing data quality for both control and test samples
- Mapping gene expression data to HUGO symbols, addressing challenges like unmapped rows and non-unique mappings
- Cleaning the data by potentially removing outliers based on a justified rationale

## Data Normalization

Normalization was carried out using methods best suited for the dataset type. The method chosen, its rationale, and the effects of normalization are thoroughly documented and demonstrated through visualizations.

## Analysis and Documentation

The final dataset consists of a dataframe where each row represents a unique HUGO gene symbol, crucial for subsequent analyses. Detailed responses are provided for:

- Description of control and test conditions
- Rationale behind dataset selection
- Challenges encountered with non-unique gene expressions and their resolutions
- Outlier analysis and handling
- Replicate management and overall dataset coverage

## Results and Discussion

This assignment provided valuable insights into the complexities of data analysis in bioinformatics. Particularly, it underscored the importance of precise data handling and the challenges of working with raw genomic data. The project not only advanced my technical skills but also deepened my understanding of the biological implications of gene expression in psychiatric conditions.

## Repository Contents

- **R Notebook (.Rmd and .html formats)**: Contains all the codes, explanations, and outputs.
- **Data/**: Folder containing processed datasets used in the notebook (note: raw data is not uploaded due to size constraints but can be downloaded directly via the notebook).
- **docs/**: Additional documentation and references.

## Further Reading and Resources

- **Quantile Normalization:** Bolstad et al. (2003)...
- **RNA-seq Analysis:** Bioconductor workflows...
- **HUGO Gene Nomenclature Committee:** Authoritative source for gene symbols...
- **Microarray Normalization Strategies:** Cheng et al. (2016)...

## How to Use This Repository

To replicate the analysis or explore the methodologies:

1. Clone this repository.
2. Run the R Notebook to understand each step of the data handling and analysis process.
3. Explore the HTML output for a comprehensive breakdown of the tasks and findings.


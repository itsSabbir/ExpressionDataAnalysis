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

Here's a detailed "Further Reading and Resources" section for your README.md, including complete references and descriptions to aid users in expanding their understanding of the methodologies and tools used in your project:


## Further Reading and Resources

- **Quantile Normalization**: Explore various normalization methods for high-density oligonucleotide array data focused on reducing variance and bias. See Bolstad et al. (2003), "A comparison of normalization methods for high-density oligonucleotide array data based on variance and bias." Bioinformatics, 19(2), 185-193. The paper discusses how different approaches can influence the outcome of data analysis. [Read the paper on PubMed](https://pubmed.ncbi.nlm.nih.gov/12538238/), [DOI link](https://doi.org/10.1093/bioinformatics/19.2.185).

- **RNA-seq Analysis**: Learn about RNA-seq differential expression analysis using Bioconductor packages such as limma, Glimma, and edgeR. These workflows provide a comprehensive guide from data preparation through to exploratory analysis and statistical testing. For a complete workflow, refer to the Bioconductor's RNA-seq workflow using **edgeR** and **DESeq2**, which are essential for gene-level exploratory analysis and differential expression.

- **HUGO Gene Nomenclature Committee**: The HGNC provides the authoritative source for human gene symbols and nomenclature, offering tools for searching synonyms, aliases, and gene families. This resource is vital for ensuring data is accurately labeled and consistent with current scientific standards. Visit the [HGNC website](https://www.genenames.org/) to explore gene symbols and their related information.

- **Microarray Normalization Strategies**: Cheng et al. (2016) introduce "CrossNorm," a novel normalization strategy for microarray data in cancer studies that aims to improve the accuracy of differential expression analysis. This paper discusses the advantages of using CrossNorm over traditional methods, particularly in the context of cancer research. The study is detailed in "CrossNorm: A Novel Normalization Strategy for Microarray Data in Cancers," Scientific Reports, 6, Article number: 18898. [Read the paper on PubMed](https://pubmed.ncbi.nlm.nih.gov/26732145/), [DOI link](https://doi.org/10.1038/srep18898).


## How to Use This Repository

To replicate the analysis or explore the methodologies:

1. Clone this repository.
2. Run the R Notebook to understand each step of the data handling and analysis process.
3. Explore the HTML output for a comprehensive breakdown of the tasks and findings.


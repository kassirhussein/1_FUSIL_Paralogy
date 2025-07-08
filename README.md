# 🔬 FUSIL and Gene Paralogy Analysis

This repository contains an R-based analysis of human genes categorized by their FUSIL (Functional Understanding of Susceptibility In Loss-of-function) classification and their paralogues. It explores relationships between gene essentiality, paralogue presence, sequence similarity, and evolutionary ancestry (MCRA).

## 📌 Key Features

- Integration of FUSIL data, BioMart paralogue annotations, and protein-coding gene lists  
- Comparative analysis of genes with and without paralogues across FUSIL bins  
- Visualization of:
  - Gene and paralogue distributions by FUSIL category  
  - Similarity thresholds and paralogue presence  
  - Alluvial plots showing FUSIL transitions (gene → paralogue)  
  - MCRA (Most Conserved Recent Ancestor) subtype enrichment  
- Reproducible and well-documented R scripts using tidyverse workflows and `ggplot2` visualizations

## 📁 Data Sources

- Ensembl BioMart (via `biomaRt`)  
- Custom TSV/CSV files:
  - FUSIL dataset  
  - Protein-coding gene list  
  - Human gene paralogues  

## 📊 Tools & Libraries

- `tidyverse`  
- `ggplot2`  
- `readxl`  
- `writexl`  
- `biomaRt`  
- `ggalluvial`  
- `reshape2`  

## 🧪 Getting Started

1. Clone this repo  
2. Open the `.R` or `.Rmd` file in RStudio  
3. Install required packages (see list above)  
4. Replace the file paths with your local file locations  
5. Run the code chunk-by-chunk or knit the R Markdown report  

## 📈 Example Outputs

> Include plots of:
> - Gene counts by FUSIL  
> - Paralogues by similarity threshold  
> - FUSIL → Paralog FUSIL Sankey plot  
> - MCRA subtype distributions  

## 🧬 Author

Your Name (or GitHub username)  
Affiliation (optional)  
Created with ❤️ in R

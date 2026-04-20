# Human eQTL Data Visualization 

This repository presents an exploratory data visualization project based on the eQTL resource from **Lappalainen et al., Nature (2013)**. The aim of the project is to explore the distribution, strength, and genomic organization of expression quantitative trait loci (eQTLs) across human chromosomes using clear statistical summaries and visual analysis.

## Project overview

Expression quantitative trait loci (eQTLs) are genomic variants associated with differences in gene expression. They provide an important link between genetic variation and molecular phenotype, helping explain how sequence differences may influence transcriptional regulation.

This project uses the eQTL dataset described in **Lappalainen et al. (2013)** and focuses on the visual exploration of chromosome-level patterns, association strength, and differences across quantitative traits. The repository was developed as a data visualization project in bioinformatics and emphasizes interpretation through informative plots rather than predictive modeling.

## Dataset

The repository is based on the eQTL file set from **Lappalainen et al., Nature 2013**.

### Sample set and sample size
- **EUR**: 373 samples
- **YRI**: 89 samples

### Quantitative traits included
- **Exon**
- **Gene**
- **Transcript ratio (trratio)**
- **Transcribed repetitive element (repeat)**
- **miRNA (mi)**

### Associations included in the dataset
The file contains:
- all associations below **false discovery rate (FDR) 5%**, or
- the **best association per gene** for exon and transcript ratio data, or
- the **best association per unit** for gene, repeat, and miRNA data.

When several variants had the same best p-value, one of them was selected randomly.

## eQTL file format

The dataset contains variant–quantitative trait association information with the following columns:

1. **SNP_ID** – Variant identifier according to dbSNP137; position-based identifier for variants not included in dbSNP  
2. **ID** – Null (`-`)  
3. **GENE_ID** – Gene identifier according to Gencode v12, miRBase v18, or repeat start-site annotation  
4. **PROBE_ID** – Quantitative trait identifier  
   - For exons: `GENEID_ExonStartPosition_ExonEndPosition`
   - For transcript ratios: transcript identifier according to Gencode v12  
5. **CHR_SNP** – Chromosome of the variant  
6. **CHR_GENE** – Chromosome of the quantitative trait  
7. **SNPpos** – Position of the variant  
8. **TSSpos** – Transcription start site of the gene or quantitative trait  
9. **Distance** – Absolute distance between SNP position and TSS  
10. **rvalue** – Spearman rank correlation rho; the sign indicates the direction of the nonreference allele effect  
11. **pvalue** – Association p-value  
12. **log10pvalue** – `-log10(pvalue)`

## Objectives

- Import and inspect the eQTL dataset
- Explore the genomic distribution of eQTLs across chromosomes
- Visualize association strength using p-value based summaries
- Compare chromosome-level patterns using different plot types
- Identify chromosomes with stronger or denser regulatory signals
- Produce a clean and interpretable visual summary of the dataset

## Main analyses

The project includes visualizations such as:

- eQTL counts by chromosome
- Density and distribution plots of association strength
- Boxplots and violin plots of `log10pvalue`
- Heatmaps summarizing chromosome-level significance patterns
- Exploratory comparisons of genomic regions with stronger regulatory enrichment

## Repository structure

```bash
.
├── README.md
├── human_eqtl_data_visualization.Rmd
├── human_eqtl_data_visualization.html
├── data/
├── figures/
```

## How to run

1. Clone the repository:

   ```bash
   git clone https://github.com/your-username/human-eqtl-data-visualization.git
   cd human-eqtl-data-visualization
   ```

2. Open the `.Rmd` file in **RStudio**.

3. Install the required packages if needed:

   ```r
   install.packages(c(
     "ggplot2", "dplyr", "tidyr", "readr", "tibble",
     "patchwork", "viridis", "RColorBrewer"
   ))
   ```

4. Place the input eQTL file inside the `data/` folder.

5. Knit the R Markdown document to HTML.

## Results and interpretation

The analysis is designed to provide an intuitive overview of how eQTL associations are distributed across the human genome. Particular attention is given to chromosomes showing higher densities of significant associations or stronger signals according to `log10pvalue`, helping highlight regions with potentially relevant regulatory activity.

## Academic context

This repository was developed as part of an academic data visualization project in bioinformatics.

## Authors

- Eira Fontanals Muñoz
- Òscar Contreras Parejo

## Reference

Lappalainen, T. et al. (2013). *Transcriptome and genome sequencing uncovers functional variation in humans*. **Nature**.

# BMI 5333: Systems Medicine: Principles and Practice

This repository contains details about the **bioinformatics project** as part of the BMI 5333: Systems Medicine: Principles and Practice course. The study, titled **Microarray-based Transcriptomic Analysis and Signature Extraction in Mice Lungs following 4 and 6 months of Cigarette Smoke (CS) Exposure**, aims to investigate the transcriptomic changes in gene expression in mouse lung tissue resulting from prolonged exposure to CS at different time intervals of four and six months, focusing on differentially expressed genes (DEGs), pathway enrichment, network analysis and hub gene identification using a publicly available microarray dataset obtained from NCBI GEO.

---

# Microarray-based Transcriptomic Analysis and Signature Extraction in Mice Lungs following 4 and 6 months of Cigarette Smoke Exposure

## Project Overview

### 1. Data Retrieval
The mouse microarray dataset with the accession ID **[GSE52509](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52509)** was accessed through the **National Center for Biotechnology Information Gene Expression Omnibus (NCBI GEO)**. 
- Data retrieval was conducted using the **GEOquery** R package.
- **Sample groups:** 4 month Control, 4 month CS-exposed samples, 6 month Control, and 6 month CS-exposed samples (n=3 per group).

---

### 2. Differential Expression Analysis
- **Tools used:** limma R package for differential expression analysis and EnhancedVolcano R package.
    
  1. 4 month CS-exposed vs. Control
  2. 6 month CS-exposed vs. Control
     
 - **Criteria for Differentially Expressed Genes (DEGs):**
   - False Discovery Rate (FDR) ≤ 0.05
   - |Fold Change| > 1.5
     
- R scripts can be found in this [folder](scripts).
  
**Figures:**
- [4 month CS-exposed vs. Control Volcano Plot](Assets/1monthCSvsAir.png)
- [6 month CS-exposed vs. Control Volcano Plot](Assets/3monthCSvsAir.png)
  
  ---

### 3. Functional Enrichment Analysis
- **Tools used:**
  - **[Metascape](https://metascape.org/gp/index.html#/main/step1)** was used for comparative analysis across experimental timepoints.
  - GO Biological Process (BP), Cellular Compartment (CC), and Molecular Function (MF) and KEGG pathway enrichments.

- **Key Parameters for Enrichment Significance:**
  - Minimum gene overlap = 3
  - Enrichment p-value cutoff = 0.05
  - Minimum enrichment factor = 1.5
  
**Figures:**
- [Comparative GO enrichment analysis-Heatmap](Assets/HeatmapSelectedGO.png)
- [Comparative KEGG enrichment analysis-Heatmap](Assets/HeatmapSelectedGO.png)
  
---

### 4. Network Analysis and Hub Gene Identification
- **Tools used:**
  - **[StringDB version 12](https://string-db.org/)** for protein-protein interaction (PPI) network visualization.
  - **[Cytoscape](https://cytoscape.org/)** for network visualization and hub gene analysis.
    - Maximal clique centrality (MCC) algorithm from the cytoHubba plugin was utilized to determine the top ten hub genes in each PPI.

**Figures:**
- [StringDB Inflammatory response PPI](Assets/inflammatoryresponse_hubgene.png)
- [Top ten hub genes in the Inflammatory response PPI](Assets/cytokineproduction_hubgene.png)
- [StringDB Cytokine-cytokine receptor interaction PPI](Assets/chemotaxis_regulation_hubgene.png)
- [Top ten hub genes in the Cytokine-cytokine receptor interaction PPI](Assets/cellmigration_hubgenes.png)

  
---

## Citation
If you use the tools or dataset mentioned in this repository in your research, please cite the following references:

- Ashburner, M., Ball, C. A., Blake, J. A., Botstein, D., Butler, H., Cherry, J. M., Davis, A. P., Dolinski, K., Dwight, S. S., Eppig, J. T., Harris, M. A., Hill, D. P., Issel-Tarver, L., Kasarskis, A., Lewis, S., Matese, J. C., Richardson, J. E., Ringwald, M., Rubin, G. M., & Sherlock, G. (2000). Gene Ontology: tool for the unification of biology. Nature Genetics 2000 25:1, 25(1), 25–29. https://doi.org/10.1038/75556
  
- Barrett, T., Wilhite, S. E., Ledoux, P., Evangelista, C., Kim, I. F., Tomashevsky, M., Marshall, K. A., Phillippy, K. H., Sherman, P. M., Holko, M., Yefanov, A., Lee, H., Zhang, N., Robertson, C. L., Serova, N., Davis, S., & Soboleva, A. (2013). NCBI GEO: archive for functional genomics data sets—update. Nucleic Acids Research, 41(D1), D991–D995. https://doi.org/10.1093/NAR/GKS1193

- Blighe K, Rana S, Lewis M (2024). EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling. R package version 1.24.0, https://github.com/kevinblighe/EnhancedVolcano

- Chin, C. H., Chen, S. H., Wu, H. H., Ho, C. W., Ko, M. T., & Lin, C. Y. (2014). cytoHubba: Identifying hub objects and sub-networks from complex interactome. BMC Systems Biology, 8(4), 1–7. https://doi.org/10.1186/1752-0509-8-S4-S11/TABLES/4

- Kanehisa, M. (2019). Toward understanding the origin and evolution of cellular organisms. Protein Science, 28(11), 1947–1951. https://doi.org/10.1002/PRO.3715
  
- Kanehisa, M., Furumichi, M., Sato, Y., Kawashima, M., & Ishiguro-Watanabe, M. (2023). KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Research, 51(D1), D587–D592. https://doi.org/10.1093/NAR/GKAC963

- Kanehisa, M., & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Research, 28(1), 27–30. https://doi.org/10.1093/NAR/28.1.27

- Li, C., & Xu, J. (2019). Feature selection with the Fisher score followed by the Maximal Clique Centrality algorithm can accurately identify the hub genes of hepatocellular carcinoma. Scientific Reports, 9(1), 1–11. https://doi.org/10.1038/s41598-019-53471-0

- Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res [Internet]. 2015 Apr 20 [cited 2023 Dec 7];43(7):e47–e47. Available from: https://dx.doi.org/10.1093/nar/gkv007

- Sean D, Meltzer PS. GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. Bioinformatics [Internet]. 2007 Jul 15 [cited 2023 Dec 7];23(14):1846–7. Available from: https://dx.doi.org/10.1093/bioinformatics/btm254

- Shannon, P., Markiel, A., Ozier, O., Baliga, N. S., Wang, J. T., Ramage, D., Amin, N., Schwikowski, B., & Ideker, T. (2003). Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome Research, 13(11), 2498–2504. https://doi.org/10.1101/GR.1239303

- Szklarczyk, D., Kirsch, R., Koutrouli, M., Nastou, K., Mehryary, F., Hachilif, R., Gable, A. L., Fang, T., Doncheva, N. T., Pyysalo, S., Bork, P., Jensen, L. J., & Von Mering, C. (2023). The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest. Nucleic Acids Research, 51(D1), D638–D646. https://doi.org/10.1093/NAR/GKAC1000

- Zhou, Y., Zhou, B., Pache, L., Chang, M., Khodabakhshi, A. H., Tanaseichuk, O., Benner, C., & Chanda, S. K. (2019). Metascape provides a biologist-oriented resource for the analysis of systems-level datasets. Nature Communications 2019 10:1, 10(1), 1–10. https://doi.org/10.1038/s41467-019-09234-6

---

For questions or issues, please contact the repository maintainer. Refer to the [final course report](Assets/Easwaran_Meena_Journal_paper.pdf) for detailed information and results.

This repository is **strictly for educational purposes and serves as a backup** for my graduate school assignments related to the **BMI 5333: Systems Medicine: Principles and Practice** course at McWilliams School of Biomedical Informatics at UTHealth Houston.

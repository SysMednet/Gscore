[![DOI](https://zenodo.org/badge/678241997.svg)](https://zenodo.org/badge/latestdoi/678241997)

## Gscore
**Gene Set Correlation Enrichment (Gscore)** is an online tool for interpreting and annotating gene expression dataset by using a dataset-derived coexpression network to measure the statistical significance of associations between your interesting differentially expressed gene (DEG) list (i.e., query) and the collection of gene sets. Based on the hypergeometric distribution with Benjamini–Hochberg correction, the tool identifies the coexpressed gene pairs between: 
  
**(1) individual DEG** in the query DEG list and all the DEGs of a certain gene set in the selected collection to determine its association significance;  
**(2) query DEG list** and all the DEGs of a certain gene set in the selected collection to determine its association significance.  
  
This online tool is free and open to all users and there is no login requirement.


## Online Service
The register-free online tool is available on: https://gscore.ibsb.nycu.edu.tw/

## Overview
![Alt text](./readme_img/1.png)  
  
To measure the statistical significance of associations between a list of selected DEGs and a specific gene set (e.g., a group of genes in a KEGG pathway), the Gscore method evaluates the enrichment of coexpressed gene pairs between all the DEGs of the selected list and all the DEGs in the collection of gene sets (e.g., gene sets for 347 KEGG human pathways) based on the hypergeometric distribution with Benjamini–Hochberg correction. The process details are as follows:  
  
**(1).** For each input gene expression dataset with samples belonging to two classes, we first identified the DEGs between control and case samples. Then, we constructed a coexpression network using the expression profiles of the case samples, in which two DEGs with a **Pearson correlation coefficient (|Pearson’s r|) ≥ c** across case samples were considered a coexpressed gene pair. Here, c can be set by the user, for example, to **0.3 (low), 0.5 (moderate), or 0.7 (high)**. Note that genes with identical expression values across all “case” samples are ignored when calculating the correlation.  
**(2).** For each DEG in the query list, Gscore first uses the coexpressed gene pairs between that DEG and all the DEGs of a gene set in the selected collection to determine the association significance for this gene set based on the hypergeometric distribution as follows:  


$$\Large p=1-\sum^{m-1}_{i=0}\frac{\binom{M}{i}\binom{N-M}{n-i}}{\binom{N}{n}}$$  

where *m* and *n* are, respectively, the numbers of coexpressed gene pairs and all possible gene pairs between each DEG in the query DEG list and all the DEGs in a specific gene set. *M* and *N* are, respectively, the total numbers of all the coexpressed gene pairs and all possible gene pairs between each DEG in the query DEG list and all the DEGs in the gene sets of the selected collection. The **FDR *q* value** for multiple hypothesis testing with the Benjamini–Hochberg method was used, and the false discovery rate was controlled at 5%. Here, the association between the DEG in the query DEG list and a certain gene set was considered statistically significant **when its *q* value was ≤ 0.05**.  

**(3).** For the query DEG list, Gscore further measured the statistical significance of association for a specific gene set based on the coexpressed gene pairs between all of the involved DEGs and all the DEGs of this gene set in the selected collection. Then, we computed the *p* value of the **hypergeometric distribution** as follows:

$$\Large p=1-\sum_{i=0}^{m_g-1}\frac{\binom{M_g}{i}\binom{N_g-M_g}{n_g-i}}{\binom{N_g}{n_g}}$$  

where $m_g$ and $n_g$ are, respectively, the numbers of coexpressed gene pairs and all possible gene pairs between all the DEGs of the query list and all the DEGs in a specific gene set. $M_g$ and $N_g$ are, respectively, the total numbers of all the coexpressed gene pairs and all possible gene pairs between all the DEGs of the query list and all the DEGs of gene sets in the selected collection. Here, the association between the query DEG list and a certain gene set was considered statistically significant **when its FDR *q* value was ≤ 0.05** (Benjamini–Hochberg correction4).
  
More detail information is available on the tutorial of the website: https://gscore.ibsb.nycu.edu.tw/tutorial.html


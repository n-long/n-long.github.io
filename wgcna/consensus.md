####Gene co-expression networks with [WGCNA R package](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/)

[Link to code used for analysis](https://github.com/n-long/n-long.github.io/blob/master/wgcna/consensus.r)

Unique counts (number of sequencing reads mapping unambiguously to gene regions) were used to represent gene expression values from 8 normal samples and 8 samples affected with neuromuscular disorder (hereafter referred to as affected). First, a consensus network was constructed to find patterns of co-expression that are unchanged between normal and affected samples. Using all samples as input, pairwise Pearson correlations between gene expression values are calculated with a transformation for negative values (0.5 + 0.5 cor(Xi,Xj), since negative expression values represent downregulation. Then, an adjacency matrix is constructed by raising the correlation matrix to a β power of 18 to approximate scale-free topology. This transformation increases the number of connections in the adjacency matrix, therefore uncovering "hub nodes" which are genes possessing high correlation to many other genes. Hub genes are more likely to involved in regulating a metabolic cascade or pathway. All genes are then assigned a network connectivity value (k) by summing the connection strengths to all other genes. By taking the shortest path between gene neighbors in the adjacency matrix, each gene is assigned to a co-expression module. With the consensus network in place, the above steps are repeated to create a network for only the affected samples. Consensus and affected-specific modules were contrasted by pairwise comparisons of gene membership (Fig. 1). Genes with unique co-expression patterns in the affected condition were tested for functional enrichment using Fisher's Exact Test (FDR < 0.05). In Tables 1-4, "# in test group" shows how many genes belong to each functional category. GO terms associated with DNA breakage and repair are consistently seen across modules, indicating a stress response. However, the intramodular connectivity values were too similar for any meaningful follow-up analysis.

![alt text](https://github.com/n-long/n-long.github.io/blob/master/wgcna/consensus.png "Consensus Network")

Figure 1 – Module membership contrasts for affected (y-axis) against consensus (x-axis) networks. Grey modules are unassigned (too low connectivity). Matrix numbers represent genes shared in both modules, while degree of red is determined by the -log(p) value of the Fisher's exact test for shared membership. Genes that significantly overlap between affected modules and the grey consensus module are inferred to have connectivity unique to the affected condition.

Table 1 - Blue module
![alt text](https://github.com/n-long/n-long.github.io/blob/master/wgcna/royalblue3.png "royalblue3")

Table 2 - Pink module
![alt text](https://github.com/n-long/n-long.github.io/blob/master/wgcna/pink4.png "pink4")

Table 3 - Yellow module
![alt text](https://github.com/n-long/n-long.github.io/blob/master/wgcna/yellow4.png "yellow4")

Table 4 - Green module
![alt text](https://github.com/n-long/n-long.github.io/blob/master/wgcna/darkolivegreen.png "darkolivegreen")

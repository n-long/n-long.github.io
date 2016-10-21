<span style="font-variant: normal">Methods</span>

<span style="font-variant: normal">Gene coexpression networks were
constructed with WGCNA \[1\] from RNA-seq data consisting of 8 normal
samples and 8 affected samples across 4 time points with replicates (day
1 and 3 of pupa stage and day 1 and 3 of adult stage). A signed
consensus network was made by calculating all pairwise correlations
between expression values across all tissue samples in parallel. An
adjacency matrix is constructed by raising the correlations to
a</span><span style="font-variant: normal"> β power of 18 to approximate
scale-free topology. Each gene is then assigned a network connectivity
value (k) by summing the connection strengths to all other genes, and
are assigned to coexpression modules of arbitrary color by similarity in
topological overlap measure, or the shortest path between gene neighbors
in the adjacency matrix. In a consensus network, these modules are
conserved in the normal and affected samples and serve as a reference
for differential network analysis \[2\]. Affected modules were created
by duplicating the parameters used for the consensus network with only
affected samples. Consensus and affected-specific modules were
contrasted by pairwise comparisons of gene membership (Fig. 1).</span>

\

<span class="sd-abs-pos"
style="position: absolute; top: 0in; left: 0in; width: 665px">![](consensus_html_83a7baf986ac76a9.png){width="665"
height="467"} </span><span style="font-variant: normal">Fig. 1 – Module
membership contrasts for affected (y-axis) against consensus (x-axis)
networks. Grey modules are unassigned (too low connectivity). Numbers
represent genes shared in both modules, while degree of red is
determined by the -log(p) value of the Fisher's exact test for shared
membership. Genes that significantly overlap between affected modules
and the grey consensus module are inferred to have affected-specific
connectivity.</span>

\

\

<span style="font-variant: normal">Affected modules that overlap with
unconnected (grey color) consensus modules were examined for
differential expression and GO terminology. Only the parts of the module
that overlap with the grey</span><span style="font-variant: normal">
</span><span style="font-variant: normal">consensus were tested with
Fisher's Exact Test (FDR &lt; 0.05) \[3\] against the entire reference
transcriptome annotations using the Blast2GO plugin \[4\] in CLC
Genomics Workbench 8 \[5\]. Only the reduced terms (i.e. most specific)
are shown for biological process.</span>

\

\

<span class="sd-abs-pos"
style="position: absolute; top: 0in; left: 0in; width: 812px">![](consensus_html_54cadf1319f01ce6.png){width="812"
height="72"} </span><span style="font-variant: normal">Fig.2 –
darkolivegreen</span>

\

<span class="sd-abs-pos"
style="position: absolute; top: 0in; left: 0in; width: 816px">![](consensus_html_6280819c0c941969.png){width="816"
height="99"} </span><span style="font-variant: normal">Fig. 3 –
yellow4</span>

\

<span class="sd-abs-pos"
style="position: absolute; top: 0in; left: 0in; width: 816px">![](consensus_html_a7298bfb839a7803.png){width="816"
height="147"} </span><span style="font-variant: normal">Fig. 3 –
pink4</span>

\

<span class="sd-abs-pos"
style="position: absolute; top: 0in; left: 0in; width: 816px">![](consensus_html_3bb74ac9834d5fbe.png){width="816"
height="217"} </span><span style="font-variant: normal">Fig. 4 –
royalblue3</span>

\

\

<span style="font-variant: normal">GO terms associated with DNA breakage
and repair appear in all modules and is consistent with free radical
damage in mitochondrial dysfunction \[6\].</span>


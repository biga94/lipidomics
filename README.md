# LOEnA
In this repository I've included all the scripts produced to perform an analysis on lipidomics profiles of 16 childrens with different obesity conditions.
At first, I had to reorder names and change the ontology, from LipidBlast ontology to LIPIDMAPS one. Then I've performed a multivariate analysis on this set of data, using a PCA method. 
Following the PCA, I decided to implement also [lipidr](https://www.lipidr.org/) package, that deal specifically with untargeted lipidomics data, in order to correctly parsing and analyze lipids and finally performing and enrichment analysis using both LSEA (a modified LSEA version, included in the package) and [LION/web](https://www.lipidontology.com/).

In a second script, called LOEnA.R, I've included all the lines of code produced to generate a dot plot in GSEA-like style.

Also, I've added my Quarto document (and corresponding pdf) where I've described all the steps performed in Descr.R

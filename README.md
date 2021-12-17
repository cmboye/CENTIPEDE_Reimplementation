# CENTIPEDE_Reimplementation

This is a re-implementation of the original CENTIPEDE package described here : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3044858/. The purpose of the CENTIPEDE package is to predict transcription factor binding based on experimental data and genomic information. We implement an EM algorithm to predict if transcription factors are actually bound. 

Our re-implementation of CENTIPEDE requires a BAM file, which contains read count data, and a peaks file, which contains genomic coordinates that indicate where peaks were called, such as in DNAse-Seq or ATAC-Seq data. The peaks file should also be scanned for the motif of interest using a program such as FIMO. Our pre-processing was done following the tutorial at https://github.com/slowkow/CENTIPEDE.tutorial. Please also see Preprocessing_prep.sh for our code describing the PWM scan we did, as well as the preparation for the optional conservation data. The peaks file is also used to obtain the distance to the nearest transcription start site (TSS). The proximity to the nearest TSS and the score from the motif scan both contribute to our prior probability. Lastly, the BAM file is used to determine the read distribution, which is later used to construct X. 

X is a matrix containing the experimental data. It is constructed such that the first half is per-base coverage from the forward strand and the second half is per-base coverage from the reverse strand. The matrix is S x L dimensions such that L is the number of motif loci from the PWM scan and S = (200 + length of motif)x2. Thus, to construct X, you need to get the per-base coverage for the motif site as well as 100bp flanking each side of the motif. We used the centipede_data() function from the tutorial at https://github.com/slowkow/CENTIPEDE.tutorial to construct our X, but this can be done in other ways as well. 

The function run_centipede() is used to run the model. It requires an input list, which is a list containing containing PWMScan, TSS data, and X. If using conservation data, this matrix must be included as a seperate argument. 

To plot the results, we provide a function called plot_cutsites(). This plots the lambda parameter per each locus, which makes it easy to visualize footprints. We also include the functions compare_results() and plot_roc2(), which we use to validate our data when compared to ENCODE ChIP-Seq data. compare_results() requires a binary vector of bound status predictions, a dataframe of coordinates corresponding to results with columns 'chr','start','end' and rows as the coordinates in "chr:start-end" format, and a bed file of ChIP-seq peaks. The output is a dataframe including precision, FDR, FPR, sensitivity, and a dataframe comparing the binding predictions with the ChIP-Seq data. This dataframe can be used as an input for plot_roc2(), which uses PRROC to plot an ROC and provide the AUC value. 

# Installation

```r
#install dependencies, these packages will prevent installation if not pre-installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("AnnotationDbi")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

#install package
install.packages("devtools")
devtools::install_github("cmboye/CENTIPEDE_Reimplementation")
```

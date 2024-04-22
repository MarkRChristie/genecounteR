Please visit this link for a [vignette](http://htmlpreview.github.io/?https://github.com/MarkRChristie/genecounteR/blob/main/vignettes/Vignette.html)
The vignette can also be found in the vignettes directory and is pasted below for convinience.

## Before starting

We recommend that users follow this vignette with the example data
before using on their own data sets. Feel free to copy and paste the
below code into a new R script when running on your own data. There are
several steps to follow which work best when followed sequentially.
Please install the libraries ‘devtools’ (for library installation) and
‘scales’ (for transparency on plots) before proceeding.

## Install the library

    devtools::install_github("MarkRChristie/genecounteR")

## Load the library

    library('genecounteR')

## Set working directory

Because you will be working with a medium sized dataset, we suggest
setting a working directory so you can later delete the file.

You can set the working directory directly or use the R package ‘here’.
For the example we will just use the Downloads directory; in practice
you may want to create a specific directory:

    setwd("~/Downloads/") # use something like this for Linux and Mac operating systems
    setwd("C:/Users/username/Downloads") # use something like this for Windows machines, where username should be substituted with your actual user name.

or you can use the pacakge ‘here’ to create relative directories. Learn
more [here](https://cloud.r-project.org/web/packages/here/index.html)

## Load the test dataset

A test dataset called “TestData.txt” with 15 samples and 120,076 “genes”
can be downloaded directly from the [genecounter github
repo](https://github.com/MarkRChristie/genecounteR/tree/main/inst/extdata)

Alternatively, you can download the file directly from Dropbox into your
current working directory:

    download.file('https://www.dropbox.com/scl/fi/zcx4ou8wahov9wivw8hxq/TestData.txt?rlkey=jdzk0qlpixsw13y3thec1zhmr&st=hcx6duqv&dl=1',
                  destfile="TestData.txt",
                  method="auto") # will fail if a file named TestData.txt already exists in your working directory

Next, we load the file into R and perform some minimal formatting. After
formatting, the first column should be a gene name or gene id, the
second column the gene lengths, and each remaining column the gene
counts per sample (one column per sample). Each row should be a
gene/feature/isoform. The test dataset comes from a standard
featurecounts pipeline.

    gene.counts <- read.table("TestData.txt", sep="\t", header=T, row.names=1, check.names=F)
    colnames(gene.counts)

    dat <- gene.counts[, -c(1:4)]  # remove all columns except gene lengths and gene counts
    dat <- cbind(1:nrow(dat), dat) # add unique gene ids; can use rownames if you want to preserve name

When importing your own dataset as “dat” make sure that each row is a
gene and that  
the columns are formatted as (col1= gene id; col2 = gene length; col3 =
sample 1 counts; col4 = sample 2 counts etc.) Run head(gene.counts) and
head(dat) after running the above code to understand required data
format.

## Zero Remover

Remove all genes that have all zeros for counts; 2 plots created

    percent.zeros.retained <- 0.8 # percentage of samples allowed to have a count of zero; applied for each gene (row) separately
    dat <- ZeroRemover(dat, percent.zeros.retained)

## Sample Standardization

Standardize gene counts by the total number of reads per sample; 1 plot
created

    dat <- SampleStandardization(dat)

## Gene Standardization

Standardize by residuals between gene length and total expression;
linear by default; creates 4 plots

    dat  <- GeneStandardization(dat)
    ref2 <- dat # for Matching Back after analyses are completed 

Users can import their own standardized data sets (using different
methods) at this step if so desired.

## Visualize Expression

Retain top x percent of genes that equal or exceed a log mean gene
count, where x is set with “threshold”; illustrated on the single plot
that is generated.

    threshold <- 0.15 # will retain top 15% of genes by log mean expression
    dat       <- GeneCountDistributionPlot(dat, threshold) 

## Rank Order Expression

Calculate rank-ordered expression within each gene; 3 plots generated

    deltas <- RankOrder(dat) 

## Boot Strap

Calculate significance of each gene; 1 plot generated, significant genes
in red

    alpha  <- 0.05
    deltas <- BootStrap(dat, deltas, alpha)
    deltas[, c(1:8, ncol(deltas))] # visual results

## Match Back

Match back to original data set, plot top genes, report gene IDs

    MatchBack(deltas, dat, ref2)

## Run All

If you would like to run all of the steps at once, you can wrap them
into a single function called RunAll and then simply call that function.
Generally not advised, but included below are all the R commands run
above without many comments, which may facilitate creating code that
works best for you!

Step 1: create function

    RunAll <- function(){
      library('genecounteR')
      setwd("~/Downloads/") 
      gene.counts <- read.table("TestData.txt", sep="\t", header=T, row.names=1, check.names=F)
      colnames(gene.counts)
      
      dat  <- gene.counts[, -c(1:4)]  # remove all columns except gene lengths and gene counts
      dat  <- cbind(1:nrow(dat), dat)
      percent.zeros.retained <- 0.8 # percentage of samples allowed to have a count of zero; applied for each gene (row) seperately
      dat  <- ZeroRemover(dat, percent.zeros.retained)
      dat  <- SampleStandardization(dat)
      dat  <- GeneStandardization(dat)
      ref2 <- dat # for Matching Back after analyses are completed 
      threshold <- 0.15
      dat    <- GeneCountDistributionPlot(dat, threshold) 
      deltas <- RankOrder(dat) 
      alpha  <- 0.05
      deltas <- BootStrap(dat, deltas, alpha)
      deltas[, c(1,3,4, ncol(deltas))] # see results table
      MatchBack(deltas, dat, ref2)
    }

Step 2: run function

    RunAll()

## Notes
**Program created by Mark Christie;
Created in R 4.4.3  
genecounteR: a simple framework for identifying constituitively expressed genes in large RNA-Seq data sets**




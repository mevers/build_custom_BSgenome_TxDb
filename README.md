# Build BSgenome and TxDb R libraries

This respository contains R scripts and a YAML config file to build

1. a custom `BSgenome` R library, and
2. a custom `TxDb` R library

based on the [Ensembl GRCh38 genome assembly/annotation](ftp://ftp.ensembl.org/pub/release-89) and the [GenBank U13369.1 canonical rDNA sequence/annotation](https://www.ncbi.nlm.nih.gov/nuccore/555853).


## Requirements

The following R libraries are required.

1. From CRAN:

    * `devtools`
    * `optparse`
    * `tidyverse`
    * `yaml`

2. From Bioconductor:

    * `BSgenome`
    * `GenomicFeatures`


You can install CRAN packages with e.g. `install.packages("tidyverse");`, and packages from Bioconductor with e.g. `source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome")`.


## How to build the R packages

**Note that this needs to be done only once. Provided building/installing was successful, the packages will be available within the R environment as libraries.**

1. Clone this respository

    ```
    git clone https://github.com/mevers/build_custom_BSgenome_TxDb
    ```

2. Change to folder that constains the scripts and config file.

    ```
    cd build_custom_BSgenome_TxDb
    ```

3. To build the R packages, there are two options.

    1. `make_packages.sh` is a convenience bash script that builds both the `BSgenome` and `TxDb` R packages. Run the script with:

        ```
        sh make_packages.sh
        ```

        Note that this will only work on Linux/Unix systems.

    2. You can build the R packages separately; `make_BSgenome_package.R` and `make_TxDb_package.R` are two R command-line scripts that build the `BSgenome` and `TxDb` packages, respectively. Run the scripts with:

        ```
        Rscript make_BSgenome_package.R -i config.yaml
        ```

        and

        ```
        Rscript make_TxDb_package.R -i config.yaml
        ```

        On a Windows systems, you might have to call `Rscript` with its full path at the Windows command prompt, e.g. `C:\R-3.4.3\bin\x64\Rscript.exe` (depending on where R is installed).

4. After the packages have been successfully build you should see two files `BSgenome.Hsapiens.mevers.hs1_1.0.0.tar.gz` and `TxDb.Hsapiens.mevers.hs1_1.0.0.tar.gz`. Install both packages in R in the usual way, e.g.

    ```
    install.packages("BSgenome.Hsapiens.mevers.hs1_1.0.0.tar.gz", repos = NULL, type = "source");
    install.packages("TxDb.Hsapiens.mevers.hs1_1.0.0.tar.gz", repos = NULL, type = "source");
    ```


## How to use the libraries

Both the `BSgenome` and `TxDb` libraries can now be loaded within R.

```
library("BSgenome.Hsapiens.mevers.hs1");
library("TxDb.Hsapiens.mevers.hs1");
```

Get a general genome sequence summary:
```
BSgenome.Hsapiens.mevers.hs1;
#Human genome:
## organism: Homo sapiens (Human)
## provider: Ensembl+NCBI
## provider version: hs1
## release date: Mar. 2018
## release name: Human masked chromosome plus MT plus rDNA sequences
## 26 sequences:
##   chr1     chr2     chr3     chr4     chr5     chr6     chr7     chr8
##   chr9     chr10    chr11    chr12    chr13    chr14    chr15    chr16
##   chr17    chr18    chr19    chr20    chr21    chr22    chrX     chrY
##   chrM     U13369.1
## (use 'seqnames()' to see all the sequence names, use the '$' or '[[' operator
## to access a given sequence)
```

Extract lengths and features of individual chromosome sequences:
```
seqinfo(BSgenome.Hsapiens.mevers.hs1);
#Seqinfo object with 26 sequences (2 circular) from hs1 genome:
#  seqnames seqlengths isCircular genome
#  chr1      248956422      FALSE    hs1
#  chr2      242193529      FALSE    hs1
#  chr3      198295559      FALSE    hs1
#  chr4      190214555      FALSE    hs1
#  chr5      181538259      FALSE    hs1
#  ...             ...        ...    ...
#  chr22      50818468      FALSE    hs1
#  chrX      156040895      FALSE    hs1
#  chrY       57227415      FALSE    hs1
#  chrM          16569       TRUE    hs1
#  U13369.1      42999       TRUE    hs1
```

Extract all U13369.1 transcripts:

```
txdb <- TxDb.Hsapiens.mevers.hs1;
gr <- transcripts(txdb, columns = c("gene_id", "tx_id", "tx_name"));
gr[seqnames(gr) == "U13369.1"];
#GRanges object with 1 range and 3 metadata columns:
#      seqnames     ranges strand |         gene_id     tx_id     tx_name
#         <Rle>  <IRanges>  <Rle> | <CharacterList> <integer> <character>
#  [1] U13369.1 [1, 13314]      + |            rDNA    199168    pre-rRNA
#  -------
#  seqinfo: 26 sequences (2 circular) from an unspecified genome; no seqlengths    
```

Or extract all U13369.1 "exons":

```
gr <- exons(txdb, columns = c("exon_id", "exon_name", "tx_name", "gene_id"));
gr[seqnames(gr) == "U13369.1"];
#GRanges object with 7 ranges and 4 metadata columns:
#      seqnames         ranges strand |   exon_id   exon_name         tx_name
#         <Rle>      <IRanges>  <Rle> | <integer> <character> <CharacterList>
#  [1] U13369.1 [    1,  3656]      + |    680542       5'ETS        pre-rRNA
#  [2] U13369.1 [ 3657,  5527]      + |    680543    18S-rRNA        pre-rRNA
#  [3] U13369.1 [ 5528,  6622]      + |    680544        ITS1        pre-rRNA
#  [4] U13369.1 [ 6623,  6779]      + |    680545   5.8S-rRNA        pre-rRNA
#  [5] U13369.1 [ 6780,  7934]      + |    680546        ITS2        pre-rRNA
#  [6] U13369.1 [ 7935, 12969]      + |    680547    28S-rRNA        pre-rRNA
#  [7] U13369.1 [12970, 13314]      + |    680548       3'ETS        pre-rRNA
#              gene_id
#      <CharacterList>
#  [1]            rDNA
#  [2]            rDNA
#  [3]            rDNA
#  [4]            rDNA
#  [5]            rDNA
#  [6]            rDNA
#  [7]            rDNA
#  -------
#  seqinfo: 26 sequences (2 circular) from an unspecified genome; no seqlengths
```


## References and details

For more details, see e.g.
* [Efficient genome searching with Biostrings and the BSgenome data packages](https://bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/GenomeSearching.pdf)
* [An overview of the Biostrings/BSgenome framework](https://www.bioconductor.org/help/course-materials/2011/BioC2011/LabStuff/BiostringsBSgenomeOverview.pdf)
* [Making and utilizing TxDb Objects](https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf)

Ensembl data are based on the

* GRCh38 genome assembly, release 89
* the GRCh38.89 GTF annotation file.


## Author

In case of questions please contact [Maurits Evers](mailto:maurits.evers@anu.edu.au).

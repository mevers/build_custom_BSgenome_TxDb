#!/usr/bin/env Rscript


# The R command line script
#  (1) downloads GRCh38 genome annotation from Ensembl,
#  (2) downloads the U13369.1 rDNA annotation,
#  (3) builds a TxDb object, and
#  (4) converts the TxDb object into an R library.
#
# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Date: 09/04/2018


## ------------------------------------------------------------------------
## Start clock
t0 <- Sys.time();


## ------------------------------------------------------------------------
## Load libraries
suppressWarnings(suppressMessages(library("optparse")));         # Python-style command line args
suppressWarnings(suppressMessages(library("yaml")));             # Parse YAML
suppressWarnings(suppressMessages(library("GenomicFeatures")));  # for TxDB
suppressWarnings(suppressMessages(library("devtools")));         # for build()
suppressWarnings(suppressMessages(library("tidyverse")));        # for tidying data


## ------------------------------------------------------------------------
## Parse command line arguments
option_list <- list(
    make_option(
        c("-i", "--input"),
        type = "character",
        default = NULL,
        help = "YAML config file",
        metavar = "character"),
#    make_option(
#        c("--annotdir"),
#        type = "character",
#        default = "annot",
#        help = "Folder where annotation files are stored [default %default]",
#        metavar = "character"),
    make_option(
        c("-f", "--forceall"),
        type = "logical",
        action = "store_true",
        default = FALSE,
        help = "Force download (this will overwrite existing files) [default %default]",
        metavar = "character")
);
opt_parser <- OptionParser(option_list = option_list);
args <- parse_args(opt_parser);
if (is.null(args$input)) {
    print_help(opt_parser);
    stop("Must give YAML config file.\n", call. = FALSE);
}


## ------------------------------------------------------------------------
## Custom function to generate timestamp
ts <- function() {
    return(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"));
}


## ------------------------------------------------------------------------
# Custom function to download file (if condition is met)
cond_download <- function(url, fn, force, id) {
    if (force | !file.exists(fn)) {
        download.file(url, fn, quiet = TRUE);
        cat(sprintf(
            "%s Finished downloading %s file %s.\n", ts(), id, fn));
    } else {
        cat(sprintf(
            "%s File %s already exists. Skipping (use -f to force download).\n",
            ts(),
            fn));
    }
}


## ------------------------------------------------------------------------
## Global variables
input <- args$input;
#annotdir <- args$annotdir;
annotdir <- "annot";
forceall <- args$forceall;
cat(sprintf("%s Parameter summary\n", ts()));
cat(sprintf(" input          = %s\n", input));
#cat(sprintf(" annotdir       = %s\n", annotdir));
cat(sprintf(" forceall       = %s\n", forceall));


## ------------------------------------------------------------------------
# Check if input files and output directory exists
if (!file.exists(input)) {
    stop(
        sprintf("Input file %s does not exists.\n", input),
        call. = FALSE);
}
if (!dir.exists(annotdir)) {
    stop(
        sprintf("Folder %s does not exists.\n", annotdir),
        call. = FALSE);
}


## ------------------------------------------------------------------------
# Parse YAML config file
cfg <- read_yaml(input);


## ------------------------------------------------------------------------
# Create folders for sequence and annotation files (if necessary)
#if (!dir.exists(annotdir)) {
#    cat(sprintf("%s Creating folder %s\n", ts(), annotdir));
#    dir.create(annotdir);
#}


## ------------------------------------------------------------------------
# Download genome gtf.gz annotation from Ensembl
cat(sprintf("%s Downloading Ensembl genome annotation file...\n", ts()));
url <- paste0(
    cfg$download$baseurl_ensembl,
    "/gtf/homo_sapiens/Homo_sapiens.GRCh38.89.gtf.gz")
fn1 <- paste0(annotdir, "/Homo_sapiens.GRCh38.89.gtf.gz");
cond_download(url, fn1, forceall, "Ensembl annotation")


## ------------------------------------------------------------------------
# Download rDNA gtf annotation from Dropbox
cat(sprintf("%s Downloading rDNA annotation file...\n", ts()));
url <- "https://www.dropbox.com/s/jnvwo8avkcc12iq/U13369.1.gtf.gz?dl=1";
fn2 <- paste0(annotdir, "/U13369.1.gtf.gz");
cond_download(url, fn2, forceall, "rDNA annotation");


## ------------------------------------------------------------------------
# Clean annotation files, merge and save in single file
cat(sprintf("%s Merging annotation files...\n", ts()));
fn <- paste0(annotdir, "/GRCh38.89+U13369.1.gtf");
if (forceall | !file.exists(fn)) {
    df1 <- read_delim(
        fn1,
        delim = "\t",
        comment = "#",
        col_names = FALSE,
        col_types = paste0(rep("c", 9), collapse = ""),
        escape_double = FALSE);
    df2 <- read_delim(
        fn2,
        delim = "\t",
        comment = "#",
        col_names = FALSE,
        col_types = paste0(rep("c", 9), collapse = ""),
        escape_double = FALSE);
    df <- df1 %>%
        filter(X1 %in% c(seq(1:22), "MT", "X", "Y")) %>%
        bind_rows(df2);
    cat(sprintf("%s Writing merged annotation files...\n", ts()));
    # Can't use write_delim because it puts double-quotes around character
    # entries, and this causes an error in makeTxDbFromGFF
    write.table(df, fn, sep = "\t", col.names = F, row.names = F, quote = F);
} else {
    cat(sprintf(
        "%s File %s already exists. Skipping (use -f to force build).\n",
        ts(),
        fn));
}


## ------------------------------------------------------------------------
# Make TxDb package
cat(sprintf("%s Preparing TxDb package...\n", ts()));
pkgname <- gsub("BSgenome", "TxDb", cfg$BSgenome$Package);
if (forceall | !dir.exists(pkgname)) {
    if (forceall & dir.exists(pkgname)) {
        cat(sprintf("%s Deleting existing folder %s...\n",
        ts(),
        pkgname));
        unlink(pkgname, recursive = TRUE);
    }
    cat(sprintf("%s Creating TxDb object from annotation file...\n", ts()));
    metadata <- data.frame(
        name = c("Resource URL"),
        value = c(""),
        stringsAsFactors = FALSE);
    txdb <- makeTxDbFromGFF(
        file = fn,
        format = "gtf",
        dataSource = cfg$BSgenome$provider,
        organism = cfg$BSgenome$organism,
        circ_seqs = eval(parse(text = cfg$BSgenome$circ_seqs)),
        metadata = metadata);
    cat(sprintf(
        "%s Creating package source files in folder %s\n",
        ts(),
        pkgname));
    makeTxDbPackage(
        txdb,
        version = cfg$BSgenome$Version,
        maintainer = cfg$BSgenome$Maintainer,
        author = cfg$BSgenome$Author,
        license = cfg$BSgenome$License,
        pkgname = pkgname,
        provider = cfg$BSgenome$provider,
        providerVersion =cfg$BSgenome$provider_version);
} else {
    cat(sprintf(
        "%s Folder %s already exists. Skipping (use -f to force rebuild).
         %12s Or manually delete target folder (e.g. 'rm -rf %s').\n",
        ts(),
        pkgname,
        "",
        pkgname));
}


### ------------------------------------------------------------------------
## Build TxDb package
cat(sprintf("%s Building TxDb package...\n", ts()));
fn <- paste0(
    gsub("BSgenome", "TxDb", cfg$BSgenome$Package),
    "_",
    cfg$BSgenome$Version,
    ".tar.gz");
if (forceall | !file.exists(fn)) {
    fn <- devtools::build(gsub("BSgenome", "TxDb", cfg$BSgenome$Package));
    cat(sprintf("%s TxDb package file %s created.\n",
                 ts(),
                 fn));
} else {
    cat(sprintf(
        "%s File %s already exists. Skipping (use -f to force rebuild).
         %12s Or manually delete target file (e.g. 'rm -f %s').\n",
        ts(),
        fn,
        "",
        fn));
}


### ------------------------------------------------------------------------
## Done
cat(sprintf("%s All done.\n", ts()));
cat(sprintf(
    "%s Install package with 'install.packages(\"%s\")'.\n",
    ts(),
    fn))

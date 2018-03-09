library(TxDb.Hsapiens.mevers.hs1);
txdb <- TxDb.Hsapiens.mevers.hs1;


seqlevels(txdb);
#[1] "1"        "2"        "3"        "4"        "5"        "6"
#[7] "7"        "8"        "9"        "10"       "11"       "12"
#[13] "13"       "14"       "15"       "16"       "17"       "18"
#[19] "19"       "20"       "21"       "22"       "X"        "Y"
#[25] "MT"       "U13369.1"


gr <- transcripts(txdb, columns = c("gene_id", "tx_id", "tx_name"));
gr[seqnames(gr) == "U13369.1"];
#GRanges object with 1 range and 3 metadata columns:
#      seqnames     ranges strand |         gene_id     tx_id     tx_name
#         <Rle>  <IRanges>  <Rle> | <CharacterList> <integer> <character>
#  [1] U13369.1 [1, 13314]      + |            rDNA    199168    pre-rRNA
#  -------
#  seqinfo: 26 sequences (2 circular) from an unspecified genome; no seqlengths


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


promoters <- promoters(txdb, upstream = 2000, downstream = 2000);
promoters[seqnames(promoters) == "U13369.1"];
#GRanges object with 1 range and 2 metadata columns:
#      seqnames        ranges strand |     tx_id     tx_name
#         <Rle>     <IRanges>  <Rle> | <integer> <character>
#  [1] U13369.1 [-1999, 2000]      + |    199168    pre-rRNA
#  -------
#  seqinfo: 26 sequences (2 circular) from an unspecified genome; no seqlengths


#library(BSgenome.Hsapiens.mevers.hs1);
#genome <- BSgenome.Hsapiens.mevers.hs1;

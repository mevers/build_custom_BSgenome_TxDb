library(TxDb.Hsapiens.mevers.hs1);
txdb <- TxDb.Hsapiens.mevers.hs1;

seqlevels(txdb);

gr <- transcripts(txdb, columns = c("gene_id", "tx_id", "tx_name"));
gr[seqnames(gr) == "U13369.1"];
#GRanges object with 8 ranges and 3 metadata columns:
#      seqnames         ranges strand |         gene_id     tx_id     tx_name
#         <Rle>      <IRanges>  <Rle> | <CharacterList> <integer> <character>
#  [1] U13369.1 [    1,  3656]      + |            rDNA    199168       5'ETS
#  [2] U13369.1 [    1, 13314]      + |            rDNA    199169    pre-rRNA
#  [3] U13369.1 [ 3657,  5527]      + |            rDNA    199170    18S-rRNA
#  [4] U13369.1 [ 5528,  6622]      + |            rDNA    199171        ITS1
#  [5] U13369.1 [ 6623,  6779]      + |            rDNA    199172   5.8S-rRNA
#  [6] U13369.1 [ 6780,  7934]      + |            rDNA    199173        ITS2
#  [7] U13369.1 [ 7935, 12969]      + |            rDNA    199174    28S-rRNA
#  [8] U13369.1 [12970, 13314]      + |            rDNA    199175       3'ETS
#  -------
#  seqinfo: 26 sequences (1 circular) from an unspecified genome; no seqlengths


gr <- exons(txdb, columns = c("exon_id", "exon_name", "tx_name", "gene_id"));
gr[seqnames(gr) == "U13369.1"];

promoters <- promoters(txdb, upstream = 2000, downstream = 2000);
promoters[seqnames(promoters) == "U13369.1"];

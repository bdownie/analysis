library(ggbio)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(biovizBase)
library(Rsamtools)
data("genesymbol", package = "biovizBase")


# Which genes to plot
geneList <- c("TET3","TTC34","ADAMTS14")
geneList <- c("KSR1P1")
annotation_file <- "/ngsvol/Genomes/Homo_sapiens/ens37.75/Homo_sapiens.GRCh37.75.gtf"

egrep.text <- paste(geneList,collapse="|")
egrep.text <- paste("egrep \"",egrep.text,"\" ",annotation_file,sep = "")

df.annotation <- read.table(pipe(egrep.text),sep="\t")

names(df.annotation) <- c("chr","source","feature","start","end","score","strand","frame","attribute")
df.annotation$gene_name = df.annotation$attribute
df.annotation$gene_id = df.annotation$attribute
df.annotation$exon_id = df.annotation$attribute
df.annotation <- df.annotation[-9]

df.annotation$gene_id <- sub(".*gene_id ","",df.annotation$gene_id,perl=TRUE)
df.annotation$gene_id <- sub(";.*","",df.annotation$gene_id,perl=TRUE)

df.annotation$gene_name <- sub(".*gene_name ","",df.annotation$gene_name,perl=TRUE)
df.annotation$gene_name <- sub(";.*","",df.annotation$gene_name,perl=TRUE)

df.annotation$exon_id <- sub(".*exon_id ","",df.annotation$exon_id,perl=TRUE)
df.annotation$exon_id <- sub(";.*","",df.annotation$exon_id,perl=TRUE)



txdb <- makeTranscriptDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="useast.ensembl.org")
#txdb <- makeTranscriptDbFromGFF(file=annotation_file,format="gtf")
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#makeTranscriptDbFromBiomart
# From http://davetang.org/muse/2013/01/02/iranges-and-genomicranges/
annotation.gene <- with(df.annotation[df.annotation$feature == "gene",], GRanges(chr, IRanges(start, end), strand, score, id=gene_id))
annotation.exon <- with(df.annotation[df.annotation$feature == "exon",], GRanges(chr, IRanges(start, end), strand, score, id=exon_id))

##
#

# Get ideogram for hg19
#p.ideo <- Ideogram(genome = "hg19")

# Limit
#p.ideo + xlim(GRanges("chr2", IRanges(1e8, 1e8+10000000)))

bfl <- c("bam/Sample_r827sCtr1_hs_chip_sr_Zeisberg_C_1.bam","bam/Sample_r827sCtr2_hs_chip_sr_Zeisberg_C_2.bam","bam/Sample_r827sTreat1_hs_chip_sr_Zeisberg_B_1.ds.bam","bam/Sample_r827sTreat2_hs_chip_sr_Zeisberg_B_2.ds.bam")
bfl <- BamFileList(bfl)



#bamfile <- BamFile(file="/ngsdata/intern/141126_SN933_0172_AC5J5GACXX/Aligned/Zeisberg_827/bam/Sample_r827sCtr1_hs_chip_sr_Zeisberg_C_1.bam")
#ctr1.bam <- BamFile(file="/ngsdata/intern/141126_SN933_0172_AC5J5GACXX/Aligned/Zeisberg_827/bam/Sample_r827sCtr1_hs_chip_sr_Zeisberg_C_1.bam")
#wh <- keepSeqlevels(txdb, "chr1")
#autoplot(ctr1.bam, chr = "chr21")

for (i in 1:length(geneList)) {
	n <- geneList[i]


	p <- ScanBamParam(which = genesymbol[n])

	ga1 <- readGAlignmentsFromBam(bfl[[1]], param=p,use.names=TRUE)
	ga2 <- readGAlignmentsFromBam(bfl[[2]], param=p,use.names=TRUE)
	ga3 <- readGAlignmentsFromBam(bfl[[3]], param=p,use.names=TRUE)
	ga4 <- readGAlignmentsFromBam(bfl[[4]], param=p,use.names=TRUE)
	p1 <- autoplot(ga1, geom = "line", stat = "coverage")
	p2 <- autoplot(ga2, geom = "line", stat = "coverage")
	p3 <- autoplot(ga3, geom = "line", stat = "coverage")
	p4 <- autoplot(ga4, geom = "line", stat = "coverage")

	

	wh = genesymbol[geneList[i]]
	wh <- range(wh, ignore.strand = TRUE)

	genes <- autoplot(Homo.sapiens, which = wh, gap.geom = "chevron")

	tracks(Control1=p1,Control2=p2,Treatment1=p3,Treatment2=p4,title=paste("Coverage of",n)) + ylim(0,100) + tracks(gene=genes)
}

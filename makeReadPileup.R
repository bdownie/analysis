library(ggbio)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(biovizBase)
library(Rsamtools)
library(org.Hs.eg.db)
library (TxDb.Hsapiens.UCSC.hg19.knownGene)
library(biomaRt)

ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="feb2014.archive.ensembl.org",path="/biomart/martservice",archive=FALSE)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)

gene_ids <- c(
	"ENSG00000186205",
	"ENSG00000133636",
	"ENSG00000105723",
	"ENSG00000148408",
	"ENSG00000131471",
	"ENSG00000092421"
	)



genes.annotation=getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id","external_gene_id"),filters=c("ensembl_gene_id"),values=gene_ids, mart= ensembl) # fuction to get  gene id's and gene name from data base
transcripts = genes.annotation$ensembl_transcript_id

genes2id <- unique(genes.annotation[,c(1,3)])
rownames(genes2id) <- genes2id[,1] 



txdb=makeTxDbFromBiomart(transcript_ids=transcripts,biomart="ENSEMBL_MART_ENSEMBL",host="feb2014.archive.ensembl.org")

genes <-        genes(txdb)
seqlevels(genes) <- sapply(seqlevels(genes),function(x) { if (any(grep("^\\d",x))) { paste("chr",x,sep="") } else { x } })
exons <-        exonsBy(txdb, "gene")

#wh.genes <- genes(txdb,vals=list(gene_id=gene_ids))
wh.transcripts <- transcriptsBy(txdb ,by="gene")


bfl <- c(
	"/ngsdata/intern/150409_SN933_0186_AC6TN4ACXX/Aligned/Boemecke_909/bam/DiseaseUntreated.bam",
	"/ngsdata/intern/150409_SN933_0186_AC6TN4ACXX/Aligned/Boemecke_909/bam/DiseaseEpi.bam",
	"/ngsdata/intern/150409_SN933_0186_AC6TN4ACXX/Aligned/Boemecke_909/bam/ControlUntreated.bam",
	"/ngsdata/intern/150409_SN933_0186_AC6TN4ACXX/Aligned/Boemecke_909/bam/ControlEpi.bam"
	)
bfn <- c(
	"Disease Untreated",
	"Disease Epi",
	"Control Untreated",
	"Control Epi"
	)

# Normalize sizes
norm <- c(
	135765722,
	134176076,
	89563374,
	85207016
)
norm <- min(norm)/norm

pdf("genes.pdf")
for (i in 1:length(gene_ids)) {
	n <- gene_ids[i]

	#p <- ScanBamParam(which = genesymbol[n])
	p <- ScanBamParam(which = genes[n])

	wh <- wh.transcripts[n][[1]]
	z <- tracks(autoplot(txdb, which = wh, gap.geom = "chevron"),xlim=reduce(ranges(wh)),scale.height=2)
	seqlevels(wh) <- sapply(seqlevels(wh),function(x) { if (any(grep("^\\d",x))) { paste("chr",x,sep="") } else { x } })
	for (j in 1:length(bfl)) {
		ga <- readGAlignments(bfl[j], param=p,use.names=TRUE)
		ga <- ga[sample(length(ga),length(ga) * norm[j])]

		ga2 <- ga[(start(ga) > start(reduce(ranges(wh)))) & (end(ga) < end(reduce(ranges(wh))))]

		x <- tracks(autoplot(ga2,  stat = "coverage",ylab=genes2id[n,2]))
		if (j == 1) { 
			y <- x
		} else { 
			y <- y + x
		}
	}
	print(y + z + theme_alignment()  )
}
dev.off()



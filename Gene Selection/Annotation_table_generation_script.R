library(biomaRt)


Species = 'Human'
if (Species == "Mouse") {
  ensembl = useEnsembl(biomart = "ensembl", 
                       dataset = "mmusculus_gene_ensembl", 
                       mirror = "useast")
}

if (Species == "Human") {
  ensembl = useEnsembl(biomart = "ensembl", 
                       dataset = "hsapiens_gene_ensembl", 
                       mirror = "useast")
}
cat("Biomart object created ! \n")
#Donwloading the information about the gene of interest 

gene_symbol= Parameters["Gene",]

cat(paste("Species studied :",Species,"\n"))
cat(paste("Gene studied :",gene_symbol,"\n"))

if (Species == "Mouse") {
  gene_sequence=getBM(mart = ensembl,attributes = c("mgi_symbol","transcript_length","percentage_gene_gc_content"))
  gene_sequence = gene_sequence[gene_sequence$mgi_symbol!="",]
  Aggregated_length = aggregate(gene_sequence$transcript_length,by=list(gene_sequence$mgi_symbol),FUN = mean)
  Aggregated_GC_content = aggregate(gene_sequence$percentage_gene_gc_content,by=list(gene_sequence$mgi_symbol),FUN = mean)
  Merged_table = data.frame(Gene= Aggregated_length$Group.1,GC_content = Aggregated_GC_content$x,Mean_transcript_length = Aggregated_length$x)
  write.table(Merged_table,"~/Gene_annotation_mouse.txt",sep="\t")
  
  }

#To modify : different attribute name for mouse and human gene symbol : hugo ?
if (Species == "Human") {
  gene_sequence=getBM(mart = ensembl,attributes = c("hgnc_symbol","transcript_length","percentage_gene_gc_content"))
  gene_sequence = gene_sequence[gene_sequence$hgnc_symbol!="",]
  Aggregated_length = aggregate(gene_sequence$transcript_length,by=list(gene_sequence$hgnc_symbol),FUN = mean)
  Aggregated_GC_content = aggregate(gene_sequence$percentage_gene_gc_content,by=list(gene_sequence$hgnc_symbol),FUN = mean)
  Merged_table = data.frame(Gene= Aggregated_length$Group.1,GC_content = Aggregated_GC_content$x,Mean_transcript_length = Aggregated_length$x)
  write.table(Merged_table,"~/Gene_annotation_human.txt",sep="\t")
}
cat("Sequences information successfuly downloaded. \n")

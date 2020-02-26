suppressMessages(library("biomaRt"))
suppressMessages(library("Biostrings"))
cat("Library loaded ! \n")
args <- commandArgs(trailingOnly = TRUE)

#Loading the text file containing all parameter values

Parameter_file = as.character(args[1])
Parameters = read.delim(Parameter_file,header = F,row.names = 1)
colnames(Parameters) = "Value"
Parameters$Value = as.character(Parameters$Value)

Output_directory = Parameters["Output_directory",]
Transcriptome_index = Parameters["Transcriptome_index",]

#Checking ouptut directory exist, if not create it
if (!dir.exists(Output_directory)) {
  cat("Output directory does not exist. Creating it... \n")
  dir.create(Output_directory)
}


#Selecting the species of interest (Mouse or human only)
Species = Parameters["Species",]

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
  gene_sequence=getBM(mart = ensembl,attributes = c("coding","mgi_symbol","ensembl_transcript_id"),
                      filters = "mgi_symbol",values = gene_symbol) 
}

#To modify : different attribute name for mouse and human gene symbol : hugo ?
if (Species == "Human") {
  gene_sequence=getBM(mart = ensembl,attributes = c("coding","hgnc_symbol","ensembl_transcript_id"),
                      filters = "hgnc_symbol",values = gene_symbol) 
}
cat("Sequences information successfuly downloaded. \n")


gene_sequence$Length=nchar(gene_sequence[,1])
gene_sequence = gene_sequence[gene_sequence$coding!="Sequence unavailable",]
gene_sequence = gene_sequence[order(gene_sequence$Length,decreasing = T),]
rownames(gene_sequence) = 1:nrow(gene_sequence)
#Checking step : how many sequences found ?
#If more than one : which one to use ?

if (nrow(gene_sequence)==0) {
  stop("No transcript sequence found for that gene. Please check the spelling of the gene !")
}

if (nrow(gene_sequence)==1) {
  cat("Only one sequence found for that gene \n")
  Isoform_of_interest = 1
}

if (nrow(gene_sequence)>1) {
  cat("More than one sequence found for that gene. Please select the transcript to use ! \n ")
  print(gene_sequence[,-1])
  cat("Select the transcript of interest (number of the row) \n")
  Isoform_of_interest = readLines('stdin',n=1)
  Isoform_of_interest = as.numeric(Isoform_of_interest)
  while (is.na(Isoform_of_interest) | Isoform_of_interest <= 0 ) {
    cat("This is not a number ! Please provide a positive integer !  \n")
    Isoform_of_interest = readLines('stdin',n=1)
    Isoform_of_interest = as.numeric(Isoform_of_interest)
  }
}

#Creating the Biostring object

gene_sequence_target =gene_sequence$coding[Isoform_of_interest]
gene_sequence_target = DNAStringSet(gene_sequence_target)
names(gene_sequence_target) = paste(gene_symbol,gene_sequence$ensembl_transcript_id[Isoform_of_interest])

#Creating the Fasta directory if it does not exist already and exporting the sequences used as fasta files

Fasta_directory = paste(Output_directory,"/Fasta/",sep = "")
if (!dir.exists(Fasta_directory)) {
  dir.create(Fasta_directory)
}

writeXStringSet(x =(gene_sequence_target), filepath = paste(Fasta_directory,gene_symbol,".fasta",sep = "") , append=FALSE,
                compress=FALSE, format="fasta",)

###Settings for probe design 

gene_sequence_target = gene_sequence_target[[1]]
l = length(gene_sequence_target)
Length_probe = as.numeric(Parameters["Length_probe",])
Probe_space = as.numeric(Parameters["Probe_space",])
Min_GC_percent = as.numeric(Parameters["Min_GC_percent",])
Max_GC_percent = as.numeric(Parameters["Max_GC_percent",])
Forbiden_sequences = DNAStringSet(x = c("AAAAA","TTTTT","GGGGG","CCCCC"))

##Designing the probes themselves

k = 1 
end_reached = F
List_probes = DNAStringSet()
List_probes_names = c()
n = 1

while (k < (l-Length_probe)  & !end_reached) {
  
  Temp_sequence = gene_sequence_target[k:(k+Length_probe-1)]
  
  ##Cheching GC percent is in the range...
  GC_percent = letterFrequency(Temp_sequence,letters = c("A","T","C","G"))
  GC_percent = (GC_percent["G"] + GC_percent["C"])/Length_probe * 100
  GC_percent_in_range = GC_percent > Min_GC_percent &  GC_percent < Max_GC_percent
  
  ##Checking for the absence of repetitive sequences 
  
  Forbiden_sequences_presence = 0
  for (i in 1:length(Forbiden_sequences)) {
    x=countPattern(Forbiden_sequences[[i]],Temp_sequence)
    Forbiden_sequences_presence = Forbiden_sequences_presence + x
  }
  Forbiden_sequences_presence = as.logical(Forbiden_sequences_presence)
  
  if (GC_percent_in_range & !Forbiden_sequences_presence) {
    #Saving the probe sequence
    List_probes[[n]] = Temp_sequence
    List_probes_names = c(List_probes_names,paste(gene_symbol,as.character(k),"-",as.character(k + Length_probe -1)))
    
    #Updating the variable
    n = n + 1 
    k = k + Length_probe + Probe_space
  } else {
    k = k + 1
  }
  
}
names(List_probes) = List_probes_names
cat("Creation of probes done ! \n")
cat(paste(as.character(length(List_probes)),"probes have been designed and passed first QC \n"))

##Exporting the sequences of the probes 
cat("Exporting the probes ...")
Subsequence_directory = paste(Output_directory,"/Subsequences/",sep = "")
if (!dir.exists(Subsequence_directory)) {
  dir.create(Subsequence_directory)
}

writeXStringSet(x =List_probes, filepath = paste(Subsequence_directory,gene_symbol,".fasta",sep = "") , append=FALSE,
                compress=FALSE, format="fasta")
cat(" done ! \n")


#Creating the alignment directory if it does not already exist

Alignment_directory = paste(Output_directory,"/Alignment/",sep = "")
if (!dir.exists(Alignment_directory)) {
  dir.create(Alignment_directory)
}

#Creating the alignment command

cat("Mapping the probes to reference transcriptome")
Bowtie2_path = Parameters["Bowtie2_path",]
Mapping_command = paste(Bowtie2_path," -x ", Transcriptome_index  ,"-f " ,paste(Subsequence_directory,gene_symbol,".fasta",sep = ""), 
                        "--no-hd -t -k 100 --very-sensitive -S", paste(Alignment_directory,gene_symbol,".sam",sep = ""),"--norc")
#Launching the command 
system(Mapping_command)


##Reading the results of the mapping 
SAM_file = read.delim(paste(Alignment_directory,gene_symbol,".sam",sep = ""),header = F,sep = "\t")
SAM_file = SAM_file[complete.cases(SAM_file),]
colnames(SAM_file) = c("Probe","Meta_score","Target","Position","Score_quality","Alignment","Ceiling","Offset","Fragment","Segment","Quality")
SAM_file$Segment = as.character(SAM_file$Segment)
n_probes = length(unique(SAM_file$Segment))


n_mapping = table(SAM_file$Segment)
order_probes = order(n_mapping,decreasing = F)

SAM_file_ordered = c()

for (j in order_probes) {
  p = names(n_mapping[j])
  u = SAM_file[SAM_file$Segment==p,]
  v = c(unique(u$Segment),n_mapping[j],paste(as.character(u$Target),collapse = " "),u$Position)
  SAM_file_ordered = rbind(SAM_file_ordered,v)
}
rownames(SAM_file_ordered) = 1:nrow(SAM_file_ordered)
colnames(SAM_file_ordered) = c("Segment","N_mapping","Mapped_genes","Position")
SAM_file_ordered = as.data.frame(SAM_file_ordered)

#If too many high quality probes : select the probes that maximise the spread over the transcript

n_unique_mapping = table(n_mapping)["1"]
if (is.na(n_unique_mapping)) {
  n_unique_mapping =1
}

Max_probes = as.numeric(Parameters["Max_probes",])

if (n_unique_mapping>=Max_probes) {
  
  cat(paste("More than",as.character(Max_probes),"higly specific probes found for gene", gene_symbol,": selecting the most spread probes on the transcript! \n"))
  SAM_file_filtered =  SAM_file[as.character(SAM_file$Probe) == gene_symbol,]
  SAM_file_filtered = SAM_file_filtered[SAM_file_filtered$Segment%in%(SAM_file_ordered$Segment[SAM_file_ordered$N_mapping==1]),]
  
  used_probes = c(which.min(SAM_file_filtered$Position),which.max(SAM_file_filtered$Position))
  d = as.matrix(dist(SAM_file_filtered$Position))^2
  
  for (i in 1:(Max_probes-2)) {
    d_prime = matrix(d[used_probes,],nrow = length(used_probes),byrow = F)
    d_prime = apply(d_prime,MARGIN = 2,FUN = min)
    used_probes = c(used_probes,which.max(d_prime))
  }
  SAM_file_final = data.frame(Segment = SAM_file_filtered$Segment[used_probes],N_mapping=1,Mapped_genes=SAM_file_filtered$Target[used_probes])
}
if (n_unique_mapping<Max_probes)  {
  SAM_file_final= SAM_file_ordered
}



#Now that were are happy we can take the reverse complement of the probes and export them :

cat("Taking the reverse complement of the probes. \n")
Sequences = SAM_file_final$Segment
Sequences = DNAStringSet(Sequences)
Sequences = reverseComplement(Sequences)
Sequences = as.character(Sequences)
SAM_file_final$Segment = Sequences


Final_directory = paste(Output_directory,"/Final_probes/",sep = "")
if (!dir.exists(Final_directory)) {
  dir.create(Final_directory)
}
write.table(SAM_file_final,file = paste(Final_directory,gene_symbol,".txt",sep = ""),sep = "\t",quote=F,row.names = F)


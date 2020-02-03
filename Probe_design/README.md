Probe design for smFISH and seqFISH
=============================

This guide has been written to efficiently design probe sets for single molecular (sm) and sequential (seq) FISH. It is inspired by probe design realised by Long Cai lab. So far it can only be used to design probes for human/mouse gene exons but can easily adapted to whole transcripts and to other species.

Installation of the required softwares, packages and database
-----------------------------------------------

This code has been tested on MacOS High Sierra and Ubuntu 18.0.4 but not on Windows. I would therefore recommend running this script on those two platforms  !
 
The first step is to install the latest version of  [R] (https://www.r-project.org/). Once this is done, two R packages have to be installed : **biomaRt** and **Biostrings**. How to install these two packages is described [there] (https://bioconductor.org/packages/release/bioc/html/biomaRt.html) and [there] (https://www.bioconductor.org/packages/release/bioc/html/Biostrings.html).

Probe design requires to map the designed probes against the whole specie transcriptome, therefore we will need a short read mapper called [Bowtie2] (http://bowtie-bio.sourceforge.net/bowtie2/). Please follow the instructions on the website for installation.

Lastly you have to download the fasta files corresponding to the whole CDS sequences of the targeted species on the [ensembl download page] (https://www.ensembl.org/info/data/ftp/index.html). The corresponding fasta file need to be transformed into a **bowtie2 index** before use. This can be done using the following command in Terminal :

```bash
bowtie2-build  /Path/to/CDS_file.fasta Index_directory/Index_CDS
```
This will generate various index files (.bt2 files) starting with the Index_CDS. 

Setting parameters values 
-----------------------------------------------

Several parameters need to be specified for probe design : the species, the targeted gene, the length of the probes... 
All those parameters can be tuned through a tab delimited text file that will be read by the R script. This text file corresponds to a two-column table with the first column containning the name of the parameter and the second the values.  A typical example of parameter file is provided with the script on the associated Github page. 

The available parameters are :

1. **Output_directory** : Path to the directory where the results of the probe design will be saved. If this directory does not already exist it will be created.
2. **Transcriptome_index** : Path to the Bowtie2 index with the suffix of the index files.
3. **Species** : Which species should be considered ? Can only be "Mouse" or "Human".
4. **Gene** : Name/symbol of the gene.
5. **Length_probe** : Length (in nucleotide) of each individual probe.
6. **Probe_space** :  Minimal distance (in nucleotide) between two neighbouring probes.
7. **Min_GC_percent** : Minimal GC content in percent (between 0 and 100).
8. **Max_GC_percent** : Maximal GC content in percent (between 0 and 100).
9. **Max_probes** : Maximal number of probes in a probeset.
10. **Bowtie2_path** : Path to Bowtie2 executable.


Running the script 
-----------------------------------------------

Once the parameters have been set the script can easily be launched through a bash command : 

```bash
Rscript Probe_design_script.R Path/to/Probe_parameter.txt
```
While in some cases you will not have to provide any additional information for the script to run, if several transcripts are found for the gene of interest you will have to specify which transcript to use.


Results of the script  
-----------------------------------------------

Four different sub-directories are created during probe design :
1. **Fasta** : contains the nucleotide sequences of the genes studied.
2. **Subsequences** : contains the sequences of the probes before any filtering.
3. **Alignment** : contains the SAM file resulting from the alignment of the probes to a reference transcriptome using bowtie2.
4. **Final_probes** : contains the final results of the analysis. Each analysed gene will produce a .txt file with one row per probe. The sequence of the probe as well as the gene to which it maps are provided.





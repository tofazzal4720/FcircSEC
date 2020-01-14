# FcircSEC is an R package for full length circRNA sequence extraction and classification

## Requirements
* R (>= 3.6.0)
* Biostrings
* seqRFLP
* stringi

## Installation
### From cran
To install the package from cran, run the command:

    install.packages("FcircSEC", dep=T)
	
### From github
To install the package from github first you need to install the package “devtools” using the following command:

    install.packages("devtools", dep=T)

The package "FcircSEC" depends on a bioconductor package "Biostrings" which cannot be installed automatically while installing "FicrcSEC" using "devtools". So, you need to install "Biostrings" manually using the following way:

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("Biostrings")

Finally, install “FcircSEC” by the following command:

    devtools::install_github("tofazzal4720/FcircSEC", dep = T)

Start analysis by typing the following command:

    library("FcircSEC")

## Extracting transcript information from the annotation file

Transcript data can be obtained using the following function:

    transcriptExtract(annotationFile, databaseName, outputfile)
Here,

`annotationFile` is the annotation file (in *gtf*, *gff* or *gff3* fromat) corresponding to the reference genome. Please use *gff* or *gff3* format for "ncbi" and *gtf* format for "ucsc" and "other". 

`databaseName` is the database name from where the annotation file was downloaded (the possible options are "ncbi", "ucsc" and "other").

`outputfile` is the name of the output file.

#### Examples
    #Loading an example annotation file and write to a file
    #Here temporary directory is created as input-output
    #directory. Please provide your own directory instead.
    out_dir<-tempdir()
    annotation_file<-data(refGenchr1)  
    annotation_file<-refGenchr1
    write.table(annotation_file, file.path(out_dir,"annotation_file.gtf"), 
         row.names=FALSE, sep="\t",quote=FALSE, col.names=FALSE)

    #Extraction of transcript information. Here, the output will be generated in file 
    #transcriptdata.txt in out_dir directory
    transcriptExtract(file.path(out_dir,"annotation_file.gtf"), "ucsc", 
    file.path(out_dir, "transcriptdata.txt"))

## Classifying circRNAs

circular RNAs can be classified using the following function:

    circClassification(transcriptdata, bedfile, outfiletxt, outfilebed)

Here,

`transcriptdata` is the transcript data extracted from the annotation file (obtained from function `transcriptExtract`).

`bedfile` is the bed file (obtained from the circRNA prediction tools) having four columns chromosome, start position, end position and strand of circRNAs.

`outfiletxt` is the output file with the detailed information of circRNA classification.

`outfilebed` is the output file with chromosome, start and end position of each circRNAs.

#### Examples
    #Loading and example transcript data and write to a file
    #Here temporary directory is created as input-output
    #directory. Please provide you own directory instead.
    out_dir<-tempdir()
    t_data<-data("transcript_data") 
    t_data<-transcript_data
    write.table(t_data, file.path(out_dir,"transcript_data.txt"), row.names=FALSE)

    #Loading an example bedfile obtained form the circRNA prediction tool and write to a file
    b_file<-data("output_CIRI")
    b_file<-output_CIRI
    write.table(b_file, file.path(out_dir,"output_CIRI.bed"), col.names=FALSE, row.names=FALSE)

    #Classification of circRNAs. Here, the output will be written in two files 
    #circRNA_class.txt and circRNA_class.bed in out_dir directory
    circClassification (file.path(out_dir,"transcript_data.txt"), 
            file.path(out_dir,"output_CIRI.bed"), file.path(out_dir, "circRNA_class.txt"), 
              file.path(out_dir, "circRNA_class.bed"))

## Generating sequences from the reference genome with specific intervals

Genomic sequences of the circRNAs is ontained from the reference genome for given circRNA boundary(start and end) using the following function:

    get.fasta(ref_genome, circ_class_bed, out_filename)

Here,

`ref_genome` is the reference genome.

`circ_class_bed` is the bed file having chromosome, start and end position of each circRNAs (obtained from function `circClassification`)

`out_filename` is the name of the output file.

#### Examples

    #Loading an example reference genome and write to a file
    #Here temporary directory is created as input-output
    #directory. Please provide you own directory instead.
    out_dir<-tempdir()
    ref_genom<-data("chr1")
    ref_genom<-chr1
    df.fasta=dataframe2fas(ref_genom, file.path(out_dir, "ref_genome.fasta"))

    #Loading an example circRNA classification bed file and write to a file
    circ_class_bed<-data("circRNA_classb")
    circ_class_bed<-circRNA_classb
    write.table(circ_class_bed, file.path(out_dir, "circ_class.bed"), 
          col.names=FALSE, row.names=FALSE)

    #Getting genomic sequences of circRNAs. The output will be 
    #generated in file circRNA_genomic_seq.fasta in out_dir directory
    get.fasta(file.path(out_dir, "ref_genome.fasta"), 
          file.path(out_dir, "circ_class.bed"), 
            file.path(out_dir, "circRNA_genomic_seq.fasta"))

## Generating full length circRNA sequences

The full length circRNA sequences are obtained using the following function:

    circSeqExt(genomic_seq, circ_class_txt, out_filename)

Here,

`genomic_seq` is the *fasta* file (obtained using function `get.fasta`) having the genomic sequences for circRNAs.

`circ_class_txt` is the circRNA classification file (obtained from function `circClassification`).

`out_filename` is the name of the output file.

#### Examples
    #Loading an example circRNA genomic sequence and write to a file
    #Here temporary directory is created as input-output
    #directory. Please provide you own directory instead.
    out_dir<-tempdir()
    circ_genomic_seq<-data("circRNA_genomic_sequence")
    circ_genomic_seq<-circRNA_genomic_sequence
    df.fasta=dataframe2fas(circ_genomic_seq, file.path(out_dir, "circ_genomic_seq.fasta"))

    #Loading an example circ_class_txt data and write to a file
    circ_class_txt<-data("circRNA_classt")
    circ_class_txt<-circRNA_classt
    write.table(circ_class_txt, file.path(out_dir, "circ_class.txt"), 
          row.names=FALSE)

    #Extracting full length circRNA sequences. Here, the output will be 
    #written in file circRNA_sequence.fasta in out_dir directory
    circSeqExt(file.path(out_dir, "circ_genomic_seq.fasta"),
         file.path(out_dir, "circ_class.txt"), file.path(out_dir, "circRNA_sequence.fasta"))


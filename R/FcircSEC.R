#'Extracting transcript information from the annotation file 
#'
#' This function extracts transcript information from the annotation file corresponding to reference genome
#' @param annotationFile The annotation file (in gtf, gff or gff3 fromat) corresponding to the reference genome
#' @param databaseName The database name from where the annotation file was downloaded (the possible options are "ncbi", "ucsc" and "other")
#' @param outputfile The name of the output file
#' @return The transcript information from the annotation file will be written in the output file 'outputfile'
#' @importFrom utils write.table
#' @importFrom utils read.table
#' @importFrom utils read.delim
#' @examples 
#' 
#' #Loading an example annotation file and write to a file
#' #Here temporary directory is created as input-output
#' #directory. Please provide your own directory instead.
#' out_dir<-tempdir()
#' annotation_file<-data(refGenchr1)  
#' annotation_file<-refGenchr1
#' write.table(annotation_file, file.path(out_dir,"annotation_file.gtf"), 
#'          row.names=FALSE, sep="\t",quote=FALSE, col.names=FALSE)
#' 
#' #Extraction of transcript information. Here, the output will be generated in file 
#' #transcriptdata.txt in out_dir directory
#' transcriptExtract(file.path(out_dir,"annotation_file.gtf"), "ucsc", 
#'     file.path(out_dir, "transcriptdata.txt"))
#' 
#' @export
#' 


transcriptExtract<-function(annotationFile, databaseName, outputfile){
  
  getFileNameExtension <- function (fn) {
    # remove a path
    splitted    <- strsplit(x=fn, split='/')[[1]]   
    # or use .Platform$file.sep in stead of '/'
    fn          <- splitted [length(splitted)]
    ext         <- ''
    splitted    <- strsplit(x=fn, split='\\.')[[1]]
    l           <-length (splitted)
    if (l > 1 && sum(splitted[1:(l-1)] != ''))  ext <-splitted [l] 
    # the extention must be the suffix of a non-empty name    
    ext
  }
  ext<-getFileNameExtension(annotationFile)
  
  if(ext=="gff" | ext=="gff3"){
    annot_file<- read.delim(annotationFile, header = FALSE, sep = '\t', skip = 8)
  } else if(ext=="gtf") {
    annot_file<-read.table(annotationFile, header = FALSE, sep = '\t')
  }else {
    print("Error: please input the annotation file in gff or gtf format as the first argument")
  }
  
  
  
  if(databaseName=="ncbi"){
    gff_exon<-annot_file[which(annot_file$V3=="exon"),]
    
    trans_id<-sub(".*transcript_id=", "", gff_exon$V9)
    if (length(which(grepl("ID=",trans_id)==T))==length(trans_id) | length(which(grepl("gene_id ",trans_id)==T))==length(trans_id)){
      print("Error: your provided annoation file is not an ncbi annotation file")
    }else{
      gene_id<-sub('.*gene=', '', gff_exon$V9)
      gene_id_final<-sub('\\;.*', '', gene_id)
      gff_exon_trans_id<-data.frame(gff_exon,trans_id, gene_id_final)
      
      wh<-which(grepl("ID=",trans_id)==T)
      
      gff_exon_trans_id_final<-gff_exon_trans_id[-wh,]
      
      
      trans_chr_id<-paste(gff_exon_trans_id_final$trans_id, gff_exon_trans_id_final$V1, sep=':')
      
      gff_exon_trans_chr_id_final<-data.frame(gff_exon_trans_id_final,trans_chr_id)
      
      gff_exon_final<-gff_exon_trans_chr_id_final
      
      gff_exon_final<-gff_exon_final[-9]
      
      gff_exon_final$V4<-gff_exon_final$V4-1
      
      colnames(gff_exon_final)<-c("V1","V2", "V3","V4","V5", "V6", "V7", "V8", "V9", "V10", "V11")
      
      trans_chr_id_1<-gff_exon_final$V11
      dup<-duplicated(trans_chr_id_1)
      trans_chr_id_uni<-trans_chr_id_1[-which(dup==T)]
      
      transcript_id=NULL
      chr=NULL
      strand=NULL
      trans_start=NULL
      trans_end=NULL
      exon_count=NULL
      exon_starts=NULL
      exon_ends=NULL
      gene=NULL
      
      for(i in 1:length(trans_chr_id_uni)){
        total_id<-gff_exon_final[which(gff_exon_final$V11==trans_chr_id_uni[i]),]
        exon_id<-data.frame(total_id$V4,total_id$V5)
        exon_id_sort<-exon_id[order(exon_id[,1]),]
        transcript_id[i]=as.character(total_id$V9[1])
        chr[i]=as.character(total_id$V1[1])
        strand[i]=as.character(total_id$V7[1])
        trans_start[i]=exon_id_sort[1,1]
        trans_end[i]=exon_id_sort[dim(exon_id_sort)[1],2]
        exon_count[i]=dim(exon_id_sort)[1]
        exon_starts[i]=paste(exon_id_sort[,1],collapse=",")
        exon_ends[i]=paste(exon_id_sort[,2],collapse=",")
        gene[i]=as.character(total_id$V10[1])
      }
      trans_data<-data.frame(transcript_id, chr, strand, trans_start, trans_end, exon_count, exon_starts, exon_ends,gene)
      
      write.table(trans_data,outputfile, sep="\t",quote=F,row.names=F)
    }
  }else if (databaseName =="ucsc"){
    gtf_exon<-annot_file[which(annot_file$V3=="exon"),]
    
    trans_id<-sub(".*transcript_id ", "", gtf_exon$V9)
    if(length(which(grepl("ID=",trans_id)==T))==length(trans_id) | length(which(grepl("gene_id ",trans_id)==T))==length(trans_id)){
      print("Error: your provided annoation file is not an ucsc annotation file")
    }else{
      trans_id_final<-sub('\\;.*', '', trans_id)
      
      gene_id<-sub('.*gene_name ', '', gtf_exon$V9)
      gene_id_final<-sub('\\;.*', '', gene_id)
      
      gtf_exon_trans_id<-data.frame(gtf_exon,trans_id_final, gene_id_final)
      
      wh<-which(grepl("gene_id",trans_id_final)==T)
      
      if (length(wh)>0)gtf_exon_trans_id_final<-gtf_exon_trans_id[-wh,] else gtf_exon_trans_id_final<-gtf_exon_trans_id
      
      trans_chr_id<-paste(gtf_exon_trans_id_final$trans_id_final, gtf_exon_trans_id_final$V1, sep=':')
      
      gtf_exon_trans_chr_id_final<-data.frame(gtf_exon_trans_id_final,trans_chr_id)
      
      gtf_exon_final<-gtf_exon_trans_chr_id_final
      
      gtf_exon_final<-gtf_exon_final[-9]
      
      gtf_exon_final$V4<-gtf_exon_final$V4-1
      
      colnames(gtf_exon_final)<-c("V1","V2", "V3","V4","V5", "V6", "V7", "V8", "V9", "V10", "V11")
      
      trans_chr_id_1<-gtf_exon_final$V11
      length(trans_chr_id_1)
      dup<-duplicated(trans_chr_id_1)
      trans_chr_id_uni<-trans_chr_id_1[-which(dup==T)]
      
      
      transcript_id=NULL
      chr=NULL
      strand=NULL
      trans_start=NULL
      trans_end=NULL
      exon_count=NULL
      exon_starts=NULL
      exon_ends=NULL
      gene=NULL
      
      for(i in 1:length(trans_chr_id_uni)){
        total_id<-gtf_exon_final[which(gtf_exon_final$V11==trans_chr_id_uni[i]),]
        exon_id<-data.frame(total_id$V4,total_id$V5)
        exon_id_sort<-exon_id[order(exon_id[,1]),]
        transcript_id[i]=as.character(total_id$V9[1])
        chr[i]=as.character(total_id$V1[1])
        strand[i]=as.character(total_id$V7[1])
        trans_start[i]=exon_id_sort[1,1]
        trans_end[i]=exon_id_sort[dim(exon_id_sort)[1],2]
        exon_count[i]=dim(exon_id_sort)[1]
        exon_starts[i]=paste(exon_id_sort[,1],collapse=",")
        exon_ends[i]=paste(exon_id_sort[,2],collapse=",")
        gene[i]=as.character(total_id$V10[1])
      }
      trans_data<-data.frame(transcript_id, chr, strand, trans_start, trans_end, exon_count, exon_starts, exon_ends,gene)
      
      write.table(trans_data,outputfile,sep="\t",quote=F,row.names=F)
    }
  }else if (databaseName =="other"){
    gtf_exon<-annot_file[which(annot_file$V3=="exon"),]
    
    trans_id<-sub(".*transcript_id ", "", gtf_exon$V9)
    if(length(which(grepl("ID=",trans_id)==T))==length(trans_id) | length(which(grepl("gene_id ",trans_id)==T))==length(trans_id)){
      print("Error: your provided annoation file is not supported by our method, please use another annotation file")
    }else{
      trans_id_final<-sub('\\;.*', '', trans_id)
      
      gene_id<-sub('.*gene_name ', '', gtf_exon$V9)
      gene_id_final<-sub('\\;.*', '', gene_id)
      
      gtf_exon_trans_id<-data.frame(gtf_exon,trans_id_final, gene_id_final)
      
      wh<-which(grepl("gene_id",trans_id_final)==T)
      
      if (length(wh)>0)gtf_exon_trans_id_final<-gtf_exon_trans_id[-wh,] else gtf_exon_trans_id_final<-gtf_exon_trans_id
      
      trans_chr_id<-paste(gtf_exon_trans_id_final$trans_id_final, gtf_exon_trans_id_final$V1, sep=':')
      
      gtf_exon_trans_chr_id_final<-data.frame(gtf_exon_trans_id_final,trans_chr_id)
      
      gtf_exon_final<-gtf_exon_trans_chr_id_final
      
      gtf_exon_final<-gtf_exon_final[-9]
      
      gtf_exon_final$V4<-gtf_exon_final$V4-1
      
      colnames(gtf_exon_final)<-c("V1","V2", "V3","V4","V5", "V6", "V7", "V8", "V9", "V10", "V11")
      
      trans_chr_id_1<-gtf_exon_final$V11
      length(trans_chr_id_1)
      dup<-duplicated(trans_chr_id_1)
      trans_chr_id_uni<-trans_chr_id_1[-which(dup==T)]
      
      
      transcript_id=NULL
      chr=NULL
      strand=NULL
      trans_start=NULL
      trans_end=NULL
      exon_count=NULL
      exon_starts=NULL
      exon_ends=NULL
      gene=NULL
      
      for(i in 1:length(trans_chr_id_uni)){
        total_id<-gtf_exon_final[which(gtf_exon_final$V11==trans_chr_id_uni[i]),]
        exon_id<-data.frame(total_id$V4,total_id$V5)
        exon_id_sort<-exon_id[order(exon_id[,1]),]
        transcript_id[i]=as.character(total_id$V9[1])
        chr[i]=as.character(total_id$V1[1])
        strand[i]=as.character(total_id$V7[1])
        trans_start[i]=exon_id_sort[1,1]
        trans_end[i]=exon_id_sort[dim(exon_id_sort)[1],2]
        exon_count[i]=dim(exon_id_sort)[1]
        exon_starts[i]=paste(exon_id_sort[,1],collapse=",")
        exon_ends[i]=paste(exon_id_sort[,2],collapse=",")
        gene[i]=as.character(total_id$V10[1])
      }
      trans_data<-data.frame(transcript_id, chr, strand, trans_start, trans_end, exon_count, exon_starts, exon_ends,gene)
      
      write.table(trans_data,outputfile,sep="\t",quote=F,row.names=F)
    }
  }else{
    print("Error: please input 'ncbi' or 'ucsc' or 'other' as the second argument")
  } 
}



#' circRNA classification using trancript information and the bed file from the circRNA prediction tools
#' 
#' This function classifies circRNAs using the transcript information obtained from annotation file and the bedfile obtained from the circRNA prediction tools 
#' @param transcriptdata The transcript data (obtained from function \code{\link[FcircSEC]{transcriptExtract}})
#' @param bedfile The bed file (obtained from the circRNA prediction tools) having four columns chromosome, circRNA start, circRNA end position and circRNA strand
#' @param outfiletxt The output file with the detailed information of circRNA classification
#' @param outfilebed The output file with chromosome, start and end position of each circRNAs
#' @return The detailed information of circRNA classification will be written in outfiletxt and only chromosome, start and end position of each circRNAs will be written in outfilebed
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @examples 
#' 
#' #Loading and example transcript data and write to a file
#' #Here temporary directory is created as input-output
#' #directory. Please provide you own directory instead.
#' out_dir<-tempdir()
#' t_data<-data("transcript_data") 
#' t_data<-transcript_data
#' write.table(t_data, file.path(out_dir,"transcript_data.txt"), row.names=FALSE)
#' 
#' #Loading an example bedfile obtained form the circRNA prediction tool and write to a file
#' b_file<-data("output_CIRI")
#' b_file<-output_CIRI
#' write.table(b_file, file.path(out_dir,"output_CIRI.bed"), col.names=FALSE, row.names=FALSE)
#' 
#' #Classification of circRNAs. Here, the output will be written in two files 
#' #circRNA_class.txt and circRNA_class.bed in out_dir directory
#' circClassification (file.path(out_dir,"transcript_data.txt"), 
#'    file.path(out_dir,"output_CIRI.bed"), file.path(out_dir, "circRNA_class.txt"), 
#'    file.path(out_dir, "circRNA_class.bed"))
#' 
#' @export
#' 


circClassification<-function(transcriptdata, bedfile, outfiletxt, outfilebed){
  transcript<-read.table(transcriptdata,header=T) #  reading transcript data file
  output_circ_tool<-read.table(bedfile) # reading bed file from the output of circRNA prediction tool#
  Splice_length=NULL
  CircRNA_type=NULL
  b_transcript=NULL
  b_strand=NULL
  b_trans_start=NULL
  b_trans_end=NULL
  b_gene=NULL
  e_count=NULL
  e_sizes=NULL
  e_offsets=NULL
  for(i in 1:dim(output_circ_tool)[1]){
    s=as.numeric(as.character(output_circ_tool$V2[i]))
    e=as.numeric(as.character(output_circ_tool$V3[i]))
    wh<-which(as.character(transcript$chr)==as.character(output_circ_tool$V1[i]) & transcript$trans_start<e & transcript$trans_end>s)
    if (length(wh)==0){
      circRNA_type="other"
      splice_length=e-s
      B_transcript="NA"
      B_strand="NA"
      B_trans_start="NA"
      B_trans_end="NA"
      B_gene="NA"
      e_count_1=1
      e_sizes_1=e-s
      e_offsets_1=0
    }else{
      trans<-transcript[wh,]
      splice_L=NULL
      dm=NULL
      convse_1=NULL
      t_length=NULL
      for(j in 1: dim(trans)[1]){
        a<-do.call("rbind", strsplit(as.character(trans$exon_starts[j]),","))
        a1<-data.frame(apply(a, 2, as.numeric))
        b<-do.call("rbind", strsplit(as.character(trans$exon_ends[j]),","))
        b1<-data.frame(apply(b, 2, as.numeric))
        colnames(a1)<-"start"
        colnames(b1)<-"end"
        exon_start_end<-cbind(a1,b1)
        
        wh1<-exon_start_end[which(exon_start_end$start<e & exon_start_end$end>s),]
        if (dim(wh1)[1]==0){
          i1<-which(exon_start_end$end<=s);  i1<-i1[length(i1)]
          i2<-which(exon_start_end$start>=e); i2<-i2[1]
          wh2<-exon_start_end[i1:i2,]
          convse_1[j]<-(abs(wh2$end[1]-s)+abs(wh2$start[2]-e))/2
          splice_L[j]<-e-s
        }else if(wh1$start[1]==s & wh1$end[dim(wh1)[1]]==e){
          convse_1[j]<-0
          splice_L[j]<-sum(wh1$end-wh1$start)
        }else{
          convse_1[j]<-(abs(wh1$start[1]-s)+abs(wh1$end[dim(wh1)[1]]-e))/2
          splice_L[j]<-e-s
        }
        dm[j]=dim(wh1)[1]
        t_length[j]=trans$trans_end[j]-trans$trans_start[j]
      }
      
      d_l<-data.frame(dm,splice_L,convse_1,t_length)
      
      if (min(d_l$convse_1)==0){
        d_2<-d_l[which(d_l$convse_1==min(d_l$convse_1)),]
        d_3<-d_2[which(d_2$splice_L==max(d_2$splice_L)),]
        d_4<-d_3[which(d_3$t_length==max(d_3$t_length)),]
      }else if (max(d_l$dm)==0) {
        d_2<-d_l[which(d_l$convse_1==min(d_l$convse_1)),]
        d_3<-d_2[which(d_2$t_length==max(d_2$t_length)),]
        d_4=d_3
      }else {
        d_2<-d_l[which(d_l$dm==max(d_l$dm)),]
        d_3<-d_2[which(d_2$t_length==max(d_2$t_length)),]
        d_4=d_3
      }
      
      tr<-as.numeric(as.character(row.names(d_4)))
      
      trans_1<-trans[tr,]
      
      trans_3=trans_1[1,]     
      B_transcript=as.character(trans_3$transcript_id)
      B_strand=as.character(trans_3$strand)
      B_trans_start=as.character(trans_3$trans_start)
      B_trans_end=as.character(trans_3$trans_end)
      B_gene=as.character(trans_3$gene)
      pp<-do.call("rbind", strsplit(as.character(trans_3$exon_starts),","))
      ppp<-data.frame(apply(pp, 2, as.numeric))
      pp1<-do.call("rbind", strsplit(as.character(trans_3$exon_ends),","))
      ppp1<-data.frame(apply(pp1, 2, as.numeric))
      colnames(ppp)<-"start"
      colnames(ppp1)<-"end"
      exon_start_end<-cbind(ppp,ppp1)
      
      wh1<-exon_start_end[which(exon_start_end$start<e & exon_start_end$end>s),]
      
      if(dim(wh1)[1]>0){
        if(wh1$start[1]<=s & wh1$end[dim(wh1)[1]]>=e){
          wh1$start[1]=s
          wh1$end[dim(wh1)[1]]=e
          circRNA_type="exonic"
          splice_length<-sum(wh1$end-wh1$start)
          e_count_1=dim(wh1)[1]
          e_sizes_1=wh1$end-wh1$start
          e_offsets_1=wh1$start-s
        }else{
          circRNA_type="other"
          splice_length<-e-s
          e_count_1=1
          e_sizes_1=e-s
          e_offsets_1=0
        }
      }else{
        circRNA_type="intronic"
        splice_length=e-s
        e_count_1=1
        e_sizes_1=e-s
        e_offsets_1=0
      }
    }
    
    Splice_length[i]=splice_length
    CircRNA_type[i]=circRNA_type
    b_transcript[i]=B_transcript
    b_strand[i]=B_strand
    b_trans_start[i]=B_trans_start
    b_trans_end[i]=B_trans_end
    b_gene[i]=B_gene
    e_count[i]=e_count_1
    e_sizes[i]=paste(e_sizes_1, collapse=",")
    e_offsets[i]=paste(e_offsets_1, collapse=",")
  }
  cRNA_id<-paste(output_circ_tool$V1,output_circ_tool$V2, sep=':')
  cRNA_id_final<-paste(cRNA_id,output_circ_tool$V3, sep='-')
  
  circ_class<-cbind(cRNA_id_final, output_circ_tool[,1:4],Splice_length, CircRNA_type, e_count, e_sizes, e_offsets, b_transcript, b_strand, b_trans_start, b_trans_end, b_gene)
  colnames(circ_class)<-c("ID", "chr", "circ_start", "circ_end", "circ_strand", "splice_L", "circ_type", "e_count", "e_sizes", "e_offsets", "b_transcript", "b_strand", "b_trans_start", "b_trans_end", "b_gene")
  
  write.table(circ_class, outfiletxt, sep="\t",quote=F,row.names=F) # writing circ_class.txt file#
  
  circ_class_bed<-circ_class[,2:4]
  write.table(circ_class_bed, outfilebed, sep="\t",quote=F,row.names=F,col.names=F) # writing circ_class.bed file#
}



#' Generating sequences from the reference genome with specific intervals
#'
#' This function can extract the sequences from the reference genome for the given intervals (start, end) of chromosomes
#' @param ref_genome The reference genome
#' @param circ_class_bed  The bed file having chromosome, start and end position of each circRNAs (obtained from function \code{\link[FcircSEC]{circClassification}})
#' @param out_filename The name of the output file
#' @return The fasta file of the sequences extracted from the reference genome for the given intervals will be written in the output file 'out_filename'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom seqRFLP dataframe2fas
#' @importFrom utils read.table
#' @examples 
#'
#' #Loading an example reference genome and write to a file
#' #Here temporary directory is created as input-output
#' #directory. Please provide you own directory instead.
#' out_dir<-tempdir()
#' ref_genom<-data("chr1")
#' ref_genom<-chr1
#' df.fasta=dataframe2fas(ref_genom, file.path(out_dir, "ref_genome.fasta"))
#' 
#' #Loading an example circRNA classification bed file and write to a file
#' circ_class_bed<-data("circRNA_classb")
#' circ_class_bed<-circRNA_classb
#' write.table(circ_class_bed, file.path(out_dir, "circ_class.bed"), 
#'     col.names=FALSE, row.names=FALSE)
#' 
#' #Getting genomic sequences of circRNAs. The output will be 
#' #generated in file circRNA_genomic_seq.fasta in out_dir directory
#' get.fasta(file.path(out_dir, "ref_genome.fasta"), 
#'    file.path(out_dir, "circ_class.bed"), 
#'    file.path(out_dir, "circRNA_genomic_seq.fasta"))
#' 
#' @export
#'

get.fasta<-function(ref_genome, circ_class_bed, out_filename){
  fastaFile <- readDNAStringSet(ref_genome) # reading reference genome
  seq_name = sub('\\ .*', '', names(fastaFile))
  sequence = paste(fastaFile)
  df <- data.frame(seq_name, sequence)
  circ_class<-read.table(circ_class_bed) # reading circ_class_bed file
  a<-as.character(df[,2])
  B=NULL
  for (i in 1:dim(circ_class)[1]){
    v<-which(as.character(circ_class$V1[i])==as.character(df[,1]))
    B[i]=substr(a[v], circ_class$V2[i]+1, circ_class$V3[i])
  }
  id_1<-paste(circ_class$V1, circ_class$V2, sep=':')
  id_2<-paste(id_1, circ_class$V3, sep='-')
  csv<-data.frame(id_2, B)
  colnames(csv)<-c("name", "seq")
  df.fasta=dataframe2fas(csv, out_filename) #writing genomic sequence as fasta file
}


#' Generating full length circRNA sequences 
#'
#' This function can extract the full length circRNA sequences from the output of the circular RNA predictions tools
#' @param genomic_seq A fasta file (obtain using function \code{\link[FcircSEC]{get.fasta}}) with the genomic sequences for circRNAs
#' @param circ_class_txt The circRNA classification file (obtained from function \code{\link[FcircSEC]{circClassification}})
#' @param out_filename The name of the output file
#' @return The fasta file containing the full length circRNA sequences will be written in the output file 'out_filename'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stringi stri_reverse
#' @importFrom seqRFLP dataframe2fas
#' @importFrom utils read.table
#' @examples 
#' 
#' #Loading an example circRNA genomic sequence and write to a file
#' #Here temporary directory is created as input-output
#' #directory. Please provide you own directory instead.
#' out_dir<-tempdir()
#' circ_genomic_seq<-data("circRNA_genomic_sequence")
#' circ_genomic_seq<-circRNA_genomic_sequence
#' df.fasta=dataframe2fas(circ_genomic_seq, file.path(out_dir, "circ_genomic_seq.fasta"))
#' 
#' #Loading an example circ_class_txt data and write to a file
#' circ_class_txt<-data("circRNA_classt")
#' circ_class_txt<-circRNA_classt
#' write.table(circ_class_txt, file.path(out_dir, "circ_class.txt"), 
#'         row.names=FALSE)
#' 
#' #Extracting full length circRNA sequences. Here, the output will be 
#' #written in file circRNA_sequence.fasta in out_dir directory
#' circSeqExt(file.path(out_dir, "circ_genomic_seq.fasta"),
#'  file.path(out_dir, "circ_class.txt"), file.path(out_dir, "circRNA_sequence.fasta"))
#' 
#' @export
#'

circSeqExt<-function (genomic_seq, circ_class_txt, out_filename){
  fastaFile1 <- readDNAStringSet(genomic_seq) # reading genomic sequence for circRNAs
  seq_name1 = names(fastaFile1)
  sequence1 = paste(fastaFile1)
  fasta_file <- data.frame(seq_name1, sequence1)
  circ_classtxt<-read.table(circ_class_txt, header=T)
  C=NULL
  for (i in 1:nrow(fasta_file)){
    a<-fasta_file[i,2]
    a<-as.character(a)
    if(circ_classtxt$circ_strand[i]=="-" & is.na(circ_classtxt$b_transcript[i])==F ) a<-chartr("acgtACGT", "tgcaTGCA", a)
    B=NULL
    pp<-do.call("rbind", strsplit(as.character(circ_classtxt$e_sizes[i]),","))
    ppp<-data.frame(apply(pp, 2, as.numeric))
    pp1<-do.call("rbind", strsplit(as.character(circ_classtxt$e_offsets[i]),","))
    ppp1<-data.frame(apply(pp1, 2, as.numeric))
    for (j in 1:circ_classtxt$e_count[i]){
      B[j]<-substr(a, ppp1[j,1]+1, ppp1[j,1]+ppp[j,1])
    }
    C[i]<-paste(B,collapse="")
    if(circ_classtxt$circ_strand[i]=="-" & is.na(circ_classtxt$b_transcript[i])==F) C[i]<-stri_reverse(C[i])
  }
  csv<-data.frame(fasta_file[,1],C)
  colnames(csv)<-c("name", "seq")
  df.fasta = dataframe2fas(csv, out_filename)
}






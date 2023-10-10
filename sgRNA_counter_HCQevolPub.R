#create count table from known Brunello sgRNA targets and FASTQ file
library(Biostrings)
library(ShortRead)
library(data.table)
library(stringr)

#####user inputs#####
use_cluster <- TRUE
fastq_input_dir <- "CRISPR_screens/HCQ/"
analysis_dir <- "CRISPR_screens/HCQ/analysis_tables/"
filename_label_split_char <- "_R1_0" #to identify name within FASTQ file label

#map filtered files
fastq_filelist <- list.files(fastq_input_dir, pattern = "fastq.gz") #R1's only
input_fastq_R1_path <- paste0(fastq_input_dir,"fastq_filelist")

#CRISPR library information
proper_guide_sequence_5prime <- "GGAAAGGACGAAACACCG"
library_target_file_path <- "CRISPR_screens/broadgpp-brunello-library-contents.txt"
library_target_gene_column <- 2
library_target_guideseq_column <- 7


#####workflow#####
sample_names <- sapply(1:length(fastq_filelist), function(fastq_file){
  str_split(fastq_filelist[fastq_file], filename_label_split_char)[[1]][1]})

library_targets <- as.data.frame(fread(library_target_file_path, stringsAsFactors = FALSE, sep = "\t"))[
  ,c(library_target_gene_column, library_target_guideseq_column)]
colnames(library_targets) <- c("Gene","sgRNA_sequence")

sgRNA_count_table <- matrix(0, nrow = nrow(library_targets)
                            ,ncol = length(fastq_filelist))

row.names(sgRNA_count_table) <- library_targets$sgRNA_sequence
colnames(sgRNA_count_table) <- sample_names

#initialize variables
sample_index <- 1
sample_name <- sample_names[1]
seq <- NULL
fastq_piece <- NULL
seq_start_vector <- NULL
seq_end_vector <- NULL
sgRNA_seqs <- NULL
if(use_cluster == TRUE){
  library(snow)
  library(parallel)
}

for(sample_index in 1:length(fastq_filelist)){ #length(fastq_filelist)
  timer <- Sys.time()
  sample_name <- sample_names[sample_index]
  
#read FASTQ file and count sgRNAs
read_counter <- 0

 streamer_link_to_fastq_file <- FastqStreamer(paste0(fastq_input_dir,fastq_filelist[sample_index])
                       ,n = 1e6)
 
 repeat { ## process chunk loop start
   fastq_piece <- yield(streamer_link_to_fastq_file)@sread
   if (length(fastq_piece) == 0)
     break
   read_counter <- as.integer(read_counter + length(fastq_piece))
   print(paste0(sample_name," : read ", prettyNum(read_counter, big.mark = ",", scientific = FALSE) ))

sgRNA_start_indexes <- vmatchPattern(
  pattern = as.character(proper_guide_sequence_5prime)
  , subject = fastq_piece
  , max.mismatch = 1
  , min.mismatch = 0
  )@ends

R1_seq_to_analyze <- which(as.numeric(unlist(sapply(1:length(sgRNA_start_indexes), function(seq){if(is.null(sgRNA_start_indexes[seq][[1]])){return(0)} else {return(1)}} ))) != 0)

seq_start_vector <- rep(1, length(fastq_piece))
seq_end_vector   <- rep(50, length(fastq_piece))
seq_start_vector <- as.numeric(unlist(sapply(1:length(sgRNA_start_indexes), function(seq){
  if(is.null(sgRNA_start_indexes[[seq]][1])){return(1)} else {
    return(sgRNA_start_indexes[[seq]])}
  } )))
seq_end_vector   <- seq_start_vector + 20
too_short_reads <- which(seq_end_vector > 50)
seq_end_vector[too_short_reads] <- 50
R1_seq_to_analyze <- R1_seq_to_analyze[which(!(R1_seq_to_analyze %in% too_short_reads))]
  
if(use_cluster == TRUE){
  library(snow)
  library(parallel)
  cluster <- makeCluster(detectCores())
  clusterExport(cluster, c("fastq_piece","seq",
                           "seq_start_vector","seq_end_vector", "sgRNA_seqs"))
  sgRNA_seqs <- parSapply(cluster,R1_seq_to_analyze, function(seq){
    library(Biostrings)
    try({
    as.character(subseq(fastq_piece[[seq]]
                        , start = seq_start_vector[seq]+1
                        , end = seq_end_vector[seq]))
    })
  })
  stopCluster(cluster)
} else {
sgRNA_seqs <- sapply(R1_seq_to_analyze, function(seq){
  try({
  as.character(subseq(fastq_piece[[seq]]
         , start = seq_start_vector[seq]+1
         , end = seq_end_vector[seq]))
  })
})
}
sgRNA_seqs2 <- sgRNA_seqs[which(sgRNA_seqs %in% row.names(sgRNA_count_table))]

sgRNA_seqs_summary <- table(sgRNA_seqs2)
sgRNA_count_table[names(sgRNA_seqs_summary),sample_name] <- 
  sgRNA_count_table[names(sgRNA_seqs_summary),sample_name] +
  sgRNA_seqs_summary


 }## process chunk loop end
write.table(sgRNA_count_table, paste0(analysis_dir, "count_matrix.tsv")
             , sep = "\t", row.names = TRUE, col.names = NA)
timer_out <- Sys.time() - timer
print(timer_out)

}

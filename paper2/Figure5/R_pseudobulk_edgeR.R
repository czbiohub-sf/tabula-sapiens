library(edgeR)

output_dir = '/hpc/mydata/siyu.he/Siyu_projects/TS_project/pseudobulk_result/'
input_dir  = '/hpc/mydata/siyu.he/Siyu_projects/TS_project/pseudobulk_input/'
# List all files in the folder
file_list <- list.files(path = input_dir)


for (file_ in file_list[1:2500]){
  
  last_four_chars <- substr(file_, nchar(file_) - 8, nchar(file_))
  if (last_four_chars=='pbulk.csv'){
    print(file_)
    count_matrix <- as.matrix(read.csv(paste0(input_dir,file_)))
    gene_names <- count_matrix[,1]
    count_matrix <- count_matrix[,-1]
    
    #count_matrix <- matrix(count_matrix)
    
    if (!is.null(dim(count_matrix))) {    count_matrix <- apply(count_matrix, 2, as.numeric)
    
    
    modified_file_ <- substr(file_, 1, nchar(file_) - 4)
    modified_file_trimmed <- trimws(modified_file_, which = "right")
    
    if (!file.exists(paste0(output_dir,paste(modified_file_trimmed, "_edgeR_LRT_marker.csv", sep = "")))){
      print(file_)
      print(file_)
      
      result_string <- paste(modified_file_trimmed, '_obs.csv', sep = '')
      
      sample_info <- read.csv(paste0(input_dir,result_string))
      
      sample_info_sex <- sample_info$sex
      
      row.names(count_matrix) <-gene_names
      
      dge <- DGEList(counts = count_matrix, group = factor(sample_info_sex))
      
      dge <- calcNormFactors(object = dge)
      
      if ((length(unique(dge$samples$group))>1) & (!file.exists(paste0(output_dir,paste(modified_file_trimmed, "_edgeR_LRT_marker.csv", sep = "")))))
      {
        print(file_)
        design.mat <- model.matrix(~ 0 + dge$samples$group)
        colnames(design.mat) <- levels(dge$samples$group)
        dge2 <- estimateGLMCommonDisp(dge,design.mat)
        dge2 <- estimateGLMTrendedDisp(dge2,design.mat, method="power")
        dge2 <- estimateGLMTagwiseDisp(dge2,design.mat)
        fit <- glmFit(dge2, design.mat)
        lrt12 <- glmLRT(fit, contrast=c(1,-1))
        result_df <- topTags(lrt12, n=Inf)$table
        
        print(modified_file_trimmed)
        write.csv(result_df, file=paste0(output_dir,paste(modified_file_trimmed, "_edgeR_LRT_marker.csv", sep = "")), row.names=TRUE)
      }
    }}
  }
}


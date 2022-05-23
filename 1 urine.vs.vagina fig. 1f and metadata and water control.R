########### library ###########
library(stringr)
library(GUniFrac)
library(vegan)
library(ggplot2)
library(ALDEx2)
library(tidyr)
library(DiscriMiner)
library(compositions)


##### prepare samples ##################################################
# get reads of vaginal samples

setwd('/Users/binzhu/secure/momspi-data/home/metadata/urine_project')
metadata_v = read.table('allVaginal.16SDB.V1-V3.uclust99.summary.txt',header = T, sep = "\t")
metadata_u = read.table('allUrine.16SDB.V1-V3.uclust99.summary.txt',header = T, sep = "\t")
metadata_r = read.delim('allRectal.16SDB.V1-V3.uclust99.summary.txt',header = T, sep = "\t")
metadata_b = read.delim('allBuccal.16SDB.V1-V3.uclust99.summary.txt',header = T, sep = "\t")

metadata = read.table('metadata.txt',header = T)
metadata2 = read.table('age.txt',header = F)
match_cv= read.table('MCKD2MV1D.txt',header = F)
match_rv= read.table('MRCD2MV1D.txt',header = F)
match_uv= read.table('MU1W2MV1D_urineChris_022520.txt',header = F)

setwd('/Users/binzhu/Desktop/Urine')
# load metadata_all variable from '1 urine.vs.vagina.R'

##### data arrangement ###########################################
colnames(metadata_v) = c('V1','V2','V3','V4','V5','V6')
colnames(metadata_u) = c('V1','V2','V3','V4','V5','V6')
colnames(metadata_r) = c('V1','V2','V3','V4','V5','V6')
colnames(metadata_b) = c('V1','V2','V3','V4','V5','V6')

metadata_v$V3 = str_replace_all(metadata_v$V3,'AT','')
metadata_u$V3 = str_replace_all(metadata_u$V3,'AT','')
metadata_r$V3 = str_replace_all(metadata_r$V3,'AT','')
metadata_b$V3 = str_replace_all(metadata_b$V3,'AT','')

metadata_v_2 = as.data.frame(str_c(metadata_v$V2,metadata_v$V3,sep = " "))
metadata_v = cbind(metadata_v$V1,metadata_v_2,metadata_v$V4)
rm(metadata_v_2)
#colnames(metadata_v) = c('V1','V2','V3')

metadata_u_2 = as.data.frame(str_c(metadata_u$V2,metadata_u$V3,sep = " "))
metadata_u = cbind(metadata_u$V1,metadata_u_2,metadata_u$V4)
rm(metadata_u_2)
#colnames(metadata_u) = c('V1','V2','V3')

metadata_r_2 = as.data.frame(str_c(metadata_r$V2,metadata_r$V3,sep = " "))
metadata_r = cbind(metadata_r$V1,metadata_r_2,metadata_r$V4)
rm(metadata_r_2)
#colnames(metadata_r) = c('V1','V2','V3')

metadata_b_2 = as.data.frame(str_c(metadata_b$V2,metadata_b$V3,sep = " "))
metadata_b = cbind(metadata_b$V1,metadata_b_2,metadata_b$V4)
rm(metadata_b_2)
#colnames(metadata_b) = c('V1','V2','V3')

# get reads table
sample_name <- unique(metadata_v$`metadata_v$V1`)
species_name <- unique(metadata_v$`str_c(metadata_v$V2, metadata_v$V3, sep = " ")`)

reads_table_v <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table_v) = species_name 
colnames(reads_table_v) = sample_name

for (a in 1: dim(metadata_v)[1]) {
   
   column_num <- which(sample_name == metadata_v[a,1])
   row_num <- which(species_name == metadata_v[a,2])
   
   reads_table_v[row_num,column_num] = as.numeric(as.character(metadata_v[a,3]))
   
}
reads_table_v = as.data.frame(reads_table_v)


sample_name <- unique(metadata_u$`metadata_u$V1`)
species_name <- unique(metadata_u$`str_c(metadata_u$V2, metadata_u$V3, sep = " ")`)

reads_table_u <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table_u) = species_name 
colnames(reads_table_u) = sample_name

for (a in 1: dim(metadata_u)[1]) {
   
   column_num <- which(sample_name == metadata_u[a,1])
   row_num <- which(species_name == metadata_u[a,2])
   
   reads_table_u[row_num,column_num] = as.numeric(as.character(metadata_u[a,3]))
}

reads_table_u = as.data.frame(reads_table_u)


sample_name <- unique(metadata_r$`metadata_r$V1`)
species_name <- unique(metadata_r$`str_c(metadata_r$V2, metadata_r$V3, sep = " ")`)

reads_table_r <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table_r) = species_name 
colnames(reads_table_r) = sample_name

for (a in 1: dim(metadata_r)[1]) {
   
   column_num <- which(sample_name == metadata_r[a,1])
   row_num <- which(species_name == metadata_r[a,2])
   
   reads_table_r[row_num,column_num] = as.numeric(as.character(metadata_r[a,3]))
}

reads_table_r = as.data.frame(reads_table_r)


sample_name <- unique(metadata_b$`metadata_b$V1`)
species_name <- unique(metadata_b$`str_c(metadata_b$V2, metadata_b$V3, sep = " ")`)

reads_table_b <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table_b) = species_name 
colnames(reads_table_b) = sample_name

for (a in 1: dim(metadata_b)[1]) {
   
   column_num <- which(sample_name == metadata_b[a,1])
   row_num <- which(species_name == metadata_b[a,2])
   
   reads_table_b[row_num,column_num] = as.numeric(as.character(metadata_b[a,3]))
}

reads_table_b = as.data.frame(reads_table_b)

reads_table_v = reads_table_v[metadata_all$sampleID]
match_uv = match_uv[order(match(match_uv$V2,metadata_all$sampleID)),]
match_uv$V1 = str_replace_all(match_uv$V1, 'MU1W','MU1D')
reads_table_u = reads_table_u[match_uv$V1]

reads_table_v_2 = reads_table_v
colnames(reads_table_v_2) = paste0('Vagina_',c(1:84))
reads_table_u_2 = reads_table_u
colnames(reads_table_u_2) = paste0('Urine_',c(1:84))
keep = colSums(reads_table_v_2) >= 10000 & colSums(reads_table_u_2) >= 10000
sum(keep)
reads_table_v_2 = reads_table_v_2[,keep]
reads_table_u_2 = reads_table_u_2[,keep]

species_list = unique(c(row.names(reads_table_v_2),row.names(reads_table_u_2)))
species_list <- species_list[order(species_list)]
reads_table_u_2 <- reads_table_u_2[order(row.names(reads_table_u_2)), ]
reads_table_v_2 <- reads_table_v_2[order(row.names(reads_table_v_2)), ]

reads_table_u_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_u_2),nrow=length(species_list)))
reads_table_v_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_u_2),nrow=length(species_list)))
row.names(reads_table_u_new) = species_list
colnames(reads_table_u_new) = colnames(reads_table_u_2)
row.names(reads_table_v_new) = species_list
colnames(reads_table_v_new) = colnames(reads_table_v_2)

keep <- species_list %in% row.names(reads_table_v_2)
reads_table_v_new[keep,] = reads_table_v_2
keep <- species_list %in% row.names(reads_table_u_2)
reads_table_u_new[keep,] = reads_table_u_2
reads_table_all <- cbind(reads_table_v_new,reads_table_u_new)
write.csv(reads_table_all, 'urine.vs.vagina_no_th_16sDB.csv')

##### get sample pair #####
sample_list = colnames(reads_table_v)
keep = (sample_list %in% match_cv$V2) & (sample_list %in% match_rv$V2) & (sample_list %in% match_uv$V2)
sum(keep)
sample_list = sample_list[keep]
match_cv = match_cv[match_cv$V2 %in% sample_list,]
match_rv = match_rv[match_rv$V2 %in% sample_list,]
match_uv = match_uv[match_uv$V2 %in% sample_list,]

sample_list = sample_list[order(sample_list)]
match_cv = match_cv[order(match_cv$V2),]
match_rv = match_rv[order(match_rv$V2),]
match_uv = match_uv[order(match_uv$V2),]

reads_table_v = reads_table_v[, colnames(reads_table_v) %in% sample_list]
reads_table_v = reads_table_v[sample_list]
reads_table_b = reads_table_b[, colnames(reads_table_b) %in% match_cv$V1]
reads_table_b = reads_table_b[match_cv$V1]
reads_table_r = reads_table_r[, colnames(reads_table_r) %in% match_rv$V1]
reads_table_r = reads_table_r[match_rv$V1]
match_uv$V1 = str_replace_all(match_uv$V1, 'MU1W','MU1D')
reads_table_u = reads_table_u[, colnames(reads_table_u) %in% match_uv$V1]
reads_table_u = reads_table_u[match_uv$V1]

##### sample total reads threshold ####################
keep = colSums(reads_table_b) >= 10000 & colSums(reads_table_r) >= 10000 & 
   colSums(reads_table_u) >= 10000 &  colSums(reads_table_v) >= 10000
sum(keep)
reads_table_b = reads_table_b[,keep]
reads_table_v = reads_table_v[,keep]
reads_table_r = reads_table_r[,keep]
reads_table_u = reads_table_u[,keep]

metadata_all$sampleID_2 = c(1:84)
x = metadata_all$sampleID_2[match(colnames(reads_table_v), metadata_all$sampleID)]
colnames(reads_table_v) = paste0('Vagina',x)
reads_table_v = reads_table_v[order(colnames(reads_table_v))]

colnames(reads_table_u) = paste0('Urine',x)
reads_table_u = reads_table_u[order(colnames(reads_table_u))]

colnames(reads_table_r) = paste0('Rectum',x)
reads_table_r = reads_table_r[order(colnames(reads_table_r))]

colnames(reads_table_b) = paste0('Cheek',x)
reads_table_b = reads_table_b[order(colnames(reads_table_b))]

metadata_all = c(colnames(reads_table_v),colnames(reads_table_u),
                 colnames(reads_table_r),colnames(reads_table_b))
   

species_list = unique(c(row.names(reads_table_v),row.names(reads_table_u),
                        row.names(reads_table_r),row.names(reads_table_b)))
species_list <- species_list[order(species_list)]
reads_table_v <- reads_table_v[order(row.names(reads_table_v)), ]
reads_table_u <- reads_table_u[order(row.names(reads_table_u)), ]
reads_table_r <- reads_table_r[order(row.names(reads_table_r)), ]
reads_table_b <- reads_table_b[order(row.names(reads_table_b)), ]

reads_table_u_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_u),nrow=length(species_list)))
reads_table_v_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_v),nrow=length(species_list)))
reads_table_r_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_r),nrow=length(species_list)))
reads_table_b_new <- as.data.frame(matrix(data=0,ncol=ncol(reads_table_b),nrow=length(species_list)))

row.names(reads_table_u_new) = species_list
colnames(reads_table_u_new) = colnames(reads_table_u)
row.names(reads_table_v_new) = species_list
colnames(reads_table_v_new) = colnames(reads_table_v)
row.names(reads_table_r_new) = species_list
colnames(reads_table_r_new) = colnames(reads_table_r)
row.names(reads_table_b_new) = species_list
colnames(reads_table_b_new) = colnames(reads_table_b)

keep <- species_list %in% row.names(reads_table_v)
reads_table_v_new[keep,] = reads_table_v
keep <- species_list %in% row.names(reads_table_u)
reads_table_u_new[keep,] = reads_table_u
keep <- species_list %in% row.names(reads_table_r)
reads_table_r_new[keep,] = reads_table_r
keep <- species_list %in% row.names(reads_table_b)
reads_table_b_new[keep,] = reads_table_b

reads_table <- cbind(reads_table_v_new,reads_table_u_new,reads_table_r_new,reads_table_b_new)
write.csv(reads_table,'four_body_sites.csv')
##### rarefaction ##################################################

reads_table_rare <- t(reads_table)
reads_table_rare = Rarefy(reads_table_rare, depth = min(rowSums(reads_table_rare)))
reads_table_rare <- reads_table_rare$otu.tab.rff
reads_table_rare <- as.data.frame(reads_table_rare)

##### diversity ##################################################
# Bray_Curtis distance
Bray_Curtis <- as.matrix(vegdist(reads_table_rare, METHODS="bray", binary=FALSE))
Bray_Curtis <- as.data.frame(Bray_Curtis)

# Running Nonmetric Multidimensional Scaling (NMDS) Ordination
NMDS <-
   metaMDS(Bray_Curtis,
           distance = "bray",
           k = 2,
           maxit = 999, 
           trymax = 50,
           wascores = TRUE)

mds_data <- as.data.frame(NMDS$points)
mds_data$Site <- metadata_all

mds_data$Site = factor(mds_data$Site, levels = c('Urine','Rectum','Vagina','Cheek'))

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Site)) +
   geom_point()+ theme_bw()
ggsave("NMDS_all.pdf", width=5.5, height=5.5)

metadata = c(rep('Vagina',39),rep('Urine',39),rep('Rectum',39),rep('Cheek',39))
metadata = as.data.frame(metadata)
result = beta_diversity(reads_table, metadata = metadata, factor_name = metadata, order = NA, NMDS_skip = T, ref_group = NA, rarefy_to = NA, pheatmap_fontsize = 50,treeheight = 10, pheatmap_y = F)
result$group_dis_sig
result$group_dis_2_p


############################### get metadata #########################
sample_name <- unique(metadata_v$`metadata_v$V1`)
species_name <- unique(metadata_v$`str_c(metadata_v$V2, metadata_v$V3, sep = " ")`)

reads_table_v <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table_v) = species_name 
colnames(reads_table_v) = sample_name

for (a in 1: dim(metadata_v)[1]) {
   
   column_num <- which(sample_name == metadata_v[a,1])
   row_num <- which(species_name == metadata_v[a,2])
   
   reads_table_v[row_num,column_num] = as.numeric(as.character(metadata_v[a,3]))
   
}
reads_table_v = as.data.frame(reads_table_v)


sample_name <- unique(metadata_u$`metadata_u$V1`)
species_name <- unique(metadata_u$`str_c(metadata_u$V2, metadata_u$V3, sep = " ")`)

reads_table_u <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table_u) = species_name 
colnames(reads_table_u) = sample_name

for (a in 1: dim(metadata_u)[1]) {
   
   column_num <- which(sample_name == metadata_u[a,1])
   row_num <- which(species_name == metadata_u[a,2])
   
   reads_table_u[row_num,column_num] = as.numeric(as.character(metadata_u[a,3]))
}

reads_table_u = as.data.frame(reads_table_u)

match_uv$V1 = str_replace_all(match_uv$V1,'MU1W','MU1D')
reads_table_u = reads_table_u[match_uv$V1]
reads_table_v = reads_table_v[match_uv$V2]

##### sample total reads threshold ####################
keep <- matrix(ncol = dim(reads_table_v)[2],)
for (a in 1: dim(reads_table_v)[2]) {
   keep[,a] = sum(reads_table_v[,a]) >=10000    # input
}
reads_table_u <- reads_table_u[,keep]
reads_table_v <- reads_table_v[,keep]

keep <- matrix(ncol = dim(reads_table_u)[2],)
for (a in 1: dim(reads_table_u)[2]) {
   keep[,a] = sum(reads_table_u[,a]) >=10000    # input
}

reads_table_u <- reads_table_u[,keep]
reads_table_v <- reads_table_v[,keep]

##### get metadata #####
keep = metadata$sampleID %in% colnames(reads_table_v)

metadata = metadata[keep,]

mean(metadata$ga_at_col)
sd(metadata$ga_at_col)


mean(metadata2$V2)
sd(metadata2$V2)



##### empty control #####
setwd('/Users/binzhu/secure/godel/gpfs_fs/momspi/16s/M16-001/stirrups')

file_list <- list.files(path = ".", pattern = '*Water*')
file_list <- file_list[str_detect(file_list,'hq_summary_97.txt')]


trials = c(1: length(file_list))

func_1 = function(trial) {
   if (file.info(file_list[trial])$size == 0) {
      return(NA)
   }
   
   reads_table = read.table(file_list[trial], sep = '|')
   
   reads_table$V3[reads_table$V3 == 'AT'] = ''
   reads_table$V3 = paste(reads_table$V2,reads_table$V3,sep='_')
   reads_table[,c(2,5,6)] = NULL

   return(reads_table)
   
}

data_3 = mclapply(trials, func_1, mc.cores = 8)

setwd('/Users/binzhu/Desktop/Urine')

data = as.data.frame(matrix(data = NA, ncol = 3, nrow = 0))
colnames(data) = c('SampleID','Taxa','No_of_Reads')

for (a in 1: length(data_3)) {
   data_2 = data_3[a]
   data_2 <- unlist(data_2)
   
   data_4 = as.data.frame(matrix(data = NA, ncol = 3, nrow = length(data_2)/3))
   colnames(data_4) = c('SampleID','Taxa','No_of_Reads')
   data_4[,1] = as.character(data_2[1:(length(data_2)/3)])
   data_4[,2] = as.character(data_2[((length(data_2)/3)+1):((length(data_2)/3)*2)])
   data_4[,3] = as.character(data_2[((length(data_2)/3)*2+1):((length(data_2)/3)*3)])

   data = rbind(data,data_4)
}

colnames(data) = c('V1','V2','V4')
sample_name <- unique(data$V1)
species_name <- unique(data$V2)

data$V2 <- as.character(data$V2)
data$V1 <- as.character(data$V1)
data$V4 <- as.numeric(as.character(data$V4))

reads_table <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table) = species_name 
colnames(reads_table) = sample_name

for (a in 1: nrow(data)) {
   
   column_num <- which(sample_name == data[a,1])
   row_num <- which(species_name == data[a,2])
   
   reads_table[row_num,column_num] =  data$V4[a]
   
}
reads_table <- as.data.frame(reads_table)
colnames(reads_table) = str_remove_all(colnames(reads_table) ,'M16-001-')
write.csv(reads_table,'Blank_control.csv')


















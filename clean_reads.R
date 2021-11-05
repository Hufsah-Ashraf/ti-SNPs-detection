#Since we see reads that don't agree to the genotype sometimes e.g. inverted reads from a sample that was '0|0'
#we need to additionally account for the genotypes to get rid of such false signals. 
#i.e. we accept 'forward' reads from a sample only if atleast one of it's haplotypes is 'REF'.
#Similarly,  we accept 'inverted' reads from a sample only if atleast one of it's haplotypes is reported to be 'INV'
library(dplyr)
library(data.table)
library(stringr)

files_link<-'per_sample_configs/'
files <- list.files(files_link,pattern = "txt")
for (i in 1:length(files)){
  files[i]<-paste0(files_link,files[i])
  
}
genotypes<-data.frame(fread('variants_freeze4inv_sv_inv_hg38_processed_arbigent_filtered_manualDotplot_filtered_PAVgenAdded_withInvCategs.tsv'))
genotypes[,15:58]<- data.frame(lapply(genotypes[,15:58], function(x) { gsub("_lowconf", "", x)  }))
genotypes[,15:58]<- data.frame(lapply(genotypes[,15:58], function(x) { gsub("noreads", "\\.\\|\\.", x)  }))
genotypes[,15:58]<- data.frame(lapply(genotypes[,15:58], function(x) { gsub("2020", "0\\|0", x)  }))
genotypes[,15:58]<- data.frame(lapply(genotypes[,15:58], function(x) { gsub("\\/", "\\|", x)  }))
#genotypes[,15:58]<- data.frame(lapply(genotypes[,15:58], function(x) { gsub("1\\/1", "1\\|1", x)  }))

colnames(genotypes)[6]<-'ID'
for (f in (1:length(files))) {
  link<-files[f]
  config_file<-data.table(fread(link))
  sample = str_match(link, "([HG|NA|GM]+[0-9]{3,5})")[,2]
  sample=ifelse(sample=='GM20509','NA20509',ifelse(sample=='GM20847','NA20847', 
                      ifelse(sample=='GM19983','NA19983',ifelse(sample=='HG002','NA24385',
                    ifelse(sample=='GM19650','NA19650', ifelse(sample=='GM19434', 'NA19434', 
                   ifelse(sample=='GM19036', 'NA19036', ifelse(sample=='GM18939', 'NA18939', 
                    ifelse(sample=='GM18534', 'NA18534', ifelse(sample=='GM12329', 'NA12329', 
                            ifelse(sample=='GM19434', 'NA19434', sample)))))))))) )
  sample_genotypes<-data.frame(genotypes$seqnames, genotypes$start, genotypes$end, genotypes$ID, genotypes[,sample])
  colnames(sample_genotypes)<-c('chrom','inv_start','inv_end','ID','GT')
  config_file<-left_join(config_file, sample_genotypes)
  config_file[['H1']]<-lapply(config_file$GT, function(x) {str_split_fixed(x, pattern = '\\|', n=Inf)[1]})
  config_file[['H2']]<-lapply(config_file$GT, function(x) {str_split_fixed(x, pattern = '\\|', n=Inf)[2]})
  config_file[['H1']]<-as.character(config_file[['H1']])
  config_file[['H2']]<-as.character(config_file[['H2']])
  config_file[['ref_fwd']]<-ifelse(config_file$H1=='0' | config_file$H2=='0', config_file$ref_fwd, 0)
  config_file[['alt_fwd']]<-ifelse(config_file$H1=='0' | config_file$H2=='0', config_file$alt_fwd, 0)
  config_file[['ref_inv']]<-ifelse(config_file$H1=='1' | config_file$H2=='1', config_file$ref_inv, 0)
  config_file[['alt_inv']]<-ifelse(config_file$H1=='1' | config_file$H2=='1', config_file$alt_inv, 0)
  config_file$H1<-NULL
  config_file$H2<-NULL
  config_file$GT<-NULL
  config_file$inv_ID<-NULL
  colnames(config_file)[length(colnames(config_file))]<-'inv_ID'
  config_file<-config_file[,c(1:4,13,5:12)]
  config_file<-config_file %>% mutate(ref_fwd = as.numeric(ref_fwd), 
                                      ref_inv = as.numeric(ref_inv),
                                      alt_fwd = as.numeric(alt_fwd),
                                      alt_inv = as.numeric(alt_inv),
                                      inv_ID = as.character(inv_ID))
  filename<-str_split_fixed(link, pattern = '/', n=Inf)[2]
  write.table(config_file, paste0('per_sample_configs_clean/',filename), col.names = T, row.names = F, quote = F)
  
}

##if needed
#files_link<-'per_sample_configs_clean/'
#files <- list.files(files_link,pattern = "txt")
#for (i in 1:length(files)){
#  files[i]<-paste0(files_link,files[i])
#  
#}
#for (f in (1:length(files))) {
#  link<-files[f]
#  config_file<-data.table(fread(link))
#  config_file<-config_file %>% mutate(ref_fwd = as.numeric(ref_fwd), 
#                                      ref_inv = as.numeric(ref_inv),
#                                      alt_fwd = as.numeric(alt_fwd),
#                                      alt_inv = as.numeric(alt_inv),
#                                      inv_ID = as.character(inv_ID))
#  write.table(config_file, link, col.names = T, row.names = F, quote = F)
#  
#}

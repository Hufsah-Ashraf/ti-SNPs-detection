#makes one table 'one_table_per_snp.txt' with all haplotype configuration counts for each snp per inversion across all samples for all inversions
library(data.table)
library(dplyr)
files_link<-'per_sample_configs_clean/'
out_link<-'one_table_per_snp.txt'
files <- list.files(files_link,pattern = "txt")
for (i in 1:length(files)){
  files[i]<-paste0(files_link,files[i])

  }

dat <- lapply(files, function(x) data.frame(fread(x)))
func <- function(...){
  df1 = list(...)[[1]]
  df2 = list(...)[[2]]
  xxx = full_join(..., by = c('chrom','inv_start', 'inv_end','inv_ID', 'snp_pos','ref_allele','alt_allele','snp_AF'))
  return(xxx)
}
print('read files')
whole<-Reduce(func, dat)
print('whole formed')
whole[(grep('ref_fwd',names(whole)))][is.na(whole[(grep('ref_fwd',names(whole)))])]<-0
whole[(grep('ref_inv',names(whole)))][is.na(whole[(grep('ref_inv',names(whole)))])]<-0
whole[(grep('alt_fwd',names(whole)))][is.na(whole[(grep('alt_fwd',names(whole)))])]<-0
whole[(grep('alt_inv',names(whole)))][is.na(whole[(grep('alt_inv',names(whole)))])]<-0
whole$ref_fwd_sum<- rowSums(whole[grep('ref_fwd',names(whole))])
whole$ref_inv_sum<- rowSums(whole[grep('ref_inv',names(whole))])
whole$alt_fwd_sum<- rowSums(whole[grep('alt_fwd',names(whole))])
whole$alt_inv_sum<- rowSums(whole[grep('alt_inv',names(whole))])  
#if whole is already written
#whole<-data.frame(fread('whole.txt'))
one_table_per_snp<-data.frame(whole$chrom, whole$inv_start, whole$inv_end, whole$inv_ID, whole$snp_pos, whole$ref_allele, whole$alt_allele, whole$snp_AF,
                              whole$ref_fwd_sum, whole$ref_inv_sum, whole$alt_fwd_sum, whole$alt_inv_sum)
colnames(one_table_per_snp)<-c('chrom','inv_start', 'inv_end','inv_ID','snp_pos','ref_allele','alt_allele','snp_AF','reads_ref_fwd','reads_ref_inv',
                               'reads_alt_fwd','reads_alt_inv')
one_table_per_snp[['minimum']]<-apply(one_table_per_snp[,9:12],1, FUN = min)
#only keep SNPs within required AF range
one_table_per_snp<-data.frame(one_table_per_snp[one_table_per_snp$snp_AF>=0.05 & one_table_per_snp$snp_AF<=0.95,])

for (i in 1:length(one_table_per_snp$chrom)) {
  one_table_per_snp[i,'signal']<-ifelse((one_table_per_snp[i,'reads_alt_fwd']>2 & one_table_per_snp[i,'reads_alt_inv']>2 & 
                                           one_table_per_snp[i,'reads_ref_fwd']>2 & one_table_per_snp[i,'reads_ref_inv']>2), 'YES', 'NO')
  
}
write.table(one_table_per_snp, out_link, row.names = F, col.names = T, quote = F)

#joining all snps for one inversion
sig_nosig_invs<-one_table_per_snp%>% group_by(chrom, inv_start, inv_end, inv_ID)%>% 
  summarize(snps_showing_signal=sum(signal=='YES'),  considered_snps= n())
sig_nosig_invs$inv_len<-(sig_nosig_invs$inv_end-sig_nosig_invs$inv_start)/1000
# seperate signalling and non-signalling inversions
sig_invs<-sig_nosig_invs[sig_nosig_invs$snps_showing_signal>0,]
non_sig_invs<-sig_nosig_invs[sig_nosig_invs$snps_showing_signal==0,]
write.table(sig_nosig_invs, 'sig_nosig_invs.txt', col.names = T, row.names = F, quote = F)
write.table(sig_invs, 'sig_invs.txt', col.names = T, row.names = F, quote = F)















































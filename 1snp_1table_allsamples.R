#makes one table with all haplotype configuration counts for each snp per inversion across all samples for all inversions
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

#minimum cell entry and sum of all cells for each snp within an inversion
#for (i in 1:length(one_table_per_snp$chrom)){
#  one_table_per_snp[i,'minimum']<-min(one_table_per_snp[i,'reads_ref_fwd'],one_table_per_snp[i,'reads_ref_inv'],
#                                      one_table_per_snp[i,'reads_alt_fwd'], one_table_per_snp[i,'reads_alt_inv'])
#  one_table_per_snp[i,'sum']<-sum(one_table_per_snp[i,'reads_ref_fwd'],one_table_per_snp[i,'reads_ref_inv'],
#                                  one_table_per_snp[i,'reads_alt_fwd'], one_table_per_snp[i,'reads_alt_inv'])
#  
#}
#one table for both signal showing and not showing inversions
write.table(one_table_per_snp, 'one_table_per_snp_full.txt', row.names = F, col.names = T, quote = F)

one_table_per_snp[['minimum']]<-apply(one_table_per_snp[,9:12],1, FUN = min)

one_table_per_snp<-data.frame(one_table_per_snp[one_table_per_snp$snp_AF>=0.05 & one_table_per_snp$snp_AF<=0.95,])

for (i in 1:length(one_table_per_snp$chrom)) {
  one_table_per_snp[i,'signal']<-ifelse((one_table_per_snp[i,'reads_alt_fwd']>2 & one_table_per_snp[i,'reads_alt_inv']>2 & 
                                           one_table_per_snp[i,'reads_ref_fwd']>2 & one_table_per_snp[i,'reads_ref_inv']>2), 'YES', 'NO')
  
}
write.table(one_table_per_snp, out_link, row.names = F, col.names = T, quote = F)

#for (i in 1:length(one_table_per_snp$chrom)) {
#  one_table_per_snp[i,'signal_v2']<-ifelse((one_table_per_snp[i,'reads_alt_fwd']>2 & one_table_per_snp[i,'reads_alt_inv']>2 & 
#                                           one_table_per_snp[i,'reads_ref_fwd']>2 & one_table_per_snp[i,'reads_ref_inv']>2)&
#                                           one_table_per_snp[i,'mapable']>0, 'YES', 'NO')
#  
#}





#POST PROCESSING
#library(dplyr)
# whole table with each snp 
#one_table_per_snp<-read.table('one_table_per_snp.txt', header=T)
#remove snps with AF not between 0.05-0.95
#one_table_per_snp<-data.frame(one_table_per_snp[one_table_per_snp$snp_AF>=0.05 & one_table_per_snp$snp_AF<=0.95,])
# table with total snps lying in each region alongwith the snps with AFbetween 20-80%
#total_snps<-read.table('total_snps.txt', header = T)
#joining all snps for one inversions
sig_nosig_invs<-one_table_per_snp%>% group_by(chrom, inv_start, inv_end, inv_ID)%>% 
  summarize(snps_showing_signal=sum(signal=='YES'),  considered_snps= n())
#add total snps to it
#sig_nosig_invs<-left_join(sig_nosig_invs, total_snps, by=c('chrom','inv_start','inv_end','inv_ID'))
sig_nosig_invs$inv_len<-(sig_nosig_invs$inv_end-sig_nosig_invs$inv_start)/1000
# seperate signalling and non-signalling inversions
sig_invs<-sig_nosig_invs[sig_nosig_invs$snps_showing_signal>0,]
non_sig_invs<-sig_nosig_invs[sig_nosig_invs$snps_showing_signal==0,]
write.table(sig_nosig_invs, 'sig_nosig_invs.txt', col.names = T, row.names = F, quote = F)
write.table(sig_invs, 'sig_invs.txt', col.names = T, row.names = F, quote = F)




###Binomial model 
#one_table_per_snp<-read.table('one_table_per_snp.txt', header=T)
#one_table_per_snp<-one_table_per_snp[one_table_per_snp$snp_AF>=0.05 & one_table_per_snp$snp_AF<=0.95,]
# with n being the total number of reads, p being 0.05 (for 5% background reads) check the probability for observing
# each minimum entry in the table
#for (i in 1:length(one_table_per_snp$chrom)) {
#  if(one_table_per_snp[i,'minimum']>0){
#    one_table_per_snp[i,'p_value']<-1-pbinom((one_table_per_snp[i,'minimum']-1), one_table_per_snp[i,'sum'],0.05)
#  }
#  else{
#    one_table_per_snp[i,'p_value']<-NA
#    
#  }
#}
#all those where p value <0.05, consider them as a signal
#for (i in 1:length(one_table_per_snp$chrom)) {
#  one_table_per_snp[i,'signal_p']<-ifelse(!is.na(one_table_per_snp[i,'p_value'])&one_table_per_snp[i,'p_value']<0.05,'YES', 'NO')
  
#}














































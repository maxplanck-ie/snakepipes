

count_comp_11wks_genomic=read.csv("/data/pospisilik/group/heyne/laura/Project_A564+A583/B6_11wks_all_genomic.coutc_countb.tsv",header = F,sep = "\t",stringsAsFactors = F)

ggplot(dat=count_comp_11wks_genomic,aes(x=(V3),y=(V6),color=V1))+geom_point()+facet_wrap(~V1)+xlab("non-collapsed UMIs per cell")+ylab("unique UMIs per cell (trancripts)")
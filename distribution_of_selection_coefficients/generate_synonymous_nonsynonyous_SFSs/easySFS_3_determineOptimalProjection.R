############ Plotting projection preview syn and mis ###################
wd="/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/easySFS/projectionPreview/"
categories=c("synonymous","missense")

for(category in categories){
preview <- read.table(paste(wd,"GoC.",category,".easySFS.projPreview.R.format.txt",sep=""),header=T,sep=",")

# get max value:
maximum <- preview %>%
  filter(snps==max(snps))


 p1 <- ggplot(preview,aes(x=projection,y=snps))+
  geom_point()+
  theme_bw()+
  geom_point(data=maximum,aes(x=projection,y=snps),color="red",shape=8,size=4)+
  scale_x_continuous(breaks=seq(2,80,2))+
  ggtitle(paste(category," sites SFS projection preview\nmaximum at ",maximum$projection," with ",maximum$snps," sites\n",sep=""))
p1
ggsave(paste(wd,category,".projectionPreview.pdf",sep=""),p1,height = 5,width=7)
}

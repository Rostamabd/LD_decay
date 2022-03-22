rm(list=ls())
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

### read the ld data for four lines
filenames = list.files(pattern="dist_*")
datalist = lapply(filenames, function(x){fread(file=x,header=T)})

### merge multiple files based on SNP_A and SNP_B
### a function to merge all files at once
all_data <- Reduce(function(x,y) {merge(x,y,by=c("SNP_A","SNP_B","CHR_A","CHR_B","BP_A","BP_B","0"))}, datalist)
colnames(all_data) <- c("SNP_A","SNP_B", "CHR_A", "CHR_B", "BP_A","BP_B","distance","line_1", "line_2")

### order based on chr and BP_A
all_data <- all_data[with(all_data, order(CHR_A, BP_A)),]
all_data <- all_data[all_data$CHR_A!=0,]

### get only the distance and r2
ld_summary <- all_data[,c("distance","line_1", "line_2")]

### mean of r2 
ld_summary$distance <- cut(ld_summary$distance,breaks=seq(from=min(ld_summary$distance)-1,to=max(ld_summary$distance)+1,by=100000))


dfr1 <- ld_summary %>% group_by(distance) %>% summarise(mean=mean(line_1),median=median(line_1))

dfr2 <- ld_summary %>% group_by(distance) %>% summarise(mean=mean(line_2),median=median(line_2))

### get the mid points for line_1

dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distance,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distance,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

### get mid points for line_2

dfr2 <- dfr2 %>% mutate(start=as.integer(str_extract(str_replace_all(distance,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distance,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))


#### combine all lines LD values
rsqure_all <- cbind(dfr1[,c("mid","mean")],dfr2[,c("mean")])

colnames(rsqure_all) <- c("dist","line_1","line_2")

### melt the data
df <- melt(rsqure_all ,  id.vars = 'dist', variable.name = "lines")


#######################################################################
#### ld decay plot, overlay plots
#######################################################################
pdf("ld_all_lines.pdf")
ggplot(df,aes(dist,value, col=lines)) + geom_line(aes(color = lines, group = lines)) +
 labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+
scale_x_continuous(breaks=c(0,1*10^6,2*10^6,3*10^6,4*10^6,5*10^6,6*10^6,7*10^6,8*10^6,9*10^6,10^7),
labels=c("0","1","2","3","4","5","6","7","8","9","10"))+
theme(legend.position="none") +theme_minimal()

dev.off()

#####################################################################

### LD for 0 to 100kb

df_100kb <- df[df$dist<=10^5,]

pdf("ld_all_lines_100kb.pdf")
ggplot(df_100kb,aes(dist,value, col=lines)) + geom_line(aes(color = lines, group = lines)) +
 labs(x="Distance (Kilobases)",y=expression(LD~(r^{2})))+
scale_x_continuous(breaks=seq(0,10^5,10000),labels=c("0","10","20","30","40","50","60","70","80","90","100"))+
theme(legend.position="none") +theme_minimal()

dev.off()

#########################################

#### now plot the LD per each population

##########################################
### plot for line_1
pdf("line_1_ld_decay.pdf")
  ggplot()+
  geom_point(data=dfr1,aes(x=start,y=mean),size=0.4,colour="grey20")+
  geom_line(data=dfr1,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
  labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0,2*10^6,4*10^6,6*10^6,8*10^6),labels=c("0","2","4","6","8"))+
  theme_bw()
dev.off()

### plot for line_2
pdf("line_2_ld_decay.pdf")
  ggplot()+
  geom_point(data=dfr2,aes(x=start,y=mean),size=0.4,colour="grey20")+
  geom_line(data=dfr1,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
  labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0,2*10^6,4*10^6,6*10^6,8*10^6),labels=c("0","2","4","6","8"))+
  theme_bw()
dev.off()


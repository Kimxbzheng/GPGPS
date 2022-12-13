setwd("/mnt/md0/wuhn/data")
### 读取表达谱 清空无关变量
coding_gene <- read.csv("protein-coding_gene.txt", header = T, sep = "\t", stringsAsFactors = F)
path <- "./EXP/"
files <- list.files(path = path, pattern = "*.csv")
for (i in 1:length(files)) {
    file_path <- paste(path, files[i], sep = "")
    temp = read.csv(file_path, header = T, stringsAsFactors = F)
    filenames <- strsplit(files[i], "[.]")[[1]][1]
    assign(filenames, temp)
}
rm(path, files, temp, file_path, filenames, i)
#setwd("/data/users/wuhaonan/ex/glioma_data")
### 读取临床信息 清空无关变量
path <- "./Clinic/"
files <- list.files(path = path, pattern = "*.csv")
for (i in 1:length(files)) {
    file_path <- paste(path, files[i], sep = "")
    temp = read.csv(file_path, header = TRUE, stringsAsFactors = F)
    filenames <- strsplit(files[i], "[.]")[[1]][1]
    assign(filenames, temp)
}

rm(path, files, temp, file_path, filenames, i)

###
library(Metrics)
library(Hmisc)
library(glmnet)
library(snowfall)

###二项分布筛选 稳定对
Stable_pair_Binomial<-function(Mat,cutoff){
  GeneID<-rownames(Mat)
  Nsample<-ncol(Mat)
  Result<-matrix(NA,ncol=2)
  for( i in 1:nrow(Mat)){
    CopyMat<-matrix(rep(Mat[i,],times=(nrow(Mat)-i)),ncol=ncol(Mat),byrow=T)   #扩充矩阵
    IF_Positive<-CopyMat-Mat[-(1:i),]  #A:第i行的基因
    # print(IF_Positive)
    IF_Positive[which(IF_Positive<0)]<- 0   #
    IF_Positive[which(IF_Positive>0)]<- 1
    Sum_AmB<-rowSums(IF_Positive)             #表达值A大于B的样本个数
    P<-matrix()
    
    for(j in 1:length(Sum_AmB)){
      P[j]<-1-pbinom(Sum_AmB[j],Nsample,0.5)
    }
    P<-as.matrix(p.adjust(P,'BH'))
    Gene<-GeneID[-(1:i)]
    
    
    TempB<-as.matrix(Gene[which(P<cutoff)])
    
    if(length(TempB)!=0){
      TempA<-as.matrix(rep(GeneID[i],times=(length(TempB))))
      R_Mat_p1<-cbind(TempA,TempB)  #B<A
    }else{
      R_Mat_p1<-matrix(NA,ncol=2)
    }      
    
    TempA<-as.matrix(Gene[which(P>(1-cutoff))])   #A>B
    
    if(length(TempA)!=0){
      TempB<-as.matrix(rep(GeneID[i],times=(length(TempA))))
      R_Mat_p2<-cbind(TempA,TempB)
    }else{
      R_Mat_p2<-matrix(NA,ncol=2)
    }   
    
    R_Mat<-rbind(R_Mat_p1,R_Mat_p2)
    Result<-rbind(Result,R_Mat)
    
  }    
  return(na.omit(Result))   #A>B
}

#########################################################################################################################
#卡阈值找稳定对
Stable_pair2<-function(Mat,cutoff){#,file1,file2,file3
  GeneID<-rownames(Mat)
  #print(head(GeneID))
  Result<-matrix(NA,ncol=2)
  #print(head(Result))
  for( i in 1:nrow(Mat)){
    CopyMat<-matrix(ncol=ncol(Mat),rep(Mat[i,],times=(nrow(Mat)-i)),byrow=T)   #扩充矩阵
    IF_Positive<-CopyMat-Mat[-(1:i),]  #A:第i行的基因
    IF_Positive[which(IF_Positive<0)]<- 0   #
    IF_Positive[which(IF_Positive>0)]<- 1
    Sum_AmB<-rowSums(IF_Positive) #表达值A大于B的样本个数
    Gene<-GeneID[-(1:i)]
    
    N<-ncol(Mat)
    
    TempB<-as.matrix(Gene[which(Sum_AmB>(N*cutoff))])
    
    if(length(TempB)!=0){
      TempA<-as.matrix(rep(GeneID[i],times=(length(TempB))))
      R_Mat_p1<-cbind(TempA,TempB)  #B<A
    }else{
      R_Mat_p1<-matrix(NA,ncol=2)
    }      
    TempA<-as.matrix(Gene[which(Sum_AmB<(N*(1-cutoff)))])   #A>B
    
    if(length(TempA)!=0){
      TempB<-as.matrix(rep(GeneID[i],times=(length(TempA))))
      R_Mat_p2<-cbind(TempA,TempB)
    }else{
      R_Mat_p2<-matrix(NA,ncol=2)
    }   
    
    R_Mat<-rbind(R_Mat_p1,R_Mat_p2)
    
    Result<-rbind(Result,R_Mat)
    # print(head(Result))
    
    
  }    
  return(na.omit(Result))   #A>B
}

#############################################################################################################################

#表达谱处理  行 列 名
GetRowNames<-function(EXP){
    

 NAMES<-EXP[,1]
  Cnames<-colnames(EXP)
  EXP<-EXP[,-1]
  rownames(EXP)<-NAMES
  colnames(EXP)<-Cnames[-1]
  return(as.matrix(EXP))
}
###########################################################################################################################

#####基因ID和numb#######
GeneID_2_Symbol<-function(Result){
    COL1<-as.matrix(coding_gene[match(Result[,1],coding_gene$entrez_id),]$symbol)
    COL2<-as.matrix(coding_gene[match(Result[,2],coding_gene$entrez_id),]$symbol)
    TABLE<-cbind(COL1,Result,COL2)
    colnames(TABLE)<-c("Gene Symbol","Gene ID","Gene ID","Gene Symbol")
    return(TABLE)
    }

########################################################################################################################
#获取秩次矩阵
 GetRank<-function(EXP_rank){ 
    
  N<-ncol(EXP_rank)
  NR<-nrow(EXP_rank)
  for(i in (1:NR)){
    Temp<-length(which(EXP_rank[i,]==0))     ##计算0值数量
    if(Temp>(N/2)){
      EXP_rank[i,]<-NA
    }
  }
  EXP_rank<-na.omit(EXP_rank)           #删除0值多的基因
  NewMat<-rownames(EXP_rank)
  #print(head(NewMat))###基因ID
  for (i in (1:N)){
    Temp1<-unique(as.numeric(EXP_rank[,i]))    #样本所以基因表达值 去重
    numb<-length(Temp1)                        #去重后还剩多少基因
    OrderSqe<-matrix(data=(1:numb),ncol=1)     #生成序列数列
    EXP_Drd<-cbind(as.matrix(sort(Temp1)),OrderSqe)    #排序
    Rank<-EXP_Drd[match(EXP_rank[,i],EXP_Drd[,1]),2]   #匹配
    NewMat<-cbind(NewMat,Rank)             #得到秩序[]
    
  }
  NewMat<-GetRowNames(NewMat)
  
  colnames(NewMat)<-colnames(EXP_rank)
  return(NewMat)   #返回表达谱的基因秩序  S*N
}
##取中位秩次
MedRank<-function(RankMat){
  Temp<-as.matrix(apply(RankMat,1,function(x){median(as.numeric(x))}))
  return(Temp)
}
###########################################################################################################################
#筛选差异基因
SelectDEG<-function(RankA,RankB,n){      ##n为选取秩序差异最大的2n个基因   rank n*2
    
    GeneID_A<-as.matrix(as.numeric(rownames(RankA)))
    GeneID_B<-as.matrix(as.numeric(rownames(RankB)))
    RankA<-cbind(GeneID_A,RankA)
    RankB<-cbind(GeneID_B,RankB)           #将行名改为第一列
    #print("n of RA")
    #print(length(RankA))
    #print(head(RankA))
    #print("n of RB")
    #print(length(RankB))
    #print(head(RankB))
    
    InterGeneID<-intersect(GeneID_A,GeneID_B)##交叠的基因
    #print("ING")
    #print(head(InterGeneID))
    RankA<-RankA[match(InterGeneID,RankA[,1]),]
    RankB<-RankB[match(InterGeneID,RankB[,1]),]          ##取交叠基因的秩序矩阵

    Def_Rank<-as.matrix(RankA[,2]-RankB[,2])             ##秩序差   ####   eg.LGG-GBM   疾病-正常
    rownames(Def_Rank)<-InterGeneID
    
    A<-names(Def_Rank[order(as.numeric(Def_Rank),decreasing = F),])[c(1:n)] 
    C<-names(Def_Rank[order(as.numeric(Def_Rank),decreasing = T),])[c(1:n)] 
   
    DEG<-list()
    if(length(intersect(A,C))==0){
      DEG[[1]]<-A
      DEG[[2]]<-C
    }else{
      Overlop<-intersect(A,C)
      DEG[[1]]<-setdiff(A,Overlop)      #去掉交叠的基因
      DEG[[2]]<-setdiff(C,Overlop)
    }
    
    DEG[[3]]<-unique(as.matrix(c(A,C)))
    print(dim(DEG[[3]]))
    names(DEG)<-c("Gene1","Gene2","All_DEG")
    #print("n of All_DEG")
    
    return(DEG)
    
  }

############################################################################################################################
    
GetDEG<-function(EXP,Clinic,n){  #clinic 统一通过lable中的0和1判断   输入的EXP是有行名
    
    TypeA<-Clinic[which(Clinic$label==1),1]
    TypeB<-Clinic[which(Clinic$label==2),1]
    print("Sample of TypeA")         #### A类名称
    print(length(TypeA))
    print("Sample of TypeB")         #### B类名称
    print(length(TypeB))
    EXP<-GetRowNames(EXP)
     
    EXP_TypeA<-EXP[,match(TypeA,colnames(EXP))]    ####  1对应于原GBM  0对应与LGG
    EXP_TypeB<-EXP[,match(TypeB,colnames(EXP))]

    Rank_TypeB<-GetRank(EXP_TypeB)
    Rank_TypeA<-GetRank(EXP_TypeA)


    MedRank_B<-MedRank(Rank_TypeB)
    MedRank_A<-MedRank(Rank_TypeA)
        
    DEG<-SelectDEG(MedRank_A,MedRank_B,n)
    Result<-list()
    Result[[1]]<-DEG[[1]]
    Result[[2]]<-DEG[[2]]
    Result[[3]]<-DEG[[3]]
    Result[[4]]<-EXP_TypeA
    Result[[5]]<-EXP_TypeB
        
    return(Result)
    }
##########################################################################################################################
## 获取差异基因表达谱
GetDEG_exp<-function(DEG,EXP,RowN = T)
{
  if(RowN==T){
    GeneID<-as.matrix(as.numeric(rownames(EXP)))
  }else{
    GeneID<-as.matrix(as.numeric(EXP[,1]))
  }
      
  EXP[match(DEG,GeneID),] 
  
}
###############################################################################################################################

ScoreMat<-function(DEG_infor,EXP){   #DEG_infor为训练出的差异基因信息矩阵
DEG_pair<-DEG_infor[[1]]
DEG<-unique(rbind(as.matrix(DEG_pair[,1]),as.matrix(DEG_pair[,2])))
EXP<-GetRowNames(EXP)
sub_EXP<-GetDEG_exp(as.matrix(as.numeric(DEG)),EXP,RowN = T)
scoreMat<-matrix(NA,nrow=nrow(DEG_pair),ncol=ncol(sub_EXP))

for (i in 1:nrow(DEG_pair)){

        EXP_A_B<-sub_EXP[match(as.matrix(DEG_pair[i,]),DEG),]
        Temp<-t(as.matrix(EXP_A_B[1,]-EXP_A_B[2,]))   #相减  A-B
        Temp[which(Temp>0)]<-1     #A>B
        Temp[which(Temp<0)]<- -1   #A<B
        scoreMat[i,]<-Temp
    }
   
    colnames(scoreMat)<-colnames(sub_EXP)
    
    return(na.omit(scoreMat))
}


##########################################################################################################################

Get_Clinic_label<-function(Clinic){
    Temp<-as.matrix(rep(NA,times=(dim(Clinic)[1])))
    #####自定义####
    Temp[which(Clinic$grade!=4&Clinic$IDH==1)]<-1
    Temp[which(Clinic$grade!=4&Clinic$IDH==0)]<-2
    

    colnames(Temp)<-"label"
    Clinic<-cbind(Clinic,Temp)
    return(Clinic)
} 
Get_Clinic_Slabel<-function(Clinic){
    Temp<-as.matrix(rep(NA,times=(dim(Clinic)[1])))
    #####自定义####
    Temp2<-median(Clinic$survival, na.rm = T)
    Temp[which(Clinic$survival<Temp2)]<-1
    Temp[which(Clinic$survival>Temp2)]<-0
    

    colnames(Temp)<-"GOOD"
    Clinic<-cbind(Clinic,Temp)
    return(Clinic)
} 
ScoreMat2<-function(DEG_pair,EXP){   #DEG_infor为训练出的差异基因信息矩阵
DEG<-unique(rbind(as.matrix(DEG_pair[,1]),as.matrix(DEG_pair[,2])))

sub_EXP<-GetDEG_exp(as.matrix(as.numeric(DEG)),EXP,RowN = T)
scoreMat<-matrix(NA,nrow=nrow(DEG_pair),ncol=ncol(sub_EXP))

for (i in 1:nrow(DEG_pair)){

        EXP_A_B<-sub_EXP[match(as.matrix(DEG_pair[i,]),DEG),]
        Temp<-t(as.matrix(EXP_A_B[1,]-EXP_A_B[2,]))   #相减  A-B
        Temp[which(Temp>0)]<-1     #A>B
        Temp[which(Temp<=0)]<-(-1)   #A<B
        scoreMat[i,]<-Temp
    }
   
    colnames(scoreMat)<-colnames(sub_EXP)
    rownames(scoreMat)<-c(1:nrow(scoreMat))
    return(na.omit(scoreMat))
}



##ihd
ls_c<-list(TCGA_LGG_clinical,CGGA_clinic,GSE16011_clinic,`E-TABM-3892_clinic`,GSE61374_clinic,GSE43113_clinic,GSE43388_clinic)
ls_e<-list(tcga_lgg_exp451,CGGA_exp,GSE16011_exp,`E-TABM-3892_exp`,GSE61374_exp,GSE43113_exp,GSE43388_exp)
name<-c("TCGA","CGGA","GSE16011","E-TABM-3892","GSE61374","GSE43113","GSE43388")
Stable_pair_fisher<-function(Mat1,Mat2,cutoff){
  GeneID<-rownames(Mat1)
  Nsample<-ncol(Mat1)
  Result<-matrix(NA,ncol=2)
  pb <- txtProgressBar(style=3)
  for( i in 1:nrow(Mat1)-1){
    CopyMat<-matrix(rep(Mat1[i,],times=(nrow(Mat1)-i)),ncol=ncol(Mat1),byrow=T)   #扩充矩阵
    IF_Positive<-CopyMat-Mat1[-(1:i),]  #A:第i行的基因
    # print(IF_Positive)
    IF_Positive[which(IF_Positive<0)]<- 0   #
    IF_Positive[which(IF_Positive>0)]<- 1
    Sum_AmB<-rowSums(IF_Positive)             #表达值A大于B的样本个数
   temp<-as.matrix(rep(dim(IF_Positive)[2],times = dim(IF_Positive)[1]))-Sum_AmB
   Sum_AmB<-cbind(Sum_AmB,temp)
    
    CopyMat<-matrix(rep(Mat2[i,],times=(nrow(Mat2)-i)),ncol=ncol(Mat2),byrow=T)   #扩充矩阵
    IF_Positive<-CopyMat-Mat2[-(1:i),]  #A:第i行的基因
    # print(IF_Positive)
    IF_Positive[which(IF_Positive<0)]<- 0   #
    IF_Positive[which(IF_Positive>0)]<- 1
    Sum_AmB2<-rowSums(IF_Positive)             #表达值A大于B的样本个数
   temp<-as.matrix(rep(dim(IF_Positive)[2],times = dim(IF_Positive)[1]))-Sum_AmB2
   Sum_AmB2<-cbind(Sum_AmB2,temp)
      
    SUM_AMB<-cbind(Sum_AmB,Sum_AmB2)
      
    P<-apply(SUM_AMB,1,function(ind){
    Temp<- matrix(c(ind[1], ind[3],ind[2],ind[4]), nrow = 2)
    Temp2<- fisher.test(Temp)
    return(Temp2[1])} )
    P<-as.matrix(p.adjust(as.matrix(unlist(P)),'BH'))
    Gene<-GeneID[-(1:i)]
   TempB<-as.matrix(Gene[which(P<cutoff)])
 if(length(TempB)!=0){
      TempA<-as.matrix(rep(GeneID[i],times=(length(TempB))))
      R_Mat_p1<-cbind(TempA,TempB)  #B<A
     TempA2<-matrix(SUM_AMB[which(P<cutoff),],ncol = 4)

if(sum(TempA2[,1]>TempA2[,2])!=0)
    {R_Mat_p1[which(TempA2[,1]>TempA2[,2]),]<-R_Mat_p1[which(TempA2[,1]>TempA2[,2]),c(2,1)]}
    }
      else{
      R_Mat_p1<-matrix(NA,ncol=2)
    }      
      Result<-rbind(Result,R_Mat_p1)
      setTxtProgressBar(pb, i/nrow(Mat1))
    } 
    return(na.omit(Result))   #A>B
}
Get_Clinic_label<-function(Clinic){
    Temp<-as.matrix(rep(NA,times=(dim(Clinic)[1])))
    #####自定义####
    Temp[which(Clinic$IDH==1)]<-1
    Temp[which(Clinic$IDH==0)]<-2
        colnames(Temp)<-"label"
    Clinic<-cbind(Clinic,Temp)
    return(Clinic)
} 
IDHROC_Polt<-function(ls_e,ls_c,name,Group_Pair){
    modelroc<-list()
     for (i in 1:length(ls_e)){    
    EXP<-GetRowNames(ls_e[[i]])
    Clinic<-ls_c[[i]]
    A<-na.omit(Get_LMtable2(Clinic,EXP,Group_Pair)[[1]])
 
    modelroc[[i]]  <- roc(A$IDH,A$LABEL1,ci=T,auc = T)
    #modelroc[[16]][2]
         names(modelroc)[i]<-name[i]
}
    #print(modelroc)

    g2 <- ggroc(modelroc,size = 1,alpha=0.6,)
print(g2+theme_classic()+ggsci::scale_color_jco()+annotate("text",x = .25, y = .30, ## 注释text的位置
           label = paste(name[1],"AUC =", round(modelroc[[1]][[16]][2], 2)))+
      annotate("text",x = .25, y = .25, ## 注释text的位置
           label = paste(name[2],"AUC =", round(modelroc[[2]][[16]][2], 2)))+
      annotate("text",x = .25, y = .20, ## 注释text的位置
           label = paste(name[3],"AUC =", round(modelroc[[3]][[16]][2], 2)))+
      annotate("text",x = .25, y = .15, ## 注释text的位置
           label = paste(name[4],"AUC =", round(modelroc[[4]][[16]][2], 2)))+
      annotate("text",x = .25, y = .10, ## 注释text的位置
           label = paste(name[5],"AUC =", round(modelroc[[5]][[16]][2], 2)))+
      annotate("text",x = .25, y = .05, ## 注释text的位置
           label = paste(name[6],"AUC =", round(modelroc[[6]][[16]][2], 2)))+
      annotate("text",x = .25, y = .00, ## 注释text的位置
           label = paste(name[7],"AUC =", round(modelroc[[7]][[16]][2], 2)))+
    theme(legend.position = 'top',axis.text.x  = element_text(size=10),axis.text.y  = element_text(size=10),legend.text =element_text(size=10)))
}
Get_LMtable2<-function(Clinic,EXP,Group_Pair){
   # EXP<-GetRowNames(EXP)
    HHH<-list()
Score<-ScoreMat2(Group_Pair,EXP)
colnames(Clinic)[1]<-"patientID"
    
Score<-t(Score)
data<-as.data.frame(Score)
colnames(data)<-1:16
Clinic<-Clinic[match(rownames(Score),Clinic[[1]]),c("patientID","survival","event","IDH","C1P19Q")]
LABEL1<-as.matrix(predict.glm(pre1,type="response",newdata=data))    ####  1为低风险，0为高风险   
Mask<-cbind(Clinic,LABEL1)
HHH[[1]]<-Mask
    return(HHH)
    }
Diff_Co_Exp<-function(EXP,Clinic,n,QQ,K){
    Clinic<-Get_Clinic_label(Clinic)
    
######## 获取差异基因######################
 
    Result<-GetDEG(EXP,Clinic,n)
   # Result[[1]]<-DEG[[1]]    A中的差异基因
   # Result[[2]]<-DEG[[2]]    B中的差异基因
   # Result[[3]]<-DEG[[3]]    AB 共有的差异基因
   # Result[[4]]<-EXP_TypeA   差异前的A表达谱
   # Result[[5]]<-EXP_TypeB   差异前的B表达谱
    
############# 获取差异表达谱###############
  
  Mat1<-GetDEG_exp(Result[[3]],Result[[4]],RowN = T)
  Mat2<-GetDEG_exp(Result[[3]],Result[[5]],RowN = T)
    

 TT2<-matrix(match(QQ[[1]],rownames( Mat1)),ncol = 2)
temp1<-as.matrix(apply(TT2,1,function(X){
    c(rcorr(t(Mat1[X,]),type=c("spearman"))[[1]][2],rcorr(t( Mat1[X,]),type=c("spearman"))$P[2])
}))
temp2<-as.matrix(apply(TT2,1,function(X){
     c(rcorr(t( Mat2[X,]),type=c("spearman"))[[1]][2],rcorr(t( Mat2[X,]),type=c("spearman"))$P[2])
}))
   
  Diff_Co<-rbind(temp1,temp2)

    text<-Diff_Co[,which( Diff_Co[2,]<K& Diff_Co[4,]<K)]
    index<-which((text[1,]>0&text[3,]>0)|(text[1,]<0&text[3,]<0))
   return(index)
}

INGENE<-as.matrix(Reduce(intersect,list(TCGA_GBM_SEQ_exp[[1]],tcga_lgg_exp451[[1]],GSE16011_exp[[1]],GSE61374_exp[[1]],`E-TABM-3892_exp`[[1]],CGGA_exp[[1]])))
tcga_lgg_exp451<-tcga_lgg_exp451[match(INGENE,tcga_lgg_exp451[[1]]),]
GSE16011_exp<-GSE16011_exp[match(INGENE,GSE16011_exp[[1]]),]
#数据整合
clinic1<-na.omit(TCGA_LGG_clinical[,c("patientID","IDH","survival","event")])
colnames(GSE16011_clinic)[1]<-"patientID"
clinic2<-na.omit(GSE16011_clinic[,c("patientID","IDH","survival","event")])
clinic3<-rbind(clinic1,clinic2)
exp1<-tcga_lgg_exp451
exp2<-GSE16011_exp[,-1]
exp3<-cbind(exp1,exp2)
clinicID<-as.matrix(as.character(gsub("_","-",TCGA_LGG_clinical[,1])))

REOF1<-function(EXP,Clinic,n,cutoff,C){
####分类###############
    Clinic<-Get_Clinic_label(Clinic)
    
######## 获取差异基因######################
 
    Result<-GetDEG(EXP,Clinic,n)
   # Result[[1]]<-DEG[[1]]    A中的差异基因
   # Result[[2]]<-DEG[[2]]    B中的差异基因
   # Result[[3]]<-DEG[[3]]    AB 共有的差异基因
   # Result[[4]]<-EXP_TypeA   差异前的A表达谱
   # Result[[5]]<-EXP_TypeB   差异前的B表达谱
    
############# 获取差异表达谱###############
  
  Mat1<-GetDEG_exp(Result[[3]],Result[[4]],RowN = T)
  Mat2<-GetDEG_exp(Result[[3]],Result[[5]],RowN = T)
    
    
############稳定对筛选方法################ 
  ############## Fisher
Stable_pair_fisher<-function(F){
  GeneID<-rownames(Mat1)
  Nsample<-ncol(Mat1)
  Result<-matrix(NA,ncol=3)
  pb <- txtProgressBar(style=3)
  for( i in 1:length(F)){
      i<-F[i]
    CopyMat<-matrix(rep(Mat1[i,],times=(nrow(Mat1)-i)),ncol=ncol(Mat1),byrow=T)   #扩充矩阵
    IF_Positive<-CopyMat-Mat1[-(1:i),]  #A:第i行的基因
    # print(IF_Positive)
    IF_Positive[which(IF_Positive<0)]<- 0   #
    IF_Positive[which(IF_Positive>0)]<- 1
    Sum_AmB<-rowSums(IF_Positive)             #表达值A大于B的样本个数
   temp<-as.matrix(rep(dim(IF_Positive)[2],times = dim(IF_Positive)[1]))-Sum_AmB
   Sum_AmB<-cbind(Sum_AmB,temp)
    
    CopyMat<-matrix(rep(Mat2[i,],times=(nrow(Mat2)-i)),ncol=ncol(Mat2),byrow=T)   #扩充矩阵
    IF_Positive<-CopyMat-Mat2[-(1:i),]  #A:第i行的基因
    # print(IF_Positive)
    IF_Positive[which(IF_Positive<0)]<- 0   #
    IF_Positive[which(IF_Positive>0)]<- 1
    Sum_AmB2<-rowSums(IF_Positive)             #表达值A大于B的样本个数
   temp<-as.matrix(rep(dim(IF_Positive)[2],times = dim(IF_Positive)[1]))-Sum_AmB2
   Sum_AmB2<-cbind(Sum_AmB2,temp)
      
    SUM_AMB<-cbind(Sum_AmB,Sum_AmB2)
      
    P<-apply(SUM_AMB,1,function(ind){
    Temp<- matrix(c(ind[1], ind[3],ind[2],ind[4]), nrow = 2)
    Temp2<- fisher.test(Temp)
    return(Temp2[1])} )
    P<-as.matrix(p.adjust(as.matrix(unlist(P)),'BH'))
    Gene<-GeneID[-(1:i)]
   TempB<-as.matrix(Gene[which(P<cutoff)])
 if(length(TempB)!=0){
      TempA<-as.matrix(rep(GeneID[i],times=(length(TempB))))
      R_Mat_p1<-cbind(cbind(TempA,TempB),P[which(P<cutoff)])  #B<A
     TempA2<-matrix(SUM_AMB[which(P<cutoff),],ncol = 4)

if(sum(TempA2[,1]>TempA2[,2])!=0)
    {R_Mat_p1[which(TempA2[,1]>TempA2[,2]),1:2]<-R_Mat_p1[which(TempA2[,1]>TempA2[,2]),2:1]}
    }
      else{
      R_Mat_p1<-matrix(NA,ncol=3)
    }      
      Result<-rbind(Result,R_Mat_p1)
      setTxtProgressBar(pb, i/nrow(Mat1))
    } 
    return(na.omit(Result))   #A>B
}
Term<-function(n,c){
    
    A<-matrix(c(1:(n-1)))
    B<-rep(c(1:c,c:1),times=(n/c+1))
    rownames(A)<-B[c(1:(n-1))]
    C<-list()
for(i in 1:c){
    C[[i]]<-c(A[which(rownames(A)==i),])
}
return(C)
}

library(snowfall)
    sfStop()
sfInit(parallel = TRUE, cpus = C)
sfExport("Mat1", "Mat2","cutoff") 
result <- sfLapply(Term(length(Result[[3]]),C), Stable_pair_fisher)
sfStop()
RES<-do.call("rbind",result)
  #############   阈值  
     
# R_A<-Stable_pair2(DEG_EXP_A,cutoff)           #Stable_pair2
# R_B<-Stable_pair2(DEG_EXP_B,cutoff)
    
#######################################################    
 
  print("The number of Reverse_Pair")
  print(length(RES)/3)
  RES<-RES[order(as.numeric(RES[,3])),]
##############返回需要的信息#################
    WHN<-list()
    WHN[[1]]<-RES
  return(WHN)  
}


QQ<-REOF1(exp3,clinic3,1500,2.2e-70,10)
#矩阵构建
CLINIC<-clinic3
EXP<-exp3
colnames(CLINIC)[1]<-"patientID"
Clinc_I<-na.omit(CLINIC[match(colnames(EXP),CLINIC[[1]]),c("patientID","IDH")])
EXP_I<-EXP[match(INGENE,EXP[[1]]),match(Clinc_I[[1]],colnames(EXP))]
SID<-c(paste("R",1:dim(Clinc_I)[1],sep = ""))
Clinc_I[[1]]<-SID
colnames(EXP_I)<-SID
EXP_I<-cbind(INGENE,EXP_I)
DEG_pair<-QQ[[1]][1:50,]
DEG<-unique(rbind(as.matrix(DEG_pair[,1]),as.matrix(DEG_pair[,2])))
EXP<-GetRowNames(EXP_I)

sub_EXP<-GetDEG_exp(as.matrix(as.numeric(DEG)),EXP,RowN = T)
Gene_A<-sub_EXP[match(DEG_pair[,1],rownames(sub_EXP)),]
Gene_B<-sub_EXP[match(DEG_pair[,2],rownames(sub_EXP)),]


rownames(Gene_A)<-NULL
rownames(Gene_B)<-NULL
Matrix01<-as.matrix(Gene_B-Gene_A)
Matrix01[which(Matrix01>0)]<-1
Matrix01[which(Matrix01<0)]<-(-1)
#colnames(Matrix01)<-NULL
ind<-Matrix01
Y<-as.matrix(Clinc_I$IDH)
FF<-list()
for(i in 1:50)
{
library(lars)
#set.seed(1)
x = t(ind)
y<-as.integer(Y)
lar1 <-lars(x,y,type = "lasso")
#plot(lar1) 
cv_sol<-cv.lars(x,y,K=80,type="lasso",mode="fraction",use.Gram=FALSE,max.steps=100)
fra1=cv_sol$index[which.min(cv_sol$cv)]
#fra
#lar1$beta[fra,]
coef2 <-coef.lars(lar1,mode="fraction",s=fra1) #s为step+1，也比图2中竖线为2的迭代次数对应，与图3中df值相等；s取值范围1-7.
#QQQ<-predict(lar1,data.frame(t(runif(dim(x)[2],min=0,max=0))),mode="fraction",s=fra1)
#QQQ$fit
#length(coef2[which(coef2!=0)])
FF[[i]]<-which(coef2!=0)
}
FF1<-Reduce(intersect,FF)
DEG_pair<-DEG_pair[FF1,]
DEG<-unique(rbind(as.matrix(DEG_pair[,1]),as.matrix(DEG_pair[,2])))
EXP<-GetRowNames(EXP_I)

sub_EXP<-GetDEG_exp(as.matrix(as.numeric(DEG)),EXP,RowN = T)
Gene_A<-sub_EXP[match(DEG_pair[,1],rownames(sub_EXP)),]
Gene_B<-sub_EXP[match(DEG_pair[,2],rownames(sub_EXP)),]

rownames(Gene_A)<-NULL
rownames(Gene_B)<-NULL
Matrix01<-as.matrix(Gene_B-Gene_A)
Matrix01[which(Matrix01>0)]<-1
Matrix01[which(Matrix01<0)]<-(-1)
#colnames(Matrix01)<-NULL
ind<-Matrix01
Y<-as.matrix(Clinc_I$IDH)
data<-as.data.frame(cbind(t(ind),Y))
colnames(data)[1:16]<-1:16
###结果保存
pre1<-glm(formula = V17~., data = data, family = binomial(link="logit"))
IDH_pair<-DEG_pair



###1p19q
setwd("/mnt/md0/wuhn/data")
### 读取临床信息 清空无关变量
path <- "./Clinic/"
files <- list.files(path = path, pattern = "*.csv")
for (i in 1:length(files)) {
    file_path <- paste(path, files[i], sep = "")
    
    temp = read.csv(file_path, header = TRUE, stringsAsFactors = F)
    filenames <- strsplit(files[i], "[.]")[[1]][1]
    assign(filenames, temp)
}
rm(path, files, temp, file_path, filenames, i)
#数据整合
clinic1<-na.omit(TCGA_LGG_clinical[,c("patientID","C1P19Q")])
colnames(GSE16011_clinic)[1]<-"patientID"
clinic2<-na.omit(GSE16011_clinic[,c("patientID","C1P19Q")])
clinic3<-rbind(clinic1,clinic2)
exp1<-tcga_lgg_exp451
exp2<-GSE16011_exp[,-1]
exp3<-cbind(exp1,exp2)
Get_Clinic_label<-function(Clinic){
    Temp<-as.matrix(rep(NA,times=(dim(Clinic)[1])))
    #####自定义####
    Temp[which(Clinic$C1P19Q==1)]<-1
    Temp[which(Clinic$C1P19Q==0)]<-2
    

    colnames(Temp)<-"label"
    Clinic<-cbind(Clinic,Temp)
    return(Clinic)
} 

QQ2<-REOF1(exp3,clinic3,1500,2.2e-70,15)
#矩阵构建
CLINIC<-clinic3
EXP<-exp3
#CLINIC<-TCGA_LGG_clinical
#EXP<-tcga_lgg_exp451
Clinc_S<-na.omit(CLINIC[match(colnames(EXP),CLINIC[[1]]),c("patientID","C1P19Q")])

EXP_S<-exp3[match(INGENE,EXP[[1]]),match(Clinc_S[[1]],colnames(EXP))]
SID<-c(paste("R",1:dim(Clinc_S)[1],sep = ""))
Clinc_S[[1]]<-SID
colnames(EXP_S)<-SID
EXP_S<-cbind(INGENE,EXP_S)
DEG_pair<-QQ2[[1]][1:50,]
DEG<-unique(rbind(as.matrix(DEG_pair[,1]),as.matrix(DEG_pair[,2])))
EXP<-GetRowNames(EXP_S)

sub_EXP<-GetDEG_exp(as.matrix(as.numeric(DEG)),EXP,RowN = T)
Gene_A<-sub_EXP[match(DEG_pair[,1],rownames(sub_EXP)),]
Gene_B<-sub_EXP[match(DEG_pair[,2],rownames(sub_EXP)),]

rownames(Gene_A)<-NULL
rownames(Gene_B)<-NULL
Matrix01<-as.matrix(Gene_B-Gene_A)
Matrix01[which(Matrix01>0)]<-1
Matrix01[which(Matrix01<0)]<-(-1)
colnames(Matrix01)<-NULL
ind<-Matrix01
Y<-as.matrix(Clinc_S$C1P19Q)
FF2<-list()
for(i in 1:50)
{library(lars)
x = t(ind)
y<-as.integer(Y)
lar2 <-lars(x,y,type = "lasso")
#plot(lar1) 
cv_sol<-cv.lars(x,y,K=80,type="lasso",mode="fraction",use.Gram=FALSE,max.steps=100)
fra2=cv_sol$index[which.min(cv_sol$cv)]
#fra
#lar1$beta[fra,]
coef2 <-coef.lars(lar2,mode="fraction",s=fra2) #s为step+1，也比图2中竖线为2的迭代次数对应，与图3中df值相等；s取值范围1-7.
#QQQ<-predict(lar2,data.frame(t(runif(dim(x)[2],min=0,max=0))),mode="fraction",s=fra2)
#QQQ$fit
FF2[[i]]<-which(coef2!=0)
}
FF2<-Reduce(intersect,FF2)
DEG_pair<-DEG_pair[FF2,]
DEG<-unique(rbind(as.matrix(DEG_pair[,1]),as.matrix(DEG_pair[,2])))
EXP<-GetRowNames(EXP_S)

sub_EXP<-GetDEG_exp(as.matrix(as.numeric(DEG)),EXP,RowN = T)
Gene_A<-sub_EXP[match(DEG_pair[,1],rownames(sub_EXP)),]
Gene_B<-sub_EXP[match(DEG_pair[,2],rownames(sub_EXP)),]

rownames(Gene_A)<-NULL
rownames(Gene_B)<-NULL
Matrix01<-as.matrix(Gene_B-Gene_A)
Matrix01[which(Matrix01>0)]<-1
Matrix01[which(Matrix01<0)]<-(-1)
colnames(Matrix01)<-NULL
ind<-Matrix01
Y<-as.matrix(Clinc_S$C1P19Q)
DEG_pair



PQROC_Polt<-function(ls_e,ls_c,name,Group_Pair){
    modelroc<-list()
     for (i in 1:length(ls_e)){    
    EXP<-GetRowNames(ls_e[[i]])
    Clinic<-ls_c[[i]]
    A<-na.omit(Get_LMtable2(Clinic,EXP,Group_Pair)[[1]])
 
    modelroc[[i]]  <- roc(A$C1P19Q,A$LABEL1,ci=T,auc = T)
    #modelroc[[16]][2]
         names(modelroc)[i]<-name[i]
}
    #print(modelroc)

    g2 <- ggroc(modelroc,size = 1,alpha=0.6,)
print(g2+theme_classic()+ggsci::scale_color_jco()+annotate("text",x = .25, y = .30, ## 注释text的位置
           label = paste(name[1],"AUC =", round(modelroc[[1]][[16]][2], 2)))+
      annotate("text",x = .25, y = .25, ## 注释text的位置
           label = paste(name[2],"AUC =", round(modelroc[[2]][[16]][2], 2)))+
      annotate("text",x = .25, y = .20, ## 注释text的位置
           label = paste(name[3],"AUC =", round(modelroc[[3]][[16]][2], 2)))+
      annotate("text",x = .25, y = .15, ## 注释text的位置
           label = paste(name[4],"AUC =", round(modelroc[[4]][[16]][2], 2)))+
      annotate("text",x = .25, y = .10, ## 注释text的位置
           label = paste(name[5],"AUC =", round(modelroc[[5]][[16]][2], 2)))+
      annotate("text",x = .25, y = .05, ## 注释text的位置
           label = paste(name[6],"AUC =", round(modelroc[[6]][[16]][2], 2)))+
      annotate("text",x = .25, y = .00, ## 注释text的位置
           label = paste(name[7],"AUC =", round(modelroc[[7]][[16]][2], 2)))+
    theme(legend.position = 'top',axis.text.x  = element_text(size=10),axis.text.y  = element_text(size=10),legend.text =element_text(size=10)))
}
Get_LMtable2<-function(Clinic,EXP,Group_Pair){
   # EXP<-GetRowNames(EXP)
    HHH<-list()
Score<-ScoreMat2(Group_Pair,EXP)
colnames(Clinic)[1]<-"patientID"
    
Score<-t(Score)
data<-as.data.frame(Score)
colnames(data)<-1:26
Clinic<-Clinic[match(rownames(Score),Clinic[[1]]),c("patientID","survival","event","IDH","C1P19Q")]
LABEL1<-as.matrix(predict.glm(pre2,type="response",newdata=data)) 
    Mask<-cbind(Clinic,LABEL1)
HHH[[1]]<-Mask
    return(HHH)
    }

data<-as.data.frame(cbind(t(ind),Y))
colnames(data)[1:26]<-1:26
##结果
pre2<-glm(formula = V27~., data = data, family = binomial(link="logit"))
PQ_pair<-DEG_pair
#PQROC_Polt(ls_e,ls_c,name,PQ_pair)



###组合
Group_Pair3<-rbind(IDH_pair,PQ_pair)
Group_Pair3<-as.data.frame(Group_Pair3,stringsAsFactors = F)
dim(Group_Pair3)
clinic1<-na.omit(TCGA_LGG_clinical[,c("patientID","survival","event","IDH","C1P19Q")])
colnames(GSE16011_clinic)[1]<-"patientID"
clinic2<-na.omit(GSE16011_clinic[,c("patientID","survival","event","IDH","C1P19Q")])
clinic3<-rbind(clinic1,clinic2)
exp1<-tcga_lgg_exp451
exp2<-GSE16011_exp[,-1]
exp3<-GetRowNames(cbind(exp1,exp2))
Clinic<-clinic3
EXP<-exp3
Score<-as.data.frame(ScoreMat2(IDH_pair,EXP))
colnames(Clinic)[1]<-"patientID"  
Score<-t(Score)
Clinic<-Clinic[match(rownames(Score),Clinic[[1]]),c("patientID","survival","event","IDH","C1P19Q")]
Mask<-Clinic
attributeListc<-as.numeric()
index<-1:dim(Group_Pair3)[1]
indexSet<-as.numeric()
oosErrorc<-as.numeric()
Mask<-Mask

for(i in index){
    attTry<-setdiff(index,attributeListc)
    errorList<-as.numeric()
    attTemp<-as.numeric()
    for(ii in attTry){
        attTemp<-append(attTemp,attributeListc)
        attTemp<-append(attTemp,ii)
        Group_Pairqqq<-Group_Pair3[attTemp,]
        #finalGroup_Pair<-Group_Pairqqq
        Score<-as.data.frame(ScoreMat2(Group_Pairqqq,EXP))
       # print(head(Score))
        LABEL3<-as.matrix(colSums(Score))
        combinationMask<-na.omit(cbind(Mask,LABEL3))
        #H<-Make_table3(ls_e,ls_c,name)
        COX<-coxph(Surv(survival, event) ~ LABEL3, data =combinationMask)
       
        #res.cut <- surv_cutpoint(combinationMask, time = "survival", event = "event",variables = c("LABEL3"),progressbar=T)
        #sur.cat <- surv_categorize(res.cut)
        #sdf1<- survdiff(Surv(survival, event) ~ LABEL3, data = Mask)
        #statistic <-log10(1 - pchisq(sdf1$chisq, length(sdf1$n) - 1))
        #statistic<-round(COX,2)
        #statistic<-summary(res.cut)[[2]]
        statistic<-summary(COX)[[14]][1]
        errorList<-append(errorList,statistic)
        attTemp<-as.numeric()
    }
    iBest<-which.max(errorList)
    attributeListc<-append(attributeListc,attTry[iBest])
    oosErrorc<-append(oosErrorc,errorList[iBest])
}

####结果
cat("Best attribute indices: ", attributeListc, "\n","Best attribute names: \n")
plot(oosErrorc,type = "l",xlab = "Number of Attributes",ylab = "statistic",main = "statistic of attributes")
finalGroup_Pair<-Group_Pair3[attributeListc[1:which.max(oosErrorc)],]
finalGroup_Pair
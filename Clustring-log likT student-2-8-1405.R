
##################
# part 1        ##
##################
# library(mvtnorm)
# library(dplyr)
# library(mclust)
#install.packages("usefun")
#library(usefun)


ClusSpaTDLik=function(M,ddf,alfa,pow){

#alfa=0.1;  # alfa is positive
bet=1    # 0<=beta<=2 
#M=140
#ddf=Inf;
#pow=1


mmad2={};adjust2={};
mmad3={};adjust3={};
mmad4={};adjust4={};
mmad5={};adjust5={};


mmad11={};adjust11={};
mmad12={};adjust12={};
mmad13={};adjust13={};
mmad14={};adjust14={};
mmad15={};adjust15={};

for(k in 1:M){

########generate sample from spatial data ######  ???? ???: ????? ???? ?? ?????? ??????? ?????   
###
s=7;
n<-s*s 
sy<-runif(n)  #### ??? ???? 
sx<-runif(n)  #### ??? ????
Loc<-cbind(sx,sy)  ### ?????? ??? ? ??? ????
dh<-as.matrix(dist(Loc))
dh=dh/max(dh)  
################# change the correlation structuer #########
m1=matrix(1:n,ncol=s)
index1= m1[col(m1)<=(s/2+1)& row(m1)<=(s/2+1)&col(m1)+ row(m1)<=s]
index2=m1[col(m1)>(s/2)& row(m1)>(s/2)&col(m1)+ row(m1)>s]
#index=c(index1,index2)
ind=c(1:n)
index3=ind[-c(index1,index2)]
######################################
rho=exp(-alfa*dh^ bet); ### ?????? ??????? ????? ???? ??

mu <- c(rep(0,n))
Sigma=rho  #### matrix correlation
Sigma=rho^pow 
# library(mvtnorm)
Data<-rmvt(1, sigma=Sigma, df=ddf)#
#Data;
############################    ???? ???: ????? ???? ?? ?? ?? ???? ?? ????? ???? ??????? ?? ?? ????? ????? 
Loc1 = matrix(Loc[,1], nrow=7)
Loc2 = matrix(Loc[,2], nrow=7)
DataMat = matrix(Data, nrow=7)
###################################
x = as.vector(Loc1)
y = as.vector(Loc2) 
Data=as.vector(DataMat)
Data[index1] = Data[index1] +rnorm(length(index1),-5, 0.01) 
Data[index2] = Data[index2] +rnorm(length(index2), 0, 0.01) 
Data[index3] = Data[index3] +rnorm(length(index3), 5, 0.01) 

#X=cbind(x,y,Data,)

lab0=c(1:49)
lab0[index1]=1
lab0[index2]=2
lab0[index3]=3
####

######################## Compute H-Clust by Vegan ###################### ???? ???: ?? ??????? ?? ????  ??-???? ???? ???? ?? ????? ?????

################ clustring by May program ############
distData1=dist(Data,method = "euclidean")
hc11=hclust(distData1,method = "complete")
hc12=hclust(distData1,method = "single")
hc13=hclust(distData1,method = "average")
hc14=hclust(distData1,method = "ward.D")
hc15=hclust(distData1,method = "ward.D2")

memb11 =cutree(hc11,k=3)
##
lab={};for(i in 1:3){
clas=mean(Data[memb11==i])
lab=c(lab,clas)
     }

memb11=recode(memb11 ,'1'=rank(lab)[1],'2'=rank(lab)[2],'3'=rank(lab)[3])
##
memb12 =cutree(hc12,k=3)
lab={};for(i in 1:3){
clas=mean(Data[memb12==i])
lab=c(lab,clas)
     }

memb12=recode(memb12 ,'1'=rank(lab)[1],'2'=rank(lab)[2],'3'=rank(lab)[3])
##
memb13 =cutree(hc13,k=3)
lab={};for(i in 1:3){
clas=mean(Data[memb13==i])
lab=c(lab,clas)
     }

memb13=recode(memb13 ,'1'=rank(lab)[1],'2'=rank(lab)[2],'3'=rank(lab)[3])

##
memb14 =cutree(hc14,k=3)
lab={};for(i in 1:3){
clas=mean(Data[memb14==i])
lab=c(lab,clas)
     }

memb14=recode(memb14 ,'1'=rank(lab)[1],'2'=rank(lab)[2],'3'=rank(lab)[3])

##
memb15 =cutree(hc15,k=3)
lab={};for(i in 1:3){
clas=mean(Data[memb15==i])
lab=c(lab,clas)
     }

memb15=recode(memb15 ,'1'=rank(lab)[1],'2'=rank(lab)[2],'3'=rank(lab)[3])

##"complete"
#mad11=1-get_percentage_of_matches(lab0,memb11)
mad11=mean(abs(lab0-memb11))
 
mmad11=c(mmad11,mad11)  

adj11=adjustedRandIndex(lab0,memb11)
adjust11=c(adjust11,adj11)
### "single"
#mad12=1-get_percentage_of_matches(lab0,memb12)
mad12= mean(abs(lab0-memb12))

mmad12=c(mmad12,mad12)  

adj12=adjustedRandIndex(lab0,memb12)
adjust12=c(adjust12,adj12)
## "average"
#mad13=1-get_percentage_of_matches(lab0,memb13)
 mad13= mean(abs(lab0-memb13))

mmad13=c(mmad13,mad13)  

adj13=adjustedRandIndex(lab0,memb13)
adjust13=c(adjust13,adj13)

## "ward.D"
#mad14=1-get_percentage_of_matches(lab0,memb14)
 mad14= mean(abs(lab0-memb14))

mmad14=c(mmad14,mad14)  

adj14=adjustedRandIndex(lab0,memb14)
adjust14=c(adjust14,adj14)
## "ward.D2"
#mad15=1-get_percentage_of_matches(lab0,memb15)
 mad15= mean(abs(lab0-memb15))

mmad15=c(mmad15,mad15)  

adj15=adjustedRandIndex(lab0,memb15)
adjust15=c(adjust15,adj15)

#############################
################# compute log-likelihood for a pair sample ############
a0=function(c1,c2){
N0 = length(c(c1,c2))
mui <- c(rep(0,N0 ))
dh[c(c1,c2),c(c1,c2)]
Sigmai=Sigma[c(c1,c2),c(c1,c2)]
Datai = Data[c(c1,c2)]
log.lik0 =dmvt(Datai,delta = mui,df=ddf, sigma = Sigmai,log=FALSE)# df=NULL means df=1
return(log.lik0)
}
#################################
ABadj=function(c1,c2){
N1 = length(c1);N2 = length(c2)
mu1 <- c(rep(0,N1 ));mu2 <- c(rep(0,N2 ))
dh1=dh[c1,c1];dh2=dh[c2,c2]
Sigma1=Sigma[c1,c1];Sigma2=Sigma[c2,c2]
Sigmai=Sigma[c(c1,c2),c(c1,c2)]
Datai = Data[c(c1,c2)]
mui=c(rep(0,length(Datai ) ))
Data1 = Data[c1];Data2 = Data[c2]
log.lik0 =dmvt(Datai,delta = mui,df=ddf, sigma = Sigmai,log=FALSE)

if(length(c1)==1&length(c2)==1){
dt(Data1, df=ddf,log = FALSE)
log.lik1=dt(Data1, df=ddf,log = FALSE)
log.lik2=dt(Data2, df=ddf,log = FALSE)
return(log.lik0/(N1+N2)/(log.lik1*log.lik2/(N1*N2)))
}

if(length(c1)>1&length(c2)==1){
log.lik1 =dmvt(Data1,delta = mu1,df=ddf, sigma = Sigma1,log=FALSE)
log.lik2=dt(Data2, df=ddf,log = FALSE)
return(log.lik0/(N1+N2)/(log.lik1*log.lik2/(N1*N2)))
}

if(length(c1)==1&length(c2)>1){
log.lik1=dt(Data1, df=ddf,log = FALSE)
log.lik2=dmvt(Data2,delta = mu2,df=ddf, sigma = Sigma2,log=FALSE)
return(log.lik0/(N1+N2)/(log.lik1*log.lik2/(N1*N2)))
}
if(length(c1)>1&length(c2)>1){
log.lik1=dmvt(Data[c1],delta =c(rep(0,length(c1))),df=ddf, sigma = Sigma[c1,c1],log=FALSE)
log.lik2=dmvt(Data2,delta = mu2,df=ddf, sigma = Sigma2,log=FALSE)
return(log.lik0/(N1+N2)/(log.lik1*log.lik2/(N1*N2)))
}
}

##################
# part 5        ##
############# clustring sample based on  loglikelihood ## ??? ???? ????? ???? 7 ? 8  ?? ??? ?????? ???? ???? ?? ?????? ??????
X=cbind(x,y,Data)
clustlog=function(X,a,nc=1){ 
n=dim(X)[1];cli={};for(i in 1:n)cli[i]=list(i)
stage=1
repeat{
fit={};ij={};t=1
for(i in 1:(length(cli)-1)){for(j in (i+1):length(cli)){
fit=c(fit,a(as.vector(cli[[i]]),as.vector(cli[[j]])))
ij[t]=list(c(i,j));t=t+1
}}
maxij=ij[which(fit==max(fit))]
cli[[maxij[[1]][1]]]=c(cli[[maxij[[1]][1]]],cli[[maxij[[1]][2]]])
cli=cli[-maxij[[1]][2]]
print(cli)
print(c("stage",stage));stage=stage+1 
if(length(cli)==nc){
lab={};cliv={}
for(tt in 1:length(cli)){lab[tt]=length(cli[[tt]]);cliv=c(cliv,cli[[tt]])}
#print(lab)
#print(cliv)
label=rep(1:length(cli),lab)
return(label[order(cliv)]);break}
}
}

memb2=clustlog(X[1:n,],a=a0,nc=3)  ## 2 number of cluster  ### ????? ???? ?????? ??
##################
lab={};for(i in 1:3){
clas=mean(Data[memb2==i])
lab=c(lab,clas)
     }

memb2=recode(memb2 ,'1'=rank(lab)[1],'2'=rank(lab)[2],'3'=rank(lab)[3])

mad2=1-get_percentage_of_matches(lab0,memb2) ## ???? ???? ???? ??? ??? ??- ????
mmad2=c(mmad2,mad2)
adj2=adjustedRandIndex(lab0,memb2)
adjust2=c(adjust2,adj2)
###################################################
memb4=clustlog(X[1:49,],a=ABadj,nc=3)  ## 3 number of cluster  ### ????? ???? ?????? ??
####
lab={};for(i in 1:3){
clas=mean(Data[memb4==i])
lab=c(lab,clas)
     }

memb4=recode(memb4 ,'1'=rank(lab)[1],'2'=rank(lab)[2],'3'=rank(lab)[3])

#mad4=mean(abs(lab0-memb4)) ## ???? ???? ???? ??? ??? ??- ????
mad4=1-get_percentage_of_matches(lab0,memb4)

mmad4=c(mmad4,mad4)
adj4=adjustedRandIndex(lab0,memb4)
adjust4=c(adjust4,adj4)
##########

}
p=0.00
matrix(c("Mean Eucl",mean(mmad11,trim=p),  mean(mmad12,trim=p),  mean(mmad13,trim=p),  mean(mmad14,trim=p),  mean(mmad15,trim=p),mean(mmad2,trim=p),
mean(mmad4,trim=p),
         "ARI  Eucl    ", mean(adjust11,trim=p),mean(adjust12,trim=p),mean(adjust13,trim=p),mean(adjust14,trim=p),mean(adjust15,trim=p),mean(adjust2,trim=p),
mean(adjust4,trim=p),
                  "Method","complete","single","average","ward.D","ward.D2","Log-Like1","Log-ABadj"),ncol=3)
}

#a2=ClusSpaTDLik(M=100,ddf=2,alfa=0.1,pow=1)

# a2
    # [,1]                   [,2]               
#[1,] "Mean Eucl"            "ARI  Eucl    "    
#[2,] "0.00102040816326531"  "0.997269980774512"
#[3,] "0.00408163265306122"  "0.994523809523809"
#[4,] "0.00408163265306122"  "0.994523809523809"
#[5,] "0.00102040816326531"  "0.997269980774512"
#[6,] "0.00142857142857143"  "0.996140948516448"
#[7,] "0.00285714285714286"  "0.992788622314912"
#[8,] "0.000612244897959183" "0.998133851707859" 
#     [,3]       
3[1,] "Method"   
#[2,] "complete" 
#[3,] "single"   
#[4,] "average"  
#[5,] "ward.D"   
#[6,] "ward.D2"  
#[7,] "Log-Like1"
#[8,] "Log-ABadj"


ptm <- proc.time()

a2=ClusSpaTDLik(M=5,ddf=1,alfa=0.1,pow=1)

 Stop the clock
proc.time() - ptm




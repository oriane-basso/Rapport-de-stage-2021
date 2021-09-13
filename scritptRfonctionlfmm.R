datedujour="2308"
#Input:
CHR="01"
datestructureresult="1708"
envcsv="lattitude1708VA"
namecsvenvoutput= "LatVA"
#variable du test:
Kl=6
iterationsnb=3
burninnb=1
repnb=4


library(LEA)
project=import.snmfProject(paste("CHR",CHR, "genoformat",datestructureresult,"_snmfProject.zip",sep=""),force=TRUE)


vcf2lfmm(paste("Chr_",CHR,".vcf",sep=""))
best = which.min(cross.entropy(project, K = Kl))
# Impute the missing genotypes

impute(project, paste("Chr_",CHR,".lfmm",sep=""),
       method = 'mode', K=Kl, run = best)

Y<- read.csv2(paste(envcsv, ".csv",sep=""))
write.env(Y,paste(namecsvenvoutput,"env",sep="."))

project.lfmm=lfmm(paste("Chr_",CHR, ".lfmm_imputed.lfmm",sep=""),paste(namecsvenvoutput,".env",sep=""),K=Kl,
                  iterations=iterationsnb,burnin=burninnb,repetitions=repnb)

export.lfmmProject(paste("Chr_",CHR,".lfmmProject",sep=""))

# compute adjusted p-values (justement j'ai pas l'impression qu'elles soient ajuster mes pvalues ici? Si elles sont ajustés, en fait y'a déjà ajustement pour confouding factor. Par contre pas ajusté avec inflation factor )
#lfmm.pvalues utilise Fisher/stouffer method pour combiner les z-scores pour des runs multiples (à creuser cf: Olivier F)
p = lfmm.pvalues(project.lfmm, Kl) #lfmm estime facteurs confondant et effect size basé sur algo MCMC et lfmm.pvalues permet identifier les SNP qui montre association avec eco gradient, en corrigant par facteurs confondant
#lfmm.pvalues returns a vector of p values computed from a combination of lfmm runs 
#NOTE: il semble que je peux directement corriger par l'inflation facteur grâce à cetet fonction si je lu ajoute lambda?
pvalues = p$pvalues
#GWAS significance test
pdf(file = "plotpv.pdf")
par(mfrow = c(2,1))
hist(pvalues, col = "lightblue")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .7)
dev.off()



#multiple testing issues
EXPFDR<-c() # expected FDR
OBSFDR<-c()
ESTTPR<-c()
i=1
c05<-c()
c1<-c()
c015<-c()
c2<-c()

# Benjamini-Hochberg's algorithm:
# estimated FDR and True Positive Rate
for (alpha in c(.05,.1,.15,.2)) {
  EXPFDR[i]<-alpha
  L = length(pvalues)
  w = which(sort(pvalues) < alpha * (1:L) / L)
  candidates = order(pvalues)[w]
  
  if (i==1) {c05=candidates}
  if (i==2) {c1=candidates}
  if (i==3) {c15=candidates}
  if (i==4) {c2=candidates}
  
  Lc = length(candidates)
  estimated.FDR = sum(candidates <= 350)/Lc
  OBSFDR[i]=round(estimated.FDR, digits = 2)
  
  estimated.TPR = sum(candidates > 350)/50
  ESTTPR[i]=round(estimated.TPR, digits = 2)
  i=i+1}


multipletestingissues<-cbind(EXPFDR,OBSFDR,ESTTPR)
max_length<-max(c(length(c05),length(c1),length(c015),length(c2)))


#Pour 0,05
pvalues2 = data.frame(pvaluesbefore=p$pvalues)
pvalues2$sort<-1:nrow(pvalues2)
cand05<-data.frame(sort=c05)
cd05<-merge(pvalues2,cand05,by="sort")
colnames(cd05)<-c("c05","pvalue")

#Pour 0,1
cand1<-data.frame(sort=c1)
cd1<-merge(pvalues2,cand1,by="sort")
colnames(cd1)<-c("c1","pvalue")

#Pour 0,15
cand15<-data.frame(sort=c15)
cd15<-merge(pvalues2,cand15,by="sort")
colnames(cd15)<-c("c15","pvalue")

#Pour 0,2
cand2<-data.frame(sort=c2)
cd2<-merge(pvalues2,cand2,by="sort")
colnames(cd2)<-c("c2","pvalue")

candidatetot<-data.frame(candidate05=c(cd05$c05,rep(NA,max_length-length(cd05$c05))),
                         pvaluebeforealgo05=c(cd05$pvalue,rep(NA,max_length-length(cd05$pvalue))),
                         candidate1=c(cd1$c1,rep(NA,max_length-length(cd1$c1))),
                         pvaluebeforealgo1=c(cd1$pvalue,rep(NA,max_length-length(cd1$pvalue))),
                         candidate15=c(cd15$c15,rep(NA,max_length-length(cd15$c15))),
                         pvaluebeforealgo15=c(cd15$pvalue,rep(NA,max_length-length(cd15$pvalue))),
                         candidate2=c(cd2$c2,rep(NA,max_length-length(cd2$c2))),
                         pvaluebeforealgo2=c(cd2$pvalue,rep(NA,max_length-length(cd2$pvalue))))

write.csv(candidatetot,paste("candidatetot.csv",datedujour,sep=""))
write.csv(multipletestingissues,paste("multipletestingissues",datedujour,".csv",sep=""))
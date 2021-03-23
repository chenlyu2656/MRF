library(betareg)
library(CompQuadForm)

source('MRF_fun.r')

# trait sample (n=300)
# y0, y1 and y2 represent methylation traits (beta values), ranging from [0, 1] with each following a beta distribution
# y0 was simulated with no causal association with genotypes
# y1 was simulated to causally associated with 10% of SNPs in one direction (positive associations)
# y2 was simulated to causally associated with 10% of SNPs in two directions (50% positive associations and 50% negative associations)
trait <-read.table("Trait.txt", header=T)

# genotype sample (n=300, nsnp=358)
# most SNPs are rare variants with MAF < 0.05
geno <- read.table("Genotype.txt", header=T)

# covariate sample
# x1 is a quantitative coviate, while x2 is a binary covariate
cov <- read.table("Covariate.txt", header=T)

# MRF
MAF <- colMeans(geno)/2
wt<-dbeta(MAF,1,25)

# Assuming methylation trait follows a beta distribution using GR kernel (other kernel functions can be applied)
p.ggrf.b0   <-MRF(trait$y0,Z=as.matrix(geno),X=as.matrix(cov),weights=(wt^2),type='beta',kernel='GR')$pvalue
p.ggrf.b1   <-MRF(trait$y1,Z=as.matrix(geno),X=as.matrix(cov),weights=(wt^2),type='beta',kernel='GR')$pvalue
p.ggrf.b2   <-MRF(trait$y2,Z=as.matrix(geno),X=as.matrix(cov),weights=(wt^2),type='beta',kernel='GR')$pvalue

# MRF can also be used assuming methylation trait follows a normal distribution after logit transformation
y.l0<-log(trait$y0/(1-trait$y0))
y.l1<-log(trait$y1/(1-trait$y1))
y.l2<-log(trait$y2/(1-trait$y2))

p.ggrf.l0   <-MRF(y.l0,Z=as.matrix(geno),X=as.matrix(cov),weights=(wt^2),type='normal',kernel='GR')$pvalue
p.ggrf.l1   <-MRF(y.l1,Z=as.matrix(geno),X=as.matrix(cov),weights=(wt^2),type='normal',kernel='GR')$pvalue
p.ggrf.l2   <-MRF(y.l2,Z=as.matrix(geno),X=as.matrix(cov),weights=(wt^2),type='normal',kernel='GR')$pvalue

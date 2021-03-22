source('MRF_fun.r')

# trait sample (n=300)
# y0, y1 and y2 ranges from [0, 1] with each following a beta distribution
# y.l0, y.l1 and y.l0 are the corresponding logit transformations
trait <-read.table("Trait.txt", header=T)

# genotype sample (n=300, nsnp=358)
# most SNPs are rare variants with MAF < 0.05
geno <- read.table("Genotype.txt", header=T)

# covariate sample
# x1 is a quantitative cov, while x2 is a binary cov
cov <- read.table("Covariate.txt", header=T)

# GGRF
MAF <- colMeans(geno)/2
wt<-dbeta(MAF,1,25)

# Assuming trait follows a beta distribution
p.ggrf.b0   <-MRF(trait$y0,Z=as.matrix(geno),X=as.matrix(cov),weights=(wt^2),type='beta',kernel='GR')$pvalue
p.ggrf.b1   <-MRF(trait$y1,Z=as.matrix(geno),X=as.matrix(cov),weights=(wt^2),type='beta',kernel='GR')$pvalue
p.ggrf.b2   <-GGRF(trait$y2,Z=as.matrix(geno),X=as.matrix(cov),weights=(wt^2),type='beta',kernel='GR')$pvalue

# Assuming trait follows a normal distribution after logit transformation
p.ggrf.l0   <-MRF(trait$y.l0,Z=as.matrix(geno),X=as.matrix(cov),weights=(wt^2),type='normal',kernel='GR')$pvalue
p.ggrf.l1   <-MRF(trait$y.l1,Z=as.matrix(geno),X=as.matrix(cov),weights=(wt^2),type='normal',kernel='GR')$pvalue
p.ggrf.l2   <-MRF(trait$y.l2,Z=as.matrix(geno),X=as.matrix(cov),weights=(wt^2),type='normal',kernel='GR')$pvalue

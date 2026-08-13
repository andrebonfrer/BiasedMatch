# Simulations for Examples (a) and (b)

# Installing packages 
# Only required for the first time of use
#install.packages("BiasedUrn"); install.packages("nortest"); 
#install.packages("EnvStats"); install.packages("coda")

setwd("/Users/andrebonfrer/Dropbox/Research/MatchModels/Analysis/BiasInference/MHAlg")
# Loading packages
# Required for every time of use
library("BiasedUrn"); 
library(nortest); 
library(EnvStats); 
library(coda)
library(xtable) 

# Defining functions 
source("MHALGfun.R")

# flags
run.MLE <- TRUE
run.plot <- TRUE

# Calculation for Example (a) in the paper
set.seed(20); 
N1=150; N2=200; N=50; m1=40; m2=30; mt<-10
Bn <- 250; Jn = 10000

system.time(resFG <- MHALG(B=Bn, J=Jn, N1vec=N1, N2vec=N2, N0vec=N, 
		      M1vec=m1, M2vec=m2, mvec=mt,  		
		Cvec=1, eta=3230, tau=10000, mu0=1.15, sig0=10^-4, 
		muinit=1.15, rinit=0, del=0.7, hf=1)
)
resFG$rar  # 0.81  OK

# check convergence
s <- 1; step <- 10
tg1<- resFG$muv[-(1:Bn)][seq(s, Jn, step)];   
tg2<- resFG$lamv[-(1:Bn)][seq(s, Jn, step)];   
tg3<- resFG$rm[-(1:Bn)][seq(s, Jn, step)];   

s1<-c(serialCorrelationTest(tg1,test="rank.von.Neumann")$p.value, serialCorrelationTest(tg2,test="rank.von.Neumann")$p.value,
      serialCorrelationTest(tg3,test="rank.von.Neumann")$p.value)
s2<-c(serialCorrelationTest(tg1,test="AR1.yw")$p.value, 
      serialCorrelationTest(tg2,test="AR1.yw")$p.value,
      serialCorrelationTest(tg3,test="AR1.yw")$p.value)
s3<-c(ad.test(tg1)$p.value, ad.test(tg2)$p.value, ad.test(tg3)$p.value)
s4<-c(cvm.test(tg1)$p.value, cvm.test(tg2)$p.value, cvm.test(tg3)$p.value)
s5<-c(lillie.test(tg1)$p.value, lillie.test(tg2)$p.value, lillie.test(tg3)$p.value)
s6<-c(min((1-pnorm(geweke.diag(tg1)$z))*2, (1-pnorm(-geweke.diag(tg1)$z))*2),
      min((1-pnorm(geweke.diag(tg2)$z))*2, (1-pnorm(-geweke.diag(tg2)$z))*2),
      min((1-pnorm(geweke.diag(tg3)$z))*2, (1-pnorm(-geweke.diag(tg3)$z))*2))


rn <- NULL; rn <- c(rn, "SerialCorrelationTest[rank.von.Neumann]",
                    "SerialCorrelationTest[AR1.yw]", "NormalityTest[Anderson-Darling]",
                    "NormalityTest[Cramer-von Mises]", "NormalityTest[Lilliefors]",
                    "ConvergenceTest[Geweke]")
estparmat<- rbind(s1,s2,s3,s4,s5,s6)
estparmat <- data.frame(estparmat,row.names = rn)
colnames(estparmat) <- c("mu","lambda","r")

# table for convergence
xtable(estparmat,digits=c(0, rep(3,3)))



if(run.plot) {
  
wvec=seq(0.1,20,0.01); 
lenwvec=length(wvec); 
L0vec=rep(NA,lenwvec) 
for(i in 1:lenwvec){ 
  wval=wvec[i]
  L0vec[i]=PROBM(m=mt, N1=N1,N2=N2,N=N,m1=m1,m2=m2,w=wval) 
  }
K0vec=L0vec*dlnorm(wvec,1.15,1.76); 
const=INTEG(wvec,K0vec); 
expt=INTEG(wvec, wvec*K0vec)/const
vart=INTEG(wvec, ((wvec-expt)^2)*K0vec)/const
wvec[which.max(K0vec)]         # posterior mode
expt                               # posterior mean
sqrt(vart)                     # posterior standard deviation
  
# Plotting Figure for trace plots
pdf("Figures/traceplots.pdf",width=10)
par(mfrow=c(2,2)); 
plot(resFG$muv,type="l",xlab=" ",ylab=""); 
mtext("m",family="HersheySymbol",side=2,padj=-2.5,font=3,cex=1.5)
mtext("iteration",family="serif",side=1,padj=3.5,font=1,cex=1.5)
plot(resFG$lamv,type="l",xlab="",ylab="");
mtext("l",family="HersheySymbol",side=2,padj=-2.5,font=3,cex=1.5)
mtext("iteration",family="serif",side=1,padj=3.5,font=1,cex=1.5)
plot(resFG$rm,type="l",xlab="",ylab=""); 
mtext("r",family="serif",side=2,padj=-2.5,font=3,cex=1.5)
mtext("iteration",family="serif",side=1,padj=3.5,font=1,cex=1.5)
plot(exp(resFG$rm),type="l",xlab="",ylab="")
mtext("w",family="serif",side=2,padj=-2.5,font=3,cex=1.5)
mtext("iteration",family="serif",side=1,padj=3.5,font=1,cex=1.5)
dev.off()

# Plotting Figure for density in main paper for Example (a)
pdf("Figures/histExample1.pdf",width=9)
hist(exp(resFG$rm),prob=T,nclass=50,xlab="",ylab="",main="",xlim=c(0,22),ylim=c(0,0.21)); 
mtext("Density",family="serif",side=2,padj=-3.5,font=3,cex=1.5)
mtext("w",family="serif",side=1,padj=3.5,font=3,cex=1.5)

den=density(exp(resFG$rm),adjust=1.5)
imode=(1:length(den$x))[den$y==max(den$y)]; 
wmode=den$x[imode]
c(imode,wmode)

lines(den,lty=2,lwd=2); 
abline(v=wmode,lty=2,lwd=2)
lines(wvec,K0vec/const,lwd=3,lty=1)
abline(v=wvec[which.max(K0vec)],lty=1,lwd=2) 
legend(13,0.15,c("Exact posterior & mode", "MC estimate & mode"),lty=c(1,2),lwd=c(4,2))
dev.off()
}


# Calculating the posterior mode and the MLE for Example (a) (reported in a table in the paper)
if(run.MLE){ 
    mvec<-seq(0,30); nm<-length(mvec); mdv<-mlev<-rep(NA,nm)
    for(j in 1:nm){
      mt<-mvec[j]
      wvec<-c(0, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2, seq(0.1, 200, 0.1), 
      seq(200, 1000, 100), 10^4, 10^5, 10^6, 10^7);   
      lenwvec=length(wvec); 
      L0vec=rep(NA,lenwvec) 
        for(i in 1:lenwvec){ 
          wval=wvec[i]
	        L0vec[i]=PROBM(m=mt, N1=N1,N2=N2,N=N,m1=m1,m2=m2,w=wval) 
	      }
    K0vec=L0vec*dlnorm(wvec,1.15,1.76); 
    mlev[j]<-wvec[which.max(L0vec)]
    mdv[j]<-wvec[which.max(K0vec)]
    }
}



###################################################################################################
# Calculation for Example (b) - set to Nk = 5 for first example, Nk = 50 for the second
N1=150; N2=200; N=50; m1=40; m2=30;
Nk = 50; K=6; n = K*Nk; Cvec=NULL; for(k in 1:K) Cvec=c(Cvec,rep(k,n/K)); Cvec

N1vec=rep(N1,n); N2vec=rep(N2,n); N0vec=rep(N,n)
M1vec=rep(m1,n); M2vec=rep(m2,n); mvec=rep(NA,n)

# suggested prior to stabilize mu
eta=3230; tau=10000; mu0=1.15; sig0=0.0001
set.seed(20)
mu<-rnorm(1,mu0,sig0); lam<-rgamma(1,eta,tau); sig<-sqrt(1/lam)

set.seed(10); 
rvec=rnorm(K,mu,sig); rvec 
wvec=exp(rvec); wvec 

set.seed(30); for(i in 1:n){	Wval=wvec[Cvec[i]]
	N1val=N1vec[i]; N2val=N2vec[i]; N0val=N0vec[i]
	M1val=M1vec[i]; M2val=M2vec[i]
	  y1samp=rWNCHypergeo(nran=1, m1=N0val, m2=N1val-N0val,
		n=M1val, odds=Wval)
	mqvecsamp=rMWNCHypergeo(nran=1, 
	    m=c(y1samp,N0val-y1samp,N2val-N0val), M2val, c(Wval,Wval,1))
	mvec[i]=mqvecsamp[1]   } 
mvec

Jn<-1000;  Bn<-250
set.seed(50);    # Change the number for different simulations

eta <- 3230; tau<- 10000;  mu0 = 1.15; sig0 = 0.0001; 
system.time(resg<-MHALG(B=Bn, J=Jn,  N1vec=N1vec,  N2vec=N2vec, N0vec=N0vec, 
	M1vec=M1vec,  M2vec=M2vec, mvec=mvec,   Cvec=Cvec, 
	eta=eta, tau=10^4,  mu0 = mu0, sig0=sig0,  muinit=1.15, 
	rinit=rvec, del=0.75*c(0.4, 0.4, 0.75, 1.2, 0.4, 0.3)  )
)
resg$rar    #   acceptance rates

s <- 4; step = 10
tg1<- resg$muv[-(1:Bn)][seq(s, Jn, step)];   
tg2<- resg$lamv[-(1:Bn)][seq(s, Jn, step)];   
tg3<- resg$rm[-(1:Bn),][seq(s, Jn, step),];   

s1<-c(serialCorrelationTest(tg1,test="rank.von.Neumann")$p.value, serialCorrelationTest(tg2,test="rank.von.Neumann")$p.value,
      apply(tg3, 2, function(x) serialCorrelationTest(x,test="rank.von.Neumann")$p.value))
s2<-c(serialCorrelationTest(tg1,test="AR1.yw")$p.value, 
      serialCorrelationTest(tg2,test="AR1.yw")$p.value,
      apply(tg3, 2, function(x) serialCorrelationTest(x,test="AR1.yw")$p.value))
s3<-c(ad.test(tg1)$p.value, ad.test(tg2)$p.value, apply(tg3, 2, function(x) ad.test(x)$p.value))
s4<-c(cvm.test(tg1)$p.value, cvm.test(tg2)$p.value, apply(tg3, 2, function(x) cvm.test(x)$p.value))
s5<-c(lillie.test(tg1)$p.value, lillie.test(tg2)$p.value, apply(tg3, 2, function(x) lillie.test(x)$p.value))
s6<-c(min((1-pnorm(geweke.diag(tg1)$z))*2, (1-pnorm(-geweke.diag(tg1)$z))*2),
      min((1-pnorm(geweke.diag(tg2)$z))*2, (1-pnorm(-geweke.diag(tg2)$z))*2),
      apply(tg3, 2, function(x) min((1-pnorm(geweke.diag(x)$z))*2, (1-pnorm(-geweke.diag(x)$z))*2))
      )


rn <- NULL; rn <- c(rn, "SerialCorrelationTest[rank.von.Neumann]",
                    "SerialCorrelationTest[AR1.yw]", "NormalityTest[Anderson-Darling]",
                    "NormalityTest[Cramer-von Mises]", "NormalityTest[Lilliefors]",
                    "ConvergenceTest[Geweke]")
estconv<- rbind(s1,s2,s3,s4,s5,s6)
estconv <- data.frame(estconv,row.names = rn)
colnames(estconv) <- c("mu","lambda",paste0("r",1:K))

# table for convergence
xtable(estconv,digits=c(0, rep(3,2+K)))

estparmat <- restable(resg,Bn=250,sim=TRUE)

xtable(estparmat,digits = c(0, rep(3,3+K)))




# Use SMRMMHALG

Bn <- 500
Jn <- 1000
X<-model.matrix(~factor(Cvec))
resg <- SMRMMHALG(B=Bn, J=Jn,  N1vec=N1vec,  N2vec=N2vec, N0vec=N0vec, 
	M1vec=M1vec,  M2vec=M2vec, mvec=mvec, X = X, 
	B0 = matrix(c(1,rep(0,dim(X)[2]-1)),nrow=1),
	V0 = 1*diag(dim(X)[2]),
	eta=eta, tau=10^4,  mu0=mu0, sig0=10^-4,  muinit=1.15, 
	rinit=matrix(rep(1, nrow(X)), nrow=1), del = c(0.4, 0.4, 0.35, 0.7, 0.4, 0.3))



True  value
3.2067
0.1137
4.3082
-1.2070
0.3602
2.1576
Posterior Mean
2.9450
0.0974
4.3590
-2.3064
0.6531
1.9242
Posterior Mode
2.9201
0.1269
4.2021
-2.1257
0.6355
1.9052
95% Credible Interval
(2.5264, 3.3992)
(-0.3296, 0.5307)
(3.6575, 5.2244)
(-4.3004, -0.9793)
(0.2904, 1.0403)
(1.5923, 2.2733)
Batch means
2.9490
0.1008
4.3636
-2.3092
0.6529
1.9239
All the true values are contained in the 95% credible intervals.


# Calculating the posterior means via the batch means approach
install.packages("batchmeans")
library(batchmeans)
bm(resg$rm[Bn:Jn,1]); bm(resg$rm[Bn:Jn,2]); bm(resg$rm[Bn:Jn,3])
bm(resg$rm[Bn:Jn,4]); bm(resg$rm[Bn:Jn,5]); bm(resg$rm[Bn:Jn,6]) 


# Calculating mean of means via resampling
RS<-100                   # Numbers of times for resampling
SD<-sample.int(1000, RS)     # Randomly generating values of seeds
muvr<-lamvr<-sigvr<-r1mnvr<-r2mnvr<-r3mnvr<-r4mnvr<-r5mnvr<-r6mnvr<-rep(NA,RS)
r1mdvr<-r2mdvr<-r3mdvr<-r4mdvr<-r5mdvr<-r6mdvr<-rep(NA,RS)
r1lbvr<-r2lbvr<-r3lbvr<-r4lbvr<-r5lbvr<-r6lbvr<-rep(NA,RS)
r1ubvr<-r2ubvr<-r3ubvr<-r4ubvr<-r5ubvr<-r6ubvr<-rep(NA,RS)

# CAUSION: EXTREMELY LONG TIME CALCULATION ( > 60 hours )
Jn<-1000;  Bn<-50
for(j in 1:RS){
set.seed(SD[j]);    # Change the number for different simulations
resg<-MHALG(B=Bn, J=Jn,  N1vec=N1vec,  N2vec=N2vec, N0vec=N0vec, 
	M1vec=M1vec,  M2vec=M2vec, mvec=mvec,   Cvec=Cvec, 
	eta=eta, tau=10^4,  mu0=mu0, sig0=10^-4,  muinit=1.15, 
	rinit=rvec, del=c(0.4, 0.4, 0.75, 1.2, 0.4, 0.3)  )
muvr[j]<- mean(resg$muv[Bn:Jn])
lamvr[j]<- mean(resg$lamv[Bn:Jn])
sigvr[j]<- mean(1/sqrt(resg$lamv[Bn:Jn]))
r1mnvr[j]<- mean(resg$rm[Bn:Jn,1]); r2mnvr[j]<- mean(resg$rm[Bn:Jn,2])
r3mnvr[j]<- mean(resg$rm[Bn:Jn,3]); r4mnvr[j]<- mean(resg$rm[Bn:Jn,4])
r5mnvr[j]<- mean(resg$rm[Bn:Jn,5]); r6mnvr[j]<- mean(resg$rm[Bn:Jn,6])
den=density( (resg$rm[Bn:Jn, 1]),adjust=1.5)
imode=(1:length(den$x))[den$y==max(den$y)]; r1mdvr[j]=den$x[imode]
r1lbvr=quantile(resg$rm[Bn:Jn, 1],0.025); r1ubvr=quantile(resg$rm[Bn:Jn, 1],0.975)
den=density( (resg$rm[Bn:Jn, 2]),adjust=1.5)
imode=(1:length(den$x))[den$y==max(den$y)]; r2mdvr[j]=den$x[imode]
r2lbvr=quantile(resg$rm[Bn:Jn, 2],0.025); r2ubvr=quantile(resg$rm[Bn:Jn, 2],0.975)
den=density( (resg$rm[Bn:Jn, 3]),adjust=1.5)
imode=(1:length(den$x))[den$y==max(den$y)]; r3mdvr[j]=den$x[imode]
r3lbvr=quantile(resg$rm[Bn:Jn, 3],0.025); r3ubvr=quantile(resg$rm[Bn:Jn, 3],0.975)
den=density( (resg$rm[Bn:Jn, 4]),adjust=1.5)
imode=(1:length(den$x))[den$y==max(den$y)]; r4mdvr[j]=den$x[imode]
r4lbvr=quantile(resg$rm[Bn:Jn, 4],0.025); r4ubvr=quantile(resg$rm[Bn:Jn, 4],0.975)
den=density( (resg$rm[Bn:Jn, 5]),adjust=1.5)
imode=(1:length(den$x))[den$y==max(den$y)]; r5mdvr[j]=den$x[imode]
r5lbvr=quantile(resg$rm[Bn:Jn, 5],0.025); r5ubvr=quantile(resg$rm[Bn:Jn, 5],0.975)
den=density( (resg$rm[Bn:Jn, 6]),adjust=1.5)
imode=(1:length(den$x))[den$y==max(den$y)]; r6mdvr[j]=den$x[imode]
r6lbvr=quantile(resg$rm[Bn:Jn, 6],0.025); r6ubvr=quantile(resg$rm[Bn:Jn, 6],0.975)
}

mean(muvr); mean(lamvr); mean(sigvr); mean(r1mnvr); mean(r2mnvr); mean(r3mnvr)
mean(r4mnvr); mean(r5mnvr); mean(r6mnvr)

# Plotting Figure 3
x11(9,5.5)
par(mfrow=c(2,2));
plot(resg$muv,type="l",xlab=" ",ylab=""); 
mtext("m",family="HersheySymbol",side=2,padj=-3.5,font=3)
mtext("iteration",family="serif",side=1,padj=3.5,font=1)
plot(resg$lamv,type="l",xlab="",ylab="");
mtext("l",family="HersheySymbol",side=2,padj=-3.5,font=3)
mtext("iteration",family="serif",side=1,padj=3.5,font=1)
plot(resg$rm[,1],type="l",xlab="",ylab=""); 
mtext("r",family="serif",side=2,padj=-3.5,font=3)
mtext("1",side=2,padj=-5.5,adj=0.55,cex=0.48)
mtext("iteration",family="serif",side=1,padj=3.5,font=1)
plot(resg$rm[,2],type="l",xlab="",ylab=""); 
mtext("r",family="serif",side=2,padj=-3.5,font=3)
mtext("2",side=2,padj=-5.5,adj=0.55,cex=0.48)
mtext("iteration",family="serif",side=1,padj=3.5,font=1)

par(mfrow=c(2,2));
plot(resg$rm[,3],type="l",xlab="",ylab=""); 
mtext("r",family="serif",side=2,padj=-3.5,font=3)
mtext("3",side=2,padj=-5.5,adj=0.55,cex=0.48)
mtext("iteration",family="serif",side=1,padj=3.5,font=1)
plot(resg$rm[,4],type="l",xlab="",ylab=""); 
mtext("r",family="serif",side=2,padj=-3.5,font=3)
mtext("4",side=2,padj=-5.5,adj=0.55,cex=0.48)
mtext("iteration",family="serif",side=1,padj=3.5,font=1)
plot(resg$rm[,5],type="l",xlab="",ylab=""); 
mtext("r",family="serif",side=2,padj=-3.5,font=3)
mtext("5",side=2,padj=-5.5,adj=0.55,cex=0.48)
mtext("iteration",family="serif",side=1,padj=3.5,font=1)
plot(resg$rm[,6],type="l",xlab="",ylab=""); 
mtext("r",family="serif",side=2,padj=-3.5,font=3)
mtext("6",side=2,padj=-5.5,adj=0.55,cex=0.48)
mtext("iteration",family="serif",side=1,padj=3.5,font=1)



# Juice 
N1 <- 168; N2 <- 125; N <- 104; m1 <- m2 <- 25
MA <- seq(1,N1,1)
wl = 2.41; wu = 3.85; wm = 3.11


# median m2
mmvec1 <- mlvec1 <- muvec1 <- rep(0,length(MA))
for(ma in MA) mlvec1[ma] <- DISTM(N1,N2,N,m1=ma, m2=m2, w= wl)$Em
for(ma in MA) mmvec1[ma] <- DISTM(N1,N2,N,m1=ma, m2=m2, w= wm)$Em
for(ma in MA) muvec1[ma] <- DISTM(N1,N2,N,m1=ma, m2=m2, w= wu)$Em

# min m2 (12)
m2 <- 12
mmvec2 <- mlvec2 <- muvec2 <- rep(0,length(MA))
for(ma in MA) mlvec2[ma] <- DISTM(N1,N2,N,m1=ma, m2=m2, w= wl)$Em
for(ma in MA) mmvec2[ma] <- DISTM(N1,N2,N,m1=ma, m2=m2, w= wm)$Em
for(ma in MA) muvec2[ma] <- DISTM(N1,N2,N,m1=ma, m2=m2, w= wu)$Em


# max m2 (42)
m2 <- 42
mmvec3 <- mlvec3 <- muvec3 <- rep(0,length(MA))
for(ma in MA) mlvec3[ma] <- DISTM(N1,N2,N,m1=ma, m2=m2, w= wl)$Em
for(ma in MA) mmvec3[ma] <- DISTM(N1,N2,N,m1=ma, m2=m2, w= wm)$Em
for(ma in MA) muvec3[ma] <- DISTM(N1,N2,N,m1=ma, m2=m2, w= wu)$Em

# no bias, mean m2
m2 <- 25
mvec0 <- rep(0,length(MA))
for(ma in MA) mvec0[ma] <- DISTM(N1,N2,N,m1=ma, m2=m2, w= 1)$Em

fontfamily = "Times"

pdf(file="Figures/EmEV.pdf", family=fontfamily)

par("font.lab" = 1, "font.axis" = 3, mgp = c(2.5,0.85,0))
plot(MA, muvec3, type='l', col='red', cex.lab = 1.25, ylab = expression("Expected matches ("*italic("E[m])")), xlab = expression("Number of items promoted ("*italic("M"[A])*")") )
#mtext(, side = 1, padj = 4, adj = 0.73)
lines(mmvec3,col='black',lwd=3)
lines(mlvec3,col='red')
text(85,29, expression(italic("M"[B]*"=42")))

lines(muvec1, type='l', col='red')
lines(mmvec1,col='black',lwd=3)
lines(mlvec1,col='red')
text(95,19, expression(italic("M"[B]*"=25")))
text(100,14, expression(italic("M"[B]*"=25")), col='darkgreen')

lines(muvec2, type='l', col='red')
lines(mmvec2,col='black',lwd=3)
lines(mlvec2,col='red')
text(140,12, expression(italic("M"[B]*"=12")))

lines(mvec0,col='darkgreen',lwd=2)

la <- c(expression(italic("w = 2.41,3.85")), expression(italic("w = 3.11")), expression(italic("w = 1"))) 
legend(10,38,la,lty=c(1,1,1), lwd = c(1,3,2), col=c("red","black","darkgreen"))

dev.off()



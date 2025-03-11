Nth = read.csv('Nth.csv')
price = read.csv('price.csv')
prom = read.csv('prom.csv')

## In data Nth and following analysis, h=consumer id, b = basket id, bn = basket count, t = time id
## price and prom are the category level variable for 50 weeks and 10 categories


require(MASS)
require(flexmix)
require(reshape2)
N = 5 #number of category

T = 50 #number of weeks

H = 100 #number of consumers


## prepare the explanatory variables for the poisson regression
## check Appendix A
##create all the possible baskets combinations
x = 1:N
B = 2^N-1
Z = matrix(0,B,N)
n = 0
for (i in 1:N){
  c =t(combn(x,i))
  J = nrow(c)
  for (j in 1:J){
    Z[n+j,c[j,]]=1
  }
  n = n+J
}

## Each row of Z corresponds to one basket content zb

anam <- paste("alpha", 1:N, sep="")

Zb = as.data.frame(Z)
colnames(Zb) = anam
Zb$b = 1:B



thetab = matrix(0,B,N*(N-1)/2)

for(i in 1:B){
  Zi = matrix(c(Z[i,],rep(0,N)),N,2)
  Zi = Zi%*%t(Zi)
  thetab[i,] = Zi[lower.tri(Zi)]
}

thetab = as.data.frame(thetab)
thetab$b = 1:B





## bundle price and prom at each week t
priceb = Z%*%t(price)
priceb = as.data.frame(priceb)
colnames(priceb)=1:T
priceb$b = 1:B
priceb = melt(priceb,id.vars='b')
colnames(priceb)=c('b','t','price')

promb = Z%*%t(prom)
colnames(promb)=1:T
promb = as.data.frame(promb)
promb$b = 1:B
promb = melt(promb,id.vars='b')
colnames(promb)=c('b','t','prom')








#### one-step approach
## short run regression data
snhb1 = merge(Nth,Zb, by='b')
## If using a subset of basket, e.g. basket of size 1 and 2, B1 = N*(N+1)/2+2
### snhb1 = merge(Nth,Zb[1:B1,], by='b') 
snhb1 = merge(snhb1,thetab,by='b')
snhb1 = merge(snhb1,priceb,by=c('b','t'))
snhb1 = merge(snhb1,promb,by=c('b','t'))






## latent-class poisson regression

Model <- FLXMRglmfix(family = "poisson")

f7s=stepFlexmix(bn~ (alpha1 + alpha2 + alpha3 + alpha4 + alpha5
                     +V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10+price+prom)|h, 
                k=4:6, data = snhb1,  model = Model,nrep=5)
f4s =getModel(f7s,'BIC')

## In the above regression coefficients of V1...V10 represent the lower_tri elements of Theta matrix: 
### theta_21, theta_31, theta_32, theta_41, theta_42, theta_43, theta_51,theta_52, theta_53, theta_54.


## results
pis = prior(f4s) ## segment sizes
S = length(pis) ## number of segments

## if S = 6, one should re-run the latent-class analysis with higher range for k



## consumer-segment membership, posterior probability that of consumer h in segment s
tao = posterior(f4s)
tau = as.data.frame(tao)
tau$h = snhb1$h
tau = unique(tau)
tao=as.matrix(tau[,1:S])





## preferences
alphas = t(parameters(f4s)[2:(N+1),])

## cross-category dependent matrix
Theta = list()
for(s in 1:4){
  Thetas = matrix(0,N,N) 
  Thetas[lower.tri(Thetas)]= c(parameters(f4s)[(N+2):(N*(N+1)/2+1),s])
  Thetas[upper.tri(Thetas)] = t(Thetas)[upper.tri(Thetas)]
  Theta[[s]]=Thetas
}
## Theta[[s]] is the cross-category dependent matrix of segment s

## marketing mix coefficients
beta1s = parameters(f4s)['coef.price',] ## price coefficient
beta2s = parameters(f4s)['coef.prom',] ## prom coefficient

## frequency scale
Rs = parameters(f4s)['coef.(Intercept)',]














### two-step approach

## step 1
## long-run data
lnhb = aggregate(Nth$bn, by = list(Nth$h,Nth$b), FUN = sum)
colnames(lnhb) = c('h','b','bn')


## lnhb is the long-run basket data, now we merge in the independent variables

lnhb = merge(lnhb,Zb, by='b')
lnhb = merge(lnhb,thetab,by='b')


#Latent class regression to Estimate segments' long-run Alpha and Theta

Model <- FLXMRglmfix(family = "poisson")


f7 <- stepFlexmix(bn~ (alpha1 + alpha2 + alpha3 + alpha4 + alpha5+V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10)|h, 
                  k=4:6, data = lnhb,  model = Model,nrep=5)

## In the above regression coefficients of V1...V10 represent the lower_tri elements of Theta matrix: 
### theta_21, theta_31, theta_32, theta_41, theta_42, theta_43, theta_51,theta_52, theta_53, theta_54.


f4 =getModel(f7,'BIC')

pi0 = prior(f4)

S = length(pi0) ## number of segments

## if S = 6, one should re-run the latent-class analysis with higher range for k



##step 2

## resume long-run data and add it in short-run model


## consumer-segment id: posterior distribution of consumer h in segment s
tao = posterior(f4)
tau = as.data.frame(tao)
tau$h = lnhb$h
tau = unique(tau)
tao=as.matrix(tau[,1:S])



## Long-run Theta matrix of segments
Theta1 = list()
for(s in 1:4){
  Theta1s = matrix(0,N,N) 
  Theta1s[lower.tri(Theta1s)]= c(parameters(f4)[(N+2):(N*(N+1)/2+1),s])
  Theta1s[upper.tri(Theta1s)] = t(Theta1s)[upper.tri(Theta1s)]
  Theta1[[s]]=Theta1s
}




## compute 0.5*zb'(Theta)zb for each basket b
Thetab = matrix(0,B,S)
for(i in 1:B){
  for(s in 1:S){
    Thetab[i,s] = .5*Z[i,]%*%Theta1[[s]]%*%Z[i,]
  }
}


## Compute the value of 0.5*zb'(Theta)zb for each basket b and each consumer h
Thetab = tao%*%t(Thetab)
colnames(Thetab)=1:B
Thetab = as.data.frame(Thetab)
Thetab$h = 1:H 
Thetab = melt(Thetab,id.vars='h')
colnames(Thetab)=c('h','b','Thetab')


## long-run alpha value of segments
alpha1 = parameters(f4)[2:(N+1),] 


## Compute the value of alpha'zb for each basket b and each consumer h
alphab = tao%*%t(alpha1)
alphab = alphab%*%t(Z)
colnames(alphab)=1:B
alphab = as.data.frame(alphab)
alphab$h = 1:H 
alphab = melt(alphab,id.vars='h')
colnames(alphab)=c('h','b','alphab') 


## long-run frequency scale
R1 = parameters(f4)[1,]- log(T)

## long-run frequency scale for each consumer
Rt = tao%*%R1
Rt = Rt 
Rt = data.frame(R1=Rt)
Rt$h = 1:H






## add price and prom to short-run data
meanprice= colMeans(price)
meanprice = Z%*%meanprice
meanprice = data.frame(meanprice = meanprice)
meanprice$b = 1:B

meanprom= colMeans(prom)
meanprom = Z%*%meanprom
meanprom = data.frame(meanprom = meanprom)
meanprom$b = 1:B

priceb = Z%*%t(price)

priceb = as.data.frame(priceb)
colnames(priceb)=1:T
priceb$b = 1:B
priceb = melt(priceb,id.vars='b')
colnames(priceb)=c('b','t','price')

promb = Z%*%t(prom)
colnames(promb)=1:T
promb = as.data.frame(promb)
promb$b = 1:B
promb = melt(promb,id.vars='b')
colnames(promb)=c('b','t','prom')

## snhb is the short-run basket data in step 2

snhb = merge(Nth,alphab,by=c('h','b'))
snhb = merge(snhb,priceb,by=c('b','t'))
snhb = merge(snhb,promb,by=c('b','t'))
snhb = merge(snhb,meanprice,by='b')
snhb = merge(snhb,meanprom,by='b')
snhb = merge(snhb,Thetab,by=c('h','b'))

snhb$price = snhb$price - snhb$meanprice
snhb$prom = snhb$prom-snhb$meanprom

snhb = merge(snhb,Rt,by='h')



snhb = merge(snhb,tau,by='h')

snhb$price1= snhb$price*snhb$V1
snhb$price2= snhb$price*snhb$V2
snhb$price3= snhb$price*snhb$V3
snhb$price4= snhb$price*snhb$V4


snhb$prom1= snhb$prom*snhb$V1
snhb$prom2= snhb$prom*snhb$V2
snhb$prom3= snhb$prom*snhb$V3
snhb$prom4= snhb$prom*snhb$V4





## final estimation


gl = glm(bn ~ 0+R1+alphab+Thetab
         +price1+price2+price3+price4+
           prom1+prom2+prom3+prom4,
         snhb,family = 'poisson')

step2_para = summary(gl)$coefficients


## results

## marketing mix coefficients
beta1s = step2_para[c('price1','price2','price3','price4'),1] ## price coefficients
beta2s = step2_para[c('prom1','prom2','prom3','prom4'), 1] ## prom coefficients

## frequency scale
R = R1*step2_para['R1',1]

## preferences
price_bias = colMeans(price)%*%t(beta1s)
prom_bias = colMeans(prom)%*%t(beta2s)
alpha0 = alpha1*step2_para['alphab',1] - price_bias-prom_bias

alphas = t(alpha0)-colMeans(alpha0)



## cross-category dependency matrix
Theta = list()
for(s in 1:S){
  Theta[[s]] = Theta1[[s]]*step2_para['Thetab',1]
}





X1 <- rnorm(100)
X2 <- rnorm(100,5)
X3 <- rnorm(100,0,8)
X4 <- rnorm(100,10)
X5 <- rnorm(100,10,3)
mat <- cbind(X1,X2,X3,X4,X5)
par(mfrow=c(1,1))
boxplot(mat,notch=TRUE)

nearestDist <- matrix(1:50,nrow = 50, ncol = 1)
for (i in 1:50){
  repetitive_nn[i] <- TSPsolve(couts, 'repetitive_nn')
  nearest_insertion[i] <- TSPsolve(couts, 'nearest_insertion')
  two_opt[i] <- TSPsolve(couts, 'two_opt')
  nearestDist[i] <-TSPnearest(couts)$longueur
  branchDist[i] <- TSPbranch(couts)
}

mat_boxplot <- cbind(repetitive_nn,nearest_insertion,two_opt,nearestDist,branchDist)
par(mfrow=c(1,1))
boxplot(mat_boxplot,notch=TRUE)

t.test(nearestDist - branchDist, paired=FALSE, mu=0)


mean(branchDist)
mean(nearestDist)


results <- vector(length = 250)
methods <- vector(length = 250)

results <- c(nearestDist, branchDist, nearest_insertion,repetitive_nn, two_opt)

for (i in 1:250){
  if (i < 51){
    methods[i] <- 'nearest'
  }
  else if (i < 101){
    methods[i] <- 'Branch&Bound'
  }
  else if (i < 151) {
    methods[i] <- 'nearest_insertion'
  }
  else if (i < 201) {
    methods[i] <- 'repetitive_nn'
  }
  else {
    methods[i] <- 'two_opt'
  }
}


pairwise.t.test(results,methods,adjust.method='bonferroni')

microbenchmark(TSPsolve(couts, 'repetitive_nn'), TSPsolve(couts, 'nearest_insertion'),
               TSPsolve(couts, 'two_opt'), TSPnearest(couts)$longueur, TSPbranch(couts),
               times=20, 
               setup={
                n <- 10
               sommets <- data.frame(x = runif(n), y = runif(n))
               couts <- distance(sommets)
               })

seqn <- seq(4,20,1)
temps <- matrix(ncol=10, nrow=length(seqn))
for (i in 1:length(seqn)){
  temps[i,] <- microbenchmark(TSPsolve(couts, method = 'branch'),
                             times = 10,
                             setup = { n <- seqn[i]
                             couts <- distance(cbind(x = runif(n), y = runif(n)))}
  )$time
}

par(mfrow=c(1,2)) # 2 graphiques sur 1 ligne
matplot(seqn, temps, xlab='n', ylab='temps')
matplot(seqn, log(temps)^2, xlab='n', ylab=expression(log(temps)^2))

vect_temps <- log(as.vector(temps))^2
vect_dim <- rep(seqn,times=10)
temps.lm <- lm(vect_temps ~ vect_dim)
summary(temps.lm)
par(mfrow=c(2,2))
plot(temps.lm)
shapiro.test(residuals(temps.lm))

temps.moy <- rowMeans(temps)
vect_temps.moy <- log(as.vector(temps.moy))^2
temps.moy.lm <- lm(vect_temps.moy ~ seqn)
summary(temps.moy.lm)
par(mfrow=c(2,2))
plot(temps.moy.lm)
shapiro.test(residuals(temps.moy.lm))


data.graph <- read.csv('DonneesTSP.csv',header=TRUE,dec='.',sep=',',quote="\"")

dataset <- data.frame(matrix(ncol = 7, nrow = 70))
colnames(dataset) <- c("sqrtdim", "mean.long","mean.dist", 
                       "sd.dist","mean.deg", "sd.deg", "diameter")
dataset$sqrtdim = sqrt(data.graph$dim)
dataset$mean.long = data.graph$mean.long
dataset$mean.dist = data.graph$mean.dist
dataset$sd.dist = data.graph$sd.dist
dataset$mean.deg = data.graph$mean.deg
dataset$sd.deg = data.graph$sd.deg
dataset$diameter = data.graph$diameter

vect_tps <- log(as.vector(data.graph$tps))^2
tps.lm <- lm(vect_tps ~., data = dataset)
summary(tps.lm)
par(mfrow=c(2,2))
plot(tps.lm)

step(tps.lm)

datasetAfterAIC <- data.frame(matrix(ncol = 6, nrow = 70))
colnames(datasetAfterAIC) <- c("sqrtdim", "mean.long", 
                       "sd.dist","mean.deg", "sd.deg", "diameter")
datasetAfterAIC$sqrtdim = sqrt(data.graph$dim)
datasetAfterAIC$mean.long = data.graph$mean.long
datasetAfterAIC$sd.dist = data.graph$sd.dist
datasetAfterAIC$mean.deg = data.graph$mean.deg
datasetAfterAIC$sd.deg = data.graph$sd.deg
datasetAfterAIC$diameter = data.graph$diameter

tps.lm_aic <- lm(vect_tps ~., data = datasetAfterAIC)
summary(tps.lm_aic)
par(mfrow=c(2,2))
plot(tps.lm_aic)

step(tps.lm_aic)

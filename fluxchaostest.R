
library (nonlinearTseries)
library (tseriesChaos)
library (plot3D)
library (deSolve) #ODE solver package
library (mblm)
library (fields)
library (doParallel)
library (snow)
#### Load functions ####
date.fun <- function(date.fluxnet){
  return(as.numeric(substring(as.character(date.fluxnet), c(1,5,7,9,11), c(4,6,8,10,12))))
}

eladio e gordo
embed_udf <-function(x, mmax=6, d=NA){
  if(is.na(d)){
    #Step 1: Calculate embedding delay
    #Step 1a. Compute average mutual information (AMI) function
    mutual.out<-mutual(x, lag.max=100, plot=F) #mutual(tseriesChaos)
    #Step 1b: Use embedding approach to calculate delay at which AMI hits
    #its first minimum
    mutual.em<-embedd(mutual.out,2,1) #embedd(tseriesChaos)
    mutual.adj.length<-length(mutual.out)-1 #lose 1 observation to delay
    mutual.hold<-matrix(0,mutual.adj.length,1)
    for (i in 1:mutual.adj.length){ #loop to compute successfive value differences
      mutual.test<-if(mutual.em[i,1]>mutual.em[i,2])TRUE else break
      mutual.hold[i,1]<-mutual.test #"TRUE' = 1 in R, so min delay occurs at sum(TRUE)
    } #end iloop
    mutual.hold.sum<-sum(mutual.hold)
    #Step 1c: Estimate embeddeding delay (d). If the mutual information function is
    #decreasing across all 100 delays, set delay at max = 100
    d<-if(mutual.hold.sum<100)mutual.hold.sum else 100 #Embedding delay
  }
  #Step 2: Compute Theiler window parameter (tw) required for false nearest
  # neighbours test
  #Step 2a: Autocorrelation function
  lag.max=100
  acf.run<-acf(x,lag.max, plot=F)
  acf.out<-acf.run$acf #array of acf values
  #Step 2b: Use embedding approach to calculate delay at which AMI hits
  #its first minimum.
  acf.em<-embedd(acf.out,2,1) #embedd(tseriesChaos)
  acf.adj.length<-length(acf.out)-1 #lose 1 observation to lag
  acf.hold<-matrix(0,acf.adj.length,1)
  for (i in 1:acf.adj.length){ #loop to compute successfive value differences
    acf.test<-if(acf.em[i,1]>acf.em[i,2])TRUE else break
    acf.hold[i,1]<-acf.test #"TRUE' = 1 in R, so min delay occurs at sum(TRUE)
  } #end iloop
  #Step 2c: Estimate Theiler window (tw)
  tw<-sum(acf.hold)
  #Step 3: Embedding dimension (m)
  #Step 3a: False nearest neighbours function
  m.max<-mmax #maximum number of embedding dimensions to consider
  fn.out<-false.nearest(x,m.max,d,tw) #false.nearest(tseriesChaos)
  fn.out[is.na(fn.out)] <- 0 #set NA in fn.out to zero
  #plot(fn.out)
  #Step 3b: Find delay at which false nearest neighbours decrease below set tolerance
  #Output vector of fnn percentages from fn.out
  fnp<-c(fn.out[1],fn.out[3],fn.out[5],fn.out[7],fn.out[9],fn.out[11])
  fnp.tol<-fnp>0.15 #If fnp greater than tolerance of 15%, T entered into fnp.tol
  fnp.tol.sum<-sum(fnp.tol) #sum up number of T's
  m<-if(fnp.tol.sum<m.max)fnp.tol.sum+1 else m.max #Embedding dimension
  #Step 4: Embed time series (Mx)
  #If m=1, embedd routine crashes due to 'subscript out of bounds' error--need to
  #guarantee an embedding dimension of at least two:
  if(m<=1){m<-2} else {m}
  Mx<-embedd(x,m,d) #embedd(tseriesChaos)
  #Results
  results.embed_udf<-list(d,m,tw,Mx)
  return(results.embed_udf)
} # modified from NLTS analysis book
attract2d <- function(var, lag, type="p", cex=0.5, col="red", burn=0){
  if(burn!=0){
    aux <- embedd(var[-c(1:burn)],3,lag)} else {aux <- embedd(var,3,lag)}
  plot(x=aux[,1],y=aux[,2], type=type, color = col, ylab=paste("Var at lag", lag), xlab="Var", cex.symbols=cex)
}
attract <- function(var, lag, type="p", bty="b2", ..., phi = 10, theta = 20, add=F, pch=NULL, lwd=NULL, cex=0.5, col="red", alpha=1, burn=0, xlim=NULL, ylim=NULL, zlim=NULL){
  if(burn!=0){
    aux <- embedd(var[-c(1:burn)],3,lag)} else {aux <- embedd(var,3,lag)}
  # scatterplot3d(x=aux[,1],y=aux[,2],z=aux[,3], type=type, color = col, ylab=paste("Var at lag", lag), xlab="Var", zlab=paste("Var at lag", 2*lag), cex.symbols=cex, xlim=xlim, ylim=ylim, zlim=zlim)
  scatter3D(x=aux[,1],y=aux[,2],z=aux[,3], ..., add=add, alpha = alpha, ticktype = "detailed", type=type, bty= bty, lwd=lwd, cex=cex, pch=pch, phi = phi, theta = theta, col = col, ylab=paste("Var at lag", lag), xlab="Var", zlab=paste("Var at lag", 2*lag), xlim=xlim, ylim=ylim, zlim=zlim)
}
attract.or <- function(var, type="p", bty="b2", ..., phi = 10, theta = 20, add=F, pch=NULL, lwd=NULL, cex=0.5, col="red", alpha=1, burn=0, xlim=NULL, ylim=NULL, zlim=NULL){
  if(burn!=0){
    var <- var[-c(1:burn), ]}
  # scatterplot3d(x=aux[,1],y=aux[,2],z=aux[,3], type=type, color = col, ylab=paste("Var at lag", lag), xlab="Var", zlab=paste("Var at lag", 2*lag), cex.symbols=cex, xlim=xlim, ylim=ylim, zlim=zlim)
  scatter3D(x=var[,1],y=var[,2],z=var[,3], ..., add=add, alpha = alpha, ticktype = "detailed", type=type, bty= bty, lwd=lwd, cex=cex, pch=pch, phi = phi, theta = theta, col = col, ylab=paste("PC2"), xlab="PC1", zlab=paste("PC3"), xlim=xlim, ylim=ylim, zlim=zlim)
}

lorenz_udf<-function(x0,y0,z0,sigma,beta,rho,t.end,delta) {
  #Initial conditions and parameters
  state<-c(x=x0,y=y0,z=z0) #initial conditions
  parameters<-c(sigma,beta,rho) #parameters giving chaotic dynamics
  #ODE model
  model<-function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      dx<- sigma*(y-x) #Lorenz equations
      dy<- x*(rho-z)-y
      dz<- x*y-beta*z
      list(c(dx,dy,dz))
    }) #end with (as.list...
  } #end model function
  #Solution
  times<-seq(0,t.end,by=delta) #integration step
  out<-ode(y=state,times=times,func=model,parms=parameters,method="lsoda")
  #Solution variables
  x<-out[,"x"];y<-out[,"y"];z<-out[,"z"]
  results<-list(x,y,z)
  return(results)
} #end user-defined function
plot_c <- function(M,y1= 0.00000001,y2 = 1, xlim=NULL, ylab=NULL, xlab=NULL, main=NULL ){

  plot(M[,1],M[,2],log="xy",ylim=c(y1,y2),type="l", xlim=xlim, ylab=ylab, xlab=xlab, main=main, las=1)
  for(i in c(2:ncol(M))) lines(M[,1],M[,i])

}
plot_b <- function(M,y1= 0.00000001,y2 = 1, xlim=NULL, ylab=NULL, xlab=NULL, main=NULL ){

  plot(M[,1],M[,2],ylim=c(y1,y2),type="l", xlim=xlim, ylab=ylab, xlab=xlab, main=main, las=1)
  for(i in c(2:ncol(M))) lines(M[,1],M[,i])

}

## Information dimension function ##
pv.fun1 <- function(var1){
  var <- na.omit(var)
  if(var1[1]==var1[2]){return(0)}else{
    return((1-(min(var1[1],var1[2])/max(var1[1],var1[2]))))}
}
pv.fun <- function(var){
  var <- na.omit(var)
  if(length(var)<2){return(NA)} else {
    return(mean(as.numeric(combn(var, 2, FUN = pv.fun1, simplify = F))))}}
id.calc <- function(d, add=0, rang=20, prop.rang.max=1, prop.rang.min=1,
                    ttail=5, cut=0.99, cut.low=0.01, w.var=0.5, dims=NA, correction="loess"){
  # initialise
  # add <- 1
  # rang <- 20
  # prop.rang.max <- 1
  # prop.rang.min <- 1
  # ttail <- 5
  # cut <- 0.99 # cut values of density higher than...
  # cut.low <- 0.01 # cut values of density lower than...
  # w.var <- 0.35 # weigh of variability vs. Rsq, 0.5 is equal contribution
  # dims <- 1:30
  # set up M matrix for analysis
  # plot_c(M,y1 = 1e-02, y2=14)
  # plot_c(cbind(M[,1],1/exp(M[,-c(1:6)])), y1 = 0.0000001, y2=1)
  # plot_c(M[,-c(2:6)], y1=8, y2=15, xlim=c(0.05,1))
  temp <- t(as.matrix(d$log.radius))
  if(!is.numeric(dims[1])){dims <- 1:(ncol(temp)-1)}
  p <- as.numeric(row.names(temp))
  M <- matrix(nrow=length(temp[,1]),ncol=length(temp[1,])+2)
  M[,1] <- p;M[,2] <- 0
  for(j in c(3:length(M[1,]))) M[,j] <- 10^(temp[,j-2])
  # plot_c(cbind(M[,1],1/M[,2:ncol(M)]), y1 = 0.1, y2=10)
  # calculate p at which variability betwen dimensions is lower #
  M <- M[,c(1,(dims+1))]
  M <- subset(M, (M[,1]<cut)&(M[,1]>cut.low))
  # var.row <- 1-apply(M[,-c(1:3)], MARGIN = 1, FUN=pv.fun)
  # var.row <- var.row/max(var.row,na.rm=T)
  # plot(var.row)
  ### starts running
  rang <- round(rang,digit=0)
  ttail <- round(ttail, digit=0)
  prop.rang <- seq(prop.rang.max, prop.rang.min, length=(ncol(M)-1))
  sel <- list()
  cdims <- data.frame(embdim=c(1:(ncol(M)-1)),id=rep(NA,times=(ncol(M)-1)), eps=NA, rsq=NA, locus=NA)
  b <- matrix(ncol=(ncol(M)-1), nrow=nrow(M))
  c <- b # check correlations
  bbb <- b # matrix for slopes
  ccc <- b
  M[M==0] <- NA
  for (i in 1:(ncol(M)-1)){
    for (j in (rang+1):(nrow(M)-rang)){
      if (length(which(is.na(M[(j-rang):(j+rang),(i+1)])))>=(rang/2)){b [j,i] <- NA} else {
        # b [j,i] = summary(lm(log(M[(j-rang):(j+rang),(i+1)])~ log(M[(j-rang):(j+rang),1])))$r.squared # to use local Rsq
        bbb [j,i] = 1-summary(lm(log(M[(j-rang):(j+rang),(i+1)])~ log(M[(j-rang):(j+rang),1])))$coefficients[2,1]
        b [j,i] <- bbb [j,i]} # to use local slope
    }
  }

  # calculate variability between rows (1-pv= stability)
  var.row <- 1-apply(b[,-c(1:3)], MARGIN = 1, FUN=pv.fun) # estimate stability (1-PV)
  var.row [is.infinite(var.row)] <- NA
  var.row <- as.numeric((scale(var.row)+(abs(min(scale(var.row), na.rm=T))+1))/max(scale(var.row)+(abs(min(scale(var.row), na.rm=T))+1),na.rm=T))

  for (i in 1:(ncol(M)-1)){
    c [,i] <- b[,i]
    ccc [,i] <- bbb[,i]
    # b [,i] <- b[,i]*(1-w.var) * (var.row*w.var) # multiplicative weighting
    b [,i] <- b[,i]*(1-w.var) + (var.row*w.var) # additive weighting
    if (length(which(!is.na(b[,i])))<5){} else {
      cdims$rsq [i] <- max(b[,i],na.rm=T) # maximum rsq
      cdims$locus[i]<- (which(b[,i]==cdims$rsq[i])) # place max rsq
      cdims$eps [i] <- M[cdims$locus[i],1] # epsilon of max rsq
      sel[[i]] <- c((cdims$locus[i]-(round(rang*prop.rang[i],0))):(cdims$locus[i]+(round(rang*prop.rang[i],0))))
      cdims$id [i] <- 1/as.numeric(lm(log(M[sel[[i]],(i+1)]) ~ log(M[sel[[i]],1]))$coef[2])
    }
  }


  ## correcting selection ##
  sel.c <- list()
  if(correction=="none"){cdims$locus.c <- cdims$locus} # no correction
  if(correction=="loess"){cdims$locus.c <- cdims$locus
    cdims$locus.c [!is.na(cdims$id)] <- round(predict(loess(locus ~ embdim, data=cdims, span=1)),0)} # loess correction
  if(correction=="constant"){cdims$locus.c <- cdims$locus
    cdims$locus.c [(nrow(cdims)-ttail+1):nrow(cdims)] <- cdims$locus.c [(nrow(cdims)-ttail+1)]} # constant epsilon correction
  if(correction=="increasing"){cdims$locus.c <- cdims$locus
    cdims$locus.c [(nrow(cdims)-ttail+1):nrow(cdims)] <- round(seq(cdims$locus.c [(nrow(cdims)-ttail+1)], cdims$locus.c [(nrow(cdims)-ttail+1)] + add, length=ttail),0)} # constant increasing/decreasing epsilon correction
  for (i in 1:(nrow(cdims))){
    if(is.na(cdims$locus.c[i])){
      sel.c[[i]]<-NA
      cdims$rsq.c [i] <- NA
      cdims$eps.c [i] <- NA
      cdims$id.c [i] <- NA} else {
      cdims$rsq.c [i] <- b[cdims$locus.c [i],i]
      cdims$eps.c [i] <- M[cdims$locus.c [i],1]
      sel.c[[i]] <- c((cdims$locus.c[i]-(round(rang*prop.rang[i],0))):(cdims$locus.c[i]+(round(rang*prop.rang[i],0))))
      cdims$id.c [i] <- 1/as.numeric(lm(log(M[sel.c[[i]],(i+1)]) ~ log(M[sel.c[[i]],1]))$coef[2])
    }
  }


  par(mfrow=c(4,1), mar=c(3,4.5,1,1))
  # plot_c(cbind(M[,1],c*100), y1 = 85, y2=(max(c*100,na.rm=t)+5), ylab="R-squared", xlim=c(0.05,1))
  # plot_c(cbind(M[,1],b*100), y1 = 85, y2=(max(b*100,na.rm=t)+5), ylab="R-squared*Var", xlim=c(0.05,1))
  plot_c(cbind(M[,1],c*100), y1 = 80, y2=(max(c*100,na.rm=t)+5), ylab="Ln-Ln Slope", xlim=c(0.03,1))
  plot_c(cbind(M[,1],b*100), y1 = 80, y2=100, ylab="Ln-Ln Slope*Var", xlim=c(0.03,1))
  for (i in 2:nrow(cdims)){
    lines(rep(cdims$rsq.c[i]*100, times=1+2*prop.rang[i]*rang) ~ M[c((cdims$locus.c[i]-prop.rang[i]*rang):(cdims$locus.c[i]+prop.rang[i]*rang)),1], col="red", lwd=1)
  }
  points(cdims$rsq*100 ~ cdims$eps, bg="blue", pch=21)
  points(cdims$rsq.c*100 ~ cdims$eps.c, bg="red", pch=21)
  plot(id~embdim, data=cdims,las=1, bg="blue", pch=21, ylab="Info. dimension (id)", xlab="Embeddng dimension")
  plot(id.c~embdim, data=cdims,las=1, bg="red", pch=21, ylab="Info. dimension (id)", xlab="Embeddng dimension")

  mm <-summary.lm(mblm(id.c ~ embdim, data=cdims[(nrow(cdims)-ttail+1):nrow(cdims),], repeated=T))$coefficients[2,c(1,2,4)]
  id <- data.frame(id=mean(cdims$id.c[(nrow(cdims)-ttail+1):nrow(cdims)]),id.se=sd(cdims$id.c[(nrow(cdims)-ttail+1):nrow(cdims)])/sqrt(ttail))

  lines(rep(id$id,ttail) ~ c((nrow(cdims)-ttail+1):nrow(cdims)), lwd=4, col="blue")

  if(mm[3]>=0.05){
    cond <- paste("Plateau reached at ", round(id$id,2), " +/- ", round(id$id.se,2), "; Pval=", round(mm[3],5), sep="")
  } else {
    cond <- paste("Plateau was not reached, sorry... ", round(id$id,2), " +/- ", round(id$id.se,2), "; Pval=", round(mm[3],5), sep="")
  }
  print(cond)


  # based on density function
  aux <- numeric()
  for (z in 2:ncol(ccc)){
    # plot(density(na.omit(c[,z]), bw = 0.02), xlim=c(0.5,1))
    a <- density(na.omit(ccc[,z]), bw = 0.02)
    a$y <- a$y[a$x>(0.6)&a$x<(1)]
    a$x <- a$x[a$x>(0.6)&a$x<(1)]
    aa <- data.frame(x=a$x, y=a$y)
    aa <- aa[order(aa$x, decreasing = T),]
    aaa <- data.frame(embedd(aa$y, d=1, m=2))
    aaa [(nrow(aaa)+1), c(1,2)] <- NA
    aaa$res <- aaa$V1.1 -aaa$V1.0
    aa <- cbind(aa, aaa)
    site <- aa$x[which(aa$res<0)][1]
    if (class(site)=="integer"){site <- aa$x[which(aa$y==max(aa$y))]}
    aux [z] <- 1/(1-site)
  }

  aabb <- data.frame(a=aux[(length(aux)-4):length(aux)], b=c(1:5))
  mm <-summary.lm(mblm(a ~ b,dataframe = aabb, repeated=T))$coefficients[2,c(1,2,4)]
  aux2 <- data.frame(id=mean(aux[(length(aux)-4):length(aux)]),se=sd(aux[(length(aux)-4):length(aux)])/sqrt(5), pval=as.numeric(mm[3]))
  par(mfrow=c(2,1), mar=c(3,4.5,1,1))
  plot(aux, pch=21, bg="red", las=1, ylab="Info. Dimension (D1)", xlab="Embedding dimension")
  lines(rep(aux2$id,5) ~ c((max(dims)-4):max(dims)), col="blue", lwd=2)
  plot_c(cbind(M[,1],(1/(1-ccc[,-1]))), y1 = 0.5, y2=(max((1/(1-ccc[,-1])),na.rm=t)+0.5), ylab="Ln-Ln Slope", xlim=c(0.03,1))
  return(list(slopes= cbind(M[,1],(1/(1-ccc[,-1]))), rsq.method=cbind(id,pval=round(mm[3],5)), dens.method=aux2))
}
## Correlation dimension function ##
cd.calc <- function(d, exp.smooth=F, add=0, rang=100, prop.rang.max=0.5, prop.rang.min=0.5, ttail=5, dims=NA, w.var=0.5, tails=F, nontails=13:18, correction="loess"){ # correction="loess", "constant", "increasing"
  # add <- 0 # try 100 for nee
  # rang <- 10 # fix 100 for nee
  # prop.rang.max <- 0.5 # 0.5 equals rang*2 + 1
  # prop.rang.min <- 0.5
  # ttail <- 5
  # w.var <- 0.5 # variability weight, 0.5 means equal
  if(is.na(dims[1])){dims <- 1:(ncol(d)-1)}
  # tails <- F
  # nontails <- c(13:17)
  if(tails == F) {cttail <- ((length(dims)-ttail+1):length(dims))} else {cttail <- nontails}

  ## smooth scalling exponents?
  if(exp.smooth==T){
    for (i in 1:ncol(d)){
      d[,i] <- predict(loess(d[,i]~c(1:length(d[,i])),span = 0.15))
    }
  }

  ### starts running
  x <- d[,c(1,(dims+1))]
  rang <- round(rang,digit=0)
  ttail <- round(ttail, digit=0)
  prop.rang <- seq(prop.rang.max, prop.rang.min, length=(ncol(x)-1))
  sel <- list()
  cdims <- data.frame(embdim=c(1:(ncol(x)-1)),d2=rep(NA,times=(ncol(x)-1)), eps=NA, rsq=NA, locus=NA)
  b <- matrix(ncol=(ncol(x)-1), nrow=nrow(x))
  c <- b
  f <- matrix(ncol=ncol(c), nrow=nrow(c))
  f[,c(1,2)] <- 0
  x[x==0] <- NA
  for (i in 1:(ncol(x)-1)){
    for (j in (rang+1):(nrow(x)-rang)){
      if (length(which(is.na(x[(j-rang):(j+rang),(i+1)])))>=(rang/2)){b [j,i] <- NA} else {
        md <- summary(lm(log(x[(j-rang):(j+rang),(i+1)])~ log(x[(j-rang):(j+rang),1])))
        b [j,i] = md$r.squared
        c [j,i] = md$coefficients[2,1]}
    }
  }

  # for (y in 2:ncol(c)){
  #   f [,y] <- abs(c[,y]-c[,y-1])/c[,y-1]
  # } # substraction

  # calculate variability (estability=1-PV) per epsilon for ech band
  for (i in 3:ncol(c)){
    for (w in 1:nrow(c)){
      if(any(is.na(c(c[w,i],c[w,i-1],c[w,i-2])))) {f [w,i] <- 0} else
      {f [w,i] <- pv.fun(c(c[w,i],c[w,i-1],c[w,i-2]))}
    }
  } # PV 3 bands

  # # calculate stability (1-PV) per epsilon all bands together
  # for (w in 1:nrow(c)){
  #   if(any(is.na(c(c[w,])))) {f [w,1] <- NA} else
  #   {f [w,1] <- pv.fun(c(c[w,]))}
  # }# PV all bands
  # f [,-1] <- f[,1]

  # correct matrices for var.weight
  f <- (1-f)*w.var # convert to stability
  summary (f)
  b <- b*(1-w.var)

  # selection of site
  for (i in 1:(ncol(x)-1)){
    # b[!is.finite(b[,i]),i] <- NA
    # maximum of the series
    cdims$rsq [i] <- max((f[,i])+b[,i],na.rm=T) # maximum rsq
    cdims$locus[i]<- (which((f[,i])+b[,i]==cdims$rsq[i])) # place max rsq
    cdims$eps [i] <- x[cdims$locus[i],1] # epsilon of max rsq
    sel[[i]] <- c((cdims$locus[i]-(round(rang*prop.rang[i],0))):(cdims$locus[i]+(round(rang*prop.rang[i],0))))
    cdims$d2 [i] <- as.numeric(lm(log(x[sel[[i]],(i+1)]) ~ log(x[sel[[i]],1]))$coef[2])
  }

  ## correcting selection ##
  sel.c <- list()
  cdims$locus.c <- round(predict(loess(locus ~ embdim, data=cdims, span=1)),0) # smooth-only correction
  if(correction=="loess"){cdims$locus.c <- cdims$locus} # no correction
  if(correction=="constant"){cdims$locus.c [(ncol(b)-ttail+1):ncol(b)] <- cdims$locus.c [(ncol(b)-ttail+1)]} # constant epsilon correction
  if(correction=="increasing"){cdims$locus.c [(ncol(b)-ttail+1):ncol(b)] <- round(seq(cdims$locus.c [(ncol(b)-ttail+1)], cdims$locus.c [(ncol(b)-ttail+1)] + add, length=ttail),0)} # constant increasing/decreasing epsilon correction
  for (i in 1:(ncol(x)-1)){
    cdims$rsq.c [i] <- b[cdims$locus.c [i],i]+(f[cdims$locus.c [i],i])
    cdims$eps.c [i] <- x[cdims$locus.c [i],1]
    sel.c[[i]] <- c((cdims$locus.c[i]-(round(rang*prop.rang[i],0))):(cdims$locus.c[i]+(round(rang*prop.rang[i],0))))
    cdims$d2.c [i] <- as.numeric(lm(log(x[sel.c[[i]],(i+1)]) ~ log(x[sel.c[[i]],1]))$coef[2])
  }


  # D2 based on density function
  aux <- numeric()
  dims <- dims
  lim <- 0.05
  L <- c; L[L<lim] <- NA; L <- na.omit(L)
  for (z in 3:length(dims)){
    # plot(density(na.omit(c[,z]), bw = 0.08))
    a <- density(L[,(z-2):z], bw = 0.07)
    aa <- data.frame(x=a$x, y=a$y)
    # for first maximum density
    # aa <- aa[order(aa$x, decreasing = T),]
    # aaa <- data.frame(embedd(aa$y, d=2, m=2))
    # aaa [(nrow(aaa)+1):(nrow(aaa)+2), c(1,2)] <- NA
    # aaa$res <- round((aaa$V1.2 - aaa$V1.0)/aaa$V1.0, digit=3)
    # aa <- cbind(aa, aaa)
    # site <- aa$x[which(aa$res<(-0.01))][1] # chosing first maximum
    # if (class(site)=="integer"){site <- aa$x[which(aa$y==max(aa$y))]}
    # for maximum density
    site <- aa$x[which(aa$y==max(aa$y))] # slope at maximum density
    aux [z] <- site

    if (z==length(dims)){
      t <- density(L, bw = 0.08)
      tt <- data.frame(x=t$x, y=t$y)
      # for first maximum
      # tt <- tt[order(tt$x, decreasing = T),]
      # ttt <- data.frame(embedd(tt$y, d=2, m=2))
      # ttt [(nrow(ttt)+1):(nrow(ttt)+2), c(1,2)] <- NA
      # ttt$res <- round((ttt$V1.2 - ttt$V1.0)/ttt$V1.0, digit=3)
      # tt <- cbind(tt, ttt)
      # autt <- mean(tt$x[(which(tt$res<(-0.01))[1]-2):(which(tt$res<(-0.01))[1]+2)])
      # autt.se <- sd(tt$x[(which(tt$res<(-0.01))[1]-2):(which(tt$res<(-0.01))[1]+2)])/sqrt(5)
      # if (class(site)=="integer"){autt <- tt$x[which(tt$y==max(tt$y))]}

      # for maximum density
      autt <- mean(tt$x[(which(tt$y==max(tt$y))-2):(which(tt$y==max(tt$y))+2)])
      autt.se <- sd(tt$x[(which(tt$y==max(tt$y))-2):(which(tt$y==max(tt$y))+2)])/sqrt(5)
    }
  }
  aabb <- data.frame(a=aux[(length(aux)-4):length(aux)], b=c(1:5))
  mm <-summary.lm(mblm(a ~ b,dataframe = aabb, repeated=T))$coefficients[2,c(1,2,4)]
  aux2 <- data.frame(cd=c(mean(aux[(length(aux)-4):length(aux)]),autt),se=c(sd(aux[(length(aux)-4):length(aux)])/sqrt(5),autt.se), pval=c(as.numeric(mm[3]),NA))
  rownames(aux2) <- c("Single", "Overall")

  ## minimum slope (5) ##
  slataux <- NA
  plataux <- NA
  for (i in 5:length(cdims$embdim)){
    aaa <- summary.lm(mblm(d2.c ~ embdim,dataframe = cdims[c((i-4):i),], repeated=T))$coefficients[2,]
    slataux[i] <- abs(as.numeric(aaa[1]))
    plataux[i] <- as.numeric(aaa[4])
  }
  d2.minslope <- mean(cdims$d2.c[(which(slataux==min(slataux, na.rm=T))-4):which(slataux==min(slataux, na.rm=T))])
  d2.minslope.se <- sd(cdims$d2.c[(which(slataux==min(slataux, na.rm=T))-4):which(slataux==min(slataux, na.rm=T))])/sqrt(5)
  d2.minslope.p <- plataux[which(slataux==min(slataux, na.rm=T))]

  # graphic output
  par(mfrow=c(5,1), mar=c(3,4.5,1,1))
  plot_b(na.omit(cbind(x[,1],(b+f)*100)), y1 = 60, y2=100, ylab="R-squared x Stability (1-PV)")
  for (i in 2:nrow(cdims)){
    lines(rep(cdims$rsq.c[i]*100, times=1+2*prop.rang[i]*rang) ~ x[c((cdims$locus.c[i]-prop.rang[i]*rang):(cdims$locus.c[i]+prop.rang[i]*rang)),1], col="red", lwd=1)
  }
  points(cdims$rsq*100 ~ cdims$eps, bg="blue", pch=21)
  points(cdims$rsq.c*100 ~ cdims$eps.c, bg="red", pch=21)
  plot_b(na.omit(cbind(x[,1],c)), y1 = 0.01, y2=max(c, na.rm=T), ylab="Local Exponents (slope)")
  plot(d2~embdim, data=cdims,las=1, bg="blue", pch=21, ylab="Correlation dimension (D2)", xlab="Embeddng dimension")
  plot(d2.c~embdim, data=cdims,las=1, bg="red", pch=21, ylab="Corrected D2", xlab="Embeddng dimension")

  # estimate density of D2
  dn <- density(cdims$d2.c);dn <- dn$x[which(dn$y==max(dn$y))]
  # estimate CD
  mm <-summary.lm(mblm(d2.c ~ embdim, data=cdims[cttail,], repeated=T))$coefficients[2,c(1,2,4)]
  d2 <- data.frame(d2=mean(cdims$d2.c[cttail]),d2.se=sd(cdims$d2.c[cttail])/sqrt(length(cttail)))
  # plot line on graph
  lines(rep(d2$d2,length(cttail)) ~ cttail, lwd=4, col="blue")
  lines(rep(dn,length(dims)) ~ c(1:length(dims)), lwd=4, col="black")
  lines(rep(d2.minslope,5) ~ c((which(slataux==min(slataux, na.rm=T))-4):which(slataux==min(slataux, na.rm=T))), lwd=4, col="red")

  # Plot density estimation of D2
  plot(aux, pch=21, bg="red", las=1, ylab="Corr. Dimension (D2)", xlab="Embedding dimension")
  lines(rep(aux2$cd[1],5) ~ c((length(dims)-4):length(dims)), col="blue", lwd=2)
  lines(rep(aux2$cd[2],length(1:max(dims))) ~ c(1:max(dims)), col="blue", lwd=2)

  # print message
  if(mm[3]>=0.05){
    cond <- paste("End plateau reached at ", round(d2$d2,2), " +/- ", round(d2$d2.se,2), "; Pval=", round(mm[3],5), sep="")
  } else {
    cond <- paste("End plateau was not reached, sorry... ", round(d2$d2,2), " +/- ", round(d2$d2.se,2), "; Pval=", round(mm[3],5), sep="")
  }
  print(cond)

  # print message previous plateau?
  if(d2.minslope.p>=0.05){
    cond <- paste("Min slope plateau reached at ", round(d2.minslope,2), " +/- ", round(d2.minslope.se,2), "; Pval=", round(d2.minslope.p,5), sep="")
  } else {
    cond <- paste("Plateau was not reached, minimum slope at ", round(d2.minslope,2), " +/- ", round(d2.minslope.se,2), "; Pval=", round(d2.minslope.p,5), sep="")
  }
  print(cond)

  ##print result rsq method
  # cbind(d2,pval=round(mm[3],5), dens.d2=dn, pv=round(100*pv.fun(cdims$d2.c[cttail]), digit=1))
  ## print results density estimation of D2
  # cbind(aux2[1,], aux2[2,-3], pv=round(100*pv.fun(aux[(length(aux)-4):length(aux)]), digit=1))
  par(mfrow=c(1,1), mar=c(3,4.5,1,1))
  plot_b(na.omit(cbind(x[,1],c)), y1 = 0.01, y2=max(c, na.rm=T), ylab="Local Exponents (slope)")
  return(list(slopes=x,
              param=list(scal.exp=d, c=c, b=b, f=f, prop.rang=prop.rang, rang=rang, cttail=cttail,
                         dims=dims, aux=aux, aux2=aux2, slataux=slataux), d2=d2, dn=dn, cdims=cdims, d2.minslope=d2.minslope,
              rsq.method=data.frame(d2, slope=round(mm[1],3), pval=round(mm[3],5), d2.minslope=d2.minslope, d2.minslope.se=d2.minslope.se, d2.minslope.p=d2.minslope.p, dens.d2=dn, pv=round(100*pv.fun(cdims$d2.c[cttail]), digit=1)),
              dens.method=data.frame(cd.dens.single=aux2[1,1], se.single=aux2[1,2],pval=aux2[1,3], cd.dens.overall=aux2[2,1], se.overall=aux2[2,2], pv=round(100*pv.fun(aux[(length(aux)-4):length(aux)]), digit=1))))
}
cd.mfm <- function(var, tau="auto", delay.technique="ami", sel.method="first.minimum", extra.m=10, tw="auto", prop.rang.max=0.5, prop.rang.min=0.5, neps=1000, rang=100, eps.min=0.01){ # tau is embedding delay
  # calculate delay (tau) if auto
  if (tau=="auto"){
    tau <- timeLag(var, technique = delay.technique, selection.method=sel.method, lag.max = 100, do.plot = F)
    a <- print(paste("Delay (tau) =",tau))
  }
  a <- print(paste("Delay (tau) =",tau))
  attract(var [1:c(round(0.25*length(var),0))], lag=tau, type="p", main= a) # show attractor
  emb.dim <- estimateEmbeddingDim(var, time.lag = tau, max.embedding.dim = 25, do.plot=F) # calculating embeding dimension #
  if (tw=="auto"){
    #Step 2: Compute Theiler window parameter (tw) required for false nearest
    # neighbours test
    #Step 2a: Autocorrelation function
    lag.max=100
    acf.run<-acf(var,lag.max, plot=F)
    acf.out<-acf.run$acf #array of acf values
    #Step 2b: Use embedding approach to calculate delay at which AMI hits
    #its first minimum.
    acf.em<-embedd(acf.out,2,1) #embedd(tseriesChaos)
    acf.adj.length<-length(acf.out)-1 #lose 1 observation to lag
    acf.hold<-matrix(0,acf.adj.length,1)
    for (i in 1:acf.adj.length){ #loop to compute successfive value differences
      acf.test<-if(acf.em[i,1]>acf.em[i,2])TRUE else break
      acf.hold[i,1]<-acf.test #"TRUE' = 1 in R, so min delay occurs at sum(TRUE)
    } #end iloop
    #Step 2c: Estimate Theiler window (tw)
    tw<-sum(acf.hold)
  }
  cdd <- d2(var, m = emb.dim+extra.m, d = tau, t = tw, eps.min = eps.min, neps=neps) # run d2 routine
  return(cd_calc=cd.calc(cdd, dims=1:(emb.dim+extra.m), prop.rang.max = prop.rang.max, prop.rang.min = prop.rang.min, rang = rang))
}
cd.calc.plots <- function(cd.calc.obj){
  x <- cd.calc.obj$slopes
  c <- cd.calc.obj$param$c
  b <- cd.calc.obj$param$b
  f <- cd.calc.obj$param$f
  prop.rang <- cd.calc.obj$param$prop.rang
  rang <- cd.calc.obj$param$rang
  cttail <- cd.calc.obj$param$cttail
  slataux <- cd.calc.obj$param$slataux
  d2.minslope <- cd.calc.obj$d2.minslope
  aux <- cd.calc.obj$param$aux
  aux2 <- cd.calc.obj$param$aux2
  cdims <- cd.calc.obj$cdims
  dims <- cd.calc.obj$param$dims
  dn <- cd.calc.obj$dn
  d2 <- cd.calc.obj$d2
  # graphic output
  par(mfrow=c(5,1), mar=c(3,4.5,1,1))
  plot_b(na.omit(cbind(x[,1],(b+f)*100)), y1 = 60, y2=100, ylab="R-squared x Stability (1-PV)")
  for (i in 2:nrow(cdims)){
    lines(rep(cdims$rsq.c[i]*100, times=1+2*prop.rang[i]*rang) ~ x[c((cdims$locus.c[i]-prop.rang[i]*rang):(cdims$locus.c[i]+prop.rang[i]*rang)),1], col="red", lwd=1)
  }
  points(cdims$rsq*100 ~ cdims$eps, bg="blue", pch=21)
  points(cdims$rsq.c*100 ~ cdims$eps.c, bg="red", pch=21)
  plot_b(na.omit(cbind(x[,1],c)), y1 = 0.01, y2=max(c, na.rm=T), ylab="Local Exponents (slope)")
  plot(d2~embdim, data=cdims,las=1, bg="blue", pch=21, ylab="Correlation dimension (D2)", xlab="Embeddng dimension")
  plot(d2.c~embdim, data=cdims,las=1, bg="red", pch=21, ylab="Corrected D2", xlab="Embeddng dimension")
  # plot line on graph
  lines(rep(d2$d2,length(cttail)) ~ cttail, lwd=4, col="blue")
  lines(rep(dn,length(dims)) ~ c(1:length(dims)), lwd=4, col="black")
  lines(rep(d2.minslope,5) ~ c((which(slataux==min(slataux, na.rm=T))-4):which(slataux==min(slataux, na.rm=T))), lwd=4, col="red")

  # Plot density estimation of D2
  plot(aux, pch=21, bg="red", las=1, ylab="Corr. Dimension (D2)", xlab="Embedding dimension")
  lines(rep(aux2$cd[1],5) ~ c((length(dims)-4):length(dims)), col="blue", lwd=2)
  lines(rep(aux2$cd[2],length(1:max(dims))) ~ c(1:max(dims)), col="blue", lwd=2)

  ##print result rsq method
  # cbind(d2,pval=round(mm[3],5), dens.d2=dn, pv=round(100*pv.fun(cdims$d2.c[cttail]), digit=1))
  ## print results density estimation of D2
  # cbind(aux2[1,], aux2[2,-3], pv=round(100*pv.fun(aux[(length(aux)-4):length(aux)]), digit=1))
  par(mfrow=c(1,1), mar=c(3,4.5,1,1))
  plot_b(na.omit(cbind(x[,1],c)), y1 = 0.01, y2=max(c, na.rm=T), ylab="Local Exponents (slope)")

}
summary.cd.calc <- function(cd.calc.obj){
  return(data.frame(cd.calc.obj$rsq.method,
                    cd.calc.obj$dens.method))
}
#Code 7.3 Compute correlation dimension for time series using Takens estimator
cd_takens<-function(x,m_max) {
  cd.m<-matrix(0,m_max,1) #store estimated correlation dimensions
  for(j in 2:m_max) { #iterate over embedding dimensions (m)
    #Step 1: Compute embedding delay and Theiler window
    embed<-embed_udf(x)
    d<-embed[[1]] #embedding delay
    tw<-embed[[3]] #Theiler window
    #Step 2: Compute embedded data matrix
    library(tseriesChaos)
    Mx<-embedd(x,j,d) #embedded data matrices
    n<-nrow(Mx)
    #Step 3: Compute distance matrix
    library(fields)
    dist<-rdist(Mx) #Distance matrix
    dist.vec<-as.vector(dist)
    #Step 4: Adjust distance matrix to Theiler window
    #Step 4a: reduce to relevant rows and columns
    dist.1<-dist[(tw+1):n,1:(n-tw)]
    #Step 4b: reduce to relevant distances in each column
    dist.2<-matrix(0,(n-tw),(n-tw)) #storage matrix for loop
    for(i in 1:(n-tw)) { #iterate over columns of dist.1
      column<-dist.1[,i][i:(n-tw)]
      length(column)<-n-tw #make reduced columns same length for storage
      dist.2[,i]<-column
    } #end loop i
    #Step 5: Put distances in matrix into vector
    dist.vec1<-as.vector(dist.2) #distances not satisfying Theiler window
    #have NAs
    dist.vec2<-na.omit(dist.vec1) #remove NAs, only distances satisfying
    #Theiler window remain
    if(tw==0) {distances<-dist.vec} else {distances<-dist.vec2}
    #Step 6: Set upper threshold on distance
    x.thresh<-0.5*sqrt(sum(x^2)/length(x)) #1/2 root mean squared error (rms)
    distances.thresh<-distances[-c(which(distances>x.thresh))] #omit distances over threshold
    #Step 7: Takens estimation of correlation dimension
    cd<--1/mean(log(distances.thresh/x.thresh))
    cd.m[j,]<-cd #vector storing correlation dimensions over range of
    #embedding dimensions
  } #end loop j
  return(list(cd=cd.m, delay=d, tw=tw))
} #end user-defined function
cd_fertakens<-function(x, m_min=2, m_max="auto", d="ami", tw="auto", sel.meth="first.e.decay") {
  #Step 1: Compute embedding delay (d), max embedding dimension (m_max) and Theiler window
  if(d=="ami"){
    d <- timeLag(x, technique = "ami", selection.method=sel.meth, lag.max = 100, do.plot = F)
  }
  if(d=="acf"){
    d <- timeLag(x, technique = "acf", selection.method=sel.meth, lag.max = 100, do.plot = F)
  }
  # calculate embedding dimension if m_max="auto"
  if (m_max=="auto"){
    m_max <- 10+estimateEmbeddingDim(x, time.lag = d, max.embedding.dim = 25, do.plot=F) # calculating embeding dimension #
  }
  cd.m<-matrix(0,m_max,1) #store estimated correlation dimensions
  if (tw=="auto"){
    #Step 2: Compute Theiler window parameter (tw) required for false nearest
    # neighbours test
    #Step 2a: Autocorrelation function
    lag.max=100
    acf.run<-acf(x,lag.max, plot=F)
    acf.out<-acf.run$acf #array of acf values
    #Step 2b: Use embedding approach to calculate delay at which AMI hits
    #its first minimum.
    acf.em<-embedd(acf.out,2,1) #embedd(tseriesChaos)
    acf.adj.length<-length(acf.out)-1 #lose 1 observation to lag
    acf.hold<-matrix(0,acf.adj.length,1)
    for (i in 1:acf.adj.length){ #loop to compute successfive value differences
      acf.test<-if(acf.em[i,1]>acf.em[i,2])TRUE else break
      acf.hold[i,1]<-acf.test #"TRUE' = 1 in R, so min delay occurs at sum(TRUE)
    } #end iloop
    #Step 2c: Estimate Theiler window (tw)
    tw<-sum(acf.hold)
  }
  #iterate over embedding dimensions (m)
  for(j in m_min:m_max) {
    #Step 2: Compute embedded data matrix
    Mx<-embedd(x,j,d) #embedded data matrices
    n<-nrow(Mx)
    #Step 3: Compute distance matrix
    dist<-rdist(Mx) #Distance matrix
    dist.vec<-as.vector(dist)
    #Step 4: Adjust distance matrix to Theiler window
    #Step 4a: reduce to relevant rows and columns
    dist.1<-dist[(tw+1):n,1:(n-tw)]
    #Step 4b: reduce to relevant distances in each column
    dist.2<-matrix(0,(n-tw),(n-tw)) #storage matrix for loop
    for(i in 1:(n-tw)) { #iterate over columns of dist.1
      column<-dist.1[,i][i:(n-tw)]
      length(column)<-n-tw #make reduced columns same length for storage
      dist.2[,i]<-column
    } #end loop i
    #Step 5: Put distances in matrix into vector
    dist.vec1<-as.vector(dist.2) #distances not satisfying Theiler window
    #have NAs
    dist.vec2<-na.omit(dist.vec1) #remove NAs, only distances satisfying
    #Theiler window remain
    if(tw==0) {distances<-dist.vec} else {distances<-dist.vec2}
    #Step 6: Set upper threshold on distance
    x.thresh<-0.5*sqrt(sum(x^2)/length(x)) #1/2 root mean squared error (rms)
    distances.thresh<-distances[-c(which(distances>x.thresh))] #omit distances over threshold
    #Step 7: Takens estimation of correlation dimension
    cd<--1/mean(log(distances.thresh/x.thresh))
    cd.m[j,]<-cd #vector storing correlation dimensions over range of
    #embedding dimensions
  } #end loop j

  ## minimum slope (5) ##
  cdims <- data.frame(embdim=c(1:length(cd.m)), cd=as.numeric(cd.m))
  slataux <- NA
  plataux <- NA
  for (i in 5:length(cdims$embdim)){
    aaa <- summary.lm(mblm(cd ~ embdim,dataframe = cdims[c((i-4):i),], repeated=T))$coefficients[2,]
    slataux[i] <- abs(as.numeric(aaa[1]))
    plataux[i] <- as.numeric(aaa[4])
  }
  d2.minslope <- mean(cdims$cd[(which(slataux==min(slataux, na.rm=T))-4):which(slataux==min(slataux, na.rm=T))])
  d2.minslope.se <- sd(cdims$cd[(which(slataux==min(slataux, na.rm=T))-4):which(slataux==min(slataux, na.rm=T))])/sqrt(5)
  d2.minslope.p <- plataux[which(slataux==min(slataux, na.rm=T))]
  plot(cd~embdim, data=cdims,las=1, bg="red", pch=21, ylab="Correlation dimesion (D2)", xlab="Embeddng dimension")

  # estimate CD
  mm <-summary.lm(mblm(cd ~ embdim, data=tail(cdims,5), repeated=T))$coefficients[2,c(1,2,4)]
  d2 <- data.frame(d2=mean(tail(cdims$cd,5)),d2.se=sd(tail(cdims$cd,5))/sqrt(length(5)), d2.slope= as.numeric(mm[1]), pval=as.numeric(mm[3]))
  # plot line on graph
  lines(rep(d2$d2,5) ~ tail(cdims$embdim,5), lwd=4, col="blue")
  lines(rep(d2.minslope,5) ~ c((which(slataux==min(slataux, na.rm=T))-4):which(slataux==min(slataux, na.rm=T))), lwd=4, col="red")

  return(list("cds"=cdims,
              "cd"=data.frame(d2.tail=d2$d2, d2.tail.se=d2$d2.se, d2.tail.slope.pval=d2$pval,
                              d2.minslope=d2.minslope, d2.minslope.se=d2.minslope.se, d2.minslope.p=d2.minslope.p),
              "delay"=d, "tw"=tw))
} #end user-defined function
lorenz_udf<-function(x0,y0,z0,sigma,beta,rho,t.end,delta) {
  #Initial conditions and parameters
  state<-c(x=x0,y=y0,z=z0) #initial conditions
  parameters<-c(sigma,beta,rho) #parameters giving chaotic dynamics
  #ODE model
  model<-function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      dx<- sigma*(y-x) #Lorenz equations
      dy<- x*(rho-z)-y
      dz<- x*y-beta*z
      list(c(dx,dy,dz))
    }) #end with (as.list...
  } #end model function
  #Solution
  times<-seq(0,t.end,by=delta) #integration step
  library(deSolve) #ODE solver package
  out<-ode(y=state,times=times,func=model,parms=parameters,method="lsoda")
  #Solution variables
  x<-out[,"x"];y<-out[,"y"];z<-out[,"z"]
  results<-list(x,y,z)
  return(results)
} #end user-defined function
embed_udf <-function(x){
  #Step 1: Calculate embedding delay
  #Step 1a. Compute average mutual information (AMI) function
  library(tseriesChaos)
  mutual.out<-mutual(x) #mutual(tseriesChaos)
  #Step 1b: Use embedding approach to calculate delay at which AMI hits
  #its first minimum
  mutual.em<-embedd(mutual.out,2,1) #embedd(tseriesChaos)
  mutual.adj.length<-length(mutual.out)-1 #lose 1 observation to delay
  mutual.hold<-matrix(0,mutual.adj.length,1)
  for (i in 1:mutual.adj.length){ #loop to compute successfive value differences
    mutual.test<-if(mutual.em[i,1]>mutual.em[i,2])TRUE else break
    mutual.hold[i,1]<-mutual.test #"TRUE' = 1 in R, so min delay occurs at sum(TRUE)
  } #end iloop
  mutual.hold.sum<-sum(mutual.hold)
  #Step 1c: Estimate embeddeding delay (d). If the mutual information function is
  #decreasing across all 20 delays, set delay at max = 20
  d<-if(mutual.hold.sum<20)mutual.hold.sum else 20 #Embedding delay
  #Step 2: Compute Theiler window parameter (tw) required for false nearest
  # neighbours test
  #Step 2a: Autocorrelation function
  lag.max=100
  acf.run<-acf(x,lag.max)
  acf.out<-acf.run$acf #array of acf values
  #Step 2b: Use embedding approach to calculate delay at which AMI hits
  #its first minimum.
  acf.em<-embedd(acf.out,2,1) #embedd(tseriesChaos)
  acf.adj.length<-length(acf.out)-1 #lose 1 observation to lag
  acf.hold<-matrix(0,acf.adj.length,1)
  for (i in 1:acf.adj.length){ #loop to compute successfive value differences
    acf.test<-if(acf.em[i,1]>acf.em[i,2])TRUE else break
    acf.hold[i,1]<-acf.test #"TRUE' = 1 in R, so min delay occurs at sum(TRUE)
  } #end iloop
  #Step 2c: Estimate Theiler window (tw)
  tw<-sum(acf.hold)
  #Step 3: Embedding dimension (m)
  #Step 3a: False nearest neighbours function
  m.max<-6 #maximum number of embedding dimensions to consider
  fn.out<-false.nearest(x,m.max,d,tw) #false.nearest(tseriesChaos)
  fn.out[is.na(fn.out)] <- 0 #set NA in fn.out to zero
  #plot(fn.out)
  #Step 3b: Find delay at which false nearest neighbours decrease below set tolerance
  #Output vector of fnn percentages from fn.out
  fnp<-c(fn.out[1],fn.out[3],fn.out[5],fn.out[7],fn.out[9],fn.out[11])
  fnp.tol<-fnp>0.15 #If fnp greater than tolerance of 15%, T entered into fnp.tol
  fnp.tol.sum<-sum(fnp.tol) #sum up number of T's
  m<-if(fnp.tol.sum<m.max)fnp.tol.sum+1 else m.max #Embedding dimension
  #Step 4: Embed time series (Mx)
  #If m=1, embedd routine crashes due to 'subscript out of bounds' error--need to
  #guarantee an embedding dimension of at least two:
  if(m<=1){m<-2} else {m}
  Mx<-embedd(x,m,d) #embedd(tseriesChaos)
  #Results
  results.embed_udf<-list(d,m,tw,Mx)
  return(results.embed_udf)
} #end user-defined function
season.dates <-function(st.date, nyears){
  yini <- as.numeric(substring(as.character(st.date), c(1), c(4)))
  yys <- rep(yini:(yini+(nyears-1)), each=4)
  mnt.ini <- rep(c("01010000", "04010000", "07010000", "10010000"), times=nyears)
  mnt.end <- rep(c("03312330", "06302330", "09302330", "12312330"), times=nyears)
  res <- data.frame(date.ini=as.numeric(paste(yys,mnt.ini, sep="")), date.end=as.numeric(paste(yys,mnt.end, sep="")))
  return(res)
}
season.dates.yrs <-function(yys){
  mnt.ini <- rep(c("01010000", "04010000", "07010000", "10010000"), times=length(yys))
  mnt.end <- rep(c("03312330", "06302330", "09302330", "12312330"), times=length(yys))
  yys <- rep(yys, each=4)
  res <- data.frame(date.ini=as.numeric(paste(yys,mnt.ini, sep="")), date.end=as.numeric(paste(yys,mnt.end, sep="")))
  return(res)
}


date.fun <- function(date.fluxnet){
  return(as.numeric(substring(as.character(date.fluxnet), c(1,5,7,9,11), c(4,6,8,10,12))))}


#### prova ####
library (deSolve)
library (deSolve)
library (tseriesChaos)
library (fields)

## functions to perform a logistic map simulaiton
f.x<- function(x,r){
  r*x*(1-x)
}
f.temp<-function(xinit,nstep,r){ # starting function f.temp
  xt<- numeric()
  x<- xinit
  xt[1]<- x
  for(i in 2:nstep){
    y<- f.x(x,r)
    x<- y
    xt[i]<- x}

  # plot(xt,type="b",xlab="time",ylab="x(t)",
  #      cex.lab=1.7,cex.axis=1.3,lwd=2)
  return(xt)
  #xt # comment to skip iterates
} # ending function f.temp
## function to draw the time plot ##
f.temp.graph<-function(xinit,nstep,r){ # starting function f.temp
  xt<- numeric()
  x<- xinit
  xt[1]<- x
  for(i in 2:nstep){
    y<- f.x(x,r)
    x<- y
    xt[i]<- x
  }
  plot(xt,type="b",xlab="time",ylab="x(t)",
       cex.lab=1.7,cex.axis=1.3,lwd=2)
  #xt # comment to skip iterates
} # ending function f.temp
### parameters and initial conditions
# vary r: r<- 2.8, r<- 3.2, r<- 3.5, r<- 4
xinit<- 0.01
nstep<- 10000
r<- 4
f.temp.graph(xinit,25,r)
logistmap.ts <-f.temp(xinit,nstep,r) # call up the time plot
x <- logistmap.ts
# Generate Lorenz data with user-defined function 'lorenz_udf'
lorenz_results<-lorenz_udf(x0=1,y0=1,z0=1,sigma=10,beta=8/3,rho=28,
                           t.end=100,delta=0.01)
x <- lorenz_results[[1]]
# initiating routine
tau <- timeLag(x, technique = "ami", selection.method="first.e.decay", lag.max = 100, do.plot = T)
a <- print(paste("Delay (tau) =",tau))
attract(x [1:4000], lag=tau, type="p", main= a) # show attractor
emb.dim <- estimateEmbeddingDim(x, time.lag = tau, max.embedding.dim = 30) # calculating embeding dimension #

# using cd_takens
x <- logistmap.ts
attract(x [1:4000], lag=10, type="p", main= a) # show attractor
emb.dim <- estimateEmbeddingDim(x, time.lag = 10, max.embedding.dim = 30) # calculating embeding dimension #

m_max<-10 #maximum embedding dimension
cd_m<-cd_takens(x,m_max=m_max)
cd_m
a <- Sys.time()
cd_m2<-cd_fertakens(x,m_min=2, m_max="auto", d = "ami", sel.meth = "first.e.decay") # sel.meth="first.e.decay" (default) or "first.minimum" (like in cd_takens)
cd_m2
Sys.time() - a
attract(x [1:4000], lag=cd_m$delay, type="l", main= a) # show attractor
attract(x [1:4000], lag=1, type="p", main= a) # show attractor


m<-1:m_max #revised x-axis
cd.m.plot<-plot(m,cd_m$cd,xaxt="n",xlab="",ylab="correlation dimension",
                pch=15,type="b")
#Custom x-axis, side=1 is x-axis, "at" defines 'tick-mark' numbers
axis(side=1,at=1:length(cd_m))
#Label custom x-axis, "line=2.5" puts label below tick-mark numbers
mtext(side=1,"embedding dimension",line=2,cex=0.8)
#Label points in plot
round(cd_m,digits=2)
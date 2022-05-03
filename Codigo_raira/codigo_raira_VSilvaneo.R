library(Matrix)   # Pode ser substituido pela bilbioteca Matrix.
library(tidyr)
library(ggplot2)
library(plotly)
library(dplyr)
library(latex2exp)

calcula_max=function(pre_max){
  if(length(pre_max)==0 | sum(pre_max**2)<10**-20){
    pre_max=1
  }else{
    pre_max=max(pre_max)
  }
  scaled_max=log10(pre_max)
  category=scaled_max%%1
  value=10**(floor(log10(max(pre_max))))
  if(category<0.1){
    value=value/10
  }else{
    if(category<0.25){
      value=value/5
    }else{
      if(category<0.5){
        value=value/2
      }
    }
  }
  interval_size=(pre_max%/%value)+2
  max_value=value*interval_size
  
  return(list(value,interval_size,max_value))
}

poisson_gi_exp <- function(y,m0 = 0, C0 = 1, F1,G1,D1,W1, pop, IC_prob=0.95){
  
  # Definindo quantidades
  n1 <- nrow(F1)
  
  r <- 1
  T <- length(y)
  n <- n1 
  FF <- F1
  G <- G1
  
  # D1.aux <- matrix(rep(D1,n1^2), ncol = n1 )
  # D2.aux <- matrix(rep(D2,n2^2), ncol = n2 )
  # D.aux <- as.matrix(bdiag(D1.aux, D2.aux))
  D.aux <- D1
  D <- ifelse(D.aux == 0, 1, D.aux)
  
  W <- W1
  
  # Definindo objetos
  at <- matrix(0, ncol=T, nrow=n)
  mt <- matrix(0, ncol=T, nrow=n)
  ft <- matrix(0, ncol=1, nrow=T)
  qt <- matrix(0, ncol=1, nrow=T)
  alphat <- matrix(0, ncol=1, nrow=T)
  betat <- matrix(0, ncol=1, nrow=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  Rt <- array(rep(diag(n),T),dim=c(n,n,T))
  gt = pt = a = b= 0
  rep <- 1
  pred = var.pred = icl.pred = icu.pred = matrix(0, ncol=rep, nrow=T)
  media.post = var.post = icu.post = icl.post = matrix(0, ncol=rep, nrow=T)
  a.post = b.post = eqm =0
  E.th3 = E.th4 = raiz2 = raiz1= matrix(0, ncol=rep, nrow=T)
  tau0_star <- rep(NA, l = T)
  tau1_star <- rep(NA, l = T)
  tau0 <- rep(NA, l = T)
  tau1 <- rep(NA, l = T)
  fstar <- rep(NA, l = T)
  qstar <- rep(NA, l = T)
  
  norm_ic=qnorm(1-(1-IC_prob)/2)
  
  ## Algoritmo
  
  # Priori
  
  at[,1] <- G%*%m0
  Rt[,,1] <-G%*%C0%*%(t(G))*D[,,1]+W[,,1]
  
  reduc_RFF=Rt[,,1]%*%FF[,1]
  
  # Previsão 1 passo a frente
  ft[1,] <- t(FF[,1])%*%at[,1] + pop[1]
  qt[1,] <- t(FF[,1])%*%Rt[,,1]%*%FF[,1]
  
  # Compatibilizando prioris
  
  a[1] <- (1/qt[1,])
  b[1] <- (exp(-ft[1,]+0.5*qt[1,])/(qt[1,])) 
  
  # Posteriori
  
  a.post[1] <- a[1] + y[1]
  b.post[1] <- b[1] + 1
  
  gt[1] <- log(a.post[1]/b.post[1]) + 1/(2*a.post[1])
  pt[1] <- (2*a.post[1]-1)/(2*a.post[1]^2)
  
  mt[,1] <- at[,1]+Rt[,,1]%*%FF[,1]*(gt[1]-ft[1,])*(1/(qt[1,]))
  Ct[,,1] <- Rt[,,1] - (Rt[,,1]%*%FF[,1]%*%t(Rt[,,1]%*%FF[,1]))*(1 - pt[1]/qt[1,])*(1/qt[1,])
  
  media.post[1] <- a.post[1]/(b.post[1])
  var.post[1] <- a.post[1]/(b.post[1]^2)
  icu.post[1] <- media.post[1] + norm_ic*sqrt(var.post[1])
  icl.post[1] <- media.post[1] - norm_ic*sqrt(var.post[1])
  
  # Preditiva em t = 1
  
  pred[1] <- a[1]/ b[1]
  var.pred <- a[1]*(b[1]+1)/(b[1])^2
  icl.pred[1]<-qnbinom((1-IC_prob)/2, a[1], (b[1]/(b[1] +1)))
  icu.pred[1]<-qnbinom(1-(1-IC_prob)/2, a[1], (b[1]/(b[1] +1)))
  
  # Passo 2 até t
  
  start<- proc.time()
  for(t in 2:T){
    
    # Priori
    
    at[,t] <-  G%*%mt[,t-1]
    Rt[,,t] <- G%*%Ct[,,t-1]%*%(t(G))*D[,,t]+W[,,t]
    
    reduc_RFF=Rt[,,t]%*%FF[,t]
    
    # Previsão 1 passo a frente
    
    ft[t,] <- t(FF[,t])%*%at[,t] + pop[t]
    qt[t,] <- t(FF[,t])%*%reduc_RFF
    
    # Compatibilizando prioris
    
    a[t] <- (1/qt[t,])
    b[t] <- (exp(-ft[t,]+0.5*qt[t,])/(qt[t,]))
    
    # Posteriori
    
    a.post[t] <- a[t] + y[t]
    b.post[t] <- b[t] + 1
    
    gt[t] <- log(a.post[t]/b.post[t])+1/(2*a.post[t])
    pt[t] <- (2*a.post[t]-1)/(2*a.post[t]^2)
    
    mt[,t] <- at[,t]+reduc_RFF*(gt[t]-ft[t,])*(1/(qt[t,]))
    Ct[,,t] <- Rt[,,t] - (reduc_RFF%*%t(reduc_RFF))*(1 - pt[t]/qt[t,])*(1/qt[t,])
    
    media.post[t] <- a.post[t]/b.post[t]
    var.post[t] <- a.post[t]/(b.post[t]^2)
    icu.post[t] <- media.post[t] + norm_ic*sqrt(var.post[t])
    icl.post[t] <- media.post[t] - norm_ic*sqrt(var.post[t])
    
    # Preditiva em t = 1
    
    pred[t] <- a[t]/b[t]
    var.pred <- a[t]*(b[t]+1)/(b[t])^2
    icl.pred[t] <- qnbinom((1-IC_prob)/2, a[t], (b[t]/(b[t] +1)))
    icu.pred[t] <- qnbinom(1-(1-IC_prob)/2, a[t], (b[t]/(b[t] +1)))
  }
  
  mts <- mt
  Cts <- Ct
  
  var_index=matrix(apply(Ct,3,diag),n,T)!=0
  
  for(t in (T-1):1){
    var_ref=var_index[,t]
    restricted_Rt=Rt[var_ref,var_ref,t+1]
    restricted_Ct=Ct[var_ref,var_ref,t]
    simple_Rt_inv=restricted_Ct%*%t(G[var_ref,var_ref])%*%solve(restricted_Rt)
    
    mts[var_ref,t] <- mt[var_ref,t] + simple_Rt_inv%*%(mts[var_ref,t+1] - at[var_ref,t+1])
    Cts[var_ref,var_ref,t] <- restricted_Ct - simple_Rt_inv%*%(restricted_Rt - Cts[var_ref,var_ref,t+1])%*%t(simple_Rt_inv)
  }
  
  result <- list(mt,Ct,
                 ft, qt,
                 a,b,
                 a.post, b.post,
                 FF, G, D,W,
                 pred, icl.pred, icu.pred,
                 mts, Cts ,
                 IC_prob,var.pred,exp(pop),pop)
  names(result) <- c("mt",  "Ct",
                     "ft", "qt",
                     "alpha", "beta",
                     "alpha_star", "beta_star",
                     "F", "G", "D","W",
                     "pred", "icl.pred", "icu.pred",
                     "mts", "Cts",
                     "IC_prob","var.pred",'offset','log_offset')
  return(result)
  
}

poisson_testa_par <- function(y,m0 = 0, C0 = 1, F1,G1,D1,W1, pop, IC_prob=0.95){
  
  # Definindo quantidades
  n1 <- nrow(F1)
  
  r <- 1
  T <- length(y)
  n <- n1 
  FF <- F1
  G <- G1
  
  # D1.aux <- matrix(rep(D1,n1^2), ncol = n1 )
  # D2.aux <- matrix(rep(D2,n2^2), ncol = n2 )
  # D.aux <- as.matrix(bdiag(D1.aux, D2.aux))
  D.aux <- D1
  D <- ifelse(D.aux == 0, 1, D.aux)
  
  W <- W1
  
  # Definindo objetos
  at <- matrix(0, ncol=T, nrow=n)
  mt <- matrix(0, ncol=T, nrow=n)
  ft <- matrix(0, ncol=1, nrow=T)
  qt <- matrix(0, ncol=1, nrow=T)
  alphat <- matrix(0, ncol=1, nrow=T)
  betat <- matrix(0, ncol=1, nrow=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  Rt <- array(rep(diag(n),T),dim=c(n,n,T))
  gt = pt = a = b= 0
  rep <- 1
  pred = var.pred = icl.pred = icu.pred = matrix(0, ncol=rep, nrow=T)
  media.post = var.post = icu.post = icl.post = matrix(0, ncol=rep, nrow=T)
  a.post = b.post = eqm =0
  E.th3 = E.th4 = raiz2 = raiz1= matrix(0, ncol=rep, nrow=T)
  tau0_star <- rep(NA, l = T)
  tau1_star <- rep(NA, l = T)
  tau0 <- rep(NA, l = T)
  tau1 <- rep(NA, l = T)
  fstar <- rep(NA, l = T)
  qstar <- rep(NA, l = T)
  
  norm_ic=qnorm(1-(1-IC_prob)/2)
  
  ## Algoritmo
  
  # Priori
  
  at[,1] <- G%*%m0
  Rt[,,1] <-G%*%C0%*%(t(G))*D[,,1]+W[,,1]
  
  reduc_RFF=Rt[,,1]%*%FF[,1]
  
  # Previsão 1 passo a frente
  ft[1,] <- t(FF[,1])%*%at[,1] + pop[1]
  qt[1,] <- t(FF[,1])%*%reduc_RFF
  
  # Compatibilizando prioris
  
  a[1] <- (1/qt[1,])
  b[1] <- (exp(-ft[1,]+0.5*qt[1,])/(qt[1,])) 
  
  # Posteriori
  
  a.post[1] <- a[1] + y[1]
  b.post[1] <- b[1] + 1
  
  gt[1] <- log(a.post[1]/b.post[1]) + 1/(2*a.post[1])
  pt[1] <- (2*a.post[1]-1)/(2*a.post[1]^2)
  
  mt[,1] <- at[,1]+reduc_RFF*(gt[1]-ft[1,])*(1/(qt[1,]))
  Ct[,,1] <- Rt[,,1] - (reduc_RFF%*%t(reduc_RFF))*(1 - pt[1]/qt[1,])*(1/qt[1,])
  
  # Preditiva em t = 1
  
  pred[1] <- a[1]/ b[1]
  var.pred <- a[1]*(b[1]+1)/(b[1])^2
  
  # Passo 2 até t
  
  start<- proc.time()
  for(t in 2:T){
    
    # Priori
    
    at[,t] <-  G%*%mt[,t-1]
    Rt[,,t] <- G%*%Ct[,,t-1]%*%(t(G))*D[,,t]+W[,,t]
    
    reduc_RFF=Rt[,,t]%*%FF[,t]
    
    # Previsão 1 passo a frente
    
    ft[t,] <- t(FF[,t])%*%at[,t] + pop[t]
    qt[t,] <- t(FF[,t])%*%reduc_RFF
    
    # Compatibilizando prioris
    
    a[t] <- (1/qt[t,])
    b[t] <- (exp(-ft[t,]+0.5*qt[t,])/(qt[t,])) 
    
    # Posteriori
    
    a.post[t] <- a[t] + y[t]
    b.post[t] <- b[t] + 1
    gt[t] <- log(a.post[t]/b.post[t])+1/(2*a.post[t])
    pt[t] <- (2*a.post[t]-1)/(2*a.post[t]^2)
    
    mt[,t] <- at[,t]+reduc_RFF*(gt[t]-ft[t,])*(1/(qt[t,]))
    Ct[,,t] <- Rt[,,t] - (reduc_RFF%*%t(reduc_RFF))*(1 - pt[t]/qt[t,])*(1/qt[t,])
    
    # Preditiva em t = 1
    
    pred[t] <- a[t]/b[t]
  }
  
  result <- list(mt,Ct,ft, qt, a,b, FF, G, D,W, pred,var.pred,exp(pop),pop)
  names(result) <- c("mt",  "Ct","ft", "qt",
                     "alpha", "beta","F", "G", "D","W","pred","var.pred",'offset','log_offset')
  return(result)
  
}

gera_bloco_poly <- function(order,value=1,name='Var_Poly',D=1,m0=0,C0=1,W=0){
  G=diag(order)
  t=length(value)
  FF=matrix(c(value,rep(0,(order-1)*t)),order,t,byrow = TRUE)
  if(order==2){
    G[1,2]=1
  }else{if(order>2){
    diag(G[1:(order-1),2:order])=1
  }}
  if(length(m0)<order){
    m0=rep(m0,order)
  }
  names=list()
  names[[name]]=c(1:order)
  return(list('FF'=FF,'G'=G,'D'=array(1,c(order,order,t))*D,'W'=array(diag(order),c(order,order,t))*W,'m0'=m0,'C0'=diag(order)*C0,'names'=names,'order'=order,'n'=order,'t'=t))
}
gera_bloco_sazo <- function(period,value=1,name='Var_Sazo',D=1,m0=0,C0=1,W=0){
  w=2*pi/period
  G <- matrix(c(cos(w),-sin(w),sin(w),cos(w)),2,2)
  if(length(m0)<2){
    m0=rep(m0,2)
  }
  names=list()
  names[[name]]=c(1,2)
  t=max(length(value),1)
  
  FF=matrix(c(value,rep(0,t)),2,t,byrow = TRUE)
  
  return(list('FF'=FF,'G'=G,'D'=array(1,c(2,2,length(value)))*D,'W'=array(diag(2),c(2,2,length(value)))*W,'m0'=m0,'C0'=diag(2)*C0,'names'=names,'period'=period,'n'=2,'t'=length(value)))
}
concat_bloco <- function(...){
  blocks=list(...)
  
  n=0
  t=1
  names=list()
  for(block in blocks){
    ref_names=block$names
    for(name in names(ref_names)){
      ref_names[[name]]=ref_names[[name]]+n
    }
    names=c(names,ref_names)
    if(block$t>1){
      if(block$t!=t & t>1){
        stop(paste('Error: Blocks should have same length or length equal 1. Got',block$t,'and',t))
      }
      t=block$t
    }
    n=n+block$n
  }
  for(name in names(names)){
    ref_idx=which(names(names)==name)
    n_names=length(ref_idx)
    if(n_names>1){
      names(names)[ref_idx]=paste0(names(names)[ref_idx],'_',c(1:n_names))
    }
  }
  
  FF=matrix(0,n,t)
  G=matrix(0,n,n)
  D=array(0,c(n,n,t))
  W=array(0,c(n,n,t))
  m0=c()
  C0=matrix(0,n,n)
  position=1
  for(block in blocks){
    current_range=position:(position+block$n-1)
    FF[current_range,]=block$FF
    G[current_range,current_range]=block$G
    D[current_range,current_range,]=block$D
    W[current_range,current_range,]=block$W
    m0=c(m0,block$m0)
    C0[current_range,current_range]=block$C0
    position=position+block$n
  }
  return(list('FF'=FF,'G'=G,'D'=D,'W'=W,'m0'=m0,'C0'=C0,'n'=n,'t'=t,'names'=names))
}

ajusta_modelo <- function(...,data_out,kernel=poisson_gi_exp,offset=NULL,log_offset=NULL){
  structure=concat_bloco(...)
  if(!is.null(offset) & !is.null(log_offset)){
    stop('Erro: Cannot set both offset and log_offset. Choose only one.')
  }else{if(!is.null(offset)){
    log_offset=log(offset)
  }else{if(!is.null(log_offset)){
    offset=exp(log_offset)
  }else{
    offset=1
    log_offset=0
  }}}
  if(1==length(log_offset)){
    log_offset=rep(log_offset,length(data_out))
    offset=rep(offset,length(data_out))
  }
  if(length(data_out)!=length(log_offset)){
    stop('Erro: offset/log_offset does not have the same length as data_out')
  }
  
  if(structure$t==1){
    structure$D=array(structure$D,c(structure$n,structure$n,length(data_out)))
    structure$W=array(structure$W,c(structure$n,structure$n,length(data_out)))
    structure$FF=matrix(structure$FF,structure$n,length(data_out))
    structure$t=length(data_out)
  }
  if(length(data_out)!=structure$t){
    stop('Erro: data_out does not have the same length as structure.')
  }
  
  model=kernel(y=data_out,
               m0=structure$m0,
               C0=structure$C0,
               F1=structure$FF,
               G1=structure$G,
               D1=structure$D,
               W1=structure$W,
               pop=log_offset)
  
  model$names=structure$names
  model$data_out=data_out
  
  return(model)
  
}

predict=function(model,t=1,offset=NULL,log_offset=NULL,FF=NULL,D=NULL,W=NULL,plot=TRUE,IC_prob=0.95){
  n=dim(model$mts)[1]
  t_last=dim(model$mts)[2]
  
  if(t>10){
    warning('Warning: Prediction window is big, results will probabily be unreliable.')
  }
  
  #### Consistency check ####
  if(is.null(FF)){
    FF=array(model$F[,t_last],c(n,t))
  }
  if(is.null(D)){
    D=array(model$D[,,t_last],c(n,n,t))
    D[,,2:t_last]=0
  }else{if(all(dim(D)==1)){
    D=ifelse(model$D==0,0,D)
    D[,,2:t_last]=0
  }else{if(length(dim(D))==2 | (length(dim(D))==3 & dim(D)[3]==1)){
    D=array(D,c(n,n,t))
    D[,,2:t_last]=0
  }}
  }
  if(is.null(W)){
    W=array(model$W[,,t_last],c(n,n,t))
    W[,,2:t_last]=0
  }else{if(all(dim(W)==1)){
    W=array(diag(n)*W,c(n,n,t))
    W[,,2:t_last]=0
  }else{if(length(dim(W))==2 | (length(dim(W))==3 & dim(W)[3]==1)){
    W=array(W,c(n,n,t))
    W[,,2:t_last]=0
  }}
  }
  if(!is.null(offset) & !is.null(log_offset)){
    stop('Erro: Cannot set both offset and log_offset. Choose only one.')
  }else{if(!is.null(offset)){
    log_offset=log(offset)
  }else{if(!is.null(log_offset)){
    offset=exp(log_offset)
  }else{
    offset=1
    log_offset=0
  }}}
  if(1==length(log_offset)){
    log_offset=rep(log_offset,t)
    offset=rep(offset,t)
  }
  
  if(dim(FF)[2]!=t){
    stop(paste0('Error: FF should have one column for each time or exactly 1 column, got ',dim(FF)[2],'!=',t,'.'))
  }
  if(dim(FF)[1]!=n){
    stop(paste0('Error: FF should have one line for each latent variable in the model, got ',dim(FF)[1],'!=',n,'.'))
  }
  if(dim(D)[3]!=t){
    stop(paste0('Error: D should have 3º dimention equal to t or 1, got ',dim(D)[3],'.'))
  }
  if(dim(D)[1]!=n | dim(D)[2]!=n){
    stop(paste0('Error: D should have 1º and 2º dimentions equal the number of latent variables in the model, got (',dim(D)[1],',',dim(D)[2],')!=(',n,',',n,').'))
  }
  if(dim(W)[1]!=n | dim(W)[2]!=n){
    stop(paste0('Error: W should have 1º and 2º dimentions equal the number of latent variables in the model, got (',dim(W)[1],',',dim(W)[2],')!=(',n,',',n,').'))
  }
  if(dim(W)[3]!=t){
    stop(paste0('Error: W should have 3º dimention equal to t or 1, got ',dim(W)[3],'.'))
  }
  if(length(offset)!=t){
    stop(paste0('Error: Offset should have length 1 or equal to t, got ',dim(offset)[1],'!=',n,'.'))
  }
  #####
  
  G=model$G
  
  m0=model$mt[,t_last]
  C0=model$Ct[,,t_last]
  
  D <- ifelse(D == 0, 1, D)
  
  # Definindo objetos
  at <- matrix(0, ncol=t, nrow=n)
  mt <- matrix(0, ncol=t, nrow=n)
  ft <- matrix(0, ncol=1, nrow=t)
  qt <- matrix(0, ncol=1, nrow=t)
  Ct <- array(rep(diag(n),t),dim=c(n,n,t))
  Rt <- array(rep(diag(n),t),dim=c(n,n,t))
  a = b= 0
  pred = var.pred = icl.pred = icu.pred = matrix(0, ncol=1, nrow=t)
  
  ## Algoritmo
  
  # Priori
  
  at[,1] <- G%*%m0
  Rt[,,1] <-G%*%C0%*%(t(G))*D[,,1]+W[,,1]
  
  reduc_RFF=Rt[,,1]%*%FF[,1]
  
  # Previsão 1 passo a frente
  ft[1,] <- t(FF[,1])%*%at[,1] + log_offset[1]
  qt[1,] <- t(FF[,1])%*%reduc_RFF
  
  a[1] <- (1/qt[1,])
  b[1] <- (exp(-ft[1,] -0.5*qt[1,])/(qt[1,])) 
  
  # Preditiva em t = 1
  
  pred[1] <- a[1]/ b[1]
  var.pred <- a[1]*(b[1]+1)/(b[1])^2
  icl.pred[1]<-qnbinom((1-IC_prob)/2, a[1], (b[1]/(b[1] +1)))
  icu.pred[1]<-qnbinom(1-(1-IC_prob)/2, a[1], (b[1]/(b[1] +1)))
  
  for(i in c(2:t)){
    # Priori
    
    at[,i] <- G%*%at[,i-1]
    Rt[,,i] <-G%*%Rt[,,i-1]%*%(t(G))*D[,,i]+W[,,i]
    
    reduc_RFF=Rt[,,i]%*%FF[,i]
    
    # Previsão 1 passo a frente
    ft[i,] <- t(FF[,i])%*%at[,i] + log_offset[i]
    qt[i,] <- t(FF[,i])%*%reduc_RFF
    
    a[i] <- (1/qt[i,])
    b[i] <- (exp(-ft[i,] -0.5*qt[i,])/(qt[i,])) 
    
    pred[i] <- a[i]/ b[i]
    var.pred <- a[i]*(b[i]+1)/(b[i])^2
    icl.pred[i]<-qnbinom((1-IC_prob)/2, a[i], (b[i]/(b[i] +1)))
    icu.pred[i]<-qnbinom(1-(1-IC_prob)/2, a[i], (b[i]/(b[i] +1)))
  }
  if(plot){
    fill_list=c('#2596be','#2596be','black')
    names(fill_list)=c(paste0('Prediction I.C. (',(100*IC_prob) %>% round(),'%)'),'Prediction','Observed values')
    color_list=c('#2596be','black')
    names(color_list)=c('Prediction','Observed values')
    
    print(
      ggplotly(
        ggplot()+
          geom_point(aes(x=c(1:t)+t_last,y=pred,color='Prediction',fill='Prediction'))+
          geom_ribbon(aes(x=c(1:t)+t_last,ymin=icl.pred,ymax=icu.pred,fill=paste0('Prediction I.C. (',(100*IC_prob) %>% round(),'%)'),color=paste0('Prediction I.C. (',(100*IC_prob) %>% round(),'%)')),alpha=0.25)+
          geom_point(aes(x=c(1:t_last),y=model$data_out,color='Observed values',fill='Observed values'))+
          scale_fill_manual('',na.value=NA,values=fill_list)+
          scale_color_manual('',na.value=NA,values=color_list)+
          theme_bw()
      )
    )
  }
  
  return(list('pred'=pred,'var.pred'=var.pred,'icl.pred'=icl.pred,'icu.pred'=icu.pred,'at'=at,'Rt'=Rt))
}
eval_past=function(model,smooth=FALSE,t_offset=0){
  if(smooth & t_offset>0){
    t_offset=0
    warning('t_offset is only used if smooth is set to TRUE.')
  }
  n=dim(model$mt)[1]
  t_last=dim(model$mt)[2]
  
  at=array(0,c(n,t_last))
  Rt=array(0,c(n,n,t_last))
  
  ft=array(0,c(t_last))
  qt=array(0,c(t_last))
  a=array(0,c(t_last))
  b=array(0,c(t_last))
  FF=model$F
  G=array(model$G,c(n,n,t_last))
  D=model$D
  W=model$W
  
  pred=c(1:t_last)*0
  
  for(i in c(1:(t_last-t_offset))+t_offset){
    if(smooth){
      at[,i]=model$mts[,i]
      Rt[,,i]=model$Cts[,,i]
    }else{
      at[,i]=model$mt[,i-t_offset]
      Rt[,,i]=model$Ct[,,i-t_offset]
      if(t_offset>0){
        at[,i]=G[,,i-t_offset+1]%*%at[,i]
        Rt[,,i]=G[,,i-t_offset+1]%*%Rt[,,i]%*%t(G[,,i-t_offset+1])*D[,,1]+W[,,1]
        if(t_offset>1){
          multi_G=diag(n)
          for(j in c(2:t_offset)){
            multi_G=G[,,i-t_offset+j]%*%multi_G
          }
          at[,i]=multi_G%*%at[,i]
          Rt[,,i]=multi_G%*%Rt[,,i]%*%t(multi_G)
        }
      }
    }
    ft[i] <- t(FF[,i])%*%at[,i] + model$log_offset[i]
    qt[i] <- t(FF[,i])%*%Rt[,,i]%*%FF[,i]
    
    a[i] <- (1/qt[i])
    b[i] <- (exp(-ft[i] +0.5*qt[i])/(qt[i])) 
    
    pred[i] <- a[i]/ b[i]
    #icl.pred[i]<-qnbinom((1-IC_prob)/2, a[i], (b[i]/(b[i] +1)))
    #icu.pred[i]<-qnbinom(1-(1-IC_prob)/2, a[i], (b[i]/(b[i] +1)))
  }
  
  return(list('pred'=pred,'a'=a,'b'=b))
}
show_fit=function(model,IC_prob=0.95,smooth=TRUE,dinamic=TRUE,t_offset=0){
  n=dim(model$mt)[1]
  t_last=dim(model$mt)[2]
  eval=eval_past(model,smooth,t_offset)
  pred=eval$pred
  a=eval$a
  b=eval$b
  icl.pred<-qnbinom((1-IC_prob)/2, a, (b/(b +1)))
  icu.pred<-qnbinom(1-(1-IC_prob)/2, a, (b/(b +1)))
  
  labels_names=c(paste0('Prediction I.C. (',(100*IC_prob) %>% round(),'%)'),'Prediction','Observed values')
  labels_names=c('I.C. (95%)','Predição','Valores observados')
  
  fill_list=c('#2596be','#2596be','black')
  names(fill_list)=labels_names
  color_list=c('#2596be','black')
  names(color_list)=labels_names[-1]
  
  max_value=calcula_max(model$data_out-min(model$data_out))[[3]]+min(model$data_out)
  min_value=-calcula_max(-(model$data_out-max(model$data_out)))[[3]]+max(model$data_out)
  
  plt=ggplot()+
    geom_point(aes(x=c(1:t_last),y=pred,color=labels_names[2],fill=labels_names[2]))+
    geom_ribbon(aes(x=c(1:t_last),ymin=icl.pred,ymax=icu.pred,fill=labels_names[1],color=labels_names[1]),alpha=0.25)+
    geom_point(aes(x=c(1:t_last),y=model$data_out,color=labels_names[3],fill=labels_names[3]))+
    scale_fill_manual('',na.value=NA,values=fill_list)+
    scale_color_manual('',na.value=NA,values=color_list)+
    scale_y_continuous(name='$y_t$')+
    scale_x_continuous('Time')+
    theme_bw()+
    coord_cartesian(ylim=c(min_value,max_value))
  if(dinamic){
    plt=ggplotly(plt)
  }
  return(list('plot'=plt,'pred'=pred,'icl.pred'=icl.pred,'icu.pred'=icu.pred,'r'=a,'p'=(b/(b +1))))
}

plot_lat_var=function(model,var,smooth=TRUE,cut_off=10,dinamic=TRUE,exp_y=FALSE){
  if(!(var %in% names(model$names))){
    stop(paste0('Error: Invalid selected variable. Got ',var,', expected one of the following:\n',names(model$names)))
  }
  test_data=rnorm(100,0,100)
  
  indice=model$names[[var]][1]
  size=length(indice)
  t=dim(model$mts)[2]
  if(exp_y){
    ref_mt=if(smooth){model$mts}else{model$mt}
    ref_Ct=if(smooth){model$Cts}else{model$Ct}
    lambda=array(0,c(size,t))
    dummy_index=rep(1,size)
    
    for(i in c(1:t)){
      lambda[,i]=exp(ref_mt[indice,i])-1
    }
    m1=(lambda+1)*ref_mt[indice,]-log(lambda+1)*(lambda+1)+lambda
    std_mat=array(0,c(size,size,t))
    for(i in c(1:t)){
      lambda_mat=if(size>1){diag(lambda[,i]+1)}else{lambda[,i]+1}
      std_mat[,,i]=lambda_mat%*%ref_Ct[indice,indice,i]%*%t(lambda_mat)
    }
  }else{
    m1=if(smooth){model$mts[indice,]}else{model$mt[indice,]}
    std_mat=if(smooth){model$Cts[indice,indice,]}else{model$Ct[indice,indice,]}
  }
  m1=m1 %>% matrix(size,t) %>% t
  
  if(size>1){
    std_mat=std_mat %>% apply(3,diag)
  }
  std_mat=std_mat %>% sqrt %>% matrix(size,t) %>% t
  
  lim_i=m1-2*std_mat
  lim_s=m1+2*std_mat
  
  m1=as.data.frame(m1)
  lim_i=as.data.frame(lim_i)
  lim_s=as.data.frame(lim_s)
  std_mat=as.data.frame(std_mat)
  
  names(m1)=paste0('Variable ',c(1:dim(m1)[2]))
  names(lim_i)=paste0('Variable ',c(1:dim(lim_i)[2]))
  names(lim_s)=paste0('Variable ',c(1:dim(lim_s)[2]))
  names(std_mat)=paste0('Variable ',c(1:dim(std_mat)[2]))
  
  max_value=calcula_max(m1-min(m1))[[3]]+min(m1)
  min_value=-calcula_max(-(m1-max(m1)))[[3]]+max(m1)
  
  m1$time=c(1:dim(m1)[1])
  lim_i$time=c(1:dim(lim_i)[1])
  lim_s$time=c(1:dim(lim_s)[1])
  std_mat$time=c(1:dim(std_mat)[1])
  
  m1=m1 %>% pivot_longer(1:size) %>% rename(mean=value)
  lim_i=lim_i %>% pivot_longer(1:size) %>% rename(lim_i=value)
  lim_s=lim_s %>% pivot_longer(1:size) %>% rename(lim_s=value)
  std_mat=std_mat %>% pivot_longer(1:size) %>% rename(std=value)

  plot_data = m1%>%
      inner_join(std_mat,by=c('time','name')) %>%
      inner_join(lim_i,by=c('time','name')) %>%
      inner_join(lim_s,by=c('time','name'))
  
  n_var=length(unique(plot_data$name))
  color_list=rainbow(n_var,s=0.5)
  names(color_list)=paste(unique(plot_data$name),'point estimate')
  
  fill_list=rainbow(n_var,s=0.5)
  names(fill_list)=paste0(unique(plot_data$name),' I.C. (',0.95*100 %>% round(),'%)')
  
  plt=ggplot(plot_data[plot_data$time>cut_off,])+
    geom_hline(yintercept=0,linetype='dashed')+
    scale_x_continuous('Time')+
    scale_color_manual('',values=color_list,na.value=NA)+
    scale_fill_manual('',values=fill_list,na.value=NA)+
    labs(title=paste0(var,' (',ifelse(smooth,'smoothed','only filtered'),')'))+
    scale_y_continuous('Parameter value')+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90))+
    geom_ribbon(aes(x=time,ymin=lim_i,ymax=lim_s,fill=paste0(name,' I.C. (',0.95*100 %>% round(),'%)'),color=paste0(name,' I.C. (',0.95*100 %>% round(),'%)')),alpha=0.25)+
    geom_line(aes(x=time,y=mean,color=paste(name,'point estimate'),fill=paste(name,'point estimate')))+
    coord_cartesian(ylim=c(min_value,max_value))
  if(dinamic){
    plt=ggplotly(plt)
  }
  list('plot'=plt,'table'=plot_data)
}


# ### Exemplos ####
# 
# ### Importando dados ####
# 
# 
# dados=read.csv('varicela\\data\\varicela internacoes.csv')[,c(1,7:162)]
# dados[1:2,1]='00 a 04 anos'
# dados[5:8,1]='15 a 49 anos'
# dados[9:12,1]='50 anos e mais'
# dados=aggregate(.~FaixaEtaria,dados,sum)[,-1]
# 
# pre_exp=read.csv2('varicela\\data\\populacao 2000-2020.csv')[-12,c(1,10:22)]
# pre_exp[4:7,1]='15 a 49 anos'
# pre_exp[8:11,1]='50 anos e mais'
# pre_exp=aggregate(.~FaixaEtaria,pre_exp,sum)[,-1]
# 
# dummy=matrix(0,dim(pre_exp)[1],0)
# nomes=c()
# for(ano in c(2008:2020)){
#   for(mes in c(1:12)){
#     nomes=c(nomes,paste0('X',ano,'.',mes))
#     dummy=cbind(dummy,pre_exp[,ano-2007])
#   }
# }
# pre_exp=dummy
# 
# idade_indice=1
# inter=as.numeric(dados[idade_indice,])
# exp=as.data.frame(pre_exp)[idade_indice,]
# names(exp)=nomes
# exp=as.numeric(exp)
# 
# N <- dim(dados)[2]
# indice_inter=69
# data_lab=names(dados)
# 
# indic_inter=c(rep(0,indice_inter-1),rep(1,N-indice_inter+1))
# covid=c(rep(0,146),rep(1,N-146))
# var_covid=ifelse(is.na(covid),0,covid)
# 
# ### Criando modelo ####
# 
# nivel_bloc=gera_bloco_poly(2,D=1/0.9,name='Nível')
# inter_bloc=gera_bloco_poly(1,value=indic_inter,D=1/1,name='Vacina')
# covid_bloc=gera_bloco_poly(1,value=var_covid,D=1/1,name='Covid')
# sazo_bloc=gera_bloco_sazo(12,D=1/0.98,name='Sazonalidade')
# 
# estrutura=concat_bloco(nivel_bloc,inter_bloc,covid_bloc,sazo_bloc)
# 
# teste <- ajusta_modelo(estrutura,data_out=inter,offset=exp)
# 
# plot_lat_var(teste,'Vacina',smooth = T,exp_y=TRUE)
# 
# predicao=predict(teste,t=10,offset=exp[N],log_offset=NULL,FF=NULL,D=NULL,W=NULL,plot=T)

# # exp=rpois(100,10)+2
# # y=rpois(100,10*exp*c(rep(1,50),rep(2,50)))
# # D_value=log(y)/rnorm(100,2)
# # indic=c(rep(0,50),rep(1,50))
# # 
# # 
# # A=gera_bloco_poly(1,value=log(exp),m0=1,C0=0)
# # B=gera_bloco_sazo(12,D=1/0.95)
# # C=gera_bloco_poly(3,D=1/0.9)
# # D=gera_bloco_poly(1,value=D_value,D=1/0.9)
# # E=gera_bloco_poly(1,value=indic,m0=0,C0=0,D=1/0.9)
# # E$W[,,51]=1
# # 
# # estrutura=concat_bloco(A,E,D,B,C)
# # 
# # teste1=ajusta_modelo(y,structure=estrutura)
# # 
# # 
# # B=gera_bloco_sazo(12,D=1/0.95)
# # C=gera_bloco_poly(3,D=1/0.9)
# # D=gera_bloco_poly(1,value=D_value,D=1/0.9)
# # E=gera_bloco_poly(1,value=indic,m0=0,C0=0,D=1/0.9)
# # E$W[,,51]=1
# # 
# # estrutura=concat_bloco(E,D,B,C)
# # 
# # teste2=ajusta_modelo(y,structure=estrutura,offset=exp)
# # 
# # plot(teste1$pred)
# # points(teste2$pred)
# # 
# # print(mean(abs(teste1$pred-teste2$pred)/teste2$pred))
# # 
# # plot(teste1$mts[2,])
# 
# 

# # Exemplo de uso --------
# y <- c(131.7, 322.6, 285.6, 105.7,
#        80.4, 285.1, 347.8, 68.9,
#        203.3, 375.9, 415.9, 65.8,
#        177.0, 438.3, 463.2, 136.0,
#        192.2, 442.8, 509.6, 201.2,
#        196.0, 478.6, 688.6, 259.8,
#        352.5, 508.1, 701.5, 325.6,
#        305.9, 422.2, 771.0, 329.3,
#        384.0, 472.0, 852.0)
# y <- trunc(y/10)
# 
# plot(y, type="l", xlab = "")
# n <- 6
# N <- length(y)
# w  <- 2*pi/4
# 
# FF <- matrix(c(rep(1,N),rep(0,N),rep(1,N),rep(0,N),rep(1,N),rep(0,N)), ncol=N,nrow=n, byrow=T)
# 
# 
# m0 <- matrix(0,ncol=1, nrow=n)
# C0 <- diag(n)
# 
# G1 <- matrix(c(1,1,0,1), nrow=2, byrow =T)
# G2 <- matrix(c(cos(w),-sin(w),sin(w),cos(w)),2,2)
# G3 <- matrix(c(cos(2*w),-sin(2*w),sin(2*w),cos(2*w)),2,2)
# G  <- as.matrix(bdiag(G1, G2, G3))
# 
# # Matriz de desconto
# 
# D=array(0,c(n,n,N))
# 
# # Tendência
# delta.1 <- 0.9
# D[1:2,1:2,] <- 1/delta.1
# 
# # Sazonalidade
# delta.2 <- 0.9
# D[3:4,3:4,] <- 1/delta.2
# 
# # Sazonalidade
# delta.3 <- 0.9
# D[5:6,5:6,] <- 1/delta.3
# 
# 
# teste <- poisson_gi_exp(y,m0, C0, F1 = FF ,G1 = G, D1 = D, pop=log(1))
# aux <- rep(NA, 9)
# plot(c(aux,y[10:N]), pch=20, ylim=c(0,140), xlab = "", ylab = "")
# lines(c(aux,teste$pred[10:N]), col="#d8790d", lwd = 2)
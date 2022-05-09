poisson_gi_exp <- function(y,m0 = 0, C0 = 1, FF,G,D,W, pop, IC_prob=0.95){

  # Definindo quantidades

  r <- 1
  T <- length(y)
  n <- dim(FF)[1]

  # D1.aux <- matrix(rep(D1,n1^2), ncol = n1 )
  # D2.aux <- matrix(rep(D2,n2^2), ncol = n2 )
  # D.aux <- as.matrix(bdiag(D1.aux, D2.aux))
  D.aux <- D1
  D <- ifelse(D.aux == 0, 1, D.aux)

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

  reduc_RFF=Rt[,,1]%*%FF[,1,1]

  # Previsão 1 passo a frente
  ft[1,] <- t(FF[,1,1])%*%at[,1] + pop[1]
  qt[1,] <- t(FF[,1,1])%*%reduc_RFF

  # Compatibilizando prioris

  a[1] <- (1/qt[1,])
  b[1] <- (exp(-ft[1,] -0.5*qt[1,])/(qt[1,]))

  # Posteriori

  a.post[1] <- a[1] + y[1]
  b.post[1] <- b[1] + 1

  gt[1] <- log(a.post[1]/b.post[1]) + 1/(2*a.post[1])
  pt[1] <- (2*a.post[1]-1)/(2*a.post[1]^2)

  mt[,1] <- at[,1]+reduc_RFF*(gt[1]-ft[1,])*(1/(qt[1,]))
  Ct[,,1] <- Rt[,,1] - (reduc_RFF%*%t(reduc_RFF))*(1 - pt[1]/qt[1,])*(1/qt[1,])

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

    reduc_RFF=Rt[,,t]%*%FF[,1,t]

    # Previsão 1 passo a frente

    ft[t,] <- t(FF[,1,t])%*%at[,t] + pop[t]
    qt[t,] <- t(FF[,1,t])%*%reduc_RFF

    # Compatibilizando prioris

    a[t] <- (1/qt[t,])
    b[t] <- (exp(-ft[t,] )/(qt[t,]))

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

multinom_gi <- function(y,m0=0, C0=1, FF,G,D,W, pop, IC_prob=0.95){
  
  # Função a ser otimizada
  # 
  # otim1 <- function(x, f1,f2,q1,q2,q12, Omega){
  #   alpha1 <- x[1]
  #   alpha2 <- x[2]
  #   alpha3 <- x[3]
  #   
  #   eqs = rbind(
  #     f1 - digamma(alpha1) + digamma(alpha1 + alpha2 + alpha3),
  #     f2 - digamma(alpha2) + digamma(alpha1 + alpha2 + alpha3),
  #     q1 - trigamma(alpha1)  - trigamma(alpha1 + alpha2 + alpha3),
  #     q2 - trigamma(alpha2) - trigamma(alpha1 + alpha2 + alpha3),
  #     q12 - trigamma(alpha1 + alpha2 + alpha3))
  #   return(t(eqs)%*%Omega%*%eqs)
  # }
  # 
  
  model_tau0_e_tau1 <- function(x, parms){
    sub_last=digamma(x[length(x)] - sum(x[1:(length(x)-1)]))
    digamma_vec=digamma(x)
    
    f_all=f-digamma_vec[-length(x)]+sub_last
    last_guy=media.log+digamma_vec[length(x)]-sub_last
    
    f_all=c(f_all,last_guy)
    
    return(f_all)
  }
  

  f <- function(m) {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    m
  }
  
  otim2 <- function(x, f,q, Omega){
    sub_last=digamma(x[length(x)] - sum(x[1:(length(x)-1)]))
    digamma_vec=digamma(x)
    
    f_all=f-digamma_vec[-length(x)]+sub_last
    last_guy=media.log+digamma_vec[length(x)]-sub_last
    
    f_all=c(f_all,last_guy)
    
    calc_helper=1 + sum(exp(f_all))
    
    H=exp(f_all)%*%t(exp(f_all))/(calc_helper**2)
    diag(H)=-(exp(f_all)*calc_helper-(exp(f_all)**2))/(calc_helper**2)
    
    media.log = 
      -log(calc_helper) + (H%>%q) %>% diag %>% sum
    
    return(t(eqs)%*%Omega%*%eqs)
  }
  
  
  # Definindo quantidades
  T <- nrow(y)
  n <- dim(FF)[1]
  
  D.aux <- D
  D <- ifelse(D.aux == 0, 1, D.aux)
  
  r = ncol(y)-1
  m0 <- matrix(m0,n,1)
  C0 <- C0
  mt <- matrix(0,nrow=n,ncol=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  Rt <- array(rep(0,T),dim=c(n,n,T))
  Pt <- array(rep(0,T),dim=c(n,n,T))
  Wt <- array(rep(0,T),dim=c(n,n,T))
  ft <- matrix(0,nrow=T,ncol=r)
  at <- matrix(0,nrow=n,ncol=T)
  Qt <-  array(rep(0,T),dim=c(r,r,T))
  et <- matrix(0,nrow=T,ncol=r)
  At <- array(rep(0,T),dim=c(n,r,T))
  dt <- matrix(0,nrow=T,ncol=r)
  nt <- matrix(0,nrow=T,ncol=r)
  
  f_star <- matrix(0,nrow=T,ncol=r)
  Q_star <- array(0,c(r,r,T))
  
  mt <- matrix(0,nrow=n,ncol=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  
  tau <- matrix(NA,nrow=r+1,ncol=T)
  
  alpha <- matrix(NA,nrow=r+1,ncol=T)
  
  alpha_star <- matrix(NA,nrow=r+1,ncol=T)
  
  tau_star = matrix(NA,nrow=r+1,ncol=T)
  
  # #Auxiliar do desconto
  #if(is.null(dim(G)) == FALSE){
  #  matrixaux1 <- G == 0
  #  if(G[1,1] == 1 & G[1,2] == 1 & G[2,2] == 1 ) {matrixaux1[2] <- 0}
  #  tira1 <- which(matrixaux1 == 1)
  #  mantem1 <- which(matrixaux1 == 0)}
  # 
  # if(is.null(dim(G2)) == FALSE){
  #   matrixaux2 <- G2 == 0 
  #   if(G2[1,1] == 1 & G2[1,2] == 1 & G2[2,2] == 1 ) {matrixaux2[2] <- 0}
  #   tira2 <- which(matrixaux2 == 1)
  #   mantem2 <- which(matrixaux2 == 0)}
  # 
  
  D=ifelse(D==0,1,D)
  
  # Priori em t = 1
  
  at[,1]          = (G%*%m0)[,1]
  Pt            <-  G%*%C0%*%(t(G))
  Rt[,,1]       <- as.matrix(D[,,1]*Pt)+W[,,1]
  
  # Previsao em t = 1
  ft[1,]        <- (t(FF[,,1])%*%at[,1])[,1]
  Qt[,,1]       <- as.matrix(t(FF[,,1])%*%Rt[,,1]%*%FF[,,1])
  
  # minimizando...
  
  f = ft[1,]
  q = Qt[,,1]
  
  calc_helper=1 + sum(exp(f))
  
  H=exp(f)%*%t(exp(f))/(calc_helper**2)
  diag(H)=-(exp(f)*calc_helper-(exp(f)**2))/(calc_helper**2)
  
  media.log = 
    -log(calc_helper) + 0.5*(H%*%q) %>% diag %>% sum
  
  parms = c(f, media.log)
  
  ss1 <- multiroot(f = model_tau0_e_tau1 , start = c(rep(0.01,r),0.01*(r+1)), parms = parms)
  
  tau[,1] <- as.vector(ss1$root)
  
  alpha[,1]      <- tau[,1]  
  alpha[r+1,1]      <- tau[r+1,1] - sum(tau[-r-1,1]) 
  
  # alpha_star[-r-1,1] <-    alpha[-r-1,1]  +  y[1,-r-1]
  # alpha_star[r+1,1] <-    alpha[r+1,1]  +  N[1]-sum(y[1,-r-1])
  
  alpha_star[,1] <-    alpha[,1]  +  y[1,]
  
  tau_star[,1]  <-  alpha_star[,1]
  tau_star[r+1,1]  <-  sum(alpha_star[,1])
  
  # Posteriori
  f_star[1,] <-   digamma(alpha_star[-r-1,1]) -  digamma(alpha_star[r+1,1])
  Q_star[,,1] <-   trigamma(alpha_star[r+1,1])
  diag(Q_star[,,1]) <- trigamma(alpha_star[-r-1,1])+trigamma(alpha_star[r+1,1])
  
  At[,,1] <- as.matrix(Rt[,,1]%*%FF[,,1]%*%ginv(Qt[,,1]))
  mt[,1] <- at[,1] + At[,,1]%*%(f_star[1,] -ft[1,])
  Ct[,,1] <- Rt[,,1] +  At[,,1]%*%(Q_star[,,1] - Qt[,,1])%*%t(At[,,1])
  
  for(t in 2:T){
    at[,t] = as.matrix(G%*%mt[,t-1])
    Pt <- G%*%Ct[,,t-1]%*%(t(G))
    Rt[,,t] <- as.matrix(D[,,t]*Pt)+W[,,t]
    
    # Previsao em t = 1
    ft[t,]        <- (t(FF[,,t])%*%at[,t])[,1]
    Qt[,,t]       <- as.matrix(t(FF[,,t])%*%Rt[,,t]%*%FF[,,t])
    
    # minimizando...
    
    f = ft[t,]
    q = Qt[,,t]
    calc_helper=1 + sum(exp(f))
    
    H=exp(f)%*%t(exp(f))/(calc_helper**2)
    diag(H)=-(exp(f)*calc_helper-(exp(f)**2))/(calc_helper**2)
    
    media.log = 
      -log(calc_helper) + 0.5*(H%*%q) %>% diag %>% sum
    
    parms = c(f, media.log)
    
    ss1 <- multiroot(f = model_tau0_e_tau1 , start = c(rep(0.01,r),0.01*(r+1)), parms = parms)
    
    tau[,t] <- as.vector(ss1$root)
    
    alpha[,t]      <- tau[,t]  
    alpha[r+1,t]      <- tau[r+1,t] - sum(tau[-r-1,t]) 
    
    # alpha_star[-r-1,t] <-    alpha[-r-1,t]  +  y[t,-r-1]
    # alpha_star[r+1,t] <-    alpha[r+1,t]  +  N[t]-sum(y[t,-r-1])
    
    alpha_star[,t] <-    alpha[,t]  +  y[t,]
    
    tau_star[,t]  <-  alpha_star[,t]
    tau_star[r+1,t]  <-  sum(alpha_star[,t])
    
    # Posteriori
    f_star[t,] <-   digamma(alpha_star[-r-1,t]) -  digamma(alpha_star[r+1,t])
    Q_star[,,t] <-   trigamma(alpha_star[r+1,t])
    diag(Q_star[,,t]) <- trigamma(alpha_star[-r-1,t])+trigamma(alpha_star[r+1,t])
    
    At[,,t] <- as.matrix(Rt[,,t]%*%FF[,,t]%*%ginv(Qt[,,t]))
    mt[,t] <- at[,t] + At[,,t]%*%(f_star[t,] -ft[t,])
    Ct[,,t] <- Rt[,,t] +  At[,,t]%*%(Q_star[,,t] - Qt[,,t])%*%t(At[,,t])
    
  }
  
  
  mts <- matrix(0, ncol=T, nrow=n)
  Cts <- array(rep(diag(n),T),dim=c(n,n,T))
  
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
                 ft, Qt,
                 f_star, Q_star,
                 alpha,alpha_star,
                 tau,tau_star,
                 FF, G, D,W,
                 mts, Cts ,
                 exp(pop),pop)
  names(result) <- c("mt",  "Ct",
                     "ft", "Qt",
                     'f_star', 'Q_star',
                     'alpha','alpha_star',
                     'tau','tau_star',
                     "FF", "G", "D","W",
                     "mts", "Cts",
                     'offset','log_offset')
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
  b[1] <- (exp(-ft[1,] -0.5*qt[1,])/(qt[1,]))

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
    b[t] <- (exp(-ft[t,] )/(qt[t,]))

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

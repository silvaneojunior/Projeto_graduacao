normal.analise <- function(y,m01, C01,m02,C02, F1,F2,G1,G2,D1,D2){
  
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
  
  model_tau <- function(x, parms){
    
    tau0=x[1]
    tau1=x[2]
    tau2=x[3]
    tau3=x[4]
    
    Eq1=(q1 + f1^2)*exp(f2 + q2/2)
    
    par1=(tau2^2 - 4*tau1*tau3)
    n0=2*tau0+ 1
    
    P11=(n0*(tau2^2))/(2*tau1*par1)
    P12=-1/(2*tau1)
    P1=P11+P12
    
    Eq2=f1*exp(f2 + q2/2)
    P2=-tau2*n0/par1
    
    Eq3=exp(f2 + 0.5*q2)
    P3=2*tau1*n0/par1
    
    Eq4=f2
    P4=digamma(n0/2) - log(par1/(4*tau1))
    
    output=c(
    F1 =  Eq1-P1 ,
    F2 =  Eq2-P2,
    F3 =  Eq3-P3,
    F4 =  Eq4-P4)
    # print('c0')
    # print(-2*x[2])
    # print('mu0')
    # print(-x[3]/(2*x[2]))
    # print('d0/2')
    # print((x[3]**2)/(4*x[2])-x[4])
    # print('n0/2')
    # print(x[1]+0.5)
    # print('output')
    # print(output)
    # print('F1')
    # print((q1 + f1^2)*exp(f2 + q2/2))
    # print('F2')
    # print(f1*exp(f2 + q2/2))
    # print('F3')
    # print(exp(f1 + 0.5*q2))
    # print('F4')
    # print(f2)
    return(output)
  }
  
  root_init=c(0.3,-0.3,0.0,-0.3)
  x=root_init
  c(
    (x[3]^2) - 4*x[2]*x[4]
    )
  
  f <- function(m) {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    m
  }
  
  # otim2 <- function(x, f1,f2,q1,q2,q12, Omega){
  #   tau1 <- x[1]
  #   tau2 <- x[2]
  #   tau0 <- x[3]
  #   media.log = log(1/(1+exp(f1) + exp(f2))) - 
  #     (q1/2)*(exp(f1)*exp(f2) + exp(f1))/(exp(f1)+exp(f2)+ 1)^2 - 
  #     (q2/2)*(exp(f2)*exp(f1) + exp(f2))/(exp(f1)+exp(f2)+1)^2 + 
  #     q12*(exp(f1 + f2)/(1+exp(f1)+exp(f2))) 
  #   
  #   eqs = rbind(
  #     f1  - digamma(tau1) + digamma(tau0 - tau1 - tau2),
  #     f2  - digamma(tau2) + digamma(tau0 - tau1 - tau2),
  #     digamma(tau0) - digamma(tau0 - tau1 - tau2) + media.log)
  #   
  #   return(t(eqs)%*%Omega%*%eqs)
  # }
  # 
  
  
  # Definindo quantidades
  n1 <- nrow(F1)
  n2 <- nrow(F2)
  T <- length(y)
  n <- n1 + n2
  F <- as.matrix(bdiag(F1[,1], F2[,1]))
  G <- as.matrix(bdiag(G1,G2))
  
  D.aux <- as.matrix(bdiag(D1, D2))
  D <- ifelse(D.aux == 0, 1, D.aux)
  
  m0 <- as.matrix(c(m01, m02), nrow = 1)
  C0 <- as.matrix(bdiag(C01, C02))
  
  r = 2
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
  
  f1star <- matrix(0,nrow=T)
  f2star <- matrix(0,nrow=T)
  Q1star <- matrix(0,nrow=T)
  Q2star <- matrix(0,nrow=T)
  Q12star <- matrix(0,nrow=T)
  fstar <-  matrix(0,nrow=T,ncol=r)
  Qstar <-  matrix(0,nrow=T,ncol=r)
  
  m1 <- matrix(0,nrow=n1,ncol=T)
  C1 <- array(rep(diag(n1),T),dim=c(n1,n1,T))
  
  tau0 <- rep(NA,l=T)
  tau1 <- rep(NA,l=T)
  tau2 <- rep(NA,l=T)
  tau3 <- rep(NA,l=T)
  
  tau0_star = rep(NA, l = T)
  tau1_star = rep(NA, l = T)
  tau2_star = rep(NA, l = T)
  tau3_star = rep(NA, l = T)
  
  # #Auxiliar do desconto
  if(is.null(dim(G)) == FALSE){
    matrixaux1 <- G == 0
    if(G[1,1] == 1 & G[1,2] == 1 & G[2,2] == 1 ) {matrixaux1[2] <- 0}
    tira1 <- which(matrixaux1 == 1)
    mantem1 <- which(matrixaux1 == 0)}
  # 
  # if(is.null(dim(G2)) == FALSE){
  #   matrixaux2 <- G2 == 0 
  #   if(G2[1,1] == 1 & G2[1,2] == 1 & G2[2,2] == 1 ) {matrixaux2[2] <- 0}
  #   tira2 <- which(matrixaux2 == 1)
  #   mantem2 <- which(matrixaux2 == 0)}
  # 
  
  # Priori em t = 1
  
  D=ifelse(D==0,1,D)
  
  at[,1]          = G%*%m0
  Pt            <-    G%*%C0%*%(t(G))
  Rt[,,1]       <- D*Pt
  
  
  # Previsao em t = 1
  ft[1,]        <- t(F)%*%at[,1]
  Qt[,,1]       <- t(F)%*%Rt[,,1]%*%F 
  
  # minimizando...
  
  f1=ft[1,1]
  f2=ft[1,2]
  q1 = Qt[1,1,1]
  q2 = Qt[2,2,1]
  q12 = Qt[1,2,1]
  
  parms = c(f1,f2, q1, q2)
  
  print(parms)
  
  ss1 <- multiroot(f = model_tau , start = root_init, parms = parms)
  print(ss1$root)
  print(model_tau(ss1$root,parms))
  
  tau0[1] <- ss1$root[1]
  tau1[1] <- ss1$root[2]
  tau2[1] <- ss1$root[3]
  tau3[1] <- ss1$root[4]
  
  
  tau0_star[1]  <-  tau0[1] + 1/2
  tau1_star[1]  <-  tau1[1] - 1/2
  tau2_star[1]  <-  tau2[1] + y[1]
  tau3_star[1]  <-  tau3[1] - y[1]^2
  
  aux_1=digamma(tau0_star[1] + 0.5) -log(((tau2_star[1]^2)/(4*tau1_star[1])) - tau3_star[1])
  
  # Posteriori
  f1star[1,] <-   -2*tau2_star[1]/(2*tau1_star[1])
  f2star[1,] <-   aux_1
  Q1star[1,] <-   (tau2_star[1]^2)/(4*tau1_star[1]^2) - 8*(tau1_star[1]^2)*(tau0_star[1] + 1/2)/(tau2_star[1]^2 - 4*tau1_star[1]*(tau3_star[1]))
  Q2star[1,] <-   trigamma(tau0_star[1] + 0.5) + aux_1^2
  Q12star[1,] <-   (-tau2_star[1]/(2*tau1_star[1]))*aux_1
  
  fstar <- c(f1star[1,],   f2star[1,])
  Qstar <- matrix(c( Q1star[1,],  Q12star[1,],  Q12star[1,],  Q2star[1,]), byrow =F, ncol = 2)
  
  At[,,1] <- Rt[,,1]%*%F%*%ginv(Qt[,,1])
  mt[,1] <- at[,1] + At[,,1]%*%(fstar -ft[1,])
  Ct[,,1] <- Rt[,,1] +  At[,,1]%*%(Qstar - Qt[,,1])%*%t(At[,,1])
  
  for(t in 2:T){
    print('time')
    print(t)
    
    at[,t] = G%*%mt[,t-1] 
    Pt <- G%*%Ct[,,t-1]%*%(t(G))
    Rt[,,t] <- D*Pt
    
    # Previsao em t = 1
    ft[t,] <-  t(F)%*%at[,t]
    Qt[,,t] <- t(F)%*%Rt[,,t]%*%F
    
    # minimizando...
    
    f1=ft[t,1]
    f2=ft[t,2]
    q1 = Qt[1,1,t]
    q2 = Qt[2,2,t]
    q12 = Qt[1,2,t]
    
    parms = c(f1,f2, q1, q2)
    
    ss1 <- multiroot(f = model_tau , start = root_init, parms = parms)
    
    print(model_tau(ss1$root,parms))
    
    tau0[t] <- ss1$root[1]
    tau1[t] <- ss1$root[2]
    tau2[t] <- ss1$root[3]
    tau3[t] <- ss1$root[4]
    
    tau0_star[t]  <-  tau0[t] + 1/2
    tau1_star[t]  <-  tau1[t] - 1/2
    tau2_star[t]  <-  tau2[t] + y[t]
    tau3_star[t]  <-  tau3[t] - y[t]^2
    
    aux_1=digamma(tau0_star[t] + 0.5) -log(((tau2_star[t]^2)/(4*tau1_star[t])) - tau3_star[t])
    
    # Posteriori
    f1star[t,] <-   -2*tau2_star[t]/(2*tau1_star[t])
    f2star[t,] <-   aux_1
    Q1star[t,] <-   (tau2_star[t]^2)/(4*tau1_star[t]^2) - 8*(tau1_star[t]^2)*(tau0_star[t] + 1/2)/(tau2_star[t]^2 - 4*tau1_star[t]*(tau3_star[t]))
    Q2star[t,] <-   trigamma(tau0_star[t] + 0.5) + aux_1^2
    Q12star[t,] <-   (-tau2_star[t]/(2*tau1_star[t]))*aux_1
    
    
    fstar <- c(f1star[t,],   f2star[t,])
    Qstar <- matrix(c( Q1star[t,],  Q12star[t,],  Q12star[t,],  Q2star[t,]), byrow =F, ncol = 2)
    
    At[,,t] <- Rt[,,t]%*%F%*%ginv(Qt[,,t])
    mt[,t] <- at[,t] + At[,,t]%*%(fstar -ft[t,])
    Ct[,,t] <- Rt[,,t] +  At[,,t]%*%(Qstar - Qt[,,t])%*%t(At[,,t])
    
  }
  
  
  mts <- matrix(0, ncol=T, nrow=n)
  Cts <- array(rep(diag(n),T),dim=c(n,n,T))
  
  mts[,T] <- mt[, T]
  Cts[,,T] <- Ct[,,T]
  for(t in (T-1):1){
    mts[,t] <- mt[,t] + Ct[,,t]%*%t(G)%*%solve(Rt[,,t+1])%*%(mts[,t+1] - at[,t+1])
    Cts[,,t] <- Ct[,,t] - Ct[,,t]%*%t(G)%*%solve(Rt[,,t+1])%*%(Rt[,,t+1] - Cts[,,t+1])%*%t(Ct[,,t]%*%t(G)%*%solve(Rt[,,t+1]))
  }
  
  result <- list(mt,Ct,f1star, Q1star,ft, Qt,
                 mts, Cts, tau0_star, tau1_star,tau2_star ,tau3_star, F, G, D, tau0, tau1, tau2, tau3)
  names(result) <- c("mt",  "Ct",  "f1star", "Q1star","ft", "qt",
                     "mts","Cts", "tau0_star", "tau1_star",  "tau_2_star", "tau3_star",
                     "F", "G", "D", "tau0", "tau1", "tau2", "tau4")
  return(result)
}

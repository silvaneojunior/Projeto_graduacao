library(MASS)
library(ggplot2)
library(plotly)

W=matrix(10,1,1)
V_pre=matrix(c(1,0,0,1),2,2)
C0=diag(c(1,1))

theta=matrix(1,1000,2)
Y=matrix(0,1000,1)

F=matrix(c(1,0),1,2)
G=matrix(c(1,0,1,1),2,2)
V=G%*%V_pre%*%t(G)

for(i in 1:999){
  theta[i+1,]=mvrnorm(n=1,mu=G%*%theta[i,],Sigma=V)
  Y[i+1,]=mvrnorm(n=1,mu=F%*%theta[i+1,],Sigma=W)
}

ggplotly(
  ggplot()+
  geom_point(aes(x=1:1000,y=Y,color='obs.'))+
  geom_line(aes(x=1:1000,y=theta[,1],color='theta_1'))+
  geom_line(aes(x=1:1000,y=theta[,2],color='theta_2'))+
  theme_bw()
)


#-----------------------------------------------------------

init=1
end=100
theta_true=theta[init:end,]
Y=Y[init:end,]

#-----------------------------------------------------------

n=50000
theta=array(0,dim=c(n,end-init+1,2))

W_inv=solve(W)
V_inv=solve(V)
Var_inv=(t(F)%*%solve(W)%*%F+solve(V))
Var=solve(Var_inv)

pre_calc1=V_inv%*%G
pre_calc2=t(F)%*%W_inv
pre_calc3=solve(pre_calc2%*%F+V_inv)
pre_calc4=pre_calc2%*%t(Y)

theta[,1,]=mvrnorm(n=n,mu=pre_calc4[,1]%*%Var,Sigma=solve(t(F)%*%W_inv%*%F+C0))
for(i in 2:(end-init+1)){
  theta[,i,]=mvrnorm(n=n,mu=c(0,0),Sigma=pre_calc3)
  theta[,i,]=theta[,i,]+t(pre_calc1%*%t(theta[,i-1,])+pre_calc4[,i])%*%Var
}

sample_T=theta[c(1:5000)*10,,]
mean_T=colMeans(sample_T)
pred=c()
q975_1=c()
q025_1=c()
q975_2=c()
q025_2=c()
pred975=c()
pred025=c()

for(i in 1:100){
  q975_1=c(q975_1,quantile(sample_T[,i,1],0.975))
  q025_1=c(q025_1,quantile(sample_T[,i,1],0.025))
  q975_2=c(q975_2,quantile(sample_T[,i,2],0.975))
  q025_2=c(q025_2,quantile(sample_T[,i,2],0.025))
  next_pred=F%*%G%*%t(sample_T[,i,])+rnorm(n/10,0,W)
  pred975=c(pred975,quantile(next_pred,0.975))
  pred025=c(pred025,quantile(next_pred,0.025))
  pred=c(pred,mean(next_pred))
}

ggplotly(
  ggplot()+
    geom_point(aes(x=1:(end-init+1),y=Y,color='obs.',fill='obs.'))+
    geom_line(aes(x=1:(end-init+1)+1,y=pred,color='obs.',fill='obs.',linetype='prediction'))+
    geom_ribbon(aes(x=1:(end-init+1)+1,ymin=pred025,ymax=pred975,color='obs.',fill='obs.',linetype='prediction'),alpha=0.25)+
    geom_line(aes(x=1:(end-init+1),y=mean_T[,1],color='theta_1',fill='theta_1',linetype='estimated'))+
    geom_ribbon(aes(x=1:(end-init+1),ymin=q025_1,ymax=q975_1,color='theta_1',fill='theta_1',linetype='estimated'),alpha=0.25)+
    geom_line(aes(x=1:(end-init+1),y=colMeans(theta[c(1:5000)*10,,])[,2],color='theta_2',fill='theta_2',linetype='estimated'))+
    geom_ribbon(aes(x=1:(end-init+1),ymin=q025_2,ymax=q975_2,color='theta_2',fill='theta_2',linetype='estimated'),alpha=0.25)+
    geom_line(aes(x=1:(end-init+1),y=theta_true[,1],color='theta_1',fill='theta_1',linetype='true'))+
    geom_line(aes(x=1:(end-init+1),y=theta_true[,2],color='theta_2',fill='theta_2',linetype='true'))+
    theme_bw()
)

plot((Y[-1]-pred[-1000])/sqrt(W+F%*%V%*%t(F)))
print(
  var(
    (Y[-1]-pred[-1000])/sqrt(W+F%*%V%*%t(F))
    )
  )
acf((Y[-1]-pred[-1000])/sqrt(W+F%*%V%*%t(F)))

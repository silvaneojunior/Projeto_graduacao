library(ggplot2)

data=read.csv2(paste0('data/Ano/','varicela internacoes',' 2007-2021.csv'))[1:12,c(-16,-17)]

y=as.numeric(as.matrix(data)[1,-c(1,2)])
n=5000
samples=matrix(0,n,13)
beta_0=1


  last_guy=rgamma(n,shape=beta_0*0+y[1],rate=1)
  samples[,1]=last_guy
  for(j in c(2:12)){
    momentum=last_guy
    last_guy=rgamma(n,shape=beta_0*momentum+y[j],rate=beta_0+1)
    samples[,j]=last_guy
  }
  momentum=last_guy
  last_guy=rgamma(n,shape=beta_0*momentum+y[13],rate=beta_0+1)
  samples[,13]=last_guy

samples=samples[c(1:(n/10))*10,]

q025=c()
q975=c()
pred025=c()
pred975=c()
for(i in c(1:13)){
  q025=c(q025,quantile(samples[,i],0.025))
  q975=c(q975,quantile(samples[,i],0.975))
  pred_sample=rpois(n,samples[,i])
  pred025=c(pred025,quantile(pred_sample,0.025))
  pred975=c(pred975,quantile(pred_sample,0.975))
}
years=c(2008:2020)
{
ggplot()+
    geom_ribbon(aes(x=years,ymin=pred025,ymax=pred975),fill='#aaaaff',alpha=0.25)+
    geom_ribbon(aes(x=years,ymin=q025,ymax=q975),fill='#0000ff',alpha=0.25)+
    geom_line(aes(x=years,y=colMeans(samples)))+
    geom_point(aes(x=years,y=colMeans(samples)))+
    geom_point(aes(x=years,y=y),color='red')+
    scale_x_continuous('Ano',breaks=years)+
    theme_bw()
}
par(mfrow=c(4,4))
for(i in c(1:13)){
  acf(samples[,i])
}
par(mfrow=c(1,1))

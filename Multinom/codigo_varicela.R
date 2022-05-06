library(dlm)   # Pode ser substituido pela bilbioteca Matrix.
library(tidyr) # Opcional: ajuda na clareza do código da barra de loading.
library(hms)   # Opcional: serve para formatar o ETA na barra de loading.
library(beepr) # Opcional: serve para chamar a função beep.
library(feather)

# Ajustando os dados

dados=read.csv('data/varicela internacoes.csv')[,c(1,7:162)]
dados[1:2,1]='00 a 04 anos'
dados[5:8,1]='15 a 49 anos'
dados[9:11,1]='50 a 79 anos'
#dados[9:12,1]='50 anos e mais'
dados=aggregate(.~FaixaEtaria,dados,sum)[,-1]

pre_exp=read.csv2('data/populacao 2000-2020.csv')[-12,c(1,10:22)]
pre_exp[4:7,1]='15 a 49 anos'
#pre_exp[8:11,1]='50 anos e mais'
pre_exp[8:10,1]='50 a 79 anos'
pre_exp=aggregate(.~FaixaEtaria,pre_exp,sum)[,-1]

dummy=matrix(0,dim(pre_exp)[1],0)
nomes=c()
for(ano in c(2008:2020)){
  for(mes in c(1:12)){
    nomes=c(nomes,paste0('X',ano,'.',mes))
    dummy=cbind(dummy,pre_exp[,ano-2007])
  }
}
pre_exp=dummy
#idade_indice=3

t_offset=1
indice_inter=69

T_final=dim(dados)[2]

out_var=5
n_var=out_var*4+2

#### nível ####

FF_nivel=array(0,c(out_var*2,out_var,T_final))
G_nivel=diag(out_var*2)

for(i in c(1:out_var)){
  FF_nivel[(i-1)*2+1,i,]=1
  G_nivel[2*i-1,2*i]=1
}

D_nivel=matrix(1/0.95,out_var*2,out_var*2)

#### vacina ####

FF_vac=array(0,c(out_var,out_var,T_final))
for(i in out_var){
  FF_vac[i,i,(T_final-indice_inter-9+1):T_final]=1
}

G_vac=matrix(1,out_var,out_var)

D_vac=matrix(1/1,out_var,out_var)

#### covid ####

FF_cov=array(0,c(out_var,out_var,T_final))
for(i in out_var){
  FF_cov[i,i,]=c(rep(0,146),rep(1,T_final-146))
}

G_cov=matrix(1,out_var,out_var)

D_cov=matrix(1/1,out_var,out_var)

#### sazonalidade ####

FF_sazo=array(0,c(2,out_var,T_final))
FF_sazo[1,,]=1

w=2*pi/12
G_sazo=matrix(c(cos(w),sin(w),-sin(w),cos(w)),2,2)

D_sazo=matrix(1/0.95,2,2)

FF=array(0,c(n_var,out_var,T_final))
FF[1:(2*out_var),,]=FF_nivel
FF[(2*out_var+1):(3*out_var),,]=FF_vac
FF[(3*out_var+1):(4*out_var),,]=FF_cov
FF[(4*out_var+1):(4*out_var+2),,]=FF_sazo
G=bdiag(G_nivel,G_vac,G_cov,G_sazo)
D=bdiag(D_nivel,D_vac,D_cov,D_sazo)

D=array(D,c(n_var,n_var,T_final))
total=rowSums(t(dados[1:(out_var+1),]))

resultado=multinom_gi(
  y=t(dados[1:(out_var+1),]),
  m0=matrix(0,n_var,1),
  C0=diag(n_var)*50,
  FF=FF,
  G=G,
  D=D,
  W=0,
  pop=total)

ps=exp(resultado$ft)/(rowSums(exp(resultado$ft[,1:out_var]))+1)
ps=cbind(ps,1-rowSums(ps))

plot_data=dados[1:(out_var+1),]
for(i in c(1:T_final)){
  plot_data[,i]=plot_data[,i]/sum(plot_data[,i])
}

plot_data=plot_data %>%
  t %>%
  as.data.frame %>%
  mutate(t=1:T_final,) %>%
  pivot_longer(1:6) %>%
  mutate(name='V'%>%paste0(name),
         value=value)

plot_ps=ps %>%
  as.data.frame %>%
  mutate(t=1:T_final) %>%
  pivot_longer(1:6)

plot_data=plot_data %>%
  inner_join(plot_ps,c('t','name'))


ggplot(plot_data)+
  geom_point(aes(x=t,y=value.x,color=name,shape='Observado'))+
  geom_line(aes(x=t,y=value.y,color=name,linetype='Estimado'))+
  geom_vline(xintercept=indice_inter+9,linetype='dashed')+
  theme_bw()

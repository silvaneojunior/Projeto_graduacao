library(shiny)
library(ggplot2)
library(plotly)
library(stringr)
library(tidyr)
library(dplyr)
options(scipen = 9999,encoding = 'UTF-8')

mes2num=list(
  Jan='01',
  Fev='02',
  Mar='03',
  Abr='04',
  Mai='05',
  Jun='06',
  Jul='07',
  Ago='08',
  Set='09',
  Out='10',
  Nov='11',
  Dez='12'
)
get_month=function(index){mes2num[[index]]}


get_data=function(data_name){
  pre_data=read.csv(paste0('data/Mes/',data_name,'.csv'))
  pre_data[1,-1]=pre_data[1,-1]+pre_data[2,-1]
  pre_data[1,1]='00 a 04 anos'
  data=pivot_longer(pre_data[-2,],-1)
  names(data)=c('Faixa_etária','Data','Valor')
  data$Faixa_etária=factor(data$Faixa_etária)
  data$Data=as.Date(paste0(substr(data$Data,2,5),'-',sapply(substr(data$Data,7,9),get_month),'-01'))
  return(data)
}

labels_idade=c("00 a 04 anos",
               "05 a 09 anos",
               "10 a 14 anos",
               "15 a 19 anos",
               "20 a 29 anos",
               "30 a 39 anos",
               "40 a 49 anos",
               "50 a 59 anos",
               "60 a 69 anos",
               "70 a 79 anos",
               "80 anos e mais")

data=merge(get_data('varicela internacoes'),get_data('varicela obitos'),by=c("Faixa_etária",'Data'))
data$iden=1
data$mes=data$Data
data$Data=substr(data$Data,1,4)


pop=cbind(read.csv2(paste0('data/Ano/populacao 2000-2020.csv'))[1:11,-c(2:9)])
pop=pivot_longer(pop,-1)
names(pop)=c('Faixa_etária','Data','Valor')
pop$Data=substr(pop$Data,2,5)
data=merge(data,pop,by=c('Faixa_etária','Data'))
data$Data=data$mes
ref_data=data[,-6]
names(ref_data)=c('Faixa_etária','Data','Internações','Óbitos','1','Residentes')
ref_data$internacoes_por_100k=100000*ref_data$Internações/ref_data$Residentes
ref_data$Faixa_etária=as.character(ref_data$Faixa_etária)
ref_data$Faixa_etária=as.factor(ref_data$Faixa_etária)
ref_data=ref_data[order(ref_data$Faixa_etária,ref_data$Data),]

data=ref_data[,c(1,2,3)]
data=pivot_wider(data,names_from=c(1),values_from=c(3))

# Clusterização

ProporcaoEstabilidade = 0.2

# Adicionando um vetor de valores aleatórios para testar a concistência da clusteriação
#data[,1]=runif(156,1,2)
#mat_data=as.matrix(data)

mat_data=as.matrix(data[,-1])
if(F){
  mat_data[13:156,]=mat_data[13:156,]/mat_data[1:(156-12),]
  mat_data=mat_data[-c(1:12),]
  }else{
  mat_data[2:156,]=mat_data[2:156,]/mat_data[1:(156-1),]
  mat_data=mat_data[-c(1),]
  }

# Transformando dados para o formato 1,0,-1
mat_data=ifelse(mat_data>=1+ProporcaoEstabilidade,1,ifelse(mat_data<1-ProporcaoEstabilidade,-1,0))

plt_data=mat_data %>% as.data.frame %>%
  mutate(time=c(1:155)) %>%
  pivot_longer(-time)

raw_data=data[-1,] %>% as.data.frame %>%
  mutate(time=c(1:155)) %>%
  pivot_longer(-c(time,Data))

plt_data=plt_data %>% left_join(raw_data,by=c('time','name'))

ggplot(plt_data %>% filter(name %in% c('00 a 04 anos') & time<50))+
  geom_line(aes(x=time,y=mean(value.y)+sd(value.y)*2*value.x,color='Vetor estabilidade',linetype='Vetor estabilidade'))+
  geom_line(aes(x=time,y=value.y,color='Série original',linetype='Série original'))+
  scale_color_manual('',values=c('black','grey'))+
  scale_linetype_manual('',values=c('solid','dashed'))+
  scale_x_continuous('Data',breaks=c(0:5)*12+1,labels=c(0:5)+2008)+
  scale_y_continuous('')+
  theme_bw()


if(T){
  # Calculando distância com a matriz de correlação
  distancia=as.dist(1-cor(mat_data))
}else{
  # Calculando distância euclidiana entre os vetores
  distancia=dist(t(mat_data))
}
# Exibindo dendograma
plot(hclust(distancia,method="complete"),ylab='Distância',xlab='',main='')

# Gráfico com a série histórica dos dados usados na clusterização
plot_data=as.data.frame(mat_data)
plot_data$data=sort(unique(ref_data$Data))[-c(1:12)]
plot_data=pivot_longer(plot_data,c(1:11))
ggplotly(ggplot(plot_data)+
  geom_line(aes(x=data,y=value,color=name))+
    geom_hline(yintercept=1,linetype='dashed')+
    theme_bw())


# Gráfico de médias móveis
data=ref_data[,c(1,2,3)]
data=pivot_wider(data,names_from=c(1),values_from=c(3))
data=data['2011-12-01'<data$Data,]

ma_data=as.data.frame(data['2013-08-01'<data$Data & data$Data<'2020-01-01',])

for(i in c(1:dim(ma_data)[1])){
  ma_data[i,-1]=colMeans(data[0:11+i,-1])
}
for(j in c(2:dim(ma_data)[2])){
  ma_data[,j]=ma_data[,j]/ma_data[1,j]
}

plot_ma_data=pivot_longer(ma_data,2:12)

ggplotly(
  ggplot(plot_ma_data)+
    geom_line(aes(x=Data,y=value,color=name))+
    scale_x_date('Ano',date_breaks='1 year',labels=function(x){substr(x,1,4)})+
    scale_y_continuous('Quantidade de internações',
                       labels=function(x){paste0(100*x,'%')})+
    scale_color_hue('Faixa etária')+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90))
  )

ggplot(plot_ma_data)+
  geom_line(aes(x=Data,y=value))+
  scale_x_date('Ano',date_breaks='1 year',labels=function(x){substr(x,1,4)})+
  scale_y_continuous('Percentual de internações (em relação a 2012)',
                     labels=function(x){paste0(100*x,'%')},
                     breaks=c(0:6)/4,
                     limits=c(0,6)/4)+
  facet_wrap(~name)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90))


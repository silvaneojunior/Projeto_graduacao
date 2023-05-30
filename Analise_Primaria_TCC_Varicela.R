library(shiny)
library(ggplot2)
library(plotly)
library(stringr)
library(tidyr)
options(scipen = 9999,encoding = 'UTF-8')

calcula_max=function(pre_max){
  if(length(pre_max)==0){
    pre_max=10
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

set_format=function(text){
  formatC(text,0, format="d", big.mark='.',decimal.mark = ',')
}

fonte='Fonte: DATASUS, Julho de 2021' 

turn_dynamic=function(pl){return(ggplotly(pl) %>%
                                   layout(annotations = 
                                            list(x = 1, y = 0, text = fonte, 
                                                 showarrow = F, xref='paper',yref='paper',
                                                 xanchor='right', yanchor='bottom', xshift=0, yshift=0,
                                                 font=list(size=15, color="black"))
                                   ))}

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
  pre_data[1,1]='0 a 4 anos'
  data=pivot_longer(pre_data[-2,],-1)
  names(data)=c('Faixa_etária','Data','Valor')
  data$Faixa_etária=factor(data$Faixa_etária)
  data$Data=as.Date(paste0(substr(data$Data,2,5),'-',sapply(substr(data$Data,7,9),get_month),'-01'))
  return(data)
}

labels_idade=c("0 a 4 anos",
               "5 a 9 anos",
               "10 a 14 anos",
               "15 a 19 anos",
               "20 a 29 anos",
               "30 a 39 anos",
               "40 a 49 anos",
               "50 a 59 anos",
               "60 a 69 anos",
               "70 a 79 anos",
               "80 anos e mais")

data=get_data('varicela internacoes')
data=merge(data,get_data('varicela obitos'),by=c('Faixa_etária','Data'))
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
ref_data$Faixa_etária[ref_data$Faixa_etária=='5 a 9 anos']='05 a 09 anos'
ref_data$Faixa_etária=as.factor(ref_data$Faixa_etária)
ref_data=ref_data[order(ref_data$Faixa_etária,ref_data$Data),]
#2 Criando as faixas et?rias ___________________________________________
ProporcaoEstabilidade = 0.025
PEV=c(1-ProporcaoEstabilidade,1+ProporcaoEstabilidade)
Faixa1 = ref_data[1:156,7]
a=c(0)
v1=c()
for(i in 1:155){if(Faixa1[i+1]/Faixa1[i]<PEV[1]){v1[i]=-1}else if(Faixa1[i+1]/Faixa1[i]<PEV[2]){v1[i]=0}else{v1[i]=1}}
v1=c(a,v1)
plot(v1, col=2, type="b")
VetorEstabilidade1 = v1==0
VetorSubida1 = v1==1
VetorDescida1 = v1==-1
ContadorEstabilidade1 = round(sum(VetorEstabilidade1)/length(VetorEstabilidade1),4)
ContadorSubida1 = round(sum(VetorSubida1)/length(VetorSubida1),4)
ContadorDescida1 = round(sum(VetorDescida1)/length(VetorDescida1),4)
F1=c(ContadorDescida1,ContadorEstabilidade1,ContadorSubida1, sum(ContadorDescida1,ContadorEstabilidade1,ContadorSubida1))
F1
#________________________________________________________
Faixa2 = ref_data[157:312,7]
v2=c()
for(i in 1:155){if(Faixa2[i+1]/Faixa2[i]<PEV[1]){v2[i]=-1}else if(Faixa2[i+1]/Faixa2[i]<PEV[2]){v2[i]=0}else{v2[i]=1}}
v2=c(a,v2)
plot(v2, col=2, type="b")
VetorEstabilidade2 = v2==0
VetorDescida2 = v2==-1
VetorSubida2 = v2==1
ContadorEstabilidade2 = round(sum(VetorEstabilidade2)/length(VetorEstabilidade2),4)
ContadorDescida2 = round(sum(VetorDescida2)/length(VetorEstabilidade2),4)
ContadorSubida2 = round(sum(VetorSubida2)/length(VetorEstabilidade2),4)
F2 = c(ContadorDescida2,ContadorEstabilidade2,ContadorSubida2,sum(ContadorEstabilidade2,ContadorDescida2,ContadorSubida2))
F2
#_________________________________________________________
#___________________________________________________________
Faixa3 = ref_data[(156*2+1):(156*3),7]
Faixa3
v3=c()
for(i in 1:155){if(Faixa3[i+1]/Faixa3[i]<PEV[1]){v3[i]=-1}else if(Faixa3[i+1]/Faixa3[i]<PEV[2]){v3[i]=0}else{v3[i]=1}}
v3=c(a,v3)
plot(v3, col=2, type="b")
VetorEstabilidade3 = v3==0
VetorSubida3 = v3==1
VetorDescida3 = v3==-1
ContadorEstabilidade3 = round(sum(VetorEstabilidade3)/length(VetorEstabilidade3),4)
ContadorSubida3 = round(sum(VetorSubida3)/length(VetorEstabilidade3),4)
ContadorDescida3 = round(sum(VetorDescida3)/length(VetorEstabilidade3),4)
F3 = c(ContadorDescida3,ContadorEstabilidade3,ContadorSubida3,sum(ContadorDescida3,ContadorEstabilidade3,ContadorSubida3))
F3

#____________________________________________________________
Faixa4 = ref_data[(156*3+1):(156*4),7]
Faixa4
v4=c()
for(i in 1:155){if(Faixa4[i+1]/Faixa4[i]<PEV[1]){v4[i]=-1}else if(Faixa4[i+1]/Faixa4[i]<PEV[2]){v4[i]=0}else{v4[i]=1}}
v4=c(a,v4)
plot(v4, col=2, type="b")
VetorEstabilidade4 = v4==0
VetorDescida4 = v4==-1
VetorSubida4 = v4==1
ContadorEstabilidade4 = round(sum(VetorEstabilidade4)/length(VetorEstabilidade4),4)
ContadorDescida4 = round(sum(VetorDescida4)/length(VetorEstabilidade4),4)
ContadorSubida4 = round(sum(VetorSubida4)/length(VetorEstabilidade4),4)
F4=c(ContadorDescida4,ContadorEstabilidade4,ContadorSubida4,sum(ContadorDescida4,ContadorEstabilidade4,ContadorSubida4))
F4
#___________________________________________________________
Faixa5 = ref_data[(156*4+1):(156*5),7]
Faixa5
v5=c()
for(i in 1:155){if(Faixa5[i+1]/Faixa5[i]<PEV[1]){v5[i]=-1}else if(Faixa5[i+1]/Faixa5[i]<PEV[2]){v5[i]=0}else{v5[i]=1}}
v5=c(a,v5)
plot(v5, col=2, type="b")
VetorEstabilidade5 = v5==0
VetorDescida5 = v5==-1
VetorSubida5 = v5==1
ContadorEstabilidade5 = round(sum(VetorEstabilidade5)/length(VetorEstabilidade5),4)
ContadorSubida5 = round(sum(VetorSubida5)/length(VetorEstabilidade5),4)
ContadorDescida5 = round(sum(VetorDescida5)/length(VetorEstabilidade5),4)
F5=c(ContadorDescida5,ContadorEstabilidade5,ContadorSubida5,sum(ContadorDescida5,ContadorEstabilidade5,ContadorSubida5))
F5
#___________________________________________________________
Faixa6 = ref_data[(156*5+1):(156*6),7]
Faixa6
v6=c()
for(i in 1:155){if(Faixa6[i+1]/Faixa6[i]<PEV[1]){v6[i]=-1}else if(Faixa6[i+1]/Faixa6[i]<PEV[2]){v6[i]=0}else{v6[i]=1}}
v6=c(a,v6)
plot(v6, col=2, type="b")
VetorEstabilidade6 = v6==0
VetorDescida6 = v6==-1
VetorSubida6 = v6==1
ContadorEstabilidade6 = round(sum(VetorEstabilidade6)/length(VetorEstabilidade6),4)
ContadorSubida6 = round(sum(VetorSubida6)/length(VetorEstabilidade5),4)
ContadorDescida6 = round(sum(VetorDescida6)/length(VetorEstabilidade5),4)
F6=c(ContadorDescida6,ContadorEstabilidade6,ContadorSubida6,sum(ContadorDescida6,ContadorEstabilidade6,ContadorSubida6))
F6
#____________________________________________________________
Faixa7 = ref_data[(156*6+1):(156*7),7]
Faixa7
v7=c()
for(i in 1:155){if(Faixa7[i+1]/Faixa7[i]<PEV[1]){v7[i]=-1}else if(Faixa7[i+1]/Faixa7[i]<PEV[2]){v7[i]=0}else{v7[i]=1}}
v7=c(a,v7)
plot(v7, col=2, type="b")
VetorEstabilidade7 = v7==0
VetorDescida7 = v7==-1
VetorSubida7 = v7==1
ContadorEstabilidade7 = round(sum(VetorEstabilidade7)/length(VetorEstabilidade5),4)
ContadorSubida7 = round(sum(VetorSubida7)/length(VetorEstabilidade5),4)
ContadorDescida7 = round(sum(VetorDescida7)/length(VetorEstabilidade5),4)
F7=c(ContadorDescida7,ContadorEstabilidade7,ContadorSubida7,sum(ContadorDescida7,ContadorEstabilidade7,ContadorSubida7))
F7
#____________________________________________________________
Faixa8 = ref_data[(156*7+1):(156*8),7]
Faixa8
v8=c()
for(i in 1:155){if(Faixa8[i+1]/Faixa8[i]<PEV[1]){v8[i]=-1}else if(Faixa8[i+1]/Faixa8[i]<PEV[2]){v8[i]=0}else{v8[i]=1}}
v8=c(a,v8)
plot(v8, col=2, type="b")
VetorEstabilidade8 = v8==0
VetorDescida8 = v8==-1
VetorSubida8 = v8==1
ContadorEstabilidade8 = round(sum(VetorEstabilidade8)/length(VetorEstabilidade5),4)
ContadorSubida8 = round(sum(VetorSubida8)/length(VetorEstabilidade5),4)
ContadorDescida8 = round(sum(VetorDescida8)/length(VetorEstabilidade5),4)
F8=c(ContadorDescida8,ContadorEstabilidade8,ContadorSubida8,sum(ContadorDescida8,ContadorEstabilidade8,ContadorSubida8))
F8
#____________________________________________________________
Faixa9 = ref_data[(156*8+1):(156*9),7]
Faixa9
v9=c()
for(i in 1:155){if(Faixa9[i+1]/Faixa9[i]<PEV[1]){v9[i]=-1}else if(Faixa9[i+1]/Faixa9[i]<PEV[2]){v9[i]=0}else{v9[i]=1}}
v9=c(a,v9)
plot(v9, col=2, type="b")
VetorEstabilidade9 = v9==0
VetorDescida9 = v9==-1
VetorSubida9 = v9==1
ContadorEstabilidade9 = round(sum(VetorEstabilidade9)/length(VetorEstabilidade5),4)
ContadorSubida9 = round(sum(VetorSubida9)/length(VetorEstabilidade9),4)
ContadorDescida9 = round(sum(VetorDescida9)/length(VetorEstabilidade9),4)
F9=c(ContadorDescida9,ContadorEstabilidade9,ContadorSubida9,sum(ContadorDescida9,ContadorEstabilidade9,ContadorSubida9))
F9
#____________________________________________________________
Faixa10 = ref_data[(156*9+1):(156*10),7]
Faixa10
v10=c()
for(i in 1:155){if(Faixa10[i+1]/Faixa10[i]<PEV[1]){v10[i]=-1}else if(Faixa10[i+1]/Faixa10[i]<PEV[2]){v10[i]=0}else{v10[i]=1}}
v10=c(a,v10)
plot(v10, col=2, type="b")
VetorEstabilidade10 = v10==0
VetorDescida10 = v10==-1
VetorSubida10 = v10==1
ContadorEstabilidade10 = round(sum(VetorEstabilidade10)/length(VetorEstabilidade5),4)
ContadorSubida10 = round(sum(VetorSubida10)/length(VetorEstabilidade9),4)
ContadorDescida10 = round(sum(VetorDescida10)/length(VetorEstabilidade9),4)
F10=c(ContadorDescida10,ContadorEstabilidade10,ContadorSubida10,sum(ContadorDescida10,ContadorEstabilidade10,ContadorSubida10))
F10
#____________________________________________________________
Faixa11 = ref_data[(156*10+1):(156*11),7]
Faixa11
v11=c()
for(i in 1:155){if(Faixa11[i+1]/Faixa11[i]<PEV[1]){v11[i]=-1}else if(Faixa11[i+1]/Faixa11[i]<PEV[2]){v11[i]=0}else{v11[i]=1}}
v11=c(a,v11)
plot(v11, col=2, type="b")
VetorEstabilidade11 = v11==0
VetorDescida11 = v11==-1
VetorSubida11 = v11==1
ContadorEstabilidade11 = round(sum(VetorEstabilidade11)/length(VetorEstabilidade5),4)
ContadorSubida11 = round(sum(VetorSubida11)/length(VetorEstabilidade9),4)
ContadorDescida11 = round(sum(VetorDescida11)/length(VetorEstabilidade9),4)
F11=c(ContadorDescida11,ContadorEstabilidade11,ContadorSubida11,sum(ContadorDescida11,ContadorEstabilidade11,ContadorSubida11))
F11



F1.1 = F1[1:3]
F2.1 = F2[1:3]
F3.1 = F3[1:3]
F4.1 = F4[1:3]
F5.1 = F5[1:3]
F6.1 = F6[1:3]
F7.1 = F7[1:3]
F8.1 = F8[1:3]
F9.1 = F9[1:3]
F10.1 = F10[1:3]
F11.1 = F11[1:3]



MatrizEstabilidade = matrix(rbind(F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11),nrow=11,ncol=4)
rownames(MatrizEstabilidade) = c("0 a 4 anos",
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
colnames(MatrizEstabilidade) = c("Prop. Baixa","Prop. Estabilidade", "Prop. Alta","Total")
MatrizEstabilidade
pie(F1.1,labels=F1.1,radius=1,col = c("green","yellow","blue"),main=c("De 0 a 4 Anos"))
pie(F2.1,labels=F2.1,radius=1,col = c("green","yellow","blue"),main=c("De 5 a 9 Anos"))
pie(F3.1,labels=F3.1,radius=1,col = c("green","yellow","blue"),main=c("De 10 a 14 Anos"))
pie(F4.1,labels=F4.1,radius=1,col = c("green","yellow","blue"),main=c("De 15 a 19 Anos"))
pie(F5.1,labels=F5.1,radius=1,col = c("green","yellow","blue"),main=c("20 a 29 Anos"))
pie(F1.1,labels=F1.1,radius=1,col = c("green","yellow","blue"),main=c("De 30 a 39 Anos"))
pie(F2.1,labels=F2.1,radius=1,col = c("green","yellow","blue"),main=c("De 40 a 49 Anos"))
pie(F3.1,labels=F3.1,radius=1,col = c("green","yellow","blue"),main=c("De 50 a 59 Anos"))
pie(F4.1,labels=F4.1,radius=1,col = c("green","yellow","blue"),main=c("De 60 a 69 Anos"))
pie(F5.1,labels=F5.1,radius=1,col = c("green","yellow","blue"),main=c("De 70 a 79 Anos"))
pie(F5.1,labels=F5.1,radius=1,col = c("green","yellow","blue"),main=c("80 Anos e Mais"))

savehistory("~/R/history.Rhistory")

par(mfrow=c(1,1))
Variacoes = cbind(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11)
colnames(Variacoes) = c("0 a 4 anos",
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
Correlacoes = 1-cor(Variacoes)
Distancias = as.dist(Correlacoes)
Clusters = hclust(Distancias)
plot(Clusters)


vetortesteI = runif(156, 0, 1)
vetortesteF = c()
for(i in 1:156){if(vetortesteI[i]<=ProporcaoEstabilidade){vetortesteF[i]=0}else if(vetortesteI[i]<=(ProporcaoEstabilidade+((1-ProporcaoEstabilidade)/2))){vetortesteF[i]=-1}else{vetortesteF[i]=1}}
vetortesteF
sum(vetortesteF)
sum(vetortesteF==0)
VariacoesTeste = cbind(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,vetortesteF)
colnames(VariacoesTeste) = c("0 a 4 anos",
                                 "05 a 09 anos",
                                 "10 a 14 anos",
                                 "15 a 19 anos",
                                 "20 a 29 anos",
                                 "30 a 39 anos",
                                 "40 a 49 anos",
                                 "50 a 59 anos",
                                 "60 a 69 anos",
                                 "70 a 79 anos",
                                 "80 anos e mais","VetorAleatório")
CorrelacoesTeste = 1-cor(VariacoesTeste)
DistanciaTeste = as.dist(CorrelacoesTeste)
ClustersTeste = hclust(DistanciaTeste, "average")
plot(ClustersTeste)
v12 = vetortesteF
cor(VariacoesTeste)

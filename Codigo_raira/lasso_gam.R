library(dlm)   # Pode ser substituido pela bilbioteca Matrix.
library(tidyr) # Opcional: ajuda na clareza do código da barra de loading.
library(hms)   # Opcional: serve para formatar o ETA na barra de loading.
library(beepr) # Opcional: serve para chamar a função beep.
library(mgcv)
library(plsmselect)
source('C:\\Jupyter\\TCC\\Projeto_graduacao\\Codigo_raira\\codigo_raira_VSilvaneo.R')

# Ajustando os dados
raw_dados=read.csv('Wania/data/dados_wania_completo.CSV',row.names=1)
raw_dados=t(raw_dados) %>% as.data.frame
dados=raw_dados
dados$Tempo=c(1:dim(dados)[1])
dados$N_OBITOS=dados$N_OBITOS %>% as.integer
dados$log_pop=dados$Pac_dia %>% log
dados$intervencao=c(rep(0,54-1),rep(1,101-54+1))
dados$prorp_covid=dados$Entradas_COVID/dados$Entradas

new_vars=c('intervencao',
           'prorp_covid',
           'N_ERC',
           'N_PaR_carba',
           'N_ESBL',
           'Mero',
           'PoliB',
           'Zerbaxa',
           'Torgena',
           'Amica',
           'Aztreo',
           'Tige',
           'Ceftaroline',
           'Pipetazo',
           'Erta',
           'Cipro',
           'Soma_Antimicro'
           )

for(var in new_vars){
  dados[[var]]=(dados[[var]]-mean(dados[[var]]))/sd(dados[[var]])
}

dados$X=model.matrix(~Mero+PoliB+Zerbaxa+Torgena+Amica+Aztreo+Tige+Ceftaroline+Pipetazo+Erta+Cipro+intervencao+prorp_covid+N_ERC+N_PaR_carba+N_ESBL, data=dados)[,-1]
#+intervencao+prorp_covid+N_ERC+N_PaR_carba+N_ESB

suv_t=gamlasso(N_OBITOS~X+s(Tempo),
               family="poisson" ,data=dados,offset='log_pop',linear.penalty='l1')
plot(dados$Tempo,dados$N_OBITOS)
lines(dados$Tempo,suv_t$gam$fitted.values)

summary(suv_t)

#### Versão alternativa com análise de fatores e PCA ####

clust_data=dados[,names(dados) %in% new_vars]
eigen_decomp=cor(clust_data) %>% eigen
cumsum(eigen_decomp$values/sum(eigen_decomp$values))

factanal(clust_data,6)

var=new_vars[2]

plot_data=dados[,names(dados) %in% c('N_OBITOS','Tempo',new_vars)]

plot_data[['N_OBITOS']]=(plot_data[['N_OBITOS']]-mean(plot_data[['N_OBITOS']]))/sd(plot_data[['N_OBITOS']])

plot_data=plot_data %>% pivot_longer(-c(Tempo))

colors=rainbow(length(new_vars),s=0.6)
names(colors)=new_vars
colors[['N_OBITOS']]='black'

ggplotly(
ggplot(plot_data)+
  geom_point(aes(x=Tempo,y=value,color=name))+
  scale_color_manual(values=colors)+
  theme_bw()
)

cor_dados=raw_dados
cor_dados$N_OBITOS=cor_dados$N_OBITOS/cor_dados$Pac_dia
cor(cor_dados[cor_dados$Entradas_COVID==0,])[1,]

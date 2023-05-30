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

dados$outcome=dados$Soma_Antimicro-dados$Zerbaxa-dados$Torgena

suv_t=gam(outcome~intervencao+prorp_covid+s(Tempo),
               family=gaussian ,data=dados)
plot(dados$Tempo,dados$outcome)
lines(dados$Tempo,suv_t$fitted.values)

summary(suv_t)

cut_in=29
cut_out=101

dados$outcome=dados$Pipetazo

suv_t=gam(outcome[cut_in:cut_out]~intervencao[cut_in:cut_out]+prorp_covid[cut_in:cut_out]+s(Tempo[cut_in:cut_out]),
          family=gaussian ,data=dados)
plot(dados$Tempo[cut_in:cut_out],dados$outcome[cut_in:cut_out])
lines(dados$Tempo[cut_in:cut_out],suv_t$fitted.values)

summary(suv_t)


rowSums(dados[,9:20])

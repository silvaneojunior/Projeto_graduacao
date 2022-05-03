library(dlm)   # Pode ser substituido pela bilbioteca Matrix.
library(tidyr) # Opcional: ajuda na clareza do código da barra de loading.
library(hms)   # Opcional: serve para formatar o ETA na barra de loading.
library(beepr) # Opcional: serve para chamar a função beep.

source('C:\\Jupyter\\TCC\\Projeto_graduacao\\Codigo_raira\\codigo_raira_VSilvaneo.R')
print(getwd())
# Ajustando os dados
dados=read.csv2('Wania/data/dados_wania_covid.CSV',row.names=1)[-1,]

obt=as.numeric(dados[1,1:(101-3)]) # Variável alvo
exp=as.numeric(dados[2,1:(101-3)])#/1000 # Exposição
NERC=as.numeric(dados[3,1:(101-3)])
PaRCarba=as.numeric(dados[4,1:(101-3)])
ESBL=as.numeric(dados[5,1:(101-3)])
var_names=c('obt','exp','NERC','ParCarba','ESBL')

cut_data=dados[3:5,1:(101-3)]

N <- length(obt) # Tamanho da série temporal
indice_inter=54 # ìndice da intervenção (Junho de 2017)
pre_covid=as.numeric(dados[7,1:(101-3)])
indice_covid=min(c(1:(101-3))[ifelse(is.na(pre_covid),0,pre_covid)>0])

# Criando função para calcular o Error dado uma escolha de parâmetros.
roda_modelo=function(y_resp,flag_multi,delta_m1,delta_inter,delta_covid,indic_covid,usa_vel){
  n <- 3 + usa_vel # Número de variáveis
  indic_inter=c(rep(0,indice_inter-1),rep(1,N-indice_inter+1))
  
  covid=as.numeric(dados[as.numeric(indic_covid),1:(101-3)])
  var_covid=ifelse(is.na(covid),0,covid)/ifelse(indic_covid=='7',as.numeric(dados[6,1:(101-3)]),exp) # Remove NA's no dados.
  
  nivel=gera_bloco_poly(1+usa_vel,D=1/delta_m1)
  inter=gera_bloco_poly(1,value=indic_inter,D=1/delta_inter,m0=0,C0=0)
  covid=gera_bloco_poly(1,value=var_covid,D=1/delta_covid,m0=0,C0=0)

  inter$D[,,1:indice_inter] <- 1
  inter$D[,,indice_inter] <- 1/0.7
  inter$W[,,indice_inter] <- 1

  covid$D[,,1:indice_covid] <- 1
  covid$D[,,indice_covid] <- 1/0.7
  covid$W[,,indice_covid] <- 1
  
  estrutura=concat_bloco(nivel,inter,covid)

  if(flag_multi=='relativo'){
    NERC_bloc=gera_bloco_poly(1,value=NERC/exp,D=1/1)
    PaRCarba_bloc=gera_bloco_poly(1,value=PaRCarba/exp,D=1/1)
    ESBL_bloc=gera_bloco_poly(1,value=ESBL/exp,D=1/1)

    estrutura=concat_bloco(estrutura,NERC_bloc,PaRCarba_bloc,ESBL_bloc)
  }else{
    if(flag_multi=='TRUE'){
      NERC_bloc=gera_bloco_poly(1,value=NERC,D=1/1)
      PaRCarba_bloc=gera_bloco_poly(1,value=PaRCarba,D=1/1)
      ESBL_bloc=gera_bloco_poly(1,value=ESBL,D=1/1)

      estrutura=concat_bloco(estrutura,NERC_bloc,PaRCarba_bloc,ESBL_bloc)
    }
  }
  
  teste <- ajusta_modelo(data_out=y_resp,struture=estrutura,offset=exp,kernel=poisson_testa_par)
  
  # Error calculado a partir do 10º instante de tempo.
  erros=y_resp-teste$pred
  erros=abs(erros)/teste$pred
  
  erros=erros[-c(1:(10-1))]
  erro=mean(erros)
  return(erro)
}

# Valor central inicial da busca
m0=rep(0.9,3)
search_radio=0.1
initial_granul=10

initial_range=c(-initial_granul:initial_granul)/(initial_granul)

# Matriz para armazenar os valores do Error em cada combinação de parâmtros testado.
data_mat=as.data.frame(matrix(0,1,8),)
names(data_mat)=c('resp','flag_multi','Delta.M1','Delta.Inter','Delta.COV','Indic_COVID','Usa_Vel','Error')
data_mat$Error=Inf

count=0 # Contagem do progresso

init=Sys.time() # Armazenando o valor do tempo de início.
# O loop a seguir testa os valores para diversas combinações de parâmetros
for(resp in c(1)){
best_m=m0
best_e=Inf
y_resp=as.numeric(dados[resp,1:(101-3)])
multi_range=if(resp==1){c(FALSE)}else{c(FALSE)}
for(flag_multi in multi_range){
for(d_m1 in m0[1]+search_radio*initial_range){
for(d_cov in m0[2]+search_radio*initial_range){
for(d_in in m0[3]+search_radio*initial_range){
for(indic_covid in c(7)){
for(use_vel in c(TRUE,FALSE)){
  count=count+1
  
  tryCatch({
    # Tenta rodar o modelo
    erro=roda_modelo(y_resp,flag_multi,d_m1,d_in,d_cov,indic_covid,use_vel)
    erro=ifelse(is.nan(erro),Inf,erro)
  },
  error=function(e){
    # Caso ocorra um problema, o erro é fixado em "ínfinito".
    print('a)')
    erro=Inf
  },
  finally={
    data_mat=rbind(data_mat,c(var_names[resp],flag_multi,d_m1,d_in,d_cov,indic_covid,use_vel,erro))
    ##################################################################################
    ##### Este treicho do código é opcional: serve pra gerar a barra de loading  #####
    #perc=count/((1*2+3*1)*2*2*(2*initial_granul+1)**3)
    perc=count/(2*(2*initial_granul+1)**3)
    cur_time=Sys.time()
    
    qtd1=min(49,floor(49*perc))
    qtd2=49-qtd1
    
    cat(paste0('[',paste0(rep('=',qtd1),collapse = ''),'>',paste0(rep(' ',qtd2),collapse = ''),']  ',
               (100*perc) %>% round(2) %>% format(nsmall=2) %>% paste0('%'),
               ' - ETA: ',((1-perc)*difftime(cur_time, init, unit="secs")/perc) %>% as.numeric %>% round %>% hms,
               ' - Melhor: ',best_e %>% round(2),
               ' - Atual: ',erro %>% round(2),
               '\r'))
    if(best_e>erro | best_e==Inf){
      best_m=c(d_m1,d_in,d_cov)
      best_e=erro
    }
  ##################################################################################
  }
  )
}}}}}}}
# O valor central da busca é alterado ao final de cada etapa, se tornando igual ao melhor conjunto de hiper-parâmetros encontrado.
m0=best_m
beep()
print(Sys.time()-init)
# Ordenando os dados pelo Error
data_mat=data_mat[order(data_mat$Error,decreasing=FALSE),]

# Exibindo melhor conjunto de dados encontrado
print(data_mat[1,])

write.csv2(data_mat,'Wania/grid_data/MRE_par_com_multiresistente.csv',row.names = F)

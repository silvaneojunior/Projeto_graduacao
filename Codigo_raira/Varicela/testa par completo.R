library(dlm)   # Pode ser substituido pela bilbioteca Matrix.
library(tidyr) # Opcional: ajuda na clareza do código da barra de loading.
library(hms)   # Opcional: serve para formatar o ETA na barra de loading.
library(beepr) # Opcional: serve para chamar a função beep.
library(feather)

source('C:\\Jupyter\\TCC\\Projeto_graduacao\\Codigo_raira\\codigo_raira_VSilvaneo.R')
# Ajustando os dados

dados=read.csv('Varicela/data/varicela internacoes.csv')[,c(1,7:162)]
dados[1:2,1]='00 a 04 anos'
dados[5:8,1]='15 a 49 anos'
dados[9:11,1]='50 a 79 anos'
#dados[9:12,1]='50 anos e mais'
dados=aggregate(.~FaixaEtaria,dados,sum)[,-1]

pre_exp=read.csv2('Varicela/data/populacao 2000-2020.csv')[-12,c(1,10:22)]
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

# Valor central inicial da busca
m0=c(rep(0.9,4),0.9)
search_radio=0.1
initial_granul=5

initial_range=c(-initial_granul:initial_granul)/(initial_granul)

t_offset=1
name='MRE_pred_all_delay_ferias_agrup2'
indice_inter=69

# Criando função para calcular o Error dado uma escolha de parâmetros.
roda_modelo=function(delta_m1,delta_inter,delta_covid,usa_vel,primavera_flag,delta_primavera,sazo_size,delta_sazo,delay){
  true_indice_inter=indice_inter+delay
  indic_inter=c(rep(0,true_indice_inter-1),rep(1,N-true_indice_inter+1))
  covid=c(rep(0,146),rep(1,N-146))
  
  var_covid=ifelse(is.na(covid),0,covid) # Remove NA's no dados.
  
  nivel_bloc=gera_bloco_poly(1+usa_vel,D=1/delta_m1,C0=1)
  inter_bloc=gera_bloco_poly(1,value=indic_inter,D=1/delta_inter,C0=0)
  inter_bloc$D[,,1:true_indice_inter] <- 1
  inter_bloc$D[,,true_indice_inter] <- 1/0.7
  inter_bloc$W[,,true_indice_inter] <- 1
  
  covid_bloc=gera_bloco_poly(1,value=var_covid,D=1/delta_covid,C0=0)
  covid_bloc$D[,,1:146] <- 1
  covid_bloc$D[,,146] <- 1/0.8
  covid_bloc$W[,,146] <- 1
  
  indic_ferias=rep(0,N)
  indic_ferias[c(1:N)[c(1:N)%%12 %in% c(1,2,6,12)]]=1
  indic_ferias[79]=0
  indic_ferias[78]=1
  ferias_bloc=gera_bloco_poly(1,value=indic_ferias,D=1/1)
  
  estrutura=concat_bloco(nivel_bloc,inter_bloc,covid_bloc,ferias_bloc)
  #estrutura=concat_bloco(nivel_bloc,inter_bloc,covid_bloc)
  
  if(primavera_flag!='Sem'){
    if(primavera_flag=='Transição'){
      indic_primavera=((c(1:N)-9)%%12==0) %>% as.numeric
    }else{
      indic_primavera=((c(1:N)-9)%%12==0 | (c(1:N)-10)%%12==0 | (c(1:N)-11)%%12==0) %>% as.numeric
    }
    
    prima=gera_bloco_poly(1,value=indic_primavera)
    prima$D[,,indic_primavera==1] <- 1/delta_primavera
    estrutura=concat_bloco(estrutura,prima)
  }
  for(sazo in sazo_size){
    estrutura=concat_bloco(estrutura,gera_bloco_sazo(sazo,D=1/delta_sazo))
  }
  
  ################################################################################

  
  error_offset=24
  
  teste <- ajusta_modelo(data_out=inter,struture=estrutura,offset=exp,kernel=poisson_testa_par)
  t_offset=1
  pred=eval_past(teste,smooth = F,t_offset=t_offset)$pred
  objt=inter[-c(1:t_offset)]
  erros=objt-pred
  erros=abs(erros)/objt
  erros=erros[-c(1:(error_offset))]
  erro1=mean(erros,na.rm=T)
  
  t_offset=6
  pred=eval_past(teste,smooth = F,t_offset=t_offset)$pred
  objt=inter[-c(1:t_offset)]
  erros=objt-objt
  erros=abs(erros)/pred
  erros=erros[-c(1:(error_offset))]
  erro2=mean(erros,na.rm=T)
  
  t_offset=12
  pred=eval_past(teste,smooth = F,t_offset=t_offset)$pred
  objt=inter[-c(1:t_offset)]
  erros=objt-objt
  erros=abs(erros)/pred
  erros=erros[-c(1:(error_offset))]
  erro3=mean(erros,na.rm=T)
  
  # P=teste$var.pred
  # G=(inter-teste$pred)**2
  # erro=P+G
  # erro=sum(erro[-c(1:(10-1))])
  
  # erros=pnbinom(inter, teste$alpha, (teste$beta/(teste$beta +1)),log=TRUE)
  # erro=sum(erros,na.rm=T)

  return(list('erro1'=erro1,'erro2'=erro2,'erro3'=erro3))
  
}

sazo_list=list('Sem'=c(),
               'Anual'=c(12))

basic_loop=3*length(initial_range)**3
sazo_loop=(3*length(initial_range)+1)
prima_loop=(length(initial_range)+1)
total=basic_loop*sazo_loop*prima_loop

age_len=dim(dados)[1]

for(idade_indice in c(age_len:age_len)){
  inter=as.numeric(dados[idade_indice,])
  exp=as.data.frame(pre_exp)[idade_indice,]
  names(exp)=nomes
  exp=as.numeric(exp)
  
  N <- dim(dados)[2]
  data_lab=names(dados)
  
  m0=rep(0.9,5)
  # Valor inicial do melhor parâmtros e erro.
  best_m=m0
  best_e=Inf
  
  # Matriz para armazenar os valores do Error em cada combinação de parâmtros testado.
  fileConn1<-file(paste0('Varicela/grid_data/',name,'_1_',idade_indice,'.csv'))
  open(fileConn1,'w')
  writeLines(paste('Delta.M1','Delta.Inter','Delta.COV','Usa_Vel','Primavera_flag','Delta.prima','Sazo.Cycle','Delta.sazo','Delay','Error',sep=','),fileConn1)

  # Matriz para armazenar os valores do Error em cada combinação de parâmtros testado.
  fileConn2<-file(paste0('Varicela/grid_data/',name,'_6_',idade_indice,'.csv'))
  open(fileConn2,'w')
  writeLines(paste('Delta.M1','Delta.Inter','Delta.COV','Usa_Vel','Primavera_flag','Delta.prima','Sazo.Cycle','Delta.sazo','Delay','Error',sep=','),fileConn2)
  
  
  # Matriz para armazenar os valores do Error em cada combinação de parâmtros testado.
  fileConn3<-file(paste0('Varicela/grid_data/',name,'_12_',idade_indice,'.csv'))
  open(fileConn3,'w')
  writeLines(paste('Delta.M1','Delta.Inter','Delta.COV','Usa_Vel','Primavera_flag','Delta.prima','Sazo.Cycle','Delta.sazo','Delay','Error',sep=','),fileConn3)
  
  count=0 # Contagem do progresso
  init=Sys.time() # Armazenando o valor do tempo de início.
  
  total=0
  {for(delay in c(0:12)){
    for(use_vel in c(TRUE,FALSE)){
    for(d_m1 in m0[1]+search_radio*initial_range){
      for(d_cov in c(1)){
        for(d_in in c(1)){
          for(sazo_size in names(sazo_list)){
            seach_range=if(sazo_size==0){c(0)}else{search_radio*initial_range/2}
            for(d_sazo in m0[4]+0.05+seach_range){
              for(primavera_flag in c('Sem')){
                seach_range_p=if(primavera_flag=='Sem'){c(0)}else{search_radio*initial_range/2}
                for(d_prima in m0[5]+0.05+seach_range_p){
                  total=total+1
                  }}}}}}}}}}
  
  for(delay in c(0:12)){  
  for(use_vel in c(TRUE,FALSE)){
  for(d_m1 in m0[1]+search_radio*initial_range){
  for(d_cov in c(1)){
  for(d_in in c(1)){
  for(sazo_size in names(sazo_list)){
    seach_range=if(sazo_size==0){c(0)}else{search_radio*initial_range/2}
  for(d_sazo in m0[4]+0.05+seach_range){
  for(primavera_flag in c('Sem')){
    seach_range_p=if(primavera_flag=='Sem'){c(0)}else{search_radio*initial_range/2}
  for(d_prima in m0[5]+0.05+seach_range_p){
    count=count+1
    
    sazo_cicle=sazo_list[[as.character(sazo_size)]]
    erros=roda_modelo(d_m1,d_in,d_cov,use_vel,primavera_flag,d_prima,sazo_cicle,d_sazo,delay)

    writeLines(paste(d_m1,d_in,d_cov,use_vel,primavera_flag,d_prima,sazo_size,d_sazo,delay,erros$erro1,sep=','),fileConn1)
    writeLines(paste(d_m1,d_in,d_cov,use_vel,primavera_flag,d_prima,sazo_size,d_sazo,delay,erros$erro2,sep=','),fileConn2)
    writeLines(paste(d_m1,d_in,d_cov,use_vel,primavera_flag,d_prima,sazo_size,d_sazo,delay,erros$erro3,sep=','),fileConn3)
    
    perc=count/total
    cur_time=Sys.time()
    
    qtd1=min(49,floor(49*perc))
    qtd2=49-qtd1
    
    cat(paste0('[',paste0(rep('=',qtd1),collapse = ''),'>',paste0(rep(' ',qtd2),collapse = ''),']  ',
               paste(count,'/',total) %>% paste0(' ',(100*perc) %>% round(2) %>% format(nsmall=2) %>% paste0('%')),
               ' - ETA: ',((1-perc)*difftime(cur_time, init, unit="secs")/perc) %>% as.numeric %>% round %>% hms,
               '\r'))
    ##################################################################################
  }}}}}}}}}
    
  beep()

  close(fileConn1)
  data_mat=read.csv(paste0('Varicela/grid_data/',name,'_1_',idade_indice,'.csv'))
  data_mat=data_mat[order(data_mat$Error,decreasing=FALSE),]
  write.csv(data_mat,paste0('Varicela/grid_data/',name,'_1_',idade_indice,'.csv'),row.names = F)
  write_feather(data_mat,paste0('Varicela/grid_data/',name,'_1_',idade_indice,'.feather'))
  
  close(fileConn2)
  data_mat=read.csv(paste0('Varicela/grid_data/',name,'_6_',idade_indice,'.csv'))
  data_mat=data_mat[order(data_mat$Error,decreasing=FALSE),]
  write.csv(data_mat,paste0('Varicela/grid_data/',name,'_6_',idade_indice,'.csv'),row.names = F)
  write_feather(data_mat,paste0('Varicela/grid_data/',name,'_6_',idade_indice,'.feather'))
  
  close(fileConn3)
  data_mat=read.csv(paste0('Varicela/grid_data/',name,'_12_',idade_indice,'.csv'))
  data_mat=data_mat[order(data_mat$Error,decreasing=FALSE),]
  write.csv(data_mat,paste0('Varicela/grid_data/',name,'_12_',idade_indice,'.csv'),row.names = F)
  write_feather(data_mat,paste0('Varicela/grid_data/',name,'_12_',idade_indice,'.feather'))
}
# library(ggplot2)
# library(plotly)
# 
# idade_indice=1
# data_mat=read.csv2(data_mat,paste0('MSE_par_misto_idade',idade_indice,'.csv'),row.names = F)
# 
# ggplot(data_mat[data_mat$Delta.COV==1 & data_mat$Usa_Vel==0 & data_mat$Sazo.Cycle==0 & data_mat$Delta.sazo==0.8 & data_mat$Ganul==1,])+
#   geom_tile(aes(x=Delta.M1,y=Delta.Inter,fill=MSE))+
#   scale_fill_gradientn(colors=rainbow(10))+
#   theme_bw()
# 
# plot_dados=data_mat[data_mat$Ganul==1,]
# plot_dados$Sazo.Cycle=as.factor(plot_dados$Sazo.Cycle)
# ggplot(plot_dados)+
#   geom_boxplot(aes(x=Sazo.Cycle,y=MSE,group=Sazo.Cycle,fill=Sazo.Cycle))+
#   theme_bw()
# 
# plot_dados=data_mat[data_mat$Ganul==1,]
# plot_dados$Sazo.Cycle=as.factor(plot_dados$Sazo.Cycle)
# ggplotly(
#   ggplot(plot_dados)+
#   geom_density(aes(x=MSE,color=Sazo.Cycle,fill=Sazo.Cycle))+
#   theme_bw()
# )

# Ajustando os dados

# dados=read.csv('varicela internacoes.csv')[,c(1,7:162)]
# dados[1:3,1]='00 a 09 anos'
# dados[5:8,1]='15 a 49 anos'
# dados[9:12,1]='50 anos e mais'
# dados=aggregate(.~FaixaEtaria,dados,sum)[,-1]
# 
# pre_exp=read.csv2('populacao 2000-2020.csv')[-12,c(1,10:22)]
# pre_exp[1:2,1]='00 a 09 anos'
# pre_exp[4:7,1]='15 a 49 anos'
# pre_exp[8:11,1]='50 anos e mais'
# pre_exp=aggregate(.~FaixaEtaria,pre_exp,sum)[,-1]
# 
# dummy=matrix(0,dim(dados)[1],0)
# nomes=c()
# for(ano in c(2008:2020)){
#   for(mes in c(1:12)){
#     nomes=c(nomes,paste0('X',ano,'.',mes))
#     dummy=cbind(dummy,pre_exp[,ano-2007])
#   }
# }
# pre_exp=dummy
# #idade_indice=3
# for(idade_indice in c(1:dim(dados)[1])){
#   inter=as.numeric(dados[idade_indice,])
#   exp=as.data.frame(pre_exp)[idade_indice,]
#   names(exp)=nomes
#   exp=as.numeric(exp)
#   
#   N <- dim(dados)[2]
#   indice_inter=69+6
#   data_lab=names(dados)
#   
#   # Valor inicial do melhor parâmtros e erro.
#   m0=rep(0.9,5)
#   best_m=m0
#   best_e=Inf
#   
#   # Matriz para armazenar os valores do Error em cada combinação de parâmtros testado.
#   fileConn<-file(paste0('Group2\\MRE_par',idade_indice,'.csv'))
#   open(fileConn,'w')
#   writeLines(paste('Ganul','Delta.M1','Delta.Inter','Delta.COV','Usa_Vel','Primavera_flag','Delta.prima','Sazo.Cycle','Delta.sazo','Error',sep=','),fileConn)
#   
#   # Granul controla o refinamento do Grid, quanto mais fina a malha da busca.
#   for(granul in c(1:1)){
#     count=0 # Contagem do progresso
#     init=Sys.time() # Armazenando o valor do tempo de início.
#     # O loop a seguir testa os valores para diversas combinações de parâmetros
#     
#     for(use_vel in c(TRUE,FALSE)){
#       # use_vel indica se um Modelo de crescimento linear será usado
#       # d_m1 é o fator de desconto associado ao nível;
#       for(d_m1 in m0[1]+search_radio*initial_range){
#         # d_cov é o fator de desconto associado ao efeito do COVID;
#         for(d_cov in m0[2]+search_radio*initial_range){
#           # d_in é o fator de desconto associado à intervenção;
#           for(d_in in m0[3]+search_radio*initial_range){
#             # indic_covid é o índice da variável COVID;
#             for(sazo_size in names(sazo_list)){
#               if(sazo_size==0){
#                 seach_range=c(0)
#               }else{
#                 seach_range=search_radio*initial_range/2
#               }
#               for(d_sazo in m0[4]+0.05+seach_range){
#                 for(primavera_flag in c('Sem','Transição')){
#                   if(primavera_flag=='Sem'){
#                     seach_range_p=c(0)
#                   }else{
#                     seach_range_p=search_radio*initial_range/2
#                   }
#                   for(d_prima in m0[5]+0.05+seach_range_p){
#                     count=count+1
#                     
#                     d_m1=min(1,d_m1)
#                     d_in=min(1,d_in)
#                     d_cov=min(1,d_cov)
#                     d_sazo=min(1,d_sazo)
#                     d_prima=min(1,d_prima)
#                     sazo_cicle=sazo_list[[as.character(sazo_size)]]
#                     
#                     if(sazo_size>0){
#                       erro=roda_modelo(d_m1,d_in,d_cov,use_vel,primavera_flag,d_prima,sazo_cicle,d_sazo)
#                     }
#                     
#                     writeLines(paste(granul,d_m1,d_in,d_cov,use_vel,primavera_flag,d_prima,sazo_size,d_sazo,erro,sep=','),fileConn)
#                     #data_mat[count,]=c(granul,d_m1,d_in,d_cov,use_vel,sazo_size,d_sazo,erro)
#                     ##################################################################################
#                     ##### Este treicho do código é opcional: serve pra gerar a barra de loading  #####
#                     perc=count/total
#                     cur_time=Sys.time()
#                     
#                     qtd1=min(49,floor(49*perc))
#                     qtd2=49-qtd1
#                     
#                     cat(paste0('[',paste0(rep('=',qtd1),collapse = ''),'>',paste0(rep(' ',qtd2),collapse = ''),']  ',
#                                paste(count,'/',total) %>% paste0(' ',(100*perc) %>% round(2) %>% format(nsmall=2) %>% paste0('%')),
#                                ' - ETA: ',((1-perc)*difftime(cur_time, init, unit="secs")/perc) %>% as.numeric %>% round %>% hms,
#                                ' - Melhor: ',best_e %>% round(4),
#                                ' - Atual: ',erro %>% round(4),
#                                '\r'))
#                     if(best_e>erro){
#                       best_m=c(d_m1,d_in,d_cov,d_sazo,d_prima)
#                       best_e=erro
#                     }
#                     ##################################################################################
#                   }
#                 }
#               }
#             }
#           }
#         }
#       }
#     }
#     # O valor central da busca é alterado ao final de cada etapa, se tornando igual ao melhor conjunto de hiper-parâmetros encontrado.
#     m0=best_m
#     beep() # Esta linha é opcional: serve para gerar o aviso sonoro ao término da busca.
#     #search_radio=0.75*search_radio/initial_granul
#   }
#   close(fileConn)
#   data_mat=read.csv(paste0('Group2\\MRE_par',idade_indice,'.csv'))
#   
#   # Ordenando os dados pelo Error
#   data_mat=data_mat[order(data_mat$Error,decreasing=FALSE),]
#   
#   # Exibindo melhor conjunto de dados encontrado
#   print(data_mat[1,])
#   
#   write.csv(data_mat,paste0('Group2\\MRE_par',idade_indice,'.csv'),row.names = F)
# }
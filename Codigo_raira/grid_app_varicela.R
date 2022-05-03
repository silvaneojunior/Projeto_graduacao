library(dlm)
library(shiny)
library(ggplot2)
library(plotly)
library(tidyr)
library(dplyr)
library(latex2exp)
library(feather)

source('codigo_raira_VSilvaneo.R',encoding='UTF-8')

dados=read.csv('Varicela/data/varicela internacoes.csv')[,c(1,7:162)]
dados[1:2,1]='00 a 04 anos'
dados[5:8,1]='15 a 49 anos'
dados[9:12,1]='50 anos e mais'
dados=aggregate(.~FaixaEtaria,dados,sum)
faixa_list=dados$FaixaEtaria
dados=dados[,-1]

pre_exp=read.csv2('Varicela/data/populacao 2000-2020.csv')[-12,c(1,10:22)]
pre_exp[4:7,1]='15 a 49 anos'
pre_exp[8:11,1]='50 anos e mais'
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

sazo_list=list('Sem'=c(),
               'Anual'=c(12),
               'Semestral'=c(6),
               'Trimestral'=c(3),
               'Anual+Semestral'=c(12,6),
               'Anual+Trimestral'=c(12,3),
               'Semestral+Trimestral'=c(6,3))

N <- dim(dados)[2]
indice_inter=69
data_lab=names(dados)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel(""),
    
    withMathJax(),
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            radioButtons("grid_type",
                         "Versão do grid:",
                                  list(
                                      'Com delay'='MRE_pred_all_delay',
                                      'Com delay e ferias'='MRE_pred_all_delay_ferias'
                                  ),
                                  selected='MRE_pred_all_delay'),
            radioButtons("age_index",
                         "Idade:",
                         list(
                             '0 a 4 anos'=1,
                             '5 a 9 anos'=2,
                             '10 a 14 anos'=3,
                             '15 a 49 anos'=4,
                             '50 anos e mais'=5
                         ),
                         selected=1),
            radioButtons("Var1",
                         'Variável fixa',
                         choices=list(
                             'Nível'='Delta.M1',
                             'Sazonalidade'='Delta.sazo',
                             'Defasagem'='Delay'
                         ),
                         selected='Delta.sazo'),
            uiOutput('slider_Var1'),
            uiOutput('radio_vel'),
            radioButtons("pred_offset",
                        'Offset da previsão',
                        choices=list(
                            '1 mês à frente'=1,
                            '6 meses à frente'=6,
                            '1 ano à frente'=12
                        ),
                        selected=1),
            uiOutput('radio_var')
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotlyOutput("distPlot",width='800px',height ='600px'),
            actionButton("update_values",'Atualizar parâmetros'),
            plotlyOutput("plot_pred",width='800px',height ='500px'),
            fluidRow(downloadButton('save_pred','Salvar dados')),
            #fluidRow(plotlyOutput("plot_tx",width='800px',height ='500px')),
            tabsetPanel(
                tabPanel('Filtragem',
                         fluidRow(plotlyOutput("plot_var_filt",width='800px',height ='500px')),
                         fluidRow(downloadButton('save_var_filt','Salvar dados'))
                ),
                tabPanel('Suavização',
                         fluidRow(plotOutput("plot_var_suv",width='800px',height ='500px')),
                         fluidRow(downloadButton('save_var_suv','Salvar dados'))
                ),
                selected='Suavização'
            ),
            tabsetPanel(
                tabPanel('Filtragem',
                         fluidRow(plotlyOutput("plot_perc_filt",width='800px',height ='500px'))
                ),
                tabPanel('Suavização',
                         fluidRow(plotlyOutput("plot_perc_suv",width='800px',height ='500px'))
                ),
                selected='Suavização'
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    coords=reactiveValues(x_cut=0,y_cut=1)
    
    data_mat=eventReactive(c(input$age_index,input$pred_offset,input$grid_type),{
        req(input$age_index)
        data_mat=read_feather(paste0('Varicela/grid_data/',input$grid_type,'_',input$pred_offset,'_',input$age_index,'.feather')) %>%
            filter(Error<Inf) %>%
            mutate(Usa_Vel=Usa_Vel %>% as.logical) %>%
            arrange(desc(-Error))
        data_mat
    })
    expo=eventReactive(c(input$age_index),{
        expo=as.data.frame(pre_exp)[input$age_index %>% as.numeric,]
        names(expo)=nomes
        as.numeric(expo)
        })
    inter=eventReactive(c(input$age_index),{
        as.numeric(dados[input$age_index %>% as.numeric,])
    })
    
    observeEvent(c(plot_data()),{
        name_varplot=c('Delay','Delta.M1','Delta.sazo')
        name_varplot=name_varplot[name_varplot!=input$Var1]
        coords$x_cut=plot_data()[1,][[name_varplot[1]]]
        coords$y_cut=plot_data()[1,][[name_varplot[2]]]
    })
    observeEvent(event_data('plotly_click',source='heatmap', priority = "event"),{
        coords$x_cut=event_data('plotly_click',source='heatmap', priority = "event")$x
        coords$y_cut=event_data('plotly_click',source='heatmap', priority = "event")$y
        
        coords$x_cut=min(12,coords$x_cut)
        coords$x_cut=max(0,coords$x_cut)
        
        coords$y_cut=min(1,coords$y_cut)
        coords$y_cut=max(0.8,coords$y_cut)
    },ignoreNULL = TRUE)
    
    observeEvent(input$update_values,{
        plot_data=data_mat()
        
        values=plot_data[plot_data$Delay==coords$x_cut,]
        values=values[1,]
       
        updateRadioButtons(
            inputId="Usa_Vel",
            selected = values$Usa_Vel
        )
        updateSliderInput(
            inputId="Val1",
            value = values[[input$Var1]]
        )
        coords$y_cut=values$Delta.M1
        
    },ignoreNULL = TRUE)
    
    Dados_DLM=eventReactive(c(as.logical(input$Usa_Vel),coords$x_cut,coords$y_cut,input$age_index,input$Var1,input$Val1,input$grid_type),{
        req(input$Val1)
        name_varplot=c('Delay','Delta.M1','Delta.sazo')
        name_varplot=name_varplot[name_varplot!=input$Var1]
        
        plot_data=data_mat()
        min_Error=min(plot_data$Error,na.rm=TRUE)
        max_Error=max(plot_data$Error,na.rm=TRUE)
        
        plot_data=plot_data[plot_data$Usa_Vel==as.numeric(as.logical(input$Usa_Vel)) &
                                plot_data$Sazo.Cycle=='Anual' &
                                plot_data[[name_varplot[1]]]==coords$x_cut &
                                plot_data[[name_varplot[2]]]==coords$y_cut &
                                plot_data[[input$Var1]]==input$Val1,]
        ref_data=plot_data
        
        delta_m1=ref_data$Delta.M1[1]
        delta_Inter=ref_data$Delta.Inter[1]
        delta_Cov=ref_data$Delta.COV[1]
        delta_prima=ref_data$Delta.prima[1]
        delta_sazo=ref_data$Delta.sazo[1]
        indic_vel=ref_data$Usa_Vel[1]
        sazo_size=sazo_list[['Anual']]
        
        covid=c(rep(0,146),rep(1,N-146))
        Exp_covid=ifelse(is.na(covid),0,covid)
        
        indice_inter_true=indice_inter+ref_data$Delay
        
        indic_inter=c(rep(0,indice_inter_true-1),rep(1,N-indice_inter_true+1))
        var_covid=Exp_covid
        
        nivel=gera_bloco_poly(1+indic_vel,name='Nível',D=1/delta_m1,C0=1)#,C0=10**2)
        inter=gera_bloco_poly(1,value=indic_inter,name='Vacina',D=1/delta_Inter,C0=0)#,C0=0.2**2)
        inter$D[,,1:indice_inter_true] <- 1
        inter$D[,,indice_inter_true] <- 1/0.7
        inter$W[,,indice_inter_true] <- 1
        
        covid=gera_bloco_poly(1,value=var_covid,name='COVID',D=1/delta_Cov,C0=0)#,C0=0.2**2)
        covid$D[,,1:146] <- 1
        covid$D[,,146] <- 1/0.8
        covid$W[,,146] <- 1
        
        estrutura=concat_bloco(nivel,inter,covid)
        
        if(input$grid_type=='MRE_pred_all_delay_ferias'){
            indic_ferias=rep(0,N)
            indic_ferias[c(1:N)[c(1:N)%%12 %in% c(1,2,6,12)]]=1
            indic_ferias[79]=0
            indic_ferias[78]=1
            ferias_bloc=gera_bloco_poly(1,value=indic_ferias,D=1/1,name='Férias')
            
            estrutura=concat_bloco(estrutura,ferias_bloc)
        }
        
        for(sazo in sazo_size){
            estrutura=concat_bloco(estrutura,gera_bloco_sazo(sazo,D=1/delta_sazo,name='Sazonalidade'))
        }
        
        teste=ajusta_modelo(data_out=inter(),struture=estrutura,offset=expo())
        
        return(teste)
        
    })
    
    output$plot_pred <- renderPlotly({
        
        indice_inter_true=indice_inter+ifelse(input$Var1=='Delay',input$Val1,coords$x_cut)
        
        max_value=calcula_max(Dados_DLM()$data_out-min(Dados_DLM()$data_out))[[3]]+min(Dados_DLM()$data_out)
        min_value=-calcula_max(-(Dados_DLM()$data_out-max(Dados_DLM()$data_out)))[[3]]+max(Dados_DLM()$data_out)
        
        pre_labels=data_lab[c(c(1:13)*12+1)]
        pre_labels=substr(pre_labels,2,5)
        
        ggplotly(
        show_fit(Dados_DLM(),smooth = F,t_offset=as.numeric(input$pred_offset),dinamic=FALSE)$plot+
        scale_x_continuous('Data',breaks=c(c(1:13)*12+1),labels=pre_labels)+
        geom_vline(xintercept=indice_inter_true,linetype='dashed')+
        theme(axis.text.x = element_text(angle = 90))+
        coord_cartesian(xlim=c(10+as.numeric(input$pred_offset),dim(Dados_DLM()$mt)[2]),
                        ylim=c(min_value,max_value))
        )
    })
    
    output$plot_tx <- renderPlotly({
        indice_inter_true=indice_inter+ifelse(input$Var1=='Delay',input$Val1,coords$x_cut)
        
        r=show_fit(Dados_DLM(),smooth = T,dinamic=TRUE)$r
        p=show_fit(Dados_DLM(),smooth = T,dinamic=TRUE)$p
        
        mean=r*(1-p)/p
        lim_i=qnbinom(0.025,r,p)
        lim_s=qnbinom(0.975,r,p)
        
        mean=mean/Dados_DLM()$offset
        lim_i=lim_i/Dados_DLM()$offset
        lim_s=lim_s/Dados_DLM()$offset
        
        t_last=length(mean)
        
        max_value=calcula_max(mean-min(mean))[[3]]+min(mean)
        min_value=0
        
        fill_list=c('#2596be','#2596be','black')
        names(fill_list)=c(paste0('I.C. (95%)'),'Taxa de internações estimada','Taxa de internações observada')
        color_list=c('#2596be','black')
        names(color_list)=c('Taxa de internações estimada','Taxa de internações observada')
        
        plt=ggplot()+
            geom_point(aes(x=c(1:t_last),y=mean,color='Taxa de internações estimada',fill='Taxa de internações estimada'))+
            geom_ribbon(aes(x=c(1:t_last),ymin=lim_i,ymax=lim_s,fill=paste0('I.C. (95%)'),color=paste0('I.C. (95%)')),alpha=0.25)+
            geom_point(aes(x=c(1:t_last),y=Dados_DLM()$data_out/Dados_DLM()$offset,color='Taxa de internações observada',fill='Taxa de internações observada'))+
            scale_fill_manual('',na.value=NA,values=fill_list)+
            scale_color_manual('',na.value=NA,values=color_list)+
            scale_y_continuous(name='Taxa de internações')+
            theme_bw()+
            scale_x_continuous('Data',breaks=c(c(1:13)*12+1,indice_inter_true),labels=data_lab[c(c(1:13)*12+1,indice_inter_true)])+
            geom_vline(xintercept=indice_inter_true,linetype='dashed')+
            theme(axis.text.x = element_text(angle = 90))+
            coord_cartesian(xlim=c(10+as.numeric(input$pred_offset),dim(Dados_DLM()$mt)[2]),
                            ylim=c(min_value,max_value))
    })
    
    output$plot_var_filt <- renderPlotly({
        req(input$Usa_Vel)
        req(input$var_index)
        indice_inter_true=indice_inter+ifelse(input$Var1=='Delay',input$Val1,coords$x_cut)
            ggplotly(
                plot_lat_var(Dados_DLM(),input$var_index,smooth = F,dinamic=FALSE)$plot+
                    scale_x_continuous('Data',breaks=c(c(1:13)*12+1,indice_inter_true),labels=data_lab[c(c(1:13)*12+1,indice_inter_true)])+
                    geom_vline(xintercept=indice_inter_true,linetype='dashed')+
                    theme(axis.text.x = element_text(angle = 90))
            )
    })
    output$plot_var_suv <- renderPlotly({
        req(input$Usa_Vel)
        req(input$var_index)

        indice_inter_true=indice_inter+ifelse(input$Var1=='Delay',input$Val1,coords$x_cut)
        ggplotly(
            plot_lat_var(Dados_DLM(),input$var_index,smooth = T,dinamic=FALSE)$plot+
                scale_x_continuous('Data',breaks=c(c(1:13)*12+1,indice_inter_true),labels=data_lab[c(c(1:13)*12+1,indice_inter_true)])+
                geom_vline(xintercept=indice_inter_true,linetype='dashed')+
                theme(axis.text.x = element_text(angle = 90))
        )
    })
    # output$plot_var_suv <- renderPlot({
    #     req(input$Usa_Vel)
    #     req(input$var_index)
    #     req(input$Usa_sazo)
    #     
    #     dados_suv=eval_past(Dados_DLM(),smooth=TRUE)
    #     
    #     erro=dados_suv$pred-Dados_DLM()$data_out
    #     var_erro=dados_suv$a*(dados_suv$b+1)/(dados_suv$b)^2
    #     erro=erro[-c(1:12)]
    #     var_erro=var_erro[-c(1:12)]
    #     
    #     erro=erro/sqrt(var_erro)
    #     #erro=erro/sd(erro)
    #     
    #     p_obs=c(1:length(erro))/(length(erro)+1)
    #     q_teo=qnorm(p_obs)
    #     
    #     max_value=calcula_max(abs(erro))[[2]]
    #     print(max_value)
    #     
    #     ggplot()+
    #         geom_line(aes(x=c(-max_value,max_value),y=c(-max_value,max_value)))+
    #         geom_point(aes(x=q_teo,y=sort(erro)))+
    #         scale_x_continuous('Quantis observados',expand=c(0,0))+
    #         scale_y_continuous('Quantis teóricos',expand=c(0,0))+
    #         theme_bw()+
    #         coord_fixed()
    #     
    #     # hist(erro,breaks=20 ,freq =FALSE)
    #     # lines(c(-50:50)/10,dnorm(c(-50:50)/10))
    # })
    output$plot_var_suv <- renderPlot({
        req(input$Usa_Vel)
        req(input$var_index)

        indice_inter_true=indice_inter+ifelse(input$Var1=='Delay',input$Val1,coords$x_cut)

        pre_labels=data_lab[c(c(1:13)*12+1)]
        pre_labels=substr(pre_labels,2,5)

            plot_lat_var(Dados_DLM(),input$var_index,smooth = T,dinamic=FALSE)$plot+
                scale_x_continuous('Data',breaks=c(c(1:13)*12+1),labels=pre_labels,expand=c(0,0))+
                scale_y_continuous('Valor do parâmetro')+
                geom_vline(xintercept=indice_inter_true,linetype='dashed')+
                guides(color='none',fill='none')+
                theme(axis.text.x = element_text(angle = 90),text = element_text(size = 20))+
                labs(title=input$var_index)
    })
    
    output$save_pred <- downloadHandler(
        filename = function() {
            paste('prediction one step ahead ',faixa_list[input$age_index %>% as.numeric] ,".csv",sep='')
        },
        content = function(file) {
            
            eval=eval_past(Dados_DLM(),FALSE,1)
            
            prediction=eval$pred
            alpha=eval$a
            beta=eval$b
            
            r=alpha
            p=(beta/(beta +1))
            
            observed_value=Dados_DLM()$data_out
            icl.pred<-qnbinom((1-0.95)/2, r, p)
            icu.pred<-qnbinom(1-(1-0.95)/2, r, p)
            time=c(1:length(prediction))
            
            table_data=data.frame(time,observed_value,prediction,icl.pred,icu.pred,r,p)
            
            write.csv2(table_data, file, row.names = FALSE)
        }
    )
    output$save_var_suv <- downloadHandler(
        filename = function() {
            paste(input$var_index,faixa_list[input$age_index %>% as.numeric] ,"smooth.csv")
        },
        content = function(file) {
            write.csv2(plot_lat_var(Dados_DLM(),input$var_index,smooth = T,dinamic=FALSE)$table, file, row.names = FALSE)
        }
    )
    output$save_var_filt <- downloadHandler(
        filename = function() {
            paste(input$var_index,faixa_list[input$age_index %>% as.numeric],"filtered.csv")
        },
        content = function(file) {
            write.csv2(plot_lat_var(Dados_DLM(),input$var_index,smooth = F,dinamic=FALSE)$table, file, row.names = FALSE)
        }
    )
    output$plot_perc_filt <- renderPlotly({
        req(input$Usa_Vel)
        req(input$var_index)
        indice_inter_true=indice_inter+ifelse(input$Var1=='Delay',input$Val1,coords$x_cut)
        
        ggplotly(
            plot_lat_var(Dados_DLM(),input$var_index,smooth = F,exp_y=TRUE,dinamic=F)$plot+
                labs(title='Efeito percentual (' %>% paste(input$var_index,')'))+
                geom_vline(xintercept=indice_inter_true,linetype='dashed')+
                scale_x_continuous('Data',breaks=c(c(1:13)*12+1,indice_inter_true),labels=data_lab[c(c(1:13)*12+1,indice_inter_true)])+
                scale_y_continuous('$\\Delta y_t\\%$',labels=function(x){(100*x) %>% paste0('%')})
        )
    })
    output$plot_perc_suv <- renderPlotly({
        req(input$Usa_Vel)
        req(input$var_index)
        indice_inter_true=indice_inter+ifelse(input$Var1=='Delay',input$Val1,coords$x_cut)
        
        ggplotly(
            plot_lat_var(Dados_DLM(),input$var_index,smooth = T,exp_y=TRUE,dinamic=F)$plot+
                labs(title='Efeito percentual (' %>% paste(input$var_index,')'))+
                geom_vline(xintercept=indice_inter_true,linetype='dashed')+
                scale_x_continuous('Data',breaks=c(c(1:13)*12+1,indice_inter_true),labels=data_lab[c(c(1:13)*12+1,indice_inter_true)])+
                scale_y_continuous('$\\Delta y_t\\%$',labels=function(x){(100*x) %>% paste0('%')})
        )
    })
    # output$plot_perc_suv <- renderPlot({
    #     req(input$Usa_Vel)
    #     req(input$var_index)
    #     req(input$Usa_sazo)
    #     indice_inter_true=indice_inter+coords$x_cut
    #     
    #     pre_labels=data_lab[c(c(1:13)*12+1)]
    #     pre_labels=substr(pre_labels,2,5)
    # 
    #         plot_lat_var(Dados_DLM(),input$var_index,smooth = T,tranform_y=function(y){(exp(y)-1)},dinamic=F)$plot+
    #             labs(title='Efeito percentual (' %>% paste0(input$var_index,')'))+
    #             geom_vline(xintercept=indice_inter_true,linetype='dashed')+
    #             scale_x_continuous('Data',breaks=c(c(1:13)*12+1),labels=pre_labels,expand=c(0,0))+
    #             scale_y_continuous('$\\Delta y\\%$' %>% TeX,labels=function(x){(100*x) %>% paste0('%')})+
    #             guides(color='none',fill='none')+
    #             theme(axis.text.x = element_text(angle = 90),text = element_text(size = 20))
    # })
    
    
    output$slider_Var1 = renderUI({
        req(input$Usa_Vel)
        req(input$Var1)
        plot_data=data_mat()
        values=as.numeric(plot_data[plot_data$Usa_Vel==as.logical(input$Usa_Vel) & plot_data$Sazo.Cycle=='Anual',][1,])
        names(values)=names(plot_data)
        
        var_values=plot_data[[input$Var1]]
        
        minimo=min(var_values)
        minimo2=min(var_values[var_values>minimo])
        
        slider2=sliderInput('Val1',
                            input$Var1,
                            min = minimo,
                            max = max(plot_data[[input$Var1]]),
                            value = values[[input$Var1]],
                            step=round(minimo2-minimo,2)
        )
        return(fluidRow(slider2))
    })
    output$radio_var = renderUI({
        req(input$Usa_Vel)
        options=names(Dados_DLM()$names)
        
        radioButtons("var_index",
                    "Variável para o gráfico de efeitos:",
                    choices=options,
                    selected='Vacina')
    })
    output$radio_vel = renderUI({
        req(data_mat())
        radioButtons("Usa_Vel",
                     "Estrutura do nível:",
                     choices=list(
                         'Modelo de 1ª ordem'=FALSE,
                         'Modelo de crescimento linear'=TRUE
                     ),
                     selected=data_mat()$Usa_Vel[1])
    })
    
    plot_data=eventReactive(c(input$Indic_COV,input$Usa_Vel,input$Var1,input$Val1,input$age_index,input$grid_type),{
        plot_data=data_mat()
        req(input$Usa_Vel)
        req(input$Val1)

        plot_data=plot_data[plot_data$Usa_Vel==as.logical(input$Usa_Vel) & plot_data$Sazo.Cycle=='Anual' & plot_data$Primavera_flag=='Sem',]
        plot_data=plot_data[plot_data[[input$Var1]]==input$Val1,]
        return(plot_data)
    })
    output$distPlot <- renderPlotly({
        name_varplot=c('Delay','Delta.M1','Delta.sazo')
        label_varplot=c('Defasagem','Fator de desconto para o nível','Fator de desconto para a sazonalidade')
        label_varplot=label_varplot[name_varplot!=input$Var1]
        name_varplot=name_varplot[name_varplot!=input$Var1]
        
        plot_data=data_mat()
        min_Error=min(plot_data$Error,na.rm=TRUE)
        max_Error=max(plot_data$Error,na.rm=TRUE)
        
        plot_data=plot_data()
        #plot_data$MSE=ifelse(plot_data$MSE>min(plot_data$MSE),plot_data$MSE,0)
        x_cut=coords$x_cut
        y_cut=coords$y_cut
        
        ref_Error=min(plot_data[plot_data[[name_varplot[1]]]==x_cut & plot_data[[name_varplot[2]]]==y_cut,]$Error,na.rm=TRUE)
        placeholder=plot_data[plot_data[[name_varplot[1]]]==x_cut & plot_data[[name_varplot[2]]]==y_cut,]

        x_min=plot_data()[1,][[name_varplot[1]]]
        y_min=plot_data()[1,][[name_varplot[2]]]
        label_name=ifelse(x_min==x_cut & y_min==y_cut,'Mínimo: ','Error: ')
        plot_data$Error_color=plot_data$Error
        
        ggplotly(
            ggplot(plot_data)+
                geom_tile(aes_string(x=name_varplot[1],y=name_varplot[2],fill='Error_color',value='Error'))+
                geom_text(aes(x=ifelse(x_cut+0.025+0.02>1,x_cut-0.025,x_cut+0.025),y=ifelse(y_cut+0.01+0.02>1,y_cut-0.01,y_cut+0.01),label=paste0(label_name,ref_Error %>% round(2))))+
                geom_point(aes(x=x_cut,y=y_cut))+
                geom_point(aes(x=x_min,y=y_min))+
                scale_fill_gradientn('EAR',colours=c('red','yellow','green'),values=c(0,0.5,1))+#,limits=c(0.1,0.7))+
                geom_hline(yintercept=y_cut,linetype='dashed')+
                geom_vline(xintercept=x_cut,linetype='dashed')+
                #labs(title='Este gráfico é clicável')+
                scale_x_continuous(label_varplot[1],expand=c(0,0),breaks=c(c(0:12),x_cut))+
                scale_y_continuous(label_varplot[2],expand=c(0,0),breaks=c(c(16:20)/20,y_cut))+
                theme_bw(),
            source='heatmap'
        )
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

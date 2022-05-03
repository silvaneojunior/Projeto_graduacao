#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(dlm)
library(shiny)
library(shinyjs)
library(ggplot2)
library(plotly)

source('codigo_raira_VSilvaneo.R',encoding='UTF-8')

dados=read.csv2('Wania/data/dados_wania_covid.CSV',row.names=1)[-1,]

data_lab=names(dados)[1:(101-3)]
obt=as.numeric(dados[1,1:(101-3)])
expo=as.numeric(dados[2,1:(101-3)])#/1000
NERC=as.numeric(dados[3,1:(101-3)])
PaRCarba=as.numeric(dados[4,1:(101-3)])
ESBL=as.numeric(dados[5,1:(101-3)])

var_names=c('obt','expo','NERC','ParCarba','ESBL')

N <- length(obt)
inter_indic=54

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel(""),
    useShinyjs(),
    withMathJax(),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
        radioButtons("y_resp",
                         "Variável resposta:",
                         choices=list(
                             'Óbitos'=1,
                             'Nº ERC'=3,
                             'PaRCarba'=4,
                             'ESBL'=5
                         ),
                         selected=1),
        radioButtons("Multi_flag",
                     "Usar multiresistentes como regressora?",
                     choices=list(
                         'Sim (relativisado pela exposição)'='relativo',
                         'Sim (Bruto)'='TRUE',
                         'Não'='FALSE'
                     ),
                     selected='FALSE'),
        radioButtons("var1",
                         "Variável fixa:",
                         choices=list(
                             'Fator de desconto para o nível'=3,
                             'Fator de desconto para a intervenção'=4,
                             'Fator de desconto para a COVID'=5
                         ),
                         selected=5),
        uiOutput('slider_Var1'),
        radioButtons("Usa_Vel",
                     "Estrutura do nível:",
                     choices=list(
                         'Modelo de 1ª ordem'=0,
                         'Modelo de crescimento linear'=1
                     ),
                     selected=0),
        radioButtons("Indic_COV",
                     "Forma de entrada da COVID:",
                     choices=list(
                         'Entradas'=7,
                         'Pacientes-dia'=8
                     )),
        sliderInput("pred_offset",
                    'Offset da previsão',
                    min = 1,
                    max = 24,
                    value = 1,
                    step=1
        ),
        uiOutput('radio_var')
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotlyOutput("distPlot",width='675px'),
           plotlyOutput("plot_pred",width='675px'),
           fluidRow(downloadButton('save_pred','Salvar dados')),
           tabsetPanel(
               tabPanel('Filtragem',
                        fluidRow(plotlyOutput("plot_var_filt",width='675px')),
                        fluidRow(downloadButton('save_var_filt','Salvar dados'))
               ),
               tabPanel('Suavização',
                        fluidRow(plotlyOutput("plot_var_suv",width='675px')),
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
    coords=reactiveValues(x_cut=0.88,y_cut=0.92)
    
    data_mat=read.csv2('Wania/grid_data/MRE_par_com_multiresistente.csv')
    data_mat$Error=data_mat$Error %>% as.numeric
    data_mat$Delta.M1=data_mat$Delta.M1 %>% as.numeric
    data_mat$Delta.Inter=data_mat$Delta.Inter %>% as.numeric
    data_mat$Delta.COV=data_mat$Delta.COV %>% as.numeric
    data_mat=data_mat[data_mat$Error<Inf & !is.nan(data_mat$Error),]
    
    observeEvent(c(input$var1,input$val1,input$y_resp,input$Multi_flag,input$Usa_Vel),{
        if(input$y_resp==1){
            enable('Multi_flag')
        }else{
            if(input$Multi_flag=='TRUE'){
                updateRadioButtons(inputId='Multi_flag', selected = 'FALSE')
            }
            disable('Multi_flag')
        }
        name_varplot=names(data_mat)[3:5][-(as.numeric(input$var1)-2)]
        coords$x_cut=plot_data()[1,][[name_varplot[1]]]
        coords$y_cut=plot_data()[1,][[name_varplot[2]]]
        })
    
    observeEvent(event_data('plotly_click',source='heatmap', priority = "event"),{
        coords$x_cut=event_data('plotly_click',source='heatmap', priority = "event")$x
        coords$y_cut=event_data('plotly_click',source='heatmap', priority = "event")$y
        coords$x_cut=min(coords$x_cut,1)
        coords$x_cut=max(coords$x_cut,0.8)
        coords$y_cut=min(coords$y_cut,1)
        coords$y_cut=max(coords$y_cut,0.8)
    },ignoreNULL = TRUE)
    
    Dados_DLM=eventReactive(c(input$y_resp,input$Multi_flag,input$var1,input$Usa_Vel,input$Indic_COV,input$val1,coords$x_cut,coords$y_cut),{
        req(input$val1)
        req(input$var1)
        name_var1=names(data_mat)[as.numeric(input$var1)]
        name_varplot=names(data_mat)[3:5][-(as.numeric(input$var1)-2)]
        
        plot_data=data_mat[data_mat$Indic_COV==input$Indic_COV & data_mat$resp==var_names[input$y_resp %>% as.numeric] & data_mat$flag_multi==input$Multi_flag,]
        min_Error=min(plot_data$Error)
        max_Error=max(plot_data$Error)
        
        plot_data=plot_data[plot_data$Usa_Vel==ifelse(input$Usa_Vel==0,'FALSE','TRUE') & plot_data[,as.numeric(input$var1)]==input$val1,]
        ref_data=plot_data[1,-c(1,2)]
        ref_data[[name_varplot[1]]]=coords$x_cut
        ref_data[[name_varplot[2]]]=coords$y_cut
        
        delta_m1=ref_data[1,1]
        delta_Inter=ref_data[1,2]
        delta_Cov=ref_data[1,3]
        indic_vel=ref_data[1,5]
        indic_vel=ifelse(indic_vel %>% as.logical,1,0)
        
        y_resp=as.numeric(dados[input$y_resp %>% as.numeric,1:(101-3)])
        covid=as.numeric(dados[7,1:(101-3)])
        Exp_covid=ifelse(is.na(covid),0,covid)
        Exp_covid=Exp_covid/ifelse(input$Indic_COV=='7',as.numeric(dados[6,1:(101-3)]),expo)
        
        n <- 3+indic_vel
        indice_inter=inter_indic
        indic_inter=c(rep(0,indice_inter-1),rep(1,N-indice_inter+1))
        var_covid=Exp_covid
        
        nivel=gera_bloco_poly(1+indic_vel,name='Nível',D=1/delta_m1)
        inter=gera_bloco_poly(1,value=indic_inter,name='Intervenção',D=1/delta_Inter)
        covid=gera_bloco_poly(1,value=var_covid,name='COVID',D=1/delta_Cov)
        
        estrutura=concat_bloco(nivel,inter,covid)
        
        if(input$Multi_flag=='TRUE'){
            NERC_bloc=gera_bloco_poly(1,value=NERC,name='NERC',D=1/1)
            PaRCarba_bloc=gera_bloco_poly(1,value=PaRCarba,name='PaRCarba',D=1/1)
            ESBL_bloc=gera_bloco_poly(1,value=ESBL,name='ESBL',D=1/1)
            
            estrutura=concat_bloco(estrutura,NERC_bloc,PaRCarba_bloc,ESBL_bloc)
        }else{
            if(input$Multi_flag=='relativo'){
                NERC_bloc=gera_bloco_poly(1,value=NERC/expo,name='NERC',D=1/1)
                PaRCarba_bloc=gera_bloco_poly(1,value=PaRCarba/expo,name='PaRCarba',D=1/1)
                ESBL_bloc=gera_bloco_poly(1,value=ESBL/expo,name='ESBL',D=1/1)
                
                estrutura=concat_bloco(estrutura,NERC_bloc,PaRCarba_bloc,ESBL_bloc)
            }
        }
        
        # C0[n-1,n-1]=C0[n-1,n-1]/30
        # C0[n,n]=C0[n,n]/30
        
        # estrutura$D[n-1,n-1,1:inter_indic] <- 1
        # estrutura$D[n-1,n-1,inter_indic] <- 1/0.8
        # estrutura$D[n,n,1:87] <- 1
        # estrutura$D[n,n,87] = 1/0.8
        
        ajusta_modelo(data_out=y_resp,struture=estrutura,offset=expo)
    })
    
    output$plot_pred <- renderPlotly({
        max_value=calcula_max(Dados_DLM()$data_out-min(Dados_DLM()$data_out))[[3]]+min(Dados_DLM()$data_out)
        min_value=-calcula_max(-(Dados_DLM()$data_out-max(Dados_DLM()$data_out)))[[3]]+max(Dados_DLM()$data_out)
        
        
        labels_ano=data_lab[c(c(1:8)*12+1)]
        labels_ano='20' %>% paste0(labels_ano %>% substr(5,6))
        
        ggplotly(
        show_fit(Dados_DLM(),smooth = F,t_offset=input$pred_offset,dinamic=FALSE)$plot+
            scale_x_continuous('Data',breaks=c(c(1:8)*12+1),labels=labels_ano)+
            geom_vline(xintercept=inter_indic,linetype='dashed')+
            theme(axis.text.x = element_text(angle = 90))+
            coord_cartesian(xlim=c(10+input$pred_offset,dim(Dados_DLM()$mt)[2]),
                            ylim=c(max(0,min_value),max_value))
        )
    })
    output$plot_var_filt <- renderPlotly({
        req(input$Usa_Vel)
        req(input$var_index)
        
        labels_ano=data_lab[c(c(1:8)*12+1,inter_indic)]
        labels_ano='20' %>% paste0(labels_ano %>% substr(4,5))
        
        ggplotly(
        plot_lat_var(Dados_DLM(),input$var_index,smooth = F,dinamic=FALSE)$plot+
        scale_x_continuous('Data',breaks=c(c(1:8)*12+1,inter_indic),labels=labels_ano)+
        geom_vline(xintercept=inter_indic,linetype='dashed')+
        theme(axis.text.x = element_text(angle = 90))+
        coord_cartesian(xlim=c(10+input$pred_offset,dim(Dados_DLM()$mt)[2]))
        )
    })
    output$plot_var_suv <- renderPlotly({
        req(input$Usa_Vel)
        req(input$var_index)
        
        ggplotly(
        plot_lat_var(Dados_DLM(),input$var_index,smooth = T,dinamic=FALSE)$plot+
        scale_x_continuous('Data',breaks=c(c(1:8)*12+1,inter_indic),labels=data_lab[c(c(1:8)*12+1,inter_indic)])+
        geom_vline(xintercept=inter_indic,linetype='dashed')+
        theme(axis.text.x = element_text(angle = 90))+
        coord_cartesian(xlim=c(10+input$pred_offset,dim(Dados_DLM()$mt)[2]))
        )
    })
    
    output$save_pred <- downloadHandler(
        filename = function() {
            paste('prediction one step ahead ',var_names[input$y_resp %>% as.numeric],".csv",sep='')
        },
        content = function(file) {
            
            eval=eval_past(Dados_DLM(),FALSE,1)
            
            prediction=eval$pred
            alpha=eval$a
            beta=eval$b
            
            r=alpha
            p=(beta/(beta +1))
            
            icl.pred<-qnbinom((1-0.95)/2, r, p)
            icu.pred<-qnbinom(1-(1-0.95)/2, r, p)
            time=c(1:length(prediction))
            
            table_data=data.frame(time,prediction,icl.pred,icu.pred,r,p)
            
            write.csv2(table_data, file, row.names = FALSE)
        }
    )
    output$save_var_suv <- downloadHandler(
        filename = function() {
            paste(input$var_index,var_names[input$y_resp %>% as.numeric] ,"smooth.csv")
        },
        content = function(file) {
            write.csv2(plot_lat_var(Dados_DLM(),input$var_index,smooth = T,dinamic=FALSE)$table, file, row.names = FALSE)
        }
    )
    output$save_var_filt <- downloadHandler(
        filename = function() {
            paste(input$var_index,var_names[input$y_resp %>% as.numeric],"filtered.csv")
        },
        content = function(file) {
            write.csv2(plot_lat_var(Dados_DLM(),input$var_index,smooth = F,dinamic=FALSE)$table, file, row.names = FALSE)
        }
    )
    
    output$plot_perc_filt <- renderPlotly({
        req(input$Usa_Vel)
        ggplotly(
            plot_lat_var(Dados_DLM(),'Intervenção',smooth = F,tranform_y=function(y){(exp(y)-1)},dinamic=F)$plot+
                labs(title='Efeito percentual (' %>% paste('Intervenção',')'))+
                geom_vline(xintercept=inter_indic,linetype='dashed')+
                scale_x_continuous('Data',breaks=c(c(1:8)*12+1,inter_indic),labels=data_lab[c(c(1:8)*12+1,inter_indic)])+
                scale_y_continuous('$\\Delta y_t\\%$',labels=function(x){(100*x) %>% paste0('%')})
        )
    })
    output$plot_perc_suv <- renderPlotly({
        req(input$Usa_Vel)
        ggplotly(
            plot_lat_var(Dados_DLM(),'Intervenção',smooth = T,tranform_y=function(y){(exp(y)-1)},dinamic=F)$plot+
                labs(title='Efeito percentual (' %>% paste('Intervenção',')'))+
                geom_vline(xintercept=inter_indic,linetype='dashed')+
                scale_x_continuous('Data',breaks=c(c(1:8)*12+1,inter_indic),labels=data_lab[c(c(1:8)*12+1,inter_indic)])+
                scale_y_continuous('$\\Delta y_t\\%$',labels=function(x){(100*x) %>% paste0('%')})
        )
    })
    
    output$slider_Var1 = renderUI({
        plot_data=data_mat[data_mat$Indic_COV==input$Indic_COV &
                   data_mat$resp==var_names[input$y_resp %>% as.numeric] &
                   data_mat$flag_multi==input$Multi_flag &
                   data_mat$Usa_Vel==ifelse(input$Usa_Vel==0,'FALSE','TRUE'),]
        value=plot_data[1,as.numeric(input$var1)]
        print(input$Multi_flag)
        print(plot_data[1,])
        
        sliderInput("val1",
                    names(data_mat)[as.numeric(input$var1)],
                    min = 0.8,
                    max = 1,
                    value = value,
                    step=0.01
                    )
        })
    output$radio_var = renderUI({
        req(input$Usa_Vel)
        options=names(Dados_DLM()$names)
        
        radioButtons("var_index",
                     "Variável para o gráfico de efeitos:",
                     choices=options,
                     selected='Intervenção')
    })
    
    plot_data=eventReactive(c(input$y_resp,input$Multi_flag,input$var1,input$Indic_COV,input$val1,input$Usa_Vel),{
        req(input$val1)
        req(input$var1)
        plot_data=data_mat[data_mat$Indic_COV==input$Indic_COV & data_mat$resp==var_names[input$y_resp %>% as.numeric] & data_mat$flag_multi==input$Multi_flag,]
        plot_data=plot_data[plot_data$Usa_Vel==ifelse(input$Usa_Vel==0,'FALSE','TRUE') & plot_data[,as.numeric(input$var1)]==input$val1,]
        print('plot_data')
        print(plot_data[1,])
        return(plot_data)
        })

    output$distPlot <- renderPlotly({
        req(input$val1)
        req(input$var1)
        name_var1=names(data_mat)[as.numeric(input$var1)]
        name_varplot=names(data_mat)[3:5][-(as.numeric(input$var1)-2)]
        
        plot_data=data_mat[data_mat$Indic_COV==input$Indic_COV,]
        min_Error=min(plot_data$Error)
        max_Error=max(plot_data$Error)
        
        plot_data=plot_data()
        #plot_data$Error=ifelse(plot_data$Error>min(plot_data$Error),plot_data$Error,0)
        x_cut=coords$x_cut
        y_cut=coords$y_cut
        
        ref_Error=min(plot_data[plot_data[[name_varplot[1]]]==x_cut & plot_data[[name_varplot[2]]]==y_cut,8])
        
        x_min=plot_data()[1,][[name_varplot[1]]]
        y_min=plot_data()[1,][[name_varplot[2]]]
        label_name=ifelse(x_min==x_cut & y_min==y_cut,'Mínimo: ','Error: ')
        plot_data$Error_color=plot_data$Error
        
        labels_eixo=list('Delta.M1'='Fator de desconto para o nível',
                         'Delta.COV'='Fator de desconto para a COVID',
                         'Delta.Inter'='Fator de desconto para a intervenção')
        
        ggplotly(
        ggplot(plot_data)+
            geom_tile(aes_string(x=name_varplot[1],y=name_varplot[2],fill='Error_color',value='Error'))+
            geom_text(aes(x=ifelse(x_cut+0.025+0.02>1,x_cut-0.025,x_cut+0.025),y=ifelse(y_cut+0.01+0.02>1,y_cut-0.01,y_cut+0.01),label=paste0(label_name,ref_Error %>% round(2))))+
            geom_point(aes(x=x_cut,y=y_cut))+
            geom_point(aes(x=x_min,y=y_min))+
            scale_fill_gradientn(colours=c('red','green'))+
            geom_hline(yintercept=y_cut,linetype='dashed')+
            geom_vline(xintercept=x_cut,linetype='dashed')+
            #scale_colour_gradientn(values=c(0,min_Error,max_Error),colors=c('black','red','green'))+
            labs(title='Este gráfico é clicável')+
            scale_x_continuous(labels_eixo[[name_varplot[1]]],expand=c(0,0),breaks=c(c(16:20)/20,x_cut))+
            scale_y_continuous(labels_eixo[[name_varplot[2]]],expand=c(0,0),breaks=c(c(16:20)/20,y_cut))+
            theme_bw(),
        source='heatmap'
        )
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

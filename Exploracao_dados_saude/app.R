library(shiny)
library(ggplot2)
library(plotly)
library(tidyr)
options(scipen = 9999)

calcula_max=function(pre_max){
    value=ifelse(log10(max(pre_max))%%1>0.1,10**(floor(log10(max(pre_max)))),10**(floor(log10(max(pre_max)))-1))
    interval_size=((max(pre_max)%/%value)+2)
    
    value=ifelse(interval_size<4,value/2,value)
    interval_size=ifelse(interval_size<4,interval_size*2,interval_size)
    
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

get_data=function(data_name){
    pre_data_1=read.csv2(paste0('data/Ano/',data_name,' 1984-2007.csv'))[1:12,-12]
    pre_data_2=read.csv2(paste0('data/Ano/',data_name,' 2007-2021.csv'))[1:12,c(-1,-16,-17)]
    
    data_inter=cbind(pre_data_1,pre_data_2)
    names(data_inter)=c('Idade',paste0('X',c(1998:2006)),'X2007_old','X2007_new',paste0('X',c(2008:2020)))
    data_inter[2,2:25]=data_inter[1,2:25]+data_inter[2,2:25]
    data_inter=data_inter[-1,]
    data_inter$Idade=c(2,7,12,17,24.5,34.5,44.5,54.5,64.5,74.5,80)
    return(data_inter)
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

varicela_inter=get_data('varicela internacoes')
varicela_obt=get_data('varicela obitos')
gast_inter=get_data('gastroenterite internacoes')
gast_obt=get_data('gastroenterite obitos')

pop=cbind(read.csv2(paste0('data/ano/populacao 1998-1999.csv'))[1:11,],read.csv2(paste0('data/ano/populacao 2000-2020.csv'))[1:11,c(-1)])
names(pop)=c('Idade',paste0('X',c(1998:2020)))
pop$Idade=c(2,7,12,17,24.5,34.5,44.5,54.5,64.5,74.5,80)

ui <- fluidPage(

    titlePanel(fluidRow("Exposição de dados sobre gastroenterite e varicela",align='center',style='padding-bottom: 10%; padding-top: 10%')),

    mainPanel(
        tabsetPanel(
            tabPanel('Gastroenterite',
                     fluidRow(
                         column(1,offset=1,
                                fluidRow(radioButtons('eixo_gast_idade','Eixo X',list('Por ano'=F,'Por idade'=T))),
                                fluidRow(radioButtons('log_scale_gastro','Escala',list('Padrão'=F,'Logarítmica'=T))),
                                fluidRow(radioButtons('act_2007_gastro','2007',list('Valor mais antigo'='old','Valor mais recente'='new','Soma dos valores'='sum','Interpolação'='inter','Ignorar'='null'),selected='sum')),
                                fluidRow(radioButtons('denom_gastro','Denominador para taxa',list('População'=F,'Internações'=T))),
                                align='left',
                                style='padding-top: 10%;'),
                         column(8,
                                tabsetPanel(
                                    tabPanel('População residente',
                                             plotlyOutput("gastro_pop",height='600px', width = '1000px'),align="center"),
                                    tabPanel('Internações',
                                             plotlyOutput("gastro_int",height='600px', width = '1000px'),align="center"),
                                    tabPanel('Óbitos',
                                             plotlyOutput("gastro_obt",height='600px', width = '1000px'),align="center"),
                                    tabPanel('Taxa de óbitos',
                                             plotlyOutput("gastro_tx",height='600px', width = '1000px'),align="center")
                                )
                         )
                     )),
            tabPanel('Varicela',
                     fluidRow(
                         column(1,offset=1,
                                fluidRow(radioButtons('eixo_vari_idade','Eixo X',list('Por ano'=F,'Por idade'=T))),
                                fluidRow(radioButtons('log_scale_vari','Escala',list('Padrão'=F,'Logarítmica'=T))),
                                fluidRow(radioButtons('act_2007_vari','2007',list('Valor mais antigo'='old','Valor mais recente'='new','Soma dos valores'='sum','Interpolação'='inter','Ignorar'='null'),selected='sum')),
                                fluidRow(radioButtons('denom_vari','Denominador para taxa',list('População'=F,'Internações'=T))),
                                align='left',
                                style='padding-top: 10%;'),
                         column(8,
                                tabsetPanel(
                                    tabPanel('População residente',
                                             plotlyOutput("vari_pop",height='600px', width = '1000px'),align="center"),
                                    tabPanel('Internações',
                                             plotlyOutput("vari_int",height='600px', width = '1000px'),align="center"),
                                    tabPanel('Óbitos',
                                             plotlyOutput("vari_obt",height='600px', width = '1000px'),align="center"),
                                    tabPanel('Taxa de internações',
                                             plotOutput("vari_tx",height='600px', width = '1000px'),align="center")
                                )
                         )
                     ))),
            width=12
        )
)

server <- function(input, output) {
    
    ref_pop=eventReactive(input$act_2007_gastro,{
        dados=pivot_longer(pop,c(2:24))
        dados$name=as.factor(substr(dados$name,2,5))
        dados
    })
    gastro_inter=eventReactive(input$act_2007_gastro,{
        dados=pivot_longer(gast_inter,c(2:25))
        if(input$act_2007_gastro=='new'){
            dados=dados[dados$name!='X2007_old',]
            dados$name[dados$name=='X2007_new']='X2007'
            
        }else{
            if(input$act_2007_gastro=='old'){
                dados=dados[dados$name!='X2007_new',]
                dados$name[dados$name=='X2007_old']='X2007'
                
            }else{
                if(input$act_2007_gastro=='sum'){
                    placeholder=dados$value[dados$name=='X2007_old']+dados$value[dados$name=='X2007_new']
                    dados=dados[dados$name!='X2007_old',]
                    dados$name[dados$name=='X2007_new']='X2007'
                    dados$value[dados$name=='X2007']=placeholder
                    
                }else{
                    if(input$act_2007_gastro=='inter'){
                        placeholder=(dados$value[dados$name=='X2006']+dados$value[dados$name=='X2008'])/2
                        dados=dados[dados$name!='X2007_old',]
                        dados$name[dados$name=='X2007_new']='X2007'
                        dados$value[dados$name=='X2007']=placeholder
                        
                    }else{
                        dados=dados[dados$name!='X2007_old',]
                        dados$name[dados$name=='X2007_new']='X2007'
                        dados$value[dados$name=='X2007']=NA
                    }
                }
            }
        }
        dados$name=as.factor(substr(dados$name,2,5))
        dados
    })
    gastro_obt=eventReactive(input$act_2007_gastro,{
        dados=pivot_longer(gast_obt,c(2:25))
        if(input$act_2007_gastro=='new'){
            dados=dados[dados$name!='X2007_old',]
            dados$name[dados$name=='X2007_new']='X2007'
            
        }else{
            if(input$act_2007_gastro=='old'){
                dados=dados[dados$name!='X2007_new',]
                dados$name[dados$name=='X2007_old']='X2007'
                
            }else{
                if(input$act_2007_gastro=='sum'){
                    placeholder=dados$value[dados$name=='X2007_old']+dados$value[dados$name=='X2007_new']
                    dados=dados[dados$name!='X2007_old',]
                    dados$name[dados$name=='X2007_new']='X2007'
                    dados$value[dados$name=='X2007']=placeholder
                    
                }else{
                    if(input$act_2007_gastro=='inter'){
                        placeholder=(dados$value[dados$name=='X2006']+dados$value[dados$name=='X2008'])/2
                        dados=dados[dados$name!='X2007_old',]
                        dados$name[dados$name=='X2007_new']='X2007'
                        dados$value[dados$name=='X2007']=placeholder
                        
                    }else{
                        dados=dados[dados$name!='X2007_old',]
                        dados$name[dados$name=='X2007_new']='X2007'
                        dados$value[dados$name=='X2007']=NA
                    }
                }
            }
        }
        dados$name=as.factor(substr(dados$name,2,5))
        dados
    })
    
    vari_inter=eventReactive(input$act_2007_vari,{
        dados=pivot_longer(varicela_inter,c(2:25))
        if(input$act_2007_vari=='new'){
            dados=dados[dados$name!='X2007_old',]
            dados$name[dados$name=='X2007_new']='X2007'
            
        }else{
            if(input$act_2007_vari=='old'){
                dados=dados[dados$name!='X2007_new',]
                dados$name[dados$name=='X2007_old']='X2007'
                
            }else{
                if(input$act_2007_vari=='sum'){
                    placeholder=dados$value[dados$name=='X2007_old']+dados$value[dados$name=='X2007_new']
                    dados=dados[dados$name!='X2007_old',]
                    dados$name[dados$name=='X2007_new']='X2007'
                    dados$value[dados$name=='X2007']=placeholder
                    
                }else{
                    if(input$act_2007_vari=='inter'){
                        placeholder=(dados$value[dados$name=='X2006']+dados$value[dados$name=='X2008'])/2
                        dados=dados[dados$name!='X2007_old',]
                        dados$name[dados$name=='X2007_new']='X2007'
                        dados$value[dados$name=='X2007']=placeholder
                        
                    }else{
                        dados=dados[dados$name!='X2007_old',]
                        dados$name[dados$name=='X2007_new']='X2007'
                        dados$value[dados$name=='X2007']=NA
                    }
                }
            }
        }
        dados$name=as.factor(substr(dados$name,2,5))
        dados
    })
    vari_obt=eventReactive(input$act_2007_vari,{
        dados=pivot_longer(varicela_obt,c(2:25))
        if(input$act_2007_vari=='new'){
            dados=dados[dados$name!='X2007_old',]
            dados$name[dados$name=='X2007_new']='X2007'
            
        }else{
            if(input$act_2007_vari=='old'){
                dados=dados[dados$name!='X2007_new',]
                dados$name[dados$name=='X2007_old']='X2007'
                
            }else{
                if(input$act_2007_vari=='sum'){
                    placeholder=dados$value[dados$name=='X2007_old']+dados$value[dados$name=='X2007_new']
                    dados=dados[dados$name!='X2007_old',]
                    dados$name[dados$name=='X2007_new']='X2007'
                    dados$value[dados$name=='X2007']=placeholder
                    
                }else{
                    if(input$act_2007_vari=='inter'){
                        placeholder=(dados$value[dados$name=='X2006']+dados$value[dados$name=='X2008'])/2
                        dados=dados[dados$name!='X2007_old',]
                        dados$name[dados$name=='X2007_new']='X2007'
                        dados$value[dados$name=='X2007']=placeholder
                        
                    }else{
                        dados=dados[dados$name!='X2007_old',]
                        dados$name[dados$name=='X2007_new']='X2007'
                        dados$value[dados$name=='X2007']=NA
                    }
                }
            }
        }
        dados$name=as.factor(substr(dados$name,2,5))
        dados
    })
    
    output$gastro_pop <- renderPlotly({
        dados=ref_pop()
        
        maximos=calcula_max(dados$value)
        value=maximos[[1]]
        interval_size=maximos[[2]]
        max_value=maximos[[3]]
        
        if(input$eixo_gast_idade){
            linha=geom_line(aes(x=Idade,y=value,color=name))
            scale_x=scale_x_continuous('Faixa etária',breaks=unique(dados$Idade),labels = labels_idade,expand=c(0,0))
            scale_color=scale_color_hue('Ano')
        }else{
            dados$Idade=as.factor(dados$Idade)
            levels(dados$Idade)=labels_idade
            dados$name=as.numeric(as.character(dados$name))
            linha=geom_line(aes(x=name,y=value,color=Idade))
            scale_x=scale_x_continuous('Ano',breaks=unique(dados$name),expand=c(0,0))
            scale_color=scale_color_hue('Faixa etária',labels = labels_idade)
            
        }
        if(input$log_scale_gastro){
            scale_y=scale_y_log10('População residente')
        }else{
            scale_y=scale_y_continuous('População residente',
                                       breaks=c(0:interval_size)*value,
                                       limits=c(0,max_value),expand=c(0,0),
                                       labels=set_format(abs(c(0:interval_size))*value))
        }
        turn_dynamic(
            ggplot(dados)+
                linha+
                scale_x+
                scale_y+
                theme_bw()+
                labs(title='Quantidade de pessoas residentes no Brasil')+
                theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
        )
    })
    output$gastro_int <- renderPlotly({
        dados=gastro_inter()
        
        maximos=calcula_max(dados$value)
        value=maximos[[1]]
        interval_size=maximos[[2]]
        max_value=maximos[[3]]
        
        if(input$eixo_gast_idade){
            linha=geom_line(aes(x=Idade,y=value,color=name,shape=name))
            ponto=geom_point(aes(x=Idade,y=value,color=name,shape=name))
            scale_x=scale_x_continuous('Faixa etária',breaks=unique(dados$Idade),labels = labels_idade,expand=c(0,0))
            scale_color=scale_color_hue('Ano')
        }else{
            dados$Idade=as.factor(dados$Idade)
            levels(dados$Idade)=labels_idade
            dados$name=as.numeric(as.character(dados$name))
            linha=geom_line(aes(x=name,y=value,color=Idade,shape=Idade))
            ponto=geom_point(aes(x=Idade,y=value,color=Idade,shape=Idade))
            scale_x=scale_x_continuous('Ano',breaks=unique(dados$name),expand=c(0,0))
            scale_color=scale_color_hue('Faixa etária',labels = labels_idade)
                
        }
        if(input$log_scale_gastro){
            scale_y=scale_y_log10('Internações',expand=c(0,0))
        }else{
            scale_y=scale_y_continuous('Internações',
                                       breaks=c(0:interval_size)*value,
                                       limits=c(0,max_value),expand=c(0,0),
                                       labels=set_format(abs(c(0:interval_size))*value))
        }
        turn_dynamic(
            ggplot(dados)+
                linha+
                ponto+
                scale_x+
                scale_y+
                scale_color+
                theme_bw()+
                labs(title='Quantidade de internações por gastroenterite no Brasil.')+
                theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
        )
    })
    output$gastro_obt <- renderPlotly({
        dados=gastro_obt()
        
        maximos=calcula_max(dados$value)
        value=maximos[[1]]
        interval_size=maximos[[2]]
        max_value=maximos[[3]]
        
        if(input$eixo_gast_idade){
            linha=geom_line(aes(x=Idade,y=value,color=name))
            scale_x=scale_x_continuous('Faixa etária',breaks=unique(dados$Idade),labels = labels_idade,expand=c(0,0))
            scale_color=scale_color_hue('Ano')
        }else{
            dados$Idade=as.factor(dados$Idade)
            levels(dados$Idade)=labels_idade
            dados$name=as.numeric(as.character(dados$name))
            linha=geom_line(aes(x=name,y=value,color=Idade))
            scale_x=scale_x_continuous('Ano',breaks=unique(dados$name),expand=c(0,0))
            scale_color=scale_color_hue('Faixa etária',labels = labels_idade)
            
        }
        if(input$log_scale_gastro){
            scale_y=scale_y_log10('Óbitos',expand=c(0,0))
        }else{
            scale_y=scale_y_continuous('Óbitos',
                                       breaks=c(0:interval_size)*value,
                                       limits=c(0,max_value),expand=c(0,0),
                                       labels=set_format(abs(c(0:interval_size))*value))
        }
        turn_dynamic(
            ggplot(dados)+
                linha+
                scale_x+
                scale_y+
                scale_color+
                theme_bw()+
                labs(title='Quantidade de óbitos por gastroenterite no Brasil')+
                theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
        )
    })
    output$gastro_tx <- renderPlotly({
        if(input$denom_gastro){
            dados_exp=gastro_inter()
        }else{
            dados_exp=ref_pop()
        }
        dados_obt=gastro_obt()
        
        dados=dados_exp
        dados$value=ifelse(dados_exp$value==0,NA,dados_obt$value/dados_exp$value)
        
        maximos=calcula_max(dados$value)
        value=maximos[[1]]
        interval_size=maximos[[2]]
        max_value=maximos[[3]]
        
        if(input$eixo_gast_idade){
            linha=geom_line(aes(x=Idade,y=value,color=name))
            scale_x=scale_x_continuous('Faixa etária',breaks=unique(dados$Idade),labels = labels_idade,expand=c(0,0))
            scale_color=scale_color_hue('Ano')
        }else{
            dados$Idade=as.factor(dados$Idade)
            levels(dados$Idade)=labels_idade
            dados$name=as.numeric(as.character(dados$name))
            linha=geom_line(aes(x=name,y=value,color=Idade))
            scale_x=scale_x_continuous('Ano',breaks=unique(dados$name),expand=c(0,0))
            scale_color=scale_color_hue('Faixa etária',labels = labels_idade)
            
        }
        if(input$log_scale_gastro){
            scale_y=scale_y_log10('Taxa de mortalidade',
                                  breaks=10**c(-8:0),
                                  limits=10**c(-8,0),expand=c(0,0))
        }else{
            scale_y=scale_y_continuous('Taxa de mortalidade',expand=c(0,0))
        }
        turn_dynamic(
            ggplot(dados)+
                linha+
                scale_x+
                scale_y+
                scale_color+
                theme_bw()+
                labs(title=ifelse(input$denom_gastro,'Taxa de óbitos por gastroenterite (a cada internação)','Taxa de óbitos por gastroenterite (a cada residente)'))+
                theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
        )
    })
    
    output$vari_pop <- renderPlotly({
        dados=ref_pop()
        
        maximos=calcula_max(dados$value)
        value=maximos[[1]]
        interval_size=maximos[[2]]
        max_value=maximos[[3]]
        
        if(input$eixo_vari_idade){
            linha=geom_line(aes(x=Idade,y=value,color=name))
            scale_x=scale_x_continuous('Faixa etária',breaks=unique(dados$Idade),labels = labels_idade,expand=c(0,0))
            scale_color=scale_color_hue('Ano')
        }else{
            dados$Idade=as.factor(dados$Idade)
            levels(dados$Idade)=labels_idade
            dados$name=as.numeric(as.character(dados$name))
            linha=geom_line(aes(x=name,y=value,color=Idade))
            scale_x=scale_x_continuous('Ano',breaks=unique(dados$name),expand=c(0,0))
            scale_color=scale_color_hue('Faixa etária',labels = labels_idade)
            
        }
        if(input$log_scale_vari){
            scale_y=scale_y_log10('População residente',expand=c(0,0))
        }else{
            scale_y=scale_y_continuous('População residente',
                                       breaks=c(0:interval_size)*value,
                                       limits=c(0,max_value),expand=c(0,0),
                                       labels=set_format(abs(c(0:interval_size))*value))
        }
        turn_dynamic(
            ggplot(dados)+
                linha+
                scale_x+
                scale_y+
                scale_color+
                theme_bw()+
                labs(title='Quantidade de pessoas residentes no Brasil')+
                theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
        )
    })
    output$vari_int <- renderPlotly({
        dados=vari_inter()
        
        maximos=calcula_max(dados$value)
        value=maximos[[1]]
        interval_size=maximos[[2]]
        max_value=maximos[[3]]
        
        if(input$eixo_vari_idade){
            linha=geom_line(aes(x=Idade,y=value,color=name))
            ponto=geom_point(aes(x=Idade,y=value,color=name,shape=name))
            scale_x=scale_x_continuous('Faixa etária',breaks=unique(dados$Idade),labels = labels_idade,expand=c(0,0))
            scale_color=scale_color_hue('Ano')
        }else{
            dados$Idade=as.factor(dados$Idade)
            levels(dados$Idade)=labels_idade
            dados$name=as.numeric(as.character(dados$name))
            linha=geom_line(aes(x=name,y=value,color=Idade))
            ponto=geom_point(aes(x=name,y=value,color=Idade,shape=Idade))
            scale_x=scale_x_continuous('Ano',breaks=unique(dados$name),expand=c(0,0))
            scale_color=scale_color_hue('Faixa etária',labels = labels_idade)
            
        }
        if(input$log_scale_vari){
            scale_y=scale_y_log10('Internações',expand=c(0,0))
        }else{
            scale_y=scale_y_continuous('Internações',
                                       breaks=c(0:interval_size)*value,
                                       limits=c(0,max_value),expand=c(0,0),
                                       labels=set_format(abs(c(0:interval_size))*value))
        }
        turn_dynamic(
            ggplot(dados)+
                linha+
                ponto+
                scale_x+
                scale_y+
                scale_color+
                theme_bw()+
                labs(title='Quantidade de internações por varicela no Brasil')+
                theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
        )
    })
    output$vari_obt <- renderPlotly({
        dados=vari_obt()
        
        maximos=calcula_max(dados$value)
        value=maximos[[1]]
        interval_size=maximos[[2]]
        max_value=maximos[[3]]
        
        if(input$eixo_gast_idade){
            linha=geom_line(aes(x=Idade,y=value,color=name))
            scale_x=scale_x_continuous('Faixa etária',breaks=unique(dados$Idade),labels = labels_idade,expand=c(0,0))
            scale_color=scale_color_hue('Ano')
        }else{
            dados$Idade=as.factor(dados$Idade)
            levels(dados$Idade)=labels_idade
            dados$name=as.numeric(as.character(dados$name))
            linha=geom_line(aes(x=name,y=value,color=Idade))
            scale_x=scale_x_continuous('Ano',breaks=unique(dados$name),expand=c(0,0))
            scale_color=scale_color_hue('Faixa etária',labels = labels_idade)
            
        }
        if(input$log_scale_vari){
            scale_y=scale_y_log10('Óbitos',expand=c(0,0))
        }else{
            scale_y=scale_y_continuous('Óbitos',
                                       breaks=c(0:interval_size)*value,
                                       limits=c(0,max_value),expand=c(0,0),
                                       labels=set_format(abs(c(0:interval_size))*value))
        }
        turn_dynamic(
            ggplot(dados)+
                linha+
                scale_x+
                scale_y+
                scale_color+
                theme_bw()+
                labs(title='Quantidade de óbitos por varicela no Brasil')+
                theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
        )
    })
    output$vari_tx <- renderPlot({
        if(input$denom_vari){
            dados_exp=vari_inter()
        }else{
            dados_exp=ref_pop()
        }
        dados_obt=vari_inter()
        
        dados=dados_exp
        dados$value=ifelse(dados_exp$value==0,NA,dados_obt$value/dados_exp$value)
        
        maximos=calcula_max(dados$value)
        value=maximos[[1]]
        interval_size=maximos[[2]]
        max_value=maximos[[3]]
        
        if(input$eixo_vari_idade){
            linha=geom_line(aes(x=Idade,y=100000*value/12,color=name,shape=name))
            ponto=geom_point(aes(x=Idade,y=100000*value/12,color=name,shape=name))
            scale_x=scale_x_continuous('Faixa etária',breaks=unique(dados$Idade),labels = labels_idade,expand=c(0,0))
            scale_color=scale_color_hue('Ano')
        }else{
            dados$Idade=as.factor(dados$Idade)
            levels(dados$Idade)=labels_idade
            dados$name=as.numeric(as.character(dados$name))
            linha=geom_line(aes(x=name,y=100000*value/12))
            ponto=geom_point(aes(x=name,y=100000*value/12))
            scale_x=scale_x_continuous('Ano',breaks=c(0:10)*2+2000,expand=c(0,0))
            scale_color=scale_color_hue('Faixa etária',labels = labels_idade)
            
        }
        if(input$log_scale_vari){
            scale_y=scale_y_log10('Taxa de internações',
                                  breaks=10**c(-8:0),
                                  limits=10**c(-8,0),expand=c(0,0))
        }else{
            scale_y=scale_y_continuous('Internações a cada 100.000 residente',expand=c(0,0),limits=c(0,8))
        }
            ggplot(dados)+
                linha+
                ponto+
                scale_x+
                scale_y+
                scale_color+
                theme_bw()+
                scale_shape_manual(values=c(0:10))+
                facet_wrap(~Idade)+
                geom_vline(xintercept=2008,linetype='dashed')+
                labs(title=ifelse(input$denom_gastro,'Média mensal de internações por varicela','Média mensal de internações por varicela'))+
                theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    })
}

shinyApp(ui = ui, server = server)

library(shiny)
library(ggplot2)
library(plotly)
library(stringr)
library(tidyr)
options(scipen = 9999,encoding = 'UTF-8')
print(sessionInfo())

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
    pre_data[,1]=as.character(pre_data[,1])
    pre_data[1,-1]=pre_data[1,-1]+pre_data[2,-1]
    pre_data[1,1]='00 a 04 anos'
    pre_data[,1]=as.factor(pre_data[,1])
    print(pre_data)
    data=pivot_longer(pre_data[-2,],-1)
    names(data)=c('Faixa_etaria','Data','Valor')
    data$Faixa_etaria=factor(data$Faixa_etaria)
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

data=get_data('varicela internacoes')
data=merge(data,get_data('varicela obitos'),by=c('Faixa_etaria','Data'))
data$iden=1
data$mes=data$Data
data$Data=substr(data$Data,1,4)


pop=cbind(read.csv2(paste0('data/Ano/populacao 2000-2020.csv'))[1:11,-c(3:8)])
pop=pivot_longer(pop,-1)
names(pop)=c('Faixa_etaria','Data','Valor')
pop$Data=substr(pop$Data,2,5)
data=merge(data,pop,by=c('Faixa_etaria','Data'))
data$Data=data$mes
ref_data=data[,-6]
names(ref_data)=c('Faixa_etaria','Data','Internações','Óbitos','1','Residentes')
#pop$Idade=c(2,7,12,17,24.5,34.5,44.5,54.5,64.5,74.5,80)

ui <- fluidPage(
    
    titlePanel(fluidRow("Exposição de dados sobre varicela",align='center',style='padding-bottom: 10%; padding-top: 10%')),
    
    fluidRow(
        column(1,offset=1,
               fluidRow(radioButtons('eixo_idade','Eixo X',list('Por ano'=F,'Por idade'=T))),
               fluidRow(radioButtons('group','Agrupamento',list('DATASUS'=0,
                                                                'Reunião 30/07'=1,
                                                                'Reunião 25/08\n(5 grupos)'=2,
                                                                'Reunião 25/08\n(4 grupos)'=3))),
               fluidRow(radioButtons('log_scale','Escala',list('Padrão'=F,'Logarítmica'=T))),
               fluidRow(radioButtons('numerador','Numerador',list('Internações'='Internações','Óbitos'='Óbitos','População residente'='Residentes','1'='1'),selected='Internações')),
               fluidRow(radioButtons('denominador','Denominador',list('Internações'='Internações','Óbitos'='Óbitos','População residente'='Residentes','1'='1'),selected='Residentes')),
               align='left',
               style='padding-top: 5%;'),
        column(8,
               plotlyOutput("absol",height='600px', width = '1000px')),
        align='center'
    )
)

server <- function(input, output) {
    
    output$absol <- renderPlotly({
        data=ref_data
        if(input$group<1){
            values_idade=c(2,7,12,17,24.5,34.5,44.5,54.5,64.5,74.5,80)
        }else{if(input$group<2){
            data$Faixa_etaria=as.character(data$Faixa_etaria)
            data$Faixa_etaria[data$Faixa_etaria %in% c("15 a 19 anos","20 a 29 anos","30 a 39 anos")]='15 a 39 anos'
            data$Faixa_etaria[data$Faixa_etaria %in% c("40 a 49 anos","50 a 59 anos")]='40 a 59 anos'
            data$Faixa_etaria[data$Faixa_etaria %in% c("60 a 69 anos","70 a 79 anos","80 anos e mais")]='60 anos e mais'
            data$Faixa_etaria=as.factor(data$Faixa_etaria)
            data=aggregate(.~Faixa_etaria+Data,data,sum)
            data[['1']]=1
            values_idade=c(2,7,12,27,44.5,60)
        }else{if(input$group<3){
            data$Faixa_etaria=as.character(data$Faixa_etaria)
            data$Faixa_etaria[data$Faixa_etaria %in% c('15 a 19 anos',"20 a 29 anos","30 a 39 anos",'40 a 49 anos')]='15 a 49 anos'
            data$Faixa_etaria[data$Faixa_etaria %in% c("50 a 59 anos","60 a 69 anos",'70 a 79 anos','80 anos e mais')]='50 anos e mais'
            data$Faixa_etaria=as.factor(data$Faixa_etaria)
            data=aggregate(.~Faixa_etaria+Data,data,sum)
            data[['1']]=1
            values_idade=c(2,7,12,32,50)
        }else{
            data$Faixa_etaria=as.character(data$Faixa_etaria)
            data$Faixa_etaria[data$Faixa_etaria %in% c("05 a 09 anos","10 a 14 anos")]='05 a 14 anos'
            data$Faixa_etaria[data$Faixa_etaria %in% c('15 a 19 anos',"20 a 29 anos","30 a 39 anos",'40 a 49 anos')]='15 a 49 anos'
            data$Faixa_etaria[data$Faixa_etaria %in% c("50 a 59 anos","60 a 69 anos",'70 a 79 anos','80 anos e mais')]='50 anos e mais'
            data$Faixa_etaria=as.factor(data$Faixa_etaria)
            data=aggregate(.~Faixa_etaria+Data,data,sum)
            data[['1']]=1
            values_idade=c(2,9.5,32,50)
        }
        }
        }
    data$valor=ifelse(data[[input$denominador]]==0,0,data[[input$numerador]]/data[[input$denominador]])
    
    if(input$eixo_idade){
        levels(data$Faixa_etaria)=values_idade
        data$Faixa_etaria=as.numeric(as.character(data$Faixa_etaria))
        linha=geom_line(aes(x=Faixa_etaria,y=valor,color=Data))
        scale_x=scale_x_continuous('Faixa etária',breaks=unique(data$Faixa_etaria),labels = labels_idade,expand=c(0,0))
        scale_color=scale_color_hue('Ano')
    }else{
        linha=geom_line(aes(x=Data,y=valor))#,color=Faixa_etaria))
        scale_x=scale_x_date('',date_breaks='1 years',expand=c(0,0),date_labels = "%Y")
        scale_color=scale_color_hue('Faixa etária',labels = labels_idade)
        
    }
    if(input$log_scale){
        scale_y=scale_y_log10(input$numerador)
    }else{
        if(input$denominador!='1'){
            label=paste(input$numerador,'(a cada 100,000 ',str_to_lower(input$denominador),')')
            data$valor=data$valor*100000
        }else{
            label=paste(input$numerador)
        }
        maximos=calcula_max(data$valor)
        value=maximos[[1]]
        interval_size=maximos[[2]]
        max_value=maximos[[3]]
        scale_y=scale_y_continuous(label,
                                   breaks=c(0:interval_size)*value,
                                   limits=c(0,max_value),expand=c(0,0),
                                   labels=set_format(abs(c(0:interval_size))*value))
    }
    turn_dynamic(
        ggplot(data)+
            linha+
            scale_x+
            scale_y+
            theme_bw()+
            #labs(title='Quantidade de pessoas residentes no Brasil')+
            facet_wrap(~Faixa_etaria)+
            theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    )
    })
}

shinyApp(ui = ui, server = server)

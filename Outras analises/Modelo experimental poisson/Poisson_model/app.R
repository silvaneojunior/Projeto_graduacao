#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(latex2exp)

data=read.csv2(paste0('data/Ano/','varicela internacoes',' 2007-2021.csv'))[1:12,c(-16,-17)]
n=5000
beta_0=1

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Ajuste do modelo temporal Poisson"),
    withMathJax(),

    # Sidebar with a slider input for number of bins 
    fluidRow(
        sidebarLayout(
        sidebarPanel(
            sliderInput("Indice",
                        "Indice da faixa etária:",
                        min = 1,
                        max = 12,
                        value = 1),
            sliderInput("beta",
                        "\\( \\beta \\):",
                        min = 0.0,
                        max = 2.0,
                        value = 1.0,
                        step = 0.01),
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
           fluidRow(plotOutput("distPlot",width = "800px",height = "600px",)),
           fluidRow(plotOutput("distPlot2",width = "800px",height = "600px",)),
           fluidRow(plotOutput("distPlot3",width = "800px",height = "600px",))
        )
        )
    ),
    fluidRow(
        sidebarPanel(
            sliderInput("Indice_ano",
                        "Indice do ano:",
                        min = 2008,
                        max = 2020,
                        value = 2008),
            sliderInput("n",
                        "n:",
                        min = 1,
                        max = 100,
                        value = 10),
            sliderInput("b",
                        "Burn-in:",
                        min = 1,
                        max = 1000,
                        value = 10)
            
        ),
    mainPanel(
        fluidRow(plotOutput("distPlot4",width = "800px",height = "600px",)),
        fluidRow(plotOutput("distPlot5",width = "800px",height = "600px",)),
        fluidRow(plotOutput("distPlot6",width = "800px",height = "600px",)),
        fluidRow(plotOutput("distPlot7",width = "800px",height = "600px",)),
        fluidRow(plotOutput("distPlot8",width = "800px",height = "600px",))
    )
    )
)

server <- function(input, output) {
    values=reactiveValues(y=0,samples=0,alpha=0,prorp=1)
        
    observeEvent(c(input$Indice,input$beta),{
        values$y=as.numeric(as.matrix(data)[input$Indice,-c(1,2)])
        
        values$samples=matrix(1/13,n+1,13)
        beta_0=input$beta
        values$prorp=c(1)
        values$alpha=rgamma(n+1,sum(values$y[7:13]),1)
        for(i in 1:n+1){
            if(i%%100==0){
                print(i)
            }
            placeholder=c()
            prorp=rgamma(1,sum(values$samples[i-1,1:12]),sum(values$samples[i-1,2:13]))
            alpha=values$alpha[i]/sum(values$samples[i-1,7:13])
            values$prorp=c(values$prorp,prorp)
            values$alpha[i]=alpha
            last_guy=rgamma(1,shape=values$y[1],rate=1)
            placeholder=c(placeholder,last_guy)
            for(j in c(2:6)){
                last_guy=rgamma(1,shape=beta_0*last_guy+values$y[j],rate=prorp*beta_0+1)
                placeholder=c(placeholder,last_guy)
            }
            for(j in c(7:13)){
                last_guy=rgamma(1,shape=beta_0*last_guy+values$y[j],rate=prorp*beta_0+alpha)
                placeholder=c(placeholder,last_guy)
            }
            values$samples[i,1:13]=placeholder
        }
    })
    
    output$distPlot <- renderPlot({
        y=values$y
        samples=values$samples[c(input$b:floor(n/input$n))*input$n,]
        alpha=values$alpha[c(input$b:floor(n/input$n))*input$n]
        
        samples[,7:13]=samples[,7:13]*alpha
        
        years=2007+c(1:length(y))
        
        q025=c()
        q975=c()
        pred025=c()
        pred975=c()
        for(i in c(1:13)){
            q025=c(q025,quantile(samples[,i],0.025))
            q975=c(q975,quantile(samples[,i],0.975))
            pred_sample=rpois(n,samples[,i])
            pred025=c(pred025,quantile(pred_sample,0.025))
            pred975=c(pred975,quantile(pred_sample,0.975))
        }
        ggplot()+
            geom_ribbon(aes(x=years,ymin=pred025,ymax=pred975,fill='I.C. para $y_i$'),alpha=0.5)+
            geom_ribbon(aes(x=years,ymin=q025,ymax=q975,fill='I.C. para $\\lambda_i$'),alpha=0.25)+
            geom_line(aes(x=years,y=colMeans(samples)))+
            geom_point(aes(x=years,y=colMeans(samples)))+
            geom_point(aes(x=years,y=y),color='red')+
            scale_x_continuous('Ano',breaks=years)+
            scale_y_continuous('Quantidade de internações')+
            scale_fill_manual('',values=c('#0000ff','#aaaaff'),labels=TeX)+
            theme_bw()+theme(legend.position = 'bottom')
    })
    
    output$distPlot2 <- renderPlot({
        hist(values$alpha[c(input$b:floor(n/input$n))*input$n])
    })
    
    output$distPlot3 <- renderPlot({
        y=values$y
        samples=values$samples[c(input$b:floor(n/input$n))*input$n,]
        
        years=2007+c(1:length(y))
        
        q025=c()
        q975=c()
        pred025=c()
        pred975=c()
        for(i in c(1:13)){
            q025=c(q025,quantile(samples[,i],0.025))
            q975=c(q975,quantile(samples[,i],0.975))
        }
        ggplot()+
            geom_ribbon(aes(x=years,ymin=q025,ymax=q975,fill='I.C. para $\\lambda_i$'),alpha=0.25)+
            geom_line(aes(x=years,y=colMeans(samples)))+
            geom_point(aes(x=years,y=colMeans(samples)))+
            scale_x_continuous('Ano',breaks=years)+
            scale_y_continuous('Quantidade de internações')+
            scale_fill_manual('',values=c('#0000ff'),labels=TeX)+
            theme_bw()+theme(legend.position = 'bottom')
    })
    
    output$distPlot4 <- renderPlot({
        acf(values$samples[c(input$b:floor(n/input$n))*input$n,input$Indice_ano-2007])
    })
    output$distPlot5 <- renderPlot({
        print('alpha')
        print(summary(values$alpha))
        acf(values$alpha[c(input$b:floor(n/input$n))*input$n])
    })
    output$distPlot6 <- renderPlot({
        plot(values$alpha[c(input$b:floor(n/input$n))*input$n])
    })
    output$distPlot7 <- renderPlot({
        print('prorp')
        print(summary(values$prorp))
        acf(values$prorp[c(input$b:floor(n/input$n))*input$n])
    })
    output$distPlot8 <- renderPlot({
        plot(values$prorp[c(input$b:floor(n/input$n))*input$n])
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

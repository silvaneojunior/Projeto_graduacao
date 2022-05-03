library(shiny)
library(shinyjs)

final_ui=fluidPage(
  tags$head(tags$style(".modal-dialog{ width:1000px}")),
  tags$head(tags$style(".modal-body{ min-height:800px}")),
  withMathJax(),
  uiOutput('main_ui'))

final_server=function(input, output) {
  output$main_ui=renderUI({
    if(ifelse(is.null(input$varicela),0,input$varicela)>0){
      removeModal()
      source('grid_app_varicela.R',encoding='UTF-8')
      server(input, output)
      return(ui)}else{
    if(ifelse(is.null(input$wania),0,input$wania)>0){
      removeModal()
      source('grid_app_Wania.R',encoding='UTF-8')
      server(input, output)
      return(ui)}else{
      
      return(
          showModal(
            modalDialog(fluidRow(
                tags$style(HTML("
                  .btn {
                    height: 800px;
                    width: 450px;
                    border: 10px transparent gray;
                    }

                    ")),
              actionButton('varicela','Dados de varicela',style='font-size:20pt'),
              actionButton('wania','Dados da intervenção da Wania',style='font-size:20pt'),
              align='center'),
              footer=NULL,
              size='l'
              )
            )
          )
      }}
  })
}

shinyApp(ui = final_ui, server = final_server)
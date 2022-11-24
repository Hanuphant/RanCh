library(httr)
library(jsonlite)
library(dplyr)
library(stringr)
library(tidyr)
library(modules)
api <- use("ChIPRmod/ENCODEAPI.R")

library(shiny)
library(shinythemes)

# Define UI for application that takes in RNA-Seq featurecounts and ChIP-Seq narrowPeak files
ui <- fluidPage(theme = shinytheme("cerulean"),
                titlePanel("ChIPR"),
                tabPanel("RNA-Seq",
                         sidebarLayout(
                             sidebarPanel(
                                 fileInput("rna_seq", "Select RNA-Seq featurecounts file",
                                           accept = c("text/csv", "text/comma-separated-values,text/plain",
                                                      ".csv")),
                                 textInput("tcp_name", "Enter transcription factor name"),
                                    selectInput("biosample_classification", "Select biosample classification",
                                                choices = api$initial_search(input$tcp_name)),
                                    selectInput("biosample", "Select biosample",
                                                choices = api$search_based_biosample_classifications(tcp_name = input$tcp_name, biosample_classifications = input$biosample_classification)),

                             ),
                         ),

                                mainPanel(
                                    textOutput(outputId = "chip_seq_biosample_classification"),
                                )
                ),


server <- function(input, output) {

    output$chip_seq_biosample_classification <- renderText({
        paste("The input transcription factor is", input$tcp_name,'and the biosample classification is', input$biosample_classification)
    })

shinyApp(ui = ui, server = server)

library(httr)
library(jsonlite)
library(dplyr)
library(stringr)
library(tidyr)

# Get intial search results for the transcription factor
initial_search <- function(tcp_name){

    # This function takes a transcription factor name and returns the initial search results
    #
    # Args:
    #     tcp_name: The name of the transcription factor to search for
    # Returns:
    #     biosample_classifications: A list of biosample classifications for the transcription factor


    # Set up the url for the API call
    experiment_type <- 'TF+ChIP-seq'
    target_label <- tcp_name
    url <- paste0('https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false&assay_title=',experiment_type,'&target.label=',target_label,'&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&files.file_type=bed+narrowPeak&format=json')

    # Make the API call
    response <- GET(url)
    if (response$status_code != 200) {
    return(paste0("Error: ", response$status_code))
    }

    # Parse the response
    response <- fromJSON(content(response, "text", encoding = "UTF-8"))

    # Get the biosample classifications
    biosample_classifications <- response$facets$terms[[9]]$key

    return(unlist(biosample_classifications))
}


# Now with the biosample_classifications the user will have the option to select which biological sample does the RNA-Seq belong to.

#TODO: User input for biological sample here

search_based_biosample_classifications <- function(tcp_name, biosample_classification){

    # Function to search for the TF-ChipSeq experiments based on the biosample classification
    #
    # Args:
    #     tcp_name: The name of the transcription factor to search for
    #     biosample_classification: The biosample classification to search for
    # Returns:
    #     biosample: Returns the biosample cell line info


    # Set up the url for the API call
    experiment_type <- 'TF+ChIP-seq'
    target_label <- tcp_name
    biosample_classification <- str_replace_all(biosample_classification, ' ', '+')

    url <- paste0('https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false&assay_title=',experiment_type,'&target.label=',target_label,'&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&files.file_type=bed+narrowPeak&format=json&biosample_ontology.classification=',biosample_classification)

    # Make the API call
    response <- GET(url)
    if (response$status_code != 200) {
    return(paste0("Error: ", response$status_code))
    }

    # Parse the response
    response <- fromJSON(content(response, "text", encoding = "UTF-8"))
    biosamples <- response$facets$terms[[10]]$key

    return(unlist(biosamples))
}

#TODO: User input for biological sample here

search_based_biosample <- function(biosample_classification, tcp_name, biosample, jobid){

    # Function to search for the TF-ChipSeq experiments based on the biosample classification
    #
    # Args:
    #     tcp_name: The name of the transcription factor to search for
    #     biosample_classification: The biosample classification to search for
    #     biosample: The biosample cell line to search for
    # Returns:
    #     experiment_accession: Returns the experiment accession number

    target_label <- tcp_name
    biosample_classification <- str_replace_all(biosample_classification, ' ', '+')
    biosample <- str_replace_all(biosample, ' ', '+')

    url <- paste0('https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false&assay_title=TF+ChIP-seq&target.label=',target_label,'&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&files.file_type=bed+narrowPeak&biosample_ontology.classification=',biosample_classification,'&biosample_ontology.term_name=',biosample,'&format=json')
    # print(url)
    # Make the API call
    response <- GET(url)
    if (response$status_code != 200) {
        return(paste0("Error: ", response$status_code))
    }

    # Parse the response
    response <- fromJSON(content(response, "text", encoding = "UTF-8"))
    experiment_accession <- response$`@graph`$accession[1]

    # Response from the experiment accession number
    download_narrowpeak_files(accession = experiment_accession, jobid = jobid)

    return (experiment_accession)
}

download_narrowpeak_files <- function (accession, jobid){
    # Function to download the narrowpeak file
    #
    # Args:
    #     accession: The accession number of the experiment
    # Returns:
    #     None: Get the ids of the narrowpeak files

    # Set up the url for the API call
    url <- paste0('https://www.encodeproject.org/experiments/',accession,'/?format=json')
    # print(url)

    response <- GET(url)
    if (response$status_code != 200) {
    return(paste0("Error: ", response$status_code))
    }

    # Parse the response
    response <- fromJSON(content(response, "text", encoding = "UTF-8"))

    # Get mask for the narrowpeak file format type
    mask <- response$files$file_format_type == 'narrowPeak'
    mask <- response$files$preferred_default[mask]
    hrefs <- na.omit(response$files$href[mask])

    # Create directory for the narrowpeak files based on the jobid
    dir.create(paste0('jobs/',jobid), recursive = TRUE)

    # Download these narrowpeak files
    for (href in hrefs){
        print(href)
        filename <- unlist(str_split(href, '/'))[length(unlist(str_split(href, '/')))]
        download.file(paste0('https://www.encodeproject.org',href), destfile = paste0('jobs/',jobid,'/',filename))
    }
}




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
                                 textInput("tcp_name", "Enter transcription factor name", value = "CTCF"),
                                    selectInput("biosample_classification", "Select biosample classification",
                                                choices = c('dummytissue', 'dummycell', 'dummyorganism')),
                                    selectInput("biosample", "Select biosample",
                                                choices = c('dummytissue', 'dummycell', 'dummyorganism')),
                                    actionButton("submit", "Submit")

                             ),
                         

                                mainPanel(
                                    textOutput(outputId = "chip_seq_biosample_classification"),
                                )
                         ),
                ),
)


server <- function(input, output, session) {

    observe({
        new_choices <- initial_search(input$tcp_name)
        updateSelectInput(
            session = session, 
            inputId = "biosample_classification",
            choices = new_choices)

    })

    observe({
        new_choices <- search_based_biosample_classifications(input$tcp_name, input$biosample_classification)
        updateSelectInput(
            session = session, 
            inputId = "biosample",
            choices = new_choices)

    })

    output$chip_seq_biosample_classification <- renderText({
        paste0("The input transcription factor is ", input$tcp_name, " and the biosample classification is ", input$biosample_classification, " and the biosample is ", input$biosample)
    })

    observeEvent(input$submit, {
        jobid <- 'sampleid'
        search_based_biosample(input$biosample_classification, input$tcp_name, input$biosample, jobid = jobid)
    })



}

shinyApp(ui = ui, server = server)

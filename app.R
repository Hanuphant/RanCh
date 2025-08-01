library(httr)
library(jsonlite)
library(dplyr)
library(stringr)
library(tidyr)
# library(markdown)


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")

library(DESeq2)


# BETApath <- "/projects/team5/conda/envs/beta/bin/BETA"
BETApath <- "~/miniconda3/envs/beta/bin/BETA"

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
    response <- fromJSON(httr::content(response, "text", encoding = "UTF-8"))

    # Get the biosample classifications
    biosample_classifications <- response$facets$terms[[9]]$key

    return(unlist(biosample_classifications))
}


# Now with the biosample_classifications the user will have the option to select which biological sample does the RNA-Seq belong to.

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
    response <- fromJSON(httr::content(response, "text", encoding = "UTF-8"))
    biosamples <- response$facets$terms[[10]]$key

    return(unlist(biosamples))
}

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
    response <- fromJSON(httr::content(response, "text", encoding = "UTF-8"))
    experiment_accession <- response$`@graph`$accession[1]

    # Response from the experiment accession number
    listedreturn <- download_narrowpeak_files(accession = experiment_accession, jobid = jobid)
    return(listedreturn)
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
    response <- fromJSON(httr::content(response, "text", encoding = "UTF-8"))

    # Get mask for the narrowpeak file format type
    narrowPeakmask <- response$files$file_format_type == 'narrowPeak'
    filetypemask <- response$files$file_type == 'bed narrowPeak'
    defaultmask <- response$files$preferred_default
    hrefs <- na.omit(response$files$href[narrowPeakmask & filetypemask & defaultmask])

    assembly <- na.omit(response$files$assembly[narrowPeakmask & filetypemask & defaultmask])

    # If hrefs is more than one, then we need to download the first one
    if (length(hrefs) > 1){
        hrefs <- hrefs[1]
    }

    # If assembly is more than one, retrieve the first one

    if (length(assembly) > 1){
        assembly <- assembly[1]
    }

    # For GRCh38 and hg38 return hg38
    if (assembly == 'GRCh38' || assembly == 'hg38'){
        assembly <- 'hg38'
    }

    # For GRCh37 and hg19 return hg19
    if (assembly == 'GRCh37' || assembly == 'hg19'){
        assembly <- 'hg19'
    }

    # Create directory for the narrowpeak files based on the jobid
    dir.create(paste0('jobs/',jobid), recursive = TRUE)

    # Download these narrowpeak files

    filename <- unlist(str_split(hrefs, '/'))[length(unlist(str_split(hrefs, '/')))]
    download.file(paste0('https://www.encodeproject.org',hrefs), destfile = paste0('jobs/',jobid,'/',filename))

    return (list(hrefs = filename, assembly = assembly))
}

generateJobTitle <- function(n = 5000) {
  a <- do.call(paste0, replicate(8, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}


library(shiny)
library(shinythemes)

# Define UI for application that takes in RNA-Seq featurecounts and ChIP-Seq narrowPeak files
ui <- fluidPage(theme = shinytheme("cerulean"),
                shinyjs::useShinyjs(),
                titlePanel("RanCh: RNA-Seq and ChIP-Seq Analysis"),
                tabPanel("Tool",
                         sidebarLayout(
                             sidebarPanel(
                                fileInput("rna_seq_counts", "Select RNA-Seq counts matrix file",
                                           accept = c("text/csv", "text/comma-separated-values,text/plain",".csv", '.tabular', '.tsv')),
                                fileInput("rna_seq_meta", "Select RNA-Seq metadata file",
                                          accept = c("text/csv", "text/comma-separated-values,text/plain",".csv", '.tabular', '.tsv')), 
                                selectInput("factor", "Select the factor to test",
                                            choices = c("None")),
                                 textInput("tcp_name", "Enter transcription factor name", value = "CTCF"), # nolint
                                    selectInput("biosample_classification", "Select biosample classification", # nolint # nolint
                                                choices = c('dummytissue', 'dummycell', 'dummyorganism')),
                                    selectInput("biosample", "Select biosample",
                                                choices = c('dummytissue', 'dummycell', 'dummyorganism')),
                                    actionButton("submit", "Submit")

                             ),
                         

                                mainPanel(
                                    tabsetPanel(
                                        tabPanel("Tool",                                     
                                            h3(textOutput(outputId = "chip_seq_biosample_classification")),

                                            textOutput(outputId = "betaLog"),

                                            h2(textOutput(outputId = "UP")),

                                            # Render a table
                                            fluidRow(
                                                column(7,
                                                    dataTableOutput(outputId = "upreg")

                                                )
                                            ),
                                            h3(textOutput(outputId = "DOWN")),
                                            
                                            fluidRow(
                                                column(7,
                                                    dataTableOutput(outputId = "downreg")
                                                )
                                            ), 

                                            # Download button to download the results
                                            downloadButton("downloadData", "Download Output Data")
                                        ),
                                        tabPanel("How it works", 
                                            fluidRow(
                                                column(12, includeHTML("howitworks.html"))
                                            ),
                                            # plotOutput(outputId = "pipeline")
                                            tags$img(src='Pipeline_RanCH.png', style="display: block; margin-left: auto; margin-right: auto;", width = "70%", height = "70%"),
                                            fluidRow(
                                                column(12, includeHTML("howitworks2.html"))
                                            )


                                        )
                                        )                                    
                                )
                         ),
                )
)


server <- function(input, output, session) {

    # output$pipeline <- renderImage({
    #     src <- "/home/shrey/ENCODE2/team5/Pipeline_RanCH.png"
    #     paste0("<img src='", src, "' style='width:100%;height:100%;'/>")
    # })
    

    observeEvent(input$rna_seq_meta,{
        print(input$rna_seq_meta)
        observe_metadata <- read.table(input$rna_seq_meta$datapath, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
        new_choices <- colnames(observe_metadata)
        updateSelectInput(
            session = session, 
            inputId = "factor", 
            choices = new_choices)
    })

    observe({
        jobid <<- generateJobTitle(1)
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

        listedreturn <- search_based_biosample(input$biosample_classification, input$tcp_name, input$biosample, jobid = jobid)
        narrowPeakFile <- listedreturn$hrefs
        assembly <- listedreturn$assembly
        
        withProgress(message = 'Running RanCh', value = 0, {
        n <- 3
        tryCatch({
        # Read the counts matrix
        countData <- read.table(input$rna_seq_counts$datapath, sep = '\t', header = TRUE, row.names = 1)

        # Read the metadata
        metadata <- read.table(input$rna_seq_meta$datapath, sep = '\t', header = TRUE, row.names = 1)
        }, error = function(e){
            print(e)
            showNotification(paste0("Error in reading the counts matrix or metadata files. Please check them and try again. : ", e$message), type = "error", duration = 10000)
        })
        incProgress(1/n, detail = paste("Running DESeq2"))

        # Ordering the columns of the counts matrix based on the metadata
        metadata <- metadata[match(colnames(countData), row.names(metadata)), ]

        # Delete the rna_seqs files
        if (file.exists(input$rna_seq_counts$datapath)){
            unlink(input$rna_seq_counts$datapath)
        }

        if (file.exists(input$rna_seq_meta$datapath)){
            unlink(input$rna_seq_meta$datapath)
        }

        # Delete the files
        countData <- as.data.frame(countData)
        metadata <- as.data.frame(metadata)
        print(type(input$factor))
        print(colnames(metadata[c(input$factor)]))

        # Catch DESeq2 errors
        tryCatch({
            dds <- DESeqDataSetFromMatrix(countData = countData, colData = metadata, design = formula(paste("~", input$factor)))
        }, error = function(e) {
            print(e)
            showNotification(paste0("Error in DESeq2. Please check the metadata file and try again. : ", e$message), type = "error", duration = 10000)
        })
        # dds <- DESeqDataSetFromMatrix(countData = countData, colData = metadata, design = formula(paste("~", input$factor)))
        dds <- DESeq(dds)



        res <- results(dds)
        resdf <- as.data.frame(res)
        colnames(resdf) <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
        resdf <- resdf[c("log2FoldChange", "pvalue")]

        # Remove the rows with NA values
        resdf <- na.omit(resdf)

        # Filter the results based on the pvalue and log2FoldChange
        resdf <- resdf[which(resdf$pvalue < 0.05), ]
        resdf <- resdf[which(abs(resdf$log2FoldChange) > 1), ]
        incProgress(1/n, detail = paste("Running BETA"))
        # Write the results to a file
        write.csv(resdf, file = paste0('jobs/',jobid,'/results.csv'))

        print("Running name change")
        res <- read.csv(paste0('jobs/',jobid,'/results.csv'), header = TRUE, sep = ',')
        colnames(res) <- c("GeneName", "log2FoldChange", "pvalue")
        write.table(res, file = paste0('jobs/',jobid,'/results2.tsv'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

        print("Extracting the narrowPeak file")
        R.utils::gunzip(paste0('jobs/',jobid,'/',narrowPeakFile), paste0('jobs/',jobid,'/chip_seq.bed'))

        # Change directory where shell commands can be run
        wd_before <- getwd()
        setwd(paste0('jobs/',jobid, '/'))

        print("Running BETA")
        system(paste0(BETApath, ' basic -p chip_seq.bed -e results2.tsv -g ',assembly,' -k BSF -n ',jobid,' -o . --gname2 -c 0.05 > beta.logs'), intern = TRUE)   
        
        # Change directory back to original
        setwd(wd_before)

        if (file.exists(paste0('jobs/',jobid,'/',jobid,'_downtarget.txt')) || file.exists(paste0('jobs/',jobid,'/',jobid,'_uptarget.txt'))){

        
        if (file.exists(paste0('jobs/',jobid,'/',jobid,'_downtarget.txt'))) {   
            output$DOWN <- renderText({
            paste0("DOWNREGULATED GENES")
        })

            downreg = read.table(paste0('jobs/',jobid,'/',jobid,'_downtarget.txt'), sep = '\t', header = TRUE)
            output$downreg <- renderDataTable(downreg)
        }

        if (file.exists(paste0('jobs/',jobid,'/',jobid,'_uptarget.txt'))) {

            output$UP <- renderText({
            paste0("UPREGULATED GENES")
        })

            upreg = read.table(paste0('jobs/',jobid,'/',jobid,'_uptarget.txt'), sep = '\t', header = TRUE)
            output$upreg <- renderDataTable(upreg)
        }

        if (file.exists(paste0('jobs/',jobid,'/',jobid,'_downtarget.txt')) || file.exists(paste0('jobs/',jobid,'/',jobid,'_uptarget.txt'))){
        shinyjs::show("downloadData")
        }

        } else {
            showNotification("Error in BETA. Please try other transcription factor or check other parameters.", type = "error", duration = 10000)
            betaLog <- readLines(paste0('jobs/',jobid,'/beta.logs'))
            output$betaLog <- renderText({
                paste0("Please check the beta logs and try again.\n\n", betaLog[length(betaLog)])
            })
        }
        incProgress(1/n, detail = "RanCh finished. Check output.")
        })       
    })

    observe({
        shinyjs::hide("downloadData")

        if (file.exists(paste0('jobs/',jobid))){
            if (file.exists(paste0('jobs/',jobid,'/',jobid,'_downtarget.txt')) || file.exists(paste0('jobs/',jobid,'/',jobid,'_uptarget.txt'))){
                shinyjs::show("downloadData")
            }
        }
    })



    output$downloadData <- downloadHandler(
        filename = function() {
            paste0(jobid, "_results.tar")
        },
        content = function(file) {
            tar(file, paste0('jobs/',jobid,'/'))
        }
    )

    session$onSessionEnded(function() {
        unlink(paste0('jobs/'), recursive = TRUE)
    })
      
}

shinyApp(ui = ui, server = server)


# Import library for API calls
library(httr)
library(jsonlite)
library(dplyr)
library(stringr)
library(tidyr)

export('initial_search', 'search_based_biosample_classifications', 'search_based_biosample')

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
    stop("Error: ", response$status_code)
    }

    # Parse the response
    response <- fromJSON(content(response, "text", encoding = "UTF-8"))

    # Get the biosample classifications
    biosample_classifications <- response$facets$terms[[9]]$key

    return(biosample_classifications)
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
    print(url)
    # Make the API call
    response <- GET(url)
    if (response$status_code != 200) {
    stop("Error: ", response$status_code)
    }

    # Parse the response
    response <- fromJSON(content(response, "text", encoding = "UTF-8"))
    biosamples <- response$facets$terms[[10]]$key

    return(biosamples)
}

#TODO: User input for biological sample here

search_based_biosample <- function(biosample_classification, tcp_name, biosample){

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

    url <- paste0('https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false&assay_title=TF+ChIP-seq&target.label=',target_label,'&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&files.file_type=bed+narrowPeak&biosample_ontology.classification=',biosample_classification,'&biosample_ontology.term_name=',biosample,'&format=json')
    # Make the API call
    response <- GET(url)
    if (response$status_code != 200) {
    stop("Error: ", response$status_code)
    }

    # Parse the response
    response <- fromJSON(content(response, "text", encoding = "UTF-8"))
    experiment_accession <- response$`@graph`$accession[1]

    # Response from the experiment accession number
    download_narrowpeak_files(experiment_accession)

    return (experiment_accession)
}

download_narrowpeak_files <- function (accession){
    # Function to download the narrowpeak file
    #
    # Args:
    #     accession: The accession number of the experiment
    # Returns:
    #     None: Get the ids of the narrowpeak files

    # Set up the url for the API call
    url <- paste0('https://www.encodeproject.org/experiments/',accession,'/?format=json')

    response <- GET(url)
    if (response$status_code != 200) {
    stop("Error: ", response$status_code)
    }

    # Parse the response
    response <- fromJSON(content(response, "text", encoding = "UTF-8"))

    # Get mask for the narrowpeak file format type
    mask <- response$files$file_format_type == 'narrowPeak'
    hrefs <- na.omit(response$files$href[mask])

    # Download these narrowpeak files
    for (href in hrefs){
        filename <- unlist(str_split(href, '/'))[length(unlist(str_split(href, '/')))]
        download.file(paste0('https://www.encodeproject.org',href), destfile = paste0(filename))
    }
}



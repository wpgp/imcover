
#' Download WUENIC data files
#'
#' Retrieve the input data files and WUENIC coverage estimates with further
#' processing of \code{ic.df} objects.
#' @param url a character string naming the URL of the resource to download.
#'   Default is WHO immunization data website.
#' @param destfile Optional. Character string with the name of where the
#'   downloaded file is saved.
#' @param quiet Should status messages and progress bars be hidden? Default is
#'   \code{FALSE}.
#' @param attempts integer of the number of attempts to download the file.
#' @param mode a character string specifying the write mode. 'wb' is most common
#'   to write binary files.
#' @param ... additional arguments to be passed on to \code{download.file}.
#' @details The attributes of the ic data are used to match columns. Column
#'   names are taken from the first object.
#' @return \code{download_wuenic} returns a \code{data.frame} while
#'   \code{download_coverage} and \code{downlaod_survey} return objects of type
#'   \code{ic.df}.
#' @name download
#' @export
download_wuenic <- function(url, destfile,
                            quiet = FALSE, attempts = 3, mode = 'wb', ...){
  if(missing(url)){
    url <- 'https://immunizationdata.who.int/assets/additional-data/wuenic_input_to_pdf.xlsx'
  }

  if(missing(destfile)){
    destfile <- tempfile(fileext = '.xlsx')
  }

  tries <- 1
  retval <- 1

  while(retval != 0L && tries <= attempts){
    retval <- download.file(url,
                            desfile,
                            method = 'auto',
                            quiet = quiet,
                            mode = mode)
  }

  dat <- data.frame(readxl::read_excel(destfile, sheet = 2))
  return(dat)
}


#' @name download
#' @export
download_coveraage <- function(){

}


#' @name download
#' @export
download_survey <- function(){

}

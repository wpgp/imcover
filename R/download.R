
#' Download WUENIC data files
#'
#' Retrieve the input data files and WUENIC coverage estimates with further
#' processing of \code{ic.df} objects.
#' @param destfile Optional. Character string with the name of where the
#' downloaded file is saved.
#' @param url Optional. A character string naming the URL of the resource to
#'   download. Default is the .xlsx from the WHO immunization data website.
#' @param quiet Should status messages and progress bars be hidden? Default is
#'   \code{FALSE}.
#' @param attempts integer of the number of attempts to download the file.
#'   Default is 3.
#' @param mode a character string specifying the write mode. 'wb' is most common
#'   to write binary files.
#' @param return_ic Should an object of type \code{ic.df} be returned? Default
#'   is \code{TRUE}.
#' @param ... additional arguments to be passed on to \code{download.file}.
#' @return The functions return a \code{data.frame} while
#'   \code{download_coverage} and \code{downlaod_survey} return objects of type
#'   \code{ic.df} unless the option \code{return_ic} is \code{FALSE}.
#' @seealso \code{\link[imcover]{ic_data}}
#' @name download
#' @export
download_wuenic <- function(destfile, url,
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
                            destfile,
                            method = 'auto',
                            quiet = quiet,
                            mode = mode,
                            ...)
  }

  if(retval != 0){
    strop('Error downloading file.')
  }

  dat <- data.frame(readxl::read_excel(destfile, sheet = 2))
  return(dat)
}


#' @name download
#' @export
download_coverage <- function(destfile, url,
                               quiet = FALSE, attempts = 3, mode = 'wb',
                               return_ic = TRUE, ...){
  if(missing(url)){
    url <- 'https://whowiise.blob.core.windows.net/upload/coverage--2020.xlsx'
  }

  if(missing(destfile)){
    destfile <- tempfile(fileext = '.xlsx')
  }

  tries <- 1
  retval <- 1

  while(retval != 0L && tries <= attempts){
    retval <- download.file(url,
                            destfile,
                            method = 'auto',
                            quiet = quiet,
                            mode = mode,
                            ...)
  }

  if(retval != 0){
    stop('Error downloading file.')
  }

  # data cleaning
  dat <- data.frame(readxl::read_excel(destfile, sheet = 1))
  names(dat) <- tolower(names(dat))

  dat <- dat[dat$group == 'Countries', ]  # drop aggregations
  dat$group <- NULL  # drop from data.frame

  # drop empty records
  dat <- dat[!is.na(dat$coverage) | (!is.na(dat$target_number) & !is.na(dat$doses)), ]
  # add regional ID
  dat$region <- get_region(dat$code, type = 'm49')

  if(return_ic){
    dat <- ic_data(dat)
  } else{
    return(dat)
  }

}


#' @name download
#' @export
download_survey <- function(){

}

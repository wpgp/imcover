
#' Download WUENIC data files
#'
#' Retrieve the input data files and WUENIC coverage estimates with further
#' processing of \code{ic.df} objects.
#' @param destfile Optional. Character string with the name of where the
#' downloaded file is saved.
#' @param url Optional. A character string naming the URL of the resource to
#'   download. Default is the .xlsx from the WHO immunization data website.
#' @param use_cache Logical. Look for an already downloaded file is the
#'   \code{destfile} location? Default is \code{TRUE}.
#' @param quiet Should status messages and progress bars be hidden? Default is
#'   \code{FALSE}.
#' @param attempts integer of the number of attempts to download the file.
#'   Default is 3.
#' @param mode a character string specifying the write mode. 'wb' is most common
#'   to write binary files.
#' @param return_ic Should an object of type \code{ic.df} be returned? Default
#'   is \code{TRUE}.
#' @param reduce For surveys only, should records with insufficient coverage
#'   be dropped? Default is \code{TRUE}.
#' @param minSample For surveys only, what is the minimum sample size to keep
#'   survey records? Default is 300.
#' @param adjust For survey only, should a recall-bias adjustment be applied?
#'   Default is \code{TRUE}.
#' @param adjVacc Character vector of vaccines to adjust when \code{adjust} is
#'   \code{TRUE}. Default is "DTP" and "PCV".
#' @param add_region Optional. Character string for the type of region code to
#'   add based on the ISO3 code of the country. Options are: 'who' or 'm49'.
#' @param ... additional arguments to be passed on to \code{download.file}.
#' @return The functions return a \code{data.frame} while
#'   \code{download_coverage} and \code{downlaod_survey} return objects of type
#'   \code{ic.df} unless the option \code{return_ic} is \code{FALSE}.
#' @seealso \code{\link[imcover]{ic_data}}
#' @name download
#' @export
download_wuenic <- function(destfile, url, use_cache = TRUE,
                            quiet = FALSE, attempts = 3, mode = 'wb',
                            return_ic = TRUE,
                            add_region = 'who', ...){
  if(missing(url)){
    url <- 'https://immunizationdata.who.int/docs/librariesprovider21/additional-datasets/wuenic-input-to-pdf.xlsx'
  }

  if(missing(destfile)){
    destfile <- file.path(tempdir(), 'wuenic.xlsx') # tempfile(fileext = '.xlsx')
  }

  tries <- 1
  retval <- 1

  if(use_cache){
    if(file.exists(destfile)) retval <- 0
  }

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

  dat <- data.frame(readxl::read_excel(destfile, sheet = 2, na = c("", "NA")))

  if(return_ic){
    admin <- dat[, c("ISOCountryCode", "Year", "Vaccine", "AdministrativeCoverage")]
    admin$ISOCountryCode <- toupper(admin$ISOCountryCode)
    names(admin)[4] <- 'coverage'
    admin$source <- "admin"

    offic <- dat[, c("ISOCountryCode", "Year", "Vaccine", "GovernmentEstimate")]
    offic$ISOCountryCode <- toupper(offic$ISOCountryCode)
    names(offic)[4] <- 'coverage'
    offic$source <- "official"

    if(!is.null(add_region)){
      stopifnot(length(add_region) == 1L && add_region %in% c('who', 'm49'))
      admin$region <- get_region(admin$ISOCountryCode, type = add_region)
      offic$region <- get_region(offic$ISOCountryCode, type = add_region)
    }

    dat <- rbind(admin, offic)
    # drop missing
    dat <- dat[!is.na(dat$coverage), ]

    dat <- ic_data(dat,
                   region = 'region',
                   country = 'ISOCountryCode',
                   time = 'Year',
                   vaccine = 'Vaccine',
                   coverage = 'coverage',
                   source = 'source')
  }

  return(dat)
}


#' @name download
#' @export
download_coverage <- function(destfile, url, use_cache = TRUE,
                              quiet = FALSE, attempts = 3, mode = 'wb',
                              return_ic = TRUE,
                              add_region = 'who', ...){
  if(missing(url)){
    url <- 'https://whowiise.blob.core.windows.net/upload/coverage--2020.xlsx'
  }

  if(missing(destfile)){
    destfile <- file.path(tempdir(), 'coverage.xlsx')  # tempfile(fileext = '.xlsx')
  }

  tries <- 1
  retval <- 1

  if(use_cache){
    if(file.exists(destfile)) retval <- 0
  }

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
  if(!is.null(add_region)){
    stopifnot(length(add_region) == 1L && add_region %in% c('who', 'm49'))
    dat$region <- get_region(dat$code, type = add_region)
  }

  if(return_ic){
    dat <- ic_data(dat)
  } else{
    return(dat)
  }

}


#' @name download
#' @export
download_survey <- function(destfile, url, use_cache = TRUE,
                            quiet = FALSE, attempts = 3, mode = 'wb',
                            return_ic = TRUE,
                            reduce = TRUE, minSample = 300,
                            adjust = TRUE, adjVacc = c("DTP", "PCV"),
                            add_region = 'who', ...){

  if(missing(url)){
    url <- 'https://immunizationdata.who.int/docs/librariesprovider21/additional-datasets/coverage-survey-data.xlsx'
  }

  if(missing(destfile)){
    destfile <- file.path(tempdir(), 'survey.xlsx')  # tempfile(fileext = '.xls')
  }

  tries <- 1
  retval <- 1

  if(use_cache){
    if(file.exists(destfile)) retval <- 0
  }

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
  dat <- data.frame(readxl::read_excel(destfile, skip = 1, sheet = 2))

  # drop empty records
  dat <- dat[!is.na(dat$coverage), ]
  # add source
  dat$source <- 'survey'
  # add regional ID
  if(!is.null(add_region)){
    stopifnot(length(add_region) == 1L && add_region %in% c('who', 'm49'))
    dat$region <- get_region(dat$ISO3, type = add_region)
  }

  if(return_ic){
    dat <- ic_survey(dat, region = 'region', country = 'ISO3',
                     time = 'cohortYear', vaccine = 'vaccine',
                     coverage = 'coverage', source = 'source',
                     reduce = reduce, minSample = minSample)

    dat <- survey_adjust(dat, adjVacc = adjVacc)
  } else{
    return(dat)
  }
}


#' @name download
#' @export
download_denom <- function(destfile, url, use_cache = TRUE,
                           quiet = FALSE, attempts = 3, mode = 'wb', ...){
  if(missing(url)){
    url <- 'https://immunizationdata.who.int/docs/librariesprovider21/additional-datasets/wuenic-input-to-pdf.xlsx'
  }

  if(missing(destfile)){
    destfile <- file.path(tempdir(), 'wuenic.xlsx') # tempfile(fileext = '.xlsx')
  }

  tries <- 1
  retval <- 1

  if(use_cache){
    if(file.exists(destfile)) retval <- 0
  }

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

  dat <- data.frame(readxl::read_excel(destfile, sheet = 2, na = c("", "NA")))

  denom <- dat[, c("ISOCountryCode", "Year", "Vaccine", "SurvivingInfantsUNPD")]
  denom$country <- toupper(denom$ISOCountryCode)
  denom$vaccine <- toupper(denom$Vaccine)
  denom$time <- denom$Year
  names(denom)[which(names(denom) == "SurvivingInfantsUNPD")] <- 'population'

  return(denom[, c('country', 'time', 'vaccine', 'population')])
}


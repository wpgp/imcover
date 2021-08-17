## code to prepare `yovi` dataset goes here
# making internal data
library(readxl)

path <- "data-raw"

## Older yovi data
# sh <- excel_sheets(file.path(getwd(), path, "year_vaccine_introduction.xlsx"))
# sh <- sh[-1]
#
# res <- lapply(sh, function(s){
#   print(s)
#   sheet <- read_excel(file.path(getwd(), path, "year_vaccine_introduction.xlsx"),
#                       sheet = s, na = "n/a")
#   sheet$vaccine <- s
#   return(sheet)
# })
# # match names
# nms <- lapply(res, function(r) names(r))
# nms <- unique(unlist(nms))
#
# res <- lapply(res, function(r){
#   r[, setdiff(nms, names(r))] <- NA
#   return(r)
# })
# # combine
# res <- do.call(rbind, res)
#
# # clean up
# res <- as.data.frame(res)
#
# names(res)[1] <- "ISO3code"
# names(res)[3] <- "WHOregion"
# names(res)[5] <- "year_introduced"
#
# yovi <- res[, c("ISO3code", "WHOregion", "vaccine", "year_introduced")]
# yovi[yovi$year_introduced %in% c("No", "Not introduced", "n/d"), "year_introduced"] <- NA
# yovi$year_introduced <- gsub("Prior or = ", "", yovi$year_introduced, fixed = T)
# yovi$year_introduced <- gsub("prior or = ", "", yovi$year_introduced, fixed = T)
# yovi$year_introduced <- gsub("prior or= ", "", yovi$year_introduced, fixed = T)
# yovi$year_introduced <- gsub("prior to", "", yovi$year_introduced, fixed = T)
# yovi$year_introduced <- gsub("to ", "", yovi$year_introduced, fixed = T)
#
# yovi$year_introduced <- as.numeric(yovi$year_introduced)
# table(yovi$year_introduced)

## import yovi data from 2021
ls <- list.files(path = path,
                 pattern = glob2rx("Introduction*.xlsx"),
                 full.names = T)

yovi <- lapply(ls, function(l){
  dat <- readxl::read_excel(l)
  dat_l <- reshape2::melt(dat,
                          id.vars = 'Country / Region',
                          variable.name = 'year')
  dat_l$year <- as.numeric(levels(dat_l$year)[as.numeric(dat_l$year)])

  dat_l[dat_l$value != 'Yes', 'value'] <- 'No'
  dat_l <- dat_l[dat_l$value == 'Yes', ]

  dat_agg <- aggregate(list('year_introduced' = dat_l$year),
                       by = list('country' = dat_l$`Country / Region`),
                       min)

  if(grepl("Measles", l)){
    dat_agg$vaccine <- "MCV2"
  } else if(grepl("PCV", l)){
    dat_agg$vaccine <- "PCV3"
  }

  dat_agg$ISO3Code <- countrycode::countryname(dat_agg$country, destination = 'iso3c')
  return(dat_agg)
})

yovi <- do.call(rbind.data.frame, yovi)
yovi <- yovi[, c("ISO3Code", "vaccine", "year_introduced")]
yovi <- yovi[order(yovi$ISO3Code, yovi$vaccine, yovi$year_introduced), ]
row.names(yovi) <- 1:nrow(yovi)

## code to create region codes list goes here

whoregion <- unique(yovi[, c("ISO3code", "WHOregion")])
names(whoregion) <- c("ISO3code", "region")

# download data from: https://unstats.un.org/unsd/methodology/m49/overview/

path <- "data-raw"

dat <-read.csv(file.path(path, 'UNSD_m49.csv'), stringsAsFactors = F)
names(dat) <- tolower(names(dat))

m49region <- dat[, c("iso.alpha3.code", "region.name", "sub.region.name")]
names(m49region) <- c("ISO3code", "region", "subregion")

## OUTPUT
usethis::use_data(yovi, m49region, whoregion, overwrite = TRUE, internal = TRUE)


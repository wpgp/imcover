## code to prepare `yovi` dataset goes here
# making internal data
library(readxl)

path <- "data-raw"

sh <- excel_sheets(file.path(getwd(), path, "year_vaccine_introduction.xlsx"))
sh <- sh[-1]

res <- lapply(sh, function(s){
  print(s)
  sheet <- read_excel(file.path(getwd(), path, "year_vaccine_introduction.xlsx"),
                      sheet = s, na = "n/a")
  sheet$vaccine <- s
  return(sheet)
})
# match names
nms <- lapply(res, function(r) names(r))
nms <- unique(unlist(nms))

res <- lapply(res, function(r){
  r[, setdiff(nms, names(r))] <- NA
  return(r)
})
# combine
res <- do.call(rbind, res)

# clean up
res <- as.data.frame(res)

names(res)[1] <- "ISO3code"
names(res)[3] <- "WHOregion"
names(res)[5] <- "year_introduced"

yovi <- res[, c("ISO3code", "WHOregion", "vaccine", "year_introduced")]
yovi[yovi$year_introduced %in% c("No", "Not introduced", "n/d"), "year_introduced"] <- NA
yovi$year_introduced <- gsub("Prior or = ", "", yovi$year_introduced, fixed = T)
yovi$year_introduced <- gsub("prior or = ", "", yovi$year_introduced, fixed = T)
yovi$year_introduced <- gsub("prior or= ", "", yovi$year_introduced, fixed = T)
yovi$year_introduced <- gsub("prior to", "", yovi$year_introduced, fixed = T)
yovi$year_introduced <- gsub("to ", "", yovi$year_introduced, fixed = T)

yovi$year_introduced <- as.numeric(yovi$year_introduced)
table(yovi$year_introduced)


## code to create region codes list goes here
# download data from: https://unstats.un.org/unsd/methodology/m49/overview/

path <- "data-raw"

dat <-read.csv(file.path(path, 'UNSD_m49.csv'), stringsAsFactors = F)
names(dat) <- tolower(names(dat))

region <- dat

## OUTPUT
usethis::use_data(yovi, region, overwrite = TRUE, internal = TRUE)


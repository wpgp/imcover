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


usethis::use_data(yovi, internal = TRUE)


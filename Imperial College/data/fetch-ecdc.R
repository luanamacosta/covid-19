library(lubridate)
# library(readxl)

# date_offset <- 0
url <- "https://opendata.ecdc.europa.eu/covid19/casedistribution/csv"
# date_iso <- as.character(Sys.Date() - date_offset)
# url <- sprintf(url_string, date_iso)

url_page <- "https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide"
tryCatch({
  code <- download.file(url, "data/COVID-19-up-to-date.csv")
  if (code != 0) {
    stop("Erro ao baixar arquivo")
  }
},
error = function(e) {
  stop(sprintf("Erro ao baixar arquivo '%s': %s, por favor verifique %s",
               url, e$message, url_page))
})


d <- read.csv("data/COVID-19-up-to-date.csv", stringsAsFactors = FALSE)
d$t <- lubridate::decimal_date(as.Date(d$dateRep, format = "%d/%m/%Y"))
d <- d[order(d$'countriesAndTerritories', d$t, decreasing = FALSE), ]
names(d)[names(d) == "countriesAndTerritories"] <- "Paises.e.territorios"
names(d)[names(d) == "deaths"] <- "Mortes"
names(d)[names(d) == "cases"] <- "Casos"
names(d)[names(d) == "dateRep"] <- "DataRep"
saveRDS(d, "data/COVID-19-up-to-date.rds")

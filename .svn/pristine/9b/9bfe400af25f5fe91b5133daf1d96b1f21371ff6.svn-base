library(quantmod)
getSymbols(c("MSFT", "SAP", "ADS", "^GSPC", "^DJI"))
MSFT <- diff(log(MSFT[,6]))[-1]
SAP <- diff(log(SAP[,6]))[-1]
ADS <- diff(log(ADS[,6]))[-1]
GSPC <- diff(log(GSPC[,6]))[-1]
DJI <- diff(log(DJI[,6]))[-1]
stocks <- data.frame(MSFT = as.vector(MSFT), SAP = as.vector(SAP), ADS = as.vector(ADS), GSPC = as.vector(GSPC), DJI = as.vector(DJI))
rownames(stocks) <- index(MSFT)
devtools::use_data(stocks, overwrite = TRUE)

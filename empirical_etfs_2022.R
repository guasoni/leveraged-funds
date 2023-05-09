#library(fredr)
library(tseries)
library(tidyverse)
library(zoo)


# Replace with your own local directory

setwd("C:\\Users\\Paolo\\SynologyDrive\\Dropbox\\research\\meanvar\\Leveraged_Funds\\QF\\Empirics")
ustreas = read_csv("yield-curve-rates-1990-2021.csv", col_types = cols(Date = col_date(format = "%m/%d/%y")) )
ustreas2022 = read_csv("yield-curve-rates-2022.csv")

# Alternative command that downloads from Treasury website and writes locally
# ustreas2022 = read_csv("https://home.treasury.gov/resource-center/data-chart-center/interest-rates/daily-treasury-rates.csv/2022/all?type=daily_treasury_yield_curve&field_tdr_date_value=2022&page&_format=csv", col_types = cols(Date = col_date(format = "%m/%d/%Y")))
# write_csv(ustreas2022,"yield-curve-rates-2022.csv")

ustreas = rbind(ustreas,ustreas2022) %>% arrange(desc(Date))

allticks = read_csv("usfunds.csv")

# Download data from Yahoo - optional

funds = setdiff(sort(unique(allticks$Ticker)),"QQQQ")
download_fund = function (ticker) {
  zts = get.hist.quote(ticker,compression = "d", quote = "AdjClose")
  df = fortify.zoo(zts)
  df$Ticker = ticker
  df
}
listprices = lapply(funds, download_fund)
allprices = do.call("rbind",listprices) 
#write_csv(allprices,"allprices_2022.csv.bz2")

# Process data

allprices = read_csv("allprices_2022.csv.bz2")

allprices = allprices %>% 
  rename(Date = Index, Price = Adjusted) %>% 
  group_by(Ticker) %>% 
  mutate(SimpleRet = Price/lag(Price) - 1) %>%
  filter(!is.na(SimpleRet))

mindate = allprices %>% left_join(allticks) %>% 
  filter(Factor != 1) %>% 
  ungroup() %>% 
  summarise(mindate = min(Date)) %>% .$mindate

ratedf = data.frame(Date = sort(unique(allprices$Date))) %>%
  filter(Date >= mindate) %>%
  left_join(ustreas) %>% rename(Rate = `1 Mo`) %>% select(Date, Rate)
ratedf[,-1] = na.locf(ratedf[,-1])

allprices = allprices %>% 
  left_join(ratedf) %>%
  left_join(allticks) 

allprices = allprices %>%
  mutate(Return = (1 + SimpleRet) /(1 + Rate/100/252) - 1,
         RetExp = (1 + Return) * (1 + Expenses/252) - 1)

indexwide = allprices %>% filter(Factor == 1) %>%
  group_by(Date, Underlying) %>%
  summarise(IndexRet = median(RetExp))

retswidenet = allprices %>% ungroup() %>% 
  select(Date, Underlying, Factor, Return) %>% 
  arrange(Underlying,Date, Factor) %>%
  group_by(Date, Underlying) %>%
  pivot_wider(names_from = Factor, values_from = Return, names_sort = TRUE) %>%
  left_join(indexwide)

retswidegro = allprices %>% ungroup() %>% 
  select(Date, Underlying, Factor, RetExp) %>% 
  arrange(Underlying,Date, Factor) %>%
  group_by(Date, Underlying) %>%
  pivot_wider(names_from = Factor, values_from = RetExp, names_sort = TRUE) %>%
  left_join(indexwide)


retswidenet = retswidenet %>% filter(Date < as.Date("2023-01-01"))
retswidenet = retswidenet %>% filter(!(Underlying %in% c("Dow Jones Industrial","S&P Midcap 400")) || !(Date %in% as.Date(c("2021-05-05","2021-05-06"))))
retswidenet = retswidenet[rowSums(!is.na(retswidenet[,-2:-1]))>=2,]

retswidegro = retswidegro %>% filter(Date < as.Date("2023-01-01"))
retswidegro = retswidegro %>% filter(!(Underlying %in% c("Dow Jones Industrial","S&P Midcap 400")) || !(Date %in% as.Date(c("2021-05-05","2021-05-06"))))
retswidegro = retswidegro[rowSums(!is.na(retswidegro[,-2:-1]))>=2,]

lisfact = c(-3,-2,-1,1,2,3)
underlist = sort(unique(allticks$Underlying))
underlist

# Selects one underlying to perform analysis
#thisund = underlist[3]

report = function (thisund) {
  undfun = subset(retswidenet,(Underlying == thisund) )
  undgro = subset(retswidegro,(Underlying == thisund) )
  
  commondate = as.Date(max(sapply(undfun[,3:(3+length(lisfact)-1)], function (x) min(undfun$Date[!is.na(x)]))),origin = "1970-01-01")
  undfun = subset(undfun,Date >= commondate)
  undgro = subset(undgro,Date >= commondate)
  
  # Keeps only factors with at least one year of data
  keepfact = lisfact[colSums(!is.na(undfun[,as.character(lisfact)]))>100]
  
  nyears = as.numeric(max(undfun$Date)-min(undfun$Date))/365.25
  busdays = nrow(undfun)/nyears
  volat = sd(undfun$IndexRet,na.rm = TRUE)*sqrt(busdays)
  nobs = apply(undfun[,as.character(keepfact)],2,function (x) sum(!is.na(x)))
  
  indret = undfun$IndexRet
  # indretsq = undfun$IndexRet^2
  # indretgro = undgro$IndexRet
  # indretsqgro = undgro$IndexRet^2
  
  # regs   = apply(undfun[,as.character(keepfact)],2,function (x) lm(x ~ indret + indretsq))
  # reggro = apply(undgro[,as.character(keepfact)],2,function (x) lm(x ~ indretgro + indretsqgro))
  regs   = apply(undfun[,as.character(keepfact)],2,function (x) lm(x ~ indret ))
  reggro = apply(undgro[,as.character(keepfact)],2,function (x) lm(x ~ indret ))
  alphas   = sapply(regs,function (x) summary(x)$coefficients[,1])[1,]*busdays*100
  alphagro = sapply(reggro,function (x) summary(x)$coefficients[,1])[1,]*busdays*100
  betas =    sapply(regs,function (x) summary(x)$coefficients[,1])[2,]
  # gammas =   sapply(regs,function (x) summary(x)$coefficients[,1])[3,]
  betstat = (sapply(regs,function (x) (summary(x)$coefficients[2,1]))-keepfact)/sapply(regs,function (x) summary(x)$coefficients[2,2])
  trerrs =   sapply(regs,function (x) summary(x)$coefficients[,2])[1,]*sqrt(busdays)*100
  rsqrd =    sapply(regs,function (x) summary(x)$r.squared)*100
  vols =      apply(undfun[,as.character(keepfact)],2,function (x) sd(x,na.rm = TRUE))*sqrt(busdays)*100
  alphadj    = -6.9282 * alphas   * trerrs /(keepfact^2*(1-keepfact)^2*volat^3)
  alphadjgro = -6.9282 * alphagro * trerrs /(keepfact^2*(1-keepfact)^2*volat^3)
  listicks = subset(allticks, Underlying == thisund)
  listicks = as.character(listicks$Ticker[match(lisfact,listicks$Factor)])
  dataset = rbind(lisfact,100*trerrs,alphas,alphagro,alphadj,alphadjgro,betas,betstat,vols,nobs/busdays)
  rownames(dataset) = c("Factor","Tracking Error","Tracking Diff. (Net)","Tracking Diff. (Gross)", "ImpSpread","ImpSpreadGross","Beta","T Stat","Volatility","Years")
  colnames(dataset) = keepfact
  
  output = as.data.frame(round(t(dataset),2))
  output1 = output[,1]
  output = sapply(output[,-1],function(x) formatC(x,2,format = 'f',digits = 2))
  output = cbind (listicks,output1,output)
  rbind(c(as.character(thisund),rep("",10)),output)
}

output =  do.call("rbind",lapply(underlist,report))

write.table(cbind(rep(" ",nrow(output)),output), "outtable.csv",sep = "&\t",eol = "\\\\ \n",row.names = FALSE,quote=FALSE)


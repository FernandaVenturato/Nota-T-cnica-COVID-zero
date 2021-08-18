# Análise dos dados do Programa COVID-0 - UFJF-GV
# Fernanda Venturato Roquim
#================================================================================
rm(list=ls())

# Chamando pacotes necessários
library(readr)
library(gamlss)
library(gamlss.util)
library(moments)

#================================================================================
# Funcao Centil (centiles.() modificada)

centil <- function (obj, xvar, cent = c(0.4, 2, 10, 25, 50, 75, 90, 98, 
                                        99.6), legend = TRUE, ylab = "y", xlab = "x", main = NULL, 
                    main.gsub = "@", xleg = min(xvar), yleg = max(obj$y), xlim = range(xvar), 
                    ylim = range(obj$y), save = FALSE, plot = TRUE, points = TRUE, 
                    pch = 15, cex = 0.5, col = gray(0.7), col.centiles = 1:length(cent) + 
                      2, lty.centiles = 1, lwd.centiles = 1, ...) 
{
  if (!is.gamlss(obj)) 
    stop(paste("This is not an gamlss object", "\n", ""))
  if (missing(xvar)) {
    xvar <- all.vars(obj$call$formula)[[2]]
    if (any(grepl("data", names(obj$call)))) {
      DaTa <- eval(obj$call[["data"]])
      xvar <- get(xvar, envir = as.environment(DaTa))
    }
  }
  xvarO <- deparse(substitute(xvar))
  xvar <- try(xvar, silent = TRUE)
  if (any(class(xvar) %in% "try-error")) {
    DaTa <- eval(obj$call[["data"]])
    xvar <- get(xvarO, envir = as.environment(DaTa))
  }
  fname <- obj$family[1]
  qfun <- paste("q", fname, sep = "")
  Title <- paste("Centile curves using", fname, sep = " ")
  main <- if (is.null(main)) 
    paste("Centile curves using", fname, sep = " ")
  else gsub(main.gsub, Title, main)
  oxvar <- xvar[order(xvar)]
  oyvar <- obj$y[order(xvar)]
  if (is.matrix(obj$y)) {
    oyvar <- obj$y[, 1][order(xvar)]
    ylim <- range(obj$y[, 1])
    yleg = max(obj$y[, 1])
  }
  if (plot) {
    lty.centiles <- rep(lty.centiles, length(cent))
    lwd.centiles <- rep(lwd.centiles, length(cent))
    col.centiles <- rep(col.centiles, length(cent))
    if (points == TRUE) {
      plot(oxvar, oyvar, type = "h", col = col, pch = pch, 
           cex = cex, xlab = xlab, ylab = ylab, xlim = xlim, 
           ylim, ...)
    }
    else {
      plot(oxvar, oyvar, type = "n", col = col, pch = pch, 
           xlab = xlab, ylab = ylab, xlim = xlim, ylim, 
           ...)
    }
    title(main)
  }
  col <- 3
  lpar <- length(obj$parameters)
  ii <- 0
  per <- rep(0, length(cent))
  for (var in cent) {
    if (lpar == 1) {
      newcall <- call(qfun, var/100, mu = fitted(obj, "mu")[order(xvar)])
    }
    else if (lpar == 2) {
      newcall <- call(qfun, var/100, mu = fitted(obj, "mu")[order(xvar)], 
                      sigma = fitted(obj, "sigma")[order(xvar)])
    }
    else if (lpar == 3) {
      newcall <- call(qfun, var/100, mu = fitted(obj, "mu")[order(xvar)], 
                      sigma = fitted(obj, "sigma")[order(xvar)], nu = fitted(obj, 
                                                                             "nu")[order(xvar)])
    }
    else {
      newcall <- call(qfun, var/100, mu = fitted(obj, "mu")[order(xvar)], 
                      sigma = fitted(obj, "sigma")[order(xvar)], nu = fitted(obj, 
                                                                             "nu")[order(xvar)], tau = fitted(obj, "tau")[order(xvar)])
    }
    ii <- ii + 1
    ll <- eval(newcall)
    if (plot) {
      lines(oxvar, ll, col = col.centiles[ii], lty = lty.centiles[ii], 
            lwd = lwd.centiles[ii], ...)
    }
    per[ii] <- (1 - sum(oyvar > ll)/length(oyvar)) * 100
    if (!save) 
      cat("% of cases below ", var, "centile is ", per[ii], 
          "\n")
  }
  if (plot) {
    if (legend == TRUE) 
      legend(list(x = xleg, y = yleg), legend = cent, col = col.centiles, 
             lty = lty.centiles, lwd = lwd.centiles, ncol = 1, 
             ...)
  }
  if (save) {
    return(cbind(cent, per))
  }
}

#================================================================================

#Lendo os dados

DadosGV <- read_csv("Área de Trabalho/Programa Covid/Modelagem/GV.csv", 
                    col_types = cols(casos = col_integer(), 
                                     index = col_integer(),
                                     data = col_date(format = "%d/%m/%Y"), 
                                     obitos = col_integer()))

DadosIPA <- read_csv("Área de Trabalho/Programa Covid/Modelagem/IPA.csv", 
                     col_types = cols(casos = col_integer(), 
                                      index = col_integer(),
                                      data = col_date(format = "%d/%m/%Y"), 
                                      obitos = col_integer()))

DadosJF <- read_csv("Área de Trabalho/Programa Covid/Modelagem/JF.csv", 
                    col_types = cols(casos = col_integer(), 
                                     index = col_integer(),
                                     data = col_date(format = "%d/%m/%Y"), 
                                     obitos = col_integer()))

DadosBH <- read_csv("Área de Trabalho/Programa Covid/Modelagem/BH.csv", 
                    col_types = cols(casos = col_integer(), 
                                     index = col_integer(),
                                     data = col_date(format = "%d/%m/%Y"), 
                                     obitos = col_integer()))


#================================================================================

# Limpando os dados (remoção de observações em branco)

GV <- na.omit(DadosGV)
IPA <- na.omit(DadosIPA)
JF <- na.omit(DadosJF)
BH <- na.omit(DadosBH)

#================================================================================

# Obtendo algumas medidas descritivas dos dados

tudo <- data.frame(DadosGV$casos, DadosGV$obitos, DadosIPA$casos, DadosIPA$obitos, 
                   DadosJF$casos, DadosJF$obitos, DadosBH$casos, DadosBH$obitos)

# Minimo
for(i in 1:8){
  print(min(na.omit(tudo[,i])))
}

# Maximo
for(i in 1:8){
  print(max(na.omit(tudo[,i])))
}

# Media
for(i in 1:8){
  print(mean(na.omit(tudo[,i])))
}

# Mediana
for(i in 1:8){
  print(median(na.omit(tudo[,i])))
}

# Variancia
for(i in 1:8){
  print(var(na.omit(tudo[,i])))
}

# Desvio padrao
for(i in 1:8){
  print(sd(na.omit(tudo[,i])))
}

# Coeficiente (taxa) de variação das variáveis
for(i in 1:8){
  print(100*sd(na.omit(tudo[,i]))/mean(na.omit(tudo[,i])))
}

# Assimetria
for(i in 1:8){
  print(skewness(na.omit(tudo[,i])))
}

# Curtose
for(i in 1:8){
  print(kurtosis(na.omit(tudo[,i])))
}
#================================================================================

dataX <- c("Jul/2020", "Ago/2020", "Set/2020", "Out/2020", "Nov/2020", "Dez/2020",
           "Jan/2021", "Fev/2021", "Mar/2021", "Abr/2021", "Mai/2021", "Jun/2021",
           "Jul/2021", "Ago/2021")

# Obtendo os graficos de caixas

#Gráfico de caixas dos casos
par(mfrow=c(2,2))
boxplot(GV$casos, main = "Governador Valadares", ylab = "Casos")
boxplot(IPA$casos, main = "Ipatinga", ylab = "Casos")
boxplot(JF$casos, main = "Juiz de Fora", ylab = "Casos")
boxplot(BH$casos, main = "Belo Horizonte", ylab = "Casos")

#Gráfico de caixas dos obitos
boxplot(GV$obitos, main = "Governador Valadares", ylab = "Óbitos")
boxplot(IPA$obitos, main = "Ipatinga", ylab = "Óbitos")
boxplot(JF$obitos, main = "Juiz de Fora", ylab = "Óbitos")
boxplot(BH$obitos, main = "Belo Horizonte", ylab = "Óbitos")

#================================================================================
# Intervalo de Predicao

novasdatas <- read_csv("Área de Trabalho/Programa Covid/Modelagem/PRED.csv", 
                       col_types = cols(data = col_date(format = "%d/%m/%Y")))


#================================================================================
# Analise dos novos casos confirmados em Governador Valadares

#Modelo ZINBI
CGV <- gamlss(formula = casos ~ pb(as.numeric(data)), sigma.formula = ~pb(as.numeric(data)),  
             family = ZINBI, data = GV, trace = T, method=mixed(30,100))

#Análise de Resíduos
rqres.plot(CGV, howmany =8, cex=.5, pch=20, col='black', ylim.all =2, plot.type="all"
           , ylim = c(-2, 2))
plot(CGV)

#Gráfico
centil(CGV, xvar=GV$data, cent = c(5,95), col.cent=c(4,2), 
       main = "Modelo ZINBI para predição de novos casos de COVID-19 em 
       Governador Valadares", ylab = "Novos Casos", xlab = "Data", 
       legend = TRUE, xlim = as.Date(c("2020-07-01", "2021-08-10")), xaxt= "n")

PredCGV <- predict(CGV, newdata = novasdatas, type = "response")
PredCGV

fittedCGV <- fitted(CGV)
fittedCGV

prediCGV <- predictAll(CGV, newdata=novasdatas)
prediCGV

LI_CGV <- qZINBI(0.05, mu=prediCGV$mu, sigma = prediCGV$sigma, nu=prediCGV$nu)
LS_CGV <- qZINBI(0.95, mu=prediCGV$mu, sigma = prediCGV$sigma, nu=prediCGV$nu)
             
points(GV$data, fittedCGV, col="black", type = "l")
points(novasdatas$data,PredCGV,col="black", type = "l",  lty=2)
points(novasdatas$data, LI_CGV, type = "l", col = "blue",  lty=2)
points(novasdatas$data, LS_CGV, type = "l", col = "red",  lty=2)
axis(1, at = pretty(GV$data, n=14), labels = dataX , las = 1)

#================================================================================
# Analise dos novos casos confirmados em Ipatinga

#Modelo ZINBI
CIPA <- gamlss(formula = casos ~ pb(as.numeric(data)), sigma.formula = ~pb(as.numeric(data)),  
              family = ZINBI, data = IPA, trace = T, method=mixed(30,100))

#Análise de Resíduos
rqres.plot(CIPA, howmany =8, cex=.5, pch=20, col='black', ylim.all =2, plot.type="all", ylim = c(-2, 2))
plot(CIPA)

#Gráfico
centil(CIPA, xvar=IPA$data, cent = c(5,95), col.cent=c(4,2), 
       main = "Modelo ZINBI para predição de novos casos de COVID-19 em 
       Ipatinga", ylab = "Novos Casos", xlab = "Data", 
       legend = TRUE, xlim = as.Date(c("2020-07-01", "2021-08-10")), xaxt= "n")

PredCIPA <- predict(CIPA, newdata = novasdatas, type = "response")
PredCIPA

fittedCIPA <- fitted(CIPA)
fittedCIPA

prediCIPA <- predictAll(CIPA, newdata=novasdatas)
prediCIPA

LI_CIPA <- qZINBI(0.05, mu=prediCIPA$mu, sigma = prediCIPA$sigma, nu=prediCIPA$nu)
LS_CIPA <- qZINBI(0.95, mu=prediCIPA$mu, sigma = prediCIPA$sigma, nu=prediCIPA$nu)

points(IPA$data, fittedCIPA, col="black", type = "l")
points(novasdatas$data,PredCIPA,col="black", type = "l",  lty=2)
points(novasdatas$data, LI_CIPA, type = "l", col = "blue",  lty=2)
points(novasdatas$data, LS_CIPA, type = "l", col = "red",  lty=2)
axis(1, at = pretty(IPA$data, n=14), labels = dataX , las = 1)

#================================================================================
# Analise dos novos casos confirmados em Juiz de Fora

#Modelo ZINBI
CJF <- gamlss(formula = casos ~ pb(as.numeric(data)), sigma.formula = ~pb(as.numeric(data)),  
              family = ZINBI, data = JF, trace = T, method=mixed(30,100))

#Análise de Resíduos
rqres.plot(CJF, howmany =8, cex=.5, pch=20, col='black', ylim.all =2, plot.type="all", ylim = c(-2, 2))
plot(CJF)

#Gráfico
centil(CJF, xvar=JF$data, cent = c(5,95), col.cent=c(4,2), 
       main = "Modelo ZINBI para predição de novos casos de COVID-19 em 
       Juiz de Fora", ylab = "Novos Casos", xlab = "Data", 
       legend = TRUE, xlim = as.Date(c("2020-07-01", "2021-08-10")), xaxt= "n")

PredCJF <- predict(CJF, newdata = novasdatas, type = "response")
PredCJF

fittedCJF <- fitted(CJF)
fittedCJF

prediCJF <- predictAll(CJF, newdata=novasdatas)
prediCJF

LI_CJF <- qZINBI(0.05, mu=prediCJF$mu, sigma = prediCJF$sigma, nu=prediCJF$nu)
LS_CJF <- qZINBI(0.95, mu=prediCJF$mu, sigma = prediCJF$sigma, nu=prediCJF$nu)

points(JF$data, fittedCJF, col="black", type = "l")
points(novasdatas$data,PredCJF,col="black", type = "l",  lty=2)
points(novasdatas$data, LI_CJF, type = "l", col = "blue",  lty=2)
points(novasdatas$data, LS_CJF, type = "l", col = "red",  lty=2)
axis(1, at = pretty(JF$data, n=14), labels = dataX , las = 1)

#================================================================================
# Analise dos novos casos confirmados em Belo Horizonte

#Modelo ZINBI
CBH <- gamlss(formula = casos ~ pb(as.numeric(data)), sigma.formula = ~pb(as.numeric(data)),  
              family = ZINBI, data = BH, trace = T, method=mixed(30,100))

#Análise de Resíduos
rqres.plot(CBH, howmany =8, cex=.5, pch=20, col='black', ylim.all =2, plot.type="all", ylim = c(-2, 2))
plot(CBH)

#Gráfico
centil(CBH, xvar=BH$data, cent = c(5,95), col.cent=c(4,2), 
       main = "Modelo ZINBI para predição de novos casos de COVID-19 em 
       Belo Horizonte", ylab = "Novos Casos", xlab = "Data", 
       legend = TRUE, xlim = as.Date(c("2020-07-01", "2021-08-10")), xaxt= "n")

PredCBH <- predict(CBH, newdata = novasdatas, type = "response")
PredCBH

fittedCBH <- fitted(CBH)
fittedCBH

prediCBH <- predictAll(CBH, newdata=novasdatas)
prediCBH

LI_CBH <- qZINBI(0.05, mu=prediCBH$mu, sigma = prediCBH$sigma, nu=prediCBH$nu)
LS_CBH <- qZINBI(0.95, mu=prediCBH$mu, sigma = prediCBH$sigma, nu=prediCBH$nu)

points(BH$data, fittedCBH, col="black", type = "l")
points(novasdatas$data,PredCBH,col="black", type = "l",  lty=2)
points(novasdatas$data, LI_CBH, type = "l", col = "blue",  lty=2)
points(novasdatas$data, LS_CBH, type = "l", col = "red",  lty=2)
axis(1, at = pretty(BH$data, n=14), labels = dataX , las = 1)

#================================================================================
# Analise dos novos obitos confirmados em Governador Valadares

#Modelo ZINBI
OGV <- gamlss(formula = obitos ~ pb(as.numeric(data)), sigma.formula = ~pb(as.numeric(data)),  
              family = ZINBI, data = GV, trace = T, method=mixed(30,100))

#Análise de Resíduos
rqres.plot(OGV, howmany =8, cex=.5, pch=20, col='black', ylim.all =2, plot.type="all", ylim = c(-2, 2))
plot(OGV)

#Gráfico
centil(OGV, xvar=GV$data, cent = c(5,95), col.cent=c(4,2), 
       main = "Modelo ZINBI para predição de novos óbitos de COVID-19 em 
       Governador Valadares", ylab = "Novos Óbitos", xlab = "Data", 
       legend = TRUE, xlim = as.Date(c("2020-07-01", "2021-08-10")), xaxt= "n")

PredOGV <- predict(OGV, newdata = novasdatas, type = "response")
PredOGV

fittedOGV <- fitted(OGV)
fittedOGV

prediOGV <- predictAll(OGV, newdata=novasdatas)
prediOGV

LI_OGV <- qZINBI(0.05, mu=prediOGV$mu, sigma = prediOGV$sigma, nu=prediOGV$nu)
LS_OGV <- qZINBI(0.95, mu=prediOGV$mu, sigma = prediOGV$sigma, nu=prediOGV$nu)

points(GV$data, fittedOGV, col="black", type = "l")
points(novasdatas$data,PredOGV,col="black", type = "l",  lty=2)
points(novasdatas$data, LI_OGV, type = "l", col = "blue",  lty=2)
points(novasdatas$data, LS_OGV, type = "l", col = "red",  lty=2)
axis(1, at = pretty(GV$data, n=14), labels = dataX , las = 1)

#================================================================================
# Analise dos novos obitos confirmados em Ipatinga

#Modelo ZINBI
OIPA <- gamlss(formula = obitos ~ pb(as.numeric(data)), sigma.formula = ~pb(as.numeric(data)),  
               family = ZINBI, data = IPA, trace = T, method=mixed(30,100))

#Análise de Resíduos
rqres.plot(OIPA, howmany =8, cex=.5, pch=20, col='black', ylim.all =2, plot.type="all", ylim = c(-2, 2))
plot(OIPA)

#Gráfico
centil(OIPA, xvar=IPA$data, cent = c(5,95), col.cent=c(4,2), 
       main = "Modelo ZINBI para predição de novos óbitos de COVID-19 em 
       Ipatinga", ylab = "Novos Óbitos", xlab = "Data", 
       legend = TRUE, xlim = as.Date(c("2020-07-01", "2021-08-10")), xaxt= "n")

PredOIPA <- predict(OIPA, newdata = novasdatas, type = "response")
PredOIPA

fittedOIPA <- fitted(OIPA)
fittedOIPA

prediOIPA <- predictAll(OIPA, newdata=novasdatas)
prediOIPA

LI_OIPA <- qZINBI(0.05, mu=prediOIPA$mu, sigma = prediOIPA$sigma, nu=prediOIPA$nu)
LS_OIPA <- qZINBI(0.95, mu=prediOIPA$mu, sigma = prediOIPA$sigma, nu=prediOIPA$nu)

points(IPA$data, fittedOIPA, col="black", type = "l")
points(novasdatas$data,PredOIPA,col="black", type = "l",  lty=2)
points(novasdatas$data, LI_OIPA, type = "l", col = "blue",  lty=2)
points(novasdatas$data, LS_OIPA, type = "l", col = "red",  lty=2)
axis(1, at = pretty(IPA$data, n=14), labels = dataX , las = 1)

#================================================================================
# Analise dos novos obitos confirmados em Juiz de Fora

#Modelo ZINBI
OJF <- gamlss(formula = obitos ~ pb(as.numeric(data)), sigma.formula = ~pb(as.numeric(data)),  
              family = ZINBI, data = JF, trace = T, method=mixed(30,100))

#Análise de Resíduos
rqres.plot(OJF, howmany =8, cex=.5, pch=20, col='black', ylim.all =2, plot.type="all", ylim = c(-2, 2))
plot(OJF)

#Gráfico
centil(OJF, xvar=JF$data, cent = c(5,95), col.cent=c(4,2), 
       main = "Modelo ZINBI para predição de novos óbitos de COVID-19 em 
       Juiz de Fora", ylab = "Novos Óbitos", xlab = "Data", 
       legend = TRUE, xlim = as.Date(c("2020-07-01", "2021-08-10")), xaxt= "n")

PredOJF <- predict(OJF, newdata = novasdatas, type = "response")
PredOJF

fittedOJF <- fitted(OJF)
fittedOJF

prediOJF <- predictAll(OJF, newdata=novasdatas)
prediOJF

LI_OJF <- qZINBI(0.05, mu=prediOJF$mu, sigma = prediOJF$sigma, nu=prediOJF$nu)
LS_OJF <- qZINBI(0.95, mu=prediOJF$mu, sigma = prediOJF$sigma, nu=prediOJF$nu)

points(JF$data, fittedOJF, col="black", type = "l")
points(novasdatas$data,PredOJF,col="black", type = "l",  lty=2)
points(novasdatas$data, LI_OJF, type = "l", col = "blue",  lty=2)
points(novasdatas$data, LS_OJF, type = "l", col = "red",  lty=2)
axis(1, at = pretty(JF$data, n=14), labels = dataX , las = 1)

#================================================================================
# Analise dos novos obitos confirmados em Belo Horizonte

#Modelo ZINBI
OBH <- gamlss(formula = obitos ~ pb(as.numeric(data)), sigma.formula = ~pb(as.numeric(data)),  
              family = ZINBI, data = BH, trace = T, method=mixed(30,100))

#Análise de Resíduos
rqres.plot(OBH, howmany =8, cex=.5, pch=20, col='black', ylim.all =2, plot.type="all", ylim = c(-2, 2))
plot(OBH)

#Gráfico
centil(OBH, xvar=BH$data, cent = c(5,95), col.cent=c(4,2), 
       main = "Modelo ZINBI para predição de novos óbitos de COVID-19 em 
       Belo Horizonte", ylab = "Novos Óbitos", xlab = "Data", 
       legend = TRUE, xlim = as.Date(c("2020-07-01", "2021-08-10")), xaxt= "n")

PredOBH <- predict(OBH, newdata = novasdatas, type = "response")
PredOBH

fittedOBH <- fitted(OBH)
fittedOBH

prediOBH <- predictAll(OBH, newdata=novasdatas)
prediOBH

LI_OBH <- qZINBI(0.05, mu=prediOBH$mu, sigma = prediOBH$sigma, nu=prediOBH$nu)
LS_OBH <- qZINBI(0.95, mu=prediOBH$mu, sigma = prediOBH$sigma, nu=prediOBH$nu)

points(BH$data, fittedOBH, col="black", type = "l")
points(novasdatas$data,PredOBH,col="black", type = "l",  lty=2)
points(novasdatas$data, LI_OBH, type = "l", col = "blue",  lty=2)
points(novasdatas$data, LS_OBH, type = "l", col = "red",  lty=2)
axis(1, at = pretty(BH$data, n=14), labels = dataX , las = 1)

#================================================================================
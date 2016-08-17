############## CELL CYCLE DISTRIBUTION ESTIMATION FROM DAPI INTENSITIES ################
# 
#   licensed under the MIT License:
#     
#   Copyright (c) 2016 Andreas Stengl, David Hoerl, Heinrich Leonhardt and Jonas Helma
#   
#   Permission is hereby granted, free of charge, to any person obtaining a copy of this
#   software and associated documentation files (the "Software"), to deal in the Software
#   without restriction, including without limitation the rights to use, copy, modify, merge,
#   publish, distribute, sublicense, and/or sell copies of the Software, and to permit 
#   persons to whom the Software is furnished to do so, subject to the following conditions:
#     
#   The above copyright notice and this permission notice shall be 
#   included in all copies or substantial portions of the Software.
#   
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
#   TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULARPURPOSE AND NONINFRINGEMENT. IN NO EVENT
#   SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLEFOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
#   IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
#   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#   
########################################################################################


library(GenSA)

################################## Function Defs ###############################

# sum of three Gaussian functions
threegauss = function(x, params) {
  res = (params[3]*exp(-1/2*(x-params[1])^2/params[2]^2) +
           params[6]*exp(-1/2*(x-params[4])^2/params[5]^2) +
           params[9]*exp(-1/2*(x-params[7])^2/params[8]^2))
  res  
}

# single Gaussian
agauss = function(x, params) {
  res = (params[3]*exp(-1/2*(x-params[1])^2/params[2]^2))
  res
}

# plateau with Gaussian fade-out
gaussline = function(x, params){
  res = c()
  for(xi in x){
    if (xi < params[1]){
      y = (params[3]*exp(-1/2*(xi-params[1])^2/params[4]^2))
    }
    else if (xi > params[2]) {
      y=  (params[3]*exp(-1/2*(xi-params[2])^2/params[4]^2))
    }
    else{
      y =( params[3])
    }
    res = c(res, y)
  }
  res
}

# sum of two Gaussian peaks plus a constant plateau with Gaussian fade-out
twogaussplusbroadendgauss = function(x, params) {
  res = (params[3]*exp(-1/2*(x-params[1])^2/params[2]^2) +
           params[6]*exp(-1/2*(x-params[4])^2/params[5]^2) +
           gaussline(x, params[7:10]))
  res
}

# guess plataeu height as minimum of vals between indices x1 and x2
guessplateau = function(vals, x1,x2) {
  pg = vals[x1-1 + which.min(vals[c(x1):c(x2)])]
  dif = vals[x2] - pg
  
  # x2 is the minimum, take a slightly lower value as the guess
  if(dif <= 0){
    pg = vals[x2]*0.99999999
  }
  pg
}

# integrate peak area of G1, S, G2 and sum of them
integratepeaks = function(xs, params){
  G1 = integrate(function(x) {agauss(x, params[1:3])}, lower=min(xs), upper=max(xs))
  G2 = integrate(function(x) {agauss(x, params[4:6])}, lower=min(xs), upper=max(xs))
  S = integrate(function(x) {gaussline(x, params[7:10])}, lower=min(xs), upper=,max(xs))
  Pop = G1$value+S$value+G2$value
  res = c(G1$value,S$value,G2$value,Pop)
  res
}



################ Fitting script #########################

# for saving the results of multiple runs
# do not re-initialize every time
curveareas = data.frame()


########### START ##################

## 1) the integrated DAPI intensities per cell should be assigned to DAPIIntensities

DAPIIntensities = c()

## in our case:
# platerow = 3
# platecolumn = 3
# DAPIIntensities = subset(SKBR3,Row %in% platerow & Column == platecolumn)$IDapiAbs


## 2) get xs and ys
# xs = histogram bin mids
# ys = density curve
h = hist(DAPIIntensities, breaks = 100, probability = T)
xs = h$mids
ys = h$density



# 3) Guess parameters
# Guesses for peak positions
xG1 = xs[which.max(smooth(ys))]
xp1g = which.max(ys)
pgxg = xp1g - 1 + which.min(ys[xp1g:(which.max(ys)*2.5)])
xp2g = pgxg - 1 + which.max(ys[pgxg:100])
xG2 = xs[xp2g]

# Guesses for peak and plateau heights
peak1guess = ys[which.max(smooth(ys))]
plateauguess <- guessplateau(ys,xp1g,xp2g)
peak2guess = ys[xp2g]

# Guess for all 10 params
guess = c(xG1, xG1/6, (peak1guess-plateauguess),
          xG2, xG1/6, (peak2guess-plateauguess),
          xG1, xG2, plateauguess, xG1/6)

# lower and upper bounds
guessupper = guess + guess*c(0.5, 0.9, 0.05,
                             0.1, 0.5, 0.5,
                             0.25, 0.15, 0.1, 0.9)
guesslower = guess - guess*c(0.5, 0.1, 0.3,
                             0.05, 0.5, 0.1,
                             0.5, 0.05, 0.9, 0.1)


## 4) do optimization
res = GenSA(NULL,function(guess) {sum((ys-twogaussplusbroadendgauss(xs, guess))^2)},
             upper = guessupper, lower = guesslower )


## 5) plot histogram + fitted curves

plot(xs, ys, type = "h")
  lines(xs, agauss(xs, res$par[1:3]), col="red")
  lines(xs, agauss(xs, res$par[4:6]), col="green")
  lines(xs, gaussline(xs, res$par[7:10]), col="blue")
  lines(xs, twogaussplusbroadendgauss(xs, res$par), col="cyan")
  
  # plot guesses
  abline(v = xG1, col='red')
  abline(v = xG2,col ='green')
  abline(h = peak1guess, col='red')
  abline(h = plateauguess, col='blue')
  abline(h = peak2guess, col='green')
  


## 6) append resulting areas under curves to accumulated results

### name for the condition of the histogram
# e.g.
 
condition = '16nM Trastuzumab'

curveareas = rbind.data.frame(curveareas, c(condition, as.numeric(integratepeaks(xs, res$par))), stringsAsFactors = FALSE)
colnames(curveareas) = c('condition', 'G1', 'S', 'G2', 'total')


################ Accumulated results plot/saving ###################

# normalize to area=1
normIntegrals = apply(curveareas[,2:5], 2, function(x){as.numeric(as.character(x))}) /
  as.numeric(as.character(curveareas[,5]))
# plot stacked bars
barplot(t(normIntegrals[,1:3]), names.arg = curveareas$condition)

# save cell cycle fractions
write.table(normIntegrals, "Z:/normIntegrals.txt", sep="\t") 



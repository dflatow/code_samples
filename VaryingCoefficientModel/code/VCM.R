
require(stats); require(graphics)
library(splines)

# helper function to repeat row
rep.row<-function(x,n){matrix(rep(x,each=n),nrow=n)}

# helper function to add confidence bands to plots
add_confidence_band = function(index, beta, cov, x, color){
  scaled_cov = x %*% cov %*% t(x)
  upper = beta + 2 * sqrt(diag(scaled_cov))
  lower = beta - 2 * sqrt(diag(scaled_cov))
  polygon( c(index, rev(index)) , c(upper , rev(lower)), 
           col=color, border=NA)
}

# degrees of freedom for natural cubic spline bases
df = 7
n = 500

data = read.table('../data/vcdata.csv', sep=',', header=TRUE)
" We estimate a linear model in the natural cubic spline basis
  of time multiplied with each of our variables. Since we have 
  3 variables (intercept, x1, x2) we will have (df * 3) basis
  functions where df is the degrees of freedom of each NCS basis.

  For example the second term in the linear model below is the
  NCS basis (df functions of t) of time each multiplied by x1. 
  This term (expanded) will result in df terms in the linear 
  model and we will have to estimate df coefficients for each of
  these bases."

t_basis =  ns(data$t, df=df, intercept = T) 
m = lm(y ~ 0 + t_basis + t_basis:x1 + t_basis:x2, data=data)
results = summary(m)

# coefficients for the intercept bases
int = results$coefficients[1:df, 1]

# coefficients for the beta_1 bases
b1 =  results$coefficients[(df+1):(2*df), 1]

# coefficients for the beta_2 bases
b2 =  results$coefficients[(2*df + 1):(3*df), 1]

# back out betas as a function of time
intercept = t_basis %*% int
beta_1 = t_basis %*% b1
beta_2 = t_basis %*% b2

# plot coefficient time series estimates
matplot(1:n, intercept, type="l",col="red", 
        ylim=range(-1.5, 3.5), 
        main="Time Varying Coefficients",
        xlab="t", ylab="coefficient value")

matplot(1:n, beta_1, type="l",col="blue", add=T)
matplot(1:n, beta_2, type="l",col="green", add=T)
legend("bottomleft", inset=.05, legend=c("intercept", "b1", "b2"), 
       pch="-", col=c("red", "blue", "green"), horiz=TRUE)


# add confidence band to intercept coefficient time series
int_cov = vcov(results)[1:df, 1:df]
add_confidence_band(1:n, intercept, int_cov, 
                    t_basis, rgb(1, 0, 0, 0.3))

# add confidence band to x1 coefficient time series
beta_1_cov = vcov(results)[(df + 1):(2*df),
                                  (df + 1):(2*df)]
add_confidence_band(1:n, beta_1, int_cov, 
                    t_basis, rgb(0, 0, 1, 0.3))

# add confidence band to x2 coefficient time series
beta_2_cov = vcov(results)[(2*df + 1):(3*df),
                                  (2*df + 1):(3*df)]
add_confidence_band(1:n, beta_2, beta_2_cov, 
                    t_basis, rgb(0, 1, 0, 0.3))
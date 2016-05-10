\name{sim_slpb}
\alias{sim_slpb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Social learning and parasitic behavior simulations}
\description{
  This function runs the simulations for both models from the paper. These simulations use infinite populations and stochastic environmental changes.
}
\usage{
sim_slpb(tmax=200,u=0.01,b=2,c=1,v=0.5,w=2,s=1,w0=10,mu=1e-4,
    p=c(0,0,0.99),q=0,r=0.01,d=0.5,e=0,h=0.5, 
    sto=TRUE , evolve=TRUE , no_r_before=0 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tmax}{Number of generations to run the simulation}
  \item{u}{Probability of environmental change}
  \item{b}{Fitness benefit of adaptive behavior}
  \item{c}{Fitness cost of innovation}
  \item{v}{Fitness cost of parasitic behavior}
  \item{w}{Transmission advantage of parasitic behavior}
  \item{s}{Probability of successful innovation}
  \item{w0}{Baseline fitness}
  \item{mu}{Probability innovation results in parasitic behavior}
  \item{p}{Vector of initial frequencies for heritable learning strategies: pure individual, pure social, critical social (in that order).}
  \item{q}{Initial frequency of adaptive behavior}
  \item{r}{Initial frequency of parasitic behavior}
  \item{d}{Probability critical social learner detects non-adaptive behavior}
  \item{e}{Probability critical social learner rejects adaptive behavior}
  \item{h}{Fitness cost of critical social learning}
  \item{sto}{Logical flag for stochastic environmental transitions}
  \item{evolve}{Logical flag to enable evolutionary dynamics of learning strategies}
  \item{no_r_before}{Integer for first generation in which parasitic behavior is allowed to appear via innovation rate \code{mu}.}
}
\details{
  This command runs a simulation of the parasitic behavior model that allows all three learning strategies: pure individual, pure social, and critical social learning. The code uses the analytical recursions derived in the paper. However, it uses stochastic environmental transitions, and so it is sufficient to explore geometric mean fitness effects that the analytical solutions might miss.

  See the examples further down this page for code sufficient to reproduce each figure from the paper.
}
\value{
    Returns a \code{list} with the following named slots.
    \item{x}{Array of frequencies of each combination of adaptive and parasitic behavior in each generation. \code{x[1,1,t]} is frequency of adaptive/parasitic behavior in generation \code{t}.}
    \item{pars}{Parameters passed to initialize simulation}
    \item{p}{Vector of frequencies of critical social learning in each generation}
    \item{pS}{Vector of frequencies of pure social learning in each generation}
    \item{q}{Vector of frequencies of adaptive behavior in each generation}
    \item{r}{Vector of frequencies of parasitic behavior in each generation}
    \item{k}{Vector of covariances of adaptive and parasitic behavior in each generation}
    \item{W_hist}{Vector of mean population fitness in each generation}
}
\references{}
\examples{
###########
# Figure 1

# frequency plot
set.seed(1)
x <- sim_slpb( tmax=200 , u=0.1 , s=0.6 , d=0 , e=0 , h=0 , v=0.5 , w=2 , mu=1e-3 , b=3 , c=1 , r=0 , w0=10 , evolve=TRUE , sto=TRUE , p=c(0.09,0.91,0) , q=0.3 )
plotsim(x , show_alleles=c(NA,"black",NA) , showk=FALSE , yaxp=c(0,1,2) )

# correlation plot
rho <- x$k / sqrt( x$r*(1-x$r) * x$q * (1-x$q) )
plot( rho , type="l" , ylab="correlation" , bty="n" , xaxt="n" , xlab="" , xlim=c(1,floor(200*1.05)) , ylim=c(min(rho,na.rm=TRUE),0.01) , yaxt="n" )
axis( 2 , at=c(0,-0.6) , labels=c(0,-0.6) )
abline( h=0 , lty=2 , lwd=0.5 )

###########
# Figure 2

# frequency plot
set.seed(1)
x <- sim_slpb( tmax=200 , u=0.2 , s=0.5 , d=0.5 , e=0 , h=1/20 , v=1/2 , w=1.2 , mu=1e-2 , b=3 , c=1 , r=0 , w0=10 , evolve=TRUE , sto=TRUE , p=c(0,0,0.99) )
plotsim(x , showk=FALSE , yaxp=c(0,1,2) )

sols <- comp_sols(x,1)
points( floor(200*1.025) , sols$q , pch=16 , col="blue" )
points( floor(200*1.025) , sols$r , pch=16 , col="red" )

# fitness plot
WI <- with(x$pars,{(w0+s*b-c)})
plot( x$W_hist/WI , type="l" , ylab="fitness" , bty="n" , xaxt="n" , xlab="" , xlim=c(1,floor(200*1.05)) , yaxt="n" )
axis( 2 , at=c(2,1,0.5) , labels=c(2,1,0.5) )
abline( h=1 , lwd=0.5 , lty=2 )

###########
# Figure 3

# frequency plot
set.seed(1)
x <- sim_slpb( tmax=200 , u=0.2 , s=0.5 , d=0.3 , e=0 , h=1/20 , v=1/2 , w=2 , mu=1e-4 , b=2 , c=1 , r=0 , w0=10 , evolve=TRUE , sto=TRUE , p=c(0,0,0.99) )
plotsim(x , showk=FALSE , yaxp=c(0,1,2) )

sols <- comp_sols(x,2)
points( floor(200*1.025) , sols$q , pch=16 , col="blue" )
points( floor(200*1.025) , sols$r , pch=16 , col="red" )

# fitness plot
WI <- with(x$pars,{(w0+s*b-c)})
plot( x$W_hist/WI , type="l" , ylab="fitness" , bty="n" , xaxt="n" , xlab="" , xlim=c(1,floor(200*1.05)) , yaxt="n" )
axis( 2 , at=c(2,1,0.5) , labels=c(2,1,0.5) )
abline( h=1 , lwd=0.5 , lty=2 )

# correlation plot
rho <- x$k / sqrt( x$r*(1-x$r) * x$q * (1-x$q) )
plot( rho , type="l" , ylab="correlation" , bty="n" , xaxt="n" , xlab="" , xlim=c(1,floor(200*1.05)) , ylim=c(min(rho,na.rm=TRUE),0.01) , yaxt="n" )
axis( 2 , at=c(0,-0.6) , labels=c(0,-0.6) )
abline( h=0 , lty=2 , lwd=0.5 )
Erho <- sols$rho
points( floor(200*1.025) , Erho , pch=16 )

###########
# Figure 4

# frequency plot
set.seed(1)
x <- sim_slpb( tmax=200 , u=0.05 , s=0.5 , d=0.75 , e=0.01 , h=1/10 , v=5 , w=2 , mu=1e-4*1 , b=10 , c=1 , r=0 , w0=10 , evolve=TRUE , sto=TRUE , p=c(0,0,0.99) , no_r_before=20 )
plotsim(x , showk=FALSE , yaxp=c(0,1,2) )

# fitness plot
WI <- with(x$pars,{(w0+s*b-c)})
plot( x$W_hist/WI , type="l" , ylab="fitness" , bty="n" , xaxt="n" , xlab="" , xlim=c(1,floor(200*1.05)) , yaxt="n" )
axis( 2 , at=c(2,1,0.5) , labels=c(2,1,0.5) )
abline( h=1 , lwd=0.5 , lty=2 )

# correlation plot
rho <- x$k / sqrt( x$r*(1-x$r) * x$q * (1-x$q) )
blank(w=1.66,h=0.4)
plot( rho , type="l" , ylab="correlation" , bty="n" , xaxt="n" , xlab="" , xlim=c(1,floor(200*1.05)) )
abline( h=0 , lty=2 , lwd=0.5 )

###########
# Figure 5

# frequency plot
set.seed(1)
x <- sim_slpb( tmax=200 , u=0.1 , s=0.1 , d=0.75 , e=0.1 , h=1/10 , v=2 , w=2 , mu=1e-4*1 , b=20 , c=1 , r=0 , w0=10 , evolve=TRUE , sto=TRUE , p=c(0,0,0.99) , no_r_before=50 )
plotsim(x , showk=FALSE , yaxp=c(0,1,2) )

# fitness plot
WI <- with(x$pars,{(w0+s*b-c)})
plot( x$W_hist/WI , type="l" , ylab="fitness" , bty="n" , xaxt="n" , xlab="" , xlim=c(1,floor(200*1.05)) , yaxt="n" )
axis( 2 , at=c(2,1,0.5) , labels=c(2,1,0.5) )
abline( h=1 , lwd=0.5 , lty=2 )

# correlation plot
rho <- x$k / sqrt( x$r*(1-x$r) * x$q * (1-x$q) )
blank(w=1.66,h=0.4)
plot( rho , type="l" , ylab="correlation" , bty="n" , xaxt="n" , xlab="" , xlim=c(1,floor(200*1.05)) , las=1 , yaxp=c(0,0.8,1) )
abline( h=0 , lty=2 , lwd=0.5 )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }


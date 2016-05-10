Parasitic behavior simulation code
==========

This R package contains the code necessary to replicate the simulations in the paper ``Parasitic behavior may hurt or help social learners''. These are stochastic, numerical simulations of the analytical model presented in the paper's appendix.

The package can be quickly installed by using ``dev_tools``:
```
install.packages('devtools')
library(devtools)
install_github("rmcelreath/parasiticbehaviorsim")
```

Load the package and bring up the overview with:
```
library(parasiticbehaviorsim)
?parasiticbehaviorsim
```

See the R help page ``?sim_slpb`` for complete code necessary to replicate each of the figures in the paper. For example, the frequency plot portion of Figure 1 can be replicated with:
```
set.seed(1)
x <- sim_slpb( tmax=200 , u=0.1 , s=0.6 , d=0 , e=0 , h=0 , v=0.5 , w=2 , mu=1e-3 , b=3 , c=1 , r=0 , w0=10 , evolve=TRUE , sto=TRUE , p=c(0.09,0.91,0) , q=0.3 )
plotsim(x , show_alleles=c(NA,"black",NA) , showk=FALSE , yaxp=c(0,1,2) )
```

# master simulation for both models

sim_slpb <- function(tmax=200,u=0.01,b=2,c=1,v=0.5,w=2,s=1,w0=10,mu=1e-4,p=c(0,0,0.99),q=0,r=0.01,
    d=0.5,e=0,h=0.5, 
    sto=TRUE , evolve=TRUE , no_r_before=0 ) {
    x <- array(NA,dim=c(tmax,2,2))

    # take p input argument and translate to allele frequencies
    p_hist <- rep(NA,tmax)
    pS_hist <- rep(NA,tmax)
    p_input <- p/sum(p)
    p <- p_input[3]
    pS <- p_input[2]
    p_hist[1] <- p
    pS_hist[1] <- pS
    
    k <- 0
    x[1,1,1] <- q*r + k
    x[1,1,2] <- q*(1-r) - k
    x[1,2,1] <- (1-q)*r - k
    x[1,2,2] <- (1-q)*(1-r) + k
    W_hist <- rep(NA,tmax)
    W_hist[1] <- w0
    mu_save <- mu
    for ( i in 2:tmax ) {
        if ( i < no_r_before ) mu <- 0 else mu <- mu_save
        ut <- sample(c(0,1),size=1,prob=c(1-u,u))
        if ( sto==FALSE ) ut <- u
        r <- x[i-1,1,1] + x[i-1,2,1] # freq parasites
        R <- r*w/(r*w+1-r)
        Q11 <- (1-ut)*ifelse( r>0 , R*x[i-1,1,1]/r , 0 )
        Q10 <- (1-ut)*ifelse( r==1 , 0 , (1-R)*x[i-1,1,2]/(1-r) )
        Q01 <- (1-ut)*ifelse( r>0 , R*x[i-1,2,1]/r , 0 ) + ut*R
        Q00 <- 1 - Q11 - Q01 - Q10
        Q <- Q11 + Q10
        x[i,1,1] <- (1-p-pS)*s*mu + p*((1-e)*Q11 + Q*e*s*mu + (1-Q)*d*s*mu) + pS*Q11
        x[i,1,2] <- (1-p-pS)*s*(1-mu) + p*((1-e)*Q10 + Q*e*s*(1-mu) + (1-Q)*d*s*(1-mu)) + pS*Q10
        x[i,2,1] <- (1-p-pS)*(1-s)*mu + p*((1-d)*Q01 + Q*e*(1-s)*mu + (1-Q)*d*(1-s)*mu) + pS*Q01
        x[i,2,2] <- 1 - x[i,1,1] - x[i,1,2] - x[i,2,1]
        WI <- w0 + s*b - mu*v - c;
        WS <- w0 + Q10*b + Q11*(b-v) - Q01*v
        WC <- w0 + Q10*b + Q11*(b-v) + Q00*(d*(s*b-c-v*mu)) + Q01*(d*(s*b-c-v*mu)-(1-d)*v) - h
        if ( evolve==TRUE ) {
            WBAR <- (p*WC + (1-p-pS)*WI + pS*WS)
            p <- p * WC / WBAR
            pS <- pS * WS / WBAR
            p_hist[i] <- p
            pS_hist[i] <- pS
            W_hist[i] <- WBAR
        }
    }
    q_hist <- sapply(1:tmax,function(i) x[i,1,1]+x[i,1,2] )
    r_hist <- sapply(1:tmax,function(i) x[i,1,1]+x[i,2,1] )
    k_hist <- sapply( 1:tmax , function(i) x[i,1,1]*x[i,2,2] - x[i,2,1]*x[i,1,2] )

    pars <- list(
        b=b,
        c=c,
        s=s,
        u=u,
        v=v,
        w=w,
        d=d,
        e=e,
        h=h,
        w0=w0,
        mu=mu,
        start_r=r,
        sto=sto
    )

    return(list(x=x,pars=pars,p=p_hist,pS=pS_hist,q=q_hist,r=r_hist,k=k_hist,W_hist=W_hist))
}

comp_sols <- function(x,sol=3) {
    if ( sol==1 ) {
        qhat <- with( x$pars , {(d*s)/(u + d*(s - s*u))} )
        rhat <- 0
        khat <- 0
        eigenvals <- with( x$pars , { c(
                (1-u)*(1-e*(1-s)-d*s),
                w*(1-d),
                w*(1-u)*(1-e)
            ) } )
    }
    if ( sol==2 ) {
        qhat <- with( x$pars , {(((-1 + d)*d*s*w)/((-1 + d*s)*(-1 + u) + (-1 + d)*w))} )
        rhat <- with( x$pars , {(-1 + w - d*w)/(-1 + w)} )
        khat <- with( x$pars , {((-1 + d)*d*s*w*(1 + (-1 + d)*w))/((-1 + w)*((-1 + d*s)*(-1 + u) + (-1 + d)*w))} )
        eigenvals <- with( x$pars , { c(
                1/(w*(1-d)),
                (1-u)*(1-e*(1-s)-d*s)/(w*(1-d)),
                (1-u)*(1-e)/(1-d)
            ) } )
    }
    if ( sol==3 ) {
        qhat <- with( x$pars , {(-(-1 + d*s)*(d - u) + (-1 + u)*(d - u + d*s*u)*w)/(d*(-1 + u)*(-1 + d*s + w))} )
        rhat <- with( x$pars , {(-1 + w - u*w)/(-1 + w)} )
        khat <- with( x$pars , {((d*(1 + s*(-1 + u)) - u)*u*w*(1 + (-1 + u)*w))/(d*(-1 + u)*(-1 + w)*(-1 + d*s + w))} )
        eigenvals <- with( x$pars , { c(
                1/(w*(1-e-u*(1-e))),
                (1-e-d*s+e*s)/(w*(1-e)),
                (1-d)/((1-e)*(1-u))
            ) } )
    }
    rho <- khat / sqrt( qhat*(1-qhat)*rhat*(1-rhat) )
    return(list(q=qhat,r=rhat,k=khat,rho=rho,stab=!any(eigenvals>1),ev=eigenvals))
}

plotsim <- function(x,showtext=FALSE,showk=TRUE,show_alleles=c(NA,NA,NA),...) {
    decn <- 3
    n <- length(x$p)
    offy <- 0.035
    plot( NULL , ylim=c(0,1) , ylab="frequency" , xlab="generation" , xlim=c(1,floor(n*1.05)) , xaxt="n" , bty="l" , ... )
    axis( 1 , at=c(1,floor(n/2),n) , labels=c(1,floor(n/2),n) )
    
    if ( showk==TRUE ) abline( h=0.5 , lwd=0.5 , col="gray" )
    lines( 1:n , x$q , col="blue" )
    lines( 1:n , x$r , col="red" )
    abline( h=1 , lty=2 , lwd=0.5 )
    abline( h=0 , lty=2 , lwd=0.5 )
    if ( !is.na(show_alleles[3]) ) lines( 1:n , x$p , col=show_alleles[3] )
    if ( !is.na(show_alleles[2]) ) lines( 1:n , x$pS , col=show_alleles[2] ) 
    if ( !is.na(show_alleles[1]) ) lines( 1:n , 1-x$p-x$pS , col=show_alleles[1] ) 
    if ( showk==TRUE ) lines( 1:n , x$k+0.5 , lwd=0.5 )
}


if ( FALSE ) {

##########################
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

############################
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

############################
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

############################
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

############################
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
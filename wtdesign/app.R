
library(shiny)

ui<-(navbarPage("Early-phase design based on toxicity and activity endpoints",
                #title
                tabPanel("Simulation", 
                         h3("Web Application for simulating operating characteristics of the Phase I/II method of Wages and Tait (2015)"),
                         #h4("Nolan A. Wages and Gina R. Petroni"),
                         h4("Division of Translational Research & Applied Statistics, University of Virginia;", a("nwages@virginia.edu", href="mailto:nwages@virginia.edu")),
                         br(),
                         h5("1. Enter an assumed set of true DLT probabilities, separated by commas.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         #input: boxes for truth input values OR text field
                         textInput("p0","True DLT probabilities",value="0.08,0.12,0.20,0.30,0.40",width=NULL),
                         
                         br(),
                         h5("2. Enter an assumed set of true activity probabilities, separated by commas.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         #input: tul
                         textInput("q0","True activity probabilities",value="0.20,0.40,0.60,0.60,0.60",width=NULL),
                         
                         br(),
                         h5("3. Enter the target DLT probability that defines the MTD for the study."),
                         #input: tul
                         numericInput("tul","Target DLT rate",0.30,0,1,0.01,NULL),
                         
                         br(),
                         h5("4. Enter the lower limit for the probability of activity that defines the threshold for futility"),
                         #input: ell
                         numericInput("ell","Minimum acceptable activity rate",0.20,0,1,0.01,NULL),
                         
                         br(),
                         h5("5. Enter the cohort size required before the next model-based update. Cohort size may be 1, 2, or 3 patients."),
                         #input: cohortsize
                         numericInput("cohortsize","Cohort size",2,1,6,1,NULL),  
                         
                         br(),
                         h5("4. Enter the maximum sample size for the study. This number should be a multiple of the cohort size entered above."),
                         #input: ncohort
                         numericInput("ssize","Maximum number of patients",30,0,NA,1,NULL),
                         
                         br(),
                         h5("7. Enter the number of simulations. A minimum of 1000 is recommended."),
                         #input: ntrial
                         numericInput("ntrial","Number of simulations", 10, 0, 100000, NA, NULL),
                         
                         br(),
                         p("8. Enter the index of the starting dose level.", strong("Note:"), "Index of", span("lowest", style = "color:blue"), "dose level is always 1. If the design allows for", em("'minus'"), "
                           dose levels (i.e. -2, -1, etc.), then the index of the starting dose should account for these lower levels (i.e. if -1 dose level allowed, index of starting dose is 2.)"),
                         #input: start.comb
                         numericInput("start.comb","Index of starting dose level",1,1,NA,1,NULL),
                         
                         
                         br(),
                         h5("9. Enter the total number of patients treated on any dose required to stop the trial. At any point in the trial, if the recommendation is to assign
                            the next cohort to a dose that already has the entered number of patients treated on the dose, the study
                            is stopped and the recommended dose is declared the optimal dose. If the entered number is larger than the maximum sample size, each trial will accrue to
                            the maximum sample size."),
                         #input: n
                         numericInput("n.stop","Number of patients needed on one dose to STOP",50,0,NA,1,NULL),
                         
                         br(),
                         h5("10. Set the seed of the random number generator."),
                         #input: rseed
                         numericInput("rseed","Random Seed",580,0,NA,1,NULL),
                         
                         #input: submit button
                         submitButton("Run simulation study",icon("flask"),NULL),
                         
                         #output: simbpocrm(p0,q0,tul,ell,psi,cohortsize,ncohort,ntrial,start.comb,n.stop,rseed)
                         verbatimTextOutput("simulation"),
                         
                         mainPanel(
                           p(strong("This application simulates operating characteristics for the early-phase method of Wages and Tait [1] with the following practical modifications.")),
                           p("1.", strong("Dose-finding algorithm:"), "The trial does not skip over an untried dose level."),
                           #p("2.", strong("Safety:"), "For the specified target DLT rate and total number of dose levels, the skeleton of power model d^exp(a) is generated according to Lee and Cheung (2009) [2] using a prior MTD located at the median dose level and a spacing measure of delta=0.05."),
                           #p("3.", strong("Prior:"), "The prior distribution on the parameter", strong("a"), "is the least informative, mean zero normal prior according to Lee and Cheung (2011) [3]."),
                           #p("4.", strong("Accuracy Index:"), "The Accuracy Index is given in Equation 6.1 in Cheung [3]."),
                           p("2.", strong("Safety Stopping Rule:"), "Stop the trial for safety if the lower limit of a 90% probability interval exceeds the target DLT rate."),
                           p("3.", strong("Futility Stopping Rule:"), "Stop the trial for futility if the upper limit of a 90% probability interval is lower than the minimum acceptable efficacy rate."),
                           strong("References:"),
                           p("[1] Wages NA, Tait C (2015). Seamless phase I/II adaptive design for oncology trials of molecularly targeted agents, ", em("Journal of Biopharmaceutical Statistics;"), strong("25:"),"903-920.")
                          
                         ),
                         tags$style(type="text/css",
                                    ".shiny-output-error { visibility: visible; }",
                                    ".shiny-output-error:before { visibility: visible; }"
                         )
                ),
                #title
                tabPanel("Implementation", 
                         h3("Web Application for implementation of the Phase I/II method of Wages and Tait (2015)"),
                         # h4("Nolan A. Wages and Gina R. Petroni"),
                         h4("Division of Translational Research & Applied Statistics, University of Virginia;", a("nwages@virginia.edu", href="mailto:nwages@virginia.edu")),
                         #h5("1. Enter an assumed set of true DLT probabilities, separated by commas. The length of this set should be equal to the number of dose levels."),
                         #input: start
                         br(),
                         h4(strong("Design / Protocol Information")),
                         
                         #h5("1. Enter the maximum sample size for the study. This number should be a multiple of the cohort size entered above."),
                         ##input: maximum sample size
                         #numericInput("maxn","Maximum number of patients",30,0,NA,1,NULL),
                         
                         #br(),
                         h5("1. Enter the target DLT rate probability that defines the MTD for the study."),
                         #input: theta
                         numericInput("theta","Target DLT rate",0.30,0,1,0.01,NULL),
                         
                         br(),
                         h5("2. Enter the lower limit for the probability of activity that defines the threshold for futility"),
                         #input: phi
                         numericInput("phi","Minimum acceptable activity rate",0.20,0,1,0.01,NULL),
                         
                         
                         br(),
                         h4(strong("Observed Trial Data (do not count 'replaced' patients)")),
                         #br(),
                         #p("1. Enter the index of the starting dose level.", strong("Note:"), "Index of", span("lowest", style = "color:blue"), "dose level is always 1. If the design allows for", em("'minus'"), "
                        #   dose levels (i.e. -2, -1, etc.), then the index of the starting dose should account for these lower levels (i.e. if -1 dose level allowed, index of starting dose is 2.)"),
                        # #input: start.comb
                        # numericInput("begin","Index of starting dose level",1,1,NA,1,NULL),
                         
                        # br(),
                         #input: boxes for truth input values OR text field
                         h5("1. Enter number of observed DLTs at each dose level. If none have been observed or a dose level has not yet been tried, enter '0'.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         textInput("y","Number of observed DLTs at each dose level (if none, enter '0')",value="0,0,1,1,2",width=NULL),
                         
                         br(),
                         h5("2. Enter the number of patients evaluated for DLT at each dose level. If a dose level has not yet been tried, enter '0'.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         #input: boxes for truth input values OR text field
                         textInput("ny","Number of patients evaluated for DLT on each dose level (if none, enter '0')",value="3,3,3,3,3",width=NULL),
                         
                         br(),
                         h5("3. Enter number of observed activity responses at each dose level. If none have been observed or a dose level has not yet been tried, enter '0'.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         #input: boxes for truth input values OR text field
                         textInput("z","Number of observed activity responses at each dose level (if none, enter '0')",value="0,0,1,1,2",width=NULL),
                         
                         br(),
                         h5("4. Enter the number of patients evaluated for response at each dose level. If a dose level has not yet been tried, enter '0'.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         #input: boxes for truth input values OR text field
                         textInput("nz","Number of patients evaluated for activity on each dose level (if none, enter '0')",value="3,3,3,3,3",width=NULL),
                         
                         br(),
                         h5("5. Enter the most recent dose level administered in the study."),
                         #input: comb.curr
                         numericInput("dose.curr","Current dose level",3,1,NA,1,NULL),
                         
                         #input: submit button
                         submitButton("Get updated recommended dose level",icon("flask"),NULL),
                         
                         #output: impcrm(y,n,target)
                         verbatimTextOutput("implementation"),
                         
                         mainPanel(
                           p(strong("This application computes a recommended dose level for the next patient in an early-phase trial according to the method of Wages and Tait [1] with the following practical modifications.")),
                           p("1.", strong("Dose-finding algorithm:"), "The trial does not skip over an untried dose level."),
                           #p("2.", strong("Safety:"), "For the specified target DLT rate and total number of dose levels, the skeleton of power model d^exp(a) is generated according to Lee and Cheung (2009) [2] using a prior MTD located at the median dose level and a spacing measure of delta=0.05."),
                           #p("3.", strong("Prior:"), "The prior distribution on the parameter", strong("a"), "is the least informative, mean zero normal prior according to Lee and Cheung (2011) [3]."),
                           #p("4.", strong("Accuracy Index:"), "The Accuracy Index is given in Equation 6.1 in Cheung [3]."),
                           p("2.", strong("Safety Stopping Rule:"), "Stop the trial for safety if the lower limit of a 90% probability interval exceeds the target DLT rate."),
                           p("3.", strong("Futility Stopping Rule:"), "Stop the trial for futility if the upper limit of a 90% probability interval is lower than the minimum acceptable efficacy rate."),
                           strong("References:"),
                           p("[1] Wages NA, Tait C (2015). Seamless phase I/II adaptive design for oncology trials of molecularly targeted agents, ", em("Journal of Biopharmaceutical Statistics;"), strong("25:"),"903-920.")
                           
                         ),
                        
                        tags$style(type="text/css",
                                   ".shiny-output-error { visibility: visible; }",
                                   ".shiny-output-error:before { visibility: visible; }"
                        )
                         
                         
                )
)

)

server <- function(input, output) {
  simbpocrm<-function(p0,q0,tul,ell,cohortsize,ssize,ntrial,start.comb,n.stop,rseed){
    set.seed(rseed)
    library(dfcrm)
    library(nnet)
    d=ncomb=length(p0)
    #delta=0.04
    p.skel<-getprior(halfwidth=0.05, target=tul, nu=ifelse(d %% 2 == 0,d/2,(d+1)/2), nlevel=d, model="empiric")
      #getprior(halfwidth=0.05, target=tul, start.comb, nlevel=d, model="empiric")
      
    ncohort=ssize/cohortsize
    q.skel<-matrix(nrow=2*d-1,ncol=d)
    q.skel[1,]<-getprior(0.04,0.70,d,d)
    reverse=sort(q.skel[1,],decreasing=TRUE)
    for(i in 2:d){
      q.skel[i,]<-c(tail(q.skel[1,],d-i+1),reverse[2:i])
    }
    for(k in 2:d){
      q.skel[k+d-1,]<-c(tail(q.skel[1,],d-k+1),rep(0.7,k-1))
    }
    
    
   # q.skel<-matrix(nrow=d,ncol=d)
   #  q.skel[1,]<-seq(seq(0,1,length.out=d)[2],seq(0,1,length.out=d)[d-1],length.out=d)
   # for(k in 2:d){
   #    q.skel[k,]<-c(tail(q.skel[1,],d-k+1),rep(max(q.skel[1,]),k-1))
   # }


    ###Load the function 'bpocrm' 
    bpocrm<-function(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,n.stop,start.comb){
      
      
      
      
      # if a single ordering is inputed as a vector, convert it to a matrix
      if(is.vector(q.skel)) q.skel=t(as.matrix(q.skel));
      
      nord.eff = nrow(q.skel);
      mprior.eff = rep(1/nord.eff, nord.eff); # prior for each efficacy ordering
      
      
      avar=1.34
      bcrmh<-function(a,p,y,n){
        s2=avar
        lik=exp(-0.5*a*a/s2)
        for(j in 1:length(p)){
          pj=p[j]**exp(a)
          lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
        }
        return(lik);
      }
      
      bcrmht<-function(a,p,y,n){
        s2=avar
        lik=a*exp(-0.5*a*a/s2)
        for(j in 1:length(p)){
          pj=p[j]**exp(a)
          lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
        }
        return(lik);
      }
      
      bcrmht2<-function(a,p,y,n){
        s2=avar
        lik=a^2*exp(-0.5*a*a/s2)
        for(j in 1:length(p)){
          pj=p[j]**exp(a)
          lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
        }
        return(lik);
      }
      
      
      
      ### run a trial 	
      ncomb = length(p.skel);   #number of combos
      y=rep(0,ncomb);
      z=rep(0,ncomb);
      n=rep(0,ncomb);  #number of toxicity/responses/patients at each dose level
      comb.curr = 1#start.comb;  # current dose level	 
      ptox.hat = numeric(ncomb); # estimate of toxicity/efficacy prob
      peff.hat = numeric(ncomb); # estimate of efficacy prob
      p.safety=numeric(ncomb)
      p.futility=numeric(ncomb)
      dose.select=rep(0,ncomb); # a vector of indicators for dose selection
      stop=stops=stopf=0; #indicate if trial stops early
      i=1	
      while(i <= ncohort)
      {
        # generate data for a new cohort of patients
        y[comb.curr] = y[comb.curr] + rbinom(1,cohortsize,p0[comb.curr]);
        z[comb.curr] = z[comb.curr] + rbinom(1,cohortsize,q0[comb.curr]);
        n[comb.curr] = n[comb.curr] + cohortsize;		
        
        if(any(n>n.stop)){
          stop<-0
          break
        }
        
        
        marginal= integrate(bcrmh,lower=-Inf,upper=Inf, p=p.skel, y=y,n=n,abs.tol = 0)$value;
        est=integrate(bcrmht,lower=-10,upper=10, p.skel, y, n,abs.tol = 0)$value/marginal
        e2=integrate(bcrmht2,lower=-10,upper=10, p.skel, y, n,abs.tol = 0)$value/marginal
        
        ptox.hat=p.skel**exp(est)
        post.var=e2-(est)^2
        crit=qnorm(0.5+0.9/2)
        ub=est+crit*sqrt(post.var)
        ptoxL=p.skel**exp(ub)
       
        if(ptoxL[1]>tul){
          stops=1
          break
        }
        
        
        marginal.eff = est.eff=e2.eff=rep(0, nord.eff);
        for(k in 1:nord.eff)
        {
          marginal.eff[k] = integrate(bcrmh,lower=-Inf,upper=Inf, p=q.skel[k,], y=z,n=n,abs.tol = 0)$value;
          est.eff[k]=integrate(bcrmht,lower=-10,upper=10, q.skel[k,], z, n,abs.tol = 0)$value/marginal.eff[k]
          e2.eff[k]=integrate(bcrmht2,lower=-10,upper=10, q.skel[k,], z, n,abs.tol = 0)$value/marginal.eff[k]
        }		
        postprob.eff = (marginal.eff*mprior.eff)/sum(marginal.eff*mprior.eff);
        # toxicity model selection, identify the model with the highest posterior prob
        if(nord.eff>1){ 
          meff.sel = which.is.max(postprob.eff); 
        } else{
          meff.sel = 1;
        }
        peff.hat=q.skel[meff.sel,]**exp(est.eff[meff.sel])
        post.var.eff=e2.eff[meff.sel]-(est.eff[meff.sel])^2
        crit.eff=qnorm(0.5+0.9/2)
        lb.eff=est.eff[meff.sel]-crit.eff*sqrt(post.var.eff)
        peffU=q.skel[meff.sel,]**exp(lb.eff)
      
        
        
        ##dose-finding algorithm
        try=length(n[n>0])	 #highest dose tried
        mtd=which.is.max(-abs(ptox.hat-tul))	
        utility=peff.hat*as.numeric(ptox.hat<=ptox.hat[mtd])
        best=which.max(utility)#ifelse(sum(n)<(2*ncohort*cohortsize/3),sample(1:ncomb,1,prob=utility/sum(utility)),which.max(utility))
        if(peffU[which.max(utility)]<ell){
          stopf=1
          break
        }
        if(comb.curr==ncomb){
          comb.curr=ifelse(best==ncomb,min(mtd,comb.curr),min(mtd,best))
        } else{
          if((best==try)||(best>comb.curr)){comb.curr=min(mtd,comb.curr+1)}
          else if (best<comb.curr) {comb.curr=min(mtd,comb.curr-1)}
          else {comb.curr=min(mtd,comb.curr)}
        }
        
     
        
        i=i+1
        
      }
      if(stop==0 & stops==0 & stopf==0){
        dose.select[comb.curr]=dose.select[comb.curr]+1;
      }
      return(list(dose.select=dose.select,tox.data=y,eff.data=z,pt.allocation=n,stops=stops,stopf=stopf))
    }
    ##########'bpocrm' end here
    
    comb.select<-tox<-eff<-pts<-matrix(nrow=ntrial,ncol=ncomb)
    nstopf=nstops=0
    
    for(i in 1:ntrial){
      result<-bpocrm(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,n.stop,start.comb)
      comb.select[i,]=result$dose.select
      tox[i,]=result$tox.data
      eff[i,]=result$eff.data
      pts[i,]=result$pt.allocation
      nstopf=nstopf+result$stopf
      nstops=nstops+result$stops
    }
    cat("Simulation results: \n\n");
    cat("True DLT probability:                        ", round(p0,3), sep="\t",  "\n");
    cat("True activity probability:                   ", round(q0,3), sep="\t",  "\n");
    cat("Optimal dose selection percentage:           ", formatC(colMeans(comb.select)*100, digits=1, format="f"), sep="\t",  "\n");
    cat("Average number of DLTs:                      ", formatC(colMeans(tox), digits=1, format="f"), sep="\t",   "\n");
    cat("Average number of responses:                 ", formatC(colMeans(eff), digits=1, format="f"), sep="\t",   "\n");
    cat("Average number of patients treated:          ", formatC(colMeans(pts), digits=1, format="f"), sep="\t",   "\n");
    cat("Percentage of trials stopped for safety:     ", nstops/ntrial*100, "\n");
    cat("Percentage of trials stopped for futility:   ", nstopf/ntrial*100, "\n\n\n");
    cat("Design specifications: \n\n");
    #cat("Size of safety run-in phase:\n");
    #cat(ncomb+2, sep="\t",  "\n");
    #cat("Size of adaptive randomization phase:\n");
    #cat(floor(n.ar), sep="\t",  "\n");
    cat("Prior on model parameters:\n");
    cat("N(0,1.34) \n");
    cat("Toxicity skeleton:\n");
    cat(round(p.skel,2), sep="\t",  "\n");
    cat("Activity skeletons:\n");
    round(q.skel,2);#cat(round(q.skel,3), sep="\t",  "\n");
  }
  
  impbpocrm<-function(y,ny,z,nz,theta,phi,dose.curr){
    library(dfcrm)
    library(nnet)
    
    d=ncomb=length(y)
    ptox.hat = numeric(ncomb); # estimate of toxicity prob
    peff.hat = numeric(ncomb); # estimate of efficacy prob
    
    delta=0.05
    #begin=min(which(ny>0))
    p.skel=getprior(halfwidth=0.05, target=theta, nu=ifelse(d %% 2 == 0,d/2,(d+1)/2), nlevel=d, model="empiric")
    toxskel=p.skel
        #getprior(halfwidth=0.04, target=theta, begin, nlevel=d, model="empiric")
      #getprior(halfwidth=delta, target=theta, nu=ifelse(d %% 2 == 0,d/2,(d+1)/2), nlevel=d, model="empiric")
    q.skel<-matrix(nrow=2*d-1,ncol=d)
    q.skel[1,]<-getprior(0.04,0.70,d,d)
    reverse=sort(q.skel[1,],decreasing=TRUE)
    for(i in 2:d){
      q.skel[i,]<-c(tail(q.skel[1,],d-i+1),reverse[2:i])
    }
    for(k in 2:d){
      q.skel[k+d-1,]<-c(tail(q.skel[1,],d-k+1),rep(0.7,k-1))
    }
    effskel=q.skel
    

    
    # if a single ordering is inputed as a vector, convert it to a matrix
    if(is.vector(q.skel)) q.skel=t(as.matrix(q.skel));
    
    nord.eff = nrow(q.skel);
    mprior.eff = rep(1/nord.eff, nord.eff); # prior for each efficacy ordering
    
    
    avar=1.34
    bcrmh<-function(a,p,y,n){
      s2=avar
      lik=exp(-0.5*a*a/s2)
      for(j in 1:length(p)){
        pj=p[j]**exp(a)
        lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
      }
      return(lik);
    }
    
    bcrmht<-function(a,p,y,n){
      s2=avar
      lik=a*exp(-0.5*a*a/s2)
      for(j in 1:length(p)){
        pj=p[j]**exp(a)
        lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
      }
      return(lik);
    }
    
    bcrmht2<-function(a,p,y,n){
      s2=avar
      lik=a^2*exp(-0.5*a*a/s2)
      for(j in 1:length(p)){
        pj=p[j]**exp(a)
        lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
      }
      return(lik);
    }
    
    
    marginal= integrate(bcrmh,lower=-Inf,upper=Inf, p=p.skel, y=y,n=ny,abs.tol = 0)$value;
    est=integrate(bcrmht,lower=-10,upper=10, p.skel, y, ny,abs.tol = 0)$value/marginal
    e2=integrate(bcrmht2,lower=-10,upper=10, p.skel, y, ny,abs.tol = 0)$value/marginal
    
    ptox.hat=p.skel**exp(est)
    post.var=e2-(est)^2
    crit=qnorm(0.5+0.9/2)
    ub=est+crit*sqrt(post.var)
    ptoxL=p.skel**exp(ub)
    
   # if(ptoxL[1]>theta){
   #   stops=1
   #   break
   # }
    
    
    marginal.eff = est.eff=e2.eff=rep(0, nord.eff);
    for(k in 1:nord.eff)
    {
      marginal.eff[k] = integrate(bcrmh,lower=-Inf,upper=Inf, p=q.skel[k,], y=z,n=nz,abs.tol = 0)$value;
      est.eff[k]=integrate(bcrmht,lower=-10,upper=10, q.skel[k,], z, nz,abs.tol = 0)$value/marginal.eff[k]
      e2.eff[k]=integrate(bcrmht2,lower=-10,upper=10, q.skel[k,], z, nz,abs.tol = 0)$value/marginal.eff[k]
    }		
    postprob.eff = (marginal.eff*mprior.eff)/sum(marginal.eff*mprior.eff);
    # toxicity model selection, identify the model with the highest posterior prob
    if(nord.eff>1){ 
      meff.sel = which.is.max(postprob.eff); 
    } else{
      meff.sel = 1;
    }
    peff.hat=q.skel[meff.sel,]**exp(est.eff[meff.sel])
    post.var.eff=e2.eff[meff.sel]-(est.eff[meff.sel])^2
    crit.eff=qnorm(0.5+0.9/2)
    lb.eff=est.eff[meff.sel]-crit.eff*sqrt(post.var.eff)
    peffU=q.skel[meff.sel,]**exp(lb.eff)
    
    
    
    ##dose-finding algorithm
    try=length(ny[ny>0])	 #highest dose tried
    mtd=which.is.max(-abs(ptox.hat-theta))	
    utility=peff.hat*as.numeric(ptox.hat<=ptox.hat[mtd])
    best=which.max(utility)#ifelse(sum(ny)<(2*maxn/3),sample(1:ncomb,1,prob=utility/sum(utility)),which.max(utility))
   
    #if(peffU[best]<phi){
    #  stopf=1
    #  break
    #}
    if(dose.curr==ncomb){
      rec=ifelse(best==ncomb,min(mtd,dose.curr),min(mtd,best))
    } else{
      if((best==try)||(best>dose.curr)){rec=min(mtd,dose.curr+1)}
      else if (best<dose.curr) {rec=min(mtd,dose.curr-1)}
      else {rec=min(mtd,dose.curr)}
    }
    ifelse(ptoxL[1]>theta || peffU[which.max(utility)]<phi,dose.curr<-"STOP STUDY",dose.curr<-rec)
    cat("Results of model-based estimation: \n\n");
    cat("Date and time:                           ", format(Sys.time()), sep="\t",  "\n"); 
    cat("Number of DLTs:        	                ", y, sep="\t",   "\n");
    cat("Number of patients evaluated for DLT:    ", ny, sep="\t",   "\n");
    cat("Number of activity responses:            ", z, sep="\t",   "\n");
    cat("Number of patients evaluated for Eff:    ", nz, sep="\t",   "\n");
    cat("Estimated DLT probabilities:             ", round(ptox.hat,2), sep="\t",  "\n");
    cat("Estimated activity probabilities:        ", round(peff.hat,2), sep="\t",  "\n");
    cat("Maximum tolerated dose:                  ",  mtd, sep="\t",  "\n");
    cat("Minimum effective and safe dose:             ",  best, sep="\t",  "\n");
    cat("Recommended dose level:                   ", dose.curr, sep="\t",  "\n\n\n");
    cat("Design specifications: \n\n");
    #cat("Size of safety run-in phase:\n");
    #cat(ncomb+2, sep="\t",  "\n");
    #cat("Size of adaptive randomization phase:\n");
    #cat(floor(nr), sep="\t",  "\n");
    cat("Prior on model parameters:\n");
    cat("N(0,1.34) \n");
    cat("Toxicity skeleton:\n");
    cat(round(toxskel,2), sep="\t",  "\n");
    cat("Activity skeletons:\n");
    round(effskel,2);#cat(round(q.skel,3), sep="\t",  "\n");
  }
  
  
  #reactive function to use inputs and produce output of functions
  output$simulation <- renderPrint(
    {
      validate(
        need(expr=input$ssize/input$cohortsize == round(input$ssize/input$cohortsize), message='Total number of patients is not multiple of cohort size!')
      )
      simbpocrm(as.numeric(unlist(strsplit(input$p0,",",TRUE))),as.numeric(unlist(strsplit(input$q0,",",TRUE))),input$tul,input$ell,input$cohortsize,input$ssize,input$ntrial,input$start.comb,input$n.stop,input$rseed)
    }                                                                                                  #p0,q0,tul,ell,psi,cohortsize,ncohort,ntrial,n.stop,start.comb,rseed
  )
  
  output$implementation <- renderPrint(
    {
      impbpocrm(as.numeric(unlist(strsplit(input$y,",",TRUE))),as.numeric(unlist(strsplit(input$ny,",",TRUE))),as.numeric(unlist(strsplit(input$z,",",TRUE))),as.numeric(unlist(strsplit(input$nz,",",TRUE))),input$theta,input$phi,input$dose.curr)
    }
  )
}

shinyApp(ui = ui, server = server)
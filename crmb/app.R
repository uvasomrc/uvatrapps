
library(shiny)


ui<-(navbarPage("Bayesian Continual Reassessment Method for Phase I Clinical Trials",
                #title
                tabPanel("Simulation",
                         h3("Web Application for simulating operating characteristics of the Bayesian CRM"),
                         #h4("Nolan A. Wages and Gina R. Petroni"),
                         h4("Division of Translational Research & Applied Statistics, University of Virginia;", a("nwages@virginia.edu", href="mailto:nwages@virginia.edu")),
                         br(),
                         h5("1. Enter an assumed set of true DLT probabilities, separated by commas.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         #input: boxes for truth input values OR text field
                         textInput("truth","True DLT probability at each dose level",value="0.04,0.11,0.25,0.40,0.55",width=NULL),
                         
                         br(),
                         h5("2. Enter the target DLT probability that defines the MTD for the study."),
                         #input: tul
                         numericInput("target","Target DLT rate",0.25,0,1,0.01,NULL),
                         
                         br(),
                         h5("3. Enter the cohort size required before the next model-based update. Cohort size may be 1, 2, or 3 patients."),
                         #input: cohortsize
                         numericInput("cohortsize","Cohort size",1,1,3,1,NULL),  
                         
                         br(),
                         h5("4. Enter the maximum sample size for the study. This number should be a multiple of the cohort size entered above."),
                         #input: ncohort
                         numericInput("ssize","Maximum number of patients",24,0,NA,1,NULL),
                         
                         br(),
                         h5("5. Enter the total number of patients treated on any dose required to stop the trial. At any point in the trial, if the recommendation is to assign
                            the next cohort to a dose that already has the entered number of patients treated on the dose, the study
                            is stopped and the recommended dose is declared the MTD. If the entered number is larger than the maximum sample size, each trial will accrue to
                            the maximum sample size."),
                         #input: ncohort
                         numericInput("n.stop","Number of patients needed on one dose to stop",25,0,NA,1,NULL),
                         
                         br(),
                         h5("6. Enter the number of simulations. A minimum of 1000 is recommended."),
                         #input: ntrial
                         numericInput("ntrial","Number of simulated trials", 10, 0, 100000, NA, NULL),
                         
                         br(),
                         p("7. Enter the index of the starting dose level.", strong("Note:"), "Index of", span("lowest", style = "color:blue"), "dose level is always 1. If the design allows for", em("'minus'"), "
                           dose levels (i.e. -2, -1, etc.), then the index of the starting dose should account for these lower levels (i.e. if -1 dose level allowed, starting dose is 2.)"),
                         #Even  if lowest dose level is", em("labeled"), "-1, index of lowest dose level is 1. 
                         #input: start
                         numericInput("x1","Index of starting dose level",1,1,NA,1,NULL),
                         
                         br(),
                         h5("8. Set the seed of the random number generator."),
                         #input: rseed
                         numericInput("rseed","Random seed",580,0,NA,1,NULL),
                         
                         br(),
                         h5("9. Specify the confidence level for safety stopping rule at the lowest study dose level."),
                         #input: tul
                         numericInput("cl","Confidence level for safety stopping",0.90,0,1,0.01,NULL),
                         
                         #input: submit button
                         submitButton("Run simulation study",icon("flask"),NULL),
                         
                         #output: simbcrm(truth,target,cohortsize,ncohort,ntrial,start,n.stop,rseed)
                         verbatimTextOutput("simulation"),
                         
                         mainPanel(
                           p(strong("This application simulates operating characteristics for the Bayesian continual reassessment method [1] with the following specifications.")),
                           p("1.", strong("Skipping Restriction:"), "The trial is not allowed to skip dose levels when escalating."),
                           p("2.", strong("Skeleton:"), "For the specified target DLT rate and total number of dose levels, the skeleton of power model d^exp(a) is generated according to Lee and Cheung (2009) [2] using a prior MTD located at the median dose level and a spacing measure of delta=0.05."),
                           p("3.", strong("Prior:"), "The prior distribution on the parameter", strong("a"), "is a mean zero normal distribution with the least informative prior variance [3]."),
                           #p("4.", strong("Accuracy Index:"), "The Accuracy Index is given in Equation 6.1 in Cheung [3]."),
                           p("4.", strong("Safety Stopping Rule:"), "Stop the trial for safety if the lower limit of an Agresti-Coull binomial confidence interval [4] for the lowest study dose level exceeds the target DLT rate"),
                           # p("Click",a("here",target="_blank",href="http://faculty.virginia.edu/model-based_dose-finding/nonparametric%20benchmark.R"),"to view the R Code for the benchmark."),
                           #p("Click",a("here", target="_blank", href="https://github.com/graham-wheeler/AplusB"),"to download the AplusB application to your computer from GitHub (set-up instructions provided)."),
                           strong("References:"),
                           p("[1] O'Quigley J, Pepe M, Fisher L (1990). Continual reassessment method: a practical design for phase I clinical trials in cancer, ", em("Biometrics;"), strong("46"),"(1): 33-48."),
                           p("[2] Lee and Cheung (2009). Model calibration in the continual reassessment method, ", em("Clinical Trials;"), strong("6"),"(3): 227-238."),
                           p("[3] Lee and Cheung (2011). Calibration of prior variance in the bayesian continual reassessment method, ", em("Statistics in Medicine;"), strong("30"),"(17): 2081-2089."),
                           p("[3] Agresti A, Coull BA (1998). Approximate is better than 'exact' for interval estimation of binomial proportions, ", em("American Statistician;"), strong("52"),": 119-126.")
                           #p("[3] Cheung (2013). Sample size formulae for the continual reassessment method, ", em("Clinical Trials;"), strong("10"),"(6): 852-861"),
                           #p("[3] Cheung YK (2011). ", em("Dose-finding by the continual reassessment method."), "Chapman and Hall/CRC press: New York.")
                           #p("[4] O'Quigley J, Shen LZ (1996). Continual reassessment method: a likelihood appraoch, ", em("Biometrics;"), strong("52"),"(2): 673-684")
                           
                         ),
                         
                         
                         tags$style(type="text/css",
                                    ".shiny-output-error { visibility: visible; }",
                                    ".shiny-output-error:before { visibility: visible; }"
                         )
                         ),
                #title
                tabPanel("Implementation",
                         h3("Web Application for implementation of the Bayesian CRM"),
                         # h4("Nolan A. Wages and Gina R. Petroni"),
                         h4("Division of Translational Research & Applied Statistics, University of Virginia;", a("nwages@virginia.edu", href="mailto:nwages@virginia.edu")),
                         #h5("1. Enter an assumed set of true DLT probabilities, separated by commas. The length of this set should be equal to the number of dose levels."),
                         #input: start
                         br(),
                         h4(strong("Design / Protocol Information")),
                         #br(),
                         p("1.Enter the index of the starting dose level.", strong("Note:"), "Index of", span("lowest", style = "color:blue"), "dose level is always 1. If the design allows for", em("'minus'"), "
                           dose levels (i.e. -2, -1, etc.), then the index of the starting dose should account for these lower levels (i.e. if -1 dose level allowed, starting dose is 2.)"),
                         #input: start
                         numericInput("begin","Index of starting dose level",1,1,NA,1,NULL),
                         
                         #br(),
                         h5("1. Enter the target DLT rate probability that defines the MTD for the study."),
                         #input: theta
                         numericInput("theta","Target DLT rate",0.25,0,1,0.01,NULL),
                         
                         br(),
                         h4(strong("Observed Trial Data (do not count 'replaced' patients)")),
                         #br(),
                         #input: boxes for truth input values OR text field
                         h5("2. Enter number of observed DLTs at each dose level. If none have been observed or a dose level has not yet been tried, enter '0'.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         textInput("y","Number of observed DLTs at each dose level",value="0,0,0,0,0",width=NULL),
                         
                         br(),
                         h5("3. Enter the number of patients evaluated for DLT at each dose level. If a dose level has not yet been tried, enter '0'.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         #input: boxes for truth input values OR text field
                         textInput("n","Number of patients evaluated for DLT at each dose level",value="1,0,0,0,0",width=NULL),
                         
                         br(),
                         h5("4. Enter the most recent dose level administered in the study."),
                         #input: dose.curr
                         numericInput("dose.curr","Current dose level",1,1,NA,1,NULL),
                         
                         
                         br(),
                         h5("5. Specify the confidence level for safety stopping rule at the lowest study dose level.."),
                         #input: tul
                         numericInput("cs","Confidence level used for safety stopping",0.90,0,1,0.01,NULL),
                         
                         
                         #input: submit button
                         submitButton("Get updated recommended dose level",icon("flask"),NULL),
                         
                         #output: impcrm(y,n,target)
                         verbatimTextOutput("implementation"),
                         
                         
                         mainPanel(
                           p(strong("This application computes a recommended dose level for the next patient in a phase I trial according to the Bayesian continual reassessment method [1] with the following specifications.")),
                           p("1.", strong("Skipping Restriction:"), "The trial is not allowed to skip dose levels when escalating."),
                           p("2.", strong("Skeleton:"), "For the specified target DLT rate and total number of dose levels, the  skeleton of power model d^exp(a) is generated according to Lee and Cheung (2009) [2] using a prior MTD located at the median dose level and a spacing measure of delta=0.05."),
                           p("3.", strong("Prior:"), "The prior distribution on the parameter", strong("a"), "is a mean zero normal distribution with the least informative prior variance [3]."),
                           p("4.", strong("Safety Stopping Rule:"), "Stop the trial for safety if the lower limit of an Agresti-Coull binomial confidence interval [4] for the lowest study dose level exceeds the target DLT rate"),
                           strong("References:"),
                           p("[1] O'Quigley J, Pepe M, Fisher L (1990). Continual reassessment method: a practical design for phase I clinical trials in cancer, ", em("Biometrics;"), strong("46"),"(1): 33-48."),
                           p("[2] Lee and Cheung (2009). Model calibration in the continual reassessment method, ", em("Clinical Trials;"), strong("6"),"(3): 227-238."),
                           p("[3] Lee and Cheung (2011). Calibration of prior variance in the bayesian continual reassessment method, ", em("Statistics in Medicine;"), strong("30"),"(17): 2081-2089."),
                           p("[4] Agresti A, Coull BA (1998). Approximate is better than 'exact' for interval estimation of binomial proportions, ", em("American Statistician;"), strong("52"),": 119-126.")
                           #p("[3] Cheung YK (2011). ", em("Dose-finding by the continual reassessment method."), "Chapman and Hall/CRC press: New York.")
                           #p("[4] O'Quigley J, Shen LZ (1996). Continual reassessment method: a likelihood appraoch, ", em("Biometrics;"), strong("52"),"(2): 673-684")
                           
                         ),
                         
                         
                         tags$style(type="text/css",
                                    ".shiny-output-error { visibility: visible; }",
                                    ".shiny-output-error:before { visibility: visible; }"
                         )
                         
                         
                         ),
                tabPanel("Safety stopping bounds",
                         h3("Web Application for generating safety stopping bounds at the lowest study dose level based on Agresti-Coull binomial confidence interval estimation."),
                         # h4("Nolan A. Wages and Gina R. Petroni"),
                         h4("Division of Translational Research & Applied Statistics, University of Virginia;", a("nwages@virginia.edu", href="mailto:nwages@virginia.edu")),
                         #h5("1. Enter an assumed set of true DLT probabilities, separated by commas. The length of this set should be equal to the number of dose levels."),
                         #input: start
                         br(),
                         
                         
                         
                         # br(),
                         h5("1. Specify the confidence level for safety stopping rule at the lowest study dose level."),
                         #input: tul
                         numericInput("clevel","Confidence level used for safety stopping",0.90,0,1,0.01,NULL),
                         
                         #br(),
                         h5("2. Enter the target DLT rate probability that defines the MTD for the study."),
                         #input: theta
                         numericInput("tul","Target DLT rate",0.25,0,1,0.01,NULL),
                         
                         
                         
                         #br(),
                         h5("3. Enter the maximum sample size for the study."),
                         #input: ncohort
                         numericInput("trialsize","Maximum number of patients",24,0,NA,1,NULL),
                         
                         
                         #input: submit button
                         submitButton("Get safety stopping bounds",icon("flask"),NULL),
                         
                         #output: impcrm(y,n,target)
                         verbatimTextOutput("stopping"),
                         
                         
                         mainPanel(
                           p(strong("This application computes the safety stopping bounds for the lowest study dose level based on Agresti-Coull binomial confidence interval estimation [1].")),
                           
                           strong("References:"),
                           p("[1] Agresti A, Coull BA (1998). Approximate is better than 'exact' for interval estimation of binomial proportions, ", em("American Statistician;"), strong("52"),": 119-126.")
                           #p("[2] Lee and Cheung (2009). Model calibration in the continual reassessment method, ", em("Clinical Trials;"), strong("6"),"(3): 227-238."),
                           #p("[3] Lee and Cheung (2011). Calibration of prior variance in the bayesian continual reassessment method, ", em("Statistics in Medicine;"), strong("30"),"(17): 2081-2089.")
                           #p("[3] Cheung YK (2011). ", em("Dose-finding by the continual reassessment method."), "Chapman and Hall/CRC press: New York.")
                           #p("[4] O'Quigley J, Shen LZ (1996). Continual reassessment method: a likelihood appraoch, ", em("Biometrics;"), strong("52"),"(2): 673-684")
                           
                         ),
                         tags$style(type="text/css",
                                    ".shiny-output-error { visibility: visible; }",
                                    ".shiny-output-error:before { visibility: visible;}"
                         )
                         
                         
                )
                )
     
)



server <- function(input, output) {
  simbcrm<-function(truth,target,cohortsize,ssize,n.stop,ntrial,x1,rseed,cl){
    
    library(dfcrm)
    library(binom)
    set.seed(rseed)    
    
    
    d=ndose=length(truth)
    skeleton=getprior(halfwidth=0.05, target=target, nu=x1, nlevel=d, model="empiric")
    #skeleton=round(getprior(halfwidth=0.05, target=target, nu=ifelse(d %% 2 == 0,d/2,(d+1)/2), nlevel=d, model="empiric"),2)
    
    
    foo <- crmsens(skeleton, target, model="empiric", detail=TRUE)
    sigma2<-seq(from = 0.1, to = 2, by = 0.001)
    kld<-rep(0,length(sigma2))
    for(j in 1:length(sigma2)){
      avar=sigma2[j]
      muj<-rep(0,length(skeleton ))
      for(i in 1:length(muj)){
        if(i==1){
          muj[i]<-pnorm(foo$Hset[i,2]/avar)-pnorm(-10/avar)
        } else {
          muj[i]<-pnorm(foo$Hset[i,2]/avar)-pnorm(foo$Hset[(i-1),2]/avar)
        }
      }
      kld[j]<-sum(muj*log(length(skeleton )*muj))
    }
    #plot(sigma2,kld)	
    
    ###'sli' is the least infomative prior standard deviation given 'skeleton'
    sli<-sigma2[which.min(kld)]
    
    
    ncohort=ssize/cohortsize
    
    
    
    ###Load the function 'bayescrm' 
    bayescrm<-function(truth,skeleton,target,cohortsize,ncohort,x1,n.stop,cl){
      
      bcrmh<-function(a,p,y,n){
        s2=sli^2
        lik=exp(-0.5*a*a/s2)
        for(j in 1:length(p)){
          pj=p[j]**exp(a)
          lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
        }
        return(lik);
      }
      
      bcrmht<-function(a,p,y,n){
        s2=sli^2
        lik=a*exp(-0.5*a*a/s2)
        for(j in 1:length(p)){
          pj=p[j]**exp(a)
          lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
        }
        return(lik);
      }
      
      bcrmht2<-function(a,p,y,n){
        s2=sli^2
        lik=a^2*exp(-0.5*a*a/s2)
        for(j in 1:length(p)){
          pj=p[j]**exp(a)
          lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
        }
        return(lik);
      }
      
      
      
      
      ### run a trial 	
      ndose = length(skeleton);   #number of combos
      y=n=rep(0,ndose);  #number of toxicities and number of treated patients at each dose level
      curr = x1;  # current dose level	 
      ptox.hat = numeric(ndose); # estimate of toxicity prob
      dose.select=rep(0,ndose); # a vector of indicators for dose selection
      stop=0; #indicate if trial stops early
      i=1	
      while(i <= ncohort)
      {
        
        
        
        # generate data for a new cohort of patients
        lasttox=rbinom(1,cohortsize,truth[curr])
        y[curr] = y[curr] + lasttox;
        n[curr] = n[curr] + cohortsize;
        
        
        
        
        marginal= integrate(bcrmh,lower=-Inf,upper=Inf, p=skeleton, y=y,n=n,abs.tol = 0)$value;
        est=integrate(bcrmht,lower=-10,upper=10, skeleton, y, n,abs.tol = 0)$value/marginal
        e2=integrate(bcrmht2,lower=-10,upper=10, skeleton, y, n,abs.tol = 0)$value/marginal
        ptox.hat=skeleton**exp(est)
        safety=ifelse(n[1]>1,binom.confint(y[1],n[1],conf.level=cl,methods="agresti-coull")$lower,0)
        
        if(safety>target){
          stop=1
          break
        }
        
        #distance=abs(ptox.hat-target)
        if (all(ptox.hat <= target)) {
          best <- length(skeleton)
        }
        else if (all(ptox.hat >= target)) {
          best <- 1
        }
        else {
          best <- order(abs(ptox.hat - target))[1]
        }
        
        # best=order(abs(ptox.hat - target))[1] #which.is.max(-distance)
        if(sum(lasttox)/cohortsize>=target){
          curr<-min(best, curr)
        }
        else {
          curr<-min(best,curr+1)
        }
        if(n[curr]>=n.stop){
          stop<-0
          break
        }
        i=i+1
      }
      if(stop==0){
        dose.select[curr]=dose.select[curr]+1;
      }
      return(list(dose.select=dose.select,tox.data=y,pt.allocation=n,stop=stop))
    }
    ##########'bpocrm' end here
    
    d.select<-tox<-pts<-matrix(nrow=ntrial,ncol=ndose)
    nstop=0
    
    for(i in 1:ntrial){
      result<-bayescrm(truth,skeleton,target,cohortsize,ncohort,x1,n.stop,cl)
      d.select[i,]=result$dose.select
      tox[i,]=result$tox.data
      pts[i,]=result$pt.allocation
      nstop=nstop+result$stop
    }
    #cat("Prior on model parameter:               N( 0,",round(sli**2,2),")","\n");
    # cat("Skeleton of working model:        ", round(skeleton,2), sep="\t",  "\n");
    cat("True DLT probability:             ", round(truth,3), sep="\t",  "\n");
    cat("MTD selection percentage:         ", formatC(colMeans(d.select)*100, digits=1, format="f"), sep="\t",  "\n");
    cat("Average number of DLTs:           ", formatC(colMeans(tox), digits=1, format="f"), sep="\t",   "\n");
    cat("Average number of patients:       ", formatC(colMeans(pts), digits=2, format="f"), sep="\t",   "\n");
    #cat("Accuracy index:                   ", round(1-length(truth)*(sum(abs(truth-target)*(colMeans(d.select)))/sum(abs(truth-target))),4), sep="\t",  "\n");
    cat("Percentage stopped for safety:    ", round(nstop/ntrial*100,3),sep="\t", "\n\n\n");
    cat("Design specifications: \n");
    cat("Prior on model parameter:\n");
    cat("N( 0,",round(sli**2,2),")","\n");
    cat("Skeleton of working model:\n");
    cat(round(skeleton,2), sep="\t",  "\n");
  }
  
  
  impbcrm<-function(y,n,theta,begin,dose.curr,cs){
    library(dfcrm)
    library(binom)
    library(nnet)
    d=ndose=length(y)
    ptox.hat = numeric(ndose);
    # skeleton=round(getprior(halfwidth=0.05, target=theta, nu=ifelse(d %% 2 == 0,d/2,(d+1)/2), nlevel=d, model="empiric"),2)
    skeleton=getprior(halfwidth=0.05, target=theta, nu=begin, nlevel=d, model="empiric")
    
    foo <- crmsens(skeleton, theta, model="empiric", detail=TRUE)
    sigma2<-seq(from = 0.1, to = 2, by = 0.001)
    kld<-rep(0,length(sigma2))
    for(j in 1:length(sigma2)){
      avar=sigma2[j]
      muj<-rep(0,length(skeleton ))
      for(i in 1:length(muj)){
        if(i==1){
          muj[i]<-pnorm(foo$Hset[i,2]/avar)-pnorm(-10/avar)
        } else {
          muj[i]<-pnorm(foo$Hset[i,2]/avar)-pnorm(foo$Hset[(i-1),2]/avar)
        }
      }
      kld[j]<-sum(muj*log(length(skeleton )*muj))
    }
    #plot(sigma2,kld)	
    
    ###'sli' is the least infomative prior standard deviation given 'skeleton'
    sli<-sigma2[which.min(kld)]
    
    
    bcrmh<-function(a,p,y,n){
      s2=sli^2
      lik=exp(-0.5*a*a/s2)
      for(j in 1:length(p)){
        pj=p[j]**exp(a)
        lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
      }
      return(lik);
    }
    
    bcrmht<-function(a,p,y,n){
      s2=sli^2
      lik=a*exp(-0.5*a*a/s2)
      for(j in 1:length(p)){
        pj=p[j]**exp(a)
        lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
      }
      return(lik);
    }
    
    bcrmht2<-function(a,p,y,n){
      s2=sli^2
      lik=a^2*exp(-0.5*a*a/s2)
      for(j in 1:length(p)){
        pj=p[j]**exp(a)
        lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
      }
      return(lik);
    }
    
    marginal= integrate(bcrmh,lower=-Inf,upper=Inf, p=skeleton, y=y,n=n,abs.tol = 0)$value;
    est=integrate(bcrmht,lower=-10,upper=10, skeleton, y, n,abs.tol = 0)$value/marginal
    e2=integrate(bcrmht2,lower=-10,upper=10, skeleton, y, n,abs.tol = 0)$value/marginal
    ptox.hat=skeleton**exp(est)
    
    safety=ifelse(n[1]>1,binom.confint(y[1],n[1],conf.level=cs,methods="agresti-coull")$lower,0)
    
    
    distance=abs(ptox.hat-theta)
    best=which.is.max(-distance)
    
    
    ifelse(safety>theta,dose.curr<-"STOP STUDY FOR SAFETY",dose.curr<-min(best,dose.curr+1))
    
    cat("Date and time:                           ", format(Sys.time()), sep="\t",  "\n");
    #cat("Prior on model parameter:                       N( 0,",round(sli**2,2),")","\n");
    # cat("Skeleton of working model:               ", round(skeleton,2), sep="\t",  "\n");
    cat("Number of DLTs:        	                ", y, sep="\t",   "\n");
    cat("Number of patients evaluated for DLT:    ", n, sep="\t",   "\n");
    cat("Estimated DLT probabilities:             ", round(ptox.hat,2), sep="\t",  "\n");
    cat("Target DLT rate:                         ", theta, sep="\t",  "\n");
    cat("Recommended dose level:                  ", dose.curr, sep="\t",  "\n\n");
    cat("Design specifications: \n");
    cat("Prior on model parameter:\n");
    cat("N( 0,",round(sli**2,2),")","\n");
    cat("Skeleton of working model:\n");
    cat(round(skeleton,2), sep="\t",  "\n");
  }
  
  toxmonitoring<-function(clevel,tul,trialsize){
    library(binom)
    mintox<-rep(0,trialsize)
    for(i in 1:trialsize){
      lv=rep(0,i)
      for(k in 1:i){
        lv[k]<-as.numeric(binom.confint(k,i,conf.level=clevel,methods="agresti-coull")$lower>tul)
      }
      mintox[i]<-ifelse(all(lv==0),0,min(which(lv==1)))
    }
    df=data.frame(mintox[2:trialsize],2:trialsize)
    names(df)=c("#DLTs","#pts")
    cat("Stop the study for safety if the observed DLT rate at lowest study dose level >= #DLTs out of #pts treated at lowest study dose level. \n\n");
    #cat(mat, sep="\t","\n")
    print(df,row.names=FALSE);
    #cat("Stop the study if number of DLTs at lowest study dose level >=:\n")
    # cat(mintox[2:trialsize], sep="\t","\n")
  }
  
  #reactive function to use inputs and produce output of functions
  output$simulation <- renderPrint(
    {
      validate(
        need(expr=input$ssize/input$cohortsize == round(input$ssize/input$cohortsize), message='Total number of patients is not multiple of cohort size!')
      )
      simbcrm(as.numeric(unlist(strsplit(input$truth,",",TRUE))),input$target,input$cohortsize,input$ssize,input$n.stop,input$ntrial,input$x1,input$rseed,input$cl)
    }                                                                                                  #p0,q0,tul,ell,psi,cohortsize,ncohort,ntrial,n.stop,start.comb,rseed
  )
  
  output$implementation <- renderPrint(
    {
      #validate(
      # need(expr=length(input$y) == length(input$n), message='Different number of dose levels specified in observed data')
      #)
      impbcrm(as.numeric(unlist(strsplit(input$y,",",TRUE))),as.numeric(unlist(strsplit(input$n,",",TRUE))),input$theta,input$begin,input$dose.curr,input$cs)
    }
  )
  output$stopping <- renderPrint(
    {
      #validate(
      # need(expr=length(input$y) == length(input$n), message='Different number of dose levels specified in observed data')
      #)
      toxmonitoring(input$clevel,input$tul,input$trialsize)
    }
  )
}

shinyApp(ui = ui, server = server)
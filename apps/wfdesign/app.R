
library(shiny)

ui<-(navbarPage("Adaptive dose-finding based on safety and feasibility in early-phase clinical trials of adoptive cell immunotherapy",
                #title
                tabPanel("Simulation", 
                         h3("Web Application for simulating operating characteristics of the Wages and Fadul (2019) design"),
                         #h4("Nolan A. Wages and Gina R. Petroni"),
                         h4("Division of Translational Research & Applied Statistics, University of Virginia;", a("nwages@virginia.edu", href="mailto:nwages@virginia.edu")),
                         br(),
                         h5("1. Enter an assumed set of true DLT probabilities, separated by commas.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         #input: boxes for truth input values OR text field
                         textInput("truth","True DLT probabilities",value="0.10,0.30,0.50,0.70,0.80",width=NULL),
                         
                         br(),
                         h5("2. Enter an assumed set of true feasibility probabilities, separated by commas.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         #input: tul
                         textInput("fprob","True feasibility probabilities",value="0.90,0.75,0.50,0.25,0.05",width=NULL),
                         
                         #br(),
                        # h5("3. Enter a set of prior means for feasibility probabilities, separated by commas.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         ##input: tul
                        # textInput("q.skel","Prior means for feasibility probabilities",value="0.90,0.85,0.80,0.75,0.70",width=NULL),
                         
                         br(),
                         h5("3. Enter the target DLT probability that defines the MTD for the study."),
                         #input: tul
                         numericInput("tul","Target DLT rate",0.30,0,1,0.01,NULL),
                         
                         br(),
                         h5("4. Enter the minimum acceptable probability of feasibility that defines the threshold for desireable feasibility"),
                         #input: ell
                         numericInput("ell","Minimum acceptable feasibility rate",0.50,0,1,0.01,NULL),
                         
                         br(),
                         h5("5. Specify the probability cutoff for defining the set of feasible dose levels."),
                         #input: tul
                         numericInput("puf","Upper probability cutoff for defining feasible doses",0.70,0,1,0.01,NULL),
                         
                         br(),
                         h5("6. Specify the probability cutoff for defining the safety stopping rule based on an unacceptable high DLT rate."),
                         #input: tul
                         numericInput("put","Upper probability cutoff for defining safety stopping rule",0.70,0,1,0.01,NULL),
                         
                         br(),
                         h5("7. Enter the maximum number of participants who will be infused at some dose level and evaluated for toxicity."),
                         #input: cohortsize
                         numericInput("Nmax","Maximum number of participants evaluated for DLT",24,0,NA,1,NULL),  
                         
                         br(),
                         h5("8. Enter the maximum target accrual for the study. This is the number of participants who will have cells extracted."),
                         #input: ncohort
                         numericInput("Mmax","Total maximum number of planned participants",30,0,NA,1,NULL),
                         
                         br(),
                         h5("9. Enter the number of simulations. A minimum of 1000 is recommended."),
                         #input: ntrial
                         numericInput("ntrial","Number of simulations", 1000, 0, 100000, NA, NULL),
                         
                         br(),
                         p("10. Enter the index of the starting dose level.", strong("Note:"), "Index of", span("lowest", style = "color:blue"), "dose level is always 1. If the design allows for", em("'minus'"), "
                           dose levels (i.e. -2, -1, etc.), then the index of the starting dose should account for these lower levels (i.e. if -1 dose level allowed, index of starting dose is 2.)"),
                         #input: start.comb
                         numericInput("start","Index of starting dose level",1,1,NA,1,NULL),
                         
                         br(),
                         h5("11. Set the seed of the random number generator."),
                         #input: rseed
                         numericInput("rseed","Random Seed",47579,0,NA,1,NULL),
                         
                         #input: submit button
                         submitButton("Run simulation study",icon("flask"),NULL),
                         
                         #output: simbpocrm(p0,q0,tul,ell,psi,cohortsize,ncohort,ntrial,start.comb,n.stop,rseed)
                         verbatimTextOutput("simulation"),
                         
                         mainPanel(
                           p(strong("This application simulates operating characteristics for the early-phase method of Wages and Fadul [1].")),
                           #p("1.", strong("Dose-finding algorithm:"), "The trial does not skip over an untried dose level."),
                           #p("2.", strong("Safety:"), "For the specified target DLT rate and total number of dose levels, the skeleton of power model d^exp(a) is generated according to Lee and Cheung (2009) [2] using a prior MTD located at the median dose level and a spacing measure of delta=0.05."),
                           #p("3.", strong("Prior:"), "The prior distribution on the parameter", strong("a"), "is the least informative, mean zero normal prior according to Lee and Cheung (2011) [3]."),
                           #p("4.", strong("Accuracy Index:"), "The Accuracy Index is given in Equation 6.1 in Cheung [3]."),
                          # p("2.", strong("Safety Stopping Rule:"), "Stop the trial for safety if the lower limit of a 90% probability interval exceeds the target DLT rate."),
                           #p("3.", strong("Futility Stopping Rule:"), "Stop the trial for futility if the upper limit of a 90% probability interval is lower than the minimum acceptable efficacy rate."),
                           strong("References:"),
                           p("[1] Wages NA, Fadul CE (2019). Adaptive dose finding based on safety and feasibility in early-phase clinical trials of adoptive cell immunotherapy.", em("Clin Trials;"), strong("in press"))
                          
                         ),
                         tags$style(type="text/css",
                                    ".shiny-output-error { visibility: visible; }",
                                    ".shiny-output-error:before { visibility: visible; }"
                         )
                ),
                #title
                tabPanel("Implementation",
                         h3("Web Application for implementation of the Wages and Fadul (2019) method"),
                         # h4("Nolan A. Wages and Gina R. Petroni"),
                         h4("Division of Translational Research & Applied Statistics, University of Virginia;", a("nwages@virginia.edu", href="mailto:nwages@virginia.edu")),
                         #h5("1. Enter an assumed set of true DLT probabilities, separated by commas. The length of this set should be equal to the number of dose levels."),
                         #input: start
                         br(),
                         h4(strong("Design / Protocol Information")),
                         #br(),
                         #p("1.Enter the index of the starting dose level.", strong("Note:"), "Index of", span("lowest", style = "color:blue"), "dose level is always 1. If the design allows for", em("'minus'"), "
                         #dose levels (i.e. -2, -1, etc.), then the index of the starting dose should account for these lower levels (i.e. if -1 dose level allowed, starting dose is 2.)"),
                         ##input: start
                         #numericInput("begin","Index of starting dose level",1,1,NA,1,NULL),
                         
                         #br(),
                         h5("1. Enter the target DLT rate probability that defines the MTD for the study."),
                         #input: theta
                         numericInput("theta","Target DLT rate",0.20,0,1,0.01,NULL),
                         
                         br(),
                         h4(strong("Observed Trial Data (do not count 'replaced' patients)")),
                         #br(),
                         #input: boxes for truth input values OR text field
                         h5("2. Enter number of observed DLTs at each dose level. If none have been observed or a dose level has not yet been tried, enter '0'.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         textInput("y","Number of observed DLTs at each dose level",value="0,1,0,0,0",width=NULL),
                         
                         br(),
                         h5("3. Enter the number of patients evaluated for DLT at each dose level. If a dose level has not yet been tried, enter '0'.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         #input: boxes for truth input values OR text field
                         textInput("n","Number of patients evaluated for DLT at each dose level",value="2,4,0,0,0",width=NULL),
                         
                         br(),
                         h5("4. Enter the most recent dose level administered in the study."),
                         #input: dose.curr
                         numericInput("dose.curr","Current dose level",1,1,NA,1,NULL),
                         
                         br(),
                         h5("5. Specify the threshold for defining the safety stopping rule based on an unacceptable high DLT rate at the lowest study dose level."),
                         #input: tul
                         numericInput("cs","Threshold for defining safety stopping rule",0.90,0,1,0.01,NULL),
                         
                         #input: submit button
                         submitButton("Get updated recommended dose level",icon("flask"),NULL),
                         
                         #output: impcrm(y,n,target)
                         verbatimTextOutput("implementation"),
                         
                         
                         mainPanel(
                           p(strong("This application computes an updated planned dose level for the next patient in a phase I trial according to the Wages and Fadul method [1].")),
                           #p("1.", strong("Skipping Restriction:"), "The trial is not allowed to skip dose levels when escalating."),
                           # p("2.", strong("Skeleton:"), "For the specified target DLT rate and total number of dose levels, the  skeleton of power model d^exp(a) is generated according to Lee and Cheung (2009) [2] using a prior MTD located at the median dose level and a spacing measure of delta=0.05."),
                           #  p("1.", strong("Prior:"), "The prior distribution is beta(a,b) chosen according to Conaway, Dunbar, and Peddada [1]."),
                           #p("4.", strong("Accuracy Index:"), "The Accuracy Index is given in Equation 6.1 in Cheung [3]."),
                           #p("2.", strong("Safety Stopping Rule:"), "Stop the trial for safety if the Pr(DLT rate at dose level 1 > target | data) > 0.95"),
                           #  p("3.", strong("Starting dose level:"), "Each trial should begin at the lowest study dose level."),
                           # p("Click",a("here",target="_blank",href="http://faculty.virginia.edu/model-based_dose-finding/nonparametric%20b
                           strong("References:"),
                           p("[1] Wages NA, Fadul CE (2019). Adaptive dose finding based on safety and feasibility in early-phase clinical trials of adoptive cell immunotherapy.", em("Clin Trials;"), strong("in press"))
                           #p("[2] Lee and Cheung (2009). Model calibration in the continual reassessment method, ", em("Clinical Trials;"), strong("6"),"(3): 227-238."),
                           #p("[3] Lee and Cheung (2011). Calibration of prior variance in the bayesian continual reassessment method, ", em("Statistics in Medicine;"), strong("30"),"(17): 2081-2089.")
                           #p("[3] Cheung YK (2011). ", em("Dose-finding by the continual reassessment method."), "Chapman and Hall/CRC press: New York.")
                           #p("[4] O'Quigley J, Shen LZ (1996). Continual reassessment method: a likelihood appraoch, ", em("Biometrics;"), strong("52"),"(2): 673-684")
                           
                         ),
                         
                         
                         tags$style(type="text/css",
                                    ".shiny-output-error { visibility: visible; }",
                                    ".shiny-output-error:before { visibility: visible; }"
                         )
                         
                         
                ),
                tabPanel("Safety stopping bounds",
                         h3("Web Application for generating safety stopping bounds at the lowest study dose level computed from the posterior distribution of the DLT rate at the lowest dose level based on a beta prior distribution."),
                         # h4("Nolan A. Wages and Gina R. Petroni"),
                         h4("Division of Translational Research & Applied Statistics, University of Virginia;", a("nwages@virginia.edu", href="mailto:nwages@virginia.edu")),
                         #h5("1. Enter an assumed set of true DLT probabilities, separated by commas. The length of this set should be equal to the number of dose levels."),
                         #input: start
                         br(),
                         
                         
                         
                         # br(),
                         h5("1. Specify the threshold for defining the safety stopping rule based on an unacceptable high DLT rate at the lowest study dose level."),
                         #input: tul
                         numericInput("clevel","Threshold for defining safety stopping rule",0.90,0,1,0.01,NULL),
                         
                         #br(),
                         h5("2. Enter the target DLT rate probability that defines the MTD for the study."),
                         #input: theta
                         numericInput("atul","Target DLT rate",0.30,0,1,0.01,NULL),
                         
                         
                         
                         #br(),
                         h5("3. Enter the maximum sample size for the study."),
                         #input: ncohort
                         numericInput("trialsize","Maximum number of patients",24,0,NA,1,NULL),
                         
                         
                         #input: submit button
                         submitButton("Get safety stopping bounds",icon("flask"),NULL),
                         
                         #output: impcrm(y,n,target)
                         verbatimTextOutput("stopping"),
                         
                         
                         mainPanel(
                           p(strong("This application computes the safety stopping bounds for the lowest study dose level computed from the posterior distribution of the DLT rate at the lowest dose level based on a beta prior distribution."))
                           
                           # strong("References:"),
                           # p("[1] Agresti A, Coull BA (1998). Approximate is better than 'exact' for interval estimation of binomial proportions, ", em("American Statistician;"), strong("52"),": 119-126.")
                           #p("[2] Lee and Cheung (2009). Model calibration in the continual reassessment method, ", em("Clinical Trials;"), strong("6"),"(3): 227-238."),
                           #p("[3] Lee and Cheung (2011). Calibration of prior variance in the bayesian continual reassessment method, ", em("Statistics in Medicine;"), strong("30"),"(17): 2081-2089.")
                           #p("[3] Cheung YK (2011). ", em("Dose-finding by the continual reassessment method."), "Chapman and Hall/CRC press: New York.")
                           #p("[4] O'Quigley J, Shen LZ (1996). Continual reassessment method: a likelihood appraoch, ", em("Biometrics;"), strong("52"),"(2): 673-684")
                           
                         ),
                         tags$style(type="text/css",
                                    ".shiny-output-error { visibility: visible; }",
                                    ".shiny-output-error:before { visibility: visible;}"
                         )
                         
                         
                ),
                tabPanel("Feasibility stopping bounds",
                         h3("Web Application for generating feasibility stopping bounds at the lowest study dose level computed from the posterior distribution of the feasibility rate at the lowest dose level based on a beta prior distribution."),
                         # h4("Nolan A. Wages and Gina R. Petroni"),
                         h4("Division of Translational Research & Applied Statistics, University of Virginia;", a("nwages@virginia.edu", href="mailto:nwages@virginia.edu")),
                         #h5("1. Enter an assumed set of true DLT probabilities, separated by commas. The length of this set should be equal to the number of dose levels."),
                         #input: start
                         br(),
                         
                         #qmeans,flevel,aell,msize
                         
                        # br(),
                        # h5("1. Enter a prior mean for the feasibility probability at the lowest dose level."),
                        # #input: tul
                        # numericInput("qmeans","Prior mean for feasibility probability at lowest dose",0.90,0,1,0.01,NULL),
                         
                         
                         
                         br(),
                         h5("1. Specify the probability cutoff for defining the set of feasible dose levels."),
                         #input: tul
                         numericInput("flevel","Upper probability cutoff for defining feasible doses",0.70,0,1,0.01,NULL),
                         
                         #br(),
                         h5("2. Enter the minimum acceptable probability of feasibility that defines the threshold for desireable feasibility"),
                         #input: ell
                         numericInput("aell","Minimum acceptable feasibility rate",0.50,0,1,0.01,NULL),
                         
                         
                         #br(),
                         h5("3. Enter the maximum target sample size for the study."),
                         #input: ncohort
                         numericInput("msize","Maximum number of patients",30,0,NA,1,NULL),
                         
                         
                         #input: submit button
                         submitButton("Get feasibility stopping bounds",icon("flask"),NULL),
                         
                         #output: impcrm(y,n,target)
                         verbatimTextOutput("fstopping"),
                         
                         
                         mainPanel(
                           p(strong("This application computes the feasibility stopping bounds computed from the posterior distribution of the feasbility rate at the lowest dose level based on a beta prior distribution."))
                           
                           # strong("References:"),
                           # p("[1] Agresti A, Coull BA (1998). Approximate is better than 'exact' for interval estimation of binomial proportions, ", em("American Statistician;"), strong("52"),": 119-126.")
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
  WFdesign.sim<-function(truth,fprob,tul,ell,puf,put,Nmax,Mmax,ntrial,start,rseed){
    set.seed(rseed)
    library(Iso)
    library(nnet)
    ndose=length(truth)
    
    ###calculate prior for unfusibility
    q.skel<-seq(0.90,0.90-(0.05*(ndose-1)),length.out=ndose)	  ##prior means for infusibility probabilities
    aplus=rep(0,length(q.skel))	          ##prior sample size at each dose level
    for(i in 1:length(q.skel)){
      x=q.skel[i]/2
      mu=q.skel[i]
      f<-function(b){
        1-pbeta(x,mu*b/(1-mu),b)-0.95
      }
      bF=uniroot(f,c(0.0001,100))$root
      aF=mu*bF/(1-mu)
      aplus[i]=round(aF+bF,2)
    }

    ###calculate prior for toxicity
    x<-2*tul
    mu<-tul
    u<-0.95
    f<-function(b){
      pbeta(x,mu*b/(1-mu),b)-u
    }
    b0<-uniroot(f,c(0.0001,100))$root
    a0<-mu*b0/(1-mu)
    round(c(a0,b0),2)

  #  put=0.9		          ##upper probability cutoff for safety
   # puf=0.7             ##upper probability cutoff for feasibility 
    ###Load the function 'WFdesign' 
    WFdesign<-function(truth,fprob,a0,b0,q.skel,aplus,tul,ell,Nmax,Mmax,start,put,puf){
      
      ### run a trial 	
      #ndose = length(truth);   #number of combos
      y=n=z=M=dose.select=ptox.hat=peff.hat=p.infeasible=rep(0,ndose); 
      jstar = start;   #current dose level	 
      stopf=stops=untx=0; #indicate if trial stops early
      i=1	
      A=q.skel*aplus
      B=aplus-A
      
      while(i <= Mmax)
      {
        L=runif(1)
        fdoses=which(fprob>L)
        Y=ifelse(length(fdoses)==0,0,max(fdoses))
        if(Y==0){untx=untx+1}
        z[fdoses]=z[fdoses]+1
        M=M+1
        
        for(j in 1:ndose){
          p.infeasible[j] = pbeta(ell, z[j] + A[j], M[j] - z[j] + B[j]);
        }
        
        if(p.infeasible[1]> puf){
          stopf=1
          break
        }
        
        
        if(Y>0){
          curr=min(jstar,Y)
          y[curr] = y[curr] + rbinom(1,1,truth[curr]);
          n[curr] = n[curr] + 1;
          tried=which(n>0)
          
          if(1 - pbeta(tul, y[1] + a0, n[1] - y[1] + b0) > put){
            stops=1
            break
          }
          
          if(length(tried)<ndose & all(y==0)){
            jstar<-jstar+1
          } else {
            u=(y[tried]+a0)/(n[tried]+a0+b0)
            pipost=pava(u,w=n[tried])
            lossvec=ifelse(pipost > tul, (1-0.5)*(pipost-tul), 0.5*(tul-pipost))  
            T=lossvec==min(lossvec)
            poss=which(T)
            if(sum(T)==1){
              sugglev=poss
            } else {
              if(all(pipost[poss]>tul)){
                sugglev=min(poss)
              } else {
                sugglev=max(poss)
              }
            }
            if(pipost[sugglev]<tul & length(tried)<ndose){
              jstar=ifelse(n[sugglev+1]==0,sugglev+1,sugglev)					
            } else {
              jstar=sugglev
            }
          }
        }
        
        fset=which(p.infeasible<puf)
        jstar=min(max(fset),jstar)
        
        if(sum(n)>=Nmax){
          stops=0
          break
        }
        i=i+1
      }
      if(stops==0 & stopf==0){
        fmtd=min(max(fset),jstar)
        dose.select[fmtd]=dose.select[fmtd]+1;
      }
      return(list(dose.select=dose.select,tox.data=y,pt.allocation=n,stopf=stopf,stops=stops,unt=untx))
    }
    ##########'WFdesign' end here
    
    
    nuntx=rep(0,ntrial)
    dose.select<-y<-z<-n<-naf<-matrix(nrow=ntrial,ncol=ndose)
    nstopf=nstops=0
    
    for(i in 1:ntrial){
      result<-WFdesign(truth,fprob,a0,b0,q.skel,aplus,tul,ell,Nmax,Mmax,start,put,puf)
      dose.select[i,]=result$dose.select
      y[i,]=result$tox.data
      n[i,]=result$pt.allocation
      #z[i,]=result$feas.data
      #naf[i,]=result$feas.eval
      nuntx[i]=result$unt
      nstopf=nstopf+result$stopf
      nstops=nstops+result$stops
    }
    cat("True DLT probability:\n");
    cat(round(truth,3), sep="\t",  "\n");
    cat("True feasibility probability:\n");
    cat(round(fprob,3), sep="\t",  "\n");
    cat("Feasible maximum tolerated dose (FMTD) selection percentage:\n");
    cat(formatC(colMeans(dose.select)*100, digits=1, format="f"), sep="\t",  "\n");
    cat("Average nmber of DLTs:\n");
    cat(formatC(colMeans(y), digits=1, format="f"), sep="\t",   "\n");
    cat("Average number of patients infused:\n");
    cat(formatC(colMeans(n), digits=1, format="f"), sep="\t",   "\n");
    #cat("number feasible:   		    			", formatC(colMeans(z), digits=1, format="f"), sep="\t",   "\n");
    #cat("number evaluated feasibility:      			", formatC(colMeans(M), digits=1, format="f"), sep="\t",   "\n");
    cat("Average number of patients not treated:\n");
    cat(formatC(mean(nuntx), digits=1, format="f"), sep="\t",   "\n");
    cat("percentage of stop (safety):\n");
    cat(nstops/ntrial*100, "\n");
    cat("percentage of stop (feasibility):\n");
    cat(nstopf/ntrial*100, "\n\n");
    cat("Design specifications: \n");
    cat("Prior distribution on DLT rate at each dose level:\n");
    cat("Beta(",round(a0,2),",",round(b0,2),")","\n");
    cat("Prior distribution on feasibility rate is Beta(a,b) with:\n");
    cat("a=(",round(q.skel*aplus,2),") and b=(",round(aplus-(q.skel*aplus),2),")","\n");
    cat("Safety stopping rule:\n");
    cat("Stop the trial if the Pr(DLT rate at the lowest study dose level > target DLT rate | data) > ",put,"\n");
    cat("Feasibility stopping rule:\n");
    cat("Stop the trial if the Pr(feasibility rate at the lowest study dose level < minimum feasbility rate | data) > ",puf,"\n");
  } 
  ##########'WFdesign.sim' end here
  
  impbcdp<-function(y,n,theta,dose.curr,cs){
    library(Iso)
    d=ndose=length(y)
    ptox.hat = numeric(ndose);
    
    ###calculate prior
    x<-2*theta
    mu<-theta
    u<-0.95
    
    f<-function(b){
      pbeta(x,mu*b/(1-mu),b)-u
    }
    b0<-uniroot(f,c(0.0001,100))$root
    a0<-mu*b0/(1-mu)
    round(c(a0,b0),2)
    
    tried=which(n>0)
    
    u=(y[tried]+a0)/(n[tried]+a0+b0)
    pipost=pava(u,w=n)
    currpi=pipost[dose.curr]
    lossvec=abs(pipost-theta)  
    T=lossvec==min(lossvec)
    poss=which(T)
    if(sum(T)==1){
      sugglev=poss
    } else {
      if(all(pipost[poss]>theta)){
        sugglev=min(poss)
      } else {
        sugglev=max(poss)
      }
    }
    
    
    if(length(tried)<ndose){
      if(pipost[sugglev]<theta){
        dose.curr=ifelse(n[sugglev+1]==0,sugglev+1,sugglev)					
      } else {
        dose.curr=sugglev
      }
    } else {
      dose.curr=sugglev
    }
    
    safety=ifelse(n[1]>1,1 - pbeta(theta, y[1] + a0, n[1] - y[1] + b0),0)
    
    
    ifelse(safety>cs,dose.curr<-"STOP STUDY FOR SAFETY",dose.curr<-dose.curr)
    
    cat("Date and time:                                 ", format(Sys.time()), sep="\t",  "\n");
    cat("Observed number of DLTs:        	      ", y, sep="\t",   "\n");
    cat("Number of patients evaluated for DLT:          ", n, sep="\t",   "\n");
    cat("Estimated DLT probability at current dose:     ", round(currpi,4), sep="\t",  "\n");
    cat("Target DLT rate:                               ", theta, sep="\t",  "\n");
    cat("Recommended next dose level:                   ", dose.curr, sep="\t",  "\n\n");
    cat("Design specifications: \n");
    cat("Prior distribution on DLT rate at each dose level:\n");
    cat("beta(",round(a0,2),",",round(b0,2),")","\n");
    cat("Safety stopping rule:\n");
    cat("Stop the trial if the Pr(DLT rate at the lowest study dose level > target | data) > ",cs,"\n");
  }
  toxmonitoring<-function(clevel,atul,trialsize){
    ###calculate prior
    x<-2*atul
    mu<-atul
    u<-0.95
    
    f<-function(b){
      pbeta(x,mu*b/(1-mu),b)-u
    }
    b0<-uniroot(f,c(0.0001,100))$root
    a0<-mu*b0/(1-mu)
    
    mintox<-rep(0,trialsize)
    for(i in 1:trialsize){
      lv=rep(0,i)
      for(k in 1:i){
        lv[k]<-lv[k]<-as.numeric(1-pbeta(atul,a0+k,b0+i-k)>clevel)
      }
      mintox[i]<-ifelse(all(lv==0),"n/a",min(which(lv==1)))
    }
    df=data.frame(mintox[2:trialsize],2:trialsize)
    names(df)=c("#DLTs","#pts")
    cat("Design specifications: \n");
    cat("Prior distribution on DLT rate at each dose level:\n");
    cat("beta(",round(a0,2),",",round(b0,2),")","\n\n");
    cat("Stop the study for safety if the observed DLT rate at lowest study dose level >= #DLTs out of #pts treated at lowest study dose level. \n\n");
    #cat(mat, sep="\t","\n")
    print(df,row.names=FALSE);
    #cat("Stop the study if number of DLTs at lowest study dose level >=:\n")
    # cat(mintox[2:trialsize], sep="\t","\n")
  }
  feasmonitoring<-function(flevel,aell,msize){
    ###calculate prior
    aplus=0	      ##prior sample size at each dose level
    
    qmeans=0.90
        ###calculate prior for unfusibility
    x=qmeans/2
    mu=qmeans
    f<-function(b){
      1-pbeta(x,mu*b/(1-mu),b)-0.95
    }
    bF=uniroot(f,c(0.0001,100))$root
    aF=mu*bF/(1-mu)
    aplus=round(aF+bF,2)
    
    
    A=qmeans*aplus
    B=aplus-A
    
    
    minfeas<-rep(0,msize)
    for(i in 1:msize){
      lv=rep(0,i)
      for(k in 1:i){
        lv[k]<-as.numeric(pbeta(aell,A+k,B+i-k)>flevel)
      }
      minfeas[i]<-ifelse(all(lv==0),"n/a",i-max(which(lv==1)))
    }
    df=data.frame(minfeas[3:msize],3:msize)
    names(df)=c("#Not feasible","#pts")
    cat("Design specifications: \n");
    cat("Prior distribution on feasibiltiy rate at the lowest dose level:\n");
    cat("beta(",round(A,2),",",round(B,2),")","\n\n");
    cat("Stop the study for infeasibility if the observed infeasibility rate at lowest study dose level >= #Not feasible out of #pts treated at lowest study dose level. \n\n");
    #cat(mat, sep="\t","\n")
    print(df,row.names=FALSE);
    #cat("Stop the study if number of DLTs at lowest study dose level >=:\n")
    # cat(minfeas[2:msize], sep="\t","\n")
  }
  
  
  #reactive function to use inputs and produce output of functions
  output$simulation <- renderPrint(
    {
      WFdesign.sim(as.numeric(unlist(strsplit(input$truth,",",TRUE))),as.numeric(unlist(strsplit(input$fprob,",",TRUE))),input$tul,input$ell,input$puf,input$put,input$Nmax,input$Mmax,input$ntrial,input$start,input$rseed)
    }                                                                                                  
  )
  
  output$implementation <- renderPrint(
    {
     # validate(
     #   need(expr=sum(as.numeric(unlist(strsplit(input$y,",",TRUE))))>0, message='Need at least one DLT to generate DLT probability estimates.'),
     #   need(expr=as.numeric(unlist(strsplit(input$n,",",TRUE)))[1]>0, message='Each trial should begin at the lowest study dose level.')
     # )
      impbcdp(as.numeric(unlist(strsplit(input$y,",",TRUE))),as.numeric(unlist(strsplit(input$n,",",TRUE))),input$theta,input$dose.curr,input$cs)
    }
  )
  output$stopping <- renderPrint(
    {
      #validate(
      # need(expr=length(input$y) == length(input$n), message='Different number of dose levels specified in observed data')
      #)
      toxmonitoring(input$clevel,input$atul,input$trialsize)
    }
  )
  output$fstopping <- renderPrint(
    {
      #validate(
      # need(expr=length(input$y) == length(input$n), message='Different number of dose levels specified in observed data')
      #)
      feasmonitoring(input$flevel,input$aell,input$msize) #qmeans,flevel,aell,msize
    }
  )
}

shinyApp(ui = ui, server = server)


library(shiny)


ui<-(navbarPage("Design of Conaway, Dunbar, and Peddada (2004) for single-agent phase I trials",
                #title
                tabPanel("Simulation",
                         h4("Web Application for simulating operating characteristics of the CDP method for single-agents"),
                         #h4("Nolan A. Wages and Gina R. Petroni"),
                         h5("Division of Translational Research & Applied Statistics, University of Virginia;", a("nwages@virginia.edu", href="mailto:nwages@virginia.edu")),
                         br(),
                         h5("1. Enter an assumed set of true DLT probabilities, separated by commas.", strong("Note:"), "The length of this set should be equal to the number of possible study dose levels."),
                         #input: boxes for truth input values OR text field
                         textInput("p0","True DLT probability at each dose level",value="0.10,0.20,0.30,0.40,0.55",width=NULL),
                         
                         br(),
                         h5("2. Enter the target DLT probability that defines the MTD for the study."),
                         #input: tul
                         numericInput("target","Target DLT rate",0.20,0,1,0.01,NULL),
                         
                         br(),
                         h5("3. Enter the cohort size required before the next model-based update. Cohort size may be 1, 2, or 3 patients."),
                         #input: cohortsize
                         numericInput("cohortsize","Cohort size",3,1,3,1,NULL),  
                         
                         br(),
                         h5("4. Enter the maximum sample size for the study. This number should be a multiple of the cohort size entered above."),
                         #input: ncohort
                         numericInput("ssize","Maximum number of patients",30,0,NA,1,NULL),
                         
                         br(),
                         h5("5. Enter the total number of patients treated on any dose required to stop the trial. At any point in the trial, if the recommendation is to assign
                            the next cohort to a dose that already has the entered number of patients treated on the dose, the study
                            is stopped and the recommended dose is declared the MTD. If the entered number is larger than the maximum sample size, each trial will accrue to
                            the maximum sample size."),
                         #input: ncohort
                         numericInput("n.stop","Number of patients needed on one dose to stop",31,0,NA,1,NULL),
                         
                         br(),
                         h5("6. Enter the number of simulations. A minimum of 1000 is recommended."),
                         #input: ntrial
                         numericInput("ntrial","Number of simulated trials", 1000, 0, 100000, NA, NULL),
                         
                       #  br(),
                       #  p("7. Enter the index of the starting dose level.", strong("Note:"), "Index of", span("lowest", style = "color:blue"), "dose level is always 1. If the design allows for", em("'minus'"), "
                       #    dose levels (i.e. -2, -1, etc.), then the index of the starting dose should account for these lower levels (i.e. if -1 dose level allowed, starting dose is 2.)"),
                        # #Even  if lowest dose level is", em("labeled"), "-1, index of lowest dose level is 1. 
                       #  #input: start
                       #  numericInput("start","Index of starting dose level",1,1,NA,1,NULL),
                         
                         br(),
                         h5("7. Set the seed of the random number generator."),
                         #input: rseed
                         numericInput("rseed","Random seed",342425,0,NA,1,NULL),
                       
                       br(),
                       h5("8. Specify the threshold for defining the safety stopping rule based on an unacceptable high DLT rate at the lowest study dose level."),
                       #input: tul
                       numericInput("cl","Threshold for defining safety stopping rule",0.90,0,1,0.01,NULL),
                         
                         #input: submit button
                         submitButton("Run simulation study",icon("flask"),NULL),
                         
                         #output: simbcrm(truth,target,cohortsize,ncohort,ntrial,start,n.stop,rseed)
                         verbatimTextOutput("simulation"),
                         
                         mainPanel(
                           p(strong("This application simulates operating characteristics for the CDP method [1] with each simulated trial beginning at the lowest study dose level.")),
                          #p("1.", strong("Skipping Restriction:"), "The trial is not allowed to skip dose levels when escalating."),
                           #p("2.", strong("Skeleton:"), "For the specified target DLT rate and total number of dose levels, the skeleton of power model d^exp(a) is generated according to Lee and Cheung (2009) [2] using a prior MTD located at the median dose level and a spacing measure of delta=0.05."),
                           #p("1.", strong("Prior:"), "The prior distribution is beta(a,b) chosen according to Conaway, Dunbar, and Peddada [1]."),
                           #p("4.", strong("Accuracy Index:"), "The Accuracy Index is given in Equation 6.1 in Cheung [3]."),
                           #p("2.", strong("Safety Stopping Rule:"), "Stop the trial for safety if the Pr(DLT rate at dose level 1 > target | data) > 0.95"),
                          #p("1.", strong("Starting dose level:"), "Each simulated trial begins at the lowest study dose level."),
                           # p("Click",a("here",target="_blank",href="http://faculty.virginia.edu/model-based_dose-finding/nonparametric%20benchmark.R"),"to view the R Code for the benchmark."),
                           #p("Click",a("here", target="_blank", href="https://github.com/graham-wheeler/AplusB"),"to download the AplusB application to your computer from GitHub (set-up instructions provided)."),
                           strong("References:"),
                           p("[1] Conaway MR, Dunbar S, Peddada S (2004). Designs for single- and multiple-agent phase I trials, ", em("Biometrics;"), strong("60"),": 661-669.")
                           #p("[2] Lee and Cheung (2009). Model calibration in the continual reassessment method, ", em("Clinical Trials;"), strong("6"),"(3): 227-238."),
                           #p("[3] Lee and Cheung (2011). Calibration of prior variance in the bayesian continual reassessment method, ", em("Statistics in Medicine;"), strong("30"),"(17): 2081-2089.")
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
                         h3("Web Application for implementation of the CDP method"),
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
                           p(strong("This application computes a recommended dose level for the next patient in a phase I trial according to the CDP method [1] for a trial beginning at the lowest study dose level.")),
                           #p("1.", strong("Skipping Restriction:"), "The trial is not allowed to skip dose levels when escalating."),
                          # p("2.", strong("Skeleton:"), "For the specified target DLT rate and total number of dose levels, the  skeleton of power model d^exp(a) is generated according to Lee and Cheung (2009) [2] using a prior MTD located at the median dose level and a spacing measure of delta=0.05."),
                        #  p("1.", strong("Prior:"), "The prior distribution is beta(a,b) chosen according to Conaway, Dunbar, and Peddada [1]."),
                          #p("4.", strong("Accuracy Index:"), "The Accuracy Index is given in Equation 6.1 in Cheung [3]."),
                          #p("2.", strong("Safety Stopping Rule:"), "Stop the trial for safety if the Pr(DLT rate at dose level 1 > target | data) > 0.95"),
                        #  p("3.", strong("Starting dose level:"), "Each trial should begin at the lowest study dose level."),
                          # p("Click",a("here",target="_blank",href="http://faculty.virginia.edu/model-based_dose-finding/nonparametric%20b
                           strong("References:"),
                          p("[1] Conaway MR, Dunbar S, Peddada S (2004). Designs for single- and multiple-agent phase I trials, ", em("Biometrics;"), strong("60"),": 661-669.")
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
                         
                         
                )
                )
     
)



server <- function(input, output) {
  simcdp<-function (p0,target,cohortsize,ssize,n.stop,ntrial,rseed,cl){
    library(Iso)
    set.seed(rseed)    
    ndose=length(p0)
    
    ###calculate prior
    x<-2*target
    mu<-target
    u<-0.95
    
    f<-function(b){
      pbeta(x,mu*b/(1-mu),b)-u
    }
    b0<-uniroot(f,c(0.0001,100))$root
    a0<-mu*b0/(1-mu)
    ncohort=ssize/cohortsize
    
    
    newcdp<-function(p0,target,a0,b0,cohortsize,ncohort,n.stop,cl){
     
      #truth,skeleton,target,cohortsize,ncohort,x1,n.stop
      
      
      
      ### run a trial 	
      ndose = length(p0);   #number of doses
      y=numeric(ndose);  #number of toxicity/responses at each dose level
      n=numeric(ndose);  #number of treated patients at each dose level
      dose.select=rep(0,ndose); # a vector of indicators for dose selection
      stop=0; #indicate if trial stops early
      curr = 1
      i=1
      while(i <= ncohort){
        y[curr] = y[curr] + rbinom(1,cohortsize,p0[curr]);
        n[curr] = n[curr] + cohortsize;
        tried=which(n>0)
        
     
        safety=ifelse(n[1]>1,1 - pbeta(target, y[1] + a0, n[1] - y[1] + b0),0)
        
        if(safety>cl){
          stop=1
          break
        }
        
        u=(y[tried]+a0)/(n[tried]+a0+b0)
        pipost=pava(u,w=n)
        lossvec=abs(pipost-target)  
        T=lossvec==min(lossvec)
        poss=which(T)
        if(sum(T)==1){
          sugglev=poss
        } else {
          if(all(pipost[poss]>target)){
            sugglev=min(poss)
          } else {
            sugglev=max(poss)
          }
        }
        
        
        if(length(tried)<ndose){
          if(pipost[sugglev]<target){
            curr=ifelse(n[sugglev+1]==0,sugglev+1,sugglev)					
          } else {
            curr=sugglev
          }
        } else {
          curr=sugglev
        }
        
        if(n[curr]>=n.stop){
          stop<-0
          break
        }
        
        i<-i+1
      }		
      if(stop==0){
        dose.select[curr]=dose.select[curr]+1;
      }
      return(list(dose.select=dose.select,tox.data=y,pt.allocation=n,stop=stop))
    }
    d.select<-tox<-pts<-matrix(nrow=ntrial,ncol=ndose)
    nstop=0
    
    for(i in 1:ntrial){
      result<-newcdp(p0,target,a0,b0,cohortsize,ncohort,n.stop,cl)
      d.select[i,]=result$dose.select
      tox[i,]=result$tox.data
      pts[i,]=result$pt.allocation
      nstop=nstop+result$stop
    }
    cat("True DLT probability:             ", round(p0,3), sep="\t",  "\n");
    cat("MTD selection percentage:         ", formatC(colMeans(d.select)*100, digits=1, format="f"), sep="\t",  "\n");
    cat("Average number of DLTs:           ", formatC(colMeans(tox), digits=1, format="f"), sep="\t",   "\n");
    cat("Average number of patients:       ", formatC(colMeans(pts), digits=2, format="f"), sep="\t",   "\n");
    #cat("Accuracy index:                   ", round(1-length(truth)*(sum(abs(truth-target)*(colMeans(d.select)))/sum(abs(truth-target))),4), sep="\t",  "\n");
    cat("Percentage stopped for safety:    ", round(nstop/ntrial*100,3),sep="\t", "\n\n");
    cat("Design specifications: \n");
    cat("Prior distribution on DLT rate at each dose level:\n");
    cat("beta(",round(a0,2),",",round(b0,2),")","\n");
    cat("Safety stopping rule:\n");
    cat("Stop the trial if the Pr(DLT rate at the lowest study dose level > target | data) > ",cl,"\n");
    
} #end of MCdesign
  
  
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
    
    cat("Date and time:                           ", format(Sys.time()), sep="\t",  "\n");
    cat("Number of DLTs:        	                ", y, sep="\t",   "\n");
    cat("Number of patients evaluated for DLT:    ", n, sep="\t",   "\n");
    cat("Estimated DLT probabilities:             ", round(pipost[tried],2), sep="\t",  "\n");
    cat("Target DLT rate:                         ", theta, sep="\t",  "\n");
    cat("Recommended dose level:                  ", dose.curr, sep="\t",  "\n\n");
    cat("Design specifications: \n");
    cat("Prior distribution on DLT rate at each dose level:\n");
    cat("beta(",round(a0,2),",",round(b0,2),")","\n");
    cat("Safety stopping rule:\n");
    cat("Stop the trial if the Pr(DLT rate at the lowest study dose level > target | data) > ",cs,"\n");
  }
  
  toxmonitoring<-function(clevel,tul,trialsize){
    ###calculate prior
    x<-2*tul
    mu<-tul
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
        lv[k]<-lv[k]<-as.numeric(1-pbeta(tul,a0+k,b0+i-k)>clevel)
      }
      mintox[i]<-ifelse(all(lv==0),0,min(which(lv==1)))
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
  
  #reactive function to use inputs and produce output of functions
  output$simulation <- renderPrint(
    {
      validate(
        need(expr=input$ssize/input$cohortsize == round(input$ssize/input$cohortsize), message='Total number of patients is not multiple of cohort size!')
      )
            simcdp(as.numeric(unlist(strsplit(input$p0,",",TRUE))),input$target,input$cohortsize,input$ssize,input$n.stop,input$ntrial,input$rseed,input$cl)
    }                                                                             #p0,target,cohortsize,ssize,n.stop,ntrial,start,rseed
  )
  
  output$implementation <- renderPrint(
    {
      validate(
        need(expr=sum(as.numeric(unlist(strsplit(input$y,",",TRUE))))>0, message='Need at least one DLT to generate DLT probability estimates.'),
        need(expr=as.numeric(unlist(strsplit(input$n,",",TRUE)))[1]>0, message='Each trial should begin at the lowest study dose level.')
      )
      impbcdp(as.numeric(unlist(strsplit(input$y,",",TRUE))),as.numeric(unlist(strsplit(input$n,",",TRUE))),input$theta,input$dose.curr,input$cs)
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
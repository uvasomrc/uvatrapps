
library(shiny)
library(nnet)

ui <- fluidPage(
  #title
  titlePanel("Application for evaluating Phase I methods using a non-parametric optimal benchmark ", windowTitle="Non-parametric Optimal Design"),
  h4("Nolan A. Wages and Nikole Varhegyi"),
  h4("Division of Translational Research & Applied Statistics, University of Virginia;", a("nwages@virginia.edu", href="mailto:nwages@virginia.edu")),
  
  br(),
  h5("1. Enter a set of assumed true DLT probabilities, separated by commas. The length of this set should be equal to the number of dose levels."),
  #input: boxes for r input values OR text field
  textInput("r","True DLT probability at each dose level",value="0.05,0.07,0.20,0.35,0.55,0.70",width=NULL),
  
  br(),
  h5("2. Enter the target DLT rate for the study."),
  #input: theta
  numericInput("theta","Target DLT rate",0.20,0,1,0.01,NULL),
  
  br(),
  h5("3. Enter the sample size for the study."),
  #input: n
  numericInput("n","Sample size",20,0,NA,1,NULL),
  
  br(),
  h5("4. Enter the number of simulated trials to be generated. A minimum of 1000 is recommended."),
  #input: B
  numericInput("B","Number of simulated trials", 2000, 0, 100000, NA, NULL),
  
  br(),
  h5("5. Set the seed of the random number generator."),
  #input: rseed
  numericInput("rseed","Random seed",580,0,NA,1,NULL),
  
  #input: submit button
  submitButton("Run simulation",icon("flask"),NULL),
  
  #output: npoptimal(r,theta,n,B)
  verbatimTextOutput("result"),
  
  mainPanel(
    p(strong("This application provides operating characteristics for accuracy of the non-parametric optimal benchmark [1].")),
    p("The Accuracy Index is given in Equation 6.1 in Cheung [2]."),
    p("Click",a("here",target="_blank",href="http://faculty.virginia.edu/model-based_dose-finding/nonparametric%20benchmark.R"),"to view the R Code for the benchmark."),
    #p("Click",a("here", target="_blank", href="https://github.com/graham-wheeler/AplusB"),"to download the AplusB application to your computer from GitHub (set-up instructions provided)."),
    strong("References:"),
    p("[1] O'Quigley J, Paoletti X, Maccario J. (2002). Non-parametric optimal design in dose-finding studies, ", em("Biostatistics;"), strong("3"),"(1): 51-56."),
    p("[2] Cheung YK (2011). ", em("Dose-finding by the continual reassessment method."), "Chapman and Hall/CRC press: New York.")
  )
)

server <- function(input, output) {
  npoptimal <- function(r, theta,n,B,rseed){
    
    set.seed(rseed) 
    select=rep(0,length(r))
    j=1
    while(j<=B){
      latent=runif(n,0,1)
      tox=matrix(nrow=n,ncol=length(r))
      i=1
      while(i<=n){
        tox[i,]<-as.numeric(latent[i]<=r)
        i=i+1
      }
      loss=abs(colMeans(tox)-theta)
      sugglev=which(loss==min(loss))
      mtd=ifelse(length(sugglev)==1,sugglev,sample(sugglev,1))
      select[mtd]=select[mtd]+1
      j=j+1
    }
    cat("True DLT probability:     ", format(round(r,3), nsmall = 2), sep="\t",  "\n");
    cat("MTD selection percentage: ", formatC((select/B)*100, digits=1, format="f"), sep="\t",  "\n");
    cat("Accuracy Index:           ", round(1-length(r)*(sum(abs(r-theta)*(select/B))/sum(abs(r-theta))),4), sep="\t",  "\n");
  }
  
  #reactive function to use inputs and produce output of npoptimal
  output$result <- renderPrint(
    {
      npoptimal(as.numeric(unlist(strsplit(input$r,",",TRUE))),input$theta,input$n,input$B,input$rseed)
    }
  )
}

shinyApp(ui = ui, server = server)
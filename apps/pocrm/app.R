library(shiny)
library(rhandsontable)
library(shinyjqui)
library(shinyjs)
library(pocrm)
library(shinyBS)

sim_jsResetCode <- "shinyjs.sim_reset = function() {history.go(0)}"
imp_jsResetCode <- "shinyjs.imp_reset = function() {history.go(0)}"

server <- function(input, output, session) {
  
  ###########################
  ## help text
  addPopover(session=session, 
             id="info_sim_orders", 
             title="Orders", 
             content="<p>All possible orderings of the combinations with regards to their toxicity probabilities. Used to generate the matrix of skeleton values corresponding to the possible orderings of the toxicity probabilities specified by orders.<button type='button' id='close' class='close'onclick='$(&quot;#info_sim_orders&quot;).popover(&quot;hide&quot;);'>&times;</button></p>",
             placement = "bottom",
             trigger = "click", 
             options = NULL)
  
  ###########################
  ## help text
  addPopover(session=session, 
             id="info_sim_inputs", 
             title="Simulation Details", 
             content="<p>
             <b>Target DLT Rate</b>
             <br>
             Acceptable DLT probability that defines the study MTD
             <br>
             <b>Acceptable Toxicity Range</b>
             <br>
             A single numeric value used to define a range of 'acceptable' DLT rates. The simulation results will report the percentage of simulated trials that recommended a combination within +/- tox.range of the target rate.
             <br>
             <b>Maximum Sample Size</b>
             <br>
             The maximum number of participants to be accrued
             <br>
             <b>Observe Stopping Rule</b>
             <br>
             Whether or not to the trial will stop once a certain number of participants have been treated on any combination
             <br>
             <b>Cohort Size</b>
             <br>
             Stage 1 escalation proceeds according to the first ordering entered into the set of possible orderings in cohorts of 1, 2, or 3 patients as specified by the user. This scheme will only be used until the first DLT is observed, at which time the design switches to Stage 2 model-based allocation. In Stage 2, the cohort size is 1.
             <br>
             <b>Number of Simulations</b>
             <br>
             The number of simulated trials. Recommended minimum is 1000.
             <br>
             <b>Safety stopping rule</b>
             <br>
             Stop the trial for safety if the lower limit of an 80% Agresti-Coull (1998) binomial confidence interval for the lowest combination exceeds the target DLT rate.
             <button type='button' id='close' class='close'onclick='$(&quot;#info_sim_inputs&quot;).popover(&quot;hide&quot;);'>&times;</button></p>",
             placement = "bottom",
             trigger = "click", 
             options = NULL)
  
  ###########################
  ## help text
  addPopover(session=session, 
             id="info_imp_combo", 
             title="Combinations", 
             content="<p>Labels for the possible combinations of drugs and doses.<button type='button' id='close' class='close'onclick='$(&quot;#info_imp_combo&quot;).popover(&quot;hide&quot;);'>&times;</button></p>",
             placement = "bottom",
             trigger = "click", 
             options = NULL)
  
  ###########################
  ## help text
  addPopover(session=session, 
             id="info_imp_orders", 
             title="Orders", 
             content="<p>All possible orderings of the combinations with regards to their toxicity probabilities. Used to generate the matrix of skeleton values corresponding to the possible orderings of the toxicity probabilities specified by orders.<button type='button' id='close' class='close'onclick='$(&quot;#info_imp_orders&quot;).popover(&quot;hide&quot;);'>&times;</button></p>",
             placement = "bottom",
             trigger = "click", 
             options = NULL)
  
  ###########################
  ## help text
  addPopover(session=session, 
             id="info_imp_inputs", 
             title="Implementation Details", 
             content="<p>
             <b>Target DLT Rate</b>
             <br>
             Acceptable DLT probability that defines the study MTD
             <br>
             <b>Observed Trial Data</b>
             <br>
             Comma-separated value (csv) file with observed trial data. The file must have two columns: the first should contain the labels for combinations (see 'Combination Labels') that have been tried to this point in the study, and the second should have indicators of DLT outcomes (Yes=1,No=0).
             <br>
             Note: please include a header in the first row and begin data entry on the second row.
             <button type='button' id='close' class='close'onclick='$(&quot;#info_imp_inputs&quot;).popover(&quot;hide&quot;);'>&times;</button></p>",
             placement = "bottom",
             trigger = "click", 
             options = NULL)
  
  ###########################
  ## simulation
  sim_vals <- 
    reactiveValues(
      probs = NULL,
      orders_data = NULL)
  
  observeEvent(input$sim_reset, 
               {js$sim_reset()})
  
  observeEvent(input$simulate,
               {shinyjs::show("sim-results-label")})
  
  observeEvent(input$sim_dims,   
               # need to wrap in isolate so only triggers once
               isolate(
                 insertUI(
                   selector = '#initial_simorder',
                   ## wrap element in a div with id for ease of removal
                   ui = tags$div(
                     orderInput("simorder", 
                                "", 
                                items = 1:(input$sim_ndoses_a*input$sim_ndoses_b), 
                                width = "auto"),
                     tags$br(),
                     id = "initial_simorder")
                 ) 
               )
  )
  
  observeEvent(input$sim_insertBtn, {
    btn <- input$sim_insertBtn
    id <- paste0('simorder', btn)
    insertUI(
      selector = '#sim_placeholder',
      ## wrap element in a div with id for ease of removal
      ui = tags$div(
        orderInput(id, "", 
                   items = 1:(input$sim_ndoses_a*input$sim_ndoses_b), 
                   width = "auto"), 
        tags$br(),
        id = id
      )
    )
    
  })
  
  observe({
    
    # capture all reactive values and look for ones that are order inputs
    orderlist <- reactiveValuesToList(input)
    
    # ind <- grepl("_order", names(orderlist))
    ind <- grepl("simorder[0-9]*_order", names(orderlist))
    
    orderlist <- orderlist[ind]
    
    # bind the values together
    orders <- do.call(rbind, orderlist)
    
    # make sure no duplicated orders are used
    orders <-
      orders[!duplicated(orders),]
    
    # for some reasone subsetting above coerces single set of orders to vector
    # make sure this is treated as a matrix ...
    # ... otherwise getwm() will not work
    if(is.vector(orders)) {
      
      orders <- matrix(orders,  ncol = length(orders))
      
    }
    
    sim_vals$orders_data <- orders
    
  })
  
  sim_ndoses <- eventReactive(input$sim_dims,{
    
    list(a = input$sim_ndoses_a,
         b = input$sim_ndoses_b)
    
  })
  
  # disable ndoses inputs after triggering dims
  
  observeEvent(input$sim_dims, {
    
    shinyjs::disable("sim_ndoses_a")
    shinyjs::disable("sim_ndoses_b")
    shinyjs::show("sim_other_inputs")
    shinyjs::show("sim_headers_help")
    shinyjs::hide("sim_dims")
    shinyjs::show("sim_reset")
    
  })
  
  # handle showing stopn
  
  observe({
    if(input$stop) {
      shinyjs::show("stopn")
    } else {
      shinyjs::hide("stopn")
    }
  })
  
  output$sim_hot <- renderRHandsontable({
    
    DF <- matrix(0, 
                 ncol = sim_ndoses()$b,
                 nrow = sim_ndoses()$a)
    
    row.names(DF) <- paste0("A", sim_ndoses()$a:1)
    colnames(DF) <- paste0("B", 1:sim_ndoses()$b)
    
    rhandsontable(DF)
    
  })
  
  output$sim_matrix <- renderTable({
    
    len <- sim_ndoses()$b * sim_ndoses()$a
    
    DF <- matrix(1:len,
                 ncol = sim_ndoses()$b, 
                 nrow = sim_ndoses()$a,
                 byrow = TRUE,
                 dimnames = 
                   list(
                     paste0("A", 1:sim_ndoses()$a),
                     paste0("B", 1:sim_ndoses()$b)
                   ))
    
    DF[rev(1:nrow(DF)),]
    
  },
  rownames = TRUE,
  colnames = TRUE)
  
  output$sim_sum <- renderPrint({
    
    req(input$simulate)
    
    out <-
      paste0(
        "Overall DLT Proportion: ",
        sim_res()$fit$percent.DLT,
        "\n",
        "Proportion of trials stopped for safety: ",
        sim_res()$fit$stop.safety,
        "\n",
        "Average Sample Size: ",
        sim_res()$fit$mean.n,
        "\n",
        "Proportion Acceptable MTD Selection: ",
        sim_res()$fit$acceptable,
        "\n",
        "Skeleton used to construct working models: ",
        paste(round(sim_res()$skeleton, 3), collapse = ","),
        "\n",
        "Safety stopping bounds:\n "
        # print(sim_res()$fit$stop.bounds)
      )
    
    
    cat(out)
    print(sim_res()$fit$stop.bounds)
    
  })
  
  output$sim_res_df <- renderTable({
    
    req(input$simulate)
    
    df <- sim_res()$df
    
    names(df) <- c("Combinations",
                   "True DLT\nProbability",
                   "MTD Selection\nProportion",
                   "Patient Allocation\nProportion")
    
    df
    
  })
  
  sim_res <- eventReactive(input$simulate,{
    
    # parsing dlt probs
    # need to get data (vector) from handsontable ...
    # then convert to matrix ...
    # ensure matrix is reversed from data input (so that vector is in order)... 
    # transpose (to get the vector in row order) ...
    # convert to vector
    r <- unlist(isolate(input$sim_hot)$data)
    r <- matrix(r,
                ncol = input$sim_ndoses_b, 
                nrow = input$sim_ndoses_a, 
                byrow = TRUE)
    r <- r[rev(1:nrow(r)),]
    r <- as.vector(t(r))
    
    # get orders generated and stored in reactive value (sim_vals)
    orders <-
      sim_vals$orders_data
    
    # make sure orders are ordered in sequence they were input
    if(!is.null(rownames(orders))) {
      orders <-
        orders[order(rownames(orders)),]
    }
    
    # double check that orders are captured as numeric
    class(orders) <- "numeric"
    
    # how many combinations of dosages are there?
    nlevel <- sim_ndoses()$a*sim_ndoses()$b
    
    #Specify the skeleton values.
    skeleton <- getprior(halfwidth = 0.05,
                         target = input$sim_target,
                         nu = floor(nlevel/2),
                         nlevel = nlevel)
    
    #Initial guesses of toxicity probabilities for each ordering.
    alpha <- getwm(orders,skeleton)
    
    #We consider all orders to be equally likely prior to the study.
    norders <- nrow(orders)
    prior.o <- rep(1/norders,norders)
    
    # set up cohort length using initial order
    vec <- orders[1,]
    x0 <- unlist(lapply(vec, function(x) rep(x, input$cohort)))
    x0 <- as.numeric(x0)
    
    # maximum sample size
    n <- input$samplesize
    
    # number of patients used to define stopping rule
    
    if(input$stop) {
      
      stop <- input$stopn
      
    } else {
      
      stop <- n + 1
      
    }
    
    # stop <- n + 1
    
    #The target toxicity rate
    theta <- input$sim_target
    
    #Number of simulations
    nsim <- input$nsim
    
    #Definition of acceptable DLT rates
    tox.range <- input$toxrange
    
    set.seed(input$seed)
    
    pocrm.sim<-function(r,alpha,prior.o,x0,stop,n,theta,nsim,tox.range){
      
      library(nnet)
      library(binom)
      
      toxmonitoring<-function(theta,n){
        mintox<-rep(0,n)
        for(i in 1:n){
          lv=rep(0,i)
          for(k in 1:i){
            lv[k]<-as.numeric(binom.confint(k,i,conf.level=0.8,methods="agresti-coull")$lower>theta)
          }
          mintox[i]<-ifelse(all(lv==0),0,min(which(lv==1)))
        }
        df=data.frame(mintox[2:n],2:n)
        names(df)=c("#DLTs","#pts")
        # cat("Stop the study for safety if the observed DLT rate at lowest study dose level >= #DLTs out of #pts treated at lowest study dose level. \n\n");
        #print(df,row.names=FALSE);
        df
      }
      
      ###Load the function 'lpocrm' 
      lpocrm<-function(r,alpha,prior.o,x0,stop,n,theta){
        
        # if a single ordering is inputed as a vector, convert it to a matrix
        if(is.vector(alpha)) alpha=t(as.matrix(alpha));
        
        nord.tox = nrow(alpha);
        mprior.tox = prior.o;  # prior for each toxicity ordering
        
        bcrml<-function(a,p1,y,n){
          lik=0
          for(j in 1:length(p1)){
            lik=lik+y[j]*a*log(p1[j])+(n[j]-y[j])*log((1-p1[j]**a));
          }
          return(lik);
        }
        
        ### run a trial 	
        ncomb = ncol(alpha);   #number of combos
        y=npts=ptox.hat=comb.select=numeric(ncomb);  
        comb.curr = x0[1];  # current dose level	 
        stoprule=0; #indicate if trial stops early
        i=1
        
        stage1<-c(x0,rep(ncol(alpha),n-length(x0)))
        
        ##Stage 1
        while(i <= n){
          y[comb.curr] = y[comb.curr] + rbinom(1,1,r[comb.curr]);
          npts[comb.curr] = npts[comb.curr] + 1;
          
          if(sum(y)==sum(npts)){
            safety=ifelse(npts[1]>1,binom.confint(y[1],npts[1],conf.level=.9,methods="agresti-coull")$lower,0)
            
            if(safety>theta){
              stoprule=1
              break
            }
            
            comb.curr<-ifelse(comb.curr==1,comb.curr,comb.curr-1)
          } else if(sum(y)==0){
            comb.curr<-ifelse(comb.curr==ncomb,comb.curr,stage1[i+1])
          } else {
            break
          }
          if(any(npts>stop)){
            stoprule<-0
            break
          }
          i=i+1
        }
        
        #Stage 2
        while(sum(npts) <= n)
        {
          if(sum(y)==0){
            stop=0
            break
          } else{
            like.tox= est.tox=rep(0, nord.tox);
            for(k in 1:nord.tox)
            {
              est.tox[k]<-optimize(f=bcrml,interval=c(0,100),p1=alpha[k,],y=y,n=npts,maximum=T)$maximum
              like.tox[k]<-optimize(f=bcrml,interval=c(0,100),p1=alpha[k,],y=y,n=npts,maximum=T)$objective
            }		
            
            postprob.tox = (exp(like.tox)*mprior.tox)/sum(exp(like.tox)*mprior.tox);
            # toxicity model selection, identify the model with the highest posterior prob
            if(nord.tox>1){ 
              mtox.sel = which.is.max(postprob.tox); 
            } else{
              mtox.sel = 1;
            }
            
            ptox.hat=alpha[mtox.sel,]**est.tox[mtox.sel]
            
            safety=ifelse(npts[1]>1,binom.confint(y[1],npts[1],conf.level=.9,methods="agresti-coull")$lower,0)
            
            if(safety>theta){
              stoprule=1
              break
            }
            
            
            loss=abs(ptox.hat-theta)
            comb.curr=which.is.max(-loss)
            if(npts[comb.curr]==stop){
              stoprule<-0
              break
            }
            
            if(sum(npts)==n){
              stoprule=0
              break
            } else{
              # generate data for a new cohort of patients
              y[comb.curr] = y[comb.curr] + rbinom(1,1,r[comb.curr]);
              npts[comb.curr] = npts[comb.curr] + 1;
            }
          }
        }
        if(stoprule==0){
          comb.select[comb.curr]=comb.select[comb.curr]+1;
        }
        return(list(MTD.selection=comb.select,tox.data=y,patient.allocation=npts,stoprule=stoprule))
      }
      ##########'lpocrm' end here
      
      
      
      ###Load the function 'lpocrm.sim' 
      lpocrm.sim<-function(nsim){
        ncomb=length(r)
        
        comb.select<-y<-npts<-matrix(nrow=nsim,ncol=ncomb)
        trialsize<-rep(0,nsim)
        
        for(i in 1:nsim){
          result<-lpocrm(r,alpha,prior.o,x0,stop,n,theta)
          comb.select[i,]=result$MTD.selection
          y[i,]=result$tox.data
          npts[i,]=result$patient.allocation
          trialsize[i]=sum(result$patient.allocation)
        }
        return(list(true.prob=r,MTD.selection=round(colMeans(comb.select),3),patient.allocation=round(colMeans(npts)/mean(trialsize),3),percent.DLT=round(sum(colMeans(y))/mean(trialsize),3),mean.n=round(mean(trialsize),3),acceptable=sum(colMeans(comb.select)[which(round(abs(r-theta),2)<=tox.range)]),
                    stop.safety=1-sum(round(colMeans(comb.select),2)),stop.bounds=toxmonitoring(theta,n)))
      }
      lpocrm.sim(nsim)
      
    }
    
    
    
    
    fit <- pocrm.sim(r,
                     alpha,
                     prior.o,
                     x0,
                     stop,
                     n,
                     theta,
                     nsim,
                     tox.range)
    
    df <-
      data.frame(
        levels = 1:length(r),
        DLT = r,
        MTD = fit$MTD.selection,
        allocation = fit$patient.allocation
      )
    
    list(fit = fit,
         orders = orders,
         df = df,
         skeleton = skeleton)
    
    
  })
  
  ###########################
  ## implementation
  imp_vals <- 
    reactiveValues(
      probs = NULL,
      orders_data = NULL)
  
  observeEvent(input$imp_reset, 
               {js$imp_reset()})
  
  observeEvent(input$implement,
               {shinyjs::show("imp-results-label")})
  
  observeEvent(input$imp_dims,   
               # need to wrap in isolate so only triggers once
               isolate(
                 insertUI(
                   selector = "#initial_imporder",
                   ## wrap element in a div with id for ease of removal
                   ui = tags$div(
                     orderInput("imporder", 
                                "", 
                                items = 1:(input$imp_ndoses_a*input$imp_ndoses_b), 
                                width = "auto"),
                     tags$br(),
                     id = "initial_imporder")
                 ) 
               )
  )
  
  observeEvent(input$imp_insertBtn, {
    btn <- input$imp_insertBtn
    id <- paste0('imporder', btn)
    insertUI(
      selector = '#imp_placeholder',
      ## wrap element in a div with id for ease of removal
      ui = tags$div(
        orderInput(id, "", 
                   items = 1:(input$imp_ndoses_a*input$imp_ndoses_b), 
                   width = "auto"), 
        tags$br(),
        id = id
      )
    )
    
    
  })
  
  observe({
    
    # capture all reactive values and look for ones that are order inputs
    orderlist <- reactiveValuesToList(input)
    
    ind <- grepl("imporder[0-9]*_order", names(orderlist))
    
    
    orderlist <- orderlist[ind]
    
    # bind the values together
    orders <- do.call(rbind, orderlist)
    
    # make sure no duplicated orders are used
    orders <-
      orders[!duplicated(orders),]
    
    # for some reasone subsetting above coerces single set of orders to vector
    # make sure this is treated as a matrix ...
    # ... otherwise getwm() will not work
    if(is.vector(orders)) {
      
      orders <- matrix(orders,  ncol = length(orders))
      
    }
    
    imp_vals$orders_data <- orders
    
  })
  
  imp_ndoses <- eventReactive(input$imp_dims,{
    
    list(a = input$imp_ndoses_a,
         b = input$imp_ndoses_b)
    
  })
  
  # disable ndoses inputs after triggering dims
  
  observeEvent(input$imp_dims, {
    
    shinyjs::disable("imp_ndoses_a")
    shinyjs::disable("imp_ndoses_b")
    shinyjs::show("imp_other_inputs")
    shinyjs::hide("imp_dims")
    shinyjs::show("imp_reset")
    shinyjs::show("imp_headers_help")
    
    
  })
  
  output$imp_matrix <- renderTable({
    
    len <- imp_ndoses()$b * imp_ndoses()$a
    
    DF <- matrix(1:len,
                 ncol = imp_ndoses()$b, 
                 nrow = imp_ndoses()$a,
                 byrow = TRUE,
                 dimnames = 
                   list(
                     paste0("A", 1:imp_ndoses()$a),
                     paste0("B", 1:imp_ndoses()$b)
                   ))
    
    DF[rev(1:nrow(DF)),]
    
  },
  rownames = TRUE,
  colnames = TRUE)
  
  output$imp_sum <- renderPrint({
    
    req(input$implement)
    
    out <- 
      paste0(
        "Ordering Probabilities: ",
        paste0(imp_res()$fit$ord.prob, collapse = ","),
        "\n",
        "Estimated Ordering: ",
        imp_res()$fit$order.est,
        "\n",
        "Model Parameter Estimate: ",
        imp_res()$fit$a.est,
        "\n",
        "Skeleton used to construct working models: ",
        paste(round(imp_res()$skeleton, 3), collapse = ","),
        "\n\n",
        "Combination Recommendation: ",
        imp_res()$fit$dose.rec
        
      )
    
    cat(out)
    
  })
  
  output$imp_res_df <- renderTable({
    
    req(input$implement)
    
    df <- imp_res()$df
    
    names(df) <- c("Combinations",
                   #"Skeleton",
                   "Estimated DLT\nProbability")
    
    df
    
  })
  
  imp_res  <- eventReactive(input$implement,{
    
    # get orders generated and stored in reactive value (imp_vals)
    orders <-
      imp_vals$orders_data
    
    # make sure orders are ordered in sequence they were input
    if(!is.null(rownames(orders))) {
      orders <-
        orders[order(rownames(orders)),]
    }
    
    # double check that orders are captured as numeric
    class(orders) <- "numeric"
    
    # how many combinations of dosages are there?
    nlevel <- imp_ndoses()$a*imp_ndoses()$b
    
    #Specify the skeleton values.
    skeleton <- getprior(halfwidth = 0.05,
                         target = input$imp_target,
                         nu = floor(nlevel/2),
                         nlevel = nlevel)
    
    
    #Initial guesses of toxicity probabilities for each ordering.
    alpha <- getwm(orders,skeleton)
    
    #We consider all orders to be equally likely prior to the study.
    norders <- nrow(orders)
    prior.o <- rep(1/norders,norders)
    
    #The target toxicity rate
    theta <- input$imp_target
    
    # get combos from csv
    
    combos_outcomes <-
      read.csv(input$imp_combos$datapath)
    
    #Combinations tried on the first 11 patients in Table 5.
    combos <- combos_outcomes[,1]
    
    #Toxicity outcomes on the first 11 patients in Table 5.
    y <- combos_outcomes[,2]
    
    fit <- 
      pocrm.imp(
        alpha,
        prior.o,
        theta,
        y,
        combos)
    
    df <-
      data.frame(
        levels = 1:length(orders[1,]),
        #skeleton = skeleton,
        `Estimated DLT Probability` = fit$ptox.est
      )
    
    list(fit = fit,
         orders = orders,
         skeleton = skeleton,
         df = df)
    
  })
  
  
}

ui <- 
  navbarPage(fluid = TRUE, title = "Partial Order Continual Reassessment Method",
             tabPanel(title = "Simulation",
                      useShinyjs(),
                      extendShinyjs(text = sim_jsResetCode),
                      
                      # fluidPage(
                      fluidRow(
                        column(4,
                               numericInput("sim_ndoses_a", 
                                            "Number of Dose Levels Drug A", 
                                            value = 2, 
                                            min = 1)
                               
                        ),
                        column(4,
                               numericInput("sim_ndoses_b", 
                                            "Number of Dose Levels Drug B", 
                                            value = 2, 
                                            min = 1)
                        ),
                        column(4,
                               tags$br(),
                               actionButton("sim_dims","Set Dimensions"),
                               shinyjs::hidden(actionButton("sim_reset", "Reset", icon = icon("refresh"))))
                      ),
                      shinyjs::hidden(fluidRow(
                        shinyjs::hidden(bsButton("renderButton", "Render")),
                        column(4,
                               tags$h3("True DLT probability at each combination")),
                        column(4,
                               tags$div(h3("Possible orderings of DLT probabilities",
                                           actionLink("info_sim_orders", "",
                                                      icon = icon("question-circle-o"))))),
                        column(4,
                               tags$div(h3("Specifications for simulation",
                                           actionLink("info_sim_inputs", "", 
                                                      icon = icon("question-circle-o"))))),
                        id = "sim_headers_help"
                      )),
                      shinyjs::hidden(fluidRow(
                        column(4,
                               rHandsontableOutput("sim_hot"),
                               tags$h3("Combination Labels"),
                               tableOutput("sim_matrix")
                        ),
                        column(4,
                               tags$div(id = "initial_simorder"),
                               tags$div(id = "sim_placeholder"),
                               tags$div(
                                 actionButton("sim_insertBtn", "Add orders", icon = icon("plus")))
                        ),
                        column(2,
                               sliderInput("sim_target", "Target DLT Rate", 
                                           min = 0, max = 1, value = 0.3),
                               sliderInput("toxrange", "Acceptable Toxicity Range", 
                                           min = 0, max = 1, value = 0.05),
                               actionButton("simulate", "Simulate", icon = icon("flask"))
                        ),
                        column(2,
                               numericInput("samplesize", "Maximum Sample Size", min = 1, max = 100000, value=60),
                               checkboxInput("stop", "Observe stopping rule", value = TRUE),
                               shinyjs::hidden(
                                 numericInput("stopn", "Number treated on any combination to stop", min = 1, max = 100000, value = 10)
                               ),
                               radioButtons("cohort", "Stage 1 Cohort Size", choices = 1:3, inline = TRUE),
                               numericInput("nsim", "Number of Simulations", min = 1, max = 100000, value = 10),
                               numericInput("seed", "Random Seed", value = sample(1:1e5, 1))),
                        id = "sim_other_inputs")
                      ),
                      shinyjs::hidden(
                        fluidRow(
                          tags$hr(),
                          column(6,
                                 tags$h3("Results")
                          ),
                          id = "sim-results-label"
                        )
                      ),
                      fluidRow(
                        column(6,
                               tableOutput("sim_res_df")
                        ),
                        column(6,
                               verbatimTextOutput("sim_sum"))
                      )
             ),
             tabPanel(title = "Implementation",
                      useShinyjs(),
                      extendShinyjs(text = imp_jsResetCode),
                      fluidRow(
                        column(4,
                               numericInput("imp_ndoses_a", 
                                            "Number of Dose Levels Drug A", 
                                            value = 2, 
                                            min = 1)
                               
                        ),
                        column(4,
                               numericInput("imp_ndoses_b", 
                                            "Number of Dose Levels Drug B", 
                                            value = 2, 
                                            min = 1)
                        ),
                        column(4,
                               tags$br(),
                               actionButton("imp_dims","Set Dimensions"),
                               shinyjs::hidden(actionButton("imp_reset", "Reset", icon = icon("refresh"))))
                      ),
                      shinyjs::hidden(fluidRow(
                        shinyjs::hidden(bsButton("renderButton", "Render")),
                        column(4,
                               tags$div(h3("Combination Labels",
                                           actionLink("info_imp_combo", "",
                                                      icon = icon("question-circle-o"))))),
                        column(4,
                               tags$div(h3("Possible orderings of DLT probabilities",
                                           actionLink("info_imp_orders", "",
                                                      icon = icon("question-circle-o"))))),
                        column(4,
                               tags$div(h3("Specifications for implementation",
                                           actionLink("info_imp_inputs", "", 
                                                      icon = icon("question-circle-o"))))),
                        id = "imp_headers_help"
                      )),
                      shinyjs::hidden(fluidRow(
                        column(4,
                               # tags$h3("Combination Labels"),
                               tableOutput("imp_matrix")
                        ),
                        column(4,
                               # tags$h3("Possible orderings of DLT probabilities"),
                               tags$div(id = "initial_imporder"),
                               tags$div(id = "imp_placeholder"),
                               tags$div(
                                 actionButton("imp_insertBtn", "Add orders", icon = icon("plus")))
                        ),
                        column(4,
                               sliderInput("imp_target", "Target DLT Rate",
                                           min = 0, max = 1, value = 0.3),
                               fileInput("imp_combos", "Observed Trial Data"),
                               actionButton("implement", "Get updated recommendation", icon = icon("flask"))
                        ),
                        id = "imp_other_inputs")
                      ),
                      shinyjs::hidden(
                        fluidRow(
                          tags$hr(),
                          column(6,
                                 tags$h3("Results")
                          ),
                          id = "imp-results-label"
                        )
                      ),
                      fluidRow(
                        column(6,
                               tableOutput("imp_res_df")),
                        column(6,
                               verbatimTextOutput("imp_sum"))
                        
                      )
             ),
             tabPanel("About",
                      includeMarkdown("about.md")
             )
             
  )

shinyApp(ui, server)
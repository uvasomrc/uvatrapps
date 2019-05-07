library(shiny)
library(rhandsontable)
library(shinyjqui)
library(shinyjs)
library(pocrm)

sim_jsResetCode <- "shinyjs.sim_reset = function() {history.go(0)}"

imp_jsResetCode <- "shinyjs.imp_reset = function() {history.go(0)}"

server <- function(input, output, session) {
  
  ###########################
  ## simulation
  sim_vals <- 
    reactiveValues(
      probs = NULL,
      orders_data = NULL)
  
  observeEvent(input$sim_reset, 
               {js$sim_reset()})
  
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
                                width = "200px"),
                     tags$br(),
                     id = "initial_simorder")
                 ) 
               )
  )
  ## keep track of elements inserted and not yet removed
  # sim_inserted <- c()
  
  observeEvent(input$sim_insertBtn, {
    btn <- input$sim_insertBtn
    id <- paste0('simorder', btn)
    insertUI(
      selector = '#sim_placeholder',
      ## wrap element in a div with id for ease of removal
      ui = tags$div(
        orderInput(id, "", 
                   items = 1:(input$sim_ndoses_a*input$sim_ndoses_b), 
                   width = "200px"), 
        tags$br(),
        id = id
      )
    )
    # sim_inserted <<- c(id, sim_inserted)
    
    
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
    shinyjs::hide("sim_dims")
    shinyjs::show("sim_reset")
    
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
    
    # DF <- matrix(len:1, 
    #              ncol = sim_ndoses()$b, 
    #              nrow = sim_ndoses()$a)
    
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
    
    mat <- unlist(isolate(input$sim_hot)$data)
    dlt <- matrix(mat, ncol = input$sim_ndoses_b, nrow = input$sim_ndoses_a, byrow = TRUE)
    list(dlt,sim_res()$orders,sim_res()$fit)
    
  })
  
  # output$sim_debug <- renderPrint({
  #   
  #   req(input$simulate)
  #   
  #   sim_vals$orders_data
  #   
  # })
  
  output$sim_res_df <- renderTable({
    
    req(input$simulate)
    
    sim_res()$df
    
  })
  
  sim_res <- eventReactive(input$simulate,{
    
    # get dlt probabilities (true toxicity rates) from handsontable
    # r <- unlist(isolate(input$sim_hot)$data)
    
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
    
    # if(input$stop) {
    #   
    #   stop <- input$stopn
    #   
    # } else {
    #   
    #   stop <- n + 1
    #   
    # }
    
    stop <- n + 1
    
    #The target toxicity rate
    theta <- input$sim_target
    
    #Number of simulations
    nsim <- input$nsim
    
    #Definition of acceptable DLT rates
    tox.range <- input$toxrange
    
    set.seed(input$seed)
    
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
         df = df)
    
  })
  
  
  
  ###########################
  ## implementation
  imp_vals <- 
    reactiveValues(
      probs = NULL,
      orders_data = NULL)
  
  observeEvent(input$imp_reset, 
               {js$imp_reset()})
  
  observeEvent(input$imp_dims,   
               # need to wrap in isolate so only triggers once
               isolate(
                 insertUI(
                   selector = "#initial_imporder",
                   ## wrap element in a div with id for ease of removal
                   ui = tags$div(
                     #orderInput("order",
                     orderInput("imporder", 
                                "", 
                                items = 1:(input$imp_ndoses_a*input$imp_ndoses_b), 
                                width = "200px"),
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
                   width = "200px"), 
        tags$br(),
        id = id
      )
    )
    
    
  })
  
  observe({
    
    # capture all reactive values and look for ones that are order inputs
    orderlist <- reactiveValuesToList(input)
    
    # ind <- grepl("_order", names(orderlist))
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
    
  })
  
  # output$imp_hot <- renderRHandsontable({
  #   
  #   DF <- matrix(0, 
  #                ncol = imp_ndoses()$b,
  #                nrow = imp_ndoses()$a)
  #   
  #   row.names(DF) <- paste0("A", imp_ndoses()$a:1)
  #   colnames(DF) <- paste0("B", imp_ndoses()$b:1)
  #   
  #   rhandsontable(DF)
  #   
  # })
  
  output$imp_matrix <- renderTable({
    
    len <- imp_ndoses()$b * imp_ndoses()$a
    
    # DF <- matrix(len:1,
    #              ncol = imp_ndoses()$b,
    #              nrow = imp_ndoses()$a,
    #              byrow = FALSE)
    # 
    # DF
    
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
    
    # mat <- unlist(isolate(input$imp_hot)$data)
    # dlt <- matrix(mat, ncol = input$imp_ndoses_b, nrow = input$imp_ndoses_a)
    list(imp_res()$orders,imp_res()$fit,imp_res()$skeleton)
    
  })
  
  output$imp_debug <- renderPrint({
    
    req(input$implement)
    
    imp_vals$orders_data
    
  })
  
  imp_res  <- eventReactive(input$implement,{
    
    # get orders generated and stored in reactive value (imp_vals)
    orders <-
      imp_vals$orders_data
    
    #Specify the skeleton values provided in Table 4.
    # skeleton <- unlist(isolate(input$imp_hot)$data)
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
    
    list(fit = fit,
         orders = orders,
         skeleton = skeleton)
    
  })
  
  
}

ui <- 
  navbarPage(fluid = TRUE, title = "POCRM",
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
                        column(4,
                               tags$h3("True DLT probability at each combination"),
                               rHandsontableOutput("sim_hot"),
                               tags$h3("Combination Labels"),
                               tableOutput("sim_matrix")
                        ),
                        column(4,
                               tags$h3("Possible orderings of DLT probabilities"),
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
                               radioButtons("cohort", "Cohort Size", choices = 1:3, inline = TRUE),
                               numericInput("nsim", "Number of Simulations", min = 1, max = 100000, value = 10),
                               numericInput("seed", "Random Seed", value = sample(1:1e5, 1))),
                        id = "sim_other_inputs")
                      ),
                      tableOutput("sim_res_df"),
                      verbatimTextOutput("sim_debug"),
                      verbatimTextOutput("sim_sum")
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
                        column(4,
                               tags$h3("Combination Labels"),
                               tableOutput("imp_matrix")
                        ),
                        column(4,
                               tags$h3("Possible orderings of DLT probabilities"),
                               tags$div(id = "initial_imporder"),
                               tags$div(id = "imp_placeholder"),
                               tags$div(
                                 actionButton("imp_insertBtn", "Add orders", icon = icon("plus")))
                        ),
                        column(2,
                               sliderInput("imp_target", "Target DLT Rate",
                                           min = 0, max = 1, value = 0.3),
                               actionButton("implement", "Get updated recommendation", icon = icon("flask"))
                        ),
                        column(2,
                               fileInput("imp_combos", "Observed Trial Data")),
                        id = "imp_other_inputs")
                      ),
                      verbatimTextOutput("imp_debug"),
                      verbatimTextOutput("imp_sum")
             )
             
  )

shinyApp(ui, server)
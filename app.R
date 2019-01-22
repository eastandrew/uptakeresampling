#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Resampling distribution of dietary uptake in birds"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("runs",
                     "Number of runs:",
                     min = 10,
                     max = 10000,
                     value = 100),
         sliderInput("mean",
                     "Lognormal Mean Concentration:",
                     min = 0.001,
                     max = 100,
                     value = 10)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$distPlot <- renderPlot({
     
     soildata <- rlnorm(n=100,meanlog=input$mean,sdlog=1)  # create distribution of soil concentrations (fake but log normal for the moment)
     amountexposed <- seq(0.1,0.9,0.1)             # percent of diet based on contaminated soil
     massdist <- seq(125,175,5)                    # birds of a range of masses
     seedbsaf <- c(0.11,0.45)                      # vector of seed bsaf data from D'Hollander
     insbsaf <- c(60,35,68,7.1)                    # vector of invert bsaf data from D'Hollander
     birdsoilresample <- c()                       # create vector for the final bird uptake data "resampling method"
     massselect <- c()                             # create other vectors for intermediate calculations and sampling
     firselect <- c()
     ratiobug2 <- c()
     ratioins2 <- c()
     ratioseed2 <- c()
     dirt <- c()
     for(i in 1:input$runs)                              # repeat resampling 1000 (or other obvs) times
     {
       ratiobug2[i] <- sample(c(0.25,0.5,0.75),1)        # sample from a proportion of diet as inverts
       ratioins2[i] <- ratiobug2[i]*(1-0.1)              # remove the soil ingestion rate from the percent of diet as bugs/seeds.  SIR is 0.1 here
       ratioseed2[i] <- (1-ratiobug2[i])*(1-0.1)         # from remainder of (1-SIR)-bug% calculate percent as seeds
       massselect[i] <- sample(massdist,1)               # pick a random bird's mass
       firselect[i] <- 0.301*massselect[i]^0.751         # calculate food ingestion rate
       dirt[i] <- sample(soildata,1)                     # pick a random concentration
       # below= picking random bsaf values and percent on contaminated soil and then combine all values to estimate uptake
       birdsoilresample[i] <- (((sample(seedbsaf,1)*dirt[i]*firselect[i]*ratioseed2[i])+(sample(insbsaf,1)*dirt[i]*firselect[i]*ratioins2[i])+(dirt[i]*firselect[i]*0.1))*sample(amountexposed,1))/massselect[i]
     }
     plot(density(log10(birdsoilresample)))       # plot log10 distribution of uptake values
     abline(v=median(log10(birdsoilresample)),lty=2)
     text(median(log10(birdsoilresample))+1,0.2,bquote(paste("Median=",.(round(median(log10(birdsoilresample)),3)))))
     
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

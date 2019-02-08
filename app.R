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
                  value = 10),
      sliderInput("sd",
                  "Lognormal Standard Deviation Concentration:",
                  min = 0.0001,
                  max = 10,
                  value = 1),
      sliderInput("Numsoil",
                  "Number of Soil Samples:",
                  min=2,
                  max=10000,
                  value=100),
      sliderInput("time",
                  "% Time on Cont. Site:",
                  min=1,
                  max=100,
                  value=50),
      sliderInput("maxmass",
                  "Largest Bird Mass:",
                  min=10,
                  max=250,
                  value=150),
      sliderInput("minmass",
                  "Smallest Bird Mass:",
                  min=10,
                  max=250,
                  value=100),
      sliderInput("mininv",
                  "Minimum % of Diet as Inverts:",
                  min=0,
                  max=100,
                  value=25),
      sliderInput("maxinv",
                  "Maximum % of Diet as Inverts:",
                  min=0,
                  max=100,
                  value=75)
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
    soildata <- rlnorm(n=input$Numsoil,meanlog=input$mean,sdlog=input$sd)  # create distribution of soil concentrations (fake but log normal for the moment)
    amountexposed <- input$time/100#seq(0.1,0.9,0.1)             # percent of diet based on contaminated soil
    massdist <- seq(input$minmass,input$maxmass,length.out=10)                    # birds of a range of masses
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
      ratiobug2[i] <- sample(seq(input$mininv,input$maxinv,length.out=10),1)/100        # sample from a proportion of diet as inverts
      ratioins2[i] <- ratiobug2[i]*(1-0.1)              # remove the soil ingestion rate from the percent of diet as bugs/seeds.  SIR is 0.1 here
      ratioseed2[i] <- (1-ratiobug2[i])*(1-0.1)         # from remainder of (1-SIR)-bug% calculate percent as seeds
      massselect[i] <- sample(massdist,1)               # pick a random bird's mass
      firselect[i] <- 0.301*massselect[i]^0.751         # calculate food ingestion rate
      dirt[i] <- sample(soildata,1)                     # pick a random concentration
      # below= picking random bsaf values and percent on contaminated soil and then combine all values to estimate uptake
      birdsoilresample[i] <- (((sample(seedbsaf,1)*dirt[i]*firselect[i]*ratioseed2[i])+(sample(insbsaf,1)*dirt[i]*firselect[i]*ratioins2[i])+(dirt[i]*firselect[i]*0.1))*sample(amountexposed,1))/massselect[i]
    }
    sdat <- summary(log10(birdsoilresample))
    summStr <- paste(names(sdat), format(sdat, digits = 2), collapse ="\n")
    sdatsoi <- summary(log10(soildata))
    summStrsoi <- paste(names(sdatsoi), format(sdatsoi, digits = 2), collapse ="\n")
    logbird <- log10(birdsoilresample)
    logdirt <- log10(dirt)
    lm1 <- lm(logbird~logdirt)
    sumlm1 <- summary(lm1)
    newx <- seq(min(logdirt),max(logdirt), length.out=100)
    conf_interval <- predict(lm1, newdata=data.frame(logdirt=newx), interval="confidence",level=0.95)
    
    
    par(mfrow=c(2,2),mai=c(0.7,0.7,0.1,0.1))
    plot(density(log10(birdsoilresample)),ylim=c(-0.2,1.2),main="")       # plot log10 distribution of uptake values
    segments(median(log10(birdsoilresample)),0,median(log10(birdsoilresample)),0.2,lty=2)
    text(median(log10(birdsoilresample)),-0.1,bquote(paste("Median=",.(round(median(log10(birdsoilresample)),3)))))
    text(min(log10(birdsoilresample)),0.9,summStr,cex=0.75)
    text(max(log10(birdsoilresample)),0.9,"log10 Uptake\nas mg/kg/day",cex=1)
    plot(density(log10(soildata)),ylim=c(-0.2,1.2),main="")
    segments(median(log10(soildata)),0,median(log10(soildata)),0.2,lty=2)
    text(median(log10(soildata)),-0.1,bquote(paste("Median=",.(round(median(log10(soildata)),3)))))
    text(min(log10(soildata)),0.9,summStrsoi,cex=0.75)
    text(max(log10(soildata)),0.9,"log10 Soil\nConcentrations\nSampled",cex=1)
    plot(logbird~logdirt,main="",ylab="log10 Uptake",xlab="log10 Soil Concentration", pch=16, cex=0.75)
    abline(lm1, col="lightblue", lwd=3)
    lines(newx, conf_interval[,2], col="blue", lty=2)
    lines(newx, conf_interval[,3], col="blue", lty=2)
    mtext(bquote(paste("r^2=",.(round(sumlm1$r.squared,2))," Slope=",.(round(sumlm1$coefficients[2,1],2))," Int=",.(round(sumlm1$coefficients[1,1],2)))), side=1, line=2.2,cex=0.75)
    plot(logbird~ratiobug2,pch=16,xlab="Proportion of diet as invertebrates",ylab="log10 Uptake")
    
  },height=700,width=750)
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(ggdark)

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Negative Binomial distribution"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        "n",
        "Number of items:",
        min = 1,
        max = 100,
        value = 20
      ),
      sliderInput(
        "size",
        "No of succefully trial:",
        min = 1,
        max = 50,
        value = 5
      ),
      sliderInput(
        "mu",
        "Probability",
        min = 0,
        max = 1,
        value = 0.5
      )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(plotOutput("distPlot"))
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$distPlot <- renderPlot({
    negbinom.dist <-
      data.frame(y = dnbinom(1:input$n, size = input$size, prob = input$mu),
                 x= seq(1:input$n))
      ggplot(negbinom.dist) + geom_col(aes(y=y, x = x),
                                           col = 'white',
                                           fill = 'white') +
        ylab("Probabilities") + xlab("No of Trials") + dark_mode() 
      
  })
  
  
}

# Run the application
shinyApp(ui = ui, server = server)

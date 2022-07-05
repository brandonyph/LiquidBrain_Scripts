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
    titlePanel("Poisson distribution"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("n",
                        "Number of items:",
                        min = 1,
                        max = 100,
                        value = 50),
          sliderInput("lambda",
                      "Success Rate:",
                      min = 1,
                      max = 50,
                      value = 30)
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
        pois.dist <- dpois(1:input$n, lambda=input$lambda)
        barplot(pois.dist, col = 'red', border = 'white')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

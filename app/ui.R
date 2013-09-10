# shiny UI

library(shiny)

# Define UI for random distribution application 
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Transcriptome Assembly Comparison"),
  
  # Sidebar with controls 
  sidebarPanel(
    # Choose an analysis
    selectInput("analysis", "Choose an analysis to display:", 
                 choices = list("macro"=1, "micro"=2, "gene"=3, "reads"=4), 
                 selected = 1),
    br(),
    br(),
    # Simple integer interval
    sliderInput("cutoff", "Target coverage cutoff", 
                min=50, max=100, value=80, step=1),
    br(),
    br(),
    numericInput("ntop", "Only show top x assemblies\n\n(0 shows all):", 10),
    br(),
    br(),
    submitButton("Update View")
  ),
  
  # Show a tabset to switch between different levels of metrics
  mainPanel(
    conditionalPanel(condition="input.analysis==1",
      tabsetPanel(
        tabPanel("Reference coverage", plotOutput("reference", width = "100%", height = "500px")), 
        tabPanel("Contig count", plotOutput("contig", width = "100%", height = "500px")), 
        tabPanel("Assembly Coverage", plotOutput("assembly", width = "100%", height = "500px"))
    )),
    conditionalPanel(condition="input.analysis==2",
      tabsetPanel(
        tabPanel("Length", plotOutput("length", width = "100%", height = "500px")), 
        tabPanel("Expression", plotOutput("expression", width = "100%", height = "500px")), 
        tabPanel("GC%", plotOutput("gc", width = "100%", height = "500px")),
        tabPanel("Shannon entropy", plotOutput("shannon", width = "100%", height = "500px"))
    )),
    conditionalPanel(condition="input.analysis==3",
      tabsetPanel(
        tabPanel("Gene list", tableOutput("genelist")), 
        tabPanel("Alignment", htmlOutput("alignment"))
    ))
  )
))

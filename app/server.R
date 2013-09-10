# shiny server

library(shiny)
library(ggplot2)
library(grid)
library(reshape2)
library(directlabels)

# Load macro data
counts <- read.csv("./totalcounts.tsv", sep="\t", head=FALSE)
globalrefcov <- read.csv("./globalreferencecoverage.tsv", sep="\t", head=FALSE)
globalasscov <- read.csv("./globalassemblycoverage.tsv", sep="\t", head=FALSE)
# Load the compiled assembly analytics (micro data)
metrics <- read.csv("allsamplesdata.csv", as.is=TRUE)
numtests <- abs((dim(metrics)[2] - 5))
assemblycov <- metrics[,c(6:(6+numtests-1))]

binsformetric <- function(metric) {
  # bin into deciles (10-percentiles)
  bins <- quantile(metric, seq(0, 1.0, by=0.1))
  return(bins)
}

gplotcovcountpc <- function(metricarray, title='', cutoff=80) {
  x <- seq(10, 100, by=10)
  mbins <- binsformetric(metricarray)
  mhistfn <- function(mlist) {
    h <- hist(metricarray[mlist >= cutoff], breaks=mbins, plot=FALSE)
    return(append(h$counts, 0, after=0))
  }
  mhists <- as.data.frame(sapply(assemblycov, mhistfn))
  mhists <- mhists[-1,]
  mhists$percentile <- x
  mhists <- melt(mhists, id.vars=c('percentile'), variable.name='assembly', value.name='count')
  nrows <- ceiling(numtests / 2)
  p <- ggplot(mhists) +
    geom_line(aes(x=percentile, y=count, group=assembly, colour=assembly)) +
    geom_dl(aes(label=assembly, x=percentile, y=count, colour=assembly), method="top.bumptwice") +
    scale_colour_hue() +
    scale_x_continuous("Percentile", breaks=x) +
    scale_y_continuous(paste("No. reconstructed >= ", cutoff, "%")) +
    ggtitle(title) +
    theme(plot.margin = unit(c(3.8, 1, 1, 1), "cm"), 
          axis.title.x = element_text(hjust=0.5, vjust=0),
          axis.title.y = element_text(hjust=0.5, vjust=0),
          title = element_text(vjust=12, face='bold', colour='#555555', size=14),
          legend.direction = "horizontal", 
          legend.position = c(0.5, 1.10)) +
    guides(col = guide_legend(nrow=nrows,
                              keywidth=3,
                              title.theme=element_text(size=15, face="bold", angle=0)))
  return(p)
}

shinyServer(function(input, output) {
  
  # plot bar chart
  gbar <- function(data, title, xlab, ylab) {
    sorteddata <- data[with(data, order(-V2)),]
    sorteddata$V1 <- factor(sorteddata$V1, levels=sorteddata$V1)
    p <- ggplot(data=sorteddata, aes(x=V1, y=V2)) +
      geom_bar(colour='#006363', stat="identity", fill="#009999") +
      scale_fill_brewer("Assembly", type="div") +
      scale_y_continuous(expand = c(0.01, 0.01)) +
      ylab(ylab) +
      xlab(xlab) +
      ggtitle(title) +
      theme(
        axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        axis.title.x = element_text(hjust=0.5, vjust=0),
        axis.title.y = element_text(hjust=0.5, vjust=0),
        panel.grid.major.x = element_blank(),
        title = element_text(vjust=1, face='bold', colour='#555555', size=14),
        panel.background = element_rect(fill="#F5F5F5"),
        panel.margin = unit(c(0,1,1,1), "cm"),
        plot.margin = unit(c(1,1,1,1), "cm"))
    return(p)
  }
  
  ##############################
  ## Macro coverage section   ##
  ##############################
  
  # Contig counts
  renderContigPlot <- reactive({
    title <- "Total contig counts for each assembly"
    p <- gbar(counts, title, "Assembly", "No. contigs")
    print(p)
  })
  output$contig <- renderPlot({ renderContigPlot() })
  
  # Reference coverage
  renderReferencePlot <- reactive({
    title <- "Proportion of reference transcripts\nrepresented in each assembly"
    p <- gbar(globalrefcov, title, "Assembly", "Coverage")
    print(p)
  })
  output$reference <- renderPlot({ renderReferencePlot() })
  
  # Assembly coverage
  renderAssemblyPlot <- reactive({ 
    title <- "Proportion of assembly contigs\nrepresented in reference (global)"
    p <- gbar(globalasscov, title, "Assembly", "Coverage")
    print(p)
  })
  output$assembly <- renderPlot({ renderAssemblyPlot() })
  
  ##############################
  ## Micro coverage section   ##
  ##############################
  
  # Length
  renderLengthPlot <- reactive({
    title <- paste("Count of genes reconstructed to at least ", input$cutoff, "% vs. length percentile")
    p <- gplotcovcountpc(metrics$Length, title=title, cutoff=input$cutoff)
    print(p)
  })
  output$length <- renderPlot({ renderLengthPlot() })
  
  # Expression
  renderExpressionPlot <- reactive({
    title <- "Count of genes reconstructed to at least 80% vs. expression percentile"
    p <- gplotcovcountpc(metrics$FPKM[metrics$FPKM > 0], title=title, cutoff=input$cutoff)
    print(p)
  })
  output$expression <- renderPlot({ renderExpressionPlot() })
  
  # GC %
  renderGCPlot <- reactive({
    metrics$GC. <-  sapply(metrics$GC., function(x) eval(parse(text=x)))
    title <- "Count of genes reconstructed to at least 80% vs. GC content percentile"
      p <- gplotcovcountpc(metrics$GC., title=title, cutoff=input$cutoff)
      print(p)
  })
  output$gc <- renderPlot({ renderGCPlot() })
  
  # Shannon Entropy
  renderShannonPlot <- reactive({
    title <- "Count of genes reconstructed to at least 80% vs. complexity percentile"
    p <- gplotcovcountpc(metrics$Complexity, title=title, cutoff=input$cutoff)
    print(p)
  })
  output$shannon <- renderPlot({ renderShannonPlot() })

  ##############################
  ## Gene coverage section    ##
  ##############################
  
  # Gene list coverage table
  renderGeneTable <- reactive({
    
  })
  output$genelist <- renderTable({ renderGeneTable() })
  
  # Alignment
  renderAlignment <- reactive({
    tags$iframe(src="../align.html", style="width: 95%; height: 90%; border: 0")
  })
  output$alignment <- renderUI({ renderAlignment() })
})
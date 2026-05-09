library(shiny)
library(ggplot2)
library(plotly)
library(viridis)
library(ggrepel)
library(dplyr)

# Sample Data
data <- read.table("dataset.txt", header = FALSE, 
                   sep = "\t", stringsAsFactors = FALSE)
colnames(data) <- c("SNP_ID", "ID", "GENE_ID", 
                    "PROBE_ID", "CHR_SNP", "CHR_GENE", 
                    "SNPpos", "TSSpos", "Distance", "rvalue", "pvalue", "log10pvalue")
# Remove 'ID' column
data <- subset(data, select = -ID)

data <- data %>%
  mutate(
    SNPpos = as.integer(SNPpos),
    Distance = as.integer(Distance),
    CHR_SNP = as.factor(CHR_SNP),
    CHR_GENE = as.factor(CHR_GENE)
  )

# Define significance threshold for Manhattan plot
threshold <- 100

# Define UI
ui <- fluidPage(
  titlePanel("Transcriptome and genome sequencing uncovers functional variation in humans"),
  strong("Authors: Eira Fontanals Muñoz, Marc Jordi Torres Valero & Òscar Contreras Parejo"), br(),
  em(" Subject: Data Visualitzation"),
  br(), br(),
  
  sliderInput("log10pvalueInput", "Select log10pvalue:", min = min(data$log10pvalue), 
              max = max(data$log10pvalue), value = c(min(data$log10pvalue), max(data$log10pvalue)), step = 1),br(),
  
  
  plotlyOutput("heatmap"), br(),
    p("Since the log10pvalue data covers a wide range and most of the important information 
      is concentrated in smaller values, a slider input allows users to select a specific 
      range. This interactive heatmap enables dynamic filtering of SNP log10pvalue across chromosomes, 
      improving the visualization of significant regions."), br(),
  
  plotlyOutput("manhattan_plot"),br(),
  verbatimTextOutput("snp_info"),br(),
  p("The built-in interactive plotly_click feature in the Manhattan plot allows users to identify 
    the most significant SNPs across all chromosomes simply by clicking on them. With the threshold 
    line added to the plot, it becomes easier for users to spot which SNPs are interesting to investigate."),
  br(),
  selectInput("chromosomeInput", "Select Chromosome:",
              choices = levels(data$CHR_SNP), selected = levels(data$CHR_SNP)[1]),br(),
  plotlyOutput("scatter_plot"), br(),
  p("By analyzing the previous plots, we can identify specific chromosomes to explore further. The select 
    input interaction allows users to plot the significance of the distance only for the chosen chromosome. 
    This enhances the visualization of the scatter plot, making it more focused and informative."),
  br(),
  
  
)

# Define Server
server <- function(input, output, session) {
  data$log10p_bin <- cut(data$log10pvalue, breaks = seq(0, max(data$log10pvalue), by = 1), include.lowest = TRUE)
  
  # Render Heatmap with log10pvalue Filter
  output$heatmap <- renderPlotly({
    filtered_data <- data[data$log10pvalue >= input$log10pvalueInput[1] & 
                            data$log10pvalue <= input$log10pvalueInput[2], ]
    
    p <- ggplot(filtered_data, aes(x = CHR_SNP, y = log10pvalue)) +
      geom_bin2d(bins = 50, aes(fill = after_stat(count / sum(count)))) +
      scale_fill_viridis_c(name = "Density") +
      theme_minimal() +
      labs(title = "Heatmap of log10pvalue Distributions Across Chromosomes",
           x = "Chromosome", y = "log10pvalue") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    
    ggplotly(p, tooltip = "fill") %>% config(displayModeBar = FALSE)
  })
  
  # Render Manhattan Plot
  output$manhattan_plot <- renderPlotly({
    p <- ggplot(data, aes(x = CHR_SNP, y = log10pvalue)) + 
      geom_point(aes(color = log10pvalue), alpha = 0.8, size = 1.5) +
      scale_color_gradient(low = "#355573", high = "#bed93b")+
      labs(title = "Manhattan Plot of log10(p-values) Across Chromosomes",
           x = "Chromosome", y = "log10(p-value)") +  
      theme_minimal() +
      geom_hline(yintercept = threshold, color = "black", linetype = "dashed", linewidth = 0.2) +
      geom_label_repel(data = data[data$log10pvalue > threshold, ], 
                       aes(label = SNP_ID, alpha=0.7), 
                       size = 3, color = "black", 
                       force = 1.7, 
                       box.padding = 0.5) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  
            axis.title = element_text(size = 14),
            legend.position = "none",  
            panel.grid.major = element_blank(),  
            panel.grid.minor = element_blank())
    
    ggplotly(p, tooltip = "text") %>% config(displayModeBar = FALSE) %>%
      event_register("plotly_click")
  })
  
  # Render Scatter Plot for Selected Chromosome
  output$scatter_plot <- renderPlotly({
    # Filter data dynamically for the chromosomeInput
    filtered_data <- data[data$CHR_SNP ==input$chromosomeInput, ]
    
    p <- ggplot(filtered_data, aes(x = Distance, y = log10pvalue, color = log10pvalue)) +
      geom_point(alpha = 0.5) +
      scale_color_viridis_c() +  
      geom_label_repel(aes(label = SNP_ID), 
                       size = 3, color = "black", 
                       force = 1.5, 
                       box.padding = 0.5) + 
      theme_minimal() +
      labs(title = paste("Scatter Plot for Chromosome", input$chromosomeInput),
           x = "Distance", y = "log10pvalue", color = "log10pvalue") +
      theme(legend.position = "right",
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    
    ggplotly(p, tooltip = "text")  %>% config(displayModeBar = FALSE)
  })
  
  # Capture Click Event and Return SNP Name
  output$snp_info <- renderPrint({
    click_data <- event_data("plotly_click")  # Capture click event
    
    if (is.null(click_data)) {
      return("Click on a point to see SNP name")
    } else {
      selected_snp <- data[data$CHR_SNP == click_data$x & data$log10pvalue == click_data$y, ]$SNP_ID
      paste0("SNP: ", selected_snp)
    }
  })
}


shinyApp(ui, server)

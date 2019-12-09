########
##Shiny app to create 
library(shiny)
library(DT)
library(tidyverse)
library(phyloseq)
library(plotly)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(viridis)
library(vegan)
options(shiny.maxRequestSize = 30*1024^2)



# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("Taxonomy Plots"),
    
    
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel("Upload Files",
                     
            fileInput("tax_file", "Choose *.taxonomy File",
                      multiple = TRUE,
                      accept = c(".taxonomy")),
            fileInput("shared_file", "Choose *.shared File",
                      multiple = TRUE,
                      accept = c(".shared")),
            fileInput("meta_file", "Choose Metadata File (*.csv)",
                      multiple = TRUE,
                      accept = c(".csv")),
            checkboxInput("rarefaction", "Calculate rarefaction curve", FALSE),
            ##since the taxonomy table take a lot of time to load, I'm giving an option
            #to make a smaller table just to see the top 50 otus
            selectInput("tax_plot_num", "Choose number of OTUs to include in taxonomy table",
                        choices = c("50", "100", "200")),
            uiOutput("order_list")
            
            # downloadButton('downloadData', 'Download Data'),
            # downloadButton('downloadPlot', 'Download Plot')
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("Taxonomy Table", 
                         dataTableOutput("tax_table")),
                tabPanel("Stacked Barplot", 
                         plotOutput("taxonomy_plot",width = "1000px", height = "850px")),
                tabPanel("Order-level Abundance", 
                         plotOutput("order_plot", width = "1000px", height = "850px")),
                tabPanel("Rarefaction Curve", plotOutput("rarefaction_plot", width = "500px", height = "500px"))
                
    )
     )
    )
)


server <- function(input, output) {
    ###making dataframe out of files: .taxonomy and .shared files
    abundance_data <- reactive({ 
        
        sharedfile <- input$shared_file
        taxfile <- input$tax_file
        meta.file <- input$meta_file
        
        #read in metadata file
        metafile <- read.csv(meta.file$datapath)
        #assign rownames to be sample.id names
        rownames(metafile) <- metafile$sample.id
        map_phy <- sample_data(metafile)
        
        #import mothur files, make a phyloseq object
       phyloseq_import <- import_mothur(mothur_shared_file = sharedfile$datapath, mothur_constaxonomy_file = taxfile$datapath)

        #added metadata to phyloseq object
        moth_merge <- merge_phyloseq(phyloseq_import,map_phy)
        
        #rename columns of the tax table
        colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
        
        moth_merge

        
    })
    
    output$rarefaction_plot = renderPlot({
## make rarefaction curve optional to make loadiing faster
      if (input$rarefaction == FALSE) {
        return(NULL)
      } else {
        
      withProgress(message = "Progress bar"
                   , detail ='Rendering plot...'
                   , value = 0,{
                     
                     # pausa
                     Sys.sleep(1)               
                     # update progress bar
                    
                     col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
                      lty <- c("solid", "dashed", "longdash", "dotdash")
                      pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
                     
                      setProgress(0.5)
                     
                     
                      rarefaction <- with(pars[1:3,],
                            rarecurve(t(otu_table(abundance_data())), col = col, lty = lty, step=50, cex=0.5))
      
                     # update progress to complete status
                     setProgress(1)
                   })
        
        }
      })   
    
    #taxonomy table listing top 50 taxa (to make loading faster...)
    output$tax_table = renderDataTable({
      data <- tax_glom(abundance_data(), taxrank = "Genus")
      
        tax_table <- tax_table(data)
        otu_table <- otu_table(data)
        
        new_table <- cbind(otu_table, tax_table)
        
        new_table[1:input$tax_plot_num,]
        
    }
        
    )  
        
   ###stacked bar chart of Phyla
    output$taxonomy_plot = renderPlot({
       
          moth_phylum <- abundance_data() %>%
            tax_glom(taxrank = "Phylum") %>% #agglomerate at phylum level
            transform_sample_counts(function(x) {100 * x/sum(x)}) %>% #transform to rel abund
            psmelt() %>%
            #filter(Abundance > 1.0) %>%
            arrange(desc(Phylum))
                      
            moth_phylum$Phylum[moth_phylum$Abundance < 1.0] <- "< 1% abund." 
                      
                      moth_phylum
                      
    
                     
          ggplot(moth_phylum, aes(x = new.name, y = Abundance, fill = Phylum)) +
            geom_bar(stat = "identity") +
            scale_fill_viridis(discrete = TRUE, option="D") +
            theme_minimal() +
            geom_text(aes(label = Phylum), color = "white", size = 5, hjust = 0.5, vjust = 3, position =     "stack") +
            ggtitle(label = "Relative Abundance of Taxa by Phylum") +
            theme(axis.title.x = element_blank(), legend.position = "none")
      
      }
        )
    
   
    
    output$order_list = renderUI({
      taxa <- tax_table(abundance_data())
      orders <- taxa[,4]
      
      selectInput(inputId = "order",
                  label = "Choose Order",
                  choices = orders,
                  selected = "")
      
      
    })
    
    
     ####THIS ISN'T WORKING YET.... Goal: to make a plot of order abundance in each sample with a drop
    ####down menu to select desired order (working... just the plot isn't working)
    
    output$order_plot = renderPlot({
      # ps.rarefied = rarefy_even_depth(abundance_data(), rngseed=1, sample.size=0.9*min(sample_sums(abundance_data())), replace=F)
      # 
      # moth_order <- ps.rarefied %>%
      #   tax_glom(taxrank = "Order", NArm =FALSE) %>% #agglomerate at phylum level
      #   transform_sample_counts(function(x) {100*x/sum(x)}) %>% #transform to rel abund
      #   psmelt() %>%
      #   filter(Abundance >= 1.0) %>%
      #   filter(Order == "input$select_order")
      x <- tax_table()
      
      as.data.frame(x)
      
      tax <- x %>%
        filter(Order == "input$order_list")
      
      
      ggplot(tax, aes(x = new.name, y = Abundance)) +
        geom_bar(stat = "identity") +
        scale_fill_viridis(discrete=TRUE, option="C") +
        theme_minimal() +
        ggtitle(label = "Relative Abundance") 
      
    })
    
    ####THIS ALSO ISN'T WORKING YET, want to add functions to download the plots and taxonomy tables
    
    # output$downloadData <- downloadHandler(
    #   filename = function() { paste(input$table_plot, '.csv', sep='') },
    #   content = function(file) {
    #     write.csv(datatasetInput(), file)
    #   }
    # )
    # output$downloadPlot <- downloadHandler(
    #   filename = function() { paste(input$table_plot, '.png', sep='') },
    #   content = function(file) {
    #     ggsave(file,plotInput())
    #   }
    # )
    # 
}

# Run the application 
shinyApp(ui = ui, server = server)
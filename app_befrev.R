# Load required libraries
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(readr)
library(psych) 
library(DT)
library(purrr)
library(forcats)
library(shinycssloaders)

# Define UI
library(shiny)
library(shinythemes)

ui <- fluidPage(
  theme = shinythemes::shinytheme("united"),
  
  # Intro and Logo Section
  fluidRow(
    column(2,
           tags$img(src = "logo.png", height = "550px", width = "300px", alt = "GeSAT Logo")
    ),
    column(10,
           div(style = "padding: 30px; background-color: #f9f9f9; border-radius: 12px; box-shadow: 0 4px 12px rgba(0,0,0,0.1);",
               tags$h2("Welcome to GeSAT", style = "color: #A020F0; font-weight: bold;"),
               tags$p("Your ultimate tool for identifying stable reference genes for qPCR data normalization",
                      style = "font-size: 18px; color: #555;"),
               
               tags$h3("Why Choose GeSAT?", style = "margin-top: 30px; color: #A020F0;"),
               tags$p("GeSAT empowers researchers to identify the most stable reference genes for qPCR data normalization using advanced, peer-reviewed methods. Upload your gene expression data and get a comprehensive stability ranking in minutes.",
                      style = "font-size: 16px; color: #666;"),
               
               tags$h4("✅ Comprehensive Analysis", style = "margin-top: 20px; color: #3e4741;"),
               tags$p("Utilizes multiple established methods including Delta Ct, NormFinder, geNorm, BestKeeper, and Mixed Model for robust results.",
                      style = "font-size: 15px; color: #555;"),
               
               tags$h4("✅ User-Friendly Interface", style = "margin-top: 20px; color: #3e4741;"),
               tags$p("Simply upload your data and explore intuitive visualizations of gene stability rankings.",
                      style = "font-size: 15px; color: #555;"),
               
               tags$h4("✅ Accurate and Reliable", style = "margin-top: 20px; color: #3e4741;"),
               tags$p("Combines results from multiple algorithms to provide a cumulative ranking, ensuring dependable reference gene selection.",
                      style = "font-size: 15px; color: #555;")
           )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("std_curve_file", "Upload Standard Curve File (CSV)", accept = ".csv"),
      tags$hr(),
      tags$h5("Download Example Standard Dilution Ct value Files"),
      tags$a(href = "std_curve_file.csv", "Download File", download = NA, target = "_blank", class = "btn btn-info"),
      fileInput("ct_file", "Upload Ct Values File (CSV)", accept = ".csv"),
      tags$hr(),
      tags$h5("Download Example Ct value Files"),
      tags$a(href = "ct_file.csv", "Download File", download = NA, target = "_blank", class = "btn btn-info"),
      selectInput("genes", "Select Genes", choices = NULL, multiple = TRUE),
      helpText("Upload a wide-format CSV file with columns: Group, Sample, and gene names as other columns.")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Standard Curves", 
                 plotOutput("std_curve_plot", height = "700px"),
                 downloadButton("download_std_curve_plot", "Download Standard Curve Plot", class = "btn-success"),
                 downloadButton("download_efficiency", "Download PCR Efficiency CSV", class = "btn-success")
        ),
        tabPanel("Ct Distribution", 
                 plotOutput("ct_dist_plot", height = "700px"),
                 downloadButton("download_ct_dist_plot", "Download Ct Distribution Plot", class = "btn-success")
        ),
        tabPanel("BestKeeper", 
                 plotOutput("bestkeeper_plot", height = "500px"),
                 downloadButton("download_bestkeeper_plot", "Download BestKeeper Plot", class = "btn-success"),
                 tableOutput("bestkeeper_stats_table"),
                 downloadButton("download_bestkeeper_stats", "Download Statistics Table", class = "btn-success"),
                 tableOutput("bestkeeper_pairwise_corr_table"),
                 downloadButton("download_bestkeeper_corr", "Download Pairwise Correlation Table", class = "btn-success"),
                 tableOutput("bestkeeper_gene_vs_bki_table"),
                 downloadButton("download_bestkeeper_vs_bki", "Download Gene vs BKI Table", class = "btn-success"),
                 tableOutput("ranking_table"),
                 downloadButton("download_bestkeeper_ranking", "Download Detailed Ranking Table", class = "btn-success")
        ),
        tabPanel("Mixed Model Stability",
                 h4("Best Reference Gene Combination (Based on ICC Lower Bound and Non-significant LRT)"),
                 textOutput("bestMixedModelCombination"),
                 br(),
                 h4("Ranked Gene Combinations (Filtered by p > 0.05, Sorted by ICC Lower Bound)"),
                 tableOutput("mixedModelTable"),
                 downloadButton("download_ranked_mixed_model", "Download Ranked Results", class = "btn-success"),
                 br(), br(),
                 h4("Visualization of Gene Combination Ranking"),
                 plotOutput("iccRankPlot"),
                 downloadButton("download_icc_rank_plot", "Download Plot", class = "btn-success")
        ),
        tabPanel("GeNorm", 
                 # --- CSS for strong flashing effect ---
                 tags$style(HTML("
    #use_ctrlgene.blink {
      animation: flash 0.8s infinite;
    }
    @keyframes flash {
      0%, 100% { background-color: #ffeb3b; box-shadow: 0 0 10px #ffeb3b; }
      50% { background-color: #ff0000; box-shadow: 0 0 15px #ff0000; }
    }
  ")),
                 
                 # --- The checkbox ---
                 checkboxInput("use_ctrlgene", "Use ctrlGene::geNorm2()", value = FALSE),
                 
                 # --- JS: blink until checked ---
                 tags$script(HTML("
    $(document).ready(function() {
      $('#use_ctrlgene').addClass('blink');
      $('#use_ctrlgene').on('change', function() {
        if (this.checked) {
          $(this).removeClass('blink');
        } else {
          $(this).addClass('blink');
        }
      });
    });
  ")),
                 
                 # --- Your existing tabsetPanel ---
                 tabsetPanel(
                   tabPanel("Stability Ranking",
                            tags$div(class = "result-container",
                                     plotOutput("genorm_plot", height = "500px"),
                                     downloadButton("download_genorm_plot", "Download Plot", class = "btn-success"),
                                     tags$br(), tags$br(),
                                     tableOutput("genorm_ranking_table"),
                                     downloadButton("download_genorm_ranking_table", "Download Table", class = "btn-success")
                            )
                   ),
                   
                   tabPanel("Pairwise Variation",
                            tags$div(class = "result-container",
                                     plotOutput("genorm_pairwise_plot", height = "500px"),
                                     tags$div(class = "alert alert-info",
                                              "The red dashed line indicates the recommended cutoff (V < 0.15)",
                                              tags$br(),
                                              "Optimal number of reference genes is below this threshold"
                                     ),
                                     downloadButton("download_genorm_pairwise_plot", "Download Plot", class = "btn-success")
                            )
                   )
                 )
        )
        
        ,
        tabPanel("Delta Ct", 
                 fluidRow(
                   column(12, 
                          h4("ΔCt Plot (Pairwise ΔCt between Reference Genes)"), 
                          plotOutput("deltaCtPlot", height = "700px"),
                          downloadButton("download_delta_ct_plot", "Download ΔCt Plot", class = "btn-success")
                   )
                 ),
                 fluidRow(
                   column(12, 
                          h4("Gene Ranking Based on Mean of ΔCt Standard Deviations"),
                          p("Lower values of Mean_of_StdDev indicate more stable reference genes."),
                          DT::dataTableOutput("summaryTable"),
                          downloadButton("download_delta_ct_summary", "Download Summary Table", class = "btn-success")
                   )
                 )
        ),
        tabPanel("NormFinder",
                 h4("Stability Ranking (Ordered)"),
                 DT::dataTableOutput("norm_ordered"),
                 downloadButton("download_norm_ordered", "Download Ordered Table", class = "btn-success"),
                 
                 h4("Unordered Stability Scores"),
                 DT::dataTableOutput("norm_unordered"),
                 downloadButton("download_norm_unordered", "Download Unordered Table", class = "btn-success"),
                 
                 h4("Top Stable Gene Pairs"),
                 DT::dataTableOutput("norm_pairs"),
                 downloadButton("download_norm_pairs", "Download Gene Pairs Table", class = "btn-success"),
                 
                 h4("Ranked Genes by Stability and Group SD"),
                 DT::dataTableOutput("norm_ranked_table"),
                 downloadButton("download_norm_ranked_table", "Download Ranked Table", class = "btn-success")
        ),
        tabPanel("Cumulative Ranking",
                 fluidRow(
                   column(
                     width = 6,
                     h4("Download Table"),
                     downloadButton("download_cumulative_table", "Detailed Table", class = "btn-success")
                   ),
                   column(
                     width = 6,
                     h4("Download Plot"),
                     downloadButton("download_cumulative_plot", "Download Plot", class = "btn-success")
                   )
                 ),
                 br(),
                 fluidRow(
                   column(
                     width = 12,
                     h4("Optional: Raw Ranked Data"),
                     downloadButton("download_cumulative_raw", "Raw Ranked Data", class = "btn-success")
                   )
                 ),
                 br(), br(),
                 fluidRow(
                   column(
                     width = 12,
                     h4("Preview: Stability Ranking Table"),
                     DTOutput("cumulative_table")
                   )
                 ),
                 br(),
                 fluidRow(
                   column(
                     width = 12,
                     h4("Preview: Stability Plot"),
                     plotOutput("cumulative_plot", height = "600px")
                   )
                 )
        )
      )
    )
  ),
  
  tags$hr(),
  tags$footer(
    style = "text-align:center; color: #888; padding: 10px; font-size: 90%;",
    HTML(
      "GeSAT | Developed by Dr.Ankita Gurao at the Buffalo Genomics Lab, led by Dr. Ranjit Singh Kataria & Dr. Mahesh Shivanand Dige, ICAR-National Bureau of Animal Genetic Resources | 2025<br><br>
    
    <b>Contact:</b><br>
    Email: <a href='mailto:ranjit.kataria@icar.org.in'>ranjit.kataria@icar.org.in</a><br>
    <b>References:</b><br>
    <b>NormFinder:</b> Andersen CL, Jensen JL, Ørntoft TF. Cancer Res. 2004;64(15):5245-5250.<br>
    <b>Delta Ct:</b> Silver N, Best S, Jiang J, Thein SL. BMC Mol Biol. 2006;7:33.<br>
    <b>BestKeeper:</b> Pfaffl MW, Tichopad A, Prgomet C, Neuvians TP. Biotechnol Lett. 2004;26(6):509-515.<br>
    <b>GeNorm:</b> Vandesompele J, De Preter K, Pattyn F, et al. Genome Biol. 2002;3(7):research0034.<br>
    <b>Mixed Model Stability:</b> Dai H, Charnigo R, Vyhlidal CA, et al. Stat Med. 2013;32(18):3115-25.<br>
    <b>ctrlGene:</b> Assess the Stability of Candidate Housekeeping Genes. (CRAN R package). Retrieved from https://github.com/cran/ctrlGene.
      "
    )
  )
)
# Define Server
server <- function(input, output, session) {
  
  # Load and process standard curve data
  std_curve_data <- reactive({
    req(input$std_curve_file)
    read_csv(input$std_curve_file$datapath)
  })
  
  # Load and process Ct data
  ct_data <- reactive({
    req(input$ct_file)
    read_csv(input$ct_file$datapath)
  })

#Bestkeeper analysis  
  bestKeeperAnalysis = function(ct_data, ctVal = TRUE) {
    # Remove Group and Sample columns (assumed to be first two columns)
    ct_data = ct_data[, -c(1, 2)]
    
    # Ensure all columns are numeric (handles character/factor columns too)
    ct_data <- data.frame(lapply(ct_data, function(x) as.numeric(as.character(x))))
    
    if (!ctVal) {
      ct_data = log2(ct_data)
    }
    
    # ---- CT Statistics (formerly pSta) ----
    N = rep(nrow(ct_data), times = ncol(ct_data))
    GM_CT = apply(ct_data, 2, psych::geometric.mean)
    AM_CT = apply(ct_data, 2, mean)
    Min_CT = apply(ct_data, 2, min)
    Max_CT = apply(ct_data, 2, max)
    
    AVEDEV = function(x) {
      sum(abs(x - mean(x))) / length(x)
    }
    
    SD_CT = apply(ct_data, 2, AVEDEV)
    CV_CT = 100 * SD_CT / AM_CT
    Min_x_fold = -1 / 2^(Min_CT - GM_CT)
    Max_x_fold = 2^(Max_CT - GM_CT)
    SD_x_fold = 2^SD_CT
    
    CT_statistics = matrix(
      c(N, GM_CT, AM_CT, Min_CT, Max_CT, SD_CT, CV_CT, Min_x_fold, Max_x_fold, SD_x_fold),
      ncol = length(GM_CT), nrow = 10, byrow = TRUE
    )
    
    # Safely assign row and column names
    if (ncol(CT_statistics) == length(GM_CT)) {
      colnames(CT_statistics) = colnames(ct_data)
      rownames(CT_statistics) = c(
        "N", "GM[CT]", "AM[CT]", "Min[CT]", "Max[CT]", "SD[+/- CT]",
        "CV[%CT]", "Min[x-fold]", "Max[x-fold]", "SD[+/- x-fold]"
      )
    } else {
      stop("Mismatch between number of rows/columns and the assigned names.")
    }
    
    CT_statistics = round(CT_statistics, 3)
    
    # ---- Pairwise Pearson Correlation ----
    n = ncol(ct_data)
    pairwise_corr = matrix(ncol = n, nrow = 2 * (n - 1))
    colnames(pairwise_corr) = colnames(ct_data)
    rowname_rz = list()
    
    for (i in 1:n) {
      for (j in 2:n) {
        # Ensure both columns are numeric
        x = as.numeric(ct_data[, i])
        y = as.numeric(ct_data[, j])
        
        # Check for any NA values after conversion
        if(any(is.na(x)) | any(is.na(y))) {
          next
        }
        
        COT_TEST = cor.test(x, y)
        pairwise_corr[2 * (j - 1) - 1, i] = COT_TEST$estimate
        pairwise_corr[2 * (j - 1), i] = COT_TEST$p.value
      }
      if (i + 1 > n) break
      rowname_rz = c(rowname_rz, colnames(ct_data[i + 1]), 'p-value')
    }
    
    rownames(pairwise_corr) = unlist(rowname_rz)
    pairwise_corr = round(pairwise_corr, 3)
    
    # ---- Gene vs BestKeeper Index ----
    BKI = apply(ct_data, 1, psych::geometric.mean)
    
    # Ensure BKI is numeric
    BKI = as.numeric(BKI)
    
    gene_vs_BKI = matrix()
    
    r_values = numeric(n)
    p_values = numeric(n)
    
    for (i in 1:n) {
      # Ensure both ct_data and BKI are numeric
      x = as.numeric(ct_data[, i])
      
      # Check for any NA values after conversion
      if(any(is.na(x)) | any(is.na(BKI))) {
        next
      }
      
      COT_TEST = cor.test(x, BKI)
      r = COT_TEST$estimate
      p = COT_TEST$p.value
      r_values[i] = r
      p_values[i] = p
      
      r2 = r^2
      count = floor((length(x) + length(BKI)) / 2)
      
      SC_xy = sum(x * BKI) - sum(x) * sum(BKI) / count
      SS_x = sum(BKI^2) - sum(BKI)^2 / length(BKI)
      SS_y = sum(x^2) - sum(x)^2 / length(x)
      slope = SC_xy / SS_x
      intercept = mean(x) - slope * mean(BKI)
      SE = sqrt(SS_y * (1 - r2) / (count - 2))
      Power_x_fold = 2^slope
      
      buf = matrix(c(r, r2, intercept, slope, SE, p, Power_x_fold))
      
      if (length(gene_vs_BKI) == 1) {
        gene_vs_BKI = buf
      } else {
        gene_vs_BKI = cbind(gene_vs_BKI, buf)
      }
    }
    
    colnames(gene_vs_BKI) = colnames(ct_data)
    rownames(gene_vs_BKI) = c(
      "coeff. of corr. [r]", "coeff. of det. [r^2]", "intercept [CT]",
      "slope [CT]", "SE [CT]", "p-value", "Power [x-fold]"
    )
    gene_vs_BKI = round(gene_vs_BKI, 3)
    
    # ---- Ranking Logic ----
    # Rank based on SD_CT (least SD_CT gets highest rank)
    SD_CT_rank = rank(SD_CT, ties.method = "min") # Smallest SD gets rank 1
    
    # Rank based on r values with significant p-values
    sig_idx <- which(p_values <= 0.05)
    r_rank <- rep(0, length(r_values))
    if(length(sig_idx) > 0) {
      r_rank[sig_idx] <- rank(-r_values[sig_idx], ties.method = "min")
    }
    
    # Combine ranks
    ranks = data.frame(
      Gene = colnames(ct_data),
      SD_CT_Rank = SD_CT_rank,
      Correlation_Rank = r_rank
    )
    
    # Print ranks to check
    print("Ranks:")
    print(ranks)
    
    # ---- Final Output ----
    result = list(
      CT.statistics = CT_statistics,
      Pairwise.Correlation = pairwise_corr,
      Gene.vs.BKI = gene_vs_BKI,
      Gene.Ranks = ranks
    )
    
    return(result)
  }
  
  # Compute BestKeeper results
  bestKeeper_results <- reactive({
    req(ct_data())  # Ensure ct_data is available
    bestKeeperAnalysis(ct_data())
  })
  
  # --- BestKeeper Results Tables ---
  
  output$bestkeeper_stats_table <- renderTable({
    req(bestKeeper_results())
    as.data.frame(bestKeeper_results()$CT.statistics) %>% 
      tibble::rownames_to_column("Descriptive Statistics")
  })
  
  output$bestkeeper_pairwise_corr_table <- renderTable({
    req(bestKeeper_results())
    as.data.frame(bestKeeper_results()$Pairwise.Correlation) %>% 
      tibble::rownames_to_column("Correlation: pairwise")
  })
  
output$ranking_table <- renderTable({
    req(bestKeeper_results())  
    bestKeeper_results()$Gene.Ranks %>%
      dplyr::rename(
        'SD[+/- CT] based ranking' = SD_CT_Rank,
        'coeff. of corr. [r] based ranking' = Correlation_Rank
      )})
  
 
output$bestkeeper_plot <- renderPlot({
    req(bestKeeper_results())  # Ensure results are available
    gene_vs_bki <- bestKeeper_results()$Gene.vs.BKI
    
    output$ranking_table <- renderTable({
      req(bestKeeper_results())
      
      # Create the ranking data frame
      ranks <- bestKeeper_results()$Gene.Ranks
      
      # Include SD_CT and r_values columns in the ranks table
      ranks_with_values <- cbind(
        ranks,
        SD_CT = round(bestKeeper_results()$CT.statistics["SD[+/- CT]", ], 3),  # Add SD_CT values
        r_values = round(bestKeeper_results()$Gene.vs.BKI["coeff. of corr. [r]", ], 3)  # Add r_values
      )
      

      # Rename the columns appropriately
      ranks_with_values <- ranks_with_values %>%
        dplyr::rename(
          'SD[+/- CT] based ranking' = SD_CT_Rank,
          'coeff. of corr. [r] based ranking' = Correlation_Rank,
          'SD[+/- CT]' = SD_CT,
          'coeff. of corr. [r]' = r_values
        )
      
      # Return the final ranks table
      ranks_with_values
    })
    
    output$bestkeeper_gene_vs_bki_table <- renderTable({
      req(bestKeeper_results())
      
      # Convert Gene vs BKI results to a data frame and add row names as a column
      gene_vs_bki_df <- as.data.frame(bestKeeper_results()$Gene.vs.BKI) %>%
        tibble::rownames_to_column("Correlation: Gene vs BKI")
      
      # Return the table
      gene_vs_bki_df
    })
    # Transpose and convert to data frame
    df <- as.data.frame(t(gene_vs_bki))
    df$Gene <- rownames(df)
    
    # Add significance column
    df$significance <- ifelse(df$`p-value` < 0.05, "Significant", "Not significant")
    
    # Plot
    output$bestkeeper_plot <- renderPlot({
      req(bestKeeper_results())
      
      gene_vs_bki <- bestKeeper_results()$Gene.vs.BKI
      df <- as.data.frame(t(gene_vs_bki))
      df$Gene <- rownames(df)
      df$significance <- ifelse(df$`p-value` < 0.05, "Significant", "Not significant")
      
      ggplot(df, aes(x = Gene, y = `coeff. of corr. [r]`, fill = significance)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = ifelse(significance == "Not significant", "ns", "")),
                  vjust = -0.5, size = 4) +
        scale_fill_manual(values = c("Significant" = "green", "Not significant" = "black")) +
        theme_minimal() +
        labs(x = "Gene", y = "Correlation with BestKeeper", title = "Gene vs BestKeeper Index") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
    
    # Statistics Table
    output$download_bestkeeper_stats <- downloadHandler(
      filename = function() paste0("BestKeeper_Statistics_", Sys.Date(), ".csv"),
      content = function(file) {
        write.csv(as.data.frame(bestKeeper_results()$CT.statistics), file)
      }
    )
    
    # Pairwise Correlation Table
    output$download_bestkeeper_corr <- downloadHandler(
      filename = function() paste0("BestKeeper_Pairwise_Correlation_", Sys.Date(), ".csv"),
      content = function(file) {
        write.csv(as.data.frame(bestKeeper_results()$Pairwise.Correlation), file)
      }
    )
    
    # Gene vs BKI Table
    output$download_bestkeeper_vs_bki <- downloadHandler(
      filename = function() paste0("BestKeeper_Gene_vs_BKI_", Sys.Date(), ".csv"),
      content = function(file) {
        write.csv(as.data.frame(bestKeeper_results()$Gene.vs.BKI), file)
      }
    )
    
    # Ranking Table
    output$download_bestkeeper_ranking <- downloadHandler(
      filename = function() paste0("BestKeeper_Ranking_", Sys.Date(), ".csv"),
      content = function(file) {
        write.csv(as.data.frame(bestKeeper_results()$Gene.Ranks), file)
      }
    )
    
    # Detailed Ranking Table (if different)
    output$download_bestkeeper_ranking <- downloadHandler(
      filename = function() paste0("BestKeeper_Detailed_Ranking_", Sys.Date(), ".csv"),
      content = function(file) {
        write.csv(as.data.frame(bestKeeper_results()$Gene.Ranks), file)
      }
    )   
    output$download_bestkeeper_plot <- downloadHandler(
      filename = function() {
        paste0("BestKeeper_plot_", Sys.Date(), ".png")
      },
      content = function(file) {
        gene_vs_bki <- bestKeeper_results()$Gene.vs.BKI
        df <- as.data.frame(t(gene_vs_bki))
        df$Gene <- rownames(df)
        df$significance <- ifelse(df$`p-value` < 0.05, "Significant", "Not significant")
        
        ggsave(file,
               plot = ggplot(df, aes(x = Gene, y = `coeff. of corr. [r]`, fill = significance)) +
                 geom_bar(stat = "identity") +
                 geom_text(aes(label = ifelse(significance == "Not significant", "ns", "")),
                           vjust = -0.5, size = 4) +
                 scale_fill_manual(values = c("Significant" = "green", "Not significant" = "black")) +
                 theme_minimal() +
                 labs(x = "Gene", y = "Correlation with BestKeeper", title = "Gene vs BestKeeper Index") +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1),
                       width = 10, height = 8, dpi = 300)
        )
      }
    )
})

# -----GENORM------------------
library(ctrlGene)
library(dplyr)

geNorm_result <- reactive({
  req(ct_data(), input$genes)
  ct_data <- ct_data()
  genes <- input$genes
  
  # Validate gene selection
  if (length(genes) < 2) {
    showNotification("Please select at least 2 genes", type = "error")
    return(NULL)
  }
  
  # Extract Ct values
  ct_values <- ct_data[, genes, drop = FALSE]
  
  # Convert to numeric matrix
  ct_matrix <- as.matrix(ct_values)
  ct_matrix <- apply(ct_matrix, 2, as.numeric)
  
  # Check for non-numeric values
  if (any(is.na(ct_matrix))) {
    invalid_count <- sum(is.na(ct_matrix))
    showNotification(
      paste("Warning:", invalid_count, "non-numeric values converted to NA"),
      type = "warning"
    )
  }
  
  # Remove rows with missing values
  if (anyNA(ct_matrix)) {
    ct_matrix <- na.omit(ct_matrix)
    showNotification(
      "Rows with missing values removed for geNorm analysis",
      type = "warning"
    )
  }
  
  # Validate matrix after cleaning
  if (nrow(ct_matrix) < 3) {
    showNotification(
      "Insufficient data after cleaning (need ≥ 3 samples)",
      type = "error"
    )
    return(NULL)
  }
  
  # Run geNorm with error handling
  tryCatch({
    if (input$use_ctrlgene) {
      ctrlGene::geNorm2(ct_matrix)
    } else {
      # Your custom geNorm implementation here
    }
  }, error = function(e) {
    showNotification(paste("geNorm Error:", e$message), type = "error")
    return(NULL)
  })
})

geNorm_table <- reactive({
  gnrm <- geNorm_result()
  
  # Properly handle M-values for the two most stable genes
  n <- nrow(gnrm)
  if(n >= 2) {
    last_m <- gnrm$Avg.M[n-1]  # Use the last calculated M-value
    gnrm$Avg.M[c(n-1, n)] <- last_m  # Assign to both stable genes
  }
  
  names(gnrm)[1] <- "Target"
  gnrm$Avg.M <- round(gnrm$Avg.M, 3)
  
  # Correct ranking (1 = most stable)
  gnrm$Rank <- dense_rank(gnrm$Avg.M)
  gnrm
})

pairwise_variation <- reactive({
  gnrm <- geNorm_table()
  rel_expr <- ct_data()
  n_genes <- nrow(gnrm)
  
  # Validate sufficient genes
  if(n_genes < 3) {
    return(data.frame(Variation = character(0), Value = numeric(0)))
  }
  
  # Order from most to least stable
  ordered_genes <- gnrm$Target[order(gnrm$Rank)]
  V_values <- numeric(0)
  V_labels <- character(0)
  
  # Calculate pairwise variations
  for (i in 2:(n_genes - 1)) {
    top_genes <- ordered_genes[1:(i+1)]
    
    # Corrected apply calls (removed extra parenthesis and added drop=FALSE)
    NF_i <- apply(
      rel_expr[, top_genes[1:i], drop = FALSE], 
      1, 
      function(x) exp(mean(log(x), na.rm = TRUE))
    )
    NF_i1 <- apply(
      rel_expr[, top_genes[1:(i+1)], drop = FALSE], 
      1, 
      function(x) exp(mean(log(x), na.rm = TRUE))
    )
    
    log_ratios <- log2(NF_i / NF_i1)
    V <- sd(log_ratios, na.rm = TRUE)
    
    V_values <- c(V_values, V)
    V_labels <- c(V_labels, paste0("V", i, "/", i+1))
  }
  
  data.frame(Variation = V_labels, Value = round(V_values, 4))
})

output$genorm_plot <- renderPlot({
  gnrm <- geNorm_table()
  # Order from least to most stable for plotting
  plot_data <- gnrm[order(-gnrm$Avg.M), ]
  
  plot(
    plot_data$Avg.M, 
    type = "o", 
    pch = 16, 
    col = "blue",
    xaxt = "n", 
    ylab = "Average Expression Stability (M)", 
    xlab = "Genes Ranked from Least to Most Stable",
    main = "geNorm: Expression Stability",
    las = 1,
    cex.main = 1.2,
    cex.lab = 1.1
  )
  axis(1, at = 1:nrow(plot_data), labels = plot_data$Target, las = 2)
  grid(nx = NA, ny = NULL)
})

output$genorm_pairwise_plot <- renderPlot({
  V_df <- pairwise_variation()
  
  if(nrow(V_df) == 0) {
    plot(0, 0, type = "n", xlab = "", ylab = "", axes = FALSE)
    text(0, 0, "Insufficient genes for pairwise variation", cex = 1.2)
    return()
  }
  
  bp <- barplot(V_df$Value, 
                names.arg = V_df$Variation, 
                col = "green",
                border = NA,
                main = "Optimal Reference Genes Determination",
                ylab = "Pairwise Variation (V)",
                ylim = c(0, max(V_df$Value) * 20),
                cex.main = 1.2,
                cex.lab = 1.1)
  text(bp, V_df$Value + 0.02 * max(V_df$Value), 
       labels = round(V_df$Value, 3),
       cex = 0.9)
  abline(h = 0.15, lty = 2, col = "red")
})

output$genorm_ranking_table <- renderTable({
  geNorm_table() %>% 
    arrange(Rank) %>%  # Show most stable genes first
    select(Target, `Stability (M)` = Avg.M, Rank)
})

# Download handlers (updated with improved plots)
output$download_genorm_plot <- downloadHandler(
  filename = function() {
    paste("geNorm_stability_", Sys.Date(), ".png", sep = "")
  },
  content = function(file) {
    gnrm <- geNorm_table()
    plot_data <- gnrm[order(-gnrm$Avg.M), ]
    
    png(file, width = 1200, height = 800, res = 100)
    par(mar = c(8, 4, 4, 2) + 0.1)  # Adjust bottom margin
    plot(plot_data$Avg.M, 
         type = "o", pch = 16, col = "blue",
         xaxt = "n", ylab = "Average Expression Stability (M)", 
         xlab = "", main = "geNorm: Expression Stability", las = 1)
    axis(1, at = 1:nrow(plot_data), labels = plot_data$Target, las = 2)
    grid(nx = NA, ny = NULL)
    dev.off()
  }
)

output$download_genorm_pairwise_plot <- downloadHandler(
  filename = function() {
    paste("geNorm_pairwise_variation_", Sys.Date(), ".png", sep = "")
  },
  content = function(file) {
    V_df <- pairwise_variation()
    if(nrow(V_df) == 0) return()
    
    png(file, width = 1000, height = 800, res = 100)
    bp <- barplot(V_df$Value, names.arg = V_df$Variation, 
                  col = "#1b9e77", border = NA,
                  main = "Optimal Reference Genes Determination",
                  ylab = "Pairwise Variation (V)",
                  ylim = c(0, max(V_df$Value) * 1.2))
    text(bp, V_df$Value + 0.02 * max(V_df$Value), 
         labels = round(V_df$Value, 3))
    abline(h = 0.15, lty = 2, col = "red")
    dev.off()
  }
)

output$download_genorm_ranking_table <- downloadHandler(
  filename = function() {
    paste("geNorm_ranking_", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    write.csv(
      geNorm_table() %>% arrange(Rank),
      file, 
      row.names = FALSE
    )
  }
)

###mixed model stability
library(shiny)
library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(boot)
library(readr)
library(car)
library(purrr)

mixedmodel_stability <- function(data, 
                                 target, 
                                 response,
                                 fixed.effects,
                                 random.effect,
                                 form, 
                                 reduced.model, 
                                 icc.model = NULL, 
                                 hypothesis.test = "LRT", 
                                 LRT.type = "global",
                                 p.threshold = 0.05,
                                 critical.t = 1.959964, 
                                 icc.interval = 0.95, 
                                 icc.type = "norm", 
                                 n.genes = 2, 
                                 n.sims = 500, 
                                 progress = TRUE) {
  
  results <- list()
  data <- as.data.frame(data)
  
  form <- as.formula(form)
  reduced.model <- as.formula(reduced.model)
  if (is.null(icc.model)) icc.model <- reduced.model else icc.model <- as.formula(icc.model)
  
  for (n in seq_along(n.genes)) {
    
    combinations <- combn(unique(data[[target]]), n.genes[n])
    
    icc.calc <- function(model) {
      vcomp <- as.numeric(data.frame(lme4::VarCorr(model))[, 4])
      icc <- vcomp[1] / sum(vcomp)
      return(icc)
    }
    
    results.icc <- data.frame(
      gene.combination = rep(NA, ncol(combinations)),
      icc = NA, icc.l = NA, icc.u = NA
    )
    
    if (hypothesis.test == "LRT") {
      if (LRT.type == "add") {
        terms <- c(fixed.effects, paste(target, fixed.effects, sep = ":"))
        results.lrt <- data.frame(matrix(NA, nrow = ncol(combinations), ncol = length(terms)))
        colnames(results.lrt) <- terms
      } else if (LRT.type == "global") {
        results.lrt <- data.frame(p.val = rep(NA, ncol(combinations)))
      } else if (LRT.type == "interactions") {
        results.lrt <- data.frame(syst.pval = NA, inter.pval = NA)
      }
    } else if (hypothesis.test == "wald.t") {
      results.wald.t <- data.frame(coef = NA, t.val = NA, p.val = NA)
    } else if (hypothesis.test == "z.as.t") {
      results.z.as.t <- data.frame(coef = NA, t.val = NA)
    }
    
    if (progress) pb <- txtProgressBar(min = 0, max = ncol(combinations), style = 3)
    
    for (i in 1:ncol(combinations)) {
      genes <- combinations[, i]
      results.icc[i, 1] <- paste(genes, collapse = ":")
      subset.data <- subset(data, data[[target]] %in% genes)
      
      h.test <- FALSE  # default
      
      if (!(hypothesis.test %in% c("LRT", "wald.t", "z.as.t", "anova", "wald.chi")))
        stop("Invalid hypothesis.test specified.")
      
      if (hypothesis.test == "LRT") {
        response.target <- paste0(response, "~", target)
        systematic <- paste(fixed.effects, collapse = "+")
        interactions <- paste(target, fixed.effects, sep = ":", collapse = "+")
        full.f <- as.formula(paste(response.target, systematic, interactions, random.effect, sep = "+"))
        systematic.f <- as.formula(paste(response.target, systematic, random.effect, sep = "+"))
        reduced.f <- as.formula(paste(response.target, random.effect, sep = "+"))
        
        if (LRT.type == "interactions") {
          reduced.mod <- lmer(reduced.f, data = subset.data, REML = FALSE)
          systematic.mod <- lmer(systematic.f, data = subset.data, REML = FALSE)
          full.mod <- lmer(full.f, data = subset.data, REML = FALSE)
          
          results.lrt[i, ] <- c(
            anova(reduced.mod, systematic.mod)[2, 8],
            anova(systematic.mod, full.mod)[2, 8]
          )
        } else if (LRT.type == "add") {
          reduced.mod <- lmer(reduced.f, data = subset.data, REML = FALSE)
          results.lrt[i, ] <- add1(reduced.mod, scope = terms, test = "Chi")[-1, 4]
        } else if (LRT.type == "global") {
          reduced.mod <- lmer(reduced.f, data = subset.data, REML = FALSE)
          full.mod <- lmer(full.f, data = subset.data, REML = FALSE)
          results.lrt[i, 1] <- anova(full.mod, reduced.mod)[2, 8]
        }
        
        h.test <- any(results.lrt[i, ] < p.threshold, na.rm = TRUE)
      }
      
      if (hypothesis.test == "wald.t") {
        full.model <- lmerTest::lmer(form, data = subset.data, REML = FALSE)
        coef.table <- data.frame(summary(full.model)$coef)
        
        if (ncol(coef.table) == 5) {
          coef.table$coefs <- rownames(coef.table)
          coef.table <- coef.table[!(coef.table$coefs %in% c("(Intercept)", paste0(target, genes))), ]
          results.wald.t[i, ] <- coef.table[which.min(coef.table[, 5]), c(6, 4, 5)]
          h.test <- any(results.wald.t[i, 3] < p.threshold, na.rm = TRUE)
        } else {
          warning(paste("Fallback to critical t =", critical.t))
        }
      }
      
      if (hypothesis.test == "z.as.t") {
        full.model <- lmer(form, data = subset.data, REML = TRUE)
        coef.table <- data.frame(summary(full.model)$coef)
        coef.table$coefs <- rownames(coef.table)
        coef.table <- coef.table[!(coef.table$coefs %in% c("(Intercept)", paste0(target, genes))), ]
        results.z.as.t[i, ] <- coef.table[which.max(abs(coef.table[, 3])), c(4, 3)]
        h.test <- abs(results.z.as.t[i, 2]) > critical.t
      }
      
      if (hypothesis.test == "anova") {
        full.model <- lmerTest::lmer(form, data = subset.data, REML = FALSE)
        coef.table <- anova(full.model)
        p.vals <- coef.table[!(rownames(coef.table) %in% c(target)), "Pr(>F)"]
        h.test <- any(p.vals < p.threshold, na.rm = TRUE)
      }
      
      if (hypothesis.test == "wald.chi") {
        full.model <- lmer(form, data = subset.data, REML = FALSE)
        coef.table <- car::Anova(full.model, type = 2)
        p.vals <- coef.table[!(rownames(coef.table) %in% c("(Intercept)", target)), "Pr(>Chisq)"]
        h.test <- any(p.vals < p.threshold, na.rm = TRUE)
      }
      
      if (!h.test) {
        icc.mod <- lmer(icc.model, data = subset.data, REML = TRUE)
        b <- bootMer(icc.mod, icc.calc, nsim = n.sims, use.u = FALSE)
        ci <- boot::boot.ci(b, conf = icc.interval, type = icc.type)
        
        if (icc.type == "norm") {
          cis <- ci$normal[2:3]
        } else if (icc.type == "basic") {
          cis <- ci$basic[4:5]
        } else if (icc.type == "perc") {
          cis <- ci$percent[4:5]
        } else stop("Unsupported ICC CI type.")
        
        results.icc[i, 2:4] <- c(b$t0, cis)
      }
      
      if (progress) setTxtProgressBar(pb, i)
    }
    
    if (progress) close(pb)
    
    r <- results.icc
    if (exists("results.lrt")) r <- cbind(r, results.lrt)
    if (exists("results.wald.t")) r <- cbind(r, results.wald.t)
    if (exists("results.z.as.t")) r <- cbind(r, results.z.as.t)
    r$n.genes <- n.genes[n]
    results[[n]] <- r
  }
  
  bind_rows(results)
}

ct_data <- reactive({
  req(input$ct_file)
  read_csv(input$ct_file$datapath)
})

df_long <- reactive({
  req(ct_data())
  ct_data() %>%
    pivot_longer(cols = -c(Group, Sample), names_to = "Gene", values_to = "Expression") %>%
    mutate(Expression = as.numeric(Expression)) %>%
    filter(!is.na(Expression))
})

mixed_model_results <- reactive({
  req(df_long())
  
  # Use a UI progress bar because the internal function uses txtProgressBar
  withProgress(message = "Running Mixed Model stability (this may take a while)", value = 0, {
    incProgress(0.02, detail = "Preparing combinations...")
    Sys.sleep(0.01)
    
    # call the mixedmodel_stability with progress = FALSE to avoid console txtProgressBar
    res <- mixedmodel_stability(
      data = df_long(),
      target = "Gene",
      response = "Expression",
      fixed.effects = "Group",
      random.effect = "(1|Sample)",
      form = "Expression ~ Group * Gene + (1|Sample)",
      reduced.model = "Expression ~ Gene + (1|Sample)",
      hypothesis.test = "LRT",
      LRT.type = "global",
      p.threshold = 0.05,
      n.genes = 2,
      n.sims = 500,
      progress = FALSE  # important: don't use txtProgressBar
    )
    
    # If you want to show incremental updates you can call incProgress repeatedly,
    # but since mixedmodel_stability itself does heavy work, we give a few steps:
    incProgress(0.70, detail = "Bootstrapping ICCs (inside mixed model)...")
    Sys.sleep(0.01)
    incProgress(0.25, detail = "Finalizing mixed-model results...")
    Sys.sleep(0.01)
    
    res
  })
})

ranked_mixed_model_results <- reactive({
  req(mixed_model_results())
  df <- mixed_model_results()
  
  if (!"icc.l" %in% names(df) || !"p.val" %in% names(df)) {
    return(data.frame(Message = "Missing required columns 'icc.l' or 'p.val'"))
  }
  
  df_filtered <- df %>%
    filter(p.val > 0.05) %>%
    arrange(desc(icc.l)) %>%
    mutate(Rank = row_number()) %>%
    select(Rank, everything())
  
  df_filtered
})

# ADDED: Mixed Model Individual Gene Ranking
mixedModelGeneRanks <- reactive({
  req(ranked_mixed_model_results())
  
  # Calculate gene participation in combinations and their ICC values
  gene_icc_scores <- ranked_mixed_model_results() %>%
    separate_rows(gene.combination, sep = ":") %>%
    group_by(Gene = gene.combination) %>%
    summarise(
      MixedModel_Score = mean(icc.l, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      MixedModel_rank = dense_rank(desc(MixedModel_Score))  # Fixed syntax
    )
  
  # Handle genes not in any valid combinations
  all_genes <- unique(ct_data() %>% select(-Group, -Sample) %>% names())
  missing_genes <- setdiff(all_genes, gene_icc_scores$Gene)
  
  if (length(missing_genes) > 0) {
    gene_icc_scores <- bind_rows(
      gene_icc_scores,
      tibble(Gene = missing_genes, 
             MixedModel_Score = NA_real_,
             MixedModel_rank = NA_integer_)
    )
  }
  
  gene_icc_scores
})

output$mixedModelTable <- renderTable({
  req(mixed_model_results())
  head(mixed_model_results(), 100)
})

output$mixedModelTable <- renderTable({
  req(ranked_mixed_model_results())
  head(ranked_mixed_model_results(), 100)
})
  
output$bestMixedModelCombination <- renderText({
  top <- head(ranked_mixed_model_results(), 1)
  if (nrow(top) == 0) {
    return("No suitable gene combinations found (all p-values ≤ 0.05).")
  }
  paste0("Best reference gene combination: ", top$gene.combination,
         " (ICC lower bound = ", round(top$icc.l, 3), ")")
})

output$download_ranked_mixed_model <- downloadHandler(
  filename = function() {
    paste("ranked_mixed_model_results_", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    write.csv(ranked_mixed_model_results(), file, row.names = FALSE)
  }
)

output$iccRankPlot <- renderPlot({
  req(ranked_mixed_model_results())
  df <- ranked_mixed_model_results()
  if (nrow(df) == 0) {
    plot.new(); text(0.5, 0.5, "No gene combinations to plot", cex = 1.2)
    return()
  }
  
  # Prepare data for plotting (ensure Rank exists)
  df <- df %>%
    mutate(Rank = row_number()) %>%
    arrange(Rank)
  
  # Build plot object (assign to variable p_icc_rank)
  p_icc_rank <- ggplot(df, aes(x = Rank, y = icc.l, label = gene.combination)) +
    geom_point(aes(color = p.val), size = 3) +
    geom_text(nudge_y = 0.02, hjust = 0, size = 3, check_overlap = TRUE) +
    scale_color_gradient(low = "blue", high = "red", name = "LRT p-value") +
    labs(
      x = "Rank",
      y = "ICC Lower Bound (95% CI)",
      title = "Gene Combination Ranking by ICC Lower Bound",
      subtitle = "Combinations with LRT p-value > 0.05 shown"
    ) +
    theme_minimal()
  
  print(p_icc_rank)   # render the plot
})

output$download_icc_rank_plot <- downloadHandler(
  filename = function() paste0("MixedModel_ICC_Rank_Plot_", Sys.Date(), ".png"),
  content = function(file) {
    req(ranked_mixed_model_results())
    df <- ranked_mixed_model_results()
    if (nrow(df) == 0) {
      # create a simple empty png with message
      png(file, width = 1200, height = 800, res = 150)
      plot.new(); text(0.5, 0.5, "No gene combinations to plot", cex = 1.2)
      dev.off()
      return()
    }
    
    df <- df %>%
      mutate(Rank = row_number()) %>%
      arrange(Rank)
    
    p_icc_rank <- ggplot(df, aes(x = Rank, y = icc.l, label = gene.combination)) +
      geom_point(aes(color = p.val), size = 3) +
      geom_text(nudge_y = 0.02, hjust = 0, size = 3, check_overlap = TRUE) +
      scale_color_gradient(low = "blue", high = "red", name = "LRT p-value") +
      labs(
        x = "Rank",
        y = "ICC Lower Bound (95% CI)",
        title = "Gene Combination Ranking by ICC Lower Bound",
        subtitle = "Combinations with LRT p-value > 0.05 shown"
      ) +
      theme_minimal()
    
    # Save explicitly using ggsave (uses the object p_icc_rank, not last_plot())
    ggsave(filename = file, plot = p_icc_rank, width = 10, height = 8, dpi = 300)
  }
)

df_long <- reactive({
  req(ct_data())
  data_long <- ct_data() %>%
    pivot_longer(
      cols = -c(Group, Sample), 
      names_to = "Gene", 
      values_to = "Expression"
    )
  print(colnames(data_long))  # Debug column names
  data_long
})

computeICCandRanks <- function(df) {
  tryCatch({
    # 1. Validate input structure
    req(df)
    if (!is.data.frame(df)) 
      stop("Input must be a data frame")
    
    # 2. Standardize column names (case-insensitive)
    colnames(df) <- tolower(colnames(df))
    required_cols <- c("sample", "gene", "expression")
    missing_cols <- setdiff(required_cols, colnames(df))
    
    if (length(missing_cols) > 0)
      stop("Missing columns: ", paste(missing_cols, collapse = ", "))
    
    # 3. Clean and format data
    clean_df <- df %>%
      rename(
        sample = "sample",
        gene = "gene",
        expression = "expression"
      ) %>%
      filter(
        !is.na(expression),
        !is.na(sample),
        !is.na(gene)
      ) %>%
      mutate(
        gene = as.character(gene),
        sample = as.character(sample)
      )
    
    # 4. Calculate ICC for each gene
    icc_results <- clean_df %>%
      group_by(gene) %>%
      group_modify(~{
        gene_data <- .x
        
        # Skip genes with insufficient data
        if (n_distinct(gene_data$sample) < 2 || nrow(gene_data) < 3) {
          return(tibble(
            gene = unique(gene_data$gene),
            ICC_Score = NA_real_,
            error = "Insufficient data"
          ))
        }
        
        # Build appropriate model formula
        model_formula <- if ("group" %in% colnames(gene_data)) {
          expression ~ (1 | sample) + (1 | group)
        } else {
          expression ~ (1 | sample)
        }
        
        # Fit model with error handling
        model <- tryCatch({
          lme4::lmer(
            model_formula,
            data = gene_data,
            control = lmerControl(calc.derivs = FALSE)
          )
        }, error = function(e) NULL)
        
        # Calculate ICC if model succeeded
        if (!is.null(model)) {
          icc_value <- tryCatch({
            performance::icc(model)$ICC_adjusted
          }, error = function(e) NA_real_)
          
          tibble(
            ICC_Score = icc_value,
            error = NA_character_
          )
        } else {
          tibble(
            ICC_Score = NA_real_,
            error = "Model convergence failed"
          )
        }
      }) %>%
      ungroup()
    
    # 5. Calculate ranks (handle NAs)
    icc_results %>%
      mutate(
        ICC_rank = dense_rank(desc(ICC_Score)),
        ICC_rank = ifelse(is.na(ICC_Score), NA_integer_, ICC_rank)
      ) %>%
      select(
        Gene = gene,
        ICC_Score,
        ICC_rank,
        error
      )
    
  }, error = function(e) {
    message("Error in computeICCandRanks: ", e$message)
    return(tibble(
      Gene = character(),
      ICC_Score = numeric(),
      ICC_rank = integer(),
      error = character()
    ))
  })
}

mixedModelResults <- reactive({
  req(df_long())
  results <- computeICCandRanks(df_long())  # From your mixed model script
  tibble::tibble(
    Gene               = results$Gene,
    ICC_Score          = results$ICC,
    ICC_StabilityRank  = results$Rank
  )
})

iccResults <- reactive({
  req(df_long())
  
  results <- computeICCandRanks(df_long())
  
  # Optional: Show warnings to user
  if (any(!is.na(results$error))) {
    errors <- results %>%
      filter(!is.na(error)) %>%
      distinct(Gene, error)
    
    showNotification(
      paste("ICC warnings:", 
            paste(errors$Gene, ":", errors$error, collapse = "; ")),
      type = "warning"
    )
  }
  
  select(results, -error)  # Remove error column for downstream use
})

##Normfinder
#NormFinder implementation
  Normfinder = function(filename, Groups = TRUE, ctVal = TRUE, pStabLim = 0.25) {
    dat0 = read.table(filename, header = TRUE, row.names = 1, colClasses = "character")
    ntotal = dim(dat0)[2]
    k0 = dim(dat0)[1]
    
    if (Groups) {
      ngenes = k0 - 1
      genenames = rownames(dat0)[-k0]
      grId = dat0[k0, ]
      dat0 = dat0[-k0, ]
    } else {
      ngenes = k0
      genenames = rownames(dat0)
      grId = rep(1, ntotal)
    }
    
    dat = matrix(as.numeric(unlist(dat0)), ngenes, ntotal)
    if (!ctVal) dat = log2(dat)
    samplenames = colnames(dat0)
    grId = factor(unlist(grId))
    groupnames = levels(grId)
    ngr = length(groupnames)
    
    nsamples = sapply(groupnames, function(g) sum(grId == g))
    
    MakeStab = function(da) {
      sampleavg = colMeans(da)
      genegroupavg = sapply(1:ngr, function(group) rowMeans(da[, grId == groupnames[group], drop = FALSE]))
      groupavg = sapply(1:ngr, function(group) mean(da[, grId == groupnames[group]]))
      
      GGvar = matrix(0, ngenes, ngr)
      for (group in 1:ngr) {
        grset = (grId == groupnames[group])
        GGvar[, group] = sapply(1:ngenes, function(gene)
          sum((da[gene, grset] - genegroupavg[gene, group] - sampleavg[grset] + groupavg[group])^2) /
            (nsamples[group] - 1))
        GGvar[, group] = (GGvar[, group] - mean(GGvar[, group])) / (1 - 2 / ngenes)
      }
      
      genegroupMinvar = matrix(0, ngenes, ngr)
      for (group in 1:ngr) {
        z = da[, grId == groupnames[group]]
        for (gene in 1:ngenes) {
          varpair = sapply(1:ngenes, function(g1) var(z[gene, ] - z[g1, ]))
          genegroupMinvar[gene, group] = min(varpair[-gene]) / 4
        }
      }
      
      GGvar = ifelse(GGvar < 0, genegroupMinvar, GGvar)
      
      dif = genegroupavg
      difgeneavg = rowMeans(dif)
      difgroupavg = colMeans(dif)
      difavg = mean(dif)
      for (gene in 1:ngenes) {
        for (group in 1:ngr) {
          dif[gene, group] = dif[gene, group] - difgeneavg[gene] - difgroupavg[group] + difavg
        }
      }
      
      nsampMatrix = matrix(rep(nsamples, ngenes), ngenes, ngr, byrow = TRUE)
      vardif = GGvar / nsampMatrix
      gamma = sum(dif * dif) / ((ngr - 1) * (ngenes - 1)) - sum(vardif) / (ngenes * ngr)
      gamma = max(0, gamma)
      
      difnew = dif * gamma / (gamma + vardif)
      varnew = vardif + gamma * vardif / (gamma + vardif)
      Ostab = rowMeans(abs(difnew) + sqrt(varnew))
      
      mud = apply(dif, 1, function(x) 2 * max(abs(x)))
      genevar = rowSums((nsamples - 1) * GGvar) / (sum(nsamples) - ngr)
      Gsd = sqrt(genevar)
      
      cbind(mud, Gsd, Ostab, rep(gamma, ngenes), GGvar, dif)
    }
    
    MakeComb2 = function(g1, g2, res) {
      gam = res[1, 4]
      d1 = res[g1, (4 + ngr + 1):(4 + ngr + ngr)]
      d2 = res[g2, (4 + ngr + 1):(4 + ngr + ngr)]
      s1 = res[g1, (4 + 1):(4 + ngr)]
      s2 = res[g2, (4 + 1):(4 + ngr)]
      rho = abs(gam * d1 / (gam + s1 / nsamples) + gam * d2 / (gam + s2 / nsamples)) *
        sqrt(ngenes / (ngenes - 2)) / 2
      rho = rho + sqrt(s1 / nsamples + gam * s1 / (nsamples * gam + s1) +
                         s2 / nsamples + gam * s2 / (nsamples * gam + s2)) / 2
      mean(rho)
    }
    
    MakeStabOne = function(da) {
      sampleavg = colMeans(da)
      geneavg = rowMeans(da)
      totalavg = mean(da)
      
      genevar0 = sapply(1:ngenes, function(gene)
        sum((dat[gene, ] - geneavg[gene] - sampleavg + totalavg)^2) /
          ((ntotal - 1) * (1 - 2 / ngenes)))
      genevar = genevar0 - mean(genevar0) / (ngenes - 1)
      
      geneMinvar = sapply(1:ngenes, function(gene) {
        varpair = sapply(1:ngenes, function(g1) var(da[gene, ] - da[g1, ]))
        min(varpair[-gene]) / 4
      })
      
      ifelse(genevar < 0, geneMinvar, genevar)
    }
    
    if (ngr > 1) {
      res = MakeStab(dat)
      gcand = which(res[, 3] < pStabLim)
      if (length(gcand) < 4) {
        if (ngenes > 3) {
          li = sort(res[, 3])[4]
          gcand = which(res[, 3] <= li)
        } else {
          gcand = 1:ngenes
        }
      }
      
      vv2 = do.call(rbind, combn(gcand, 2, simplify = FALSE, FUN = function(pair) {
        c(pair, MakeComb2(pair[1], pair[2], res))
      }))
      
      ord = order(res[, 3])
      list(
        Ordered = data.frame(GroupDif = round(res[ord, 1], 2), GroupSD = round(res[ord, 2], 2),
                             Stability = round(res[ord, 3], 2), row.names = genenames[ord]),
        UnOrdered = data.frame(GroupDif = round(res[, 1], 2), GroupSD = round(res[, 2], 2),
                               Stability = round(res[, 3], 2),
                               IGroupSD = round(sqrt(res[, (4 + 1):(4 + ngr)]), 2),
                               IGroupDif = round(res[, (4 + ngr + 1):(4 + ngr + ngr)], 2),
                               row.names = genenames),
        PairOfGenes = data.frame(Gene1 = genenames[vv2[, 1]], Gene2 = genenames[vv2[, 2]], 
                                 Stability = round(vv2[, 3], 2))
      )
    } else {
      sigma = sqrt(MakeStabOne(dat))
      siglim = min(sigma) + 0.1
      gcand = which(sigma < siglim)
      
      if (length(gcand) >= 2 && ngenes > 3) {
        vv2 = do.call(rbind, combn(gcand, 2, simplify = FALSE, FUN = function(pair) {
          dat1 = rbind(dat[-pair, ], colMeans(dat[pair, ]))
          c(pair, sqrt(MakeStabOne(dat1))[ngenes - 1])
        }))
        ord = order(sigma)
        list(
          Ordered = data.frame(GroupSD = round(sigma[ord], 2), row.names = genenames[ord]),
          PairOfGenes = data.frame(Gene1 = genenames[vv2[, 1]], Gene2 = genenames[vv2[, 2]], 
                                   GroupSD = round(vv2[, 3], 2))
        )
      } else {
        ord = order(sigma)
        list(Ordered = data.frame(GroupSD = round(sigma[ord], 2), row.names = genenames[ord]))
      }
    }
  }
  
  # Wrapper function to run NormFinder with prepared Ct data
  runNormFinder <- function(ct_data) {
    sample_info <- ct_data[, c("Sample", "Group")]
    expr_data <- ct_data[, -c(1, 2)]
    gene_names <- colnames(expr_data)
    
    expr_data <- as.data.frame(lapply(expr_data, as.numeric))
    expr_matrix <- t(expr_data)
    rownames(expr_matrix) <- gene_names
    colnames(expr_matrix) <- sample_info$Sample
    
    expr_df <- as.data.frame(expr_matrix)
    expr_df <- rbind(expr_df, Group = as.character(sample_info$Group))
    
    temp_file <- tempfile(fileext = ".txt")
    write.table(expr_df, file = temp_file, sep = "\t", quote = FALSE, col.names = NA)
    
    Normfinder(temp_file, Groups = TRUE, ctVal = TRUE)
  }
  
  # render
  # Reactive expression to store NormFinder result
  norm_result <- reactive({
    req(ct_data())
    runNormFinder(ct_data())
  })
  
  # ---- Render Tables ----
  
  output$norm_ordered <- DT::renderDataTable({
    req(norm_result())
    norm_result()$Ordered
  }, options = list(pageLength = 10))
  
  output$norm_unordered <- DT::renderDataTable({
    req(norm_result())
    norm_result()$UnOrdered
  }, options = list(pageLength = 10))
  
  output$norm_pairs <- DT::renderDataTable({
    req(norm_result())
    norm_result()$PairOfGenes
  }, options = list(pageLength = 10))
  
  output$norm_ranked_table <- DT::renderDataTable({
    req(norm_result())
    df <- norm_result()$UnOrdered
    df$Rank_IGroupSD.v1 <- rank(df$IGroupSD.V1, ties.method = "min", na.last = "keep")
    df$Rank_IGroupSD.v2 <- rank(df$IGroupSD.V2, ties.method = "min", na.last = "keep")
    df$Rank_Stability   <- rank(df$Stability, ties.method = "min", na.last = "keep")
    df[, c("Stability", "IGroupSD.V1", "IGroupSD.V2", 
           "Rank_IGroupSD.v1", "Rank_IGroupSD.v2", "Rank_Stability")]
  }, options = list(pageLength = 10))
  
  # ---- Download Handlers ----
  
  output$download_norm_ordered <- downloadHandler(
    filename = function() paste0("NormFinder_Ordered_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(norm_result()$Ordered, file, row.names = TRUE)
    }
  )
  
  output$download_norm_unordered <- downloadHandler(
    filename = function() paste0("NormFinder_Unordered_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(norm_result()$UnOrdered, file, row.names = TRUE)
    }
  )
  
  output$download_norm_pairs <- downloadHandler(
    filename = function() paste0("NormFinder_TopPairs_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(norm_result()$PairOfGenes, file, row.names = FALSE)
    }
  )
  
  output$download_norm_ranked_table <- downloadHandler(
    filename = function() paste0("NormFinder_Ranked_", Sys.Date(), ".csv"),
    content = function(file) {
      df <- norm_result()$UnOrdered
      df$Rank_IGroupSD.v1 <- rank(df$IGroupSD.V1, ties.method = "min", na.last = "keep")
      df$Rank_IGroupSD.v2 <- rank(df$IGroupSD.V2, ties.method = "min", na.last = "keep")
      df$Rank_Stability   <- rank(df$Stability, ties.method = "min", na.last = "keep")
      ranked_df <- df[, c("Stability", "IGroupSD.V1", "IGroupSD.V2", 
                          "Rank_IGroupSD.v1", "Rank_IGroupSD.v2", "Rank_Stability")]
      write.csv(ranked_df, file, row.names = TRUE)
    }
  )
  

##########deltact
# Step 1: Reactive CT data
  ct_data <- reactive({
    req(input$ct_file)
    read_csv(input$ct_file$datapath)
  })
  
# Step 2: Update gene input selection
  observe({
    req(ct_data())
    genes <- colnames(ct_data())[!colnames(ct_data()) %in% c("Sample", "Group")]
    updateSelectInput(session, "genes", choices = genes, selected = genes)
  })
  
# Step 3: Preprocess CT data
  deltact_prep <- function(ct_data) {
    sample_info <- ct_data[, c("Sample", "Group")]
    expr_data <- ct_data[, !(names(ct_data) %in% c("Sample", "Group"))] %>% mutate_all(as.numeric)
    gene_names <- colnames(expr_data)
    
    pivot_longer(cbind(sample_info, expr_data), 
                 cols = all_of(gene_names),
                 names_to = "Gene", values_to = "Ct")
  }
  
# Step 4: Compute ΔCt
  compute_delta_ct <- function(data, genes) {
    wide_data <- data %>% 
      pivot_wider(names_from = Gene, values_from = Ct)
    
    delta_ct_list <- list()
    for (i in seq_along(genes)) {
      for (j in seq_along(genes)) {
        if (i != j) {
          pair_name <- paste(genes[i], "-", genes[j], sep = "")
          delta_ct_list[[pair_name]] <- wide_data[[genes[i]]] - wide_data[[genes[j]]]
        }
      }
    }
    delta_df <- as.data.frame(delta_ct_list)
    delta_df$Group <- wide_data$Group
    delta_df
  }
  
# Step 5: Summarize and Rank
  summarize_delta_ct <- function(delta_ct_df) {
    long_df <- delta_ct_df %>% 
      pivot_longer(cols = -Group, names_to = "Gene_Pair", values_to = "ΔCt") %>%
      separate(Gene_Pair, into = c("First_Gene", "Second_Gene"), sep = "\\.", remove = FALSE)
    
    summary_stats <- long_df %>%
      group_by(First_Gene) %>%
      summarise(
        Mean_ΔCt = mean(ΔCt, na.rm = TRUE),
        StdDev_ΔCt = sd(ΔCt, na.rm = TRUE),
        .groups = "drop"
      )
    
    pairwise_stddev <- long_df %>%
      group_by(First_Gene, Second_Gene) %>%
      summarise(Pairwise_StdDev = sd(ΔCt, na.rm = TRUE), .groups = "drop")
    
    mean_stddev_per_gene <- pairwise_stddev %>%
      group_by(First_Gene) %>%
      summarise(Mean_of_StdDev = mean(Pairwise_StdDev, na.rm = TRUE), .groups = "drop")
    
    final_summary <- left_join(summary_stats, mean_stddev_per_gene, by = "First_Gene") %>%
      arrange(Mean_of_StdDev) %>%
      mutate(Rank = row_number())
    
    list(
      delta_ct_long = left_join(long_df, final_summary, by = "First_Gene"),
      summary_stats = final_summary
    )
  }
  
# Step 6: Plot ΔCt
  plot_delta_ct <- function(delta_ct_long, summary_stats, annotation_y = 12) {
    group_levels <- unique(delta_ct_long$Group)
    color_palette <- scales::hue_pal()(length(group_levels))  # generate n colors
    
    ggplot(delta_ct_long, aes(x = Gene_Pair, y = ΔCt, fill = Group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(position = position_jitter(0.2), alpha = 0.5, size = 1) +
      scale_fill_manual(values = setNames(color_palette, group_levels)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(x = "Gene Pair", y = expression(ΔC[T])) +
      geom_text(data = summary_stats,
                aes(x = First_Gene, y = annotation_y,
                    label = paste0(First_Gene,
                                   "\nx̄ΔCt: ", round(Mean_ΔCt, 2),
                                   "\nx̄SD: ", round(Mean_of_StdDev, 2))),
                inherit.aes = FALSE, size = 5, hjust = 0, color = "black")
  }
  
  output$delta_ct_plot <- renderPlot({
    plot_delta_ct(delta_ct_long, summary_stats)
  })
  
  output$download_delta_ct_plot <- downloadHandler(
    filename = function() {
      paste0("delta_ct_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      data <- processed_data()
      ggsave(file, plot = plot_delta_ct(data$delta_ct_long, data$summary_stats),
             width = 10, height = 7, dpi = 300)
    }
  )

# Step 7: Reactive pipeline
  processed_data <- reactive({
    req(ct_data(), input$genes)
    long_data <- deltact_prep(ct_data())
    delta_ct_df <- compute_delta_ct(long_data, input$genes)
    summarize_delta_ct(delta_ct_df)
  })
  
# Step 8: Output Summary Statistics Table
  output$summaryTable <- DT::renderDataTable({
    DT::datatable(
      processed_data()$summary_stats,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        dom = 'tip',
        columnDefs = list(list(className = 'dt-center', targets = "_all")),
        order = list(list(
          which(colnames(processed_data()$summary_stats) == "Rank") - 1, 'asc'
        ))
      ),
      rownames = FALSE,
      class = 'cell-border stripe'
    )
  })
  
  output$download_delta_ct_summary <- downloadHandler(
    filename = function() paste0("DeltaCt_Summary_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(processed_data()$summary_stats, file, row.names = FALSE)
    }
  )
  
# Step 9: Plot output (optional if needed)
  output$deltaCtPlot <- renderPlot({
    req(processed_data())
    plot_delta_ct(
      processed_data()$delta_ct_long,
      processed_data()$summary_stats
    )
  })


########## --- Standard Curve Plot ---
  output$std_curve_plot <- renderPlot({
    std <- std_curve_data()
    std$log_Conc <- log10(std$Concentration)
    graph_list <- list()
    
    for (gene in unique(std$Gene)) {
      gene_data <- subset(std, Gene == gene)
      model <- lm(Ct ~ log_Conc, data = gene_data)
      r_squared <- summary(model)$r.squared
      eq <- paste0("y = ", round(coef(model)[2], 4), "x + ", round(coef(model)[1], 4))
      
      g <- ggplot(gene_data, aes(x = log_Conc, y = Ct)) +
        geom_point(size = 2) +
        geom_smooth(method = "lm", color = "black") +
        labs(title = gene,
             x = expression(Log[10] ~ "cDNA concentration (" * mu * "g/" * mu * "l)"),
             y = expression(bar(x) * C[T])) +
        annotate("text", x = min(gene_data$log_Conc), y = max(gene_data$Ct) - 1, 
                 label = eq, hjust = 0, size = 4) +
        annotate("text", x = min(gene_data$log_Conc), y = max(gene_data$Ct), 
                 label = paste("R² =", round(r_squared, 4)), hjust = 0, size = 4) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 15),
              plot.title = element_text(face = "italic", size = 15))
      
      graph_list[[gene]] <- g
    }
    
    ggarrange(plotlist = graph_list, ncol = 3, nrow = ceiling(length(graph_list)/3))
  })
  
  output$download_efficiency <- downloadHandler(
    filename = function() {
      paste0("PCR_Efficiency_", Sys.Date(), ".csv")
    },
    content = function(file) {
      std <- std_curve_data()
      std$log_Conc <- log10(std$Concentration)
      
      result <- data.frame(Gene = character(), Slope = numeric(), Efficiency = numeric(), stringsAsFactors = FALSE)
      
      for (gene in unique(std$Gene)) {
        gene_data <- subset(std, Gene == gene)
        model <- lm(Ct ~ log_Conc, data = gene_data)
        slope <- coef(model)[["log_Conc"]]
        efficiency <- (10^(-1 / slope) - 1) * 100
        result <- rbind(result, data.frame(Gene = gene, Slope = slope, Efficiency = efficiency))
      }
      
      write.csv(result, file, row.names = FALSE)
    }
  )
  output$download_std_curve_plot <- downloadHandler(
    filename = function() {
      paste0("Standard_Curve_", Sys.Date(), ".png")
    },
    content = function(file) {
      std <- std_curve_data()
      std$log_Conc <- log10(std$Concentration)
      graph_list <- list()
      
      for (gene in unique(std$Gene)) {
        gene_data <- subset(std, Gene == gene)
        model <- lm(Ct ~ log_Conc, data = gene_data)
        r_squared <- summary(model)$r.squared
        eq <- paste0("y = ", round(coef(model)[2], 4), "x + ", round(coef(model)[1], 4))
        
        g <- ggplot(gene_data, aes(x = log_Conc, y = Ct)) +
          geom_point(size = 2) +
          geom_smooth(method = "lm", color = "black") +
          labs(title = gene,
               x = expression(Log[10] ~ "cDNA concentration (" * mu * "g/" * mu * "l)"),
               y = expression(bar(x) * C[T])) +
          annotate("text", x = min(gene_data$log_Conc), y = max(gene_data$Ct) - 1, 
                   label = eq, hjust = 0, size = 4) +
          annotate("text", x = min(gene_data$log_Conc), y = max(gene_data$Ct), 
                   label = paste("R² =", round(r_squared, 4)), hjust = 0, size = 4) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title = element_text(size = 15),
                plot.title = element_text(face = "italic", size = 15))
        
        graph_list[[gene]] <- g
      }
      
      final_plot <- ggarrange(plotlist = graph_list, ncol = 3, nrow = ceiling(length(graph_list)/3))
      
      ggsave(file, plot = final_plot, width = 15, height = 10, dpi = 300)
    }
  )
  
  # --- Ct Distribution Plot ---
  # --- Ct Distribution Plot (All Genes in Single Plot) ---
  output$ct_dist_plot <- renderPlot({
    req(ct_data())
    
    ct_long <- pivot_longer(ct_data(), cols = -c(Sample, Group), names_to = "Gene", values_to = "Ct")
    
    # Summary stats
    summary_data <- ct_long %>%
      group_by(Gene, Group) %>%
      summarise(mean_Ct = mean(Ct), sd_Ct = sd(Ct), .groups = "drop") %>%
      mutate(label = sprintf("%.2f ± %.2f", mean_Ct, sd_Ct))
    
    # p-values
    p_vals <- ct_long %>%
      group_by(Gene) %>%
      summarise(
        p_value = if (n_distinct(Group) == 2) t.test(Ct ~ Group)$p.value else NA_real_
      )
    
    final <- left_join(summary_data, p_vals, by = "Gene")
    
    # For p-value label positioning
    p_val_data <- ct_long %>%
      group_by(Gene) %>%
      summarise(y_pos = max(Ct, na.rm = TRUE) + 1) %>%
      left_join(p_vals, by = "Gene")
    
    # Generate dynamic color palette
    group_levels <- unique(ct_long$Group)
    palette <- setNames(RColorBrewer::brewer.pal(n = max(3, length(group_levels)), "Set2")[1:length(group_levels)],
                        group_levels)
    
    dodge_width <- 0.6  # Adjust as needed for spacing
    
    ggplot(final, aes(x = Gene, y = mean_Ct, color = Group)) +
      # Raw data points with jitter and dodge
      geom_jitter(
        data = ct_long,
        aes(y = Ct),
        position = position_jitterdodge(
          jitter.width = 0.15,
          dodge.width = dodge_width
        ),
        alpha = 0.5,
        size = 2
      ) +
      # Mean points
      geom_point(
        aes(group = Group),
        position = position_dodge(width = dodge_width),
        size = 4
      ) +
      # Error bars
      geom_errorbar(
        aes(
          ymin = mean_Ct - sd_Ct,
          ymax = mean_Ct + sd_Ct,
          group = Group
        ),
        position = position_dodge(width = dodge_width),
        width = 0.25,
        linewidth = 1
      ) +
      # Mean ± SD labels
      geom_text(
        aes(
          y = mean_Ct + sd_Ct + 0.8,
          label = label,
          group = Group
        ),
        position = position_dodge(width = dodge_width),
        size = 3.5,
        show.legend = FALSE
      ) +
      # p-value labels
      geom_text(
        data = p_val_data,
        aes(x = Gene, y = y_pos, label = ifelse(is.na(p_value), "", sprintf("p=%.3f", p_value))),
        inherit.aes = FALSE,
        size = 3.5,
        color = "black"
      ) +
      scale_color_manual(values = palette) +
      theme_minimal() +
      labs(y = "Ct Value", x = "Gene", title = "Ct Values by Gene and Group") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "top"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))  # Add space for labels
  })
  
  output$download_ct_dist_plot <- downloadHandler(
    filename = function() {
      paste("Ct_Distribution_Plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      req(ct_data())
      
      ct_long <- pivot_longer(ct_data(), cols = -c(Sample, Group), names_to = "Gene", values_to = "Ct")
      
      # Summary stats
      summary_data <- ct_long %>%
        group_by(Gene, Group) %>%
        summarise(mean_Ct = mean(Ct), sd_Ct = sd(Ct), .groups = "drop") %>%
        mutate(label = sprintf("%.2f ± %.2f", mean_Ct, sd_Ct))
      
      # p-values
      p_vals <- ct_long %>%
        group_by(Gene) %>%
        summarise(
          p_value = if (n_distinct(Group) == 2) t.test(Ct ~ Group)$p.value else NA_real_
        )
      
      final <- left_join(summary_data, p_vals, by = "Gene")
      
      p_val_data <- ct_long %>%
        group_by(Gene) %>%
        summarise(y_pos = max(Ct, na.rm = TRUE) + 1) %>%
        left_join(p_vals, by = "Gene")
      
      # Generate dynamic color palette
      group_levels <- unique(ct_long$Group)
      palette <- setNames(RColorBrewer::brewer.pal(n = max(3, length(group_levels)), "Set2")[1:length(group_levels)],
                          group_levels)
      
      dodge_width <- 0.6
      
      p <- ggplot(final, aes(x = Gene, y = mean_Ct, color = Group)) +
        geom_jitter(
          data = ct_long,
          aes(y = Ct),
          position = position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_width),
          alpha = 0.5,
          size = 2
        ) +
        geom_point(
          aes(group = Group),
          position = position_dodge(width = dodge_width),
          size = 4
        ) +
        geom_errorbar(
          aes(ymin = mean_Ct - sd_Ct, ymax = mean_Ct + sd_Ct, group = Group),
          position = position_dodge(width = dodge_width),
          width = 0.25,
          linewidth = 1
        ) +
        geom_text(
          aes(y = mean_Ct + sd_Ct + 0.8, label = label, group = Group),
          position = position_dodge(width = dodge_width),
          size = 3.5
        ) +
        geom_text(
          data = p_val_data,
          aes(x = Gene, y = y_pos, label = ifelse(is.na(p_value), "", sprintf("p=%.3f", p_value))),
          inherit.aes = FALSE,
          size = 3.5,
          color = "black"
        ) +
        scale_color_manual(values = palette) +
        theme_minimal() +
        labs(y = "Ct Value", x = "Gene", title = "Ct Values by Gene and Group") +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          legend.position = "top"
        ) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
      
      ggsave(file, plot = p, width = 10, height = 7, dpi = 300, bg = "white")
    }
  )

  # COMPLETE CUMULATIVE GEOMETRIC MEAN RANKING SCRIPT (FIXED) -------------------
  
  # 1. Stability Method Reactives ------------------------------------------------
  normResults <- reactive({
    req(ct_data())
    nf <- runNormFinder(ct_data())
    df_uno <- nf$UnOrdered
    tibble(
      Gene = rownames(df_uno),
      NF_Stability_Score = df_uno$Stability,
      NF_IGroupSD_v1_Score = df_uno$IGroupSD.V1,
      NF_IGroupSD_v2_Score = df_uno$IGroupSD.V2
    )
  })
  
  bestKeeperResults <- reactive({
    req(ct_data())
    bk <- bestKeeperAnalysis(ct_data())
    stats <- bk$CT.statistics
    corr <- bk$Gene.vs.BKI
    tibble(
      Gene = colnames(stats),
      BestKeeper_SD_Score = as.numeric(stats["SD[+/- CT]", ]),
      BestKeeper_Corr_Score = as.numeric(corr["coeff. of corr. [r]", ]),
      BestKeeper_Corr_PValue = as.numeric(corr["p-value", ])
    )
  })
  
  deltaCtResults <- reactive({
    req(processed_data())
    processed_data()$summary_stats %>%
      select(
        Gene = First_Gene,
        DeltaCtScore = Mean_of_StdDev
      )
  })
  
  geNormResults <- reactive({
    req(geNorm_table())
    geNorm_table() %>%
      select(Gene = Target, GeNormScore = Avg.M)
  })
  
  mixedModelGeneRanks <- reactive({
    req(ranked_mixed_model_results())
    gene_icc_scores <- ranked_mixed_model_results() %>%
      separate_rows(gene.combination, sep = ":") %>%
      group_by(Gene = gene.combination) %>%
      summarise(
        MixedModel_Score = mean(icc.l, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        MixedModel_rank = dense_rank(desc(MixedModel_Score))
      )
        all_genes <- unique(ct_data() %>% select(-Group, -Sample) %>% names())
        missing_genes <- setdiff(all_genes, gene_icc_scores$Gene)
        
        if(length(missing_genes) > 0) {
          gene_icc_scores <- bind_rows(
            gene_icc_scores,
            tibble(Gene = missing_genes, 
                   MixedModel_Score = NA_real_,
                   MixedModel_rank = NA_integer_)
          )
        }
        gene_icc_scores
  })
    
    # 2. Combined Ranking Calculation (FIXED) -------------------------------------
  cumulative_data <- reactive({
    req(
      normResults(), bestKeeperResults(),
      deltaCtResults(),
      geNormResults(), mixedModelGeneRanks()
    )
    
    withProgress(message = "Computing cumulative ranking", value = 0, {
      incProgress(0.05, detail = "Merging method outputs...")
      combined_data <- list(
        normResults(), bestKeeperResults(),
        deltaCtResults(),
        geNormResults(), mixedModelGeneRanks()
      ) %>% 
        reduce(full_join, by = "Gene")
      
      incProgress(0.20, detail = "Converting to numeric and ranking...")
      Sys.sleep(0.01)
      
      combined_data <- combined_data %>%
        mutate(across(ends_with("_Score"), as.numeric))
      
      incProgress(0.50, detail = "Computing per-method ranks and geometric mean...")
      Sys.sleep(0.01)
      
      # perform ranks and geomean as before
      combined_data <- combined_data %>%
        mutate(
          NF_Stab_rank = dense_rank(NF_Stability_Score),
          NF_IGroupSDv1_rank = dense_rank(NF_IGroupSD_v1_Score),
          NF_IGroupSDv2_rank = dense_rank(NF_IGroupSD_v2_Score),
          BestKeeper_SD_rank = dense_rank(BestKeeper_SD_Score),
          BestKeeper_Corr_rank = if_else(
            BestKeeper_Corr_PValue > 0.05,
            NA_real_,
            dense_rank(desc(BestKeeper_Corr_Score))
          ),
          DeltaCt_rank = dense_rank(DeltaCtScore),
          GeNorm_rank = dense_rank(GeNormScore),
          MixedModel_rank = dense_rank(desc(MixedModel_Score))
        ) %>%
        rowwise() %>%
        mutate(
          log_ranks = list(log(c(
            NF_Stab_rank, NF_IGroupSDv1_rank, NF_IGroupSDv2_rank,
            BestKeeper_SD_rank, BestKeeper_Corr_rank,
            DeltaCt_rank, GeNorm_rank,
            MixedModel_rank
          ))),
          valid_logs = sum(!is.na(log_ranks)),
          Cumulative_Geomean_Rank = ifelse(
            valid_logs > 0,
            exp(mean(log_ranks, na.rm = TRUE)),
            NA_real_
          )
        ) %>%
        ungroup() %>%
        arrange(Cumulative_Geomean_Rank) %>%
        select(-log_ranks, -valid_logs)
      
      incProgress(0.25, detail = "Done")
      combined_data
    })
  })
  
    
    # 3. Output Rendering (FIXED) -------------------------------------------------
    output$cumulative_table <- renderDT({
      req(cumulative_data())
      cumulative_data() %>%
        select(
          Gene,
          NF_Stability_Score, NF_Stab_rank,
          NF_IGroupSD_v1_Score, NF_IGroupSDv1_rank,
          NF_IGroupSD_v2_Score, NF_IGroupSDv2_rank,
          BestKeeper_SD_Score, BestKeeper_SD_rank,
          BestKeeper_Corr_Score, BestKeeper_Corr_PValue, BestKeeper_Corr_rank,
          DeltaCtScore, DeltaCt_rank,
          GeNormScore, GeNorm_rank,
          MixedModel_Score, MixedModel_rank,
          Cumulative_Geomean_Rank
        ) %>%
        rename(
          "NF Stability" = NF_Stability_Score,
          "NF Rank" = NF_Stab_rank,
          "NF IGroupSDv1" = NF_IGroupSD_v1_Score,
          "NF v1 Rank" = NF_IGroupSDv1_rank,
          "NF IGroupSDv2" = NF_IGroupSD_v2_Score,
          "NF v2 Rank" = NF_IGroupSDv2_rank,
          "BK SD" = BestKeeper_SD_Score,
          "BK SD Rank" = BestKeeper_SD_rank,
          "BK Corr (r)" = BestKeeper_Corr_Score,
          "BK p-value" = BestKeeper_Corr_PValue,
          "BK Corr Rank" = BestKeeper_Corr_rank,
          "ΔCt Score" = DeltaCtScore,
          "ΔCt Rank" = DeltaCt_rank,
          "geNorm M" = GeNormScore,
          "geNorm Rank" = GeNorm_rank,
          "Mixed Model ICC" = MixedModel_Score,
          "MM Rank" = MixedModel_rank,
          "Final Rank" = Cumulative_Geomean_Rank
        )
    }, options = list(
      pageLength = 10,
      scrollX = TRUE,
      autoWidth = TRUE,
      columnDefs = list(list(targets = "_all", className = "dt-right")), 
      rownames = FALSE))
    
    # NEW: Cumulative ranking plot
    output$cumulative_plot <- renderPlot({
      req(cumulative_data())
      topn <- cumulative_data() %>% arrange(Cumulative_Geomean_Rank) %>% slice_head(n = 20)
      
      # long per-method ranks for overlay
      overlay <- topn %>%
        select(Gene,
               NF_Stab_rank, NF_IGroupSDv1_rank, NF_IGroupSDv2_rank,
               BestKeeper_SD_rank, BestKeeper_Corr_rank,
               DeltaCt_rank, GeNorm_rank, MixedModel_rank,
               Cumulative_Geomean_Rank) %>%
        pivot_longer(cols = -c(Gene, Cumulative_Geomean_Rank),
                     names_to = "Method", values_to = "Rank")
      
      topn <- topn %>%
        mutate(Gene = fct_reorder(as.factor(Gene), Cumulative_Geomean_Rank),
               Gene = fct_rev(Gene))
      
      overlay$Gene <- factor(overlay$Gene, levels = levels(topn$Gene))
      
      ggplot() +
        geom_col(data = topn, aes(x = Gene, y = Cumulative_Geomean_Rank), fill = "grey80") +
        geom_point(data = overlay, aes(x = Gene, y = Rank, color = Method), position = position_jitter(width = 0.12, height = 0), size = 2) +
        coord_flip() +
        labs(title = "Top Stable Genes — Combined (bar = GM, points = per-method ranks)",
             y = "Rank (lower = more stable)",
             x = "Gene (top = most stable → bottom = least stable)",
             color = "Method") +
        theme_minimal() +
        theme(text = element_text(size = 11))
    })
    
    # FIXED: Download handler for detailed table
    output$download_cumulative_plot <- downloadHandler(
      filename = function() {
        paste0("Cumulative_Ranking_Plot_", Sys.Date(), ".png")
      },
      contentType = "image/png", # ensure browser treats it as PNG
      content = function(file) {
        req(cumulative_data())
        
        
        # prepare top-20 data and overlay (same semantics as your UI plot)
        topn <- cumulative_data() %>%
          dplyr::arrange(Cumulative_Geomean_Rank) %>%
          dplyr::slice_head(n = 20)
        
        
        overlay <- topn %>%
          dplyr::select(Gene,
                        NF_Stab_rank, NF_IGroupSDv1_rank, NF_IGroupSDv2_rank,
                        BestKeeper_SD_rank, BestKeeper_Corr_rank,
                        DeltaCt_rank, GeNorm_rank, MixedModel_rank,
                        Cumulative_Geomean_Rank) %>%
          tidyr::pivot_longer(cols = -c(Gene, Cumulative_Geomean_Rank),
                              names_to = "Method", values_to = "Rank")
        
        
        # ensure ordering so top = most stable
        topn <- topn %>%
          dplyr::mutate(
            Gene = forcats::fct_reorder(as.factor(Gene), Cumulative_Geomean_Rank),
            Gene = forcats::fct_rev(Gene)
          )
        
        
        # align overlay factor levels with topn
        overlay$Gene <- factor(as.character(overlay$Gene), levels = levels(topn$Gene))
        
        
        # build plot (same as your UI plot)
        p <- ggplot2::ggplot() +
          ggplot2::geom_col(data = topn, aes(x = Gene, y = Cumulative_Geomean_Rank), fill = "grey80") +
          ggplot2::geom_point(data = overlay, aes(x = Gene, y = Rank, color = Method),
                              position = ggplot2::position_jitter(width = 0.12, height = 0),
                              size = 2) +
          ggplot2::coord_flip() +
          ggplot2::labs(title = "Top Stable Genes — Combined (bar = GM, points = per-method ranks)",
                        y = "Rank (lower = more stable)",
                        x = "Gene (top = most stable → bottom = least stable)",
                        color = "Method") +
          ggplot2::theme_minimal() +
          ggplot2::theme(text = ggplot2::element_text(size = 11))
        
        
        # Write PNG reliably to the path 'file' that Shiny provides.
        # Use the PNG device so the file content and MIME type are unambiguous.
        png(filename = file, width = 10 * 300, height = 8 * 300, res = 300)
        print(p)
        dev.off()
      }
    )
  
    # FIXED: Download handler for detailed table
    output$download_cumulative_table <- downloadHandler(
      filename = function() {
        paste0("Detailed_Gene_Stability_Table_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(cumulative_data())
        write.csv(cumulative_data(), file, row.names = FALSE)
      }
    )
    
    # NEW: Download handler for raw ranked data
    output$download_cumulative_raw <- downloadHandler(
      filename = function() {
        paste0("Raw_Cumulative_Ranking_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(cumulative_data())
        cumulative_data() %>%
          select(Gene, Cumulative_Geomean_Rank) %>%
          write.csv(file, row.names = FALSE)
      }
    )
  
} 

# Run the application 
shinyApp(ui = ui, server = server)

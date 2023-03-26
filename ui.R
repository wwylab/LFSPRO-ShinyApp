library(shiny)

shinyUI(
  fluidPage(
    includeCSS("style.css"),
    # Application title
    titlePanel("Predicting mutation carrier status and future cancer risk for families with Li-Fraumeni syndrome"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(
        p(em("TP53"), "germline mutations are the main cause of Li-Fraumeni Syndrome. 
          This application is designed to estimate probabilities that: 
          1) the counselee is a ", em("TP53"), "germline mutation carrier, 
          2) the counselee develops breast cancer, sarcoma or any other cancers for their fist primary cancer diagnosis, 
          3) the counselee develops a second primary cancer in future, on the basis of his/her family cancer history. 
          The package also provides functions for identifying an LFS individual based on the Classic LFS and Chompret criteria."),
        br(),
        
        p("The online tool requires two files as input: Family Data File and Cancer Data File.
  Please check the two links below to download the example files:"),
        p(a("Family Data File", href="fam.data.csv", download="fam.data.csv"), br(),
          a("Cancer Data File", href="cancer.data.csv", download="cancer.data.csv")),
        
        p("The cancer.type column in the Cancer Data File are based on the LFS spectrum of cancers which includes adrenal cortical carcinomas (acc), brain tumors (brain), breast cancer (breast), choroid plexus carcinomas (choroid), leukemia, lung cancer (lung), osteosarcomas (ost), and soft tissue sarcomas (sts). 
          All other malignant cancer diagnoses are considered non-LFS spectrum cancer types (non.lfs) while non-malignant are labeled as benign."),
        
        br(),
        p("After running LFSPRO, there will be a figure and a table shown on the right. 
  The figure shows the pedigree structure of the input family.",
          span("The cancer history will be shown by clicking '+' in the table.", style = "color:blue")),
        
        br(),
        fileInput("file1", "Choose Family Data File",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        
        fileInput("file2", "Choose Cancer Data File",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        
        tags$div(
          title = "ID of samples selected for calculation. Leave it blank to select all samples. For example: ID1,ID3,ID4",
          textInput("txtSampleId", label = h5("Sample ID"), value = "")
        ),
        
        uiOutput('ui.action'),
        
        hr(),
        uiOutput("ui.mutation"),
        br(),
        uiOutput("ui.model"),
        br(),
        uiOutput("ui.cutoff")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        fluidRow(plotOutput("distPlot")),
        span(textOutput("warning"), style="color:red"),
        br(),
        fluidRow(verbatimTextOutput("sInfo")),
        fluidRow(DT::dataTableOutput("table")),
        fluidRow(column(12, align="center",
                        uiOutput('ui.download')))
      )
    )
  )
)

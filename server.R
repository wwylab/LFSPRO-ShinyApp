library(shiny)
library(kinship2)
library(DT)
library(data.table)
library(shinyjs)
library(shinyalert)
library(reshape2)
library(ggplot2)
library(plotly)
library(ggsci)
library(LFSPRO)

source("functions.R")
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

#Style
sketch = htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th(rowspan = 3, ' '),
      th(rowspan = 3, 'Details'),
      th(rowspan = 3, 'ID'),
      th(rowspan = 3, "Info"),
      th(rowspan = 3, 'Gene Testing'),
      th(rowspan = 3, 'ProbLFSPRO'),
      th(rowspan = 3, 'LFSPRO-carrier'),
      th(rowspan = 3, 'Chompret criteria'),
      th(rowspan = 3, 'Classic criteria'),
      th(colspan = 13, 'Cancer Risk')
    ),
    tr(
      th(colspan = 3, "Breast Cancer"),
      th(colspan = 3, "Sarcoma"),
      th(colspan = 3, "Other Cancers"),
      th(colspan = 3, "Second Primary Cancer"),
      th(rowspan = 2, "Figure")
    ),
    tr(
      lapply(rep(c('5 years', '10 years', '15 years'), 4), th)
    )
  )
))

cutoff <- 0.2
use.mutation <- TRUE
use.CS <- FALSE

shinyServer(function(input, output) {
  LFSPRO.rlt <- NULL
  myValue <- reactiveValues(idx.button = '')
  
  buttonInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  
  famdata <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath, header = TRUE, sep = ",", stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM")
  })
  
  cancerdata <- reactive({
    infile <- input$file2
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath, header = TRUE, sep = ",", stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM")
  })
  
  cid <- reactive({
    cid.tmp <- input$txtSampleId
    cid <- NULL
    if(length(cid.tmp)>0){
      cid <- strsplit(cid.tmp, ";")[[1]]
    }
    cid
  })
  
  output$ui.action <- renderUI({
    if (is.null(famdata())) return()
    if (is.null(cancerdata())) return()
    actionButton("action", "Run LFSPRO")
  })
  
  observeEvent(
    eventExpr = input$action,
    handlerExpr = {
      output$distPlot <- renderPlot({
        if (is.null(input$action)) return()
        if (input$action==0) return()
        isolate({
          fam.data <- famdata()
          idx <- which(fam.data$proband == "Y")
          fam.data <- rbind(fam.data[-idx,], fam.data[idx,])
          if (is.null(fam.data)) return(NULL)
          cancer.data <- cancerdata()
          if (is.null(cancer.data)) return(NULL)
          fam.data <- fam.update(fam.data, cancer.data)
          fam.data <- dummy.create(fam.data.process(fam.data, cancer.data))
          cancer.data <- cancer.data.process(fam.data, cancer.data)
          
          aff <- data.frame(Affected = fam.data$id %in% cancer.data$id)
          
          ped <- pedigree(id =  fam.data$id, 
                          dadid = fam.data$fid,
                          momid = fam.data$mid,
                          sex = ifelse(fam.data$gender==0, 2, 1),
                          famid = rep("fam", nrow(fam.data)),
                          status = ifelse(fam.data$vital == "A", 0, 1),
                          affected = as.matrix(aff))
          pedSel <- ped["fam"]
          id2 <- pedSel$id
          id2[fam.data$proband == "Y"] <- paste(id2[fam.data$proband == "Y"], "Proband", sep = "\n")
          mut <- !is.na(fam.data$test) & (fam.data$test == 1)
          id2[mut] <- paste(id2[mut], "Mut", sep = "\n")
          wt <- !is.na(fam.data$test) & (fam.data$test == 0)
          id2[wt] <- paste(id2[wt], "WT", sep = "\n")
          col.label <- c("Unaffected", "Proband", "Dummy")
          names(col.label) <- c("black", "red", "blue")
          legendPlot(pedSel, affected.label = c("Affected"),
                     col = ifelse(fam.data$dummy == 0, ifelse(fam.data$proband == "Y", "red", "black"), "blue"),
                     col.label = col.label, id = id2, symbolsize = 1.25, cex = 1)
        })
      })
      
      output$warning <- renderText({
        if (is.null(input$action)) return()
        if (input$action==0) return()
        isolate({
          fam.data <- famdata()
          if (is.null(fam.data)) return(NULL)
          cancer.data <- cancerdata()
          if (is.null(cancer.data)) return(NULL)
          non.informative <- (nrow(cancer.data) <= 1)
          if (non.informative) {
            warning <- "Warning: This family has limited information; The predicted risks can be unreliable."
          } else {
            warning <- ""
          }
          warning
        })
      })
      
      output$table <- DT::renderDataTable({
        if (is.null(input$action)) return(DT::datatable(NULL))
        if (input$action==0) return(DT::datatable(NULL))
        
        DT::datatable({
          isolate({
            fam.data <- famdata()
            if (is.null(fam.data)) return(NULL)
            cancer.data <- cancerdata()
            if (is.null(cancer.data)) return(NULL)
            fam.data <- fam.update(fam.data, cancer.data)
            fam.data <- dummy.create(fam.data.process(fam.data, cancer.data))
            cancer.data <- cancer.data.process(fam.data, cancer.data)
            
            fam.data$fam.id <- "fam"
            cancer.data$fam.id <- "fam"
            
            cid <- cid()
            
            if(is.null(cid) || length(cid) == 0){
              counselee.id <- data.frame(fam.id = fam.data$fam.id, id = fam.data$id)
            } else {
              idx.cid <- cid %in% fam.data$id
              if(sum(!idx.cid) > 0){
                cid.notfound <- cid[!idx.cid]
                cid <- cid[idx.cid]
                msg <- paste0("The following id(s) are not found in the input data: ",
                              paste(cid.notfound, collapse = ", "), ".")
                if(sum(idx.cid)==0){
                  msg <- paste0(msg, " All samples in the input file are selected for the calculation. ")
                } else {
                  msg <- paste0(msg, " Only the following samples are selected for the calculation: ",
                                paste(cid, collapse = ", "), ".")
                }
                shinyalert("ID Not Found", msg, type = "warning")
              }
              counselee.id <- data.frame(fam.id = "fam", id = cid)
            }
            
            ## deceased patients are marked
            id.dead <- fam.data$id[fam.data$vital == "D"]
            counselee.dead <- counselee.id[counselee.id$id %in% id.dead,]
            counselee.alive <- counselee.id[!counselee.id$id %in% id.dead,]
            idx.dead <- counselee.id$id %in% id.dead
            
            info <- NULL
            for(i in 1:length(counselee.id$id)){
              id <- counselee.id$id[i]
              idx.can <- which(cancer.data$id == id)
              idx.fam <- which(fam.data$id == id) 
              if (fam.data$dummy[idx.fam] == 0) {
                info.tmp <- ""
              } else {
                info.tmp <- "DUMMY; "
              }
              if (fam.data$vital[idx.fam] == 'D')
              {
                info.tmp <- paste0(info.tmp,'Deceased; Was ')
              }
              info.tmp <- paste0(info.tmp, fam.data$age[idx.fam], " years old ")
              info.tmp <- paste0(info.tmp, ifelse(fam.data$gender[idx.fam] == 0, "female; ", "male; "))
              if (!is.na(fam.data$test[idx.fam])){
                info.tmp <- paste0(info.tmp, ifelse(fam.data$gender[idx.fam] == 0, "Her ", "His "))
                info.tmp <- paste0(info.tmp, 'confirmed genetic testing result is ', fam.data$test[idx.fam], '; ')
              }
              if(length(idx.can)){
                for(j in idx.can){
                  info.tmp <- paste0(info.tmp, cancer.data$cancer.type[j], " at age ", cancer.data$diag.age[j], "; ")
                }
              } else {
                info.tmp <- paste0(info.tmp, "No Cancer. ")
              }
              info <- c(info, info.tmp)
            }
            
            # run LFSPRO on all selected patients
            temp.fam.data <- fam.data
            temp.fam.data$vital <- 'A'
            
            rltTmp <- runLFSPRO(temp.fam.data, cancer.data, counselee.id, mut.info = use.mutation)
            rltTmp[,7:ncol(rltTmp)][idx.dead,] <- NA
            rltTmp$info <- info
            rltTmp <- merge(rltTmp, subset(fam.data, select = c(id, vital, age)))
            LFSPRO.rlt <<- rltTmp
            
            if (use.CS == TRUE) {
              ProbLFSPRO = LFSPRO.rlt$carrier.cs
              LFSPRO = factor(ifelse(LFSPRO.rlt$carrier.cs>cutoff, "Yes", "No"), 
                              levels = c("Yes", "No"))
            } else {
              ProbLFSPRO = LFSPRO.rlt$carrier.mpc
              LFSPRO = factor(ifelse(LFSPRO.rlt$carrier.mpc>cutoff, "Yes", "No"), 
                              levels = c("Yes", "No"))
            }
            
            rlt <- data.frame(
              id = factor(LFSPRO.rlt$id, levels = LFSPRO.rlt$id),
              info = LFSPRO.rlt$info,
              test = fam.genelevel(LFSPRO.rlt$test),
              ProbLFSPRO = ProbLFSPRO,
              LFSPRO = LFSPRO,
              Chompret = factor(ifelse(LFSPRO.rlt$chompret, "Yes", "No"),
                                levels = c("Yes", "No")),
              Classic = factor(ifelse(LFSPRO.rlt$classic, "Yes", "No"),
                               levels = c("Yes", "No")),
              breast.5 = LFSPRO.rlt$breast.5,
              breast.10 = LFSPRO.rlt$breast.10,
              breast.15 = LFSPRO.rlt$breast.15,
              sarcoma.5 = LFSPRO.rlt$sarcoma.5,
              sarcoma.10 = LFSPRO.rlt$sarcoma.10,
              sarcoma.15 = LFSPRO.rlt$sarcoma.15,
              other.5 = LFSPRO.rlt$other.5,
              other.10 = LFSPRO.rlt$other.10,
              other.15 = LFSPRO.rlt$other.15,
              second.5 = LFSPRO.rlt$second.5,
              second.10 = LFSPRO.rlt$second.10,
              second.15 = LFSPRO.rlt$second.15,
              figure = buttonInput(
                FUN = actionButton,
                len = nrow(LFSPRO.rlt),
                id = "button_",
                label = "Risk Trend",
                onclick = 'Shiny.onInputChange(\"lastClick\",  this.id)'
              )
            )            
            rlt <- cbind('Details' = '&oplus;', rlt)
            rlt
          })
        }, container = sketch, filter = 'top', escape = FALSE,
        options = list(
          columnDefs = list(
            list(visible = FALSE, targets = c(0, 3)),
            list(orderable = FALSE, className = 'details-control', targets = 1)
          )
        ),
        callback = JS(
          "
        table.column(1).nodes().to$().css({cursor: 'pointer'});
        
        var format = function(d) {
          return '<div style=\"background-color:#eee; padding: .5em;\"> ' +
            d[3] + '</div>';
        };
        
        table.on('click', 'td.details-control', function() {
            var td = $(this), row = table.row(td.closest('tr'));
            if (row.child.isShown()) {
            row.child.hide();
            td.html('&oplus;');
            } else {
            row.child(format(row.data())).show();
            td.html('&CircleMinus;');
            }
        });
        "
        )) %>% 
          DT::formatRound('ProbLFSPRO', 2) %>%
          DT::formatRound(
            c("breast.5", "breast.10", "breast.15", "sarcoma.5", "sarcoma.10", "sarcoma.15", 
              "other.5", "other.10", "other.15", "second.5", "second.10", "second.15"), 2) %>%
          formatStyle(c('LFSPRO', 'Chompret', 'Classic'), color = styleEqual("Yes", 'red'))
      })
    }
  )
  
  output$ui.cutoff <- renderUI({
    if (is.null(input$action) ) return()
    if (input$action==0) return()
    
    sliderInput("cutoff", "Cutoff for probability",
                min=0, max=1, value = 0.2, step = 0.05)
  })
  
  output$ui.mutation <- renderUI({
    if (is.null(input$action) ) return()
    if (input$action==0) return()
    
    radioButtons("mutation", "Use known testing results?", 
                 choices = c("Yes","No"), selected = "Yes", inline = TRUE)
  })
  
  output$ui.model <- renderUI({
    if (is.null(input$action) ) return()
    if (input$action==0) return()
    
    radioButtons("model", "Which model do you want to use to calculate carrier probability?", 
                 choices = c("MPC","CS"), selected = "MPC", inline = TRUE)
  })
  
  output$download <- downloadHandler('LFSPRO.csv', content = function(file) {
    write.csv(LFSPRO.rlt, file, row.names = FALSE)
  })
  
  output$ui.download <- renderUI({
    if (is.null(input$action)) return()
    if (input$action==0) return()
    downloadButton('download', "Download Results")
  })
  
  observeEvent(eventExpr = input$cutoff, handlerExpr = {
    cutoff <<- input$cutoff
    
    output$table <- DT::renderDataTable({
      if (is.null(input$action)) return(DT::datatable(NULL))
      if (input$action==0) return(DT::datatable(NULL))
      
      if (use.CS == TRUE) {
        ProbLFSPRO = LFSPRO.rlt$carrier.cs
        LFSPRO = factor(ifelse(LFSPRO.rlt$carrier.cs>cutoff, "Yes", "No"), 
                        levels = c("Yes", "No"))
      } else {
        ProbLFSPRO = LFSPRO.rlt$carrier.mpc
        LFSPRO = factor(ifelse(LFSPRO.rlt$carrier.mpc>cutoff, "Yes", "No"), 
                        levels = c("Yes", "No"))
      }
      
      DT::datatable({
        isolate({
          rlt <- data.frame(
            id = factor(LFSPRO.rlt$id, levels =  LFSPRO.rlt$id),
            info = LFSPRO.rlt$info,
            test = fam.genelevel(LFSPRO.rlt$test),
            ProbLFSPRO = ProbLFSPRO,
            LFSPRO = LFSPRO,
            Chompret = factor(ifelse(LFSPRO.rlt$chompret, "Yes", "No"),
                              levels = c("Yes", "No")),
            Classic = factor(ifelse(LFSPRO.rlt$classic, "Yes", "No"),
                             levels = c("Yes", "No")),
            breast.5 = LFSPRO.rlt$breast.5,
            breast.10 = LFSPRO.rlt$breast.10,
            breast.15 = LFSPRO.rlt$breast.15,
            sarcoma.5 = LFSPRO.rlt$sarcoma.5,
            sarcoma.10 = LFSPRO.rlt$sarcoma.10,
            sarcoma.15 = LFSPRO.rlt$sarcoma.15,
            other.5 = LFSPRO.rlt$other.5,
            other.10 = LFSPRO.rlt$other.10,
            other.15 = LFSPRO.rlt$other.15,
            second.5 = LFSPRO.rlt$second.5,
            second.10 = LFSPRO.rlt$second.10,
            second.15 = LFSPRO.rlt$second.15,
            figure = buttonInput(
              FUN = actionButton,
              len = nrow(LFSPRO.rlt),
              id = "button_",
              label = "Risk Trend",
              onclick = 'Shiny.onInputChange(\"lastClick\",  this.id)'
            )
          )
          
          rlt <- cbind('Details' = '&oplus;', rlt)
          rlt
        })
      }, container = sketch, filter = 'top', escape = FALSE,
      options = list(
        columnDefs = list(
          list(visible = FALSE, targets = c(0, 3)),
          list(orderable = FALSE, className = 'details-control', targets = 1)
        )
      ),
      callback = JS(
        "
        table.column(1).nodes().to$().css({cursor: 'pointer'});
        
        var format = function(d) {
          return '<div style=\"background-color:#eee; padding: .5em;\"> ' +
            d[3] + '</div>';
        };
        
        table.on('click', 'td.details-control', function() {
            var td = $(this), row = table.row(td.closest('tr'));
            if (row.child.isShown()) {
            row.child.hide();
            td.html('&oplus;');
            } else {
            row.child(format(row.data())).show();
            td.html('&CircleMinus;');
            }
        });
        "
      )) %>% 
        DT::formatRound('ProbLFSPRO', 2) %>%
        DT::formatRound(
          c("breast.5", "breast.10", "breast.15", "sarcoma.5", "sarcoma.10", "sarcoma.15", 
            "other.5", "other.10", "other.15", "second.5", "second.10", "second.15"), 2) %>%
        formatStyle(c('LFSPRO', 'Chompret', 'Classic'), color = styleEqual("Yes", 'red'))
    })
  })
  
  observeEvent(eventExpr = input$mutation, handlerExpr = {
    use.mutation <<- (input$mutation == "Yes")
    
    output$table <- DT::renderDataTable({
      if (is.null(input$action)) return(DT::datatable(NULL))
      if (input$action==0) return(DT::datatable(NULL))
      
      DT::datatable({
        isolate({
          fam.data <- famdata()
          if (is.null(fam.data)) return(NULL)
          cancer.data <- cancerdata()
          if (is.null(cancer.data)) return(NULL)
          fam.data <- fam.update(fam.data, cancer.data)
          fam.data <- dummy.create(fam.data.process(fam.data, cancer.data))
          cancer.data <- cancer.data.process(fam.data, cancer.data)
          
          fam.data$fam.id <- "fam"
          cancer.data$fam.id <- "fam"
          
          cid <- cid()
          
          if(is.null(cid) || length(cid) == 0){
            counselee.id <- data.frame(fam.id = fam.data$fam.id, id = fam.data$id)
          } else {
            idx.cid <- cid %in% fam.data$id
            if(sum(!idx.cid) > 0){
              cid.notfound <- cid[!idx.cid]
              cid <- cid[idx.cid]
              msg <- paste0("The following id(s) are not found in the input data: ",
                            paste(cid.notfound, collapse = ", "), ".")
              if(sum(idx.cid)==0){
                msg <- paste0(msg, " All samples in the input file are selected for the calculation. ")
              } else {
                msg <- paste0(msg, " Only the following samples are selected for the calculation: ",
                              paste(cid, collapse = ", "), ".")
              }
              shinyalert("ID Not Found", msg, type = "warning")
            }
            counselee.id <- data.frame(fam.id = "fam", id = cid)
          }
          
          ## deceased patients are marked
          id.dead <- fam.data$id[fam.data$vital == "D"]
          counselee.dead <- counselee.id[counselee.id$id %in% id.dead,]
          counselee.alive <- counselee.id[!counselee.id$id %in% id.dead,]
          idx.dead <- counselee.id$id %in% id.dead
          
          info <- NULL
          for(i in 1:length(counselee.id$id)){
            id <- counselee.id$id[i]
            idx.can <- which(cancer.data$id == id)
            idx.fam <- which(fam.data$id == id) 
            if (fam.data$dummy[idx.fam] == 0) {
              info.tmp <- ""
            } else {
              info.tmp <- "DUMMY; "
            }
            if (fam.data$vital[idx.fam] == 'D')
            {
              info.tmp <- paste0(info.tmp,'Deceased; Was ')
            }
            info.tmp <- paste0(info.tmp, fam.data$age[idx.fam], " years old ")
            info.tmp <- paste0(info.tmp, ifelse(fam.data$gender[idx.fam] == 0, "female; ", "male; "))
            if (!is.na(fam.data$test[idx.fam])){
              info.tmp <- paste0(info.tmp, ifelse(fam.data$gender[idx.fam] == 0, "Her ", "His "))
              info.tmp <- paste0(info.tmp, 'confirmed genetic testing result is ', fam.data$test[idx.fam], '; ')
            }
            if(length(idx.can)){
              for(j in idx.can){
                info.tmp <- paste0(info.tmp, cancer.data$cancer.type[j], " at age ", cancer.data$diag.age[j], "; ")
              }
            } else {
              info.tmp <- paste0(info.tmp, "No Cancer. ")
            }
            info <- c(info, info.tmp)
          }
          
          # run LFSPRO on all selected patients
          temp.fam.data <- fam.data
          temp.fam.data$vital <- 'A'
          
          rltTmp <- runLFSPRO(temp.fam.data, cancer.data, counselee.id, mut.info = use.mutation)
          rltTmp[,7:ncol(rltTmp)][idx.dead,] <- NA
          rltTmp$info <- info
          rltTmp <- merge(rltTmp, subset(fam.data, select = c(id, vital, age)))
          LFSPRO.rlt <<- rltTmp
          
          if (use.CS == TRUE) {
            ProbLFSPRO = LFSPRO.rlt$carrier.cs
            LFSPRO = factor(ifelse(LFSPRO.rlt$carrier.cs>cutoff, "Yes", "No"), 
                            levels = c("Yes", "No"))
          } else {
            ProbLFSPRO = LFSPRO.rlt$carrier.mpc
            LFSPRO = factor(ifelse(LFSPRO.rlt$carrier.mpc>cutoff, "Yes", "No"), 
                            levels = c("Yes", "No"))
          }
          
          rlt <- data.frame(
            id = factor(LFSPRO.rlt$id, levels =  LFSPRO.rlt$id),
            info = LFSPRO.rlt$info,
            test = fam.genelevel(LFSPRO.rlt$test),
            ProbLFSPRO = ProbLFSPRO,
            LFSPRO = LFSPRO,
            Chompret = factor(ifelse(LFSPRO.rlt$chompret, "Yes", "No"),
                              levels = c("Yes", "No")),
            Classic = factor(ifelse(LFSPRO.rlt$classic, "Yes", "No"),
                             levels = c("Yes", "No")),
            breast.5 = LFSPRO.rlt$breast.5,
            breast.10 = LFSPRO.rlt$breast.10,
            breast.15 = LFSPRO.rlt$breast.15,
            sarcoma.5 = LFSPRO.rlt$sarcoma.5,
            sarcoma.10 = LFSPRO.rlt$sarcoma.10,
            sarcoma.15 = LFSPRO.rlt$sarcoma.15,
            other.5 = LFSPRO.rlt$other.5,
            other.10 = LFSPRO.rlt$other.10,
            other.15 = LFSPRO.rlt$other.15,
            second.5 = LFSPRO.rlt$second.5,
            second.10 = LFSPRO.rlt$second.10,
            second.15 = LFSPRO.rlt$second.15,
            figure = buttonInput(
              FUN = actionButton,
              len = nrow(LFSPRO.rlt),
              id = "button_",
              label = "Risk Trend",
              onclick = 'Shiny.onInputChange(\"lastClick\",  this.id)'
            )
          )
          
          rlt <- cbind('Details' = '&oplus;', rlt)
          rlt
        })
      }, container = sketch, filter = 'top', escape = FALSE,
      options = list(
        columnDefs = list(
          list(visible = FALSE, targets = c(0, 3)),
          list(orderable = FALSE, className = 'details-control', targets = 1)
        )
      ),
      callback = JS(
        "
        table.column(1).nodes().to$().css({cursor: 'pointer'});
        
        var format = function(d) {
          return '<div style=\"background-color:#eee; padding: .5em;\"> ' +
            d[3] + '</div>';
        };
        
        table.on('click', 'td.details-control', function() {
            var td = $(this), row = table.row(td.closest('tr'));
            if (row.child.isShown()) {
            row.child.hide();
            td.html('&oplus;');
            } else {
            row.child(format(row.data())).show();
            td.html('&CircleMinus;');
            }
        });
        "
      )) %>% 
        DT::formatRound('ProbLFSPRO', 2) %>%
        DT::formatRound(
          c("breast.5", "breast.10", "breast.15", "sarcoma.5", "sarcoma.10", "sarcoma.15", 
            "other.5", "other.10", "other.15", "second.5", "second.10", "second.15"), 2) %>%
        formatStyle(c('LFSPRO', 'Chompret', 'Classic'), color = styleEqual("Yes", 'red'))
    })
  })
  
  observeEvent(eventExpr = input$model, handlerExpr = {
    use.CS <<- (input$model == "CS")
    
    output$table <- DT::renderDataTable({
      if (is.null(input$action)) return(DT::datatable(NULL))
      if (input$action==0) return(DT::datatable(NULL))
      
      DT::datatable({
        isolate({
          fam.data <- famdata()
          if (is.null(fam.data)) return(NULL)
          cancer.data <- cancerdata()
          if (is.null(cancer.data)) return(NULL)
          fam.data <- fam.update(fam.data, cancer.data)
          fam.data <- dummy.create(fam.data.process(fam.data, cancer.data))
          cancer.data <- cancer.data.process(fam.data, cancer.data)
          
          fam.data$fam.id <- "fam"
          cancer.data$fam.id <- "fam"
          
          cid <- cid()
          
          if(is.null(cid) || length(cid) == 0){
            counselee.id <- data.frame(fam.id = fam.data$fam.id, id = fam.data$id)
          } else {
            idx.cid <- cid %in% fam.data$id
            if(sum(!idx.cid) > 0){
              cid.notfound <- cid[!idx.cid]
              cid <- cid[idx.cid]
              msg <- paste0("The following id(s) are not found in the input data: ",
                            paste(cid.notfound, collapse = ", "), ".")
              if(sum(idx.cid)==0){
                msg <- paste0(msg, " All samples in the input file are selected for the calculation. ")
              } else {
                msg <- paste0(msg, " Only the following samples are selected for the calculation: ",
                              paste(cid, collapse = ", "), ".")
              }
              shinyalert("ID Not Found", msg, type = "warning")
            }
            counselee.id <- data.frame(fam.id = "fam", id = cid)
          }
          
          ## deceased patients are marked
          id.dead <- fam.data$id[fam.data$vital == "D"]
          counselee.dead <- counselee.id[counselee.id$id %in% id.dead,]
          counselee.alive <- counselee.id[!counselee.id$id %in% id.dead,]
          idx.dead <- counselee.id$id %in% id.dead
          
          info <- NULL
          for(i in 1:length(counselee.id$id)){
            id <- counselee.id$id[i]
            idx.can <- which(cancer.data$id == id)
            idx.fam <- which(fam.data$id == id) 
            if (fam.data$dummy[idx.fam] == 0) {
              info.tmp <- ""
            } else {
              info.tmp <- "DUMMY; "
            }
            if (fam.data$vital[idx.fam] == 'D')
            {
              info.tmp <- paste0(info.tmp,'Deceased; Was ')
            }
            info.tmp <- paste0(info.tmp, fam.data$age[idx.fam], " years old ")
            info.tmp <- paste0(info.tmp, ifelse(fam.data$gender[idx.fam] == 0, "female; ", "male; "))
            if (!is.na(fam.data$test[idx.fam])){
              info.tmp <- paste0(info.tmp, ifelse(fam.data$gender[idx.fam] == 0, "Her ", "His "))
              info.tmp <- paste0(info.tmp, 'confirmed genetic testing result is ', fam.data$test[idx.fam], '; ')
            }
            if(length(idx.can)){
              for(j in idx.can){
                info.tmp <- paste0(info.tmp, cancer.data$cancer.type[j], " at age ", cancer.data$diag.age[j], "; ")
              }
            } else {
              info.tmp <- paste0(info.tmp, "No Cancer. ")
            }
            info <- c(info, info.tmp)
          }
          
          # run LFSPRO on all selected patients
          temp.fam.data <- fam.data
          temp.fam.data$vital <- 'A'
          
          rltTmp <- runLFSPRO(temp.fam.data, cancer.data, counselee.id, mut.info = use.mutation)
          rltTmp[,7:ncol(rltTmp)][idx.dead,] <- NA
          rltTmp$info <- info
          rltTmp <- merge(rltTmp, subset(fam.data, select = c(id, vital, age)))
          LFSPRO.rlt <<- rltTmp
          
          if (use.CS == TRUE) {
            ProbLFSPRO = LFSPRO.rlt$carrier.cs
            LFSPRO = factor(ifelse(LFSPRO.rlt$carrier.cs>cutoff, "Yes", "No"), 
                            levels = c("Yes", "No"))
          } else {
            ProbLFSPRO = LFSPRO.rlt$carrier.mpc
            LFSPRO = factor(ifelse(LFSPRO.rlt$carrier.mpc>cutoff, "Yes", "No"), 
                            levels = c("Yes", "No"))
          }
          
          rlt <- data.frame(
            id = factor(LFSPRO.rlt$id, levels =  LFSPRO.rlt$id),
            info = LFSPRO.rlt$info,
            test = fam.genelevel(LFSPRO.rlt$test),
            ProbLFSPRO = ProbLFSPRO,
            LFSPRO = LFSPRO,
            Chompret = factor(ifelse(LFSPRO.rlt$chompret, "Yes", "No"),
                              levels = c("Yes", "No")),
            Classic = factor(ifelse(LFSPRO.rlt$classic, "Yes", "No"),
                             levels = c("Yes", "No")),
            breast.5 = LFSPRO.rlt$breast.5,
            breast.10 = LFSPRO.rlt$breast.10,
            breast.15 = LFSPRO.rlt$breast.15,
            sarcoma.5 = LFSPRO.rlt$sarcoma.5,
            sarcoma.10 = LFSPRO.rlt$sarcoma.10,
            sarcoma.15 = LFSPRO.rlt$sarcoma.15,
            other.5 = LFSPRO.rlt$other.5,
            other.10 = LFSPRO.rlt$other.10,
            other.15 = LFSPRO.rlt$other.15,
            second.5 = LFSPRO.rlt$second.5,
            second.10 = LFSPRO.rlt$second.10,
            second.15 = LFSPRO.rlt$second.15,
            figure = buttonInput(
              FUN = actionButton,
              len = nrow(LFSPRO.rlt),
              id = "button_",
              label = "Risk Trend",
              onclick = 'Shiny.onInputChange(\"lastClick\",  this.id)'
            )
          )
          
          rlt <- cbind('Details' = '&oplus;', rlt)
          rlt
        })
      }, container = sketch, filter = 'top', escape = FALSE,
      options = list(
        columnDefs = list(
          list(visible = FALSE, targets = c(0, 3)),
          list(orderable = FALSE, className = 'details-control', targets = 1)
        )
      ),
      callback = JS(
        "
        table.column(1).nodes().to$().css({cursor: 'pointer'});
        
        var format = function(d) {
          return '<div style=\"background-color:#eee; padding: .5em;\"> ' +
            d[3] + '</div>';
        };
        
        table.on('click', 'td.details-control', function() {
            var td = $(this), row = table.row(td.closest('tr'));
            if (row.child.isShown()) {
            row.child.hide();
            td.html('&oplus;');
            } else {
            row.child(format(row.data())).show();
            td.html('&CircleMinus;');
            }
        });
        "
      )) %>% 
        DT::formatRound('ProbLFSPRO', 2) %>%
        DT::formatRound(
          c("breast.5", "breast.10", "breast.15", "sarcoma.5", "sarcoma.10", "sarcoma.15", 
            "other.5", "other.10", "other.15", "second.5", "second.10", "second.15"), 2) %>%
        formatStyle(c('LFSPRO', 'Chompret', 'Classic'), color = styleEqual("Yes", 'red'))
    })
  })
  
  observeEvent(eventExpr = input$lastClick, handlerExpr = {
    myValue$idx.button <<- as.numeric(strsplit(input$lastClick, "_")[[1]][2])
    
    showModal(modalDialog(
      title =  "Cancer Risk",
      plotlyOutput("cancerrisk"),
      fluidRow(
        column(
          6,
          verbatimTextOutput("sampleInfo")
        ),
        column(
          6,
          verbatimTextOutput("cancerInfo")
        )
      )
    ))
  })
  
  output$sampleInfo <- renderText({
    idx.button <- myValue$idx.button
    id.sel <- LFSPRO.rlt$id[idx.button]
    fam.data <- famdata()
    cancer.data <- cancerdata()
    fam.data <- fam.update(fam.data, cancer.data)
    idx.sel <- which(fam.data$id == id.sel)
    rlt <- paste0(
      "Sample id: ", id.sel,  "\n",
      ifelse(
        fam.data$vital[idx.sel] == "A",
        paste0("Age: ", fam.data$age[idx.sel], "\n"),
        paste0("Died at age of: ", fam.data$age[idx.sel], "\n")
      ),
      "Sex: ", ifelse(fam.data$gender[idx.sel] == 0, "Female", "Male"), "\n")
    rlt
  })
  
  output$cancerInfo <- renderText({
    idx.button <- myValue$idx.button
    id.sel <- LFSPRO.rlt$id[idx.button]
    cancer.data <- cancerdata()
    idx.sel <- which(cancer.data$id == id.sel)
    rlt <- paste0("Sample id: ", id.sel, "\n")
    if(length(idx.sel)==0){
      rlt <- paste0(rlt, "No Cancer")
    } else {
      for(i in 1:length(idx.sel)){
        rlt <- paste0(rlt, cancer.data$cancer.type[idx.sel[i]], "\t", cancer.data$diag.age[idx.sel[i]], "\n")
      }
    }
    rlt
  })
  
  output$cancerrisk <- renderPlotly({
    idx.button <- myValue$idx.button
    
    dplot <- data.frame(
      year = c(5, 10, 15),
      Breast = c(LFSPRO.rlt$breast.5[idx.button], LFSPRO.rlt$breast.10[idx.button], LFSPRO.rlt$breast.15[idx.button]),
      Sarcoma = c(LFSPRO.rlt$sarcoma.5[idx.button], LFSPRO.rlt$sarcoma.10[idx.button], LFSPRO.rlt$sarcoma.15[idx.button]),
      Other = c(LFSPRO.rlt$other.5[idx.button], LFSPRO.rlt$other.10[idx.button], LFSPRO.rlt$other.15[idx.button]),
      Second = c(LFSPRO.rlt$second.5[idx.button], LFSPRO.rlt$second.10[idx.button], LFSPRO.rlt$second.15[idx.button]),
      stringsAsFactors = FALSE
    )
    
    Second.check <- sum(is.na(c(dplot$Second, dplot$Breast, dplot$Sarcoma, dplot$Other)))==12
    idx.rm.col <- colSums(!is.na(dplot)) > 0
    dplot <- dplot[,idx.rm.col]
    if(is.null(ncol(dplot))){
      text <- ''
      if(LFSPRO.rlt$vital[idx.button] == 'D'){
        text = paste("This individual is deceased and is not applicable \n for future cancer risk prediction.")
      } else if(LFSPRO.rlt$age[idx.button] >= 80){
        text = paste("At this time, LFSPRO has not been validated to \n",
                     "provide a risk prediction over the age of 80 years.")
      } else if(Second.check){
        text = paste("At this time, LFSPRO has not been validated to \n",
                     "provide a risk prediction past a second primary cancer.")
      }
      fig <- plot_ly(x = 2, y = 25, mode = 'text', text = text, textfont = list(size = 14)) %>%
        layout(xaxis = list(title = "", zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE),
               yaxis = list(title = "", zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE),
               plot_bgcolor = "white")
      return(fig)
    }
    
    dplot.pop <- data.frame(
      year = c(5, 10, 15),
      Breast = c(LFSPRO.rlt$pop.breast.5[idx.button], LFSPRO.rlt$pop.breast.10[idx.button], 
                 LFSPRO.rlt$pop.breast.15[idx.button]),
      Sarcoma = c(LFSPRO.rlt$pop.sarcoma.5[idx.button], LFSPRO.rlt$pop.sarcoma.10[idx.button], 
                  LFSPRO.rlt$pop.sarcoma.15[idx.button]),
      Other = c(LFSPRO.rlt$pop.other.5[idx.button], LFSPRO.rlt$pop.other.10[idx.button], 
                LFSPRO.rlt$pop.other.15[idx.button]),
      Second = c(LFSPRO.rlt$pop.second.5[idx.button], LFSPRO.rlt$pop.second.10[idx.button], 
                 LFSPRO.rlt$pop.second.15[idx.button]),
      stringsAsFactors = FALSE
    )
    
    Second.check <- sum(is.na(c(dplot$Second, dplot$Breast, dplot$Sarcoma, dplot$Other)))==12
    idx.rm.col <- colSums(!is.na(dplot.pop)) > 0
    dplot.pop <- dplot.pop[,idx.rm.col]
    if(is.null(ncol(dplot.pop))| LFSPRO.rlt$vital[idx.button] == 'D'){
      text <- ''
      if(LFSPRO.rlt$vital[idx.button] == 'D'){
        text = paste("This individual is deceased and is not applicable \n for future cancer risk prediction.")
      } else if(LFSPRO.rlt$age[idx.button] >= 80){
        text = paste("At this time, LFSPRO has not been validated to \n",
                     "provide a risk prediction over the age of 80 years.\n")
      } else if(Second.check){
        text = paste("At this time, LFSPRO has not been validated to \n",
                     "provide a risk prediction past a second primary cancer.\n")
      }
      fig <- plot_ly(x = 2, y = 25, mode = 'text', text = text, textfont = list(size = 14)) %>%
        layout(xaxis = list(title = "", zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE),
               yaxis = list(title = "", zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE),
               plot_bgcolor = "white")
      return(fig)
    }
    
    dplot <- reshape2::melt(dplot, id.vars = "year", variable.name = "type", value.name = "risk")
    dplot$risk[is.na(dplot$risk)] <- 0
    dplot$risk <- 100 * dplot$risk
    dplot$year <- factor(dplot$year)
    levels(dplot$year) <- c("5yr", "10yr", "15yr")
    if (length(levels(dplot$type)) == 1) {
      levels(dplot$type) <- "Second primary cancer"
    }
    
    dplot.pop <- reshape2::melt(dplot.pop, id.vars = "year", variable.name = "type", value.name = "risk")
    dplot.pop$risk[is.na(dplot.pop$risk)] <- 0
    
    dplot$risk_pop <- 100 * dplot.pop$risk
    dplot$show_legend <- (dplot$type %in% c("Breast", "Second primary cancer"))
    
    fig <- dplot %>%
      group_by(type) %>%
      do(p = plot_ly(., x = ~year, y = ~risk, type = "bar", alpha = 0.7,
                     color = I("royalblue4"), name = "Personalized cancer risk", 
                     showlegend = ~show_legend[1], 
                     hovertemplate = "%{y:.2f}%") %>%
           add_trace(., x = ~year, y = ~risk_pop, type = "bar", width = 0.5, alpha = 0.7,
                     color = I("firebrick1"), name = "Population cancer risk", 
                     showlegend = ~show_legend[1], 
                     hovertemplate = "%{y:.2f}%") %>%
           layout(barmode = 'overlay', 
                  xaxis = list(title = list(text = ~type[1], font = list(size = 12), standoff = 5)), 
                  yaxis = list(title = "Probability (%)"))) %>%
      subplot(nrows = 1, shareX = TRUE, shareY = TRUE) %>%
      layout(legend = list(orientation = "h", xanchor = "center", x = 0.5), 
             hovermode = "x unified")
    
    fig
    
  })
  
})

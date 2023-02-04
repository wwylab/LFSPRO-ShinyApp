library(shiny)
library(kinship2)
library(DT)
library(data.table)
library(shinyjs)
library(shinyalert)
library(reshape2)
library(ggplot2)
library(ggsci)
library(LFSPRO)

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
      th(rowspan = 3, 'ProbLFSPRO.mpc'),
      th(rowspan = 3, 'ProbLFSPRO.cs'),
      th(rowspan = 3, 'LFSPRO-mpc-carrier'),
      th(rowspan = 3, 'LFSPRO-cs-carrier'),
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
          col.label <- c("Unaffected")
          if (sum(fam.data$proband == "Y") > 0) {
            col.label <- c(col.label, "Proband")
          }
          if (sum(fam.data$dummy == 1) > 0) {
            col.label <- c(col.label, "Dummy")
          }
          legendPlot(pedSel,
                     affected.label = c("Affected"),
                     col = ifelse(fam.data$dummy == 0, ifelse(fam.data$proband == "Y", "red", "black"), "blue"),
                     col.label = col.label, id = id2)
          
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
          id.test <- fam.data$id[!is.na(fam.data$test)]
          id.cancer <- unique(cancer.data$id)
          non.informative <- FALSE
          if ((length(id.test) <= 1) & (length(id.cancer) <= 1)) {
            if ((length(id.test) == 1) & (length(id.cancer) == 1)) {
              non.informative <- (id.test[1] == id.cancer[1])
            } else {
              non.informative <- TRUE
            }
          }
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
            
            rltTmp <- runLFSPRO(temp.fam.data, cancer.data, counselee.id, mut.info = TRUE)
            rltTmp[,7:ncol(rltTmp)][idx.dead,] <- NA
            rltTmp$info <- info
            rltTmp <- merge(rltTmp, subset(fam.data, select = c(id, vital, age)))
            LFSPRO.rlt <<- rltTmp
            
            rlt <- data.frame(
              id =factor(LFSPRO.rlt$id, levels = LFSPRO.rlt$id),
              info = LFSPRO.rlt$info,
              test = fam.genelevel(LFSPRO.rlt$test),
              ProbLFSPRO.mpc = LFSPRO.rlt$carrier.mpc,
              ProbLFSPRO.cs = LFSPRO.rlt$carrier.cs,
              LFSPRO.mpc = factor(ifelse(LFSPRO.rlt$carrier.mpc>cutoff, "Yes", "No"), 
                                  levels = c("Yes", "No")),
              LFSPRO.cs = factor(ifelse(LFSPRO.rlt$carrier.cs>cutoff, "Yes", "No"), 
                                 levels = c("Yes", "No")),
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
          DT::formatRound('ProbLFSPRO.mpc', 2) %>%
          DT::formatRound('ProbLFSPRO.cs', 2) %>%
          DT::formatRound(
            c("breast.5", "breast.10", "breast.15", "sarcoma.5", "sarcoma.10", "sarcoma.15", 
              "other.5", "other.10", "other.15", "second.5", "second.10", "second.15"), 2) %>%
          formatStyle(c('LFSPRO.mpc', 'LFSPRO.cs', 'Chompret', 'Classic'),
                      color = styleEqual("Yes", 'red'))
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
      
      DT::datatable({
        isolate({
          rlt <- data.frame(
            id = factor(LFSPRO.rlt$id, levels =  LFSPRO.rlt$id),
            info = LFSPRO.rlt$info,
            test = fam.genelevel(LFSPRO.rlt$test),
            ProbLFSPRO.mpc = LFSPRO.rlt$carrier.mpc,
            ProbLFSPRO.cs = LFSPRO.rlt$carrier.cs,
            LFSPRO.mpc = factor(ifelse(LFSPRO.rlt$carrier.mpc>cutoff, "Yes", "No"), 
                                levels = c("Yes", "No")),
            LFSPRO.cs = factor(ifelse(LFSPRO.rlt$carrier.cs>cutoff, "Yes", "No"), 
                               levels = c("Yes", "No")),
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
        DT::formatRound('ProbLFSPRO.mpc', 2) %>%
        DT::formatRound('ProbLFSPRO.cs', 2) %>%
        DT::formatRound(
          c("breast.5", "breast.10", "breast.15", "sarcoma.5", "sarcoma.10", "sarcoma.15", 
            "other.5", "other.10", "other.15", "second.5", "second.10", "second.15"), 2) %>%
        formatStyle(c('LFSPRO.mpc', 'LFSPRO.cs', 'Chompret', 'Classic'),
                    color = styleEqual("Yes", 'red'))
    })
  })
  
  observeEvent(eventExpr = input$mutation, handlerExpr = {
    mutation <<- input$mutation
    use.mutation <- (mutation == "Yes")
    
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
          
          rlt <- data.frame(
            id = factor(LFSPRO.rlt$id, levels =  LFSPRO.rlt$id),
            info = LFSPRO.rlt$info,
            test = fam.genelevel(LFSPRO.rlt$test),
            ProbLFSPRO.mpc = LFSPRO.rlt$carrier.mpc,
            ProbLFSPRO.cs = LFSPRO.rlt$carrier.cs,
            LFSPRO.mpc = factor(ifelse(LFSPRO.rlt$carrier.mpc>cutoff, "Yes", "No"), 
                                levels = c("Yes", "No")),
            LFSPRO.cs = factor(ifelse(LFSPRO.rlt$carrier.cs>cutoff, "Yes", "No"), 
                               levels = c("Yes", "No")),
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
        DT::formatRound('ProbLFSPRO.mpc', 2) %>%
        DT::formatRound('ProbLFSPRO.cs', 2) %>%
        DT::formatRound(
          c("breast.5", "breast.10", "breast.15", "sarcoma.5", "sarcoma.10", "sarcoma.15", 
            "other.5", "other.10", "other.15", "second.5", "second.10", "second.15"), 2) %>%
        formatStyle(c('LFSPRO.mpc', 'LFSPRO.cs', 'Chompret', 'Classic'),
                    color = styleEqual("Yes", 'red'))
    })
  })
  
  observeEvent(eventExpr = input$lastClick, handlerExpr = {
    myValue$idx.button <<- as.numeric(strsplit(input$lastClick, "_")[[1]][2])
    
    showModal(modalDialog(
      title =  "Cancer Risk",
      plotOutput("cancerrisk"),
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
  
  output$cancerrisk <- renderPlot({
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
      
      gp <- ggplot() + 
        annotate("text", x = 4, y = 25, size=6, label = text) + theme_void()
      return(gp)
    }
    
    
    dplot.pop <- data.frame(
      year = c(5, 10, 15),
      
      Pop.Breast = c(LFSPRO.rlt$pop.breast.5[idx.button], LFSPRO.rlt$pop.breast.10[idx.button], 
                     LFSPRO.rlt$pop.breast.15[idx.button]),
      Pop.Sarcoma = c(LFSPRO.rlt$pop.sarcoma.5[idx.button], LFSPRO.rlt$pop.sarcoma.10[idx.button], 
                      LFSPRO.rlt$pop.sarcoma.15[idx.button]),
      Pop.Other = c(LFSPRO.rlt$pop.other.5[idx.button], LFSPRO.rlt$pop.other.10[idx.button], 
                    LFSPRO.rlt$pop.other.15[idx.button]),
      Pop.Second = c(LFSPRO.rlt$pop.second.5[idx.button], LFSPRO.rlt$pop.second.10[idx.button], 
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
      gp <- ggplot() + 
        annotate("text", x = 4, y = 25, size=6, label = text) + theme_void()
      return(gp)
    }
    
    dplot.2 <- reshape2::melt(dplot, id.vars = "year", variable.name = "type", value.name = "risk")
    dplot.2$risk[is.na(dplot.2$risk)] <- 0
    dplot.2$risk_per <- 100 * dplot.2$risk
    dplot.2$risk_per2 <- round(dplot.2$risk_per, 3)
    dplot.2$year <- factor(dplot.2$year)
    levels(dplot.2$year) <- c("5yr", "10yr", "15yr")
    dplot.2$pop <- factor(0)
    
    dplot.pop.2 <- reshape2::melt(dplot.pop, id.vars = "year", variable.name = "type", value.name = "risk")
    dplot.pop.2$risk[is.na(dplot.pop.2$risk)] <- 0
    dplot.pop.2$risk_per <- 100 * dplot.pop.2$risk
    dplot.pop.2$year <- factor(dplot.pop.2$year)
    levels(dplot.pop.2$year) <- c("5yr", "10yr", "15yr")
    levels(dplot.pop.2$type) <- levels(dplot.2$type)
    dplot.pop.2$pop <- factor(1)
    
    colors <- c("royalblue4", "firebrick1") 
    labels <- c("Personalized cancer risk", "Population cancer risk")
    
    gp <- ggplot() +
      geom_col(data = dplot.2, aes(x = year, y = risk_per, fill = pop), 
               alpha = 0.7, width = 0.6, position = position_dodge(0.7)) +
      geom_col(data = dplot.pop.2, aes(x = year, y = risk_per, fill = pop), 
               width = 0.4, position = position_dodge(0.7)) +
      geom_text(data = dplot.2, vjust = -0.25, 
                aes(x = year, y = risk_per, label = ifelse(risk_per2 == 0, 'NA', risk_per2))) +
      facet_wrap(~ type, strip.position = "bottom") +
      theme_light() + 
      scale_fill_manual(values = colors, labels = labels) +
      labs(x = NULL, y = "Probability (%)") + 
      theme(legend.title = element_blank(), legend.position = "bottom") + 
      theme(text = element_text(size = 12), 
            panel.spacing = unit(0, "lines"), 
            strip.background = element_blank(),
            strip.text = element_text(color = "black"),
            strip.placement = "outside")
    
    gp
    
  })
  
})

format.pv <- function(pv){
  if(pv < 0.001){
    return("<0.001")
  } else {
    return(format(round(pv, digits = 3), nsmall = 3))
  }
}

runLFSPRO <- function(fam.data, cancer.data, counselee.id){
  rlt.classic <- lfsClassic(fam.data, cancer.data, counselee.id)
  rlt.chompret <- lfsChompret2015(fam.data, cancer.data, counselee.id)
  rlt.lfspro <- lfspro(fam.data, cancer.data, counselee.id,method="MPC")
  rlt.lfspro.cs <- lfspro(fam.data, cancer.data, counselee.id,method="CS")
  rlt.lfspro.pop <- lfspro.pop(fam.data, cancer.data, counselee.id)
  
  cs.id <- rlt.lfspro$Cancer_specific_risks$Breast_risks$counselee.id
  mpc.id <- rlt.lfspro$Multiple_primary_cancer_risks$id
  cs.idx <- which(counselee.id$id %in% cs.id)
  mpc.idx <- which(counselee.id$id %in% mpc.id)
  
  risk <- rep(NA, nrow(counselee.id))
  
  rlt <- data.frame(
    id = counselee.id$id,
    gene.testing = fam.data[fam.data$id %in% counselee.id$id,]$gene.testing,
    classic = rlt.classic$result,
    chompret = rlt.chompret$result,
    carrier.mpc = rlt.lfspro$Mutation_probability$mutation_probability,
    carrier.cs = rlt.lfspro.cs$Mutation_probability$mutation_probability,
    
    breast.5 = replace(risk, cs.idx, rlt.lfspro$Cancer_specific_risks$Breast_risks$`breast in 5 yrs`),
    breast.10 = replace(risk, cs.idx, rlt.lfspro$Cancer_specific_risks$Breast_risks$`breast in 10 yrs`),
    breast.15 = replace(risk, cs.idx, rlt.lfspro$Cancer_specific_risks$Breast_risks$`breast in 15 yrs`),
    
    sarcoma.5 = replace(risk, cs.idx, rlt.lfspro$Cancer_specific_risks$Sarcoma_risks$`sarcoma in 5 yrs`),
    sarcoma.10 = replace(risk, cs.idx, rlt.lfspro$Cancer_specific_risks$Sarcoma_risks$`sarcoma in 10 yrs`),
    sarcoma.15 = replace(risk, cs.idx, rlt.lfspro$Cancer_specific_risks$Sarcoma_risks$`sarcoma in 15 yrs`),
    
    other.5 = replace(risk, cs.idx, rlt.lfspro$Cancer_specific_risks$Other_cancers$`others in 5 yrs`),
    other.10 = replace(risk, cs.idx, rlt.lfspro$Cancer_specific_risks$Other_cancers$`others in 10 yrs`),
    other.15 = replace(risk, cs.idx, rlt.lfspro$Cancer_specific_risks$Other_cancers$`others in 15 yrs`),
    
    second.5 = replace(risk, mpc.idx, rlt.lfspro$Multiple_primary_cancer_risks$`5 years`),
    second.10 = replace(risk, mpc.idx, rlt.lfspro$Multiple_primary_cancer_risks$`10 years`),
    second.15 = replace(risk, mpc.idx, rlt.lfspro$Multiple_primary_cancer_risks$`15 years`),
    
    pop.breast.5 = replace(risk, cs.idx, rlt.lfspro.pop$Cancer_specific_risks$Breast_risks$`breast in 5 yrs`),
    pop.breast.10 = replace(risk, cs.idx, rlt.lfspro.pop$Cancer_specific_risks$Breast_risks$`breast in 10 yrs`),
    pop.breast.15 = replace(risk, cs.idx, rlt.lfspro.pop$Cancer_specific_risks$Breast_risks$`breast in 15 yrs`),
    
    pop.sarcoma.5 = replace(risk, cs.idx, rlt.lfspro.pop$Cancer_specific_risks$Sarcoma_risks$`sarcoma in 5 yrs`),
    pop.sarcoma.10 = replace(risk, cs.idx, rlt.lfspro.pop$Cancer_specific_risks$Sarcoma_risks$`sarcoma in 10 yrs`),
    pop.sarcoma.15 = replace(risk, cs.idx, rlt.lfspro.pop$Cancer_specific_risks$Sarcoma_risks$`sarcoma in 15 yrs`),
    
    pop.other.5 = replace(risk, cs.idx, rlt.lfspro.pop$Cancer_specific_risks$Other_cancers$`others in 5 yrs`),
    pop.other.10 = replace(risk, cs.idx, rlt.lfspro.pop$Cancer_specific_risks$Other_cancers$`others in 10 yrs`),
    pop.other.15 = replace(risk, cs.idx, rlt.lfspro.pop$Cancer_specific_risks$Other_cancers$`others in 15 yrs`),
    
    pop.second.5 = replace(risk, mpc.idx, rlt.lfspro.pop$Multiple_primary_cancer_risks$`5 years`),
    pop.second.10 = replace(risk, mpc.idx, rlt.lfspro.pop$Multiple_primary_cancer_risks$`10 years`),
    pop.second.15 = replace(risk, mpc.idx, rlt.lfspro.pop$Multiple_primary_cancer_risks$`15 years`),
    
    stringsAsFactors = FALSE
  )
  rlt
}

dummy.create <- function(fam.data) {
  fam <- fam.data
  
  id <- fam$id 
  fid <- fam$fid
  fid[fid == 0] <- NA
  mid <- fam$mid
  mid[mid == 0] <- NA
  gender <- fam$gender
  age <- fam$age
  gene.testing <- fam$gene.testing
  vital <- fam$vital
  proband <- fam$proband
  dummy <- rep(0, length(id))
  
  dad.id.not.found <- fid[!(fid %in% id) & !is.na(fid)]
  n <- length(dad.id.not.found)
  id <- c(id, dad.id.not.found)
  fid <- c(fid, rep(NA, n))
  mid <- c(mid, rep(NA, n))
  gender <- c(gender, rep(1, n))
  age <- c(age, rep(1, n))
  gene.testing <- c(gene.testing, rep(NA, n))
  vital <- c(vital, rep('A', n))
  proband <- c(proband, rep('N', n))
  dummy <- c(dummy, rep(1, n))
  
  mum.id.not.found <- mid[!(mid %in% id) & !is.na(mid)]
  n <- length(mum.id.not.found)
  id <- c(id, mum.id.not.found)
  fid <- c(fid, rep(NA, n))
  mid <- c(mid, rep(NA, n))
  gender <- c(gender, rep(0, n))
  age <- c(age, rep(1, n))
  gene.testing <- c(gene.testing, rep(NA, n))
  vital <- c(vital, rep('A', n))
  proband <- c(proband, rep('N', n))
  dummy <- c(dummy, rep(1, n))
  
  for (j in 1:length(id)) {
    add <- FALSE
    if (is.na(fid[j]) & !is.na(mid[j])) {
      mid.j <- mid[j]
      fid.j <- fid[mid==mid.j]
      fid.j <- unique(na.omit(fid.j))
      if (length(fid.j) == 1){
        fid[j] <- fid.j
      }else if(length(fid.j) > 1){
        fid[j] <- sample(fid.j,1)
      }else{
        fid[j] <- max(id) + 1
        gender <- c(gender, 1)
        add <- TRUE
      }
    } else if (!is.na(fid[j]) & is.na(mid[j])) {
      fid.j <- fid[j]
      mid.j <- mid[fid==fid.j]
      mid.j <- unique(na.omit(mid.j))
      if (length(mid.j) == 1){
        mid[j] <- mid.j
      }else if (length(mid.j) > 1){
        mid[j] <- sample(mid.j,1)
      }else{
        mid[j] <- max(id) + 1
        gender <- c(gender, 0)
        add <- TRUE
      }
    }
    if (add == TRUE) {
      id <- c(id, max(id) + 1)
      fid <- c(fid, NA)
      mid <- c(mid, NA)
      age <- c(age, 1)
      gene.testing <- c(gene.testing, NA)
      vital <- c(vital, 'A')
      proband <- c(proband, 'N')
      dummy <- c(dummy, 1)
    }
  }
  fam <- data.frame(id, fid, mid, gender, age, gene.testing, vital, proband, dummy)
  return(fam)
}

fam.data.process <- function(fam.data, cancer.data) {
  pedigree.notes.1 <- fam.data[, "PedigreeNotes1"]
  n.1 <- sum(!is.na(pedigree.notes.1))
  
  if (n.1 > 0) {
    age.extra.1 <- rep(NA, n.1)
    for (i in 1:n.1) {
      note <- pedigree.notes.1[i]
      if (note != "") {
        age <- as.numeric(unlist(regmatches(note, gregexpr("[[:digit:]]+", note))))
        if (length(age) > 0) {
          age <- age[(age <= 100) & (age >= 1)]
          age <- mean(age)
          month.1 <- unlist(regmatches(note, gregexpr(" m ", note)))
          month.2 <- unlist(regmatches(note, gregexpr("month|months", note)))
          if (length(c(month.1, month.2)) > 0) {
            age <- ceiling(age / 12)
          }
          age.extra.1[i] <- age
        }
      }
    }
  }
  
  #######################
  
  pedigree.notes.2 <- fam.data[, "PedigreeNotes2"]
  n.2 <- sum(!is.na(pedigree.notes.2))
  
  if (n.2 > 0) {
    age.extra.2 <- rep(NA, n.2)
    for (i in 1:n.2) {
      note <- pedigree.notes.2[i]
      if (note != "") {
        age <- as.numeric(unlist(regmatches(note, gregexpr("[[:digit:]]+", note))))
        if (length(age) > 0) {
          age <- age[(age <= 100) & (age >= 1)]
          age <- mean(age)
          month.1 <- unlist(regmatches(note, gregexpr(" m ", note)))
          month.2 <- unlist(regmatches(note, gregexpr("month|months", note)))
          if (length(c(month.1, month.2)) > 0) {
            age <- ceiling(age / 12)
          }
          age.extra.2[i] <- age
        }
      }
    }
  }
  
  #######################
  
  pedigree.notes.3 <- fam.data[, "PedigreeNotes3"]
  n.3 <- sum(!is.na(pedigree.notes.3))
  
  if (n.3 > 0) {
    age.extra.3 <- rep(NA, n.3)
    for (i in 1:n.3) {
      note <- pedigree.notes.3[i]
      if (note != "") {
        age <- as.numeric(unlist(regmatches(note, gregexpr("[[:digit:]]+", note))))
        if (length(age) > 0) {
          age <- age[(age <= 100) & (age >= 1)]
          age <- mean(age)
          month.1 <- unlist(regmatches(note, gregexpr(" m ", note)))
          month.2 <- unlist(regmatches(note, gregexpr("month|months", note)))
          if (length(c(month.1, month.2)) > 0) {
            age <- ceiling(age / 12)
          }
          age.extra.3[i] <- age
        }
      }
    }
  }
  
  #######################
  
  age.extra <- rep(NA, nrow(fam.data))
  
  if (n.1 > 0) {
    idx.temp <- (is.na(age.extra) & !is.na(age.extra.1))
    age.extra[idx.temp] <- age.extra.1[idx.temp]
  }
  
  if (n.2 > 0) {
    idx.temp <- (is.na(age.extra) & !is.na(age.extra.2))
    age.extra[idx.temp] <- age.extra.2[idx.temp]
  }
  
  if (n.3 > 0) {
    idx.temp <- (is.na(age.extra) & !is.na(age.extra.3))
    age.extra[idx.temp] <- age.extra.3[idx.temp]
  }
  
  #######################
  
  missing.age.idx <- is.na(fam.data[, "age"])
  fam.data[missing.age.idx, "age"] <- age.extra[missing.age.idx]
  fam.data[is.na(fam.data[, "age"]), "age"] <- 1
  fam.data[fam.data[, "age"] == 0, "age"] <- 1
  
  fam.data <- subset(fam.data, select = -c(PedigreeNotes1, PedigreeNotes2, PedigreeNotes3))
  
  #######################
  
  id.list <- fam.data[, "id"]
  
  for (i in 1:length(id.list)) {
    id <- id.list[i]
    age <- fam.data[i, "age"]
    age.diag <- cancer.data[cancer.data$id == id, "diag.age"]
    age.diag.last <- ifelse(length(age.diag) == 0, 0, max(age.diag))
    if (!is.na(age.diag.last)) {
      if (age.diag.last > age) {
        fam.data[i, "age"] <- age.diag.last
      }
    }
  }
  
  #######################
  
  return(fam.data)
}

cancer.data.process <- function(fam.data, cancer.data) {
  for (i in 1:nrow(cancer.data)) {
    if (is.na(cancer.data[i, "diag.age"])) {
      id <- cancer.data[i, "id"]
      idx <- which(fam.data[, "id"] == id)
      age <- fam.data[idx, "age"]
      cancer.data[i, "diag.age"] <- age
    }
  }
  cancer.data[cancer.data[, "diag.age"] == 0, "diag.age"] <- 1
  
  return(cancer.data)
}

lfspro.pop <- function(fam.data, cancer.data, counselee.id, penetrance.all=NULL,
                       allef=list(c(0.9994,0.0006)), nloci=1, mRate=0.00012){
  fam.data <- fam.data[order(fam.data$fam.id, fam.data$id),]
  cancer.data <- cancer.data[order(cancer.data$fam.id, cancer.data$id),]
  
  num.cancer <- nrow(cancer.data) 
  cancer.type.num <- rep(-1, num.cancer)
  colnames(counselee.id) <- c("fam.id", "id")
  
  for(i in 1:num.cancer){
    tmp <- lfspro.cancer.type[cancer.data$cancer.type[i]]
    if(is.na(tmp)){
      print(paste("Cannot find cancer ", cancer.data$cancer.type[i], 
                  " in the LFSPRO predefined cancer type", sep = ""))
      print("LFSPRO predefined cancer types are: ")
      print(cancer.type.all)
      print("Please check the input cancer information data.")
      
      num.counselee <- nrow(counselee.id)
      pp <- rep(-1, num.counselee)
      
      rlt <- data.frame(cbind(counselee.id, pp),check.names = FALSE, stringsAsFactors = F)
      colnames(rlt) <- c("fam.id", "id", "pp")
      return(rlt)
    }
    cancer.type.num[i] <- tmp
  }
  cancer.data$cancer.type <- cancer.type.num
  fam.cancer.data <- combinedata(fam.data, cancer.data)
  
  num.fam <- length(fam.cancer.data)
  risk.mpc.output <- NULL
  risk.cs.output <- NULL
  risk.mpc.final <- data.frame()
  invalid_counselee <- data.frame()
  pp.all <- NULL
  counselee.id_new <- data.frame()
  
  for(i in 1:num.fam){
    cid_all <- counselee.id$id[counselee.id$fam.id == fam.cancer.data[[i]]$fam.id[1]]
    cid <- cid_all[fam.cancer.data[[i]]$vital[which(fam.cancer.data[[i]]$id %in% cid_all)] == "A"]
    if (length(cid_all)>length(cid)){
      print("Some input counselee are dead. Details in Table invalid_counselee.")
      cid_invalid <- cid_all[fam.cancer.data[[i]]$vital[which(fam.cancer.data[[i]]$id %in% cid_all)] == "D"]
      invalid_counselee_tmp <- data.frame(ID=cid_invalid, 
                                          fam=rep(fam.cancer.data[[i]]$fam.id[1],
                                                  length(cid_invalid)))
      invalid_counselee <- rbind(invalid_counselee, invalid_counselee_tmp)
    }
    
    counselee.id_new_temp <- data.frame(ID=cid, 
                                        fam=rep(fam.cancer.data[[i]]$fam.id[1],
                                                length(cid)))
    counselee.id_new <- rbind(counselee.id_new, counselee.id_new_temp)
    
    if(length(cid)<1){
      print(paste("Cannot find any counselee id in family ", 
                  fam.cancer.data[[i]]$fam.id[1], sep=""))
      next
    }
    
    ## Carrier probability calculation with MPC
    if (is.null(penetrance.all)){
      penetrance.all = parameter.mpc
    }
    data.obj <- convert.data(fam.cancer.data)
    data.obj1 <- data.obj[[1]]
    data.obj2 <- data.obj[[2]]
    pp.tmp <- lfsproC.mpc(fam.cancer.data[[i]], penetrance.all, data.obj1[[i]],
                          data.obj2[[i]], cid, allef, nloci, mRate) 
    pp.tmp[,1] <- 1
    pp.tmp[,2] <- 0
    pp.tmp[,3] <- 0
    pp.all <- rbind(pp.all, pp.tmp)
    
    ## risk prediction
    cid_num.cancer <- fam.cancer.data[[i]]$num.cancer[which(fam.cancer.data[[i]]$id %in% cid)] 
    cid.na <- cid[which(cid_num.cancer==0)] #counselee without previous cancers
    cid.1 <- cid[which(cid_num.cancer>=1)] #counselee with previous primary cancer
    
    pp.na <- pp.tmp[which(cid_num.cancer==0),]
    dim(pp.na) <- c(sum(cid_num.cancer==0), 3)
    pp.1 <- pp.tmp[which(cid_num.cancer>=1),]
    dim(pp.1) <- c(sum(cid_num.cancer>=1), 3)
    if (length(cid.na)>0){
      pp.na[,1] <- 1
      pp.na[,2] <- 0
      pp.na[,3] <- 0
      risk.cs.temp <- risk.cs(fam.cancer.data[[i]], lfspenet.cs, cid.na, pp.na)
    } else {
      risk.cs.temp <- NULL
    }
    if (is.null(risk.cs.output)) {
      risk.cs.output <- c(risk.cs.output, risk.cs.temp)
    } else {
      risk.cs.output <- Map(list,risk.cs.output,risk.cs.temp)
    }
    
    if (length(cid.1)>0){
      risk.mpc.temp <- risk.mpc(fam.cancer.data[[i]], cid.1, penetrance.all)
      risk.mpc.output <- data.frame(risk.mpc.temp, stringsAsFactors = F)
      colnames(risk.mpc.output) <- c("fam.id", "id", "age", "gender", "first.cancer",
                                     "5 years (wildtype)", "10 years(wildtype)", "15 years (wildtype)", 
                                     "5 years (mutation)", "10 years(mutation)", "15 years (mutation)")
      counselee.id[which(cid_num.cancer>=1),]
      counselee.id.1 <- data.frame(fam.id=fam.cancer.data[[i]]$fam.id[1], id=cid.1)
      risk.all <- combined.risk.mpc(pp.1, risk.mpc.output, counselee.id.1)
      risk.mpc.final <- rbind(risk.mpc.final, risk.all)
    }
  }
  #browser()
  pp <- 1 - pp.all[, 1]
  rlt <- data.frame(cbind(counselee.id_new, pp), check.names = FALSE, stringsAsFactors = F)
  colnames(rlt) <- c("fam.id", "id", "mutation_probability")
  output <- list(rlt, risk.cs.output, na.omit(risk.mpc.final), invalid_counselee)
  names(output) <- c("Mutation_probability", "Cancer_specific_risks",
                     "Multiple_primary_cancer_risks", "Invalid_counselee")
  return(output)
}

##############################Updated ages/correct gender in family data

fam.update <- function(fam.data, cancer.data){
  id <- fam.data$id
  cancer.id <- unique(cancer.data$id)
  fam.data$gender[fam.data$gender==2] <- 0
  if (length(cancer.id) > 0){
    for (i in cancer.id){
      row.i <- which(id == i)
      age.i <- fam.data$age[row.i]
      diag.age.i <- cancer.data$diag.age[cancer.data$id == i]
      if (sum(is.na(diag.age.i)) < length(diag.age.i)){
        max.diag <- max(na.omit(diag.age.i))
        if (is.na(age.i) || (age.i < max.diag)){fam.data$age[row.i] = max.diag}
      }
    }
  }
  if (!("gene.testing" %in% colnames(fam.data))) {
    fam.data$gene.testing <- ""
  }
  return(fam.data)
}

fam.genelevel <- function(gene.testing){
  gene.abbrev <- c('TN', 'UN', 'VUS', 'LPV', 'PV', 'Mosaic')
  gene.test.full <- c('True Negative','Uniformative negative','VUS',
                      'Likely pathogenic','Pathogetnic','Suspected Mosaic')
  replace(gene.testing, gene.testing %in% gene.test.full, 
          gene.abbrev[na.omit(match(gene.testing, gene.test.full))])
}


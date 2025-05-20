library(shiny)
library(shinybusy)
library(ggplot2)
library(ggpubr)  # For theme_pubclean()
library(ggtext)
library(eefAnalytics)
library(lme4)
library(shinyjs)
library(colourpicker)
library(plyr)

ComparePlot5 <- function(eefAnalyticsList, group, Conditional = TRUE, ES_Total = TRUE, modelNames,
                         custom_labels = NULL, custom_colors = NULL,
                         custom_title = NULL, x_label = NULL, y_label = NULL,
                         vline_color = "black"){
  if(!is(eefAnalyticsList,"list")){stop("eefAnalyticsList is not a list.")}
  
  if(!all(unlist(lapply(eefAnalyticsList,function(x) is(x,"eefAnalytics"))))){stop("Not all list objects are a eefAnalytics class object.")}
  
  if(missing(modelNames)){stop("modelNames must be specified.")}
  if(missing(group)){stop("group must be specified.")}
  
  plotObject( analyticObject = eefAnalyticsList,
              group = group,
              Conditional = Conditional,
              ES_Total = ES_Total,
              compare = TRUE,
              modelNames = modelNames,
              custom_labels = custom_labels,
              custom_colors = custom_colors,
              custom_title = custom_title,
              x_label = x_label,
              y_label = y_label,
              vline_color = vline_color)
  
}




##################################
#      Internal plot function    #
##################################


plotObject <- function(analyticObject, group, Conditional, ES_Total, slope, compare, modelNames,
                       custom_labels = NULL, custom_colors = NULL,
                       custom_title = NULL, x_label = NULL, y_label = NULL,
                       vline_color = "black", ...){
  
  if(Conditional ==TRUE){analyticObject2=analyticObject; Condname="Conditional"}
  if(Conditional ==FALSE){analyticObject2=analyticObject$Unconditional; Condname="Unconditional"}
  if(ES_Total ==TRUE){ES_TW<-"Total"}
  if(ES_Total ==FALSE){ES_TW<-"Within"}
  
  if(compare==TRUE & !is.null(names(analyticObject))){stop("Specify the list of objects to compare")}
  if(!is.null(group)){
    trtname <-rownames(analyticObject2$ES)
    if(is.null(trtname)){trtname <-names(analyticObject2$ES)}
    trtpos <-substr(trtname, (nchar(trtname)-nchar(group)+1),  nchar(trtname))==group
    trt <- trtname[trtpos]
  }
  #bootstrap, permutation and Pprobability plot for SRT model
  #---------------------------------------------------------
  if(sum(analyticObject$Method=="LM") ==1) {
    if(is.null(group)){stop("Group must be specified.")}
    if(!is.null(group)& sum(names(analyticObject)=="ProbES")== 0){
      if(sum(names(analyticObject)=="Bootstrap"|
             names(analyticObject)=="permES")==0){stop("Only relevant for bootstrapped or permutated values")}
      if(sum(names(analyticObject)=="Bootstrap")==1){
        ntp <- nrow(as.matrix(analyticObject$ES))
        if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
        
        obs.est <- analyticObject2$ES[trt,1]
        tmp2 <- as.numeric(analyticObject2$Bootstrap[,trt])
        xlabs=paste0(Condname," Bootstrap estimates")
        hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main="")
        abline(v=obs.est,col="red",lwd=2,lty=1)
        abline(v=0,col="grey48",lwd=2,lty=1)
        legend("topright",c("Observed Estimate","Zero-Line"),col=c("red","grey48"),bty="n",lty=1,lwd=2)
      }
      
      if(sum(names(analyticObject2)=="permES")==1){
        ntp <- nrow(as.matrix(analyticObject2$ES))
        if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
        
        
        Perm.names <-names(analyticObject2$permES)
        
        obs.est <- analyticObject2$ES[trt,1]
        tmp2 <- as.numeric(analyticObject2$permES[,trt])
        pvalue <- round(mean(abs(tmp2)> abs(obs.est)),3)
        xlabs=paste0("Permutation values (PermES) based on ",Condname, " ES")
        
        hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main=paste("P(|PermES| > |ES|)=",pvalue,sep=""))
        abline(v=obs.est,col="red",lwd=2,lty=2)
        abline(v=-obs.est,col="red",lwd=2,lty=2)
        legend("topright",c("(-) Observed Estimate"),col=c("red"),bty="n",lty=2,lwd=2)
      }
    }
    if(sum(names(analyticObject)=="ProbES")> 0 ){
      if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
      tmp2 <- analyticObject2$ProbES[[which(trtpos==TRUE)]]
      
      thd0<- regmatches(rownames(tmp2),  gregexpr("[[:digit:]]+\\.*[[:digit:]]*",rownames(tmp2)))
      thd <- as.numeric(unlist(thd0))
      
      par_original <- par()[c("mar","xpd")]
      par_original0<- par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(thd, tmp2[,1], col="blue", type="b",ylim=c(0,max(tmp2)), xlab="Threshold", ylab="Posterior probability")
      on.exit(par(par_original0))
      on.exit(par(par_original))
    }
    
  }
  
  
  #bootstrap, permutation, and Pprobability plot for CRT and MST model
  #--------------------------------------------------------------------
  if(sum(analyticObject$Method=="MLM")==1){
    
    if(is.null(group)){
      
      tmp000  <- data.frame(analyticObject$SchEffects)#use analyticObject since both (un)condition has the same SchEffects object.
      if(slope==FALSE){
        tmp00 <- tmp000[,grep("Schools|Intercept|Estimate",names(tmp000))]
        mar11 <-  c(5, 4, 4, 2) + 0.1
      }
      if(slope==TRUE & dim(tmp000)[2] ==2){stop("x must be mstFREQ or mstBAyes object")}
      if(slope==TRUE & dim(tmp000)[2] >2){
        tmp00 <- tmp000[,!(names(tmp000) %in% "Intercept")]
        if(dim(tmp00)[2]==2){mar11 <- c(5, 4, 4, 2) + 0.1}
        if(dim(tmp00)[2]==3){mar11 <- c(5, 2, 4, 0) + 1.0}
        if(dim(tmp00)[2] >3){mar11 <- c(3, 2, 0, 0) + 1.0}}
      
      op <- par(mfrow = c(floor(dim(tmp00)[2]/2),round(dim(tmp00)[2]/2)),
                mar = mar11)
      
      
      for(i in 2:dim(tmp00)[2]){
        tmp <- data.frame(y=tmp00[,i],x=c(1:length(tmp00[,i])))
        tmp2 <- tmp[order(tmp$y),]
        ylabs=gsub("trt", "Intervention ",gsub("Estimate","Intercept", names(tmp00)[i]))
        barplot(tmp2$y,names.arg=tmp2$x,las=2,col="cornflowerblue",border="cornflowerblue")
        if(dim(tmp00)[2]<=2){mtext(ylabs, side = 2.5, line = 2)}
        if(dim(tmp00)[2] >2){mtext(ylabs, side = 2, line = 1.7, cex = 0.8)}
      }
      lines1=-2.5
      if(dim(tmp00)[2] >3){lines1=-1}
      title(xlab="School labels", outer = TRUE, line = lines1,cex.lab = 1.2)
      on.exit(par(op))
    }
    
    
    
    
    
    if( !is.null(group) & sum(names(analyticObject)=="Bootstrap")>0){
      ntp <- length(analyticObject2$ES)
      if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
      Boot.names <-names(analyticObject2$Bootstrap)
      obs.est <- analyticObject2$ES[[trt]][ES_TW,1]
      tmp2 <- as.numeric(analyticObject2$Bootstrap[,grep(ES_TW,grep(trt,Boot.names, ignore.case = T, value = T))])
      pvalue <- round(mean(abs(tmp2)> abs(obs.est)),3)
      xlabs=paste0("Bootstrap estimates for ",Condname, " ES_",ES_TW)
      hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main="")
      abline(v=obs.est,col="red",lwd=2,lty=1)
      abline(v=0,col="grey48",lwd=2,lty=1)
      legend("topright",c("Observed Estimate","Zero-Line"),col=c("red","grey48"),bty="n",lty=1,lwd=2)
    }
    
    
    if( !is.null(group) & sum(names(analyticObject)=="permES")>0){
      ntp <- ifelse(is.list(analyticObject$ES),length(analyticObject$ES),1)
      if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
      Perm.names <-names(analyticObject2$permES)
      obs.est <- analyticObject2$ES[[trt]][ES_TW,1]
      tmp2 <- as.numeric(analyticObject2$permES[,grep(ES_TW,grep(trt,Perm.names, ignore.case = T, value = T))])
      pvalue <- round(mean(abs(tmp2)> abs(obs.est)),3)
      xlabs=paste0("Permutation values(PermES) based on ",Condname, " ES_",ES_TW)
      hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main=paste("P(|PermES| > |ES|)=",pvalue,sep=""))
      abline(v=obs.est,col="red",lwd=2,lty=2)
      abline(v=-obs.est,col="red",lwd=2,lty=2)
      legend("topright",c("Observed Estimate"),col=c("red"),bty="n",lty=2,lwd=2)
    }
    
    
    if( !is.null(group) &sum(names(analyticObject)=="ProbES")> 0 ){
      
      if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
      tmp2 <- analyticObject2$ProbES[[which(trtpos==TRUE)]]
      tmp2.within<- tmp2[, grep("with",names(tmp2), ignore.case = TRUE)]
      tmp2.total <- tmp2[, grep("total",names(tmp2), ignore.case = TRUE)]
      thd0<- regmatches(rownames(tmp2),  gregexpr("[[:digit:]]+\\.*[[:digit:]]*",rownames(tmp2)))
      thd <- as.numeric(unlist(thd0))
      
      par_original <- par()[c("mar","xpd")]
      op<- par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(thd, tmp2.within, col="blue", type="b",ylim=c(0,max(tmp2)), xlab="Threshold", ylab="Posterior probability",...)
      lines(thd, tmp2.total, col="red", type="b", lty=2)
      legend("topright", legend=c("within", "total"), col=c("blue", "red"), lty=1:2, cex=0.8)
      on.exit(par(op))
      on.exit(par(par_original))
      
    }
    
  }
  
  # error bar for model comparing models
  #-------------------------------------
  if(is.null(names(analyticObject))){
    
    ltp <- names(analyticObject)
    if(!is.null(ltp)){stop("Specify list of eefAnalytics objects for comparison")}
    if(is.null(group)){stop("Group number must be defined")}
    ntp <- length(analyticObject)
    if(length(modelNames)!= ntp){stop("Names must be equal to the number of eefAnalytics objects")}
    
    es.mean <- es.lower <- es.upper <- p.name <- var.name <- NULL
    for(k in 1:ntp){
      tmp <- analyticObject[[k]]
      
      if(tmp$Method=="LM"){
        trtname <-rownames(tmp$ES)
        trtpos <-substr(trtname, (nchar(trtname)-nchar(group)+1),  nchar(trtname))==group
        trt <- trtname[trtpos]
        
        if(Conditional==TRUE){tmp2 <- as.matrix(tmp$ES)}
        if(Conditional==FALSE){tmp2 <- as.matrix(tmp$Unconditional$ES)}
        trtname <-rownames(tmp$ES)
        if(is.null(trtname)){trtname <-names(tmp$ES)}
        trtpos <-substr(trtname, (nchar(trtname)-nchar(group)+1),  nchar(trtname))==group
        trt <- trtname[trtpos]
        if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
        es.mean1 <-tmp2[trt,1]
        es.lower1 <-tmp2[trt,2]
        es.upper1 <-tmp2[trt,3]
        p.name1 <- rep(modelNames[k],length(es.mean1))
        var.name1 <- rep("Within",length(es.mean1))
        
      }
      
      if(tmp$Method=="MLM"){
        trtname <-names(tmp$ES)
        if(is.null(trtname)){trtname <-names(tmp$ES)}
        trtpos <-substr(trtname, (nchar(trtname)-nchar(group)+1),  nchar(trtname))==group
        trt <- trtname[trtpos]
        if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
        if(Conditional==TRUE) {tmp2 <- tmp$ES[[trt]]}
        if(Conditional==FALSE){tmp2 <- tmp$Unconditional$ES[[trt]]}
        es.mean1 <- tmp2[,1]
        es.lower1 <-tmp2[,2]
        es.upper1 <-tmp2[,3]
        p.name1 <- rep(modelNames[k],length(es.mean1))
        var.name1 <- rownames(tmp2)
        
        
      }
      
      
      es.mean <- c(es.mean,es.mean1)
      es.lower <- c(es.lower,es.lower1)
      es.upper <- c(es.upper,es.upper1)
      p.name <- c(p.name,p.name1)
      var.name <- c(var.name,var.name1)
      
    }
    
    
    MyData1 <- data.frame(
      ES = es.mean,
      LB.95 = es.lower,
      UB.95 = es.upper,
      Variance = var.name,
      Name = p.name
    )
    
    # Replace model names with custom labels (if provided)
    if (!is.null(custom_labels)) {
      old_names <- unique(MyData1$Name)
      if (length(custom_labels) == length(old_names)) {
        for (i in seq_along(old_names)) {
          MyData1$Name[MyData1$Name == old_names[i]] <- custom_labels[i]
        }
      }
    }
    
    MyData1 <- MyData1[order(MyData1$Variance, decreasing = TRUE), ]
    MyData1$Anot <- paste0(round(MyData1$ES, 2), " [", round(MyData1$LB.95, 2), ", ", round(MyData1$UB.95, 2), "]")
    MyData1$Xaxis <- max(MyData1$UB.95, na.rm = TRUE) + 0.05
    
    # Custom colors
    if (!is.null(custom_colors)) {
      color_map <- setNames(custom_colors, unique(MyData1$Name))
    } else {
      color_map <- NULL
    }
    
    # Breaks and limits
    Mybreaks <- round(c(min(MyData1$LB.95), mean(c(min(MyData1$LB.95), max(MyData1$UB.95))), max(MyData1$UB.95)), 2)
    xlimits <- c(min(MyData1$LB.95, 0), max(MyData1$UB.95) + 0.6)
    
    # Start plot
    p <- ggplot(MyData1, aes(x = ES, y = Name)) +
      geom_point(aes(color = Name)) +
      geom_errorbarh(aes(xmin = LB.95, xmax = UB.95, color = Name), height = 0.1) +
      geom_text(aes(x = Xaxis, label = Anot), hjust = 0, nudge_x = 0.04, size = 4) +
      scale_x_continuous(
        limits = xlimits,
        breaks = Mybreaks,
        name = ifelse(is.null(x_label), expression(paste("Hedges' ", italic("g"))), x_label)
      ) +
      ylab(ifelse(is.null(y_label), "Models", y_label)) +
      geom_vline(xintercept = 0, color = ifelse(is.null(vline_color), "black", vline_color), linetype = "dashed", alpha = 0.5) +
      theme_bw() +
      theme(
        axis.text.y = element_text(color = "black"),
        text = element_text(size = 16),
        panel.spacing = unit(1, "lines"),
        plot.margin = unit(c(30, 5, 5, 5), "point"),
        legend.position = "none"
      ) +
      coord_cartesian(clip = "off")
    
    # Facet if Total and Within are both present
    if (sum(unique(MyData1$Variance) %in% "Total") > 0) {
      p <- p + facet_grid(Variance ~ ., scales = "free", space = "free")
    }
    
    # Set title if provided
    if (!is.null(custom_title)) {
      p <- p + ggtitle(custom_title)
    }
    
    # Apply manual color scale if custom colors are provided
    if (!is.null(color_map)) {
      p <- p + scale_color_manual(values = color_map)
    }
  }
  return(p)
  
}


# Define the CRT and MST data simulation functions
# --- CRT Function with Covariates ---
# Updated CRT Data Simulation Function
crtdata_simulation <- function(nt, n_schools_treated, np, ns, sigma, ICC, B0, es, seed, attrition_rates, covariates) {
  set.seed(seed)
  
  # Error checking
  if (length(n_schools_treated) != nt + 1) {
    stop("Error: The length of n_schools_treated must be number of interventions + 1 (including control group).")
  }
  if (length(es) != nt) {
    stop("Error: The length of es must be equal to number of interventions.")
  }
  if (length(attrition_rates) != nt + 1) {
    stop("Error: The length of attrition_rates must be number of interventions + 1 (including control group).")
  }
  set.seed(seed)
  
  # Step 1: Create the base data structure with people and schools
  data <- expand.grid(pupils = 1:np, schools = 1:ns)
  data$pupils <- 1:nrow(data)
  
  # Step 2: Initialize treatment assignment and assign treatment groups to schools
  interventions <- "interventions"
  unique_schools <- unique(data$schools)
  
  # Initialize treatment column (0 = control)
  data[[interventions]] <- 0
  
  # Check if the total number of schools assigned matches the number of schools (ns)
  if (sum(n_schools_treated) != ns) {
    stop("Error: The total number of schools in n_schools_treated must be equal to the total number of schools (ns).")
  }
  
  # Assign schools to control group and treatment groups
  remaining_schools <- unique_schools
  
  # Step 2a: Assign schools to the control group
  control_schools <- sample(remaining_schools, n_schools_treated[1])
  remaining_schools <- setdiff(remaining_schools, control_schools)  # Remove control schools from remaining
  
  # Step 2b: Assign schools to each treatment group
  for (i in 1:nt) {
    treated_schools <- sample(remaining_schools, n_schools_treated[i + 1])
    data[[interventions]][data$schools %in% treated_schools] <- i  # Assign treatment group i
    remaining_schools <- setdiff(remaining_schools, treated_schools)  # Remove assigned schools
  }
  
  # Generate covariates (all individuals, no NA)
  for (cov in covariates) {
    if (cov$type == "continuous") {
      data[[cov$name]] <- rnorm(nrow(data), mean = 0, sd = cov$sd)
    } else if (cov$type == "categorical") {
      data[[cov$name]] <- sample(cov$levels, nrow(data), replace = TRUE, prob = cov$probs)
    }
  }
  
  # Convert categorical covariates to numeric (0 = reference)
  for (cov in covariates) {
    if (cov$type == "categorical") {
      levels_ordered <- c(cov$reference, setdiff(cov$levels, cov$reference))
      data[[cov$name]] <- match(data[[cov$name]], levels_ordered) - 1
    }
  }
  
  # Initialize post-test scores
  posts <- "posttest"
  data[[posts]] <- NA  # Initially set all to NA
  
  # Step 4: Handle attrition for control and treatment groups
  non_attrition_idx <- 1:nrow(data)  # Start by assuming no one is attrited
  
  for (i in 0:nt) {
    group_idx <- which(data$interventions == i)
    drop <- round(length(group_idx) * attrition_rates[i + 1])
    if (drop > 0) {
      drop_ids <- sample(group_idx, drop)
      data$posttest[drop_ids] <- NA
      non_attrition_idx <- setdiff(non_attrition_idx, drop_ids)
    }
  }
  
  # Step 5: Generate random individual errors and cluster-level random effects
  e_ij <- rnorm(length(non_attrition_idx), mean = 0, sd = sigma)  # Individual-level errors for non-attrited individuals
  sigmab <- sqrt(ICC / (1 - ICC) * sigma^2)  # School-level variance for all schools
  b_i_full <- rnorm(ns, mean = 0, sd = sigmab)  # Cluster-level random effects for all schools
  b_i <- b_i_full[data$schools[non_attrition_idx]]  # Get the random effects for relevant schools
  
  # Step 6: Create treatment effects for each group
  treatment_effects <- sapply(1:nt, function(i) as.numeric(data[non_attrition_idx, interventions] == i))
  
  treatment_effects_sizes <- sweep(treatment_effects, 2, es* sqrt(sigmab^2 + sigma^2), `*`)
  
  # Calculate the total treatment effect by combining individual treatment effects with corresponding effect sizes
  total_treatment_effect <- rowSums(treatment_effects_sizes)
  
  
  cov_effects <- rep(0, nrow(data))
  if (length(covariates) > 0) {
    for (cov in covariates) {
      if (cov$type == "continuous") {
        cov_effects <- cov_effects + cov$effect * data[[cov$name]]
      } else if (cov$type == "categorical") {
        for (lvl in cov$levels) {
          if (lvl != cov$reference) {
            cov_effects <- cov_effects + ifelse(data[[cov$name]] == lvl, cov$effects[[lvl]], 0)
          }
        }
      }
    }
  }
  
  # Apply only to non-attrited rows
  # Step 7: Compute post-test scores using the combined formula
  data[non_attrition_idx, posts] <- B0 + cov_effects[non_attrition_idx] + total_treatment_effect + b_i + e_ij
  
  return(data)
}


# Updated MST Data Simulation Function
mstdata_simulation <- function(nt, tpi, np, ns, sigma, sigmab0, sigmab1, B0, es, seed, attrition_rates, covariates) {
  # Error checking: ensure tpi has length nt + 1 (control + treatment groups)
  if (length(tpi) != nt + 1) {
    stop("Error: 'tpi' must have length nt + 1 (first value for control group, remaining for treatment groups).")
  }
  
  # Error checking: ensure the sum of tpi is 100
  if (sum(tpi) != 100) {
    stop("Error: The sum of 'tpi' must be 100%.")
  }
  
  # Error checking: ensure es has length equal to nt
  if (length(es) != nt) {
    stop("Error: 'es' must have length nt (one for each treatment group).")
  }
  
  # Error checking: ensure attrition_rates has length equal to nt + 1
  if (length(attrition_rates) != nt + 1) {
    stop("Error: 'attrition_rates' must have length nt + 1 (one for control group and one for each treatment group).")
  }
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Step 1: Create the base data structure with people and schools
  data <- expand.grid(pupils = 1:np, schools = 1:ns)
  data$pupils <- 1:nrow(data)
  
  # Initialize treatment column (0 = control)
  interventions <- "interventions"
  data[[interventions]] <- 0
  
  # Total number of participants
  total_participants <- nrow(data)
  
  # Step 2: Assign control group based on tpi[1] (control percentage)
  control_size <- floor(total_participants * (tpi[1] / 100))
  control_pupils <- sample(1:total_participants, control_size)
  data[[interventions]][control_pupils] <- 0  # Assign control group (0)
  
  # Step 3: Assign treatment groups to the remaining pupils
  remaining_pupils <- setdiff(1:total_participants, control_pupils)
  remaining_n <- length(remaining_pupils)
  
  if (nt == 1) {
    # Only one treatment group: assign all remaining pupils to it
    data[[interventions]][remaining_pupils] <- 1
  } else {
    # More than one treatment group
    treatment_props <- tpi[2:(nt+1)] / sum(tpi[2:(nt+1)])  # Normalize proportions
    treatment_sizes <- round(remaining_n * treatment_props)
    
    # Adjust last group to absorb rounding error
    treatment_sizes[nt] <- remaining_n - sum(treatment_sizes[1:(nt-1)])
    
    for (i in 1:nt) {
      treated_pupils <- sample(remaining_pupils, treatment_sizes[i])
      data[[interventions]][treated_pupils] <- i
      remaining_pupils <- setdiff(remaining_pupils, treated_pupils)
    }
  }
  
  
  for (cov in covariates) {
    if (cov$type == "continuous") {
      data[[cov$name]] <- rnorm(nrow(data), mean = 0, sd = cov$sd)
    } else {
      data[[cov$name]] <- sample(cov$levels, nrow(data), replace = TRUE, prob = cov$probs)
    }
  }
  
  # Convert categorical covariates to numeric (0 = reference)
  for (cov in covariates) {
    if (cov$type == "categorical") {
      levels_ordered <- c(cov$reference, setdiff(cov$levels, cov$reference))
      data[[cov$name]] <- match(data[[cov$name]], levels_ordered) - 1
    }
  }
  
  
  # Initialize post-test scores
  posts <- "posttest"
  data[[posts]] <- NA  # Initially set all to NA
  
  # Step 5: Handle attrition for control and treatment groups
  non_attrition_idx <- 1:nrow(data)  # Start by assuming no one is attrited
  
  for (i in 0:nt) {
    group_idx <- which(data[[interventions]] == i)
    attrition_rate <- attrition_rates[i + 1]  # Use i+1 because first value is for control group
    attrition_size <- round(length(group_idx) * attrition_rate)
    
    if (attrition_size > 0) {
      attrition_idx <- sample(group_idx, attrition_size)
      data[attrition_idx, posts] <- NA  # Mark attrited participants
      non_attrition_idx <- setdiff(non_attrition_idx, attrition_idx)  # Update non-attrited participants
    }
  }
  
  # Step 6: Generate random individual errors and cluster-level random effects
  e_ij <- rnorm(length(non_attrition_idx), mean = 0, sd = sigma)  # Individual-level errors for non-attrited individuals
  
  # Cluster-level random effects
  b_i_full <- rnorm(ns, mean = 0, sd = sigmab0)  # Random intercepts for schools
  b1_i_full <- rnorm(ns, mean = 0, sd = sigmab1)  # Random slopes for treatment effect
  
  # Get the random effects for the relevant schools
  b_i <- b_i_full[data$schools[non_attrition_idx]]
  b1_i <- b1_i_full[data$schools[non_attrition_idx]]
  
  # Step 7: Use model.matrix to create treatment effect dummy variables matrix
  unique_interventions <- unique(data[non_attrition_idx, interventions])
  
  
  if (length(unique_interventions) > 1) {
    treatment_matrix <- model.matrix(~ factor(data[non_attrition_idx, interventions]) - 1)
    treatment_matrix <- treatment_matrix[, -1, drop = FALSE]  # Remove control column, keep matrix structure
    
    if (is.null(dim(treatment_matrix))) {
      treatment_matrix <- matrix(treatment_matrix, ncol = 1)
    }
    
    treatment_effects <- sweep(treatment_matrix, 2, es * sqrt(sigmab0^2 + sigmab1^2 + sigma^2), `*`)
    total_treatment_effect <- rowSums(treatment_effects)
    random_slope <- treatment_matrix * b1_i
    total_random_slope_effect <- rowSums(random_slope)
  } else {
    total_treatment_effect <- rep(0, length(non_attrition_idx))
    total_random_slope_effect <- rep(0, length(non_attrition_idx))
  }
  
  
  # Calculate the total treatment effect by summing the contributions of each treatment group
  total_treatment_effect <- rowSums(treatment_effects) 
  
  random_slope <-  treatment_matrix*b1_i
  
  # Calculate the total treatment effect by summing the contributions of each treatment group
  total_random_slope_effect <- rowSums(random_slope)
  
  cov_effects <- rep(0, nrow(data))
  if (length(covariates) > 0) {
    for (cov in covariates) {
      if (cov$type == "continuous") {
        cov_effects <- cov_effects + cov$effect * data[[cov$name]]
      } else {
        for (lvl in cov$levels) {
          if (lvl != cov$reference) {
            cov_effects <- cov_effects + ifelse(data[[cov$name]] == lvl, cov$effects[[lvl]], 0)
          }
        }
      }
    }
  }
  
  # Step 8: Compute post-test scores for non-attrited participants
  data[non_attrition_idx, posts] <- B0 + 
    cov_effects[non_attrition_idx] + 
    total_treatment_effect + 
    b_i + total_random_slope_effect + 
    e_ij
  return(data)
}


# Updated SRT Data Simulation Function
srtdata_simulation <- function(nt, tpi, np, sigma, B0, es, seed, attrition_rates, covariates) {
  # Error checking: ensure tpi has length nt + 1 (control + treatment groups)
  if (length(tpi) != nt + 1) {
    stop("Error: 'tpi' must have length nt + 1 (first value for control group, remaining for treatment groups).")
  }
  
  # Error checking: ensure the sum of tpi is 100
  if (sum(tpi) != 100) {
    stop("Error: The sum of 'tpi' must be 100%.")
  }
  
  # Error checking: ensure es has length equal to nt
  if (length(es) != nt) {
    stop("Error: 'es' must have length nt (one for each treatment group).")
  }
  
  # Error checking: ensure attrition_rates has length equal to nt + 1
  if (length(attrition_rates) != nt + 1) {
    stop("Error: 'attrition_rates' must have length nt + 1 (one for control group and one for each treatment group).")
  }
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Step 1: Create the base data structure with participants
  data <- data.frame(ID = 1:np)  # Assign unique IDs
  
  # Initialize treatment column (0 = control)
  interventions <- "interventions"
  data[[interventions]] <- 0
  
  # Step 2: Assign control group based on tpi[1] (control percentage)
  control_size <- round(np * (tpi[1] / 100))
  control_pupils <- sample(1:np, control_size)
  data[[interventions]][control_pupils] <- 0  # Assign control group (0)
  
  # Step 3: Assign treatment groups
  remaining_pupils <- setdiff(1:np, control_pupils)
  remaining_n <- length(remaining_pupils)
  
  if (nt == 1) {
    # Only one treatment group: assign everyone left
    data[[interventions]][remaining_pupils] <- 1
  } else {
    # Multiple treatment groups
    treatment_props <- tpi[2:(nt+1)] / sum(tpi[2:(nt+1)])  # Normalize proportions
    treatment_sizes <- round(remaining_n * treatment_props)
    
    # Adjust last group to take any rounding error
    treatment_sizes[nt] <- remaining_n - sum(treatment_sizes[1:(nt-1)])
    
    for (i in 1:nt) {
      treated_pupils <- sample(remaining_pupils, treatment_sizes[i])
      data[[interventions]][treated_pupils] <- i
      remaining_pupils <- setdiff(remaining_pupils, treated_pupils)
    }
  }
  
  for (cov in covariates) {
    if (cov$type == "continuous") {
      data[[cov$name]] <- rnorm(nrow(data), mean = 0, sd = cov$sd)
    } else {
      data[[cov$name]] <- sample(cov$levels, nrow(data), replace = TRUE, prob = cov$probs)
    }
  }
  
  # Convert categorical covariates to numeric (0 = reference)
  for (cov in covariates) {
    if (cov$type == "categorical") {
      levels_ordered <- c(cov$reference, setdiff(cov$levels, cov$reference))
      data[[cov$name]] <- match(data[[cov$name]], levels_ordered) - 1
    }
  }
  
  
  # Initialize post-test scores
  posts <- "posttest"
  data[[posts]] <- NA  # Initially set all to NA
  
  # Step 5: Handle attrition for control and treatment groups
  non_attrition_idx <- 1:nrow(data)  # Start by assuming no one is attrited
  
  for (i in 0:nt) {
    group_idx <- which(data[[interventions]] == i)
    attrition_rate <- attrition_rates[i + 1]  # Use i+1 because first value is for control group
    attrition_size <- round(length(group_idx) * attrition_rate)
    
    if (attrition_size > 0) {
      attrition_idx <- sample(group_idx, attrition_size)
      data[attrition_idx, posts] <- NA  # Mark attrited participants
      non_attrition_idx <- setdiff(non_attrition_idx, attrition_idx)  # Update non-attrited participants
    }
  }
  
  # Step 6: Generate post-test scores for non-attrited participants
  
  # Check how many unique intervention groups remain after attrition
  unique_interventions <- unique(data[non_attrition_idx, interventions])
  
  if (length(unique_interventions) > 1) {
    # Create treatment matrix with dummy variables (no intercept)
    treatment_matrix <- model.matrix(~ factor(data[non_attrition_idx, interventions]) - 1)
    
    # Remove the first column (control group is baseline)
    if (ncol(treatment_matrix) >= 1) {
      treatment_matrix <- treatment_matrix[, -1, drop = FALSE]
    }
    
    if (is.null(dim(treatment_matrix))) {
      treatment_matrix <- matrix(treatment_matrix, ncol = 1)
    }
    
    # Apply treatment effects using only individual-level variance (sigma)
    treatment_effects <- sweep(treatment_matrix, 2, es * sigma, `*`)
    
    # Calculate the total treatment effect
    total_treatment_effect <- rowSums(treatment_effects)
  } else {
    # No treatment groups left after attrition
    total_treatment_effect <- rep(0, length(non_attrition_idx))
  }
  
  # Generate individual-level errors for non-attrited individuals
  e_ij <- rnorm(length(non_attrition_idx), mean = 0, sd = sigma)
  
  cov_effects <- rep(0, nrow(data))
  if (length(covariates) > 0) {
    for (cov in covariates) {
      if (cov$type == "continuous") {
        cov_effects <- cov_effects + cov$effect * data[[cov$name]]
      } else {
        for (lvl in cov$levels) {
          if (lvl != cov$reference) {
            cov_effects <- cov_effects + ifelse(data[[cov$name]] == lvl, cov$effects[[lvl]], 0)
          }
        }
      }
    }
  }
  
  # Compute post-test scores
  data[non_attrition_idx, posts] <- B0 + 
    cov_effects[non_attrition_idx] + 
    total_treatment_effect + 
    e_ij
  return(data)
}



crtfutility <- function(data, post_vars = "post1", intervention_column = "treatments", Random = "schls", Nsim = 2000, 
                        Threshold = 0.05, ProbThreshold = 0.8, continuous_covariates = NULL, categorical_covariates = NULL) {
  
  # Ensure categorical covariates are factors
  if (!is.null(categorical_covariates)) {
    for (cat_var in categorical_covariates) {
      data[[cat_var]] <- as.factor(data[[cat_var]])
    }
  }
  
  # Combine covariates AFTER conversion
  covariates <- c(continuous_covariates, categorical_covariates)
  
  
  interventions <- sort(unique(data[[intervention_column]][data[[intervention_column]] != 0]))
  prob_es_values <- numeric(length(interventions))
  futility_decisions <- numeric(length(interventions))  # Initialize a single futility column
  
  output <- list()
  
  for (i in seq_along(interventions)) {
    intervention <- interventions[i]
    intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == 0)
    intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
    
    # Also ensure factor conversion inside the subset
    if (!is.null(categorical_covariates)) {
      for (cat_var in categorical_covariates) {
        intervention_data[[cat_var]] <- as.factor(intervention_data[[cat_var]])
      }
    }
    
    formula_str <- if (is.null(covariates) || length(covariates) == 0) {
      paste(post_vars, "~", intervention_column)
    } else {
      paste(post_vars, "~", intervention_column, "+", paste(covariates, collapse = " + "))
    }
    
    output[[i]] <- crtBayes(
      as.formula(formula_str),
      random = Random,
      intervention = intervention_column,
      nsim = Nsim,
      data = intervention_data,
      threshold = Threshold
    )
    
    prob_es_values[i] <- as.numeric(output[[i]]$ProbES[[1]]["Total1"])
    futility_decisions[i] <- if (prob_es_values[i] < ProbThreshold) 1 else 0
  }
  
  return(data.frame(
    Treatment = interventions,
    Futility = futility_decisions,  # Single futility column
    ProbES = prob_es_values
  ))
}


mstfutility <- function(data, post_vars = "post1", intervention_column = "treatments", Random = "schls", Nsim = 2000, 
                        Threshold = 0.05, ProbThreshold = 0.8, continuous_covariates = NULL, categorical_covariates = NULL) {
  
  # Ensure categorical covariates are factors
  if (!is.null(categorical_covariates)) {
    for (cat_var in categorical_covariates) {
      data[[cat_var]] <- as.factor(data[[cat_var]])
    }
  }
  
  # Combine covariates AFTER conversion
  covariates <- c(continuous_covariates, categorical_covariates)
  
  interventions <- sort(unique(data[[intervention_column]][data[[intervention_column]] != 0]))
  prob_es_values <- numeric(length(interventions))
  futility_decisions <- numeric(length(interventions))  # Initialize a single futility column
  
  output <- list()
  
  for (i in seq_along(interventions)) {
    intervention <- interventions[i]
    intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == 0)
    intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
    
    # Also ensure factor conversion inside the subset
    if (!is.null(categorical_covariates)) {
      for (cat_var in categorical_covariates) {
        intervention_data[[cat_var]] <- as.factor(intervention_data[[cat_var]])
      }
    }
    
    formula_str <- if (is.null(covariates) || length(covariates) == 0) {
      paste(post_vars, "~", intervention_column)
    } else {
      paste(post_vars, "~", intervention_column, "+", paste(covariates, collapse = " + "))
    }
    
    output[[i]] <- mstBayes(
      as.formula(formula_str),
      random = Random,
      intervention = intervention_column,
      nsim = Nsim,
      data = intervention_data,
      threshold = Threshold
    )
    
    prob_es_values[i] <- as.numeric(output[[i]]$ProbES[[1]]["Total1"])
    futility_decisions[i] <- if (prob_es_values[i] < ProbThreshold) 1 else 0
  }
  
  return(data.frame(
    Treatment = interventions,
    Futility = futility_decisions,  # Single futility column
    ProbES = prob_es_values
  ))
}


srtfutility <- function(data, post_vars = "post1", intervention_column = "treatments", Nsim = 2000, 
                        Threshold = 0.05, ProbThreshold = 0.8, continuous_covariates = NULL, categorical_covariates = NULL) {
  
  # Ensure categorical covariates are factors
  if (!is.null(categorical_covariates)) {
    for (cat_var in categorical_covariates) {
      data[[cat_var]] <- as.factor(data[[cat_var]])
    }
  }
  
  # Combine covariates AFTER conversion
  covariates <- c(continuous_covariates, categorical_covariates)
  
  interventions <- sort(unique(data[[intervention_column]][data[[intervention_column]] != 0]))
  prob_es_values <- numeric(length(interventions))
  futility_decisions <- numeric(length(interventions))  # Initialize a single futility column
  
  output <- list()
  
  for (i in seq_along(interventions)) {
    intervention <- interventions[i]
    intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == 0)
    intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
    
    # Also ensure factor conversion inside the subset
    if (!is.null(categorical_covariates)) {
      for (cat_var in categorical_covariates) {
        intervention_data[[cat_var]] <- as.factor(intervention_data[[cat_var]])
      }
    }
    
    formula_str <- if (is.null(covariates) || length(covariates) == 0) {
      paste(post_vars, "~", intervention_column)
    } else {
      paste(post_vars, "~", intervention_column, "+", paste(covariates, collapse = " + "))
    }
    
    output[[i]] <- srtBayes(
      as.formula(formula_str),
      intervention = intervention_column,
      nsim = Nsim,
      data = intervention_data,
      threshold = Threshold
    )
    
    prob_es_values[i] <- as.numeric(output[[i]]$ProbES[[1]]["Cond"])
    futility_decisions[i] <- if (prob_es_values[i] < ProbThreshold) 1 else 0
  }
  
  return(data.frame(
    Treatment = interventions,
    Futility = futility_decisions,  # Single futility column
    ProbES = prob_es_values
  ))
}



crtSuperiority <- function(data, post_var = "post1", intervention_column = "treatments", Random = "schls", 
                           Nsim = 2000, Threshold = 0.1, reference_intervention = 1, superiority_threshold = 0.5, continuous_covariates = NULL, categorical_covariates = NULL) {
  
  # Convert categorical covariates to factors
  if (!is.null(categorical_covariates)) {
    for (cat in categorical_covariates) {
      data[[cat]] <- as.factor(data[[cat]])
    }
  }
  
  covariates <- c(continuous_covariates, categorical_covariates)
  
  
  # Get unique interventions excluding the control group (0)
  interventions <- sort(unique(data[[intervention_column]]))
  interventions <- interventions[interventions != 0 & interventions != reference_intervention]  # Exclude control and reference intervention
  
  output <- list()
  prob_es_values <- numeric(length(interventions))
  sup_decisions <- numeric(length(interventions))
  
  for (i in seq_along(interventions)) {
    intervention <- interventions[i]
    
    # Subset data for the current intervention and reference intervention
    intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == reference_intervention)
    intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
    
    # Ensure factors in subset too
    if (!is.null(categorical_covariates)) {
      for (cat in categorical_covariates) {
        intervention_data[[cat]] <- as.factor(intervention_data[[cat]])
      }
    }
    
    # Build the formula
    formula_str <- if (is.null(covariates) || length(covariates) == 0) {
      paste(post_var, "~", intervention_column)
    } else {
      paste(post_var, "~", intervention_column, "+", paste(covariates, collapse = " + "))
    }
    
    # Run the Bayesian analysis
    output[[as.character(intervention)]] <- crtBayes(
      as.formula(formula_str),
      random = Random,
      intervention = intervention_column,
      nsim = Nsim,
      data = intervention_data,
      threshold = Threshold
    )
    
    # Extract posterior probability of superiority
    prob_es_values[i] <- as.numeric(output[[as.character(intervention)]]$ProbES[[1]]["Total1"])
    sup_decisions[i] <- if (prob_es_values[i] > superiority_threshold) 1 else 0
  }
  
  # Add the reference treatment row
  reference_row <- data.frame(
    Treatment = reference_intervention,
    ProbES = NA,  # No comparison for reference intervention
    Superiority = "Reference"
  )
  
  # Create the result table
  result <- data.frame(
    Treatment = interventions,
    ProbES = prob_es_values,
    Superiority = ifelse(sup_decisions == 1, "Superior", "Not Superior")
  )
  
  # Combine the reference row with the results
  result <- rbind(reference_row, result)
  
  # Sort the table by Treatment
  result <- result[order(result$Treatment), ]
  
  return(result)
}



mstSuperiority <- function(data, post_var = "post1", intervention_column = "treatments", Random = "schls", 
                           Nsim = 2000, Threshold = 0.1, reference_intervention = 1, superiority_threshold = 0.5, continuous_covariates = NULL, categorical_covariates = NULL) {
  
  # Convert categorical covariates to factors
  if (!is.null(categorical_covariates)) {
    for (cat in categorical_covariates) {
      data[[cat]] <- as.factor(data[[cat]])
    }
  }
  
  covariates <- c(continuous_covariates, categorical_covariates)
  
  
  # Get unique interventions excluding the control group (0)
  interventions <- sort(unique(data[[intervention_column]]))
  interventions <- interventions[interventions != 0 & interventions != reference_intervention]  # Exclude control and reference intervention
  
  output <- list()
  prob_es_values <- numeric(length(interventions))
  sup_decisions <- numeric(length(interventions))
  
  for (i in seq_along(interventions)) {
    intervention <- interventions[i]
    
    # Subset data for the current intervention and reference intervention
    intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == reference_intervention)
    intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
    
    # Ensure factors in subset too
    if (!is.null(categorical_covariates)) {
      for (cat in categorical_covariates) {
        intervention_data[[cat]] <- as.factor(intervention_data[[cat]])
      }
    }
    
    # Build the formula
    formula_str <- if (is.null(covariates) || length(covariates) == 0) {
      paste(post_var, "~", intervention_column)
    } else {
      paste(post_var, "~", intervention_column, "+", paste(covariates, collapse = " + "))
    }
    
    # Run the Bayesian analysis
    output[[as.character(intervention)]] <- mstBayes(
      as.formula(formula_str),
      random = Random,
      intervention = intervention_column,
      nsim = Nsim,
      data = intervention_data,
      threshold = Threshold
    )
    
    # Extract posterior probability of superiority
    prob_es_values[i] <- as.numeric(output[[as.character(intervention)]]$ProbES[[1]]["Total1"])
    sup_decisions[i] <- if (prob_es_values[i] > superiority_threshold) 1 else 0
  }
  
  # Add the reference treatment row
  reference_row <- data.frame(
    Treatment = reference_intervention,
    ProbES = NA,  # No comparison for reference intervention
    Superiority = "Reference"
  )
  
  # Create the result table
  result <- data.frame(
    Treatment = interventions,
    ProbES = prob_es_values,
    Superiority = ifelse(sup_decisions == 1, "Superior", "Not Superior")
  )
  
  # Combine the reference row with the results
  result <- rbind(reference_row, result)
  
  # Sort the table by Treatment
  result <- result[order(result$Treatment), ]
  
  return(result)
}


srtSuperiority <- function(data, post_var = "post1", intervention_column = "treatments", 
                           Nsim = 2000, Threshold = 0.1, reference_intervention = 1, superiority_threshold = 0.5, continuous_covariates = NULL, categorical_covariates = NULL) {
  
  # Convert categorical covariates to factors
  if (!is.null(categorical_covariates)) {
    for (cat in categorical_covariates) {
      data[[cat]] <- as.factor(data[[cat]])
    }
  }
  
  covariates <- c(continuous_covariates, categorical_covariates)
  
  # Get unique interventions excluding the control group (0)
  interventions <- sort(unique(data[[intervention_column]]))
  interventions <- interventions[interventions != 0 & interventions != reference_intervention]  # Exclude control and reference intervention
  
  output <- list()
  prob_es_values <- numeric(length(interventions))
  sup_decisions <- numeric(length(interventions))
  
  for (i in seq_along(interventions)) {
    intervention <- interventions[i]
    
    # Subset data for the current intervention and reference intervention
    intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == reference_intervention)
    intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
    
    # Ensure factors in subset too
    if (!is.null(categorical_covariates)) {
      for (cat in categorical_covariates) {
        intervention_data[[cat]] <- as.factor(intervention_data[[cat]])
      }
    }
    
    # Build the formula
    formula_str <- if (is.null(covariates) || length(covariates) == 0) {
      paste(post_var, "~", intervention_column)
    } else {
      paste(post_var, "~", intervention_column, "+", paste(covariates, collapse = " + "))
    }
    
    # Run the Bayesian analysis
    output[[as.character(intervention)]] <- srtBayes(
      as.formula(formula_str),
      intervention = intervention_column,
      nsim = Nsim,
      data = intervention_data,
      threshold = Threshold
    )
    
    # Extract posterior probability of superiority
    prob_es_values[i] <- as.numeric(output[[as.character(intervention)]]$ProbES[[1]]["Cond"])
    sup_decisions[i] <- if (prob_es_values[i] > superiority_threshold) 1 else 0
  }
  
  # Add the reference treatment row
  reference_row <- data.frame(
    Treatment = reference_intervention,
    ProbES = NA,  # No comparison for reference intervention
    Superiority = "Reference"
  )
  
  # Create the result table
  result <- data.frame(
    Treatment = interventions,
    ProbES = prob_es_values,
    Superiority = ifelse(sup_decisions == 1, "Superior", "Not Superior")
  )
  
  # Combine the reference row with the results
  result <- rbind(reference_row, result)
  
  # Sort the table by Treatment
  result <- result[order(result$Treatment), ]
  
  return(result)
}


# Function to add a new treatment to CRT data
add_new_treatmentcrt <- function(existing_data, new_schools, new_pupils_per_school, 
                                 es, attrition_rate, 
                                 post_col, intervention_col, 
                                 schools_col, pupils_col, 
                                 continuous_covariates, categorical_covariates) {
  
  # Check for CRT structure
  school_intervention_check <- aggregate(existing_data[[intervention_col]], 
                                         by = list(existing_data[[schools_col]]), 
                                         FUN = function(x) length(unique(x)))
  if (any(school_intervention_check$x > 1)) {
    stop("Dataset is not a CRT. Each school must have only one unique intervention.")
  }
  
  # Fit LMM
  all_covariates <- c(continuous_covariates, categorical_covariates)
  formula <- as.formula(paste(post_col, "~", paste(all_covariates, collapse = " + "), "+ (1 |", schools_col, ")"))
  lmer_model <- lmer(formula, data = existing_data, REML = TRUE)
  
  fixed_effects <- fixef(lmer_model)
  sigma <- sigma(lmer_model)
  school_sd <- as.numeric(VarCorr(lmer_model)[[schools_col]])^0.5
  
  new_treatment_num <- max(existing_data[[intervention_col]], na.rm = TRUE) + 1
  new_school_ids <- seq(from = max(existing_data[[schools_col]], na.rm = TRUE) + 1, length.out = new_schools)
  
  new_data <- expand.grid(
    schls = new_school_ids,
    ppls = 1:new_pupils_per_school
  )
  
  max_existing_pupil_id <- max(existing_data[[pupils_col]], na.rm = TRUE)
  new_data[[pupils_col]] <- seq(from = max_existing_pupil_id + 1, length.out = nrow(new_data))
  new_data[[schools_col]] <- rep(new_school_ids, each = new_pupils_per_school)
  new_data[[intervention_col]] <- new_treatment_num
  
  # Generate categorical covariates
  for (covariate in categorical_covariates) {
    freq_table <- table(existing_data[[covariate]])
    probs <- freq_table / sum(freq_table)
    new_data[[covariate]] <- factor(sample(names(probs), size = nrow(new_data), replace = TRUE, prob = probs), levels = names(probs))
  }
  
  # Generate continuous covariates
  for (covariate in continuous_covariates) {
    mu <- mean(existing_data[[covariate]], na.rm = TRUE)
    sd_val <- sd(existing_data[[covariate]], na.rm = TRUE)
    new_data[[covariate]] <- rnorm(nrow(new_data), mean = mu, sd = sd_val)
  }
  
  # School-level and individual residuals
  new_data$school_effects <- rep(rnorm(new_schools, mean = 0, sd = school_sd), each = new_pupils_per_school)
  new_data$individual_residuals <- rnorm(nrow(new_data), mean = 0, sd = sigma)
  
  # Baseline post-test score
  new_data[[post_col]] <- fixed_effects[1]
  
  # Add categorical effects (model.matrix)
  if (length(categorical_covariates) > 0) {
    mat <- model.matrix(~ . - 1, data = new_data[, categorical_covariates, drop = FALSE])
    for (col in colnames(mat)) {
      if (col %in% names(fixed_effects)) {
        new_data[[post_col]] <- new_data[[post_col]] + fixed_effects[col] * mat[, col]
      }
    }
  }
  
  # Add continuous effects
  for (cov in continuous_covariates) {
    if (cov %in% names(fixed_effects)) {
      new_data[[post_col]] <- new_data[[post_col]] + fixed_effects[cov] * new_data[[cov]]
    }
  }
  
  # Add treatment effect
  new_data[[post_col]] <- new_data[[post_col]] +
    (new_data[[intervention_col]] == new_treatment_num) * es * sqrt(sigma^2 + school_sd^2) +
    new_data$school_effects + new_data$individual_residuals
  
  # Remove temp columns
  new_data <- new_data[, !names(new_data) %in% c("school_effects", "individual_residuals")]
  
  # Apply attrition
  idx <- sample(seq_len(nrow(new_data)), round(nrow(new_data) * attrition_rate))
  new_data[[post_col]][idx] <- NA
  
  # Align with existing data
  missing_cols <- setdiff(names(existing_data), names(new_data))
  for (col in missing_cols) new_data[[col]] <- NA
  missing_cols2 <- setdiff(names(new_data), names(existing_data))
  for (col in missing_cols2) existing_data[[col]] <- NA
  new_data <- new_data[names(existing_data)]
  
  # Combine
  combined <- rbind(existing_data, new_data)
  
  # Drop any auxiliary columns like schls and ppls if they exist
  combined <- combined[, !names(combined) %in% c("schls", "ppls")]
  
  # Reorder columns (optional)
  key_cols <- c(pupils_col, schools_col, intervention_col, post_col)
  covariate_cols <- setdiff(c(continuous_covariates, categorical_covariates), key_cols)
  final_cols <- c(key_cols, covariate_cols)
  
  # Only include columns that actually exist
  final_cols <- intersect(final_cols, names(combined))
  combined <- combined[, final_cols, drop = FALSE]
  
  return(combined)
  
}




# Function to add a new treatment to MST data
add_new_treatmentmst <- function(existing_data, new_schools, new_pupils_per_school, 
                                 es, attrition_rate, treatment_percentage,
                                 post_col, intervention_col, 
                                 schools_col, pupils_col, 
                                 continuous_covariates, categorical_covariates) {
  
  all_covariates <- c(continuous_covariates, categorical_covariates)
  formula <- as.formula(paste(post_col, "~", paste(all_covariates, collapse = " + "), "+ (1 +", intervention_col, "|", schools_col, ")"))
  model <- lmer(formula, data = existing_data, REML = TRUE)
  
  sigma <- sigma(model)
  rand_sd <- attr(VarCorr(model)[[schools_col]], "stddev")
  sigmab0 <- rand_sd[1]
  sigmab1 <- ifelse(length(rand_sd) > 1, rand_sd[2], 0)
  fixed <- fixef(model)
  
  new_treatment_num <- max(existing_data[[intervention_col]], na.rm = TRUE) + 1
  new_school_ids <- seq(max(existing_data[[schools_col]], na.rm = TRUE) + 1, length.out = new_schools)
  
  new_data <- expand.grid(schls = new_school_ids, ppls = 1:new_pupils_per_school)
  names(new_data)[1] <- schools_col
  new_data[[pupils_col]] <- seq(max(existing_data[[pupils_col]], na.rm = TRUE) + 1, length.out = nrow(new_data))
  
  # Assign intervention
  new_data[[intervention_col]] <- unlist(lapply(seq_along(new_school_ids), function(i) {
    treated <- round(new_pupils_per_school * treatment_percentage)
    c(rep(new_treatment_num, treated), rep(0, new_pupils_per_school - treated))
  }))
  
  # Generate covariates
  for (cov in categorical_covariates) {
    freq <- table(existing_data[[cov]])
    probs <- freq / sum(freq)
    new_data[[cov]] <- factor(sample(names(probs), nrow(new_data), replace = TRUE, prob = probs), levels = names(probs))
  }
  for (cov in continuous_covariates) {
    mu <- mean(existing_data[[cov]], na.rm = TRUE)
    sd_val <- sd(existing_data[[cov]], na.rm = TRUE)
    new_data[[cov]] <- rnorm(nrow(new_data), mu, sd_val)
  }
  
  # Random effects
  school_intercepts <- rnorm(new_schools, 0, sigmab0)
  slopes <- rnorm(new_schools, 0, sigmab1)
  new_data$school_effects <- rep(school_intercepts, each = new_pupils_per_school)
  new_data$treatment_slopes <- rep(slopes, each = new_pupils_per_school)
  new_data$individual_residuals <- rnorm(nrow(new_data), 0, sigma)
  
  new_data[[post_col]] <- fixed[1] + new_data$school_effects + 
    new_data$treatment_slopes * new_data[[intervention_col]] + 
    new_data$individual_residuals + 
    (new_data[[intervention_col]] == new_treatment_num) * es * sqrt(sigma^2 + sigmab0^2 + sigmab1^2)
  
  # Apply attrition
  idx <- sample(seq_len(nrow(new_data)), round(nrow(new_data) * attrition_rate))
  new_data[[post_col]][idx] <- NA
  
  # Final structure
  # Combine datasets
  # Ensure new_data has all columns in existing_data
  missing_in_new <- setdiff(names(existing_data), names(new_data))
  for (col in missing_in_new) {
    new_data[[col]] <- NA
  }
  
  # Ensure existing_data has all columns in new_data
  missing_in_existing <- setdiff(names(new_data), names(existing_data))
  for (col in missing_in_existing) {
    existing_data[[col]] <- NA
  }
  
  # Reorder columns to match
  new_data <- new_data[names(existing_data)]
  
  # Combine safely
  combined <- rbind(existing_data, new_data)
  
  # Remove temp/internal columns
  combined <- combined[, !names(combined) %in% c("school_effects", "treatment_slopes", "individual_residuals", "schls", "ppls")]
  
  # Reorder columns
  key_cols <- c(pupils_col, schools_col, intervention_col, post_col)
  covariate_cols <- setdiff(c(continuous_covariates, categorical_covariates), key_cols)
  final_cols <- intersect(c(key_cols, covariate_cols), names(combined))
  
  combined <- combined[, final_cols, drop = FALSE]
  
  return(combined)
  
}


add_new_treatmentsrt <- function(existing_data, new_pupils, es, attrition_rate, 
                                 post_col, intervention_col, pupils_col, 
                                 continuous_covariates, categorical_covariates) {
  
  sigma <- sd(existing_data[[post_col]], na.rm = TRUE)
  new_treatment_num <- max(existing_data[[intervention_col]], na.rm = TRUE) + 1
  new_ids <- seq(max(existing_data[[pupils_col]], na.rm = TRUE) + 1, length.out = new_pupils)
  
  template <- existing_data[1:new_pupils, , drop = FALSE]
  template[[pupils_col]] <- new_ids
  template[[intervention_col]] <- new_treatment_num
  
  for (cov in categorical_covariates) {
    freq <- table(existing_data[[cov]])
    probs <- freq / sum(freq)
    template[[cov]] <- factor(sample(names(probs), size = new_pupils, replace = TRUE, prob = probs), levels = names(probs))
  }
  
  for (cov in continuous_covariates) {
    mu <- mean(existing_data[[cov]], na.rm = TRUE)
    sd_val <- sd(existing_data[[cov]], na.rm = TRUE)
    template[[cov]] <- rnorm(new_pupils, mu, sd_val)
  }
  
  template[[post_col]] <- mean(existing_data[[post_col]], na.rm = TRUE) + 
    es * sigma + 
    rnorm(new_pupils, 0, sigma)
  
  idx <- sample(seq_len(nrow(template)), round(nrow(template) * attrition_rate))
  template[[post_col]][idx] <- NA
  
  # Combine datasets
  combined <- rbind(existing_data, template)
  
  # Remove internal columns (if any exist)
  combined <- combined[, !names(combined) %in% c("school_effects", "treatment_slopes", "individual_residuals", "schls", "ppls")]
  
  # Reorder columns
  key_cols <- c(pupils_col, intervention_col, post_col)  # no schools_col in SRT
  covariate_cols <- setdiff(c(continuous_covariates, categorical_covariates), key_cols)
  final_cols <- intersect(c(key_cols, covariate_cols), names(combined))
  
  combined <- combined[, final_cols, drop = FALSE]
  
  return(combined)
  
}





# Function to plot futility decision for CRT directly within Shiny
plot_posterior_probabilities_crt <- function(
    data, 
    post_vars = "post1", 
    intervention_column = "treatments", 
    Random = "schls", 
    Nsim = 2000, 
    continuous_covariates = NULL,
    categorical_covariates = NULL, 
    VerticalLine = NULL, 
    VerticalLineColor = "#0000FF",  # New
    ProbThreshold = NULL, 
    HorizontalLineColor = "#FF0000",  # New
    threshold_range = c(0, 1.0),
    plot_title = "Posterior Probabilities Across Thresholds", 
    x_label = "Threshold", 
    y_label = "Posterior Probability", 
    custom_colors = NULL, 
    custom_labels = NULL,
    x_breaks = seq(0, 1, by = 0.1),   # New
    y_breaks = seq(0, 1, by = 0.1)    # New
) {
  
  # Convert categorical covariates to factors
  if (!is.null(categorical_covariates)) {
    for (cat in categorical_covariates) {
      data[[cat]] <- as.factor(data[[cat]])
    }
  }
  
  covariates <- c(continuous_covariates, categorical_covariates)
  
  interventions <- sort(unique(data[[intervention_column]][data[[intervention_column]] != 0]))
  thresholds <- seq(threshold_range[1], threshold_range[2], by = 0.1)
  
  output <- list()
  probabilities_matrix <- matrix(NA, nrow = length(thresholds), ncol = length(interventions))
  
  for (i in seq_along(interventions)) {
    for (j in seq_along(thresholds)) {
      intervention <- interventions[i]
      intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == 0)
      intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
      
      # Apply factor conversion again for safety
      if (!is.null(categorical_covariates)) {
        for (cat in categorical_covariates) {
          intervention_data[[cat]] <- as.factor(intervention_data[[cat]])
        }
      }
      
      formula_str <- if (is.null(covariates) || length(covariates) == 0) {
        paste(post_vars, "~", intervention_column)
      } else {
        paste(post_vars, "~", intervention_column, "+", paste(covariates, collapse = " + "))
      }
      
      output[[j]] <- crtBayes(
        as.formula(formula_str),
        random = Random,
        intervention = intervention_column,
        nsim = Nsim,
        data = intervention_data,
        threshold = thresholds[j]
      )
      
      probabilities_matrix[j, i] <- as.numeric(output[[j]]$ProbES[[1]]["Total1"])
    }
  }
  
  df <- data.frame(
    Threshold = rep(thresholds, length(interventions)),
    Probability = as.vector(probabilities_matrix),
    Intervention = factor(rep(paste("Intervention", 1:length(interventions)), each = length(thresholds)))
  )
  
  # Rename interventions if labels provided
  if (!is.null(custom_labels)) {
    levels(df$Intervention) <- custom_labels
  }
  
  x_breaks <- sort(unique(c(seq(0, 1, by = 0.1), VerticalLine)))
  y_breaks <- sort(unique(c(seq(0, 1, by = 0.1), ProbThreshold)))
  
  # Tick labels (no HTML styling now)
  x_labels <- x_breaks
  y_labels <- y_breaks
  
  p <- ggplot(df, aes(x = Threshold, y = Probability, color = Intervention)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    labs(title = plot_title, x = x_label, y = y_label) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = c(0, 1), breaks = y_breaks, labels = y_labels) +
    theme_pubclean() +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10))
    )
  
  if (!is.null(VerticalLine)) {
    p <- p + geom_vline(xintercept = VerticalLine, linetype = "dashed", color = VerticalLineColor, size = 1)
  }
  
  if (!is.null(ProbThreshold)) {
    p <- p + geom_hline(yintercept = ProbThreshold, linetype = "dashed", color = HorizontalLineColor, size = 1)
  }
  
  if (!is.null(custom_colors)) {
    # Ensure names match the levels in df$Intervention
    p <- p + scale_color_manual(values = custom_colors)
  }
  
  print(p)
}



# Function to plot futility decision for CRT directly within Shiny
plot_posterior_probabilities_mst <- function(
    data, 
    post_vars = "post1", 
    intervention_column = "treatments", 
    Random = "schls", 
    Nsim = 2000, 
    continuous_covariates = NULL,
    categorical_covariates = NULL,
    VerticalLine = NULL, 
    VerticalLineColor = "#0000FF", 
    ProbThreshold = NULL, 
    HorizontalLineColor = "#FF0000", 
    threshold_range = c(0, 1.0),
    plot_title = "Posterior Probabilities Across Thresholds", 
    x_label = "Threshold", 
    y_label = "Posterior Probability", 
    custom_colors = NULL, 
    custom_labels = NULL,
    x_breaks = seq(0, 1, by = 0.1), 
    y_breaks = seq(0, 1, by = 0.1)
) {
  # Convert categorical covariates to factors
  if (!is.null(categorical_covariates)) {
    for (cat in categorical_covariates) {
      data[[cat]] <- as.factor(data[[cat]])
    }
  }
  covariates <- c(continuous_covariates, categorical_covariates)
  
  interventions <- sort(unique(data[[intervention_column]][data[[intervention_column]] != 0]))
  thresholds <- seq(threshold_range[1], threshold_range[2], by = 0.1)
  output <- list()
  probabilities_matrix <- matrix(NA, nrow = length(thresholds), ncol = length(interventions))
  
  for (i in seq_along(interventions)) {
    for (j in seq_along(thresholds)) {
      intervention <- interventions[i]
      intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == 0)
      intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
      
      # Apply factor conversion again for safety
      if (!is.null(categorical_covariates)) {
        for (cat in categorical_covariates) {
          intervention_data[[cat]] <- as.factor(intervention_data[[cat]])
        }
      }
      
      formula_str <- if (is.null(covariates) || length(covariates) == 0) {
        paste(post_vars, "~", intervention_column)
      } else {
        paste(post_vars, "~", intervention_column, "+", paste(covariates, collapse = " + "))
      }
      
      output[[j]] <- mstBayes(
        as.formula(formula_str),
        random = Random,
        intervention = intervention_column,
        nsim = Nsim,
        data = intervention_data,
        threshold = thresholds[j]
      )
      
      probabilities_matrix[j, i] <- as.numeric(output[[j]]$ProbES[[1]]["Total1"])
    }
  }
  
  df <- data.frame(
    Threshold = rep(thresholds, length(interventions)),
    Probability = as.vector(probabilities_matrix),
    Intervention = factor(rep(paste("Intervention", 1:length(interventions)), each = length(thresholds)))
  )
  
  if (!is.null(custom_labels)) {
    levels(df$Intervention) <- custom_labels
  }
  
  x_breaks <- sort(unique(c(seq(0, 1, by = 0.1), VerticalLine)))
  y_breaks <- sort(unique(c(seq(0, 1, by = 0.1), ProbThreshold)))
  
  # Tick labels (no HTML styling now)
  x_labels <- x_breaks
  y_labels <- y_breaks
  
  p <- ggplot(df, aes(x = Threshold, y = Probability, color = Intervention)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    labs(title = plot_title, x = x_label, y = y_label) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(limits = c(0, 1), breaks = y_breaks, labels = y_labels) +
    theme_pubclean() +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10))
    ) +
    coord_cartesian(clip = "off")
  
  if (!is.null(VerticalLine)) {
    p <- p + geom_vline(xintercept = VerticalLine, linetype = "dashed", color = VerticalLineColor, size = 1)
  }
  
  if (!is.null(ProbThreshold)) {
    p <- p + geom_hline(yintercept = ProbThreshold, linetype = "dashed", color = HorizontalLineColor, size = 1)
  }
  
  if (!is.null(custom_colors)) {
    p <- p + scale_color_manual(values = custom_colors)
  }
  
  print(p)
}


# Function to plot futility decision for CRT directly within Shiny
plot_posterior_probabilities_srt <- function(
    data, 
    post_vars = "post1", 
    intervention_column = "treatments", 
    Nsim = 2000, 
    continuous_covariates = NULL,
    categorical_covariates = NULL,
    VerticalLine = NULL, 
    VerticalLineColor = "#0000FF", 
    ProbThreshold = NULL, 
    HorizontalLineColor = "#FF0000", 
    threshold_range = c(0, 1.0),
    plot_title = "Posterior Probabilities Across Thresholds", 
    x_label = "Threshold", 
    y_label = "Posterior Probability", 
    custom_colors = NULL, 
    custom_labels = NULL,
    x_breaks = seq(0, 1, by = 0.1), 
    y_breaks = seq(0, 1, by = 0.1)
) {
  # Convert categorical covariates to factors
  if (!is.null(categorical_covariates)) {
    for (cat in categorical_covariates) {
      data[[cat]] <- as.factor(data[[cat]])
    }
  }
  covariates <- c(continuous_covariates, categorical_covariates)
  
  if (is.null(covariates) || length(covariates) == 0) {
    data$zero_covariate <- 0
    covariates <- "zero_covariate"
  }
  
  interventions <- sort(unique(data[[intervention_column]][data[[intervention_column]] != 0]))
  thresholds <- seq(threshold_range[1], threshold_range[2], by = 0.1)
  output <- list()
  probabilities_matrix <- matrix(NA, nrow = length(thresholds), ncol = length(interventions))
  
  for (i in seq_along(interventions)) {
    for (j in seq_along(thresholds)) {
      intervention <- interventions[i]
      intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == 0)
      intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
      
      # Apply factor conversion again for safety
      if (!is.null(categorical_covariates)) {
        for (cat in categorical_covariates) {
          intervention_data[[cat]] <- as.factor(intervention_data[[cat]])
        }
      }
      
      formula_str <- paste(post_vars, "~", intervention_column)
      if (!is.null(covariates)) {
        formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
      }
      
      interventioncall <- paste0(intervention_column, 1)
      
      output[[j]] <- srtBayes(
        as.formula(formula_str),
        intervention = intervention_column,
        nsim = Nsim,
        data = intervention_data,
        threshold = thresholds[j]
      )
      
      probabilities_matrix[j, i] <- as.numeric(output[[j]]$ProbES[[interventioncall]]["Cond"])
    }
  }
  
  df <- data.frame(
    Threshold = rep(thresholds, length(interventions)),
    Probability = as.vector(probabilities_matrix),
    Intervention = factor(rep(paste("Intervention", 1:length(interventions)), each = length(thresholds)))
  )
  
  if (!is.null(custom_labels)) {
    levels(df$Intervention) <- custom_labels
  }
  
  x_breaks <- sort(unique(c(seq(0, 1, by = 0.1), VerticalLine)))
  y_breaks <- sort(unique(c(seq(0, 1, by = 0.1), ProbThreshold)))
  
  # Tick labels (no HTML styling now)
  x_labels <- x_breaks
  y_labels <- y_breaks  
  
  p <- ggplot(df, aes(x = Threshold, y = Probability, color = Intervention)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    labs(title = plot_title, x = x_label, y = y_label) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(limits = c(0, 1), breaks = y_breaks, labels = y_labels) +
    theme_pubclean() +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10))
    ) +
    coord_cartesian(clip = "off")
  
  if (!is.null(VerticalLine)) {
    p <- p + geom_vline(xintercept = VerticalLine, linetype = "dashed", color = VerticalLineColor, size = 1)
  }
  
  if (!is.null(ProbThreshold)) {
    p <- p + geom_hline(yintercept = ProbThreshold, linetype = "dashed", color = HorizontalLineColor, size = 1)
  }
  
  if (!is.null(custom_colors)) {
    p <- p + scale_color_manual(values = custom_colors)
  }
  
  print(p)
}




run_analysis <- function(data, 
                         post_vars = "post1", 
                         intervention_column = "treatments", 
                         Random = "schls", 
                         Nsim = 2000, 
                         Threshold = 0.05, 
                         ProbThreshold = 0.8, 
                         method = "crtBayes", 
                         crtFREQoption = "Default", 
                         nPerm = NULL, 
                         nBoot = NULL, 
                         bootType = NULL, 
                         continuous_covariates = NULL,
                         categorical_covariates = NULL,
                         input = NULL) {
  
  
  
  # Convert categorical covariates to factors
  if (!is.null(categorical_covariates)) {
    for (cat in categorical_covariates) {
      data[[cat]] <- as.factor(data[[cat]])
    }
  }
  
  covariates <- c(continuous_covariates, categorical_covariates)
  
  interventions <- sort(unique(data[[intervention_column]][data[[intervention_column]] != 0]))
  output <- list()
  
  for (i in seq_along(interventions)) {
    intervention <- interventions[i]
    intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == 0)
    intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
    
    # Convert again after subsetting
    if (!is.null(categorical_covariates)) {
      for (cat in categorical_covariates) {
        intervention_data[[cat]] <- as.factor(intervention_data[[cat]])
      }
    }
    
    # Build formula
    formula_str <- if (length(covariates) == 0) {
      paste(post_vars, "~", intervention_column)
    } else {
      paste(post_vars, "~", intervention_column, "+", paste(covariates, collapse = " + "))
    }
    
    # Fit model based on method
    if (method == "crtBayes") {
      output[[i]] <- crtBayes(
        as.formula(formula_str),
        random = Random,
        intervention = intervention_column,
        nsim = Nsim,
        data = intervention_data,
        threshold = Threshold
      )
      
    } else if (method == "crtFREQ") {
      if (crtFREQoption == "Default") {
        output[[i]] <- crtFREQ(as.formula(formula_str), random = Random, intervention = intervention_column, data = intervention_data)
      } else if (crtFREQoption == "Permutation") {
        output[[i]] <- crtFREQ(as.formula(formula_str), random = Random, intervention = intervention_column, data = intervention_data, nPerm = nPerm)
      } else if (crtFREQoption == "Bootstrap") {
        output[[i]] <- crtFREQ(as.formula(formula_str), random = Random, intervention = intervention_column, data = intervention_data, nBoot = nBoot, type = bootType)
      }
      
    } else if (method == "mstBayes") {
      output[[i]] <- mstBayes(
        as.formula(formula_str),
        random = Random,
        intervention = intervention_column,
        nsim = Nsim,
        data = intervention_data,
        threshold = Threshold
      )
      
    } else if (method == "mstFREQ") {
      if (crtFREQoption == "Default") {
        output[[i]] <- mstFREQ(as.formula(formula_str), random = Random, intervention = intervention_column, data = intervention_data)
      } else if (crtFREQoption == "Permutation") {
        output[[i]] <- mstFREQ(as.formula(formula_str), random = Random, intervention = intervention_column, data = intervention_data, nPerm = nPerm)
      } else if (crtFREQoption == "Bootstrap") {
        output[[i]] <- mstFREQ(as.formula(formula_str), random = Random, intervention = intervention_column, data = intervention_data, nBoot = nBoot, type = bootType)
      }
      
    } else if (method == "srtBayes") {
      output[[i]] <- srtBayes(
        as.formula(formula_str),
        intervention = intervention_column,
        nsim = Nsim,
        data = intervention_data,
        threshold = Threshold
      )
      
    } else if (method == "srtFREQ") {
      if (crtFREQoption == "Default") {
        output[[i]] <- srtFREQ(as.formula(formula_str), intervention = intervention_column, data = intervention_data)
      } else if (crtFREQoption == "Permutation") {
        output[[i]] <- srtFREQ(as.formula(formula_str), intervention = intervention_column, data = intervention_data, nPerm = nPerm)
      } else if (crtFREQoption == "Bootstrap") {
        output[[i]] <- srtFREQ(as.formula(formula_str), intervention = intervention_column, data = intervention_data, nBoot = nBoot)
      }
    }
  }
  
  # Handle plot customizations if enabled
  custom_labels <- if (!is.null(input$enableMultilevelCustomization) && input$enableMultilevelCustomization) {
    interventions <- sort(unique(data[[intervention_column]]))
    interventions <- interventions[interventions != 0]
    sapply(seq_along(interventions), function(i) input[[paste0("label_ml_int_", i)]])
  } else NULL
  
  custom_colors <- if (!is.null(input$enableMultilevelCustomization) && input$enableMultilevelCustomization) {
    interventions <- sort(unique(data[[intervention_column]]))
    interventions <- interventions[interventions != 0]
    sapply(seq_along(interventions), function(i) input[[paste0("color_ml_int_", i)]])
  } else NULL
  
  custom_title <- if (!is.null(input$enableMultilevelCustomization)) input$custom_multilevel_title else NULL
  x_label <- if (!is.null(input$enableMultilevelCustomization)) input$custom_multilevel_xlab else NULL
  y_label <- if (!is.null(input$enableMultilevelCustomization)) input$custom_multilevel_ylab else NULL
  vline_color <- if (!is.null(input$enableMultilevelCustomization)) input$custom_vline_color else "black"
  
  # Generate final comparison plot
  result <- ComparePlot5(
    output,
    modelNames = paste("Intervention", 1:length(interventions)),
    group = 1,
    custom_labels = custom_labels,
    custom_colors = custom_colors,
    custom_title = custom_title,
    x_label = x_label,
    y_label = y_label,
    vline_color = vline_color
  )
  
  return(result)
}




# Define the UI
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .navbar-default .navbar-nav > .active > a,
      .navbar-default .navbar-nav > .active > a:focus,
      .navbar-default .navbar-nav > .active > a:hover {
        border: 2px solid black !important;
        border-radius: 6px;
        background-color: #e6e6e6 !important;
      }
    "))
  ),
  
  # Title Row
  fluidRow(
    column(
      width = 12,
      div(
        style = "text-align: left; margin-bottom: 10px;",  # Reduced margin-bottom for less spacing
        h1("Educational Platform Trials Simulator", 
           style = "font-size: 38px; font-weight: bold; color: #333; margin: 0; padding: 0;")  # Inline style with reduced margins/padding
      )
    )
  ),
  
  # Tabs Row (Placed in a new row)
  fluidRow(
    column(
      width = 12,
      navbarPage(
        "",
        tabPanel(
          "Data Simulation",
          sidebarLayout(
            sidebarPanel(
              radioButtons("simType", "Simulation Type", choices = c("Cluster-Randomised Trial", "Multisite Trial", "Simple Randomised Trial")),
              
              # Number of Interventions - Common to all
              numericInput("nt", "Number of Interventions (excluding control)", 2, min = 1),
              
              # --- SRT Inputs ---
              conditionalPanel(
                condition = "input.simType == 'Simple Randomised Trial'",
                numericInput("np_srt", "Number of Participants", 100, min = 10),
                textInput("tpi_srt", "Percentage of Participants in Each Group (control, Intervention1, ...)", "50,30,20"),
                numericInput("sigma_srt", "Residual Standard Deviation", 1),
                numericInput("B0_srt", "Intercept", 0),
                textInput("es_srt", "Effect Sizes (comma-separated for each Intervention)", "0.2,0.3"),
                                textInput("attrition_rates_srt", "Attrition Rates (comma-separated for control, Intervention1, ...)", "0.1,0.2,0.3"),
                numericInput("seed_srt", "Random Seed", 1234),
                numericInput("num_covariates_srt", "Number of Covariates", value = 0, min = 0, step = 1),
                uiOutput("covariate_inputs_srt")
              ),
              
              # --- CRT Inputs ---
              # Conditional inputs for CRT
              conditionalPanel(
                condition = "input.simType == 'Cluster-Randomised Trial'",
                                numericInput("ns", "Number of Schools", 10, min = 1),
                numericInput("np", "Number of Pupils Per School", 100, min = 10),
                textInput("n_schools_treated", "Number of schools treated (control, Intervention1, ...)", "5,3,2"),
                numericInput("sigma", "Residual Standard Deviation", 1),
                numericInput("ICC", "Intraclass Correlation Coefficient", 0.1),
                numericInput("B0", "Intercept", 0),
                textInput("es", "Effect Sizes (comma-separated for each Intervention)", "0.1,0.2"),
                                textInput("attrition_rates", "Attrition Rates (comma-separated for control, Intervention1, ...)", "0.1,0.1,0.1"),
                numericInput("seed", "Random Seed", 1234),
                
                numericInput("num_covariates_crt", "Number of Covariates", value = 0, min = 0, step = 1),
                uiOutput("covariate_inputs_crt")
              )              ,
              
              # --- MST Inputs ---
              conditionalPanel(
                condition = "input.simType == 'Multisite Trial'",
                                numericInput("ns", "Number of Schools", 10, min = 1),
                numericInput("np", "Number of Pupils Per School", 100, min = 10),
                textInput("tpi", "Percentage of pupils in each group (control, Intervention1, ...)", "50,30,20"),
                numericInput("sigma", "Residual Standard Deviation", 1),
                numericInput("sigmab0", "Random Intercept Standard Deviation", 0.5),
                numericInput("sigmab1", "Random Slope Standard Deviation", 0.5),
                numericInput("B0", "Intercept", 0),
                textInput("es", "Effect Sizes (comma-separated for each Intervention)", "0.2,0.3"),
                textInput("attrition_rates", "Attrition Rates (comma-separated for control, Intervention1, ...)", "0.1,0.1,0.1"),
                                numericInput("seed", "Random Seed", 1234),
                numericInput("num_covariates_mst", "Number of Covariates", value = 0, min = 0, step = 1),
                uiOutput("covariate_inputs_mst")
              ),
              
              actionButton("simulate", "Simulate Data")
            ),
            mainPanel(
              tabsetPanel(
                tabPanel(
                  "Data (first 10 rows)",
                  tableOutput("dataTable"),
                  downloadButton("downloadData", "Download Data")
                )
              )
            )
          )
        ),
        tabPanel(
          "Multi-arm Analysis",
          sidebarLayout(
            sidebarPanel(
              radioButtons("method", "Analysis Method", choices = c("crtBayes", "crtFREQ", "mstBayes", "mstFREQ", "srtBayes", "srtFREQ")),
              radioButtons("mlDataSource", "Data Source",
                           choices = c("Use Simulated Data", "Upload CSV File"),
                           selected = "Upload CSV File"),
              
              conditionalPanel(
                condition = "input.mlDataSource == 'Upload CSV File'",
                fileInput("multilevelDataset", "Choose CSV File", accept = ".csv")
              ),
              uiOutput("post_var_select_multilevel"),
              uiOutput("intervention_var_select_multilevel"),
              uiOutput("random_var_select_multilevel"),
              # Allow selection of multiple covariates
              uiOutput("covariate_cont_select_multilevel"),
              uiOutput("covariate_cat_select_multilevel"),  # Updated to allow multiple selections
              conditionalPanel(
                condition = "input.method == 'crtBayes' || input.method == 'mstBayes' || input.method == 'srtBayes'",
                numericInput("multilevelNsim", "Number of MCMC iterations (>= 10000 is recommended)", 10000, min = 1),
                numericInput("multilevelThreshold", "Effect Size Threshold", 0.05, min = 0, max = 1)
              ),
              conditionalPanel(
                condition = "input.method == 'crtFREQ' || input.method == 'mstFREQ' || input.method == 'srtFREQ'",
                selectInput("crtFREQoption", "Select Confidence Interval Calculation Method", choices = c("Analytic (Default)" = "Default", "Permutation", "Bootstrap")),
                conditionalPanel(
                  condition = "input.crtFREQoption == 'Permutation'",
                  numericInput("nPerm", "Number of Permutations", 1000, min = 1)
                ),
                conditionalPanel(
                  condition = "input.crtFREQoption == 'Bootstrap'",
                  numericInput("nBoot", "Number of Bootstraps", 1000, min = 1)
                ),
                
                # bootType appears only for crtFREQ and mstFREQ when Bootstrap is selected
                conditionalPanel(
                  condition = "input.crtFREQoption == 'Bootstrap' && (input.method == 'crtFREQ' || input.method == 'mstFREQ')",
                  selectInput("bootType", "Bootstrap Type", choices = c(
                    "Re-sampling at student level" = "case(1)",
                    "Re-sampling at school level" = "case(2)",
                    "Re-sampling at both levels" = "case(1,2)",
                    "Residual bootstrapping" = "residual"
                  ))
                )
              ),
              # Remove "Number of Covariates"
              # numericInput("num_covariates", "Number of Covariates", 1, min = 0), -- Removed
              
             
              checkboxInput("enableMultilevelCustomization", "Enable Plot Customization", value = FALSE),
              conditionalPanel(
                condition = "input.enableMultilevelCustomization == true",
                div(
                  style = "margin-top: 10px; padding-top: 10px; border-top: 1px solid #ccc;",
                  
                  textInput("custom_multilevel_title", "Plot Title", "Effect Size Comparison"),
                  textInput("custom_multilevel_xlab", "X-axis Title", "Hedges' g"),
                  textInput("custom_multilevel_ylab", "Y-axis Title", "Models"),
                  colourInput("custom_vline_color", "Vertical Line Color", value = "#000000"),
                  
                  uiOutput("custom_intervention_UI_multilevel"),
                  numericInput("multilevelPlotWidth", "Plot Width (inches)", value = 12, min = 4),
                  numericInput("multilevelPlotHeight", "Plot Height (inches)", value = 9, min = 4)
                )
              ),
              
              actionButton("runAnalysis", "Run Analysis")
            ),
            mainPanel(
              plotOutput("multilevelPlot"),
              
              div(
                style = "display: flex; align-items: flex-end; gap: 8px; justify-content: flex-start; margin-top: 20px;",
                
                # Format dropdown
                div(
                  selectInput(
                    inputId = "multilevelPlotFormat",
                    label = "Download format:",
                    choices = c("TIFF" = "tiff", "PDF" = "pdf", "SVG" = "svg", "EPS" = "eps"),
                    selected = "tiff",
                    width = "120px"
                  )
                ),
                
                # Download button
                div(
                  style = "margin-bottom: 15px;",
                  uiOutput("downloadPlotUI")
                )
              ) ,  
              tags$div(
                style = "margin-top: 30px; color: #555;",
                textOutput("multilevelTimingMessage")
              )
            )
          )
        )  ,
        tabPanel(
          "Futility Analysis",
          sidebarLayout(
            sidebarPanel(
              radioButtons("futSimType", "Simulation Type", choices = c("Cluster-Randomised Trial", "Multisite Trial", "Simple Randomised Trial")),
              radioButtons("futDataSource", "Data Source",
                           choices = c("Use Simulated Data", "Upload CSV File"),
                           selected = "Upload CSV File"),
              
              # ✅ This condition now works because futDataSource exists
              conditionalPanel(
                condition = "input.futDataSource == 'Upload CSV File'",
                fileInput("dataset", "Choose CSV File", accept = ".csv")
              ),
              uiOutput("post_var_select"),
              uiOutput("intervention_var_select"),
              conditionalPanel(
                condition = "input.futSimType != 'Simple Randomised Trial'",
                uiOutput("random_var_select")
              ),
              uiOutput("covariate_cont_select_fut"),
              uiOutput("covariate_cat_select_fut"),  # Updated to allow multiple selections
              numericInput("crtNsim", "Number of MCMC iterations (>= 10000 is recommended)", 10000, min = 1),
              numericInput("crtThreshold", "Effect Size Threshold", 0.05, min = 0, max = 1),
              numericInput("crtProbThreshold", "Futility Threshold", 0.8, min = 0, max = 1),
              actionButton("analyzeFutility", "Analyze Futility")
            ),
            mainPanel(
              tableOutput("futilityTable"),
              
              # Place the download button here
              downloadButton("downloadFutilityData", "Download Futility Results"),
              tags$div(
                style = "margin-top: 30px; color: #555;",
                textOutput("futilityTimingMessage")
              )
            )
          )
        ),
        tabPanel(
          "Superiority Analysis",
          sidebarLayout(
            sidebarPanel(
              radioButtons("supSimType", "Simulation Type", choices = c("Cluster-Randomised Trial", "Multisite Trial", "Simple Randomised Trial")),
              radioButtons("supDataSource", "Data Source",
                           choices = c("Use Simulated Data", "Upload CSV File"),
                           selected = "Upload CSV File"),
              
              conditionalPanel(
                condition = "input.supDataSource == 'Upload CSV File'",
                fileInput("supDataset", "Choose CSV File", accept = ".csv")
              ),
              
              uiOutput("post_var_select_sup"),      # Post-test variable (single select)
              uiOutput("intervention_var_select_sup"),  # Intervention variable (single select)
              
              conditionalPanel(
                condition = "input.supSimType != 'Simple Randomised Trial'",
                uiOutput("random_var_select_sup")
              ),
              uiOutput("covariate_cont_select_sup"),
              uiOutput("covariate_cat_select_sup"),  # Covariates (multi-select)
              uiOutput("reference_intervention_select"),
              
              numericInput("crtSupNsim", "Number of MCMC iterations (>= 10000 is recommended)", 10000, min = 1),
              numericInput("crtSupThreshold", "Effect Size Threshold", 0.05, min = 0, max = 1),
              
              # New Superiority Threshold input
              numericInput("crtSupSuperiorThreshold", "Superiority Threshold", 0.8, min = 0, max = 1),
              
              
              
              actionButton("analyzeSuperiority", "Analyze Superiority")
            ),
            mainPanel(
              tableOutput("superiorityTable"),
              
              # Place the download button here
              downloadButton("downloadSuperiorityData", "Download Superiority Results"),
              tags$div(
                style = "margin-top: 30px; color: #555;",
                textOutput("superiorityTimingMessage")
              )
            )
          )
        ),
        tabPanel(
          "Add New Intervention",
          sidebarLayout(
            sidebarPanel(
              # Choose CRT, MST, or SRT option
              radioButtons("newSimType", "Simulation Type", choices = c("Cluster-Randomised Trial", "Multisite Trial", "Simple Randomised Trial")),
              
              # File input to upload the existing dataset
              radioButtons("newDataSource", "Data Source",
                           choices = c("Use Simulated Data", "Upload CSV File"),
                           selected = "Upload CSV File"),
              
              conditionalPanel(
                condition = "input.newDataSource == 'Upload CSV File'",
                fileInput("newTreatmentDataset", "Choose Existing Dataset (CSV)", accept = ".csv")
              ),
              
              # Show only for CRT and MST (Not for SRT)
              conditionalPanel(
                condition = "input.newSimType != 'Simple Randomised Trial'",
                numericInput("newSchools", "Number of New Schools", 2, min = 1),
                numericInput("newPupilsPerSchool", "Pupils per New School", 100, min = 10)
              ),
              
              # Common Inputs for all types (CRT, MST, SRT)
              numericInput("newEffectSize", "Effect Size for New Intervention", 0.3),
              numericInput("newAttrition", "Attrition Rate for New Intervention", 0.1, min = 0, max = 1),
              
              # MST-specific input
              conditionalPanel(
                condition = "input.newSimType == 'Multisite Trial'",
                sliderInput("treatmentPercentageMST", "Percentage of Pupils in New Intervention (MST)", 
                            min = 0, max = 1, value = 0.5, step = 0.1)
              ),
              
              # User selections for columns
              uiOutput("post_var_select_new"),
              uiOutput("intervention_var_select_new"),
              uiOutput("schools_var_select_new"),
              uiOutput("pupils_var_select_new"),
              
              # Covariate selection
              uiOutput("covariate_cont_select_new"),
              uiOutput("covariate_cat_select_new"),
              
              # Button to initiate the addition of a new treatment
              actionButton("addNewTreatment", "Add New Intervention")
            ),
            
            mainPanel(
              # Display the modified dataset
              tableOutput("newTreatmentTable"),
              downloadButton("downloadNewTreatmentData", "Download New Dataset")
            )
          )
        ),
        
        tabPanel(
          "Plot Posterior Probabilities",
          sidebarLayout(
            sidebarPanel(
              # --- Core Inputs ---
              radioButtons("plotSimType", "Simulation Type", choices = c("Cluster-Randomised Trial", "Multisite Trial", "Simple Randomised Trial")),
              radioButtons("plotDataSource", "Data Source",
                           choices = c("Use Simulated Data", "Upload CSV File"),
                           selected = "Upload CSV File"),
              
              conditionalPanel(
                condition = "input.plotDataSource == 'Upload CSV File'",
                fileInput("plotDataset", "Choose CSV File", accept = ".csv")
              ),
              
              uiOutput("post_var_select_plot"),
              uiOutput("intervention_var_select_plot"),
              
              conditionalPanel(
                condition = "input.plotSimType != 'Simple Randomised Trial'",
                uiOutput("random_var_select_plot")
              ),
              
              uiOutput("covariate_cont_select_plot"),
              uiOutput("covariate_cat_select_plot"),
              
              numericInput("crtPlotNsim", "Number of MCMC iterations (>= 10000 is recommended)", 10000, min = 1),
              
              checkboxInput("addVerticalLine", "Add a vertical line", value = FALSE),
              conditionalPanel(
                condition = "input.addVerticalLine == true",
                numericInput("crtVerticalLine", "Value for Vertical Line", value = 0.05, min = 0, max = 1)
              ),
              
              checkboxInput("addHorizontalLine", "Add a horizontal line", value = FALSE),
              conditionalPanel(
                condition = "input.addHorizontalLine == true",
                numericInput("crtPlotProbThreshold", "Value for Horizontal Line", value = 0.8, min = 0, max = 1)
              ),
              
              sliderInput("plot_threshold_range", "Plot Threshold Range", min = 0, max = 1, value = c(0, 1), step = 0.1),
              
              # ✅ NEW: Enable Plot Customization
              checkboxInput("enablePlotCustomization", "Enable Plot Customization", value = FALSE),
              
              # ✅ NEW: Shown only when checkbox is enabled
              conditionalPanel(
                condition = "input.enablePlotCustomization == true",
                div(
                  style = "margin-top: 10px; padding-top: 10px; border-top: 1px solid #ccc;",
                  textInput("custom_plot_title", "Main Plot Title", value = "Posterior Probabilities Across Thresholds"),
                  textInput("custom_x_label", "X-axis Title", value = "Threshold"),
                  textInput("custom_y_label", "Y-axis Title", value = "Posterior Probability"),
                  
                  sliderInput("x_breaks_custom", "X-axis Scale (Tick Interval)", min = 0.01, max = 0.5, value = 0.1, step = 0.01),
                  sliderInput("y_breaks_custom", "Y-axis Scale (Tick Interval)", min = 0.01, max = 0.5, value = 0.1, step = 0.01),
                  
                  colourInput("verticalLineColor", "Vertical Line Color", value = "#0000FF"),
                  colourInput("horizontalLineColor", "Horizontal Line Color", value = "#FF0000"),
                  
                  uiOutput("interventionCustomUI"),
                  # Plot size inputs
                  numericInput("posteriorPlotWidth", "Plot Width (inches)", value = 7, min = 4),
                  numericInput("posteriorPlotHeight", "Plot Height (inches)", value = 5, min = 4),
                )
              ),
              
              actionButton("plotPosterior", "Plot Posterior Probabilities")
            ),
            
            mainPanel(
              plotOutput("posteriorPlot"),
              
              div(
                style = "display: flex; align-items: flex-end; gap: 8px; justify-content: flex-start; margin-top: 20px;",
                div(
                  selectInput(
                    inputId = "posteriorPlotFormat",
                    label = "Download format:",
                    choices = c("TIFF" = "tiff", "PDF" = "pdf", "SVG" = "svg", "EPS" = "eps"),
                    selected = "tiff",
                    width = "120px"
                  )
                ),
                div(
                  style = "margin-bottom: 15px;",
                  uiOutput("downloadPosteriorPlotUI")
                )
              ),
              
              tags$div(
                style = "margin-top: 30px; color: #555;",
                textOutput("posteriorPlotTimingMessage")
              )
            )
          )
        ),
        
        # Adding the "User Manual" tab
        tabPanel(
          "User Manual",
          fluidPage(
            tags$head(
              tags$style(HTML("
          /* Increase font size for the entire User Manual tab */
          #user-manual-content {
            font-size: 16px;  /* Increase the font size */
            line-height: 1.8;  /* Adjust line height for readability */
          }
          h1 {
            font-size: 32px;  /* Increase font size for main headers */
          }
          h2 {
            font-size: 26px;  /* Increase font size for section headers */
          }
          ul {
            margin-bottom: 17px;  /* Add some spacing between lists */
          }
        "))
            ),
            div(
              id = "user-manual-content",  # Apply custom styles to this div
              HTML("
          <h1>Educational Platform Trials Simulator (EPTS) Manual</h1>
          <h2>Summary</h2>
          <p>The EPTS simulator is a web application developed using RShiny, a package within the R and RStudio statistical software environment. The software application includes functions for cluster-randomised, multisite, and simple randomised platform trial simulation, multilevel analysis, as well as futility and superiority analysis outputs. The software's rules are based on calculating Bayesian posterior probabilities of superiority and futility. To run the simulator, you need to input some data into the input bar manually. The software enables you to save and load simulation outputs. For all tabs, the input bar is located on the left side of the browser window, and the outputs are displayed on the right side. EPTS consists of seven tabs: data simulation, multilevel analysis, futility analysis, superiority analysis, add new intervention, plot posterior probabilities and User Manual.</p>

          <h2>Data Simulation</h2>
          <p>The Data Simulation tab allows the user to simulate the Cluster-Randomised Trial (CRT), Multisite Trial (MST) and Simple Randomised Trial (SRT).The simulation requires the user to input the following parameters: </p>
          <ul>
            <li><strong>Simulation Type (Cluster-Randomised Trial, Multisite Trial or Simple Randomised Trial):</strong> Select the trial design by choosing from Cluster-Randomised, Multisite or Simple Randomised Trial.</li>
            <li><strong>Number of Interventions:</strong> The number of intervention arms in the trial, excluding the control arm. <br/>Example: 2 (1 control group and 2 intervention arms)</li>
            <li><strong>Number of Schools (for CRT/MST):</strong> Total number of schools (or clusters) participating in the trial.</li>
            <li><strong>Number of Pupils Per School (for CRT/MST):</strong> The total number of pupils (or individuals) assigned to each school (or cluster).</li>
            <li><strong>Number of Schools Treated (for CRT):</strong> A comma-separated list indicating the number of schools assigned to each group (control, Intervention 1, Intervention 2, etc.). This input applies only to CRT simulation.<br/>Example: '5, 3, 2' (5 control schools, 3 in Intervention 1, 2 in Intervention 2).</li>
            <li><strong>Percentage of Pupils in Each Group (for MST):</strong> A comma-separated list specifying the percentage of pupils in each group (control, Intervention 1, Intervention 2, etc.). This input applies only to MST simulation.<br/>Example: '50, 30, 20' (50% control, 30% in Intervention 1, 20% in Intervention 2).</li>
            <li><strong>Residual Standard Deviation:</strong> The standard deviation of residual errors.</li>
            <li><strong>Intraclass Correlation Coefficient (for CRT):</strong> The ICC measures the proportion of total variance attributed to differences between clusters (e.g., schools).</li>
            <li><strong>Random Intercept Standard Deviations (for MST):</strong> In MST simulations, it specifies the standard deviations for the random intercept, quantifying baseline heterogeneity between schools.</li>
            <li><strong>Random Slope Standard Deviations (for MST):</strong> In MST simulations, users must specify the standard deviations for the random slope, quantifying differential effects of the intervention across schools through school-by-intervention interactions.</li>
            <li><strong>Standard Deviation of Pre-test Scores:</strong> The standard deviation of pre-test scores for individuals.</li>
            <li><strong>Intercept (B0):</strong> The intercept (B0) represents the regression coefficient for intercept.</li>
            <li><strong>Effect Sizes (es):</strong> A comma-separated list of effect sizes for each Intervention group. Example: '0.2, 0.3' (0.2 effect size for Intervention 1, 0.3 for Intervention 2).</li>
            <li><strong>Random Seed:</strong> A seed value used to ensure reproducibility of the simulation.</li>
            <li><strong>Attrition Rates:</strong> A comma-separated list specifying the attrition rate for the control group and each Intervention group. Attrition refers to participants dropping out of the trial.<br/>Example: '0.1, 0.2, 0.3' (10% attrition in control, 20% in Intervention 1, and 30% in Intervention 2).</li>
          
          <li><strong>Number of Covariates:</strong> The total number of covariates (independent variables) included in the model.</li>
<li><strong>Covariate Name:</strong> A label identifying the first covariate. <br/>Example: 'pretest'</li>
<li><strong>Type:</strong> Indicates whether the covariate is continuous or categorical. <br/>Example: 'continuous' or 'categorical'</li>
<li><strong>Standard Deviation (for continuous):</strong> The standard deviation of the continuous covariate. <br/>Example: 1</li>
<li><strong>Coefficient (for continuous):</strong> The coefficient of the continuous covariate in the model. <br/>Example: 0.5</li>
<li><strong>Levels (for categorical):</strong> A comma-separated list of the possible levels for a categorical covariate. <br/>Example: 'A,B,C'</li>
<li><strong>Probabilities (for categorical):</strong> A comma-separated list of the probabilities for each level. Must sum to 1. <br/>Example: '0.3,0.3,0.4'</li>
<li><strong>Reference Category (for categorical):</strong> The reference category used for comparison in the model. <br/>Example: 'A'</li>
<li><strong>Coefficient for category (for categorical):</strong> The coefficients for categories excluding the reference category. <br/>Example: 0.2</li>

          </ul>
          
          <h2>Multi-arm Analysis Tab</h2>
          <p>The Multilevel Analysis tab allows the user to analyze the CRT, MST, and SRT datasets using Frequentist and Bayesian multilevel models.</p>
          <ul>
            <li><strong>Analysis Method:</strong> Select the type of analysis to be performed. Available options include:</li>
            <li><strong>Use Simulated Data:</strong> Select this option to use simulated data generated from Data Simulaion tab. <br/>No file upload is required.</li>
<li><strong>Upload CSV File:</strong> To begin the analysis, upload a CSV file containing your trial data. <br/>The dataset should include post-test outcomes, intervention group assignments, and any covariates you wish to include in the analysis.</li>

            <ul>
              <li>crtBayes: Bayesian analysis of cluster-randomized trials using vague priors.</li>
              <li>crtFREQ: Frequentist analysis of cluster-randomized trials using multilevel models.</li>
              <li>mstBayes: Bayesian analysis of multisite randomized trials using vague priors.</li>
              <li>mstFREQ: Frequentist analysis of multisite randomized trials using multilevel models.</li>
              <li>srtBayes: Bayesian analysis of simple randomized trials using vague priors.</li>
              <li>srtFREQ: Frequentist analysis of simple randomized trials.</li>
            </ul>
            <li><strong>Select Post-test Outcome:</strong> Select the variable from your dataset that represents the outcome or dependent variable.</li>
            <li><strong>Select Intervention Variables:</strong> Select the variable from your dataset that represents the intervention or Intervention groups.</li>
            <li><strong>Select Clustering Variable (for CRT/MST):</strong> Select the variable from your dataset that represents the clustering structure (e.g., schools, hospitals). This applies only for multilevel models (CRT and MST).</li>
            <li><strong>Number of MCMC iterations (Bayesian Methods):</strong> Set the number of MCMC iterations per chain (>= 10000 is recommended). A minimum of 10,000 is recommended
to ensure convergence.</li>
            <li><strong>Effect Size Threshold:</strong> A scalar pre-specified threshold for estimating Bayesian posterior probability that the observed effect size is greater than or equal to the threshold.</li>
            <li><strong>Frequentist Method Options (crtFREQ, mstFREQ, srtFREQ):</strong></li>
            <li><strong>Select Confidence Interval Calculation Method:</strong></li>
            <ul> 
              <li>Analytic (Default): Perform standard frequentist analysis.</li>
              <li>Permutation: Specify the number of permutations required to generate a permuted p-value.</li>
              <li>Bootstrap: Perform bootstrap analysis by specifying the number of bootstraps required to generate bootstrap confidence intervals and the bootstrap type.</li>
            </ul>
            <li><strong>Select Continuous Covariates:</strong> Choose one or more continuous covariates from the dataset that should be included in the model <br/>Example: 'pretest', 'age'</li>
<li><strong>Select Categorical Covariates:</strong> Choose one or more categorical covariates from the dataset that should be included in the model. <br/>Example: 'gender', 'ethnicity'</li>
     
     <li><strong>Enable Plot Customization:</strong> Toggle this option to customize the appearance of the effect size comparison plot. Example: enables title, axis labels, colors, and dimensions.</li>
<li><strong>Plot Title:</strong> The title displayed at the top of the plot. Example: 'Effect Size Comparison'</li>
<li><strong>X-axis Title:</strong> The label for the X-axis of the plot. Example: 'Hedges' g'</li>
<li><strong>Y-axis Title:</strong> The label for the Y-axis of the plot. Example: 'Interventions'</li>
<li><strong>Vertical Line Color:</strong> The color used for the reference line (e.g., at zero) in the plot. Example: '#000000'</li>
<li><strong>Plot Width (inches):</strong> The total width of the plot in inches. Example: 12</li>
<li><strong>Plot Height (inches):</strong> The total height of the plot in inches. Example: 9</li>

           <li><strong>Run Analysis:</strong> After specifying all necessary parameters, click the 'Run Analysis' button to perform the chosen analysis method.</li>
            <li><strong>Results and Plot:</strong> Once the analysis is complete, a forest plot will be created to display the effect sizes using within and total variances. The plot can be downloaded in TIFF, PDF, SVG, or EPS formats.</li>
          </ul>
          
          <h2>Futility Analysis Tab</h2>
          <p>The Futility Analysis tab allows users to assess futility in platform trials. Users must either use simulated data or upload a CSV dataset, then configure the simulation settings.</p>
          <ul>
            <li><strong>Simulation Type (CRT, MST, SRT):</strong> Select the type of trial simulation being analyzed for futility: CRT, MST, or SRT.</li>
            <li><strong>Effect Size Threshold:</strong> The threshold for estimating the Bayesian Posterior Probability.</li>
            <li><strong>Futility Threshold:</strong> The threshold for the Bayesian Posterior Probability below which an intervention is considered futile.</li>
<li><strong>Select Continuous Covariates:</strong> Choose one or more continuous covariates from the dataset that should be included in the model <br/>Example: 'pretest', 'age'</li>
<li><strong>Select Categorical Covariates:</strong> Choose one or more categorical covariates from the dataset that should be included in the model. <br/>Example: 'gender', 'ethnicity'</li>            <li><strong>Analyze Futility:</strong> Click the 'Analyze Futility' button to run the futility analysis based on the selected parameters.</li>
            <li><strong>Futility Decision:</strong> The futility analysis results will present the probability that an intervention's effect size exceeds a specified threshold, denoted as P(Effect Size > Threshold). Additionally, the analysis will determine and indicate whether each intervention is considered futile based on futility threshold.</li>
          </ul>

          <h2>Superiority Analysis Tab</h2>
          <p>The Superiority analysis tab is designed to compare the efficacy of a an Intervention against a reference Intervention. Users must select simulated data or upload a CSV dataset.</p>
          <ul>
           <li><strong>Effect Size Threshold:</strong> The threshold for estimating the Bayesian Posterior Probability.</li>
            <li><strong>Superiority Threshold:</strong> The threshold for the Bayesian Posterior Probability below which an intervention is considered futile.</li>
            
            <li><strong>Reference Intervention:</strong> Select the reference Intervention against which all other Interventions will be compared.</li>
<li><strong>Select Continuous Covariates:</strong> Choose one or more continuous covariates from the dataset that should be included in the model <br/>Example: 'pretest', 'age'</li>
<li><strong>Select Categorical Covariates:</strong> Choose one or more categorical covariates from the dataset that should be included in the model. <br/>Example: 'gender', 'ethnicity'</li>            <li><strong>Analyze Superiority:</strong> Click the 'Analyze Superiority' button to run the superiority analysis.</li>
            <li><strong>Superiority Decision:</strong> The superiority analysis results will present the probability that an intervention's effect size exceeds that of the reference intervention, denoted as P(Effect Size > Threshold). Additionally, the analysis will determine and indicate whether each intervention is considered superior based on superiority threshold.</li>
          </ul>
          
          <h2>Add New Intervention</h2>
        <p>The Add New Intervention tab allows users to introduce a new intervention into an existing dataset. Users must either use simulated data or upload a CSV dataset, then select the type of data (CRT, MST and SRT):</p>

        <li><strong>select required variables:</strong></li>
          <ul>
              <li><strong>Select Post-test Outcome:</strong> The outcome variable for analysis.</li>
              <li><strong>Select Intervention Variable:</strong> The column containing treatment assignments.</li>
              <li><strong>Select Pupils Variable (for CRT/MST):</strong> The pupil level identifier.</li>
             <li><strong>Select ID Variable (for SRT):</strong> The individual-level identifier.</li>

              <li><strong>Select Schools Variable (for CRT/MST):</strong> The cluster-level identifier..</li>
<li><strong>Select Continuous Covariates:</strong> Choose one or more continuous covariates from the dataset that should be included in the model <br/>Example: 'pretest', 'age'</li>
<li><strong>Select Categorical Covariates:</strong> Choose one or more categorical covariates from the dataset that should be included in the model. <br/>Example: 'gender', 'ethnicity'</li>            </ul>
          </li>
          <li><strong>Specify New Intervention Parameters:</strong></li>
          <ul>
            <li><strong>Number of New Schools (for CRT/MST):</strong> Specify how many new schools will be introduced.</li>
            <li><strong>Number of Pupils Per School (for CRT/MST):</strong> The number of pupils per new school.</li>
             <li><strong>Number of Participants (for SRT):</strong> The number of participants in new intervention.</li>

            <li><strong>Effect Size for New Intervention:</strong> The effect size of the new intervention.</li>
            <li><strong>Attrition Rate for New Intervention:</strong> The proportion of participants expected to drop out of the new intervention. Default: 0.1 (10% attrition in outcome)</li>
            <li><strong>Percentage of Pupils in New Treatment (for MST):</strong> The proportion of students in new intervention groups.</li>
          </ul>

          
          
          
          <h2>Plot Posterior Probabilities Tab</h2>
          <p>This tab allows users to visualize posterior probabilities across different thresholds for multiple interventions.</p>
          <ul>
            <li><strong>Add Vertical Line:</strong> Option to add a vertical reference line to the plot, typically representing a pre-specified threshold for estimating Bayesian posterior probability.</li>
            <li><strong>Add Horizontal Line:</strong> Option to add a horizontal reference line to the plot, typically representing a threshold of Bayesian posterior probability.</li>
            <li><strong>Plot Threshold Range:</strong> Adjust the range of threshold values for the posterior probability plot.</li>
            <li><strong>Plot Posterior Probabilities:</strong> Click the 'Plot Posterior Probabilities' button to generate the plot. The generated plot will display posterior probabilities across thresholds for each intervention group. You can also download the plot in TIFF format.</li>
          </ul>
          
          <h2>Plot Posterior Probabilities Tab</h2>
  <p>This tab allows users to visualize posterior probabilities across different thresholds for multiple interventions.</p>
  <ul>
    <li><strong>Add Vertical Line:</strong> Option to add a vertical reference line to the plot, typically representing a pre-specified threshold for estimating Bayesian posterior probability.</li>
    <li><strong>Add Horizontal Line:</strong> Option to add a horizontal reference line to the plot, typically representing a threshold of Bayesian posterior probability.</li>
    <li><strong>Plot Threshold Range:</strong> Adjust the range of threshold values for the posterior probability plot.</li>
    <li><strong>Plot Posterior Probabilities:</strong> Click the 'Plot Posterior Probabilities' button to generate the plot. The generated plot will display posterior probabilities across thresholds for each intervention group. You can also download the plot in TIFF, PDF, SVG, or EPS formats.</li>

 <li><strong>Enable Plot Customization:</strong> Toggle this option to customize the appearance of the output plot. Includes titles, axis labels, colors, and dimensions.</li>
<li><strong>Main Plot Title:</strong> The title displayed at the top of the plot. Example: 'Posterior Probabilities Across Thresholds'</li>
<li><strong>X-axis Title:</strong> The label for the X-axis of the plot. Example: 'Threshold'</li>
<li><strong>Y-axis Title:</strong> The label for the Y-axis of the plot. Example: 'Posterior Probability'</li>
<li><strong>X-axis Scale (Tick Interval):</strong> The spacing between tick marks on the x-axis. Example: 0.1</li>
<li><strong>Y-axis Scale (Tick Interval):</strong> The spacing between tick marks on the y-axis. Example: 0.05</li>
<li><strong>Vertical Line Color:</strong> The color used for vertical reference line in the plot. Example: '#0000FF'</li>
<li><strong>Horizontal Line Color:</strong> The color used for horizontal reference line in the plot. Example: '#FF0000'</li>
<li><strong>Rename Intervention:</strong> Custom label for Interventions displayed in the plot legend. Example: 'Intervention 1'</li>
<li><strong>Color for Intervention:</strong> Color used to represent Interventions in the plot. Example: '#1F77B4'</li>
<li><strong>Plot Width (inches):</strong> The total width of the plot in inches. Example: 7</li>
<li><strong>Plot Height (inches):</strong> The total height of the plot in inches. Example: 5</li>

  </ul>

  <h2>Performance Considerations</h2>
  <p>The computational time required to simulate data or run analyses depends on several factors:</p>
  <ul>
    <li><strong>Number of Individuals:</strong> Larger sample sizes increase computation time, particularly in multilevel models.</li>
    <li><strong>Number of Schools (Clusters):</strong> More clusters add complexity to hierarchical models like CRTs and MSTs.</li>
    <li><strong>Number of Interventions:</strong> More intervention arms increase simulation and model-fitting time.</li>
    <li><strong>Model Type:</strong> Bayesian models typically take longer than frequentist ones. CRTs and MSTs are slower than SRTs due to their multilevel structure.</li>
    <li><strong>Number of MCMC Iterations (Bayesian):</strong> In Bayesian analyses, more iterations increase computational time.</li>
      <li><strong>Number of Permutations and Bootstraps (Frequentist):</strong> In frequentist methods, a higher number of resampling iterations increases processing time.</li>
  </ul>

  <h4>Estimated Performance Benchmarks</h4>
 <p>Here are example scenarios to help you gauge expected computation times for multilevel (Bayesian), futility, and superiority analysis:</p>
  </ul>
    <li> CRT with 10 schools, 1000 pupils per school, 3 arms, Nsim = 10,000 → ~40–50 seconds</li>
    <li> MST with 10 schools, 100 pupils per school, 3 arms, Nsim = 10,000 → ~40–60 seconds</li>
    <li> SRT with 3 arms, 1000 participants, Nsim = 10,000 → ~5–10 seconds</li>
  </ul>

  <p><strong>Note:</strong> Posterior probability plotting can take longer when using a wider <em>Plot Threshold Range</em>. Data simulation and adding a new treatment typically take a few seconds.A progress spinner is shown during computations for better user experience. After an analysis is complete, the total duration is displayed for reference.</p>
  

          <h2>Associated Paper</h2>
          <p>The associated paper of this R Shiny application can be accessed via: <a href='#'>Link to Paper</a>.</p>
     ")
            )
          )
        )
      )
    )
  )
)

# Define the server
server <- function(input, output, session) {
  simulatedData <- eventReactive(input$simulate, {
    tryCatch({
      nt <- input$nt
      covariate_spec <- switch(input$simType,
                               "Cluster-Randomised Trial" = extract_covariates("crt_cov", input$num_covariates_crt),
                               "Multisite Trial" = extract_covariates("mst_cov", input$num_covariates_mst),
                               "Simple Randomised Trial" = extract_covariates("srt_cov", input$num_covariates_srt)
      )
      
      
      if (input$simType == "Cluster-Randomised Trial") {
        np <- input$np
        ns <- input$ns
        n_schools_treated <- as.numeric(unlist(strsplit(as.character(input$n_schools_treated), ",")))
        
        # ✅ Validation: Ensure the number of treated schools matches the total number of schools
        if (sum(n_schools_treated) != ns) {
          stop("❌ The total number of treated schools must equal the total number of schools.")
        }
        
        sigma <- input$sigma
        ICC <- input$ICC
        B0 <- input$B0
        es <- as.numeric(unlist(strsplit(as.character(input$es), ",")))
        seed <- input$seed
        attrition_rates <- as.numeric(unlist(strsplit(as.character(input$attrition_rates), ",")))
        
        crtdata_simulation(nt, n_schools_treated, np, ns, sigma, ICC, B0, es, seed, attrition_rates, covariate_spec)
        
        
        
      } else if (input$simType == "Multisite Trial") {
        np <- input$np
        ns <- input$ns
        tpi <- as.numeric(unlist(strsplit(as.character(input$tpi), ",")))
        sigma <- input$sigma
        sigmab0 <- input$sigmab0
        sigmab1 <- input$sigmab1
        B0 <- input$B0
        es <- as.numeric(unlist(strsplit(as.character(input$es), ",")))
        seed <- input$seed
        attrition_rates <- as.numeric(unlist(strsplit(as.character(input$attrition_rates), ",")))
        
        mstdata_simulation(nt, tpi, np, ns, sigma, sigmab0, sigmab1, B0, es, seed, attrition_rates, covariate_spec)
        
      } else if (input$simType == "Simple Randomised Trial") {
        np <- input$np_srt
        tpi <- as.numeric(unlist(strsplit(as.character(input$tpi_srt), ",")))
        sigma <- input$sigma_srt
        B0 <- input$B0_srt
        es <- as.numeric(unlist(strsplit(as.character(input$es_srt), ",")))
        seed <- input$seed_srt
        attrition_rates <- as.numeric(unlist(strsplit(as.character(input$attrition_rates_srt), ",")))
        
        srtdata_simulation(nt, tpi, np, sigma, B0, es, seed, attrition_rates, covariate_spec)
      }
    }, error = function(e) {
      showNotification(paste("❌ Simulation Error:", e$message), type = "error", duration = NULL)
      return(NULL)
    })
  })
  
  extract_covariates <- function(prefix, n) {
    if (n == 0 || is.null(n)) return(list())
    lapply(1:n, function(i) {
      ns <- function(id) paste0(prefix, i, "_", id)
      
      type <- input[[ns("type")]]
      name <- input[[ns("name")]]
      if (is.null(type) || is.null(name) || name == "") {
        stop(paste("❌ Covariate", i, "is missing a name. Please provide a name for each covariate."))
      }
      
      if (type == "continuous") {
        list(
          name = name,
          type = "continuous",
          sd = as.numeric(input[[ns("sd")]]),
          effect = as.numeric(input[[ns("effect")]])
        )
        
      } else if (type == "categorical") {
        levels_raw <- input[[ns("levels")]]
        probs_raw <- input[[ns("probs")]]
        ref <- input[[ns("ref")]]
        
        if (any(sapply(list(levels_raw, probs_raw, ref), is.null))) return(NULL)
        
        levels <- trimws(strsplit(as.character(levels_raw), ",")[[1]])
        probs <- as.numeric(strsplit(as.character(probs_raw), ",")[[1]])
        non_ref_levels <- levels[levels != ref]
        
        # Validate that probabilities and levels match
        if (length(levels) < 2) {
          stop(paste("Covariate", name, ": Must have at least two levels."))
        }
        if (!ref %in% levels) {
          stop(paste("Covariate", name, ": Reference category must be one of the levels."))
        }
        if (length(probs) != length(levels)) {
          stop(paste("Covariate", name, ": Number of probabilities must match number of levels."))
        }
        if (abs(sum(probs) - 1) > 1e-6) {
          stop(paste("Covariate", name, ": Probabilities must sum to 1."))
        }
        
        # Extract coefficients dynamically for each non-reference level
        effects <- sapply(seq_along(non_ref_levels), function(j) {
          effect_val <- input[[paste0(prefix, i, "_effect_", j)]]
          if (is.null(effect_val)) {
            stop(paste("Covariate", name, ": Missing effect for level", non_ref_levels[j]))
          }
          effect_val
        })
        
        list(
          name = name,
          type = "categorical",
          levels = levels,
          probs = probs,
          reference = ref,
          effects = setNames(effects, non_ref_levels)
        )
        
        
      } else {
        return(NULL)
      }
    }) |> Filter(Negate(is.null), x = _)
  }
  
  
  
  
  # For CRT
  output$covariate_inputs_crt <- renderUI({
    n <- input$num_covariates_crt
    if (is.null(n) || n == 0) return(NULL)
    req(n)
    lapply(1:n, function(i) {
      ns <- function(id) paste0("crt_cov", i, "_", id)
      tagList(
        textInput(ns("name"), paste("Covariate", i, "Name")),
        selectInput(ns("type"), "Type", choices = c("continuous", "categorical")),
        conditionalPanel(
          condition = paste0("input.", ns("type"), " == 'continuous'"),
          numericInput(ns("sd"), "Standard Deviation", 1),
          numericInput(ns("effect"), "Coefficient", 0)
        ),
        conditionalPanel(
          condition = paste0("input.", ns("type"), " == 'categorical'"),
          textInput(ns("levels"), "Categories (comma-separated)", "A,B,C"),
          textInput(ns("probs"), "Probabilities (comma-separated)", "0.3,0.3,0.4"),
          selectInput(ns("ref"), "Reference Category", choices = NULL),
          uiOutput(ns("dynamic_effect_inputs"))
        ),
        tags$hr()
      )
    })
  })
  
  
  # Multisite Trial covariate inputs
  output$covariate_inputs_mst <- renderUI({
    n <- input$num_covariates_mst
    if (is.null(n) || n == 0) return(NULL)
    req(n)
    lapply(1:n, function(i) {
      ns <- function(id) paste0("mst_cov", i, "_", id)
      tagList(
        textInput(ns("name"), paste("Covariate", i, "Name")),
        selectInput(ns("type"), "Type", choices = c("continuous", "categorical")),
        conditionalPanel(
          condition = paste0("input.", ns("type"), " == 'continuous'"),
          numericInput(ns("sd"), "Standard Deviation", 1),
          numericInput(ns("effect"), "Coefficient", 0)
        ),
        conditionalPanel(
          condition = paste0("input.", ns("type"), " == 'categorical'"),
          textInput(ns("levels"), "Categories (comma-separated)", "A,B,C"),
          textInput(ns("probs"), "Probabilities (comma-separated)", "0.3,0.3,0.4"),
          selectInput(ns("ref"), "Reference Category", choices = NULL),
          uiOutput(ns("dynamic_effect_inputs"))
        ),
        tags$hr()
      )
    })
  })
  
  
  # Simple Randomised Trial covariate inputs
  output$covariate_inputs_srt <- renderUI({
    n <- input$num_covariates_srt
    if (is.null(n) || n == 0) return(NULL)
    req(n)
    lapply(1:n, function(i) {
      ns <- function(id) paste0("srt_cov", i, "_", id)
      tagList(
        textInput(ns("name"), paste("Covariate", i, "Name")),
        selectInput(ns("type"), "Type", choices = c("continuous", "categorical")),
        conditionalPanel(
          condition = paste0("input.", ns("type"), " == 'continuous'"),
          numericInput(ns("sd"), "Standard Deviation", 1),
          numericInput(ns("effect"), "Coefficient", 0)
        ),
        conditionalPanel(
          condition = paste0("input.", ns("type"), " == 'categorical'"),
          textInput(ns("levels"), "Categories (comma-separated)", "A,B,C"),
          textInput(ns("probs"), "Probabilities (comma-separated)", "0.3,0.3,0.4"),
          selectInput(ns("ref"), "Reference Category", choices = NULL),
          uiOutput(ns("dynamic_effect_inputs"))
        ),
        tags$hr()
      )
    })
  })
  
  observe({
    req(!is.null(input$num_covariates_crt), !is.na(input$num_covariates_crt))
    n <- input$num_covariates_crt
    if (n <= 0) return(NULL)
    
    for (i in 1:n) {      
      local({
        idx <- i
        ns <- function(id) paste0("crt_cov", idx, "_", id)
        
        output[[ns("dynamic_effect_inputs")]] <- renderUI({
          levels_raw <- input[[ns("levels")]]
          ref <- input[[ns("ref")]]
          req(levels_raw, ref)
          
          levels <- trimws(unlist(strsplit(levels_raw, ",")))
          non_ref_levels <- levels[levels != ref]
          
          if (length(non_ref_levels) == 0) return(NULL)
          
          lapply(seq_along(non_ref_levels), function(j) {
            numericInput(
              inputId = paste0(ns("effect_"), j),
              label = paste("Coefficient for category:", non_ref_levels[j]),
              value = 0,
              step = 0.01
            )
            
          })
        })
      })
    }
  })
  
  
  # --- Observer to update ref category selections ---
  observe({
    req(!is.null(input$num_covariates_crt), !is.na(input$num_covariates_crt))
    n <- input$num_covariates_crt
    if (n <= 0) return(NULL)
    
    for (i in 1:n) {
      ns <- function(id) paste0("crt_cov", i, "_", id)
      levels_input <- input[[ns("levels")]]
      if (!is.null(levels_input)) {
        levels <- trimws(unlist(strsplit(levels_input, ",")))
        updateSelectInput(session, ns("ref"), choices = levels, selected = levels[1])
      }
    }
  })
  
  
  observe({
    req(!is.null(input$num_covariates_mst), !is.na(input$num_covariates_mst))
    n <- input$num_covariates_mst
    if (n <= 0) return(NULL)
    
    for (i in 1:n) {
      local({
        idx <- i
        ns <- function(id) paste0("mst_cov", idx, "_", id)
        
        output[[ns("dynamic_effect_inputs")]] <- renderUI({
          levels_raw <- input[[ns("levels")]]
          ref <- input[[ns("ref")]]
          req(levels_raw, ref)
          
          levels <- trimws(unlist(strsplit(levels_raw, ",")))
          non_ref_levels <- levels[levels != ref]
          
          if (length(non_ref_levels) == 0) return(NULL)
          
          lapply(seq_along(non_ref_levels), function(j) {
            numericInput(
              inputId = paste0(ns("effect_"), j),
              label = paste("Coefficient for category:", non_ref_levels[j]),
              value = 0,
              step = 0.01
            )
            
          })
        })
      })
    }
  })
  
  
  # --- Observer to update ref category selections ---
  observe({
    req(!is.null(input$num_covariates_mst), !is.na(input$num_covariates_mst))
    n <- input$num_covariates_mst
    if (n <= 0) return(NULL)
    
    for (i in 1:n) {      ns <- function(id) paste0("mst_cov", i, "_", id)
      levels_input <- input[[ns("levels")]]
      if (!is.null(levels_input)) {
        levels <- trimws(unlist(strsplit(levels_input, ",")))
        updateSelectInput(session, ns("ref"), choices = levels, selected = levels[1])
      }
    }
  })
  
  
  observe({
    req(!is.null(input$num_covariates_srt), !is.na(input$num_covariates_srt))
    n <- input$num_covariates_srt
    if (n <= 0) return(NULL)
    
    for (i in 1:n) {
      local({
        idx <- i
        ns <- function(id) paste0("srt_cov", idx, "_", id)
        
        output[[ns("dynamic_effect_inputs")]] <- renderUI({
          levels_raw <- input[[ns("levels")]]
          ref <- input[[ns("ref")]]
          req(levels_raw, ref)
          
          levels <- trimws(unlist(strsplit(levels_raw, ",")))
          non_ref_levels <- levels[levels != ref]
          
          if (length(non_ref_levels) == 0) return(NULL)
          
          lapply(seq_along(non_ref_levels), function(j) {
            numericInput(
              inputId = paste0(ns("effect_"), j),
              label = paste("Coefficient for category:", non_ref_levels[j]),
              value = 0,
              step = 0.01
            )
            
          })
        })
      })
    }
  })
  
  
  # --- Observer to update ref category selections ---
  observe({
    req(!is.null(input$num_covariates_srt), !is.na(input$num_covariates_srt))
    n <- input$num_covariates_srt
    if (n <= 0) return(NULL)
    
    for (i in 1:n) {      ns <- function(id) paste0("srt_cov", i, "_", id)
      levels_input <- input[[ns("levels")]]
      if (!is.null(levels_input)) {
        levels <- trimws(unlist(strsplit(levels_input, ",")))
        updateSelectInput(session, ns("ref"), choices = levels, selected = levels[1])
      }
    }
  })
  
  
  importedData <- reactive({
    req(input$futDataSource)  # Ensure futDataSource input exists
    
    if (input$futDataSource == "Use Simulated Data") {
      if (is.null(simulatedData())) {
        showNotification("Please simulate data first before using it for analysis.", type = "warning")
        return(NULL)
      }
      simulatedData()
    } else {
      req(input$dataset)
      read.csv(input$dataset$datapath)
    }
  })
  
  
  
  observe({
    data <- importedData()
    updateSelectInput(session, "post_var",
                    choices = names(data),
                    selected = if ("posttest" %in% names(data)) "posttest" else NULL)

    updateSelectInput(session, "intervention_var",
                    choices = names(data),
                    selected = if ("interventions" %in% names(data)) "interventions" else NULL)

    updateSelectInput(session, "random_var",
                    choices = names(data),
                    selected = if ("schools" %in% names(data)) "schools" else NULL)
  })
  
  output$post_var_select <- renderUI({
    selectInput("post_var", "Select Post-test Outcome", choices = NULL, multiple = TRUE)
  })
  
  output$intervention_var_select <- renderUI({
    selectInput("intervention_var", "Select Intervention Variables", choices = NULL, multiple = TRUE)
  })
  
  output$random_var_select <- renderUI({
    selectInput("random_var", "Select Clustering Variable", choices = NULL)
  })
  
  output$covariate_cont_select_fut <- renderUI({
    selectInput("covariates_cont_fut", "Select Continuous Covariates", 
                choices = names(importedData()), multiple = TRUE)
  })
  
  output$covariate_cat_select_fut <- renderUI({
    selectInput("covariates_cat_fut", "Select Categorical Covariates", 
                choices = names(importedData()), multiple = TRUE)
  })
  
  
  output$futilityTimingMessage <- renderText({
    futilityTimingText()
  })
  
  futilityTimingText <- reactiveVal("")
  futilityData <- eventReactive(input$analyzeFutility, {
    show_modal_spinner(spin = "circle", text = "Analyzing futility data...")
    
    tryCatch({
      start_time <- Sys.time()
      data <- importedData()  # Load the data
      post_vars <- input$post_var  # Get post-test variable
      intervention_column <- input$intervention_var  # Get intervention variable
      cont_covs <- input$covariates_cont_fut
      cat_covs <- input$covariates_cat_fut
      
      # Convert categorical covariates to factors
      if (!is.null(cat_covs)) {
        for (cov in cat_covs) {
          data[[cov]] <- as.factor(data[[cov]])
        }
      }
      
      # Combine all covariates for the model
      covariates <- c(cont_covs, cat_covs)
      
      # Ensure that only one post-test variable and one intervention variable is selected
      if (length(post_vars) != 1 || length(intervention_column) != 1) {
        stop("Please select exactly one post-test and one intervention variable.")
      }
      
      # Call the appropriate futility function based on the simulation type
      result <- if (input$futSimType == "Cluster-Randomised Trial") {
        crtfutility(
          data = data, 
          post_vars = post_vars,
          intervention_column = intervention_column, 
          Random = input$random_var, 
          Nsim = input$crtNsim, 
          Threshold = input$crtThreshold, 
          ProbThreshold = input$crtProbThreshold, 
          continuous_covariates = cont_covs,
          categorical_covariates = cat_covs
          # Pass covariates directly
        )
      } else if (input$futSimType == "Multisite Trial") {
        mstfutility(
          data = data, 
          post_vars = post_vars,
          intervention_column = intervention_column, 
          Random = input$random_var, 
          Nsim = input$crtNsim, 
          Threshold = input$crtThreshold, 
          ProbThreshold = input$crtProbThreshold, 
          continuous_covariates = cont_covs,
          categorical_covariates = cat_covs
          # Pass covariates directly
        )
      } else {
        srtfutility(
          data = data, 
          post_vars = post_vars,
          intervention_column = intervention_column, 
          Nsim = input$crtNsim, 
          Threshold = input$crtThreshold, 
          ProbThreshold = input$crtProbThreshold, 
          continuous_covariates = cont_covs,
          categorical_covariates = cat_covs
          # Pass covariates directly
        )
      }
      
      end_time <- Sys.time()
      time_taken <- round(difftime(end_time, start_time, units = "secs"), 2)
      futilityTimingText(paste("⏱️ Futility analysis completed in", time_taken, "seconds."))
      
      remove_modal_spinner()
      return(result)
      
    }, error = function(e) {
      remove_modal_spinner()
      showNotification(paste("Error: ", e$message), type = "error", duration = NULL)
      return(NULL)  # Return NULL if error occurs
    })
  })
  
  
  # Render the futility table
  output$futilityTable <- renderTable({
    futility_result <- futilityData()  # Get the data.frame returned by the futility functions
    
    if (is.null(futility_result) || nrow(futility_result) == 0) {
      showNotification("Futility data is empty or invalid.", type = "error", duration = NULL)
      return(NULL)
    }
    
    # Rename columns for better readability
    names(futility_result)[names(futility_result) == "Treatment"] <- "Intervention"
    names(futility_result)[names(futility_result) == "ProbES"] <- "P(Effect size > Threshold)"
    
    # Modify the 'Futility' column to give more descriptive labels
    futility_result$Futility <- sapply(1:nrow(futility_result), function(i) {
      if (futility_result$Futility[i] == 1) {
        paste("Intervention", futility_result$Intervention[i], "is futile.")
      } else {
        paste("Intervention", futility_result$Intervention[i], "is not futile.")
      }
    })
    
    # Reorder columns for better display
    futility_result <- futility_result[, c("Intervention", "P(Effect size > Threshold)", "Futility")]
    
    return(futility_result)
  })
  
  # Add a download button for futility results
  output$downloadFutilityData <- downloadHandler(
    filename = function() {
      paste("futility_analysis_results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      futility_result <- futilityData()
      
      if (is.null(futility_result) || nrow(futility_result) == 0) {
        showNotification("No futility analysis results available for download.", type = "error", duration = NULL)
        return(NULL)
      }
      
      # Rename columns before saving
      names(futility_result)[names(futility_result) == "Treatment"] <- "Intervention"
      names(futility_result)[names(futility_result) == "ProbES"] <- "P(Effect size > Threshold)"
      
      # Save the data as CSV
      write.csv(futility_result, file, row.names = FALSE)
    }
  )
  
  
  
  # Event to generate and preview the futility plot
  futilityPlotData <- reactiveVal(NULL)
  
  observeEvent(input$plotFutility, {
    show_modal_spinner(spin = "circle", text = "Plotting posterior probabilities...")
    
    data <- importedData()
    nt <- length(input$post_var)
    
    futilityPlotData(list(
      nt = nt,
      data = data,
      Random = input$random_var,
      Nsim = if (input$futSimType == "Cluster-Randomised Trial") input$crtNsim else input$mstNsim,
      Threshold = input$crtThreshold,
      ProbThreshold = input$crtProbThreshold,
      covariates = input$covariates_fut,  # Add covariates as input
      threshold_range = input$threshold_range
    ))
  })
  
  # Dynamically generate intervention color + label inputs
  output$interventionCustomUI <- renderUI({
    data <- if (input$plotDataSource == "Use Simulated Data") {
      simulatedData()
    } else {
      req(input$plotDataset)
      tryCatch(read.csv(input$plotDataset$datapath), error = function(e) NULL)
    }
    
    req(data)
    intervention_var <- input$intervention_var_plot
    req(intervention_var)
    
    interventions <- sort(unique(data[[intervention_var]]))
    interventions <- interventions[interventions != 0]
    if (length(interventions) == 0) return(NULL)
    
    # 🎨 Distinct color palette (expand if needed)
    default_colors <- c(
      "#1f77b4",  "#d62728", "#2ca02c", "#ff7f0e",
      "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
      "#bcbd22", "#17becf"
    )
    
    lapply(seq_along(interventions), function(i) {
      tagList(
        textInput(
          inputId = paste0("label_int_", i),
          label = paste("Rename Intervention", interventions[i]),
          value = paste("Intervention", interventions[i])
        ),
        colourInput(
          inputId = paste0("color_int_", i),
          label = paste("Color for Intervention", interventions[i]),
          value = default_colors[(i - 1) %% length(default_colors) + 1]  # Wrap around if needed
        )
      )
    })
  })
  
  
  
  output$futilityPlot <- renderPlot({
    plot_data <- futilityPlotData()
    req(plot_data)
    
    if (input$futSimType == "Cluster-Randomised Trial") {
      # Assuming input$plotSimType == "Cluster-Randomised Trial"
      plot_posterior_probabilities_crt(
        data = plot_data$data,
        post_vars = plot_data$post_var,
        intervention_column = plot_data$intervention_var,
        Random = plot_data$Random,
        Nsim = plot_data$Nsim,
        covariates = plot_data$covariates,
        VerticalLine = plot_data$VerticalLine,
        VerticalLineColor = input$verticalLineColor,
        ProbThreshold = plot_data$ProbThreshold,
        HorizontalLineColor = input$horizontalLineColor,
        threshold_range = plot_data$threshold_range,
        plot_title = input$custom_plot_title,
        x_label = input$custom_x_label,
        y_label = input$custom_y_label,
        x_breaks = seq(0, 1, by = input$x_breaks_custom),
        y_breaks = seq(0, 1, by = input$y_breaks_custom),
        custom_colors = {
          interventions <- sort(unique(plot_data$data[[plot_data$intervention_var]]))
          interventions <- interventions[interventions != 0]
          labels <- sapply(seq_along(interventions), function(i) input[[paste0("label_int_", i)]])
          colors <- sapply(seq_along(interventions), function(i) input[[paste0("color_int_", i)]])
          names(colors) <- labels
          colors
        },
        custom_labels = {
          interventions <- sort(unique(plot_data$data[[plot_data$intervention_var]]))
          interventions <- interventions[interventions != 0]
          sapply(seq_along(interventions), function(i) input[[paste0("label_int_", i)]])
        }
      )
      
    } else if (input$futSimType == "Multisite Trial") {
      plot_posterior_probabilities_mst(
        data = plot_data$data,
        post_vars = input$post_var,
        intervention_column = input$intervention_var,
        Random = plot_data$Random,
        Nsim = plot_data$Nsim,
        covariates = plot_data$covariates,  # Add covariates to function call
        VerticalLine = input$crtVerticalLine,  # Handle vertical line input
        ProbThreshold = input$crtProbThreshold,  # Handle probability threshold input
        threshold_range = plot_data$threshold_range
      )
    } else {
      plot_posterior_probabilities_srt(
        data = plot_data$data,
        post_vars = input$post_var,
        intervention_column = input$intervention_var,
        Nsim = plot_data$Nsim,
        covariates = plot_data$covariates,  # Add covariates to function call
        VerticalLine = input$crtVerticalLine,  # Handle vertical line input
        ProbThreshold = input$crtProbThreshold,  # Handle probability threshold input
        threshold_range = plot_data$threshold_range
      )
    }
    
    remove_modal_spinner()  # Ensure the spinner is removed only after the plot is rendered
  })
  
  output$downloadFutilityPlot <- downloadHandler(
    filename = function() {
      paste("futility_plot-", Sys.Date(), ".tiff", sep = "")
    },
    content = function(file) {
      plot_data <- futilityPlotData()
      req(plot_data)
      
      # Open a TIFF device
      tiff(file, width = 9, height = 7, units = "in", res = 300)
      
      # Call the appropriate plotting function based on the simulation type
      if (input$futSimType == "Cluster-Randomised Trial") {
        # Assuming input$plotSimType == "Cluster-Randomised Trial"
        plot_posterior_probabilities_crt(
          data = plot_data$data,
          post_vars = plot_data$post_var,
          intervention_column = plot_data$intervention_var,
          Random = plot_data$Random,
          Nsim = plot_data$Nsim,
          covariates = plot_data$covariates,
          VerticalLine = plot_data$VerticalLine,
          VerticalLineColor = input$verticalLineColor,
          ProbThreshold = plot_data$ProbThreshold,
          HorizontalLineColor = input$horizontalLineColor,
          threshold_range = plot_data$threshold_range,
          plot_title = input$custom_plot_title,
          x_label = input$custom_x_label,
          y_label = input$custom_y_label,
          x_breaks = seq(0, 1, by = input$x_breaks_custom),
          y_breaks = seq(0, 1, by = input$y_breaks_custom),
          custom_colors = {
            interventions <- sort(unique(plot_data$data[[plot_data$intervention_var]]))
            interventions <- interventions[interventions != 0]
            labels <- sapply(seq_along(interventions), function(i) input[[paste0("label_int_", i)]])
            colors <- sapply(seq_along(interventions), function(i) input[[paste0("color_int_", i)]])
            names(colors) <- labels
            colors
          },
          custom_labels = {
            interventions <- sort(unique(plot_data$data[[plot_data$intervention_var]]))
            interventions <- interventions[interventions != 0]
            sapply(seq_along(interventions), function(i) input[[paste0("label_int_", i)]])
          }
        )
        
      } else if (input$futSimType == "Multisite Trial") {
        plot_posterior_probabilities_mst(
          data = plot_data$data,
          post_vars = input$post_var,
          intervention_column = input$intervention_var,
          Random = plot_data$Random,
          Nsim = plot_data$Nsim,
          covariates = plot_data$covariates,  # Add covariates to function call
          VerticalLine = plot_data$VerticalLine,   # Pass VerticalLine if available
          ProbThreshold = plot_data$ProbThreshold, # Pass ProbThreshold if available
          threshold_range = plot_data$threshold_range
        )
      } else if (input$futSimType == "Simple Randomised Trial") {
        plot_posterior_probabilities_srt(
          data = plot_data$data,
          post_vars = input$post_var,
          intervention_column = input$intervention_var,
          Nsim = plot_data$Nsim,
          covariates = plot_data$covariates,  # Add covariates to function call
          VerticalLine = plot_data$VerticalLine,   # Pass VerticalLine if available
          ProbThreshold = plot_data$ProbThreshold, # Pass ProbThreshold if available
          threshold_range = plot_data$threshold_range
        )
      }
      
      # Close the TIFF device
      dev.off()
    }
  )
  
  importedMultilevelData <- reactive({
    req(input$mlDataSource)
    
    if (input$mlDataSource == "Use Simulated Data") {
      if (is.null(simulatedData())) {
        showNotification("Please simulate data first before using it for analysis.", type = "warning")
        return(NULL)
      }
      simulatedData()
    } else {
      req(input$multilevelDataset)
      read.csv(input$multilevelDataset$datapath)
    }
  })
  
  
  output$custom_intervention_UI_multilevel <- renderUI({
    req(input$intervention_var_multilevel)
    data <- importedMultilevelData()
    interventions <- sort(unique(data[[input$intervention_var_multilevel]]))
    interventions <- interventions[interventions != 0]
    
    if (length(interventions) == 0) return(NULL)
    
    default_colors <- c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd")
    
    lapply(seq_along(interventions), function(i) {
      tagList(
        textInput(paste0("label_ml_int_", i), paste("Label for Intervention", interventions[i]), value = paste("Intervention", interventions[i])),
        colourInput(paste0("color_ml_int_", i), paste("Color for Intervention", interventions[i]), value = default_colors[(i - 1) %% length(default_colors) + 1])
      )
    })
  })
  
  
  observe({
    data <- importedMultilevelData()
    updateSelectInput(session, "post_var_multilevel",
                      choices = names(data),
                      selected = if ("posttest" %in% names(data)) "posttest" else NULL)
    
    updateSelectInput(session, "intervention_var_multilevel",
                      choices = names(data),
                      selected = if ("interventions" %in% names(data)) "interventions" else NULL)
    
    updateSelectInput(session, "random_var_multilevel",
                      choices = names(data),
                      selected = if ("schools" %in% names(data)) "schools" else NULL)
  })
  
  output$post_var_select_multilevel <- renderUI({
    selectInput("post_var_multilevel", "Select Post-test Outcome", choices = NULL, multiple = TRUE)
  })
  
  output$intervention_var_select_multilevel <- renderUI({
    selectInput("intervention_var_multilevel", "Select Intervention Variables", choices = NULL, multiple = TRUE)
  })
  
  
  output$random_var_select_multilevel <- renderUI({
    conditionalPanel(
      condition = "input.method != 'srtBayes' && input.method != 'srtFREQ'",
      selectInput("random_var_multilevel", "Select Clustering Variable", choices = NULL)
    )
  })
  
  output$covariate_cont_select_multilevel <- renderUI({
    selectInput("covariates_cont_multilevel", "Select Continuous Covariates", 
                choices = names(importedMultilevelData()), multiple = TRUE)
  })
  
  output$covariate_cat_select_multilevel <- renderUI({
    selectInput("covariates_cat_multilevel", "Select Categorical Covariates", 
                choices = names(importedMultilevelData()), multiple = TRUE)
  })
  
  
  importedSupData <- reactive({
    req(input$supDataSource)
    
    if (input$supDataSource == "Use Simulated Data") {
      if (is.null(simulatedData())) {
        showNotification("Please simulate data first before using it for analysis.", type = "warning")
        return(NULL)
      }
      simulatedData()
    } else {
      req(input$supDataset)
      read.csv(input$supDataset$datapath)
    }
  })
  
  
  # Update the available columns for post_var, intervention_var, and covariates dynamically
  observe({
    data <- importedSupData()
    updateSelectInput(session, "post_var_sup",
                      choices = names(data),
                      selected = if ("posttest" %in% names(data)) "posttest" else NULL)
    
    updateSelectInput(session, "intervention_var_sup",
                      choices = names(data),
                      selected = if ("interventions" %in% names(data)) "interventions" else NULL)
    
    updateSelectInput(session, "random_var_sup",
                      choices = names(data),
                      selected = if ("schools" %in% names(data)) "schools" else NULL)
    updateSelectInput(session, "covariates_sup", choices = names(data))  # Multiple select for covariates
  })
  
  # Render UI for selecting variables dynamically
  output$post_var_select_sup <- renderUI({
    selectInput("post_var_sup", "Select Post-test Outcome", choices = NULL)
  })
  
  output$intervention_var_select_sup <- renderUI({
    selectInput("intervention_var_sup", "Select Intervention Variable", choices = NULL)
  })
  
  output$random_var_select_sup <- renderUI({
    conditionalPanel(
      condition = "input.supSimType != 'Simple Randomised Trial'",
      selectInput("random_var_sup", "Select Clustering Variable", choices = NULL)
    )
  })
  
  output$reference_intervention_select <- renderUI({
    data <- importedSupData()
    intervention_col <- input$intervention_var_sup
    
    req(data, intervention_col)
    
    interventions <- sort(unique(data[[intervention_col]]))
    interventions <- interventions[interventions != 0]  # Exclude control if 0
    
    selectInput(
      "crtSupReference",
      label = "Reference Intervention",
      choices = interventions,
      selected = interventions[1], selectize = FALSE
    )
  })
  
  
  output$covariate_cont_select_sup <- renderUI({
    selectInput("covariates_cont_sup", "Select Continuous Covariates", choices = names(importedSupData()), multiple = TRUE)
  })
  
  output$covariate_cat_select_sup <- renderUI({
    selectInput("covariates_cat_sup", "Select Categorical Covariates", choices = names(importedSupData()), multiple = TRUE)
  })
  
  
  output$superiorityTimingMessage <- renderText({
    superiorityTimingText()
  })
  
  
  
  superiorityTimingText <- reactiveVal("")
  # Perform Superiority Analysis when the button is clicked
  superiorityData <- eventReactive(input$analyzeSuperiority, {
    show_modal_spinner(spin = "circle", text = "Analyzing superiority data...")
    
    tryCatch({
      start_time <- Sys.time()
      data <- importedSupData()
      post_var <- input$post_var_sup  # Single selected post-test variable
      intervention_var <- input$intervention_var_sup  # Single selected intervention variable
      cont_covs <- input$covariates_cont_sup
      cat_covs <- input$covariates_cat_sup
      
      # Convert categorical covariates to factors
      if (!is.null(cat_covs)) {
        for (cov in cat_covs) {
          data[[cov]] <- as.factor(data[[cov]])
        }
      }
      
      covariates <- c(cont_covs, cat_covs)
      # Multiple selected covariates
      superiority_threshold <- input$crtSupSuperiorThreshold  # New input for superiority threshold
      
      # Validate that both post-test and intervention variables are selected
      if (is.null(post_var) || is.null(intervention_var)) {
        stop("Please select both a post-test variable and an intervention variable.")
      }
      
      # Superiority analysis based on the simulation type
      result <- if (input$supSimType == "Cluster-Randomised Trial") {
        crtSuperiority(
          data = data, 
          post_var = post_var, 
          intervention_column = intervention_var, 
          Random = input$random_var_sup, 
          Nsim = input$crtSupNsim, 
          Threshold = input$crtSupThreshold, 
          reference_intervention = input$crtSupReference, 
          superiority_threshold = superiority_threshold,  # New parameter
          continuous_covariates = cont_covs,
          categorical_covariates = cat_covs
        )
      } else if (input$supSimType == "Multisite Trial") {
        mstSuperiority(
          data = data, 
          post_var = post_var, 
          intervention_column = intervention_var, 
          Random = input$random_var_sup, 
          Nsim = input$crtSupNsim, 
          Threshold = input$crtSupThreshold, 
          reference_intervention = input$crtSupReference, 
          superiority_threshold = superiority_threshold,  # New parameter
          continuous_covariates = cont_covs,
          categorical_covariates = cat_covs
        )
      } else {
        srtSuperiority(
          data = data, 
          post_var = post_var, 
          intervention_column = intervention_var, 
          Nsim = input$crtSupNsim, 
          Threshold = input$crtSupThreshold, 
          reference_intervention = input$crtSupReference, 
          superiority_threshold = superiority_threshold,  # New parameter
          continuous_covariates = cont_covs,
          categorical_covariates = cat_covs
        )
      }
      end_time <- Sys.time()
      time_taken <- round(difftime(end_time, start_time, units = "secs"), 2)
      superiorityTimingText(paste("⏱️ Superiority analysis completed in", time_taken, "seconds."))
      
      remove_modal_spinner()  # Remove spinner when done
      return(result)
      
    }, error = function(e) {
      remove_modal_spinner()
      showNotification(paste("Error: ", e$message), type = "error", duration = NULL)
      return(NULL)  # Return NULL if error occurs
    })
  })
  
  
  # Display the results of Superiority Analysis
  # Render the superiority table (UI display)
  output$superiorityTable <- renderTable({
    superiority_result <- superiorityData()
    
    if (is.null(superiority_result) || nrow(superiority_result) == 0) {
      showNotification("No results available. Please check your input.", type = "error", duration = NULL)
      return(NULL)
    }
    
    # Rename the "Treatment" column to "Intervention"
    names(superiority_result)[names(superiority_result) == "Treatment"] <- "Intervention"
    
    # Rename the "ProbES" column to "P (Effect size > Threshold)"
    names(superiority_result)[names(superiority_result) == "ProbES"] <- "P(Effect size > Threshold)"
    
    # Modify the "Superiority" column:
    superiority_result$Superiority <- ifelse(
      superiority_result$Superiority == "Reference", 
      "Reference",  # Keep "Reference" unchanged
      ifelse(
        superiority_result$Superiority == "Superior",
        "Superior to the Reference Intervention",
        "Not Superior to the Reference Intervention"
      )
    )
    
    return(superiority_result)
  })
  
  # Download handler for CSV file
  output$downloadSuperiorityData <- downloadHandler(
    filename = function() {
      paste("superiority_analysis_results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      superiority_result <- superiorityData()
      
      if (is.null(superiority_result) || nrow(superiority_result) == 0) {
        showNotification("No superiority analysis results available for download.", type = "error", duration = NULL)
        return(NULL)
      }
      
      # Rename columns for CSV
      names(superiority_result)[names(superiority_result) == "Treatment"] <- "Intervention"
      names(superiority_result)[names(superiority_result) == "ProbES"] <- "P(Effect size > Threshold)"
      
      # Convert "Superiority" column to 0 or 1 for CSV (except "Reference" remains as-is)
      superiority_result$Superiority <- ifelse(
        superiority_result$Superiority == "Reference", 
        "Reference",  # Keep "Reference" unchanged
        ifelse(superiority_result$Superiority == "Superior", 1, 0)
      )
      
      # Save as CSV
      write.csv(superiority_result, file, row.names = FALSE)
    }
  )
  
  
  
  
  output$dataTable <- renderTable({
    head(simulatedData(),n=10)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("simulated_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(simulatedData(), file, row.names = FALSE)
    }
  )
  
  
  output$downloadPlotUI <- renderUI({
    req(input$multilevelPlotFormat)
    format_text <- toupper(input$multilevelPlotFormat)
    downloadButton("downloadPlot", label = paste("Download Plot as", format_text))
  })
  
  output$interventionCustomizationUI <- renderUI({
    req(plotPosteriorData())
    data <- plotPosteriorData()
    
    if (is.null(data)) return(NULL)
    
    interventions <- sort(unique(data$data[[data$intervention_var]]))
    
    lapply(seq_along(interventions), function(i) {
      fluidRow(
        column(6,
               textInput(paste0("int_name_", i), 
                         label = paste("Rename", interventions[i]), 
                         value = interventions[i])
        ),
        column(6,
               colourInput(paste0("int_color_", i), 
                           label = paste("Color for", interventions[i]), 
                           value = "#1f77b4")  # default blue
        )
      )
    })
  })
  
  
  multilevelTimingText <- reactiveVal("")  # Stores timing for multilevel analysis
  output$multilevelTimingMessage <- renderText({
    multilevelTimingText()
  })
  
  
  analysisResult <- eventReactive(input$runAnalysis, {
    show_modal_spinner(spin = "circle", text = "Analyzing data...")
    
    tryCatch({
      start_time <- Sys.time()
      data <- importedMultilevelData()
      post_vars <- input$post_var_multilevel
      intervention_column <- input$intervention_var_multilevel
      random_var <- input$random_var_multilevel
      cont_covs <- input$covariates_cont_multilevel
      cat_covs <- input$covariates_cat_multilevel
      
      # Convert categorical covariates to factors
      if (!is.null(cat_covs)) {
        for (cov in cat_covs) {
          data[[cov]] <- as.factor(data[[cov]])
        }
      }
      
      covariates <- c(cont_covs, cat_covs)
      
      # Validate that the number of selected post-test and intervention variables are the same
      if (length(post_vars) != 1 || length(intervention_column) != 1) {
        stop("Please select exactly one post-test and one intervention variable.")
      }
      
      # If no covariates are selected, set covariates to NULL
      if (is.null(covariates) || length(covariates) == 0) {
        covariates <- NULL  # No covariates provided
      }
      
      # Call the run_analysis function with user inputs
      result <- run_analysis(
        data = data,
        post_vars = post_vars,
        intervention_column = intervention_column,
        Random = random_var,
        Nsim = input$multilevelNsim,
        Threshold = input$multilevelThreshold,
        method = input$method,
        crtFREQoption = input$crtFREQoption,
        nPerm = input$nPerm,
        nBoot = input$nBoot,
        bootType = input$bootType,
        continuous_covariates = cont_covs,
        categorical_covariates = cat_covs,
        input = input  # 🔥 pass entire input to allow access to customization inputs
      )
      
      end_time <- Sys.time()
      time_taken <- round(difftime(end_time, start_time, units = "secs"), 2)
      multilevelTimingText(paste("⏱️ Multilevel analysis completed in", time_taken, "seconds."))
      
      remove_modal_spinner()
      return(result)
    }, error = function(e) {
      remove_modal_spinner()
      showNotification(paste("Error: ", e$message), type = "error", duration = NULL)
      return(NULL)  # Return NULL in case of error
    })
  })
  
  
  output$multilevelPlot <- renderPlot({
    analysisResult()
  })
  
  
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("multilevel_analysis_plot-", Sys.Date(), ".", input$multilevelPlotFormat)
    },
    content = function(file) {
      format <- input$multilevelPlotFormat
      
      # Get width and height from user inputs
      plot_width <- input$multilevelPlotWidth
      plot_height <- input$multilevelPlotHeight
      
      if (format == "pdf") {
        pdf(file, width = plot_width, height = plot_height)
      } else if (format == "svg") {
        svg(file, width = plot_width, height = plot_height)
      } else if (format == "eps") {
        postscript(file, width = plot_width, height = plot_height, horizontal = FALSE, onefile = FALSE, paper = "special")
      } else {
        tiff(file, width = plot_width, height = plot_height, units = "in", res = 300)
      }
      
      print(analysisResult())
      dev.off()
    }
  )
  
  
  
  
  plotTriggered <- reactiveVal(FALSE)
  
  # Plot Posterior Probabilities Tab
  importedPlotData <- reactive({
    req(input$plotDataSource)
    
    if (input$plotDataSource == "Use Simulated Data") {
      if (is.null(simulatedData())) {
        showNotification("Please simulate data first before using it for analysis.", type = "warning")
        return(NULL)
      }
      simulatedData()
    } else {
      req(input$plotDataset)
      read.csv(input$plotDataset$datapath)
    }
  })
  
  
  observe({
    data <- importedPlotData()
    updateSelectInput(session, "post_var_plot",
                      choices = names(data),
                      selected = if ("posttest" %in% names(data)) "posttest" else NULL)
    
    updateSelectInput(session, "intervention_var_plot",
                      choices = names(data),
                      selected = if ("interventions" %in% names(data)) "interventions" else NULL)
    
    updateSelectInput(session, "random_var_plot",
                      choices = names(data),
                      selected = if ("schools" %in% names(data)) "schools" else NULL)
  })
  
  output$post_var_select_plot <- renderUI({
    selectInput("post_var_plot", "Select Post-test Outcome", choices = NULL)
  })
  
  output$intervention_var_select_plot <- renderUI({
    selectInput("intervention_var_plot", "Select Intervention Variable", choices = NULL)
  })
  
  output$random_var_select_plot <- renderUI({
    selectInput("random_var_plot", "Select Clustering Variable", choices = NULL)
  })
  
  output$covariate_cont_select_plot <- renderUI({
    selectInput("covariates_cont_plot", "Select Continuous Covariates", 
                choices = names(importedPlotData()), multiple = TRUE)
  })
  
  output$covariate_cat_select_plot <- renderUI({
    selectInput("covariates_cat_plot", "Select Categorical Covariates", 
                choices = names(importedPlotData()), multiple = TRUE)
  })
  
  
  
  # Observe when the plotPosterior button is clicked
  observeEvent(input$plotPosterior, {
    plotTriggered(TRUE)
    show_modal_spinner(spin = "circle", text = "Plotting posterior probabilities...")
  })
  
  output$downloadPosteriorPlotUI <- renderUI({
    req(input$posteriorPlotFormat)
    format_text <- toupper(input$posteriorPlotFormat)
    downloadButton("downloadPosteriorPlot", label = paste("Download Plot as", format_text))
  })
  
  
  
  # Cache the plot data using a reactive expression
  plotPosteriorData <- eventReactive(input$plotPosterior, {
    data <- if (input$plotDataSource == "Use Simulated Data") {
      simulatedData()
    } else {
      req(input$plotDataset)
      read.csv(input$plotDataset$datapath)
    }
    
    
    # Get user selections
    cont_covs <- input$covariates_cont_plot
    cat_covs <- input$covariates_cat_plot
    
    # Ensure categorical covariates are treated as factors
    if (!is.null(cat_covs)) {
      for (cat in cat_covs) {
        data[[cat]] <- as.factor(data[[cat]])
      }
    }
    
    
    
    # Return full data + parameters
    list(
      data = data,
      post_var = input$post_var_plot,
      intervention_var = input$intervention_var_plot,
      Random = if (input$plotSimType != "Simple Randomised Trial") input$random_var_plot else NULL,
      Nsim = input$crtPlotNsim,
      VerticalLine = if (input$addVerticalLine) input$crtVerticalLine else NULL,
      ProbThreshold = if (input$addHorizontalLine) input$crtPlotProbThreshold else NULL,
      continuous_covariates = cont_covs,
      categorical_covariates = cat_covs,
      threshold_range = input$plot_threshold_range
    )
  })
  
  
  
  posteriorPlotTimingText <- reactiveVal("")
  cachedPlot <- eventReactive(input$plotPosterior, {
    show_modal_spinner(spin = "circle", text = "Plotting posterior probabilities...")
    tryCatch({
      start_time <- Sys.time()
      plot_data <- plotPosteriorData()
      req(plot_data)
      
      p <- NULL
      
      if (input$plotSimType == "Cluster-Randomised Trial") {
        p <- plot_posterior_probabilities_crt(
          data = plot_data$data,
          post_vars = plot_data$post_var,
          intervention_column = plot_data$intervention_var,
          Random = plot_data$Random,
          Nsim = plot_data$Nsim,
          continuous_covariates = plot_data$continuous_covariates,
          categorical_covariates = plot_data$categorical_covariates,
          VerticalLine = plot_data$VerticalLine,
          VerticalLineColor = if (input$enablePlotCustomization) input$verticalLineColor else "#0000FF",
          ProbThreshold = plot_data$ProbThreshold,
          HorizontalLineColor = if (input$enablePlotCustomization) input$horizontalLineColor else "#FF0000",
          threshold_range = plot_data$threshold_range,
          plot_title = if (input$enablePlotCustomization) input$custom_plot_title else "Posterior Probabilities Across Thresholds",
          x_label = if (input$enablePlotCustomization) input$custom_x_label else "Threshold",
          y_label = if (input$enablePlotCustomization) input$custom_y_label else "Posterior Probability",
          x_breaks = if (input$enablePlotCustomization) seq(0, 1, by = input$x_breaks_custom) else seq(0, 1, by = 0.1),
          y_breaks = if (input$enablePlotCustomization) seq(0, 1, by = input$y_breaks_custom) else seq(0, 1, by = 0.1),
          custom_colors = if (input$enablePlotCustomization) {
            interventions <- sort(unique(plot_data$data[[plot_data$intervention_var]]))
            interventions <- interventions[interventions != 0]
            labels <- sapply(seq_along(interventions), function(i) input[[paste0("label_int_", i)]])
            colors <- sapply(seq_along(interventions), function(i) input[[paste0("color_int_", i)]])
            names(colors) <- labels
            colors
          } else NULL,
          custom_labels = if (input$enablePlotCustomization) {
            interventions <- sort(unique(plot_data$data[[plot_data$intervention_var]]))
            interventions <- interventions[interventions != 0]
            sapply(seq_along(interventions), function(i) input[[paste0("label_int_", i)]])
          } else NULL
        )
        
      } else if (input$plotSimType == "Multisite Trial") {
        plot_posterior_probabilities_mst(
          data = plot_data$data,
          post_vars = plot_data$post_var,
          intervention_column = plot_data$intervention_var,
          Random = plot_data$Random,
          Nsim = plot_data$Nsim,
          continuous_covariates = plot_data$continuous_covariates,
          categorical_covariates = plot_data$categorical_covariates,
          VerticalLine = plot_data$VerticalLine,
          VerticalLineColor = if (input$enablePlotCustomization) input$verticalLineColor else "#0000FF",
          ProbThreshold = plot_data$ProbThreshold,
          HorizontalLineColor = if (input$enablePlotCustomization) input$horizontalLineColor else "#FF0000",
          threshold_range = plot_data$threshold_range,
          plot_title = if (input$enablePlotCustomization) input$custom_plot_title else "Posterior Probabilities Across Thresholds",
          x_label = if (input$enablePlotCustomization) input$custom_x_label else "Threshold",
          y_label = if (input$enablePlotCustomization) input$custom_y_label else "Posterior Probability",
          x_breaks = if (input$enablePlotCustomization) seq(0, 1, by = input$x_breaks_custom) else seq(0, 1, by = 0.1),
          y_breaks = if (input$enablePlotCustomization) seq(0, 1, by = input$y_breaks_custom) else seq(0, 1, by = 0.1),
          custom_colors = if (input$enablePlotCustomization) {
            interventions <- sort(unique(plot_data$data[[plot_data$intervention_var]]))
            interventions <- interventions[interventions != 0]
            labels <- sapply(seq_along(interventions), function(i) input[[paste0("label_int_", i)]])
            colors <- sapply(seq_along(interventions), function(i) input[[paste0("color_int_", i)]])
            names(colors) <- labels
            colors
          } else NULL,
          custom_labels = if (input$enablePlotCustomization) {
            interventions <- sort(unique(plot_data$data[[plot_data$intervention_var]]))
            interventions <- interventions[interventions != 0]
            sapply(seq_along(interventions), function(i) input[[paste0("label_int_", i)]])
          } else NULL
        )
      } else {
        plot_posterior_probabilities_srt(
          data = plot_data$data,
          post_vars = plot_data$post_var,
          intervention_column = plot_data$intervention_var,
          Nsim = plot_data$Nsim,
          continuous_covariates = plot_data$continuous_covariates,
          categorical_covariates = plot_data$categorical_covariates,
          VerticalLine = plot_data$VerticalLine,
          VerticalLineColor = if (input$enablePlotCustomization) input$verticalLineColor else "#0000FF",
          ProbThreshold = plot_data$ProbThreshold,
          HorizontalLineColor = if (input$enablePlotCustomization) input$horizontalLineColor else "#FF0000",
          threshold_range = plot_data$threshold_range,
          plot_title = if (input$enablePlotCustomization) input$custom_plot_title else "Posterior Probabilities Across Thresholds",
          x_label = if (input$enablePlotCustomization) input$custom_x_label else "Threshold",
          y_label = if (input$enablePlotCustomization) input$custom_y_label else "Posterior Probability",
          x_breaks = if (input$enablePlotCustomization) seq(0, 1, by = input$x_breaks_custom) else seq(0, 1, by = 0.1),
          y_breaks = if (input$enablePlotCustomization) seq(0, 1, by = input$y_breaks_custom) else seq(0, 1, by = 0.1),
          custom_colors = if (input$enablePlotCustomization) {
            interventions <- sort(unique(plot_data$data[[plot_data$intervention_var]]))
            interventions <- interventions[interventions != 0]
            labels <- sapply(seq_along(interventions), function(i) input[[paste0("label_int_", i)]])
            colors <- sapply(seq_along(interventions), function(i) input[[paste0("color_int_", i)]])
            names(colors) <- labels
            colors
          } else NULL,
          custom_labels = if (input$enablePlotCustomization) {
            interventions <- sort(unique(plot_data$data[[plot_data$intervention_var]]))
            interventions <- interventions[interventions != 0]
            sapply(seq_along(interventions), function(i) input[[paste0("label_int_", i)]])
          } else NULL
        )
      }
      end_time <- Sys.time()
      time_taken <- round(difftime(end_time, start_time, units = "secs"), 2)
      posteriorPlotTimingText(paste("⏱️ Posterior probabilities plot generated in", time_taken, "seconds."))
      
      remove_modal_spinner()  # Remove spinner after the plot is generated
      recordPlot()  # Capture the plot  
    }, error = function(e) {
      remove_modal_spinner()
      showNotification(paste("❌ Error generating plot:", e$message), type = "error", duration = NULL)
      return(NULL)
    })
  })
  
  output$posteriorPlotTimingMessage <- renderText({
    posteriorPlotTimingText()
  })
  
  
  # Render the plot using the cached plot
  output$posteriorPlot <- renderPlot({
    plot <- cachedPlot()
    remove_modal_spinner()  # Remove the spinner after the plot is rendered
    plot  # Simply render the cached plot
  })
  
  output$downloadPosteriorPlot <- downloadHandler(
    filename = function() {
      paste0("posterior_plot-", Sys.Date(), ".", input$posteriorPlotFormat)
    },
    content = function(file) {
      format <- input$posteriorPlotFormat
      
      # 🔧 Get width and height from inputs
      plot_width <- input$posteriorPlotWidth
      plot_height <- input$posteriorPlotHeight
      
      if (format == "pdf") {
        pdf(file, width = plot_width, height = plot_height)
      } else if (format == "svg") {
        svg(file, width = plot_width, height = plot_height)
      } else if (format == "eps") {
        postscript(file, width = plot_width, height = plot_height, horizontal = FALSE, onefile = FALSE, paper = "special")
      } else {
        tiff(file, width = plot_width, height = plot_height, units = "in", res = 300)
      }
      
      replayPlot(cachedPlot())
      dev.off()
    }
  )
  
  
  
  
  output$dataTable <- renderTable({
    req(simulatedData())
    head(simulatedData(), n = 10)
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("simulated_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(simulatedData(), file, row.names = FALSE)
    }
  )
  # Load data for the new treatment tab
  # Load data for the new treatment tab
  importedNewTreatmentData <- reactive({
    req(input$newDataSource)
    
    if (input$newDataSource == "Use Simulated Data") {
      if (is.null(simulatedData())) {
        showNotification("Please simulate data first before using it for analysis.", type = "warning")
        return(NULL)
      }
      simulatedData()
    } else {
      req(input$newTreatmentDataset)
      read.csv(input$newTreatmentDataset$datapath)
    }
  })
  
  
  # Update available choices for column selection when a dataset is uploaded
  # Update available choices for column selection when a dataset is uploaded
  observe({
    data <- importedNewTreatmentData()
    
    updateSelectInput(session, "post_var_new",
                      choices = names(data),
                      selected = if ("posttest" %in% names(data)) "posttest" else NULL)
    
    updateSelectInput(session, "intervention_var_new",
                      choices = names(data),
                      selected = if ("interventions" %in% names(data)) "interventions" else NULL)
    
    updateSelectInput(session, "covariates_new", choices = names(data))
    
    # Pupils Variable (ID for SRT, Pupils for CRT/MST)
    updateSelectInput(session, "pupils_var_new", choices = names(data))
    
    # Only update school selection when the simulation type is CRT or MST
    if (input$newSimType %in% c("Cluster-Randomised Trial", "Multisite Trial")) {
      updateSelectInput(session, "schools_var_new",
                        choices = names(data),
                        selected = if ("schools" %in% names(data)) "schools" else NULL)
    }
  })
  
  # Render UI for selecting relevant columns dynamically
  output$post_var_select_new <- renderUI({
    selectInput("post_var_new", "Select Post-test Outcome", choices = NULL)
  })
  
  output$intervention_var_select_new <- renderUI({
    selectInput("intervention_var_new", "Select Intervention Variable", choices = NULL)
  })
  
  # Render Pupils/ID Variable label dynamically
  output$pupils_var_select_new <- renderUI({
    label_text <- if (input$newSimType == "Simple Randomised Trial") "Select ID Variable" else "Select Pupils Variable"
    selectInput("pupils_var_new", label_text, choices = NULL)
  })
  
  output$covariate_cont_select_new <- renderUI({
    selectInput("covariates_cont_new", "Select Continuous Covariates", choices = names(importedNewTreatmentData()), multiple = TRUE)
  })
  
  output$covariate_cat_select_new <- renderUI({
    selectInput("covariates_cat_new", "Select Categorical Covariates", choices = names(importedNewTreatmentData()), multiple = TRUE)
  })
  
  
  # Conditionally render schools variable selection only for CRT and MST
  output$schools_var_select_new <- renderUI({
    if (input$newSimType %in% c("Cluster-Randomised Trial", "Multisite Trial")) {
      selectInput("schools_var_new", "Select Schools Variable", choices = NULL)
    }
  })
  
  # Function to validate if selected columns exist in the data
  validate_columns <- function(data, required_columns) {
    cat("Validating selected columns:\n")
    print(required_columns)
    
    missing_columns <- setdiff(required_columns, names(data))
    if (length(missing_columns) > 0) {
      stop(paste("The following columns are missing in the dataset:", paste(missing_columns, collapse = ", ")))
    }
  }
  
  # Event handler for adding a new treatment
  newTreatmentData <- eventReactive(input$addNewTreatment, {
    show_modal_spinner(spin = "circle", text = "Adding new treatment...")
    
    tryCatch({
      data <- importedNewTreatmentData()
      
      # Required inputs for new treatment tab
      post_col <- input$post_var_new
      intervention_col <- input$intervention_var_new
      schools_col <- input$schools_var_new
      pupils_col <- input$pupils_var_new
      cont_covs <- input$covariates_cont_new
      cat_covs <- input$covariates_cat_new
      
      data <- importedNewTreatmentData()
      
      # Convert selected categorical variables to factor
      if (!is.null(cat_covs)) {
        for (cov in cat_covs) {
          data[[cov]] <- as.factor(data[[cov]])
        }
      }
      
      covariates <- c(cont_covs, cat_covs)
      
      
      # Check if the selected columns are not NULL
      if (is.null(post_col) || is.null(intervention_col) || is.null(pupils_col)) {
        stop("Please make sure all column selections are made.")
      }
      
      # Validate that the user-selected columns exist in the dataset
      required_columns <- c(post_col, intervention_col, pupils_col, covariates)
      validate_columns(data, required_columns)
      
      # Generate new dataset depending on the simulation type
      result <- switch(input$newSimType,
                       "Cluster-Randomised Trial" = add_new_treatmentcrt(
                         existing_data = data,
                         new_schools = input$newSchools,
                         new_pupils_per_school = input$newPupilsPerSchool,
                         es = input$newEffectSize,
                         attrition_rate = input$newAttrition,
                         post_col = post_col,
                         intervention_col = intervention_col,
                         schools_col = schools_col,
                         pupils_col = pupils_col,
                         continuous_covariates = cont_covs,
                         categorical_covariates = cat_covs
                       ),
                       "Multisite Trial" = add_new_treatmentmst(
                         existing_data = data,
                         new_schools = input$newSchools,
                         new_pupils_per_school = input$newPupilsPerSchool,
                         es = input$newEffectSize,
                         attrition_rate = input$newAttrition,
                         treatment_percentage = input$treatmentPercentageMST,
                         post_col = post_col,
                         intervention_col = intervention_col,
                         schools_col = schools_col,
                         pupils_col = pupils_col,
                         continuous_covariates = cont_covs,
                         categorical_covariates = cat_covs
                       ),
                       "Simple Randomised Trial" = add_new_treatmentsrt(
                         existing_data = data,
                         new_pupils = input$newPupilsPerSchool * input$newSchools, # Assuming SRT adds all pupils at once
                         es = input$newEffectSize,
                         attrition_rate = input$newAttrition,
                         post_col = post_col,
                         intervention_col = intervention_col,
                         pupils_col = pupils_col,
                         continuous_covariates = cont_covs,
                         categorical_covariates = cat_covs
                       )
      )
      
      # Check if the resulting data is excessively large
      if (nrow(result) > 1e6) {  # Example threshold of 1 million rows
        stop("The generated dataset is too large. Please adjust input parameters.")
      }
      
      remove_modal_spinner()
      return(result)
      
    }, error = function(e) {
      remove_modal_spinner()
      showNotification(paste("Error: ", e$message), type = "error", duration = NULL)
      return(NULL)
    })
  })
  
  
  
  # Render the modified dataset table
  output$newTreatmentTable <- renderTable({
    req(newTreatmentData())
    head(newTreatmentData(), n = 10)
  })
  
  # Download handler for the modified dataset
  output$downloadNewTreatmentData <- downloadHandler(
    filename = function() {
      paste("new_treatment_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(newTreatmentData(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
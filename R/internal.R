#' Internal plotting function for comparison plots
#'
#' @keywords internal
forestPlotMultiArms <- function(eefAnalyticsList, group, Conditional = TRUE, ES_Total = TRUE, modelNames,
                                intlabels = NULL, intcolors = NULL,
                                maintitle = NULL, xlabel = NULL, ylabel = NULL,
                                vlinecolor = "black"){
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
              intlabels = intlabels,
              intcolors = intcolors,
              maintitle = maintitle,
              xlabel = xlabel,
              ylabel = ylabel,
              vlinecolor = vlinecolor)
  
}




##################################
#      Internal plot function    #
##################################


plotObject <- function(analyticObject, group, Conditional, ES_Total, slope, compare, modelNames,
                       intlabels = NULL, intcolors = NULL,
                       maintitle = NULL, xlabel = NULL, ylabel = NULL,
                       vlinecolor = "black", ...){
  
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
    if (!is.null(intlabels)) {
      old_names <- unique(MyData1$Name)
      if (length(intlabels) == length(old_names)) {
        for (i in seq_along(old_names)) {
          MyData1$Name[MyData1$Name == old_names[i]] <- intlabels[i]
        }
      }
    }
    
    MyData1 <- MyData1[order(MyData1$Variance, decreasing = TRUE), ]
    MyData1$Anot <- paste0(round(MyData1$ES, 2), " [", round(MyData1$LB.95, 2), ", ", round(MyData1$UB.95, 2), "]")
    MyData1$Xaxis <- max(MyData1$UB.95, na.rm = TRUE) + 0.05
    
    # Custom colors
    if (!is.null(intcolors)) {
      color_map <- setNames(intcolors, unique(MyData1$Name))
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
        name = ifelse(is.null(xlabel), expression(paste("Hedges' ", italic("g"))), xlabel)
      ) +
      ylab(ifelse(is.null(ylabel), "Models", ylabel)) +
      geom_vline(xintercept = 0, color = ifelse(is.null(vlinecolor), "black", vlinecolor), linetype = "dashed", alpha = 0.5) +
      theme_bw() +
      theme(
        axis.text.y = element_text(color = "black"),
        text = element_text(size = 16),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
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
    if (!is.null(maintitle)) {
      p <- p + ggtitle(maintitle)
    }
    
    # Apply manual color scale if custom colors are provided
    if (!is.null(color_map)) {
      p <- p + scale_color_manual(values = color_map)
    }
  }
  return(p)
  
}




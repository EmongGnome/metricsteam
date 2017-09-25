#IMR chart using ggplot
#dependencies: ggplot2,IQCC,gridExtra
#nextsteps:
#change default colors for points, display UCL and LC values
#include arguments for setting axis labels
#include checks for violations

IMR<-function(DataSet, DataVector, DataLabel, LowerBoundisZero = NULL, Stage = NULL)
{
  ###Initializations###
  
  plot.new()
  #convert data table DataSet to a data frame
  DataSet<-data.frame(DataSet)
  #check for number of samples
  n<-length(DataVector)
  #make DataLabels an ordered factor
  DataLabel<-factor(DataLabel,levels=DataLabel)
  
  ### I Chart plot###
  
  #compute for c4 factor
  c4<-(sqrt(2/(n-1)))*((factorial(n/2-1))/(factorial((n-1)/2-1)))
  #compute for center line CL
  CL<-mean(DataVector)
  #compute for standard deviation of estimators
  SDofEstimators<-c4*sd(DataVector)
  #compute for UCL and LCL
  UCL=CL+3*SDofEstimators
  LCL=CL-3*SDofEstimators
  
  
  #plot I-chart if there is no staging,i.e, Stage==NULL
  if(is.null(Stage)) 
  {
    
    #check for violations using Nelson's Rules
    #initialize violations index
    ViolationsI<-vector(mode="character",length=n) 
    ViolationsI[1:n]<-"No Violations"
    
    # Nelson's QC rule 1: detect values outside + or -3 sd
    nelsonr1<-which(abs((DataVector - CL)/SDofEstimators) >= 3)
    
    # Nelson's QC rule 2: detect runs of >= 9 points on the same side of the mean
    minrun<-9 
    counts <- sign(DataVector - CL)
    result <- counts
    for (runlength in 2:minrun)
      result <- result + c(counts[runlength:n], rep(0, runlength - 1))
    nelsonr2<-which(abs(result) >= minrun)
    
    # Nelson's QC rule 3: detect strict increase or decrease in >= 6 points in a row
    # Between 6 points you have 5 instances of increasing or decreasing. Therefore minrun - 1.
    minrun<-6
    signs <- sign(c(DataVector[-1], DataVector[n]) - DataVector)
    counts <- signs
    for (rl in 2:(minrun - 1)) {
      counts <- counts + c(signs[rl:n], rep(0, rl - 1))
    }
    nelsonr3<-which(abs(counts) >= minrun - 1)
    
    # Nelson's QC rule 4: 14 points in a row alternating in direction from the mean,
    # or 14 points in a row alternating in increase and decrease
    minrun<-14
    directing_from_mean<-FALSE
    if (directing_from_mean == TRUE) {
      signs <- sign(DataVector - CL)
    } else {
      signs <- sign(c(DataVector[-1],DataVector[n]) - DataVector)
    }
    counts <- signs
    fac <- -1
    for (rl in 2:minrun) {
      counts <- counts + fac * c(signs[rl:n], rep(0, rl - 1))
      fac <- -fac
    }
    counts <- abs(counts)
    nelsonr4<-which(counts >= minrun)
    
    # Nelson's QC rule 5: two out of 3 >2 sd from mean in the same direction
    minrun<-3
    pos <- 1 * ((DataVector - CL)/ SDofEstimators  > 2)
    neg <- 1 * ((DataVector - CL) / SDofEstimators < -2)
    poscounts <- pos
    negcounts <- neg
    for (rl in 2:minrun) {
      poscounts <- poscounts + c(pos[rl:n], rep(0, rl - 1))
      negcounts <- negcounts + c(neg[rl:n], rep(0, rl - 1))
    }
    counts <- apply(cbind(poscounts, negcounts), 1, max)
    nelsonr5<-which(counts >= minrun -1)
    
    # Nelson's QC rule 6: four out of five > 1 sd from mean in the same direction
    minrun<-5
    pos <- 1 * ((DataVector - CL) / SDofEstimators > 1)
    neg <- 1 * ((DataVector - CL) / SDofEstimators < -1)
    poscounts <- pos
    negcounts <- neg
    for (rl in 2:minrun) {
      poscounts <- poscounts + c(pos[rl:n], rep(0, rl - 1))
      negcounts <- negcounts + c(neg[rl:n], rep(0, rl - 1))
    }
    counts <- apply(cbind(poscounts, negcounts), 1, max)
    nelsonr6<-which(counts >= minrun -1)    
    
    # Nelson's QC rule 7: >= 15 points in a row within 1 sd from the mean
    minrun<-15
    within <- 1 * (abs((DataVector - CL) / SDofEstimators) < 1)
    counts <- within
    for (rl in 2:minrun){
      counts <- counts + c(within[rl:n], rep(0, rl - 1))
    }
    nelsonr7<-which(counts >= minrun)
    
    
    # Nelson's QC rule 8: >= 8 points in a row all outside the m + -1s range
    minrun<-8
    outofrange <- 1 * (abs((DataVector - CL) / SDofEstimators) > 1)
    counts <- outofrange
    for (rl in 2:minrun){
      counts <- counts + c(outofrange[rl:n], rep(0, rl - 1))
    }
    nelsonr8<-which(counts >= minrun)
    
    
    
    
    #Assign tags to violations
    ViolationsI[nelsonr1] <- "Outside ctrl limits"
    ViolationsI[nelsonr2] <- "9 pts on same side of CL"
    ViolationsI[nelsonr3] <- "6 pts in a row dec or inc"
    ViolationsI[nelsonr4] <- "14 pts in a row alternating inc then dec"
    ViolationsI[nelsonr5] <- "2 out of 3 pts >2SD from CL"
    ViolationsI[nelsonr6] <- "4 out of 5pts >1SD from CL in same dir"
    ViolationsI[nelsonr7] <- "15 pts in a row within 1SD from CL"
    ViolationsI[nelsonr8] <- ">=8pts outside 1SD from CL"
    
    #set ChangeDefaultColor
    ChangeDefaultColorI<-c("No Violations" = "blue",
                           "Outside control limits"= "red",
                           "9 pts on same side of CL" = "yellow",
                           "6 pts in a row dec or inc"="orange",
                           "14 pts in a row alternating inc then dec" = "magenta4",
                           "2 out of 3 pts >2SD from CL" = "indianred2",
                           "4 out of 5pts >1SD from CL in same dir" = "darkorchid",
                           "15 pts in a row within 1SD from CL" = "gold4",
                           ">=8pts outside 1SD from CL" = "darkmagenta")
    
    
    #compute for Y-axis limits  
    YAxisLimits<-c(LCL,UCL,max(DataVector))
    YLowerBound<-min(YAxisLimits)
    YUpperBound<-max(YAxisLimits)  
    
    #plot points and lines for I control chart
    IChart<-ggplot(DataSet, aes(x=DataLabel,y=DataVector))+
      geom_point(aes(color=ViolationsI))+
      geom_line(aes(group=1))+
      
      #plot center line, UCL, and LCL
      geom_abline(intercept=CL,slope=0,color="black")+
      geom_abline(intercept=UCL,slope=0,color="blue")+
      geom_abline(intercept=LCL,slope=0,color="blue")+
      
      
      #theme settings - vertical x-axis labels  
      theme(axis.text.x=element_text(angle=45,hjust=1))+
      
      coord_cartesian(ylim=c(YLowerBound:YUpperBound))+
      
      #change default colors
      scale_colour_manual(values=ChangeDefaultColorI)
    
    
    
  }#closeif
  
  #plot I-chart if there is staging 
  else
  {
    
    #get length of stage 2,since length of Stage 1 is Stage
    nStage2<-length(DataLabel[Stage:n])
    
    #create column StageEvent to be used as ggplot facet
    DataSet$Event<-"Before"
    DataSet$Event[(Stage+1):n]<-"After"
    DataSet$Event<-factor(DataSet$Event,levels=c("Before","After"))
    
    #compute c4stage1 factor
    c4Stage1<-(sqrt(2/(Stage-1)))*((factorial(Stage/2-1))/(factorial((Stage-1)/2-1)))
    #compute CLStage1
    CLStage1<-mean(DataVector[1:Stage])
    #compute SDofEstimatorsStage1
    SDofEstimatorsStage1<-c4Stage1*sd(DataVector[1:Stage])
    #compute UCLStage2 and LCLStage2
    UCLStage1<-CLStage1+3*SDofEstimatorsStage1
    LCLStage1<-CLStage1-3*SDofEstimatorsStage1
    
    
    #compute c4Stage2 factor
    m<-n-Stage
    c4Stage2<-(sqrt(2/(m-1)))*((factorial(m/2-1))/(factorial((m-1)/2-1)))
    #compute CLStage2
    CLStage2<-mean(DataVector[(Stage+1):n])
    #compute SDofEstimatorsStage2
    SDofEstimatorsStage2<-c4Stage2*sd(DataVector[(Stage+1):n])
    #compute UCLStage2 and LCLStage2
    UCLStage2<-CLStage2+3*SDofEstimatorsStage2
    LCLStage2<-CLStage2-3*SDofEstimatorsStage2
    
    #split DataSet into two data segments DataSegment1 and DataSegment2
    DataSegment1<-DataSet[1:Stage,]
    DataSegment2<-DataSet[(Stage+1):n,]
    
    #mark outliers for Staged I chart#
    i<-1
    OutlierI<-NULL
    while(i<=Stage)
    {
      if(DataVector[i]<=UCLStage1&&DataVector[i]>=LCLStage1)
      {
        OutlierI[i]<-"within"
      } #closeif
      else{OutlierI[i]<-"outlier"}
      i<-i+1
    }#closewhile
    while(i<=n)
    {
      if(DataVector[i]<=UCLStage2&&DataVector[i]>=LCLStage2)
      {
        OutlierI[i]<-"within"
      }
      else{OutlierI[i]<-"outlier"}
      i<-i+1
    }
    
    #plot points and lines for control chart based on Stage
    IChartStage<-ggplot(DataSet, aes(x=DataLabel,y=DataVector))+
      geom_point(aes(color=OutlierI))+
      geom_line(aes(group=1))+
      facet_grid(.~Event,scales="free")+
      
      #plot center line, UCL, and LCL for Stage 1, create new aes layer using data seg
      geom_segment(data=DataSegment1,aes(x=1,y=CLStage1,xend=Stage,yend=CLStage1),
                   color="black",inherit.aes=FALSE)+
      geom_segment(data=DataSegment1,aes(x=1,y=UCLStage1,xend=Stage,yend=UCLStage1),
                   color="blue",inherit.aes=FALSE)+
      geom_segment(data=DataSegment1,aes(x=1,y=LCLStage1,xend=Stage,yend=LCLStage1),
                   color="blue",inherit.aes=FALSE)+
      
      #plot CL, UCL, and LCL for Stage 2, create new aes layer using data seg
      geom_segment(data=DataSegment2,aes(x=1,y=CLStage2,xend=(n-Stage),yend=CLStage2),
                   color="black",inherit.aes=FALSE)+
      geom_segment(data=DataSegment2,aes(x=1,y=UCLStage2,xend=(n-Stage),yend=UCLStage2),
                   color="blue",inherit.aes=FALSE)+
      geom_segment(data=DataSegment2,aes(x=1,y=LCLStage2,xend=(n-Stage),yend=LCLStage2),
                   color="blue",inherit.aes=FALSE)+
      
      
      #theme settings - vertical x-axis labels  
      theme(axis.text.x=element_text(angle=45,hjust=1))
    
    #mark outliers for both Stages
  }#closeelse
  
  ### MR Chart plot###
  
  #compute for d2 and d3 factor for MR
  D2<-d2(n)
  D3<-d3(n)
  #create vector of moving ranges
  MR1Vector<-DataVector
  MR2Vector<-DataVector[2:n]
  MRDataVector[1]<-NA
  MRDataVector[2:n]<-abs(MR1Vector[1:(n-1)]-MR2Vector)
  
  
  #compute for center line CL of MR
  MRCL<-mean(MRDataVector[2:n]) 
  #compute true standard deviation
  MRSD<-MRCL/D2
  #compute for UCL and LCL of MR
  MRUCL=MRCL+3*MRSD
  MRLCL=MRCL-3*MRSD
  
  
  #plot MR chart if there is no staging,i.e, Stage==NULL
  if(is.null(Stage)) 
  {
    
    #initialize violations index
    ViolationsMR<-vector(mode="character",length=n) 
    ViolationsMR[1:n]<-"No Violations"
    MRDataVectorV<-MRDataVector[2:n]
    v<-length(MRDataVectorV)
    
    # Nelson's QC rule 1: detect values outside + or -3 sd
    nelsonr1<-which(abs((MRDataVectorV - MRCL)/MRSD) >= 3)
    
    # Nelson's QC rule 2: detect runs of >= 9 points on the same side of the mean
    minrun<-9 
    counts <- sign(MRDataVectorV - MRCL)
    result <- counts
    for (runlength in 2:minrun)
      result <- result + c(counts[runlength:v], rep(0, runlength - 1))
    nelsonr2<-which(abs(result) >= minrun)
    
    # Nelson's QC rule 3: detect strict increase or decrease in >= 6 points in a row
    # Between 6 points you have 5 instances of increasing or decreasing. Therefore minrun - 1.
    minrun<-6
    signs <- sign(c(MRDataVectorV[-1], MRDataVectorV[v]) - MRDataVectorV)
    counts <- signs
    for (rl in 2:(minrun - 1)) {
      counts <- counts + c(signs[rl:v], rep(0, rl - 1))
    }
    nelsonr3<-which(abs(counts) >= minrun - 1)
    
    # Nelson's QC rule 4: 14 points in a row alternating in direction from the mean,
    # or 14 points in a row alternating in increase and decrease
    minrun<-14
    directing_from_mean<-FALSE
    if (directing_from_mean == TRUE) {
      signs <- sign(MRDataVectorV - MRCL)
    } else {
      signs <- sign(c(MRDataVectorV[-1],MRDataVectorV[v]) - MRDataVectorV)
    }
    counts <- signs
    fac <- -1
    for (rl in 2:minrun) {
      counts <- counts + fac * c(signs[rl:v], rep(0, rl - 1))
      fac <- -fac
    }
    counts <- abs(counts)
    nelsonr4<-which(counts >= minrun)
    
    #adjust indexing
    nelsonr1<-nelsonr1+1
    nelsonr2<-nelsonr2+1
    nelsonr3<-nelsonr3+1
    nelsonr4<-nelsonr4+1
    
    
    #Assign tags to violations
    ViolationsMR[nelsonr1] <- "Outside control limits"
    ViolationsMR[nelsonr2] <- "9 pts on same side of CL"
    ViolationsMR[nelsonr3] <- "6 pts in a row dec or inc"
    ViolationsMR[nelsonr4] <- "14 pts in a row alternating inc then dec"
    
    #set ChangeDefaultColor
    ChangeDefaultColorMR<-c("No Violations" = "blue",
                            "Outside control limits"= "red",
                            "9 pts on same side of CL" = "yellow",
                            "6 pts in a row dec or inc"="orange",
                            "14 pts in a row alternating inc then dec" = "magenta4")
    
    
    #compute for Y-axis limits of MR Chart  
    MRYAxisLimits<-c(MRLCL,MRUCL,max(MRDataVector[2:n]))
    MRYLowerBound<-min(MRYAxisLimits)
    MRYUpperBound<-max(MRYAxisLimits)
    
    #create new data frame for MR chart aesthetic layer
    MRDataSet<-cbind(DataSet,MRDataVector)
    #plot points and lines for MR control chart
    MRChart<-ggplot(MRDataSet, aes(x=DataLabel,y=MRDataVector))+
      geom_point(aes(color=ViolationsMR))+
      geom_line(aes(group=1))+
      
      #plot center line, UCL, and LCL for MR Chart
      geom_abline(intercept=MRCL,slope=0,color="black")+
      geom_abline(intercept=MRUCL,slope=0,color="blue")+
      geom_abline(intercept=MRLCL,slope=0,color="blue")+
      
      
      #theme settings - vertical x-axis labels  
      theme(axis.text.x=element_text(angle=45,hjust=1))+
      
      coord_cartesian(ylim=c(MRYLowerBound:MRYUpperBound))+
      
      #change default colors for violaton points in MR
      scale_colour_manual(values=ChangeDefaultColorMR)
    
    
  }#closeif
  
  #plot MR Chart with staging
  else
  {
    #create vector of moving ranges
    MR1Vector<-DataVector
    MR2Vector<-DataVector[2:n]
    MRDataVector[1]<-NA
    MRDataVector[2:n]<-abs(MR1Vector[1:(n-1)]-MR2Vector)
    
    #create new data frame for MR chart aesthetic layer
    MRDataSet<-cbind(DataSet,MRDataVector)
    
    #get length of stage 2,since length of Stage 1 is Stage
    nMRStage2<-length(DataLabel[Stage:n])
    
    #create column StageEvent to be used as ggplot facet
    MRDataSet$Event<-"Before"
    MRDataSet$Event[(Stage+1):n]<-"After"
    MRDataSet$Event<-factor(MRDataSet$Event,levels=c("Before","After"))
    
    #compute MR stage 1 d2 and d3 values
    D2MRStage1<-d2(Stage)
    D3MRStage1<-d3 (Stage)
    #compute MR stage 1 center line
    MRCLStage1<-mean(MRDataVector[2:n])
    #compute MR stage 1 true standard deviation
    MRSDStage1<-MRCLStage1/D2MRStage1
    #compute MR UCL and LCL for Stage 1
    MRUCLStage1<-MRCLStage1+3*MRSDStage1
    MRLCLStage1<-MRCLStage1-3*MRSDStage1
    
    #compute MR stage 2 d2 and d3 values
    D2MRStage2<-d2(n-Stage)
    D3MRStage2<-d3 (n-Stage)
    #compute MR stage 1 center line
    MRCLStage2<-mean(MRDataVector[(Stage+1):n])
    #compute MR stage 1 true standard deviation
    MRSDStage2<-MRCLStage2/D2MRStage2
    #compute MR UCL and LCL for Stage 1
    MRUCLStage2<-MRCLStage2+3*MRSDStage2
    MRLCLStage2<-MRCLStage2-3*MRSDStage2
    
    
    #split DataSet into two data segments DataSegment1 and DataSegment2
    MRDataSegment1<-MRDataSet[1:Stage,]
    MRDataSegment2<-MRDataSet[(Stage+1):n,]
    
    #mark outliers for staged MR Chart
    i<-2
    OutlierMR<-NULL
    while(i<=Stage)
    {
      if(MRDataVector[i]<=MRUCLStage1&&MRDataVector[i]>=MRLCLStage1)
      {
        OutlierMR[i]<-"within"
      } #closeif
      else{OutlierMR[i]<-"outlier"}
      i<-i+1
    }#closewhile
    while(i<=n)
    {
      if(MRDataVector[i]<=MRUCLStage2&&DataVector[i]>=MRLCLStage2)
      {
        OutlierMR[i]<-"within"
      }
      else{OutlierMR[i]<-"outlier"}
      i<-i+1
    }
    
    #plot points and lines for control chart based on Stage
    MRChartStage<-ggplot(MRDataSet, aes(x=DataLabel,y=MRDataVector))+
      geom_point(aes(color=OutlierMR))+
      geom_line(aes(group=1))+
      facet_grid(.~Event,scales="free")+
      
      #plot center line, UCL, and LCL for Stage 1, create new aes layer using data seg
      geom_segment(data=MRDataSegment1,aes(x=1,y=MRCLStage1,xend=Stage,yend=MRCLStage1),
                   color="black",inherit.aes=FALSE)+
      geom_segment(data=MRDataSegment1,aes(x=1,y=MRUCLStage1,xend=Stage,yend=MRUCLStage1),
                   color="blue",inherit.aes=FALSE)+
      geom_segment(data=MRDataSegment1,aes(x=1,y=MRLCLStage1,xend=Stage,yend=MRLCLStage1),
                   color="blue",inherit.aes=FALSE)+
      
      #plot CL, UCL, and LCL for Stage 2, create new aes layer using data seg
      geom_segment(data=MRDataSegment2,aes(x=1,y=MRCLStage2,xend=(n-Stage),yend=MRCLStage2),
                   color="black",inherit.aes=FALSE)+
      geom_segment(data=MRDataSegment2,aes(x=1,y=MRUCLStage2,xend=(n-Stage),yend=MRUCLStage2),
                   color="blue",inherit.aes=FALSE)+
      geom_segment(data=MRDataSegment2,aes(x=1,y=MRLCLStage2,xend=(n-Stage),yend=MRLCLStage2),
                   color="blue",inherit.aes=FALSE)+
      
      
      #theme settings - vertical x-axis labels  
      theme(axis.text.x=element_text(angle=45,hjust=1))
    
  }
  
  
  #plot MRChart in grid if no staging,i.e Stage==NULL
  if(is.null(Stage))
  {
    grid.arrange(IChart, MRChart, nrow=2)
  }#closeif
  
  #plot MRChart in grid with staging
  else
  {
    grid.arrange(IChartStage, MRChartStage, nrow=2)
  }
  
  print(n)
  
}#closemain
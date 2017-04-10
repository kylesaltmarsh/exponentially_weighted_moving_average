#-----------------------------------------------------------------------------------------#
# Kyle Saltmarsh
#
# EWMA Rule functions
#
#-----------------------------------------------------------------------------------------#



#-----------------------------------------------------------------------------------------#
# Main


# Inputs are default values if not specified
Apply_ADS_EWMA_Rules <- function(Input_Data, 
                                 EWMA_Window_Size = 60,
                                 EWMA_Lambda = 0.3, 
                                 Control_Limit_N_Sigma = 1,
                                 Control_Target = 0, # "EWMA_Mean"
                                 Output_Tags = NULL) {
  
  # Check for last available output from Output_Tags
  # _1 EWMA mean
  # _2 EWMA variance
  # _3 window counter 
  
  # If given an Output_Tags table containing last known value
  if (!is.null(Output_Tags)) {
    EWMA_HWM <- as.numeric(Output_Tags[grep('\\_1$', Output_Tags[["OUTPUT"]])][["VALUE"]])
    EWMA_Variance_HWM <- as.numeric(Output_Tags[grep('\\_2$', Output_Tags[["OUTPUT"]])][["VALUE"]])
    EWMA_Iteration_HWM <- as.numeric(Output_Tags[grep('\\_3$', Output_Tags[["OUTPUT"]])][["VALUE"]])
    if (length(EWMA_HWM) == 0) {
      EWMA_HWM <- mean(Input_Data[[2]],na.rm = TRUE)
      EWMA_Variance_HWM <- 0
      EWMA_Iteration_HWM <- 0
    }
  # Else set HWM to the mean, Variance to 0 and the iteration to 0
  } else {
    EWMA_HWM <- mean(Input_Data[[2]],na.rm = TRUE)
    EWMA_Variance_HWM <- 0
    EWMA_Iteration_HWM <- 0
  }
  
  # Remove rows from end while dim[1] % EWMA_Window_Size != 0
  Trim_Number <- dim(Input_Data)[1]%%EWMA_Window_Size
  Input_Data <- Input_Data[1:(dim(Input_Data)[1]-Trim_Number),]
  
  # Get the input tag name
  target = Input_Data[,colnames(.SD),.SDcols=-"TAGDATE"]
  
  # Store TAGDATE
  Output_Data <- data.table(TAGDATE = Input_Data$TAGDATE)
  
  # Function list to apply EWMA functions
  function_list <- c(Window_Mean, EWMA, EWMA_Variance, EWMA_Lower_Control_Limit, EWMA_Upper_Control_Limit,EWMA_Counter)
  
  # Window mean
  Output_Data[,paste0('Output_',as.integer(1)) := function_list[[as.integer(1)]](Input_Data[[2]], 
                                                                                 EWMA_Window_Size = EWMA_Window_Size)]
  # EWMA
  Output_Data[,paste0('Output_',as.integer(2)) := function_list[[as.integer(2)]](Input_Data[[2]], 
                                                                                 EWMA_Window_Size = EWMA_Window_Size, 
                                                                                 EWMA_Lambda = EWMA_Lambda, 
                                                                                 EWMA_HWM = EWMA_HWM)]
  # Variance
  Output_Data[,paste0('Output_',as.integer(3)) := function_list[[as.integer(3)]](Input_Data[[2]], 
                                                                                 EWMA_Window_Size = EWMA_Window_Size, 
                                                                                 EWMA_Lambda = EWMA_Lambda, 
                                                                                 EWMA_HWM = EWMA_HWM,
                                                                                 EWMA_Iteration_HWM = EWMA_Iteration_HWM,
                                                                                 EWMA_Variance_HWM = EWMA_Variance_HWM)]
  # LCL
  Output_Data[,paste0('Output_',as.integer(4)) := function_list[[as.integer(4)]](Input_Data[[2]], 
                                                                                 EWMA_Window_Size = EWMA_Window_Size, 
                                                                                 EWMA_Lambda = EWMA_Lambda, 
                                                                                 EWMA_Variance_HWM = EWMA_Variance_HWM,
                                                                                 Control_Limit_N_Sigma = Control_Limit_N_Sigma, 
                                                                                 EWMA_HWM = EWMA_HWM,
                                                                                 EWMA_Iteration_HWM = EWMA_Iteration_HWM,
                                                                                 Control_Target = Control_Target)]
  # UCL
  Output_Data[,paste0('Output_',as.integer(5)) := function_list[[as.integer(5)]](Input_Data[[2]], 
                                                                                 EWMA_Window_Size = EWMA_Window_Size, 
                                                                                 EWMA_Lambda = EWMA_Lambda, 
                                                                                 EWMA_Variance_HWM = EWMA_Variance_HWM,
                                                                                 Control_Limit_N_Sigma = Control_Limit_N_Sigma, 
                                                                                 EWMA_HWM = EWMA_HWM,
                                                                                 EWMA_Iteration_HWM = EWMA_Iteration_HWM,
                                                                                 Control_Target = Control_Target)]
  # Iteration counter
  Output_Data[,paste0('Output_',as.integer(6)) := function_list[[as.integer(6)]](Input_Data[[2]], 
                                                                                 EWMA_Window_Size = EWMA_Window_Size, 
                                                                                 EWMA_Iteration_HWM = EWMA_Iteration_HWM)]

  # Return the timeseries of results but not the input
  # Output is repeated by Window_Size to match input data length
  return(Output_Data)
  
}

# Create our EWMA model
Model_Payload <- list( # Package functions and objects as a named list
  Model_Object_Function = Apply_EWMA_Rules, # Mandatory entry point function
  Model_Family_Functions = c("Window_Mean.Rds", 
                             "EWMA.Rds",
                             "EWMA_Variance.Rds",
                             "Lower_Control_Limit.Rds",
                             "Upper_Control_Limit.Rds")#,
  #Model_Object = list() # Enter a training object here if relevant to model
)

saveRDS(Model_Payload,paste0("EWMA.Rds"))# Save as .Rds file

#-----------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------------#
# Rule functions
#
#
#-----------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------------#
# Window data

Window_Mean <- function(Input_Data, 
                        EWMA_Window_Size = 60) {
  
  # Calculate centreline for windows
  Number_windows <- ceiling(length(Input_Data)/EWMA_Window_Size)
  Window_Data_Centreline <- matrix(0,Number_windows,1)
  
  # Aggregate
  # Input data should be a multiple of window size
  for (i in 1:Number_windows) {
    
    if (i == Number_windows) {
      
      Window_Data_Centreline[i] <- mean(Input_Data[(1+(i-1)*EWMA_Window_Size):length(Input_Data)],na.rm = TRUE)
      
    } else {
      
      Window_Data_Centreline[i] <- mean(Input_Data[(1+(i-1)*EWMA_Window_Size):(i*EWMA_Window_Size)],na.rm = TRUE)
      
    }
    
  }
  
  #xbar0 = pop mean or target (0)
  # Input data should be a multiple of window size
  Window_Data_Centreline <- rep(Window_Data_Centreline, each = EWMA_Window_Size)
  
  return(Window_Data_Centreline)
  
}

Rule_Window_Mean = list(Window_Mean = Window_Mean)
saveRDS(Rule_Window_Mean,"Window_Mean.Rds")

#-----------------------------------------------------------------------------------------#
# EWMA
# The Exponentially Weighted Moving Average (EWMA) is a statistic for monitoring the 
# process that averages the data in a way that gives less and less weight to data as 
# they are further removed in time. 

EWMA <- function(Input_Data, 
                 EWMA_Window_Size = 60, 
                 EWMA_Lambda = 0.9, 
                 EWMA_HWM = 0) {
  
  # Calculate centreline for windows
  Number_windows <- ceiling(length(Input_Data)/EWMA_Window_Size)
  Window_Data_Centreline <- matrix(0,Number_windows,1)
  
  # Aggregate
  # Input data should be a multiple of window size
  for (i in 1:Number_windows) {
    
    if (i == Number_windows) {
      
      Window_Data_Centreline[i] <- mean(Input_Data[(1+(i-1)*EWMA_Window_Size):length(Input_Data)],na.rm = TRUE)
      
    } else {
      
      Window_Data_Centreline[i] <- mean(Input_Data[(1+(i-1)*EWMA_Window_Size):(i*EWMA_Window_Size)],na.rm = TRUE)
      
    }
    
  }
  # Initialise EWMA matrix
  EWMA <- matrix(0,Number_windows,1)

  EWMA[1] <- EWMA_HWM #EWMA0 is 0 if not provided or last known value
  
  # Calculate EWMA
  for (i in 2:Number_windows) {
    
    if (is.na(Window_Data_Centreline[[i-1]])) {
      
      EWMA[i] <- EWMA[i-1]
      
    } else {
      
      EWMA[i] <- EWMA_Lambda*Window_Data_Centreline[[i-1]] + (1-EWMA_Lambda)*EWMA[i-1]
      
    }
    
  }
  
  EWMA <- rep(EWMA,each = EWMA_Window_Size)
  
  return(EWMA)
  
}

Rule_EWMA = list(EWMA = EWMA)
saveRDS(Rule_EWMA,"EWMA.Rds")

#-----------------------------------------------------------------------------------------#
# EWMA
# The Exponentially Weighted Moving Average (EWMA) is a statistic for monitoring the 
# process that averages the data in a way that gives less and less weight to data as 
# they are further removed in time. 

EWMA_Variance <- function(Input_Data, 
                          EWMA_Window_Size = 60, 
                          EWMA_Lambda = 0.9, 
                          EWMA_HWM = 0, 
                          EWMA_Variance_HWM = 0, 
                          EWMA_Iteration_HWM = 0) {
  
  # Calculate centreline for windows
  Number_windows <- ceiling(length(Input_Data)/EWMA_Window_Size)
  Window_Data_Centreline <- matrix(0,Number_windows,1)
  
  # Aggregate
  # Input data should be a multiple of window size
  for (i in 1:Number_windows) {
    
    if (i == Number_windows) {
      
      Window_Data_Centreline[i] <- mean(Input_Data[(1+(i-1)*EWMA_Window_Size):length(Input_Data)],na.rm = TRUE)
      
    } else {
      
      Window_Data_Centreline[i] <- mean(Input_Data[(1+(i-1)*EWMA_Window_Size):(i*EWMA_Window_Size)],na.rm = TRUE)
      
    }
    
  }
  # Initialise EWMA matrix
  EWMA <- matrix(0,Number_windows,1)
  EWMA[1] <- EWMA_HWM #EWMA0 is 0 if not provided or last known value
  
  # Calculate EWMA
  for (i in 2:Number_windows) {
    
    if (is.na(Window_Data_Centreline[[i-1]])) {
      
      EWMA[i] <- EWMA[i-1]
      
    } else {
      
      EWMA[i] <- EWMA_Lambda*Window_Data_Centreline[[i-1]] + (1-EWMA_Lambda)*EWMA[i-1]
      
    }
    
  }
  
  Variance <- matrix(0,Number_windows,1)
  
  if (EWMA_Variance_HWM == 0 & EWMA_Iteration_HWM == 0){
    # V1 = (xbar1 - z0)^2 = (xbar1 - xbar0)^2 (xbar0 = z0 = pop mean or target)
    Variance_Temp <- (EWMA_Lambda^2)*(Window_Data_Centreline[2] - EWMA[1])^2
    Variance[1] <- Variance_Temp
  } else if (EWMA_Variance_HWM != 0 & EWMA_Iteration_HWM != 0) {
    Variance_Temp <- EWMA_Variance_HWM
    Variance[1] <- Variance_Temp
  }

  
  for (i in 2:(Number_windows)) {
    
    # Missing data, new variance is same as old variance
    if (is.na(Window_Data_Centreline[[i]])) {
      
      Variance[i] <- Variance[i-1]
      
    } else {
      
      # Variance[i] <- ((1/(i+EWMA_Iteration_HWM))*((i+EWMA_Iteration_HWM)-1)*Variance[i-1]) + (1/(i+EWMA_Iteration_HWM))*((Window_Data_Centreline[i+1] - EWMA[i])^2)
      Variance[i] <- ((i+EWMA_Iteration_HWM-2)/(i+EWMA_Iteration_HWM-1))*Variance[i-1] + ((EWMA_Lambda^2)/(i+EWMA_Iteration_HWM-1))*((Window_Data_Centreline[i] - EWMA[i-1])^2)
    }

  }
  
  Variance <- rep(Variance, each  = EWMA_Window_Size)
  
  return(Variance)
  
}

Rule_EWMA_Variance = list(EWMA_Variance = EWMA_Variance)
saveRDS(Rule_EWMA_Variance,"EWMA_Variance.Rds")

#-----------------------------------------------------------------------------------------#
# LCL


EWMA_Lower_Control_Limit <- function(Input_Data, 
                                     EWMA_Window_Size = 60, 
                                     EWMA_Lambda = 0.9, 
                                     Control_Limit_N_Sigma = 1, 
                                     EWMA_HWM = 0,
                                     EWMA_Variance_HWM = 0, 
                                     EWMA_Iteration_HWM = 0,
                                     Control_Target = 0) {

  # Calculate centreline for windows
  Number_windows <- ceiling(length(Input_Data)/EWMA_Window_Size)
  Window_Data_Centreline <- matrix(0,Number_windows,1)
  
  # Aggregate
  # Input data should be a multiple of window size
  for (i in 1:Number_windows) {
    
    if (i == Number_windows) {
      
      Window_Data_Centreline[i] <- mean(Input_Data[(1+(i-1)*EWMA_Window_Size):length(Input_Data)],na.rm = TRUE)
      
    } else {
      
      Window_Data_Centreline[i] <- mean(Input_Data[(1+(i-1)*EWMA_Window_Size):(i*EWMA_Window_Size)],na.rm = TRUE)
      
    }
    
  }
  # Initialise EWMA matrix
  EWMA <- matrix(0,Number_windows,1)
  EWMA[1] <- EWMA_HWM #EWMA0 is 0 if not provided or last known value
  
  # Calculate EWMA
  for (i in 2:Number_windows) {
    
    if (is.na(Window_Data_Centreline[[i-1]])) {
      
      EWMA[i] <- EWMA[i-1]
      
    } else {
      
      EWMA[i] <- EWMA_Lambda*Window_Data_Centreline[[i-1]] + (1-EWMA_Lambda)*EWMA[i-1]
      
    }
    
  }
  
  Variance <- matrix(0,Number_windows,1)
  
  if (EWMA_Variance_HWM == 0 & EWMA_Iteration_HWM == 0){
    # V1 = (xbar1 - z0)^2 = (xbar1 - xbar0)^2 (xbar0 = z0 = pop mean or target)
    Variance_Temp <- (EWMA_Lambda^2)*(Window_Data_Centreline[2] - EWMA[1])^2
    Variance[1] <- Variance_Temp
  } else if (EWMA_Variance_HWM != 0 & EWMA_Iteration_HWM != 0) {
    Variance_Temp <- EWMA_Variance_HWM
    Variance[1] <- Variance_Temp
  }
  
  
  for (i in 2:(Number_windows)) {
    
    # Missing data, new variance is same as old variance
    if (is.na(Window_Data_Centreline[[i]])) {
      
      Variance[i] <- Variance[i-1]
      
    } else {
      
      # Variance[i] <- ((1/(i+EWMA_Iteration_HWM))*((i+EWMA_Iteration_HWM)-1)*Variance[i-1]) + (1/(i+EWMA_Iteration_HWM))*((Window_Data_Centreline[i] - EWMA[i-1])^2)
      Variance[i] <- ((i+EWMA_Iteration_HWM-2)/(i+EWMA_Iteration_HWM-1))*Variance[i-1] + ((EWMA_Lambda^2)/(i+EWMA_Iteration_HWM-1))*((Window_Data_Centreline[i] - EWMA[i-1])^2)
    }
    
  }
  
  if (Control_Target == "EWMA_Mean") {
    
    LCL <- c(0,sapply(2:length(Variance), function(i) EWMA[i-1] - Control_Limit_N_Sigma*sqrt(Variance[i])))
    
  } else if (Control_Target != "EWMA_Mean") {
    
    LCL <- c(0,sapply(2:length(Variance), function(i) as.numeric(Control_Target) - Control_Limit_N_Sigma*sqrt(Variance[i])))
    
  }
    
  LCL <- rep(LCL, each = EWMA_Window_Size)
  
  return(LCL)
  
}

Rule_EWMA_Lower_Control_Limit = list(EWMA_Lower_Control_Limit = EWMA_Lower_Control_Limit)
saveRDS(Rule_EWMA_Lower_Control_Limit,"Lower_Control_Limit.Rds")

#-----------------------------------------------------------------------------------------#
# UCL


EWMA_Upper_Control_Limit <- function(Input_Data, 
                                     EWMA_Window_Size = 60, 
                                     EWMA_Lambda = 0.9, 
                                     Control_Limit_N_Sigma = 1, 
                                     EWMA_HWM = 0,
                                     EWMA_Variance_HWM = 0, 
                                     EWMA_Iteration_HWM = 0,
                                     Control_Target = 0) {
    
  # Calculate centreline for windows
  Number_windows <- ceiling(length(Input_Data)/EWMA_Window_Size)
  Window_Data_Centreline <- matrix(0,Number_windows,1)
  
  # Aggregate
  # Input data should be a multiple of window size
  for (i in 1:Number_windows) {
    
    if (i == Number_windows) {
      
      Window_Data_Centreline[i] <- mean(Input_Data[(1+(i-1)*EWMA_Window_Size):length(Input_Data)],na.rm = TRUE)
      
    } else {
      
      Window_Data_Centreline[i] <- mean(Input_Data[(1+(i-1)*EWMA_Window_Size):(i*EWMA_Window_Size)],na.rm = TRUE)
      
    }
    
  }
  # Initialise EWMA matrix
  EWMA <- matrix(0,Number_windows,1)
  EWMA[1] <- EWMA_HWM #EWMA0 is 0 if not provided or last known value
  
  # Calculate EWMA
  for (i in 2:Number_windows) {
    
    if (is.na(Window_Data_Centreline[[i-1]])) {
      
      EWMA[i] <- EWMA[i-1]
      
    } else {
      
      EWMA[i] <- EWMA_Lambda*Window_Data_Centreline[[i-1]] + (1-EWMA_Lambda)*EWMA[i-1]
      
    }
    
  }
  
  Variance <- matrix(0,Number_windows,1)
  
  if (EWMA_Variance_HWM == 0 & EWMA_Iteration_HWM == 0){
    # V1 = (xbar1 - z0)^2 = (xbar1 - xbar0)^2 (xbar0 = z0 = pop mean or target)
    Variance_Temp <- (EWMA_Lambda^2)*(Window_Data_Centreline[2] - EWMA[1])^2
    Variance[1] <- Variance_Temp
  } else if (EWMA_Variance_HWM != 0 & EWMA_Iteration_HWM != 0) {
    Variance_Temp <- EWMA_Variance_HWM
    Variance[1] <- Variance_Temp
  }
  
  
  for (i in 2:(Number_windows)) {
    
    # Missing data, new variance is same as old variance
    if (is.na(Window_Data_Centreline[[i]])) {
      
      Variance[i] <- Variance[i-1]
      
    } else {
      
      # Variance[i] <- ((1/(i+EWMA_Iteration_HWM))*((i+EWMA_Iteration_HWM)-1)*Variance[i-1]) + (1/(i+EWMA_Iteration_HWM))*((Window_Data_Centreline[i] - EWMA[i-1])^2)
      Variance[i] <- ((i+EWMA_Iteration_HWM-2)/(i+EWMA_Iteration_HWM-1))*Variance[i-1] + ((EWMA_Lambda^2)/(i+EWMA_Iteration_HWM-1))*((Window_Data_Centreline[i] - EWMA[i-1])^2)
    }
    
  }
  
  if (Control_Target == "EWMA_Mean") {
    
    UCL <- c(0,sapply(2:length(Variance), function(i) EWMA[i-1] + Control_Limit_N_Sigma*sqrt(Variance[i])))
    
  } else if (Control_Target != "EWMA_Mean") {
    
    UCL <- c(0,sapply(2:length(Variance), function(i) as.numeric(Control_Target) + Control_Limit_N_Sigma*sqrt(Variance[i])))
    
  }

  UCL <- rep(UCL, each = EWMA_Window_Size)
  
  return(UCL)
  
}

Rule_EWMA_Upper_Control_Limit = list(EWMA_Upper_Control_Limit = EWMA_Upper_Control_Limit)
saveRDS(Rule_EWMA_Upper_Control_Limit,"Upper_Control_Limit.Rds")

#-----------------------------------------------------------------------------------------#
# Counter


EWMA_Counter <- function(Input_Data, 
                         EWMA_Window_Size = 60, 
                         EWMA_Iteration_HWM = 0) {
  
  # Calculate centreline for windows
  Number_windows <- ceiling(length(Input_Data)/EWMA_Window_Size)
  Window_Data_Centreline <- matrix(0,Number_windows,1)
  
  EWMA_Counter_Sequence <- seq(EWMA_Iteration_HWM + 1, EWMA_Iteration_HWM + Number_windows, 1)
  
  EWMA_Counter_Sequence <- rep(EWMA_Counter_Sequence, each = EWMA_Window_Size)

  return(EWMA_Counter_Sequence)
  
}

Rule_EWMA_Counter = list(EWMA_Counter = EWMA_Counter)
saveRDS(EWMA_Counter,"EWMA_Counter.Rds")


#-----------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------------#
# End of code
#
#
#-----------------------------------------------------------------------------------------#


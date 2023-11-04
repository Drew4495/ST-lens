# Initialization ----
# |||||||||||||||||||

## Variables ----
## ||||||||||||||

# Room for times
times <- list()


## Directories ----
## ||||||||||||||||

# Data
if (!file.exists(directory.results)){
  dir.create(directory.results)
  if(VERBOSE){
    cat("\n- Results directory created!\n")
  }
}

# Images
if (!file.exists(directory.images)){
  dir.create(directory.images)
  if(VERBOSE){
    cat("\n- Images directory created!\n")
  }
}

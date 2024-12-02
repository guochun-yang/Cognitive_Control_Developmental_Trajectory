# List of packages you need
packages <- c("metafor", "stringr", "performance", "car", "rms", "mgcv", "MuMIn","visreg","readxl","fastDummies","MASS") # Add any packages you need


# Function to install missing packages
install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
    message(paste("Installed:", package))
  } else {
    message(paste("Already installed:", package))
  }
}

# Apply the function to each package in the list
lapply(packages, install_if_missing)

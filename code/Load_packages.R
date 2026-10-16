#---- Install and Load Required Packages
Ensure_packages <- function(packages_list) {
  # Check for missing packages
  new_pkg <- packages_list[!(packages_list %in% installed.packages()[, "Package"])]
  # Install missing packages if any
  if(length(new_pkg)) {
    message("Installing missing packages: ", paste(new_pkg, collapse = ", "))
    install.packages(new_pkg, dependencies = TRUE)
  }
  # Load all packages
  sapply(packages_list, require, character.only = TRUE)
  message("All packages loaded successfully.")
}


# List the packages needs here
required_libs <- c(
  "base",      # Standard statistical functions (e.g., qr for the QR decomposition)
  "MASS",       # Standard statistical functions (e.g., ginv for pseudo-inverse)
  "stats",      # Base statistical algorithms (p-values via pchisq/pnorm)
  "Matrix",     # Sparse matrix classes for efficient high-dimensional algebra
  "data.table", # Fast file reading (fread) and memory-efficient data handling
  "doParallel", # Parallel backend to accelerate loops and simulations
  "mvtnorm",    # Simulating data with correlation structures (Linkage Disequilibrium)
  "BGLR",       # BayesC comparator (spike-and-slab MCMC)
  "dplyr"       # Tidy data manipulation for summarizing results
)

# Call the function
#Ensure_packages(required_libs)

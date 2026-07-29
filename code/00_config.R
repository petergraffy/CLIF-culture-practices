# Shared project configuration.

source("utils/clif_io.R")

project_root <- clif_repo_path

paths <- list(
  data = project_path("data"),
  intermediate = project_intermediate_path(),
  output = project_output_path(),
  config = project_path("config")
)

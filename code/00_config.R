# Shared project configuration.

project_root <- normalizePath(getwd(), mustWork = FALSE)

paths <- list(
  data = file.path(project_root, "data"),
  output = file.path(project_root, "output"),
  config = file.path(project_root, "config")
)

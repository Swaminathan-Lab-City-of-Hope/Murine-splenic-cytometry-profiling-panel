# Run the selected analysis scripts in the documented order.
# From the repository root: Rscript run_pipeline.R
cmd <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", cmd, value = TRUE)
if (length(script_arg) != 1L) stop("Run this launcher with Rscript")
repo <- normalizePath(dirname(sub("^--file=", "", script_arg)),
                      winslash = "/", mustWork = TRUE)
Sys.setenv(CYTOMETRY_PROJECT_ROOT = repo)
steps <- c(
  "00_samples.R",
  "01_preprocess.R",
  "02_flowsom.R",
  "03_umap.R",
  "04_heatmap.R"
)
for (step in steps) {
  message("RUN | ", step)
  sys.source(file.path(repo, "scripts", step),
             envir = new.env(parent = globalenv()))
  message("DONE | ", step)
}
message("Pipeline completed. Compare outputs with manuscript figures before release.")

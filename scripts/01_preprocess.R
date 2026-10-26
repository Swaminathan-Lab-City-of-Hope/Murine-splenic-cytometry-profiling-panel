local({
  root <- normalizePath(Sys.getenv("CYTOMETRY_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
  fcs_dir <- file.path(root, "data", "fcs")
  res <- file.path(root, "results")
  meta <- read.csv(file.path(root, "metadata.csv"), stringsAsFactors = FALSE)
  
  dir.create(res, showWarnings = FALSE)
  
  stopifnot(
    as.character(packageVersion("cyCombine")) == "0.3.0",
    as.character(packageVersion("flowCore")) == "2.14.2"
  )
  
  ff <- flowCore::read.FCS(
    file.path(fcs_dir, meta$fcs_file[1]),
    transformation = FALSE,
    alter.names = FALSE,
    truncate_max_range = FALSE
  )
  
  p <- Biobase::pData(flowCore::parameters(ff))
  p <- p[!is.na(p$desc) & nzchar(trimws(p$desc)), c("name", "desc")]
  
  channels <- as.character(p$name)
  markers <- gsub("[^A-Za-z0-9.]+", "", as.character(p$desc))
  
  stopifnot(length(markers) == 31L, !anyDuplicated(markers))
  rm(ff, p)
  
  dat <- do.call(rbind, lapply(seq_len(nrow(meta)), function(i) {
    ff <- flowCore::read.FCS(
      file.path(fcs_dir, meta$fcs_file[i]),
      transformation = FALSE,
      alter.names = FALSE,
      truncate_max_range = FALSE
    )
    
    x <- as.data.frame(
      asinh(flowCore::exprs(ff)[, channels, drop = FALSE] / 6000),
      check.names = FALSE
    )
    
    names(x) <- markers
    x$sample <- meta$sample_id[i]
    x$batch <- meta$acquisition_batch[i]
    x$condition <- meta$genotype[i]
    x
  }))
  
  stopifnot(
    nrow(dat) == sum(meta$live_cell_events),
    all(markers %in% names(dat))
  )
  
  corrected <- cyCombine::batch_correct(
    dat,
    markers = markers,
    covar = "condition",
    norm_method = "scale",
    xdim = 8,
    ydim = 8,
    rlen = 10,
    seed = 473
  )
  
  stopifnot(
    nrow(corrected) == nrow(dat),
    !anyNA(corrected[, markers])
  )
  
  saveRDS(
    list(
      data = corrected,
      markers = markers,
      marker_map = setNames(markers, channels),
      settings = list(
        transformation = "asinh(x / 6000)",
        batch = "acquisition_batch",
        covariate = "genotype",
        cyCombine = "0.3.0",
        norm_method = "scale",
        som = "8x8",
        rlen = 10L,
        seed = 473L
      )
    ),
    file.path(res, "corrected_events.rds")
  )
  
  message(
    "PASS | events=", format(nrow(corrected), big.mark = ","),
    " | markers=", length(markers),
    " | samples=", length(unique(corrected$sample))
  )
})
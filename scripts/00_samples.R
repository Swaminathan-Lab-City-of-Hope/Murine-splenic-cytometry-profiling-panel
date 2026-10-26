# Metadata for the six live-gated, spectrally unmixed FCS files.

local({
  root <- normalizePath(
    Sys.getenv("CYTOMETRY_PROJECT_ROOT", unset = getwd()),
    winslash = "/",
    mustWork = TRUE
  )
  fcs_dir <- file.path(root, "data", "fcs")

  samples <- data.frame(
    sample_id = c("NZM1", "NZM2", "NZM3", "WT1", "WT2", "WT3"),
    genotype = c("NZM", "NZM", "NZM", "WT", "WT", "WT"),
    age = c("1 year", "1 year", "1 year", "12 weeks", "23 weeks", "18 weeks"),
    acquisition_batch = c("2026-08-15", "2026-08-15", "2026-08-22",
                          "2026-08-15", "2026-08-15", "2026-08-22"),
    live_cell_events = c(285501, 360100, 401068, 284829, 273346, 256721),
    fcs_file = c(
      "export_--Splenocytes-NZM1_1year_Unmixed_Live cells.fcs",
      "export_--Splenocytes-NZM2_1year_Unmixed_Live cells.fcs",
      "export_--Splenocytes-NZM3_1year_Unmixed_Live cells.fcs",
      "export_--Splenocytes-WT1_12weeks_Unmixed_Live cells.fcs",
      "export_--Splenocytes-WT2_23weeks_Unmixed_Live cells.fcs",
      "export_--Splenocytes-WT3_18weeks_Unmixed_Live cells.fcs"
    )
  )

  stopifnot(
    nrow(samples) == 6L,
    sum(samples$live_cell_events) == 1861565L,
    all(file.exists(file.path(fcs_dir, samples$fcs_file)))
  )

  write.csv(samples, file.path(root, "metadata.csv"), row.names = FALSE)
  message("Saved metadata.csv (", nrow(samples), " samples)")
})

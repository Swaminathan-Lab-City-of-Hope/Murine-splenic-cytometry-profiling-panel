local({
  root <- normalizePath(Sys.getenv("CYTOMETRY_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
  res <- file.path(root, "results")
  
  x <- readRDS(file.path(res, "corrected_events.rds"))
  dat <- x$data
  
  stopifnot(as.character(packageVersion("FlowSOM")) == "2.10.0")
  
  clustering <- c(
    "GL7","CD3","IgD","F4.80","CD5","CD49b","B220","CD8","CD19",
    "IAIE","CD93","CXCR4CD184","XCR1","IgM","TCRgd","PDCA1CD317",
    "CD23","CXCR5CD185","CD138","CD11b","LY6G","LY6C","CD4",
    "CD21.35","CD95FAS","NKp46","CD11c"
  )
  state <- c("KI67","CD86","PD1CD279","TBET")
  markers <- c(clustering, state)
  
  stopifnot(length(clustering) == 27L, all(markers %in% names(dat)))
  
  seed <- 20260908L
  set.seed(seed)
  
  idx <- unlist(lapply(split(seq_len(nrow(dat)), dat$sample),
                       \(i) sample(i, 50000L)), use.names = FALSE)
  
  fsom <- FlowSOM::ReadInput(
    flowCore::flowFrame(as.matrix(dat[idx, markers])),
    compensate = FALSE, transform = FALSE, scale = FALSE, silent = TRUE
  )
  
  set.seed(seed)
  fsom <- FlowSOM::BuildSOM(
    fsom, colsToUse = clustering,
    xdim = 10, ydim = 10, rlen = 10, silent = TRUE
  )
  fsom <- FlowSOM::BuildMST(fsom, silent = TRUE)
  
  node <- integer(nrow(dat))
  
  for (s in unique(dat$sample)) {
    i <- which(dat$sample == s)
    mapped <- FlowSOM::NewData(
      fsom,
      flowCore::flowFrame(as.matrix(dat[i, markers])),
      compensate = FALSE, transform = FALSE, scale = FALSE, silent = TRUE
    )
    node[i] <- FlowSOM::GetClusters(mapped)
  }
  
  stopifnot(length(node) == nrow(dat), all(node %in% 1:100))
  
  node_to_meta <- FlowSOM::metaClustering_consensus(
    as.matrix(fsom$map$codes),
    k = 25,
    seed = seed
  )
  
  meta <- node_to_meta[node]
  
  stopifnot(
    length(node_to_meta) == 100L,
    length(unique(node_to_meta)) == 25L,
    length(meta) == nrow(dat)
  )
  
  saveRDS(
    list(
      model = fsom,
      node = node,
      node_to_metacluster = node_to_meta,
      cell_metacluster = meta,
      clustering_markers = clustering,
      state_markers = state,
      training_indices = idx,
      settings = list(
        training_per_mouse = 50000L,
        grid = "10x10",
        rlen = 10L,
        scale = FALSE,
        seed = seed,
        K = 25L
      )
    ),
    file.path(res, "flowsom_K25.rds")
  )
  
  message(
    "PASS | training=", format(length(idx), big.mark = ","),
    " | mapped=", format(length(node), big.mark = ","),
    " | nodes=", length(unique(node)),
    " | K=", length(unique(meta))
  )
})
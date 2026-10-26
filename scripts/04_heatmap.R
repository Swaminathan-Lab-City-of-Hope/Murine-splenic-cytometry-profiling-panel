local({
  
  root <- normalizePath(Sys.getenv("CYTOMETRY_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
  res  <- file.path(root, "results")
  
  dat <- readRDS(file.path(res, "corrected_events.rds"))$data
  fs  <- readRDS(file.path(res, "flowsom_K25.rds"))
  pal <- readRDS(file.path(res, "03_K25_UMAP.rds"))$palette
  
  mc <- fs$cell_metacluster
  
  markers <- c(
    "CD3","CD19","B220","NKp46","CD49b","CD11b","LY6G","LY6C",
    "F4.80","CD11c","IAIE","XCR1","PDCA1CD317","CD4","CD8",
    "TCRgd","CD5","IgD","IgM","CD21.35","CD23","CD93","GL7",
    "CD95FAS","CD138","CXCR5CD185","CXCR4CD184",
    "KI67","CD86","PD1CD279","TBET"
  )
  
  labels <- c(
    "CD3","CD19","B220","NKp46","CD49b","CD11b","Ly6G","Ly6C",
    "F4/80","CD11c","MHC-II","XCR1","PDCA-1","CD4","CD8",
    "TCRγδ","CD5","IgD","IgM","CD21/35","CD23","CD93","GL7",
    "CD95/FAS","CD138","CXCR5","CXCR4","Ki-67","CD86","PD-1","T-bet"
  )
  
  stopifnot(
    length(unique(mc)) == 25L,
    all(markers %in% names(dat)),
    setequal(pal$cluster, 1:25)
  )
  
  n <- tabulate(mc, 25L)
  
  m <- vapply(
    markers,
    \(x) as.numeric(rowsum(dat[[x]], factor(mc, levels = 1:25))) / n,
    numeric(25)
  )
  
  rownames(m) <- paste0("Cluster ", 1:25)
  colnames(m) <- labels
  
  # Marker-wise z-score across K25
  z <- scale(m)
  
  # Dendrogram across markers
  hc  <- hclust(dist(t(z)), method = "complete")
  ord <- hc$order
  
  hm <- z[, ord, drop = FALSE]
  hm <- pmax(pmin(hm, 2), -2)
  
  pal <- pal[match(1:25, pal$cluster), ]
  heat_pal <- colorRampPalette(c("blue", "white", "red"))(201)
  
  write.csv(
    data.frame(cluster = 1:25, m, check.names = FALSE),
    file.path(res, "04_K25_mean_expression.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(cluster = 1:25, z, check.names = FALSE),
    file.path(res, "04_K25_scaled_mean_expression.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      position = 1:31,
      marker = colnames(z)[ord]
    ),
    file.path(res, "04_K25_marker_dendrogram_order.csv"),
    row.names = FALSE
  )
  
  png(
    file.path(res, "04_K25_marker_heatmap.png"),
    width = 190,
    height = 150,
    units = "mm",
    res = 600,
    type = "cairo"
  )
  
  layout(
    matrix(c(
      0, 2, 0,
      1, 3, 4
    ), nrow = 2, byrow = TRUE),
    widths  = c(4.0, 10.0, 1.2),
    heights = c(2.0, 8.0)
  )
  
  # ============================================================
  # CLUSTER LABELS + EXACT UMAP COLORS
  # ============================================================
  
  par(mar = c(8.5, 0, 0.3, 0), family = "sans")
  
  plot.new()
  plot.window(
    xlim = c(0, 1),
    ylim = c(0.5, 25.5),
    xaxs = "i",
    yaxs = "i"
  )
  
  y <- 25:1
  
  text(
    rep(0.84, 25),
    y,
    paste("Cluster", 1:25),
    adj = 1,
    cex = 1.05
  )
  
  points(
    rep(0.91, 25),
    y,
    pch = 16,
    cex = 1.45,
    col = pal$color
  )
  
  # ============================================================
  # MARKER DENDROGRAM
  # Same horizontal extent as heatmap
  # ============================================================
  
  par(
    mar = c(0, 0, 0.4, 0),
    family = "sans",
    xaxs = "i"
  )
  
  plot(
    as.dendrogram(hc),
    leaflab = "none",
    axes = FALSE,
    xlab = "",
    ylab = "",
    xaxs = "i"
  )
  
  # ============================================================
  # HEATMAP
  # ============================================================
  
  par(
    mar = c(8.5, 0, 0, 0),
    family = "sans",
    xaxs = "i",
    yaxs = "i"
  )
  
  image(
    x = 1:31,
    y = 1:25,
    z = t(hm[25:1, , drop = FALSE]),
    col = heat_pal,
    breaks = seq(-2, 2, length.out = 202),
    axes = FALSE,
    xlab = "",
    ylab = "",
    xaxs = "i",
    yaxs = "i"
  )
  
  axis(
    1,
    at = 1:31,
    labels = colnames(hm),
    las = 2,
    tick = FALSE,
    cex.axis = 0.90,
    line = -0.15
  )
  
  abline(
    v = seq(1.5, 30.5, 1),
    h = seq(1.5, 24.5, 1),
    col = adjustcolor("white", 0.20),
    lwd = 0.30
  )
  
  box(lwd = 0.8)
  
  # ============================================================
  # COLOR SCALE
  # ============================================================
  
  par(mar = c(8.5, 0, 0, 0), family = "sans")
  
  plot.new()
  plot.window(
    xlim = c(0, 1),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i"
  )
  
  bottom <- 0.37
  top    <- 0.63
  left   <- 0.08
  right  <- 0.30
  
  yy <- seq(bottom, top, length.out = length(heat_pal) + 1)
  
  rect(
    left,
    yy[-length(yy)],
    right,
    yy[-1],
    col = heat_pal,
    border = NA
  )
  
  rect(
    left,
    bottom,
    right,
    top,
    border = "black",
    lwd = 0.7
  )
  
  ticks <- -2:2
  ty <- bottom + ((ticks + 2) / 4) * (top - bottom)
  
  segments(
    x0 = rep(right, length(ty)),
    y0 = ty,
    x1 = rep(0.37, length(ty)),
    y1 = ty,
    lwd = 0.7
  )
  
  text(
    x = rep(0.42, length(ty)),
    y = ty,
    labels = ticks,
    adj = 0,
    cex = 0.9
  )
  
  text(
    x = 0.69,
    y = mean(c(bottom, top)),
    labels = "Scaled mean marker expression",
    srt = 90,
    cex = 0.95
  )
  
  dev.off()
  
  message(
    "PASS | K25 rows | 31 marker columns",
    " | marker dendrogram=TRUE",
    " | cluster dendrogram=FALSE",
    " | output overwritten"
  )
})
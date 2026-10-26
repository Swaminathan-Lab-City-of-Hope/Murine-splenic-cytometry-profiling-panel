# UMAP and figures
local({
  root <- normalizePath(Sys.getenv("CYTOMETRY_PROJECT_ROOT", unset = getwd()),
                        winslash = "/", mustWork = TRUE)
  res <- file.path(root, "results")
  stopifnot(dir.exists(res))

  x <- readRDS(file.path(res, "corrected_events.rds"))
  fs <- readRDS(file.path(res, "flowsom_K25.rds"))
  meta <- read.csv(file.path(root, "metadata.csv"), stringsAsFactors = FALSE)
  dat <- x$data
  markers <- fs$clustering_markers
  seed <- 20260908L

  stopifnot(
    as.character(packageVersion("uwot")) == "0.2.4",
    as.character(packageVersion("ggplot2")) == "4.0.3",
    nrow(dat) == length(fs$cell_metacluster),
    length(markers) == 27L,
    all(markers %in% names(dat)),
    all(c("sample_id", "genotype") %in% names(meta))
  )

  # UMAP
  set.seed(seed)
  idx <- unlist(lapply(split(seq_len(nrow(dat)), dat$sample),
                       \(i) sample(i, 20000L)), use.names = FALSE)

  pca <- prcomp(as.matrix(dat[idx, markers]), center = TRUE, scale. = FALSE)
  set.seed(seed)
  u <- uwot::umap(pca$x, n_neighbors = 15, min_dist = 0.2,
                  metric = "euclidean", n_components = 2, verbose = TRUE)

  d <- data.frame(
    event_index = idx,
    sample = dat$sample[idx],
    cluster = factor(fs$cell_metacluster[idx], levels = 1:25),
    UMAP1 = u[, 1],
    UMAP2 = u[, 2]
  )
  stopifnot(nrow(d) == 120000L, length(unique(d$cluster)) == 25L)

  d$genotype <- meta$genotype[match(d$sample, meta$sample_id)]
  if (anyNA(d$genotype)) stop("UMAP sample IDs do not match metadata.csv")
  d$genotype <- factor(d$genotype, levels = c("WT", "NZM"))

  sample_order <- c("WT1", "WT2", "WT3", "NZM1", "NZM2", "NZM3")
  stopifnot(all(table(factor(d$sample, levels = sample_order)) == 20000L))

  # Cluster labels
  pos <- do.call(rbind, lapply(1:25, function(k) {
    z <- d[d$cluster == k, c("UMAP1", "UMAP2")]
    cx <- median(z$UMAP1)
    cy <- median(z$UMAP2)
    j <- which.min((z$UMAP1 - cx)^2 + (z$UMAP2 - cy)^2)
    data.frame(cluster = k, n = nrow(z), x = z$UMAP1[j], y = z$UMAP2[j])
  }))

  xmid <- median(d$UMAP1)
  ymid <- median(d$UMAP2)
  xr <- diff(range(d$UMAP1))
  yr <- diff(range(d$UMAP2))
  vx <- pos$x - xmid
  vy <- pos$y - ymid
  norm <- sqrt(vx^2 + vy^2)
  norm[norm == 0] <- 1
  small <- pos$n < 2500L
  pos$label_x <- pos$x + ifelse(small, vx / norm * 0.030 * xr, 0)
  pos$label_y <- pos$y + ifelse(small, vy / norm * 0.030 * yr, 0)

  # Cluster palette
  initial_colors <- grDevices::hcl.colors(25, "Dynamic")
  xy <- as.matrix(pos[order(pos$cluster), c("x", "y")])
  ds <- as.matrix(dist(xy))
  adj <- matrix(0, 25, 25)
  for (i in 1:25) {
    j <- order(ds[i, ])[2:5]
    adj[i, j] <- 1 / ds[i, j]
  }
  adj <- pmax(adj, t(adj))
  rgb <- t(grDevices::col2rgb(initial_colors)) / 255
  lab <- grDevices::convertColor(rgb, from = "sRGB", to = "Lab")
  dc <- as.matrix(dist(lab))
  score <- function(p) sum(adj[upper.tri(adj)] * dc[p, p][upper.tri(dc)])

  set.seed(seed)
  p <- seq_len(25)
  best <- score(p)
  for (z in 1:50000) {
    ij <- sample.int(25, 2)
    q <- p
    q[ij] <- q[rev(ij)]
    s <- score(q)
    if (s > best) {
      p <- q
      best <- s
    }
  }
  palette <- data.frame(cluster = 1:25, color = initial_colors[p])
  cols <- setNames(palette$color, as.character(palette$cluster))

  # Figure 3A
  set.seed(seed)
  pdat <- d[sample.int(nrow(d)), ]
  moved <- abs(pos$label_x - pos$x) > 1e-8 |
           abs(pos$label_y - pos$y) > 1e-8

  g <- ggplot2::ggplot(pdat, ggplot2::aes(UMAP1, UMAP2, color = cluster)) +
    ggplot2::geom_point(size = 0.22, alpha = 0.80, stroke = 0) +
    ggplot2::scale_color_manual(values = cols, guide = "none") +
    ggplot2::geom_segment(
      data = pos[moved, ],
      ggplot2::aes(x = x, y = y, xend = label_x, yend = label_y),
      inherit.aes = FALSE, linewidth = 0.25, color = "black"
    ) +
    ggplot2::geom_text(
      data = pos,
      ggplot2::aes(label_x, label_y, label = cluster),
      inherit.aes = FALSE, color = "black", size = 1.8,
      fontface = "bold"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "UMAP1", y = "UMAP2") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 13, face = "bold"),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.7)
    )
  ggplot2::ggsave(
    file.path(res, "03B_K25_numbered_UMAP_spatial_palette.png"),
    g, width = 105, height = 105, units = "mm", dpi = 600, bg = "white"
  )

  # Supplementary UMAP panels
  cluster_id <- as.integer(as.character(d$cluster))
  cell_cols <- unname(cols[as.character(cluster_id)])
  stopifnot(!anyNA(cell_cols))

  n_cluster <- table(cluster_id)
  small_clusters <- as.integer(names(n_cluster[n_cluster / length(cluster_id) < 0.02]))
  x_range <- range(d$UMAP1)
  y_range <- range(d$UMAP2)
  dx <- diff(x_range)
  dy <- diff(y_range)
  x_plot <- x_range + c(-0.06, 0.06) * dx
  y_plot <- y_range + c(-0.06, 0.06) * dy
  global_x <- median(d$UMAP1)
  global_y <- median(d$UMAP2)

  get_centers <- function(index) {
    z <- data.frame(UMAP1 = d$UMAP1[index], UMAP2 = d$UMAP2[index],
                    cluster = cluster_id[index])
    out <- aggregate(cbind(UMAP1, UMAP2) ~ cluster, data = z, FUN = median)
    out[order(out$cluster), ]
  }

  add_cluster_labels <- function(centers) {
    is_small <- centers$cluster %in% small_clusters
    if (any(!is_small)) {
      text(centers$UMAP1[!is_small], centers$UMAP2[!is_small],
           labels = centers$cluster[!is_small],
           cex = 0.48, font = 2, col = "black")
    }
    if (any(is_small)) {
      sm <- centers[is_small, , drop = FALSE]
      vx <- (sm$UMAP1 - global_x) / dx
      vy <- (sm$UMAP2 - global_y) / dy
      radius <- sqrt(vx^2 + vy^2)
      near_center <- radius < 0.03
      if (any(near_center)) {
        theta0 <- (sm$cluster[near_center] * 137.507764) * pi / 180
        vx[near_center] <- cos(theta0)
        vy[near_center] <- sin(theta0)
      }
      theta <- atan2(vy, vx)
      label_x <- sm$UMAP1 + 0.038 * dx * cos(theta)
      label_y <- sm$UMAP2 + 0.038 * dy * sin(theta)
      line_x <- sm$UMAP1 + 0.030 * dx * cos(theta)
      line_y <- sm$UMAP2 + 0.030 * dy * sin(theta)
      segments(sm$UMAP1, sm$UMAP2, line_x, line_y,
               col = "black", lwd = 0.45, xpd = NA)
      text(label_x, label_y, labels = sm$cluster,
           cex = 0.46, font = 2, col = "black", xpd = NA)
    }
  }

  draw_panel <- function(index, title_text) {
    centers <- get_centers(index)
    plot(d$UMAP1[index], d$UMAP2[index], pch = 16, cex = 0.14,
         col = cell_cols[index], axes = FALSE, xlab = "", ylab = "",
         xlim = x_plot, ylim = y_plot, asp = 1)
    add_cluster_labels(centers)
    title(title_text, cex.main = 0.90, font.main = 2, line = 0.20)
  }

  # Supplementary Figure 7
  grDevices::png(file.path(res, "Supplementary_Figure_7_UMAP.png"),
                 width = 210, height = 235, units = "mm", res = 600,
                 type = "cairo")
  layout(matrix(c(
    1, 1, 1, 2, 2, 2,
    3, 3, 4, 4, 5, 5,
    6, 6, 7, 7, 8, 8
  ), nrow = 3, byrow = TRUE),
  widths = rep(1, 6), heights = c(1.5, 1, 1))
  par(family = "sans", mar = c(0.15, 0.15, 1.25, 0.15),
      oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")

  draw_panel(d$genotype == "WT", "WT")
  draw_panel(d$genotype == "NZM", "NZM")
  for (sid in sample_order) draw_panel(d$sample == sid, sid)
  grDevices::dev.off()

  # UMAP results
  saveRDS(
    list(
      coordinates = d,
      palette = palette,
      label_positions = pos,
      settings = list(
        cells_per_mouse = 20000L, total_cells = 120000L,
        features = markers, PCA_center = TRUE, PCA_scale = FALSE,
        PCs = 27L, n_neighbors = 15L, min_dist = 0.2,
        metric = "euclidean", seed = seed,
        palette_optimization_swaps = 50000L,
        figure_3A_small_cluster_cutoff = 2500L,
        supplementary_small_cluster_fraction = 0.02
      )
    ),
    file.path(res, "03_K25_UMAP.rds")
  )

  message("PASS | UMAP generated once | cells=120,000 | clusters=25",
          " | Figure 3A saved | Supplementary Figure 7 saved")
})

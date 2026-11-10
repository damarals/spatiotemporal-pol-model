# Figures 1 (study area), 2 (source coverage), 3 (TROPOMI column checks) and 4
# (concentration maps), in the theme of the paper's figures: Poppins,
# theme_light panels with a black border, maps on theme_void with the plasma
# palette and a horizontal colour bar.

paper_font <- "Poppins"
paper_colours <- list(text = "#2D3748", frame = "#4A5568", reference = "#7F8A99",
  series = c(Observed = "#C05640", Replicated = "#2B6CB0"))

#' Figure 1: study area with the reference-station deployments and the
#' low-cost devices with field data.
#'
#' @param dir data directory with processed/, decisions/ and moqa/.
#' @param path PDF file.
#' @param files data_files(dir), possibly without names.
#' @return path, invisibly.
plot_study_area <- function(dir, path, files = data_files(dir)) {
  .save_pdf(.study_area_plot(.data_paths(files)), path, width = 8, height = 7)
}

#' Figure 2: months with data for each source, 2020-2025.
#'
#' @param input read_inputs().
#' @param path PDF file.
#' @return path, invisibly.
plot_source_coverage <- function(input, path) {
  .save_pdf(.source_coverage_plot(input), path, width = 180 / 25.4, height = 46 / 25.4)
}

#' Figure 3: medians and central 95 % intervals of the observed and replicated
#' TROPOMI statistics, six panels; (a, d) and (e, f) share their x scales.
#'
#' @param summary column_checks()$summary.
#' @param path PDF file.
#' @return path, invisibly.
plot_column_checks <- function(summary, path) {
  value_cols <- c("observed_median", "observed_lower", "observed_upper",
    "replicated_median", "replicated_lower", "replicated_upper")
  s <- summary[match(names(.column_panels), summary$statistic), c("statistic", value_cols)]
  if (anyNA(s$statistic) || !all(is.finite(as.matrix(s[value_cols]))))
    stop("The summary must hold finite values for the six statistics.", call. = FALSE)
  long <- function(series, prefix) data.frame(statistic = s$statistic, series = series,
    median = s[[paste0(prefix, "_median")]], lower = s[[paste0(prefix, "_lower")]],
    upper = s[[paste0(prefix, "_upper")]], stringsAsFactors = FALSE)
  d <- rbind(long("Observed", "observed"), long("Replicated", "replicated"))
  figure <- patchwork::wrap_plots(lapply(names(.column_panels), .column_panel_plot, d = d),
    nrow = 2L, ncol = 3L)
  .save_pdf(figure, path, width = 7.4, height = 4.6)
}

#' Figure 4: posterior mean of the 2020-2025 mean and the width of its central
#' 95 % interval relative to the mean, for background NO2 and PM2.5, on the
#' cells clipped to the municipality, with the 22 locations as points.
#'
#' @param cells map_cells().
#' @param sites read_inputs()$sites (x_m, y_m in EPSG:31984).
#' @param boundary path of the municipal boundary (data/processed/boundary.gpkg).
#' @param path PDF file.
#' @return path, invisibly.
plot_maps <- function(cells, sites, boundary, path) {
  if (anyDuplicated(cells[c("cell", "pollutant")])) stop("Cells must be unique per pollutant.", call. = FALSE)
  for (p in c("NO2", "PM25"))
    if (!setequal(cells$cell[cells$pollutant == p], cells$cell))
      stop("NO2 and PM2.5 must cover the same cells.", call. = FALSE)
  keys <- c("cell", "lon", "lat", "pollutant")
  values <- rbind(data.frame(cells[keys], quantity = "mean", value = cells$mean),
    data.frame(cells[keys], quantity = "relative_width", value = relative_width(cells)))
  if (!all(is.finite(values$value))) stop("Map values must be finite.", call. = FALSE)
  outline <- sf::st_transform(sf::st_union(sf::st_read(boundary, quiet = TRUE)), 4326)
  points <- sf::st_transform(sf::st_as_sf(sites[c("location_id", "x_m", "y_m")], coords = c("x_m", "y_m"),
    crs = 31984), 4326)
  squares <- .cell_squares(values, outline)
  rows <- lapply(c("NO2", "PM25"), function(p) {
    pair <- patchwork::wrap_plots(.map_panel(squares, outline, points, p, "mean"),
      .map_panel(squares, outline, points, p, "relative_width"), nrow = 1L)
    patchwork::wrap_elements(plot = pair) + ggplot2::labs(title = .pollutant_label[[p]]) +
      ggplot2::theme(plot.title = ggplot2::element_text(family = paper_font, size = 11.5, hjust = 0,
        margin = ggplot2::margin(b = 0)))
  })
  .save_pdf(patchwork::wrap_plots(rows, nrow = 1L), path, width = 7.2, height = 2.9)
}

# Panel letter, title, x-axis title, dashed reference value and the statistics
# whose intervals set the shared x scale.
.column_panels <- list(
  dispersion = list("a", "Residual dispersion", "Mean squared\nstandardized residual", 1,
    c("dispersion", "contrast_dispersion")),
  maximum = list("b", "Maximum residual", "Maximum absolute\nstandardized residual", NULL, "maximum"),
  monthly = list("c", "Monthly persistence", quote((10^-5 ~ mol ~ m^-2)^2), NULL, "monthly"),
  contrast_dispersion = list("d", "Contrast dispersion", "Mean squared\nstandardized contrast", 1,
    c("dispersion", "contrast_dispersion")),
  within_cell = list("e", "Within-cell persistence", "Standardized\ncontrast product", 0,
    c("within_cell", "neighbor")),
  neighbor = list("f", "Neighbor similarity", "Standardized\ncontrast product", 0,
    c("within_cell", "neighbor")))

.column_panel_plot <- function(d, key) {
  p <- .column_panels[[key]]
  rows <- d[d$statistic == key, ]
  group <- d[d$statistic %in% p[[5L]], ]
  span <- range(group$lower, group$upper, p[[4L]])
  lim <- span + c(-1, 1) * 0.06 * diff(span)
  ticks <- pretty(lim, n = if (key == "monthly") 3L else 5L)
  ticks <- ticks[ticks >= lim[1L] & ticks <= lim[2L]]
  out <- ggplot2::ggplot(rows, ggplot2::aes(y = series, colour = series))
  if (!is.null(p[[4L]]))
    out <- out + ggplot2::geom_vline(xintercept = p[[4L]], colour = paper_colours$reference,
      linetype = "22", linewidth = 0.45)
  out +
    ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper, yend = series), linewidth = 1.3,
      lineend = "butt") +
    ggplot2::geom_point(ggplot2::aes(x = median, shape = series), size = 2.9) +
    ggplot2::scale_x_continuous(limits = lim, breaks = ticks, expand = c(0, 0)) +
    ggplot2::scale_y_discrete(limits = c("Replicated", "Observed"), expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_colour_manual(values = paper_colours$series, guide = "none") +
    ggplot2::scale_shape_manual(values = c(Observed = 16, Replicated = 17), guide = "none") +
    ggplot2::labs(title = paste0("(", p[[1L]], ") ", p[[2L]]), x = p[[3L]], y = NULL) +
    .results_theme(9) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 8.5, colour = paper_colours$text),
      axis.title.x = ggplot2::element_text(size = 8.5, lineheight = 0.95, margin = ggplot2::margin(t = 4)),
      plot.margin = ggplot2::margin(4, 8, 4, 4))
}

.pollutant_label <- list(NO2 = quote("Background NO"[2]), PM25 = quote("PM"[2.5]))

# The 0.01 degree squares around the cell centres (planar in longitude and
# latitude), clipped to the municipality.
.cell_squares <- function(values, outline) {
  squares <- lapply(seq_len(nrow(values)), function(i) {
    x <- values$lon[i] + c(-1, 1, 1, -1, -1) * 0.005
    y <- values$lat[i] + c(-1, -1, 1, 1, -1) * 0.005
    sf::st_polygon(list(cbind(x, y)))
  })
  cells <- sf::st_sf(values, geometry = sf::st_sfc(squares, crs = 4326))
  old <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  clipped <- suppressWarnings(sf::st_intersection(cells, outline))
  if (nrow(clipped) != nrow(cells) || !identical(clipped$value, cells$value))
    stop("Clipping to the municipality dropped or reordered cells.", call. = FALSE)
  clipped
}

.map_panel <- function(squares, outline, points, pollutant, quantity) {
  d <- squares[squares$pollutant == pollutant & squares$quantity == quantity, ]
  limits <- range(d$value)
  bar <- if (quantity == "mean") bquote(.(.pollutant_label[[pollutant]]) ~ "(µg" ~ m^-3 * ")") else
    "95% interval width / mean"
  ggplot2::ggplot(d) +
    ggplot2::geom_sf(ggplot2::aes(fill = value, colour = value), linewidth = 0.05) +
    ggplot2::geom_sf(data = outline, fill = NA, colour = paper_colours$text, linewidth = 0.3,
      inherit.aes = FALSE) +
    ggplot2::geom_sf(data = points, shape = 21, size = 1.15, stroke = 0.35, fill = "white",
      colour = paper_colours$text, inherit.aes = FALSE) +
    ggplot2::scale_fill_viridis_c(option = "plasma", limits = limits, name = bar,
      guide = ggplot2::guide_colourbar(theme = ggplot2::theme(legend.key.width = grid::unit(28, "mm"),
        legend.key.height = grid::unit(2.2, "mm")))) +
    ggplot2::scale_colour_viridis_c(option = "plasma", limits = limits, guide = "none") +
    ggplot2::coord_sf(expand = FALSE, datum = NA) +
    ggplot2::labs(title = if (quantity == "mean") "Posterior mean" else "Relative 95% width") +
    .map_theme(10)
}

# Map frame of Figure 1 in EPSG:31984 (Fortaleza and about 1 km around it) and
# the arrow curvature of each relocation leg P00->P01 ... P07->P08, both set by eye.
.study_area_frame <- c(xmin = 538770.8, ymin = 9568150.8, xmax = 567473.6, ymax = 9593483.9)
.relocation_curvature <- c(0, 0, 0.15, 0.15, 0.15, 0.30, 0.15, 0)
.fortaleza_code <- 2304400

.study_area_plot <- function(files) {
  read <- function(name) utils::read.csv(files[[name]], stringsAsFactors = FALSE)
  sites <- rbind(read("moqa_sites"), read("sites"))
  deployments <- read("deployments")
  moqa <- read("moqa")
  layer <- function(name) sf::st_transform(
    sf::st_read(files[["basemap"]], layer = name, quiet = TRUE), 31984)
  munis <- layer("municipalities")
  urban <- layer("urban_areas")

  old_s2 <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  urban <- urban[munis, , op = sf::st_intersects]
  city <- munis[munis$code_muni == .fortaleza_code, ]
  frame <- sf::st_bbox(.study_area_frame, crs = sf::st_crs(31984))

  site_xy <- function(ids) {
    i <- match(ids, sites$location_id)
    if (anyNA(i)) stop("Locations missing from sites.csv: ", toString(ids[is.na(i)]), call. = FALSE)
    as.matrix(sites[i, c("x_m", "y_m")])
  }
  deployments <- deployments[order(deployments$station_id), ]
  route <- site_xy(deployments$location_id)
  points <- function(xy) sf::st_as_sf(as.data.frame(xy), coords = c("x_m", "y_m"), crs = 31984)
  reference <- points(unique(route))
  devices <- points(site_xy(unique(moqa$location_id)))
  legs <- .relocation_legs(route)

  labels <- .municipality_labels(munis, frame)
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = munis, fill = "#ffffff", colour = "#c89ab3", linewidth = 0.3) +
    ggplot2::geom_sf(data = urban, fill = "#ddc7d5", colour = NA) +
    ggplot2::geom_sf(data = city, ggplot2::aes(linetype = "Fortaleza"), fill = NA,
      colour = "grey35", linewidth = 0.8, alpha = 0.5, key_glyph = "path") +
    ggplot2::scale_linetype_manual(values = c(Fortaleza = "solid"), name = NULL) +
    ggplot2::geom_text(data = labels, ggplot2::aes(x, y, label = label),
      size = 4.4, colour = "grey50", fontface = "bold", family = paper_font) +
    lapply(split(legs, legs$curvature), function(d) ggplot2::geom_curve(
      data = d, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      colour = "#2c5282", linewidth = 0.35, alpha = 0.7, linetype = "dashed",
      curvature = d$curvature[1], arrow = grid::arrow(length = grid::unit(5, "pt"), type = "closed"))) +
    ggplot2::geom_sf(data = devices, ggplot2::aes(shape = "device"),
      colour = "#9b2226", size = 3.3) +
    ggplot2::geom_sf(data = reference, ggplot2::aes(shape = "reference"),
      colour = "#2c5282", size = 3.5) +
    ggplot2::annotate("text", x = route[1, 1] + 650, y = route[1, 2], label = "SEUMA",
      size = 4.7, colour = "grey30", family = paste(paper_font, "SemiBold"), hjust = 0) +
    ggplot2::scale_shape_manual(values = c(reference = 15, device = 16),
      breaks = c("reference", "device"),
      labels = c(reference = "Mobile reference station", device = "Selected low-cost sensor"),
      name = NULL) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(order = 1,
        override.aes = list(colour = c("#2c5282", "#9b2226"), size = 4.8)),
      linetype = ggplot2::guide_legend(order = 2,
        override.aes = list(colour = "grey35", linewidth = 0.8))) +
    ggplot2::coord_sf(crs = sf::st_crs(31984), xlim = frame[c("xmin", "xmax")],
      ylim = frame[c("ymin", "ymax")], expand = FALSE) +
    ggplot2::layer(geom = .GeomMapDecor, stat = "identity", position = "identity",
      data = data.frame(x = NA), mapping = NULL, inherit.aes = FALSE, show.legend = FALSE,
      params = list(family = paper_font)) +
    ggplot2::theme_light(base_size = 22, base_family = paper_font) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "#e4ecf1"),
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0, 0.06),
      legend.justification = c(0, 0),
      legend.box = "vertical",
      legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.margin = ggplot2::margin(2, 10, 2, 10),
      legend.key.spacing.y = grid::unit(0.2, "cm"),
      legend.spacing.y = grid::unit(0.1, "cm"),
      legend.text = ggplot2::element_text(size = 14.5),
      legend.key.height = grid::unit(0.5, "cm"),
      legend.key.width = grid::unit(0.7, "cm")
    )
}

# Consecutive deployments, shortened by about 500 m at both ends so the arrows clear the squares.
.relocation_legs <- function(route) {
  n <- nrow(route)
  if (n - 1L != length(.relocation_curvature))
    stop("Expected ", length(.relocation_curvature) + 1L, " deployments, got ", n, ".", call. = FALSE)
  dx <- route[-1, 1] - route[-n, 1]
  dy <- route[-1, 2] - route[-n, 2]
  f <- pmin(500 / sqrt(dx^2 + dy^2), 0.35)
  data.frame(x = route[-n, 1] + dx * f, y = route[-n, 2] + dy * f,
    xend = route[-1, 1] - dx * f, yend = route[-1, 2] - dy * f, curvature = .relocation_curvature)
}

# Names of the municipalities whose centroid falls inside the frame, two moved off the legend.
.municipality_labels <- function(munis, frame) {
  centre <- sf::st_centroid(sf::st_geometry(munis))
  inside <- lengths(sf::st_intersects(centre, sf::st_as_sfc(frame))) > 0
  xy <- sf::st_coordinates(centre[inside])
  out <- data.frame(x = xy[, 1], y = xy[, 2], label = toupper(munis$name_muni[inside]))
  moved <- out$label == "MARACANAÚ"
  out$x[moved] <- out$x[moved] + 650
  out$y[moved] <- out$y[moved] + 3000
  moved <- out$label == "FORTALEZA"
  out$x[moved] <- out$x[moved] - 1000
  out
}

# Scale bar (bottom left) and north arrow (top right), drawn as ggspatial's
# annotation_scale(width_hint = 0.4) and north_arrow_fancy_orienteering().
.GeomMapDecor <- ggplot2::ggproto("GeomMapDecor", ggplot2::Geom,
  draw_panel = function(data, panel_params, coord, family = "") {
    pad <- grid::unit(0.25, "cm")
    bar <- .scale_bar(panel_params)
    height <- grid::unit(0.38, "cm")
    left <- pad + grid::unit((seq_len(bar$divisions) - 1) * bar$division, "npc")
    scale <- grid::gList(
      grid::rectGrob(x = left, y = pad, width = grid::unit(bar$division, "npc"), height = height,
        hjust = 0, vjust = 0, gp = grid::gpar(fill = rep(c("black", "white"),
          length.out = bar$divisions), col = "black", lwd = 1)),
      grid::textGrob(bar$label, x = pad + grid::unit(0.22, "cm") + grid::unit(bar$width, "npc"),
        y = pad + 0.5 * height, hjust = 0, vjust = 0.5,
        gp = grid::gpar(cex = 1.3, col = "black", fontfamily = family)))
    size <- grid::unit(1.9, "cm")
    arrow <- grid::gTree(
      children = grid::gList(
        grid::circleGrob(x = 0.505, y = 0.4, r = 0.3, default.units = "npc",
          gp = grid::gpar(fill = NA, lwd = 1, col = "black")),
        grid::polygonGrob(x = c(0.25, 0.5, 0.5, 0.75, 0.5, 0.5), y = c(0.1, 0.8, 0.3, 0.1, 0.8, 0.3),
          id = c(1, 1, 1, 2, 2, 2), default.units = "npc",
          gp = grid::gpar(lwd = 1, col = "black", fill = c("white", "black"))),
        grid::textGrob("N", x = 0.5, y = 0.95,
          gp = grid::gpar(fontfamily = family, fontsize = 15, col = "black"))),
      vp = grid::viewport(x = grid::unit(1, "npc") - pad - 0.5 * size,
        y = grid::unit(1, "npc") - pad - 0.5 * size, width = size, height = size))
    grid::gTree(children = grid::gList(scale, arrow))
  }
)

# Round scale-bar length (ggspatial's rule) for 40% of the frame width at mid height, in km.
.scale_bar <- function(panel_params) {
  x <- panel_params$x_range
  y <- mean(panel_params$y_range)
  ends <- sf::st_transform(sf::st_sfc(sf::st_point(c(x[1], y)), sf::st_point(c(x[2], y)),
    crs = panel_params$crs), 4326)
  ll <- sf::st_coordinates(ends) * pi / 180
  a <- sin(diff(ll[, 2]) / 2)^2 + cos(ll[1, 2]) * cos(ll[2, 2]) * sin(diff(ll[, 1]) / 2)^2
  width_km <- 6371.009 * 2 * asin(min(1, sqrt(a)))
  power <- floor(log10(0.4 * width_km))
  lead <- floor(0.4 * width_km / 10^power)
  if (lead == 1) {
    lead <- 10
    power <- power - 1
  } else if (lead %in% c(7, 9)) {
    lead <- lead - 1
  }
  length_km <- lead * 10^power
  division <- (if (lead < 6) 1 else 2) * 10^power
  list(width = length_km / width_km, division = division / width_km,
    divisions = round(length_km / division), label = paste(length_km, "km"))
}

.source_coverage_plot <- function(input) {
  calendar <- format(seq(as.Date("2020-01-01"), as.Date("2025-12-01"), by = "month"), "%Y-%m")
  targets <- input$targets
  months <- list(
    tropomi = input$column$month,
    acag = input$acag$month,
    no2 = targets$month[targets$source == "ground_no2"],
    bam = targets$month[targets$source == "ground_pm25"],
    # field months and daily-calibration months of the low-cost network
    moqa = c(targets$month[targets$source == "moqa"], input$calibration$month)
  )
  labels <- c(tropomi = "TROPOMI", acag = "SatPM₂.₅", no2 = "Reference NO₂",
    bam = "Reference PM₂.₅", moqa = "Low-cost network")
  colours <- c(tropomi = "steelblue", acag = "tomato", no2 = "#4a5568", bam = "#4a5568",
    moqa = "#2c5282")
  outside <- setdiff(unlist(months), calendar)
  if (length(outside)) stop("Months outside 2020-2025: ", toString(unique(outside)), call. = FALSE)
  tiles <- do.call(rbind, lapply(names(months), function(s)
    data.frame(source = s, x = match(sort(unique(months[[s]])), calendar),
      y = 6L - match(s, names(months)))))

  ggplot2::ggplot(tiles, ggplot2::aes(x, y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = source), colour = "white", linewidth = 0.28,
      width = 0.94, height = 0.64) +
    ggplot2::geom_vline(xintercept = seq(12.5, 60.5, 12), colour = "#C7D0D7", linewidth = 0.35) +
    ggplot2::scale_fill_manual(values = colours, guide = "none") +
    ggplot2::scale_x_continuous(breaks = seq(6.5, 66.5, 12), labels = 2020:2025,
      limits = c(0.5, 72.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = 1:5, labels = rev(labels),
      limits = c(0.5, 5.5), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_light(base_size = 22, base_family = paper_font) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8.2, colour = "black"),
      axis.text.x = ggplot2::element_text(size = 8.4, colour = "black"),
      axis.ticks = ggplot2::element_line(colour = "black"),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.4),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 3, 4, 5)
    )
}

.results_theme <- function(base_size) {
  ggplot2::theme_light(base_size = base_size, base_family = paper_font) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = base_size + 0.5, face = "bold", hjust = 0,
        margin = ggplot2::margin(b = 5)),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1, colour = "black"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#E5E5E5", linewidth = 0.25),
      plot.margin = ggplot2::margin(6, 8, 6, 6))
}

.map_theme <- function(base_size) {
  ggplot2::theme_void(base_size = base_size, base_family = paper_font) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size - 0.5, colour = paper_colours$text,
        hjust = 0.5, margin = ggplot2::margin(b = 1.5)),
      legend.position = "bottom",
      legend.title.position = "top",
      legend.title = ggplot2::element_text(size = base_size - 1.8, hjust = 0.5),
      legend.text = ggplot2::element_text(size = base_size - 2.2, colour = paper_colours$text),
      legend.frame = ggplot2::element_rect(colour = paper_colours$frame, linewidth = 0.3),
      legend.ticks = ggplot2::element_line(colour = paper_colours$frame),
      legend.margin = ggplot2::margin(t = 2, b = 0),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(1, 1.5, 1, 1.5))
}

# Writes the figure as a PDF with embedded fonts (cairo).
.save_pdf <- function(figure, path, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(path, width = width, height = height, family = paper_font)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(figure)
  invisible(path)
}

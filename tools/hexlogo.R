# Hex logo for the multilpa package.
#
# Self-contained: base graphics only, no hexSticker/ggplot2/showtext, matching
# the package's zero-dependency philosophy. Styled as a sibling to the `tna`
# package logo: light periwinkle fill, thick black border, bold monospace
# wordmark top-left, flat colour-filled markers with black outlines. The motif
# is three latent clusters -- one of squares, one of circles, one of triangles,
# each its own colour -- the mixture that latent profile analysis recovers.
#
# Rebuild with:
#   Rscript tools/hexlogo.R
# It writes man/figures/logo.png, which the README embeds at the top.

logo_path <- file.path("man", "figures", "logo2.png")
dir.create(dirname(logo_path), showWarnings = FALSE, recursive = TRUE)

## Palette -------------------------------------------------------------------
fill    <- "#CFE0F4"   # light periwinkle hex fill
border  <- "#000000"   # thick black border
sq_col  <- "#E8A13C"   # squares cluster
ci_col  <- "#3F6DA0"   # circles cluster
tr_col  <- "#3FA46F"   # triangles cluster

## Geometry helpers ----------------------------------------------------------
hex <- function(r = 1, rot = 90) {
  a <- (seq(0, 5) * 60 + rot) * pi / 180
  list(x = r * cos(a), y = r * sin(a))
}

## Interior half-width of a pointy-top hexagon of radius r at height y.
hex_halfwidth <- function(y, r = 1) {
  flat <- r * 0.5
  side <- r * cos(pi / 6)
  ifelse(abs(y) <= flat, side, pmax(0, side * (r - abs(y)) / (r - flat)))
}

## Boundary encompassing a cluster: a tinted disc with a black outline.
boundary <- function(cx, cy, r, col) {
  a <- seq(0, 2 * pi, length.out = 120)
  graphics::polygon(cx + r * cos(a), cy + r * sin(a),
                    col = grDevices::adjustcolor(col, 0.30),
                    border = border, lwd = 8)
}

## Draw a cluster of one shape around a centre.
cluster <- function(cx, cy, dx, dy, pch, col, cex) {
  graphics::points(cx + dx, cy + dy, pch = pch, col = border, bg = col,
                   lwd = 6, cex = cex)
}

## Canvas --------------------------------------------------------------------
grDevices::png(logo_path, width = 1600, height = 1600, bg = "transparent")
op <- graphics::par(mar = c(0, 0, 0, 0), xpd = NA, family = "mono")
plot.new()
plot.window(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), asp = 1)

## Hex body with a thick black border.
h <- hex(0.98)
graphics::polygon(h$x, h$y, col = fill, border = border, lwd = 42)

## Wordmark, top-left, bold monospace, sized as large as fits the hexagon.
wm_y <- 0.60
label <- "multilpa"
wm <- list(cex = 0.5, x = -0.4)
for (cex in seq(6, 0.5, by = -0.05)) {
  w <- graphics::strwidth(label, cex = cex, font = 2)
  ht <- graphics::strheight(label, cex = cex, font = 2)
  limit <- min(hex_halfwidth(wm_y + ht / 2, r = 0.86),
               hex_halfwidth(wm_y - ht / 2, r = 0.86)) * 0.97
  if (w <= 2 * limit) { wm <- list(cex = cex, x = -limit); break }
}
graphics::text(wm$x, wm_y, label, adj = c(0, 0.5),
               col = border, font = 2, cex = wm$cex)

## Three latent clusters: squares (top), circles (lower-left), triangles
## (lower-right), echoing a mixture of profiles. Each is the same dense blob
## of eight markers, distinguished only by shape and colour.
cex_m <- 8.6
blob_dx <- c(-0.085, 0.085, -0.17, 0.00, 0.17, -0.085, 0.085)
blob_dy <- c(0.17, 0.17, 0.00, 0.00, 0.00, -0.17, -0.17)
clusters <- list(
  list(cx = 0.00, cy = 0.24, pch = 22, col = sq_col),   # squares, top
  list(cx = -0.40, cy = -0.30, pch = 21, col = ci_col),  # circles, lower-left
  list(cx = 0.40, cy = -0.30, pch = 24, col = tr_col)    # triangles, lower-right
)
r_bound <- 0.30
for (cl in clusters) boundary(cl$cx, cl$cy, r_bound, cl$col)
for (cl in clusters) {
  cluster(cl$cx, cl$cy, blob_dx, blob_dy, pch = cl$pch, col = cl$col,
          cex = cex_m)
}

graphics::par(op)
grDevices::dev.off()
message("wrote ", logo_path)

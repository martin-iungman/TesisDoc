# Figura de zoom: panel completo (A) a la izquierda y una zona ampliada (B)
# a la derecha, con un recuadro en A y lineas de zoom hasta el panel de B.
# Requiere ggplot2, dplyr, tibble, grid y cowplot.

# Posicion (pulgadas, desde abajo-izquierda) del area del panel de un ggplot
# dibujado en un rectangulo de w x h pulgadas, y rangos de sus ejes.
panel_box <- function(p, w, h) {
  g <- ggplotGrob(p)
  pan <- g$layout[g$layout$name == "panel", ]
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grid::pushViewport(grid::viewport(width = grid::unit(w, "in"), height = grid::unit(h, "in")))
  suma_w <- function(idx) if (length(idx)) sum(grid::convertWidth(g$widths[idx], "in", valueOnly = TRUE)) else 0
  suma_h <- function(idx) if (length(idx)) sum(grid::convertHeight(g$heights[idx], "in", valueOnly = TRUE)) else 0
  n_w <- length(g$widths)
  n_h <- length(g$heights)
  pp <- ggplot_build(p)$layout$panel_params[[1]]
  list(
    x0 = suma_w(seq_len(pan$l - 1)), x1 = w - suma_w(if (pan$r < n_w) (pan$r + 1):n_w else integer(0)),
    y1 = h - suma_h(seq_len(pan$t - 1)), y0 = suma_h(if (pan$b < n_h) (pan$b + 1):n_h else integer(0)),
    xr = pp$x.range, yr = pp$y.range
  )
}

# Figura con el panel completo (A) a la izquierda y el zoom (B) a la
# derecha: recuadro en A alrededor de la zona ampliada (ylim_zoom en
# coordenadas de A; en x cubre lo que muestra B) y dos lineas desde sus
# esquinas derechas hasta las esquinas izquierdas del panel de B.
zoom_figure <- function(p_a, p_b, ylim_zoom, w_a = 9, w_b = 8, h = 8, h_b = 5.5, gap = 0.6, tag_size = 26) {
  b <- panel_box(p_b, w_b, h_b)
  a0 <- panel_box(p_a, w_a, h)
  xlim_zoom <- c(max(a0$xr[1], b$xr[1]), b$xr[2])
  p_a <- p_a + annotate("rect", xmin = xlim_zoom[1], xmax = xlim_zoom[2], ymin = ylim_zoom[1], ymax = ylim_zoom[2],
    fill = NA, colour = "grey25", linewidth = 0.7)
  a <- panel_box(p_a, w_a, h)
  w <- w_a + gap + w_b
  rect_x <- a$x0 + (xlim_zoom[2] - a$xr[1]) / diff(a$xr) * (a$x1 - a$x0)
  rect_y <- a$y0 + (ylim_zoom - a$yr[1]) / diff(a$yr) * (a$y1 - a$y0)
  lineas <- tibble::tibble(
    x = rect_x, y = rect_y,
    xend = w_a + gap + b$x0, yend = c(b$y0, b$y1)
  ) %>% dplyr::mutate(dplyr::across(c(x, xend), ~ .x / w), dplyr::across(c(y, yend), ~ .x / h))
  fig <- cowplot::ggdraw() +
    cowplot::draw_plot(p_a, x = 0, y = 0, width = w_a / w, height = 1) +
    cowplot::draw_plot(p_b, x = (w_a + gap) / w, y = 0, width = w_b / w, height = h_b / h) +
    geom_segment(data = lineas, aes(x = x, y = y, xend = xend, yend = yend), colour = "grey25", linewidth = 0.5, linetype = "dashed") +
    cowplot::draw_plot_label(c("A", "B"), x = c(0, (w_a + gap) / w), y = c(1, (h_b + 0.4) / h), size = tag_size) +
    theme(plot.background = element_rect(fill = "white", colour = NA))
  list(plot = fig, width = w, height = h)
}

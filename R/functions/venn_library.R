# Venn de dos conjuntos contenidos en la library: library como circulo
# exterior y los dos conjuntos superpuestos adentro. Geometria fija (no
# proporcional al tamano de los conjuntos). Estilo comun de
# representatividad_library_venn y filtros_consistencia_replicas.
#
# Uso: venn_library(a, b, library_ids, "Réplica 1", "Réplica 2") devuelve
# un ggplot; guardarlo con ggsave(..., width = 4.3, height = 4.3). text_scale
# agranda los rotulos si se usa a mayor tamano (ej. dentro de una figura
# compuesta).

venn_library <- function(a, b, library_ids, label_a, label_b, label_library = "Library", text_scale = 1) {
  fuera <- setdiff(union(a, b), library_ids)
  if (length(fuera) > 0) warning(length(fuera), " secuencias detectadas no estan en la library")
  n <- c(
    a_solo = length(setdiff(a, b)),
    ambas = length(intersect(a, b)),
    b_solo = length(setdiff(b, a)),
    library_solo = length(setdiff(library_ids, union(a, b)))
  )
  print(n)

  circle <- function(x0, y0, r, id) {
    t <- seq(0, 2 * pi, length.out = 361)
    tibble::tibble(x = x0 + r * cos(t), y = y0 + r * sin(t), id = id)
  }
  circles <- dplyr::bind_rows(
    circle(0, 0, 1, "library"),
    circle(-0.325, -0.05, 0.525, "a"),
    circle(0.325, -0.05, 0.525, "b")
  )
  fills <- c(library = "#EFF6FB", a = "#358AAA", b = "#AADAD4")
  alphas <- c(library = 1, a = 0.3, b = 0.5)

  labels <- tibble::tribble(
    ~x,     ~y,     ~label,                       ~size,
    0,      0.79,   label_library,                5,
    0,      0.66,   format(n[["library_solo"]]),  4.2,
    -0.47,  0.15,   label_a,                      5,
    -0.47,  -0.07,  format(n[["a_solo"]]),        4.2,
    0,      -0.07,  format(n[["ambas"]]),         4.2,
    0.59,   0.15,   label_b,                      5,
    0.59,   -0.07,  format(n[["b_solo"]]),        4.2
  )

  ggplot2::ggplot() +
    purrr::map(names(fills), function(i) {
      ggplot2::geom_polygon(
        data = circles[circles$id == i, ], ggplot2::aes(x, y),
        fill = fills[[i]], alpha = alphas[[i]], colour = NA
      )
    }) +
    # bordes en una capa aparte para que la transparencia no los atenue
    ggplot2::geom_path(data = circles, ggplot2::aes(x, y, group = id), colour = "grey15", linewidth = 0.4) +
    ggplot2::geom_text(data = labels, ggplot2::aes(x, y, label = label, size = size * text_scale), colour = "grey35", lineheight = 0.9) +
    ggplot2::scale_size_identity() +
    ggplot2::coord_equal(xlim = c(-1.02, 1.02), ylim = c(-1.02, 1.02), expand = FALSE) +
    ggplot2::theme_void()
}

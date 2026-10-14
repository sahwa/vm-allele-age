library(ggplot2)

zoo_cols <- c(purple = "#3E3A6E", purple_lt = "#A9A3D2",
              teal   = "#00A9A5", teal_lt   = "#7FD3D1",
              grey   = "#6D6E71", grey_lt   = "#C8C8C8",
              red    = "#C8102E", dark_red  = "#7B1416", salmon = "#E07B6F")

## order the palette so the first few categories are the most distinct
zoo_pal <- unname(zoo_cols[c("purple", "teal", "grey", "salmon", "purple_lt",
                             "teal_lt", "grey_lt", "dark_red", "red")])

theme_zoo <- function(base_size = 14, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line         = element_line(colour = "black", linewidth = 0.5),
      axis.ticks        = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = unit(5, "pt"),
      axis.text         = element_text(colour = "black", size = rel(0.9)),
      axis.title        = element_text(colour = "black"),
      panel.grid        = element_blank(),
      legend.background = element_rect(colour = "black", fill = "white", linewidth = 0.4),
      legend.key        = element_blank(),
      legend.title      = element_blank(),
      legend.margin     = margin(4, 6, 4, 6),
      strip.background  = element_blank(),
      strip.text        = element_text(size = rel(1), margin = margin(4, 4, 4, 4)),
      plot.title        = element_text(face = "bold", hjust = 0, size = rel(1.1),
                                       margin = margin(b = 6)),
      plot.tag          = element_text(face = "bold", size = rel(1.7)),
      plot.tag.position = c(0, 1)
    )
}

## make it the default for every plot in the session
use_zoo <- function(base_size = 14, base_family = "Helvetica") {
  theme_set(theme_zoo(base_size, base_family))
  options(
    ggplot2.discrete.colour   = function(...) scale_colour_manual(..., values = zoo_pal),
    ggplot2.discrete.fill     = function(...) scale_fill_manual(..., values = zoo_pal),
    ggplot2.continuous.colour = function(...)
      scale_colour_gradient(..., low = zoo_cols[["purple_lt"]], high = zoo_cols[["purple"]]),
    ggplot2.continuous.fill   = function(...)
      scale_fill_gradient(..., low = zoo_cols[["purple_lt"]], high = zoo_cols[["purple"]])
  )
}
use_zoo()

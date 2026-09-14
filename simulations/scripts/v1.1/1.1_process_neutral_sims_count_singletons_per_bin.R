f = glue::glue
res_list = purrr::map(1:100, function(x) {
    print(x)
    varinfo = f("rep{x}/1.1_variant_info.csv")
    dat = fread(varinfo)
    dat[, .(N_sing = sum(freq == 0.000005), N = .N, rep = x), bin]
}

res = rbindlist(res_list)
res[, bin := factor(bin, levels = gtools::mixedsort(unique(bin)))]

mean_counts <- res %>%
  group_by(bin) %>%
  summarise(
    mean_N_sing = mean(N_sing),
    label_y = max(N_sing) * 1.08,
    .groups = "drop"
  )

p_singletons_per_bin <- res %>%
  ggplot(aes(x = bin, y = N_sing)) +
  geom_jitter(
    width = 0.15,
    height = 0,
    alpha = 0.35,
    size = 1.5
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.6,
    fatten = 1.5
  ) +
  geom_label(
    data = mean_counts,
    aes(x = bin, y = label_y, 
        label = paste0("Mean = ", round(mean_N_sing, 1))),
    inherit.aes = FALSE,
    size = 3
  ) +
  labs(
    x = "Bin",
    y = "N_sing"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

FIGS <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v1.1"

ggsave(
  file.path(FIGS, "mean_singleton_per_bin.png"),
  p_singletons_per_bin
)
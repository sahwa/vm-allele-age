f = glue::glue

DATA="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v4.0"
OUT="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v4.0/dmean_test.png"
K_SELS <- sort(c(0.01, 0.005, 0.002))      # K = 0 can't be recovered from S_j

all <- map(K_SELS, function(K_SEL) {
    dat <- fread(file.path(DATA, f("sing_K{K_SEL}_seed8.TICKS.tsv")))
    setnames(dat, c("tick", "age", "S_j", "AC"))

    dat[, B2 := -S_j / K_SEL]
    dmean <- dat[, .(B2_mean = mean(B2), se = sd(B2) / sqrt(.N), n = .N), by = age][order(age)]
    dmean <- dmean[n >= 200]                # drop ages with too few mutations

    out <- rbind(
        dmean[, .(age, panel = "1. mean B2",                             value = B2_mean)],
        dmean[, .(age, panel = "2. log(mean B2): \nstraight = exponential", value = log(B2_mean))],
        dmean[, .(age, panel = "3. 1 / mean B2: \nstraight = hyperbolic",   value = 1 / B2_mean)])
    out[, K := K_SEL][]
})

total <- rbindlist(all)

p1 <- ggplot(total, aes(age, value)) +
  geom_line() +
  geom_smooth(data = total[!grepl("^1", panel)], method = "lm", se = FALSE,
              colour = "red", linewidth = 0.8, linetype = "dashed") +
  facet_grid(panel ~ K, scales = "free_y", labeller = labeller(K = label_both)) +
  theme_light() +
  theme(
    axis.text = element_text(size=18),
    axis.title = element_text(size=18),
    strip.text = element_text(size=12, colour='black'),
    strip.background = element_rect(fill='white')
    )
ggsave(OUT, p1, width = 12, height = 9, dpi = 150)

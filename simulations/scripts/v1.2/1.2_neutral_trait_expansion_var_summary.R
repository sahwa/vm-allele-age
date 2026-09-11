FIGS = "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v1.2"

dat = fread("1.2_variant_info.csv")
dat[, bin := factor(bin, levels = gtools::mixedsort(unique(bin)))]
dat[, MAF := pmin(freq, 1-freq)]

p1 = dat %>% ggplot(aes(x=MAF)) + geom_histogram(colour='black', fill='lightblue')
p2 = dat %>% ggplot(aes(x=age)) + geom_histogram(colour='black', fill='lightblue')
p_singleton_bin_hist = dat %>% 
    filter(freq == 0.000025) %>% 
    group_by(bin) %>% 
    summarise(N = n()) %>% 
    ggplot(aes(x=bin, y=N)) + 
    geom_col() +
    theme(axis.text.x = element_text(angle = 10))


ggsave(file.path(FIGS, "MAF_spectrum.png"), p1)
ggsave(file.path(FIGS, "age_spectrum.png"), p2)
ggsave(file.path(FIGS, "singleton_bin_hist.png"), p_singleton_bin_hist)




C <- fread("1.2_singleton_counts.csv")
C[, ind := str_c("ind", .I)]
p_mean_singleton_per_bin = melt(C, id.vars = "ind")[, .(mean(value)), variable] %>% 
    ggplot(aes(x=variable, y=V1)) + 
    geom_col() + 
    xlab("Bin (gens)") + 
    ylab("Mean number of singletons per individual") +
    theme(axis.text.x = element_text(angle = 10))
    
ggsave(file.path(FIGS, "mean_singleton_bin_per_ind_hist.png"), p_mean_singleton_per_bin)

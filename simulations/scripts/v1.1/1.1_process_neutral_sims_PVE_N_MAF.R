get_var_info = function(x) {

    var_info = fread(f("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1/rep{x}/1.1_variant_info.csv"))

    var_info[, MAF := fifelse(freq > 0.5, 1-freq, freq)]
    setorder(var_info, MAF)

    var_info[, pva := 2 * freq * (1-freq) * (beta^2)]
    var_info[, csum_beta := cumsum(pva)]
    var_info[, csum_beta_prop := csum_beta / max(csum_beta)]

    # Cumulative SNP count
    var_info[, N_SNP := 1:.N]
    var_info[, N_SNP_prop := N_SNP / .N]
    var_info[, REP := x]
    var_info
}

all_reps = rbindlist(lapply(1:10, get_var_info))

# Plot both on one graph with dual y-axes
p1 <- ggplot(all_reps, aes(x=MAF)) +
    geom_line(aes(y=csum_beta_prop, colour="V_A", group=REP), linewidth=1) +
    geom_line(aes(y=N_SNP_prop, colour="SNP count", group=REP), linewidth=1) +
    scale_y_continuous(name="Cumulative Proportion") +
    scale_colour_manual(values=c("V_A"="#1f77b4", "SNP count"="#ff7f0e")) +
    labs(x="Minor Allele Frequency", title="Cumulative V_A vs SNP count by MAF") +
    theme_minimal() +
    theme(legend.position="top")
	
ggsave("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1/rep1/1.1_variant_info.prop_var_explained_MAF.png", p1)
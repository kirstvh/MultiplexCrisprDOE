n_gRNA_perconstruct = 1

n_targets = 20
n_gRNA_pergene = 6
n_gRNA_total = n_targets * n_gRNA_pergene
ϵ_knockout_global = 0.8

Random.seed!(1)
ρ = 2; l = 50; u = ρ*l
m = (l+u)/2; 
sd = (u-l)/2; 
p_gRNA_library = gRNA_frequency_distribution(m, sd, l, u, n_gRNA_total; normalize=true)


Random.seed!(1)
p_active = 0.9; 
high_activity = 0.9; 
low_activity = 0.1; 
sd_act = 0.01
p_gRNA_act = gRNA_activity_distribution(p_active, high_activity, low_activity, sd_act, n_gRNA_total; visualize=false);
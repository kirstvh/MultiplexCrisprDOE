r = 2 # number of gRNAs per combinatorial gRNA/Cas construct

x = 20 # number of target genes
g = 6 # number of gRNAs per target gene
n_gRNA_total = x * g  # total number of gRNAs in the experiment
ϵ_KO = 0.8  # global knockout efficiency; fraction of mutations leading to effective gene knockout

Random.seed!(1)
ρ = 2; # ratio of the frequency of the most abundant gRNA and the frequency of the least abundant gRNA
l = 50; u = ρ*l
m = (l+u)/2; 
sd = (u-l)/2; 
p_gRNA_freq = gRNA_frequency_distribution(m, sd, l, u, n_gRNA_total; normalize=true)  # generate gRNA frequency distribution in construct library


Random.seed!(1)
f_act = 0.9  # fraction of all gRNAs that is active
ϵ_edit_act = 0.9;  # average genome editing efficiency of active gRNAs
ϵ_edit_inact = 0.1; # average genome editing efficiency of inactive gRNAs
sd_act = 0.01  # standard deviation 
p_gRNA_edit = gRNA_edit_distribution(f_act, ϵ_edit_act, ϵ_edit_inact, sd_act, n_gRNA_total; visualize=false); # generate genome editing distribution of the gRNAs
include("MultiplexCrisprDOE.jl");
using Random, Distributions
using LinearAlgebra
using Combinatorics
using BioCCP
using XLSX
using DataFrames

k = 2 # order of interaction to investigate
r = 2 # number of gRNAs per combinatorial gRNA/Cas construct

x = 20 # number of target genes
g = 6 # number of gRNAs per target gene
n_gRNA_total = x * g  # total number of gRNAs in the experiment

#### (1) Simulated Distribution

Random.seed!(1)
ρ = 2; # ratio of the frequency of the most abundant gRNA and the frequency of the least abundant gRNA
l = 50; u = Int(ρ*l)
m = Int((l+u)/2); 
sd = Int((u-l)/2); 
# p_gRNA_reads_simulated = gRNA_frequency_distribution(m, sd, l, u, n_gRNA_total; normalize = false, visualize=false)
# p_gRNA_reads_normalized_simulated = gRNA_frequency_distribution(m, sd, l, u, n_gRNA_total; normalize=true)
println("Simulate Distribution")
run(`julia main.jl gfd $m $sd $l $u $n_gRNA_total`)

#### (2) upload own data
#### Make sure that all gRNAs are ranked per target gene, and that for each target gene, there are g gRNAs
filename = "example_data.xlsx"
sheet = 1
data = DataFrame(XLSX.readtable(filename, 1)...)
p_gRNA_reads = data[!,"gRNA_read"]
p_gRNA_reads_normalized = p_gRNA_reads/sum(p_gRNA_reads)  # normalize

#### (1) Simulate editing efficiency
Random.seed!(1)
f_act = 0.9  # fraction of all gRNAs that is active
ϵ_edit_act = 0.95;  # average genome editing efficiency of active gRNAs
ϵ_edit_inact = 0.1; # average genome editing efficiency of inactive gRNAs
sd_act = 0.01  # standard deviation 
# p_gRNA_edit_simulated = gRNA_edit_distribution(f_act, ϵ_edit_act, ϵ_edit_inact, sd_act, n_gRNA_total; visualize=false); # generate genome editing distribution of the gRNAs
println("Simulate editing efficiency")
run(`julia main.jl ged $f_act $ϵ_edit_act $ϵ_edit_inact $sd_act $n_gRNA_total`)

#### (2) Upload own data
#### Make sure that all gRNAs are ranked per target gene, and that for each target gene, there are g gRNAs
using DataFrames, XLSX
filename = "example_data.xlsx"
sheet = 1
data = DataFrame(XLSX.readtable(filename, 1)...)
p_gRNA_edit = data[!,"gRNA_edit_efficiency"]
# p_gRNA_edit

ϵ_KO = 0.8  # global knockout efficiency; fraction of mutations leading to effective gene knockout

println("Knockouts per plant in the experiment:")
r = 2
run(`julia main.jl sim 4 $x $g $r $n_gRNA_total $filename $filename $ϵ_KO --i 10`)

println("Simulation-based approaches:")
for r in 1:3
    t = r
    println("mode: ", t)
    println("r: ", r)
    run(`julia main.jl sim $t $x $g $r $n_gRNA_total $filename $filename $ϵ_KO --i 10`)
end

println("BioCCP-based approaches:")
for (t, r) in zip(collect(1:9), repeat([1,2,3],3))
    println("mode: ", t)
    println("r: ", r)
    run(`julia main.jl ccp $t $x 0 $g $r $n_gRNA_total example_data.xlsx example_data.xlsx $ϵ_KO --s 5 --MN 6000`)
end

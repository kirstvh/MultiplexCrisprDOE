module MultiplexCrisprDOE
 
using Base: Integer
using Plots
using Distributions
using Combinatorics
using BioCCP
export gRNA_frequency_distribution, gRNA_activity_distribution

"""
    gRNA_library_distribution(m, sd, l, u, gRNA_total; normalize = true, visualize=false)

Generates gRNA relative frequency distribution in combinatorial gRNA/Cas9 construct library.
"""
function gRNA_frequency_distribution(m, sd, l, u, gRNA_total; normalize = true, visualize=false)
    μ = log(m/sqrt(1+sd^2/m^2)) # https://discourse.julialang.org/t/lognormal-distribution-how-to-set-mu-and-sigma/7101/6
    σ = sqrt(log(1+sd^2/m^2))
    d_gRNAlibrary = truncated(LogNormal(μ,σ), l, u)
    gRNA_abundances = collect(rand(d_gRNAlibrary, gRNA_total))
     if visualize
        return histogram(gRNA_abundances, label="", xlabel="Number of reads per gene", ylabel="absolute frequency", title="Read distribution")
    else
        if normalize
            gRNA_abundances /= sum(gRNA_abundances)
        end
        return gRNA_abundances
    end
end

"""
    gRNA_activity_distribution(p_high_activity, μ_high_activity, μ_low_activity, σ_activity, n_gRNAs; visualize=false)   

Generates bimodal distribution of genome editing efficiencies for the gRNAs.
"""
function gRNA_activity_distribution(p_high_activity, μ_high_activity, μ_low_activity, σ_activity, n_gRNAs; visualize=false)   
    d_activity = Binomial(1, p_high_activity)
    d_highactivity = truncated(Normal(μ_high_activity, σ_activity), 0.01, 1)
    d_lowactivity = truncated(Normal(μ_low_activity, σ_activity), 0.01, 1)
    p_gRNA_act = zeros(n_gRNAs) # initialize
    for i in 1:n_gRNAs
        if rand(d_activity, 1) == [1]
            p_gRNA_act[i] = rand(d_highactivity, 1)[1]
        else
            p_gRNA_act[i] = rand(d_lowactivity, 1)[1]
        end
    end
    if visualize
        return histogram(p_gRNA_act, label="", xlabel="gRNA activity", ylabel="absolute frequency", title="gRNA activity distribution")

    else
        return p_gRNA_act
    end
end

"""
    simulate_RES1(n_targets, 
                n_gRNA_pergene, 
                n_gRNA_perconstruct, 
                n_gRNA_total, 
                p_gRNA_library, 
                p_gRNA_act, 
                ϵ_knockout_global; 
                iter = 500)

Simulation-based approach for calculating RES1 of a CRISPR/Cas experiment.
"""
function simulate_RES1(n_targets, 
                        n_gRNA_pergene, 
                        n_gRNA_perconstruct, 
                        n_gRNA_total, 
                        p_gRNA_library, 
                        p_gRNA_act, 
                        ϵ_knockout_global; 
                        iter=500)

    @assert n_targets * n_gRNA_pergene == n_gRNA_total
    
    T_vec = [] #stores number of plants required for each experiment
        for i in 1:iter       
            genes_vec = [] # Initialize matrix to count pairwise interactions
            T = 0
            while genes_vec != collect(1:n_targets) # check if all pairwise combinations are present
             
                T += 1 # count how many plants must be sampled to fill pairwise interaction matrix
                
                # sample combinatorial gRNA/Cas9 construct
                gRNA_indices_construct = findall((rand(Multinomial(n_gRNA_perconstruct, p_gRNA_library))) .!= 0)
                
                # execute mutations
                gRNA_indices_mutations = [gRNA for gRNA in gRNA_indices_construct if rand(Binomial(1, p_gRNA_act[gRNA])) == 1]
            
                # effective gene knockout (loss-of-function) ?
                gRNA_indices_KO = [gRNA for gRNA in gRNA_indices_mutations if rand(Binomial(1, ϵ_knockout_global)) == 1]
            
                # which genes are knocked out?
                genes_indices_KO = Int.(ceil.(gRNA_indices_KO / n_gRNA_pergene))
                 append!(genes_vec, genes_indices_KO)
                genes_vec = Int.(sort(unique(genes_vec)))
            end
            push!(T_vec, T)   
              
        end
        E = mean(T_vec); sd = std(T_vec)
    return E, sd
end

"""
    BioCCP_RES1(n_targets, 
                n_gRNA_pergene, 
                n_gRNA_perconstruct, 
                n_gRNA_total, 
                p_gRNA_library, 
                p_gRNA_act, 
                ϵ_knockout_global)

Employs the BioCCP.jl package to calculate the RES1 of a CRISPR/Cas experiment.
"""
function BioCCP_RES1(n_targets, 
                    n_gRNA_pergene, 
                    n_gRNA_perconstruct, 
                    n_gRNA_total, 
                    p_gRNA_library, 
                    p_gRNA_act, 
                    ϵ_knockout_global)

    p_gRNAs = p_gRNA_library .* p_gRNA_act * ϵ_knockout_global
    p_genes = [sum(p_gRNAs[i:i+n_gRNA_pergene-1]) for i in 1:n_gRNA_pergene:n_gRNA_total]
    return expectation_minsamplesize(n_targets; p=p_genes, r=n_gRNA_perconstruct, normalize=false), 
            std_minsamplesize(n_targets; p=p_genes, r=n_gRNA_perconstruct, normalize=false)
end

"""
    simulate_RES2(n_targets, 
                    n_gRNA_pergene, 
                    n_gRNA_perconstruct, 
                    n_gRNA_total, 
                    p_gRNA_library, 
                    p_gRNA_act, 
                    ϵ_knockout_global; 
                    iter=500)

Simulation-based approach for calculating RES2 of multiplex CRISPR/Cas experiment.
"""
function simulate_RES2(n_targets, 
                        n_gRNA_pergene, 
                        n_gRNA_perconstruct, 
                        n_gRNA_total, 
                        p_gRNA_library, 
                        p_gRNA_act, 
                        ϵ_knockout_global; 
                        iter=500)

    @assert n_targets * n_gRNA_pergene == n_gRNA_total
    
    T_vec = [] #stores number of plants required for each experiment
        for i in 1:iter     
            if i % 100 == 0
                println("Iteration: $i ... Calculating minimum number of plants ... \n")
            end   
            X_interactions_count = zeros(n_targets, n_targets) # Initialize matrix to count pairwise interactions
            T = 0
            while X_interactions_count != ones(n_targets, n_targets) # check if all pairwise combinations are present
                T += 1 # count how many plants must be sampled to fill pairwise interaction matrix
                
                # sample combinatorial gRNA/Cas9 construct
                gRNA_indices_construct = findall((rand(Multinomial(n_gRNA_perconstruct, p_gRNA_library))) .!= 0)
                
                # execute mutations
                gRNA_indices_mutations = [gRNA for gRNA in gRNA_indices_construct if rand(Binomial(1, p_gRNA_act[gRNA])) == 1]
            
                # effective gene knockout (loss-of-function) ?
                gRNA_indices_KO = [gRNA for gRNA in gRNA_indices_mutations if rand(Binomial(1, ϵ_knockout_global)) == 1]
            
                # which genes are knocked out?
                genes_indices_KO = Int.(ceil.(gRNA_indices_KO / n_gRNA_pergene))
            
                # which pairwise combinations are present?
                interactions = collect(combinations(genes_indices_KO, 2))
                
                # Store represented combinations in matrix
                for interaction in interactions
                    j = interaction[1]; k = interaction[2]
                    X_interactions_count[j,k] = 1; X_interactions_count[k,j] = 1; X_interactions_count[j,j] = 1; X_interactions_count[k,k] = 1          
                end  
            end
            push!(T_vec, T)   
              
        end
        E = mean(T_vec); sd = std(T_vec)
    return E, sd
end

"""
    BioCCP_RES2(n_targets, 
                n_gRNA_pergene, 
                n_gRNA_perconstruct, 
                n_gRNA_total, 
                p_gRNA_library, 
                p_gRNA_act, ϵ_knockout_global)

Employs the BioCCP.jl package to calculate the RES2 of a multiplex CRISPR/Cas experiment.
"""
function BioCCP_minsamplesize_pairwise(n_targets, 
                                         n_gRNA_pergene, 
                                         n_gRNA_perconstruct, 
                                         n_gRNA_total, 
                                         p_gRNA_library, 
                                         p_gRNA_act, ϵ_knockout_global)
    
    # how many pairwise combinations of gRNAs
    ind_combinations_gRNA = collect(combinations(1:n_gRNA_total, 2))
    n_combinations_gRNA = length(ind_combinations_gRNA)
    
    # calculate probability and activity of gRNA combinations
    p_combinations_gRNA_library = zeros(n_combinations_gRNA)
    p_combinations_gRNA_act = zeros(n_combinations_gRNA)
    for i in 1:n_combinations_gRNA
        p_combinations_gRNA_library[i] = p_gRNA_library[ind_combinations_gRNA[i][1]] * p_gRNA_library[ind_combinations_gRNA[i][2]]
        p_combinations_gRNA_act[i] = p_gRNA_act[ind_combinations_gRNA[i][1]] * p_gRNA_act[ind_combinations_gRNA[i][2]]
    end
    
    # normalize probability gRNA combinations
    p_combinations_gRNA_library /= sum(p_combinations_gRNA_library)

    # select pairwise gRNA combinations of which each component codes for different gene (goal is to study combinations of knockouts in different genes)
    p_combinations_gRNA_library_interest = []
    p_combinations_gRNA_act_interest = []
    ind_combinations_gRNA_interest = []
    for i in 1:n_combinations_gRNA
        if ceil(ind_combinations_gRNA[i][1]/n_gRNA_pergene) != ceil(ind_combinations_gRNA[i][2]/n_gRNA_pergene)
            push!(p_combinations_gRNA_library_interest, p_combinations_gRNA_library[i])
            push!(p_combinations_gRNA_act_interest, p_combinations_gRNA_act[i])
            push!(ind_combinations_gRNA_interest, ind_combinations_gRNA[i])
        end
    end
        
    n_combinations_gRNA_interest = length(p_combinations_gRNA_library_interest)
    p_combinations_gRNA = p_combinations_gRNA_library_interest .* p_combinations_gRNA_act_interest * ϵ_knockout_global^2

    #### INTEGREREN PER GENCOMBINATIE
    p_genes_matrix = zeros(n_targets, n_targets)
    for i in 1:n_combinations_gRNA_interest
        gene1 = Int(ceil(ind_combinations_gRNA_interest[i][1]/n_gRNA_pergene))
        gene2 = Int(ceil(ind_combinations_gRNA_interest[i][2]/n_gRNA_pergene))
        p_genes_matrix[gene1, gene2] += p_combinations_gRNA[i]
    end

    p_genes = collect([p_genes_matrix[i, j] for j in 2:size(p_genes_matrix, 1) for i in 1:j-1])  
    n_combinations_genes = length(p_genes)
    combinations_pp = length(collect(combinations(1:n_gRNA_perconstruct, 2)))
    
    return expectation_minsamplesize(n_combinations_genes; p=p_genes, r=combinations_pp, normalize=false), std_minsamplesize(n_combinations_genes; p=p_genes, r=combinations_pp, normalize=false)
end

"""
    BioCCP_Ps2(n_targets, 
                sample_size,
                n_gRNA_pergene, 
                n_gRNA_perconstruct, 
                n_gRNA_total, 
                p_gRNA_library, 
                p_gRNA_act, 
                ϵ_knockout_global)

Employs the BioCCP.jl package to calculate the probability of full coverage of all pairwise combinations of gene knockouts 
in the combinatorial gene knockout library of a multiplex CRISPR/Cas experiment 
with respect to a specified experimental scale (number of plants analyzed; ES).
"""
function BioCCP_Ps_2(n_targets, sample_size,
                                         n_gRNA_pergene, 
                                         n_gRNA_perconstruct, 
                                         n_gRNA_total, 
                                         p_gRNA_library, 
                                         p_gRNA_act, ϵ_knockout_global)
    
    # how many pairwise combinations of gRNAs
    ind_combinations_gRNA = collect(combinations(1:n_gRNA_total, 2))
    n_combinations_gRNA = length(ind_combinations_gRNA)
    
    # calculate probability and activity of gRNA combinations
    p_combinations_gRNA_library = zeros(n_combinations_gRNA)
    p_combinations_gRNA_act = zeros(n_combinations_gRNA)
    for i in 1:n_combinations_gRNA
        p_combinations_gRNA_library[i] = p_gRNA_library[ind_combinations_gRNA[i][1]] * p_gRNA_library[ind_combinations_gRNA[i][2]]
        p_combinations_gRNA_act[i] = p_gRNA_act[ind_combinations_gRNA[i][1]] * p_gRNA_act[ind_combinations_gRNA[i][2]]
    end
    
    # normalize probability gRNA combinations
    p_combinations_gRNA_library /= sum(p_combinations_gRNA_library)

    # select pairwise gRNA combinations of which each component codes for different gene (goal is to study combinations of knockouts in different genes)
    p_combinations_gRNA_library_interest = []
    p_combinations_gRNA_act_interest = []
    ind_combinations_gRNA_interest = []
    for i in 1:n_combinations_gRNA
        if ceil(ind_combinations_gRNA[i][1]/n_gRNA_pergene) != ceil(ind_combinations_gRNA[i][2]/n_gRNA_pergene)
            push!(p_combinations_gRNA_library_interest, p_combinations_gRNA_library[i])
            push!(p_combinations_gRNA_act_interest, p_combinations_gRNA_act[i])
            push!(ind_combinations_gRNA_interest, ind_combinations_gRNA[i])
        end
    end
        
    n_combinations_gRNA_interest = length(p_combinations_gRNA_library_interest)
    p_combinations_gRNA = p_combinations_gRNA_library_interest .* p_combinations_gRNA_act_interest * ϵ_knockout_global^2

    #### INTEGREREN PER GENCOMBINATIE
    p_genes_matrix = zeros(n_targets, n_targets)
    for i in 1:n_combinations_gRNA_interest
        gene1 = Int(ceil(ind_combinations_gRNA_interest[i][1]/n_gRNA_pergene))
        gene2 = Int(ceil(ind_combinations_gRNA_interest[i][2]/n_gRNA_pergene))
        p_genes_matrix[gene1, gene2] += p_combinations_gRNA[i]
    end

    p_genes = collect([p_genes_matrix[i, j] for j in 2:size(p_genes_matrix, 1) for i in 1:j-1])  
    n_combinations_genes = length(p_genes)
    combinations_pp = length(collect(combinations(1:n_gRNA_perconstruct, 2)))
    
    return success_probability(n_combinations_genes, sample_size; p=p_genes, r=combinations_pp, normalize=false)
end
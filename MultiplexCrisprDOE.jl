
"""
    gRNA_library_distribution(m, sd, l, u, gRNA_total; normalize = true, visualize=false)

Generates gRNA relative frequency distribution in combinatorial gRNA/Cas9 construct library.
"""
function gRNA_frequency_distribution(m, sd, l, u, gRNA_total; normalize = true, visualize=false)
    d_gRNAlibrary = truncated(Normal(m, sd), l, u)
    gRNA_abundances = collect(rand(d_gRNAlibrary, gRNA_total))
     if visualize
        return histogram(gRNA_abundances, label="", 
            xlabel="Number of reads per gene", 
            ylabel="absolute frequency", title="Read distribution")
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
    simulate_Nₓ₁(n_targets, 
                n_gRNA_pergene, 
                n_gRNA_perconstruct, 
                n_gRNA_total, 
                p_gRNA_library, 
                p_gRNA_act, 
                ϵ_knockout_global; 
                iter = 500)

Simulation-based approach for calculating Nₓ₁ of a CRISPR/Cas experiment.
"""
function simulate_Nₓ₁(n_targets, 
                                         n_gRNA_pergene, 
                                         n_gRNA_perconstruct, 
                                         n_gRNA_total, 
                                         p_gRNA_library, 
                                         p_gRNA_act, ϵ_knockout_global; iter=500)
    """ 
    INPUT
  
    
    OUTPUT
    E: expected minimum number of plants to gRNA_read_distributionsee each pairwise combination at least once
    sd: standard deviation on the minimum number of plants
    """
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
    BioCCP_Nₓ₁(n_targets, 
                n_gRNA_pergene, 
                n_gRNA_perconstruct, 
                n_gRNA_total, 
                p_gRNA_library, 
                p_gRNA_act, 
                ϵ_knockout_global)

Employs the BioCCP.jl package to calculate the Nₓ₁ of a CRISPR/Cas experiment.
"""
function BioCCP_Nₓ₁(n_targets, 
                                         n_gRNA_pergene, 
                                         n_gRNA_perconstruct, 
                                         n_gRNA_total, 
                                         p_gRNA_library, 
                                         p_gRNA_act, ϵ_knockout_global)
    
    p_gRNAs = p_gRNA_library .* p_gRNA_act * ϵ_knockout_global
    p_genes = [sum(p_gRNAs[i:i+n_gRNA_pergene-1]) for i in 1:n_gRNA_pergene:n_gRNA_total]
    return expectation_minsamplesize(n_targets; p=p_genes, r=n_gRNA_perconstruct, normalize=false), 
    std_minsamplesize(n_targets; p=p_genes, r=n_gRNA_perconstruct, normalize=false)
end


"""
    simulate_Nₓ₂(n_targets, 
                    n_gRNA_pergene, 
                    n_gRNA_perconstruct, 
                    n_gRNA_total, 
                    p_gRNA_library, 
                    p_gRNA_act, 
                    ϵ_knockout_global; 
                    iter=500)

Simulation-based approach for calculating Nₓ₂ of multiplex CRISPR/Cas experiment.
"""
function simulate_Nₓ₂(n_targets, 
                                         n_gRNA_pergene, 
                                         n_gRNA_perconstruct, 
                                         n_gRNA_total, 
                                         p_gRNA_library, 
                                         p_gRNA_act, ϵ_knockout_global; iter=500)
    """ 
    INPUT
  
    
    OUTPUT
    E: expected minimum number of plants to gRNA_read_distributionsee each pairwise combination at least once
    sd: standard deviation on the minimum number of plants
    """
    @assert n_targets * n_gRNA_pergene == n_gRNA_total
#     @assert sum(p_gRNA_library) == 1
    
    T_vec = [] #stores number of plants required for each experiment
        for i in 1:iter     
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
    BioCCP_Nₓ₂(n_targets, 
                n_gRNA_pergene, 
                n_gRNA_perconstruct, 
                n_gRNA_total, 
                p_gRNA_library, 
                p_gRNA_act, ϵ_knockout_global)

Employs the BioCCP.jl package to calculate the Nₓ₂ of a multiplex CRISPR/Cas experiment.
"""
function BioCCP_Nₓ₂(n_targets, 
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

function simulate_Nₓ₂_countKOs(n_targets, 
                                         n_gRNA_pergene, 
                                         n_gRNA_perconstruct, 
                                         n_gRNA_total, 
                                         p_gRNA_library, 
                                         p_gRNA_act, ϵ_knockout_global; iter=500)
     
    @assert n_targets * n_gRNA_pergene == n_gRNA_total
    
            n_KOs = []
       
            for j in 1:100000
                               
                # sample combinatorial gRNA/Cas9 construct
                gRNA_indices_construct = findall((rand(Multinomial(n_gRNA_perconstruct, p_gRNA_library))) .!= 0)
                
                # execute mutations
                gRNA_indices_mutations = [gRNA for gRNA in gRNA_indices_construct if rand(Binomial(1, p_gRNA_act[gRNA])) == 1]
            
                # effective gene knockout (loss-of-function) ?
                gRNA_indices_KO = [gRNA for gRNA in gRNA_indices_mutations if rand(Binomial(1, ϵ_knockout_global)) == 1]
            
                # which genes are knocked out?
                genes_indices_KO = Int.(ceil.(gRNA_indices_KO / n_gRNA_pergene))
            
               push!(n_KOs, length(unique((genes_indices_KO))))
            end  
 
    return n_KOs
end

"""
    simulate_Nₓ₃(n_targets, 
                    n_gRNA_pergene, 
                    n_gRNA_perconstruct, 
                    n_gRNA_total, 
                    p_gRNA_library, 
                    p_gRNA_act, 
                    ϵ_knockout_global; 
                    iter=500)

Simulation-based approach for calculating Nₓ₃ of multiplex CRISPR/Cas experiment.
"""
function simulate_Nₓ₃(n_targets, 
                                         n_gRNA_pergene, 
                                         n_gRNA_perconstruct, 
                                         n_gRNA_total, 
                                         p_gRNA_library, 
                                         p_gRNA_act, ϵ_knockout_global; iter=500)
    """ 
    INPUT
  
    
    OUTPUT
    E: expected minimum number of plants to gRNA_read_distributionsee each pairwise combination at least once
    sd: standard deviation on the minimum number of plants
    """
    @assert n_targets * n_gRNA_pergene == n_gRNA_total
#     @assert sum(p_gRNA_library) == 1
    
    T_vec = [] #stores number of plants required for each experiment
        for i in 1:iter       
            X_interactions_count = zeros(n_targets, n_targets, n_targets) # Initialize matrix to count triple interactions
            T = 0
            while X_interactions_count != ones(n_targets, n_targets, n_targets) # check if all triple combinations are present
                T += 1 # count how many plants must be sampled to fill triple interaction matrix
                
                # sample combinatorial gRNA/Cas9 construct
                gRNA_indices_construct = findall((rand(Multinomial(n_gRNA_perconstruct, p_gRNA_library))) .!= 0)
                
                # execute mutations
                gRNA_indices_mutations = [gRNA for gRNA in gRNA_indices_construct if rand(Binomial(1, p_gRNA_act[gRNA])) == 1]
            
                # effective gene knockout (loss-of-function) ?
                gRNA_indices_KO = [gRNA for gRNA in gRNA_indices_mutations if rand(Binomial(1, ϵ_knockout_global)) == 1]
            
                # which genes are knocked out?
                genes_indices_KO = Int.(ceil.(gRNA_indices_KO / n_gRNA_pergene))
            
                # which triple combinations are present?
                interactions = collect(combinations(genes_indices_KO, 3))
                
                # Store represented combinations in matrix
                for interaction in interactions
                    j = interaction[1]
                    k = interaction[2]
                    l = interaction[3]
                    X_interactions_count[j,k,l] = 1
                    X_interactions_count[k,j,l] = 1
                    X_interactions_count[l,j,k] = 1
                    X_interactions_count[l,k,j] = 1
                    X_interactions_count[j,l,k] = 1
                    X_interactions_count[k,l,j] = 1
        
                    X_interactions_count[:,l,l] .= 1
                    X_interactions_count[:,k,k] .= 1
                    X_interactions_count[:,j,j] .= 1
                    X_interactions_count[l,:,l] .= 1
                    X_interactions_count[k,:,k] .= 1
                    X_interactions_count[j,:,j] .= 1
        
                    X_interactions_count[j,j,:] .= 1
                    X_interactions_count[k,k,:] .= 1
                    X_interactions_count[l,l,:] .= 1
                end  
            end
            push!(T_vec, T)   
              
        end
        E = mean(T_vec); sd = std(T_vec)
    return E, sd
end
   

"""
    BioCCP_Nₓ₃(n_targets, 
            n_gRNA_pergene, 
            n_gRNA_perconstruct, 
            n_gRNA_total, 
            p_gRNA_library, 
            p_gRNA_act, 
            ϵ_knockout_global)

Employs the BioCCP.jl package to calculate the RES3 of a multiplex CRISPR/Cas experiment.
"""
function BioCCP_Nₓ₃(n_targets, 
                                         n_gRNA_pergene, 
                                         n_gRNA_perconstruct, 
                                         n_gRNA_total, 
                                         p_gRNA_library, 
                                         p_gRNA_act, ϵ_knockout_global)
    
    # how many pairwise combinations of gRNAs
    ind_combinations_gRNA = collect(combinations(1:n_gRNA_total, 3))
    n_combinations_gRNA = length(ind_combinations_gRNA)
    
    # calculate probability and activity of triple gRNA combinations
    p_combinations_gRNA_library = zeros(n_combinations_gRNA)
    p_combinations_gRNA_act = zeros(n_combinations_gRNA)
    for i in 1:n_combinations_gRNA
        p_combinations_gRNA_library[i] = p_gRNA_library[ind_combinations_gRNA[i][1]] * p_gRNA_library[ind_combinations_gRNA[i][2]] * p_gRNA_library[ind_combinations_gRNA[i][3]]
        p_combinations_gRNA_act[i] = p_gRNA_act[ind_combinations_gRNA[i][1]] * p_gRNA_act[ind_combinations_gRNA[i][2]] * p_gRNA_act[ind_combinations_gRNA[i][3]]
    end
    
    # normalize probability gRNA combinations
    p_combinations_gRNA_library /= sum(p_combinations_gRNA_library)

    # select triple gRNA combinations of which each component codes for different gene (goal is to study combinations of knockouts in different genes)
    p_combinations_gRNA_library_interest = []
    p_combinations_gRNA_act_interest = []
    ind_combinations_gRNA_interest = []
    for i in 1:n_combinations_gRNA
        if ceil(ind_combinations_gRNA[i][1]/n_gRNA_pergene) != ceil(ind_combinations_gRNA[i][2]/n_gRNA_pergene) && ceil(ind_combinations_gRNA[i][1]/n_gRNA_pergene) != ceil(ind_combinations_gRNA[i][3]/n_gRNA_pergene) && ceil(ind_combinations_gRNA[i][3]/n_gRNA_pergene) != ceil(ind_combinations_gRNA[i][2]/n_gRNA_pergene)
            push!(p_combinations_gRNA_library_interest, p_combinations_gRNA_library[i])
            push!(p_combinations_gRNA_act_interest, p_combinations_gRNA_act[i])
            push!(ind_combinations_gRNA_interest, ind_combinations_gRNA[i])
        end
    end
        
    n_combinations_gRNA_interest = length(p_combinations_gRNA_library_interest)
    p_combinations_gRNA = p_combinations_gRNA_library_interest .* p_combinations_gRNA_act_interest * ϵ_knockout_global^3

    #### INTEGREREN PER GENCOMBINATIE
    p_genes_matrix = zeros(n_targets, n_targets, n_targets)
    for i in 1:n_combinations_gRNA_interest
        gene1 = Int(ceil(ind_combinations_gRNA_interest[i][1]/n_gRNA_pergene))
        gene2 = Int(ceil(ind_combinations_gRNA_interest[i][2]/n_gRNA_pergene))
        gene3 = Int(ceil(ind_combinations_gRNA_interest[i][3]/n_gRNA_pergene))
        p_genes_matrix[gene1, gene2, gene3] += p_combinations_gRNA[i]
    end
    
    combinations_genes = collect(combinations(1:n_targets, 3))
    p_genes = []
        for combination in combinations_genes
            push!(p_genes, p_genes_matrix[combination[1], combination[2], combination[3]])
        end
        
    n_combinations_genes = length(p_genes)
    combinations_pp = length(collect(combinations(1:n_gRNA_perconstruct, 3)))
    
    return expectation_minsamplesize(n_combinations_genes; p=p_genes, r=combinations_pp, normalize=false), std_minsamplesize(n_combinations_genes; p=p_genes, r=combinations_pp, normalize=false)
end

"""
    BioCCP_Pₓ₁(n_targets, 
                sample_size,
                n_gRNA_pergene, 
                n_gRNA_perconstruct, 
                n_gRNA_total, 
                p_gRNA_library, 
                p_gRNA_act, 
                ϵ_knockout_global)
Employs the BioCCP.jl package to calculate the probability of full coverage of all individual gene knockouts 
in the gene knockout library of a CRISPR/Cas experiment 
with respect to a specified plant library size (number of plants analyzed).
"""
function BioCCP_Pₓ₁(n_targets, sample_size,
                                         n_gRNA_pergene, 
                                         n_gRNA_perconstruct, 
                                         n_gRNA_total, 
                                         p_gRNA_library, 
                                         p_gRNA_act, ϵ_knockout_global)
    
    p_gRNAs = p_gRNA_library .* p_gRNA_act * ϵ_knockout_global
    p_genes = [sum(p_gRNAs[i:i+n_gRNA_pergene-1]) for i in 1:n_gRNA_pergene:n_gRNA_total]
    return success_probability(n_targets, sample_size; p=p_genes, r=n_gRNA_perconstruct, normalize=false) 
end

"""
BioCCP_γₓ₁(n_targets, 
            sample_size,
            n_gRNA_pergene, 
            n_gRNA_perconstruct, 
            n_gRNA_total, 
            p_gRNA_library, 
            p_gRNA_act, 
            ϵ_knockout_global)

Employs the BioCCP.jl package to calculate the expected represented fraction of the design space, 
here consisting of all individual gene knockouts in a gene knockout library of a CRISPR/Cas experiment, 
with respect to a specified experimental scale (number of plants analyzed; ES).
"""
function BioCCP_γₓ₁(n_targets, sample_size,
                                         n_gRNA_pergene, 
                                         n_gRNA_perconstruct, 
                                         n_gRNA_total, 
                                         p_gRNA_library, 
                                         p_gRNA_act, ϵ_knockout_global)
    
    p_gRNAs = p_gRNA_library .* p_gRNA_act * ϵ_knockout_global
    p_genes = [sum(p_gRNAs[i:i+n_gRNA_pergene-1]) for i in 1:n_gRNA_pergene:n_gRNA_total]
    return expectation_fraction_collected(n_targets, sample_size; p=p_genes, r=n_gRNA_perconstruct, normalize=false) 
end


"""
    BioCCP_Pₓ₂(n_targets, 
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
function BioCCP_Pₓ₂(n_targets, sample_size,
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

"""
    BioCCP_γₓ₂(n_targets, 
                sample_size,
                n_gRNA_pergene, 
                n_gRNA_perconstruct, 
                n_gRNA_total, 
                p_gRNA_library, 
                p_gRNA_act, 
                ϵ_knockout_global)

Employs the BioCCP.jl package to calculate the expected represented fraction of the design space, here consisting of all pairwise combinations of gene knockouts 
in the combinatorial gene knockout library of a multiplex CRISPR/Cas experiment, 
with respect to a specified plant library size (number of plants analyzed).
"""
function BioCCP_γₓ₂(n_targets, sample_size,
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
    
    return expectation_fraction_collected(n_combinations_genes, sample_size; p=p_genes, r=combinations_pp, normalize=false)
end
      
"""
    gRNA_frequency_distribution(m, sd, l, u, n_gRNA_total; normalize = true, visualize=false)

Generates frequency distribution of the gRNAs in the combinatorial gRNA/Cas9 construct library.

m: the average abundance of the gRNAs (in terms of absolute or relative frequency)
sd: the standard deviation on the gRNA abundances (in terms of absolute or relative frequency)
l: minimal gRNA abundance (in terms of absolute or relative frequency)
u: maximal gRNA abundance (in terms of absolute or relative frequency)
n_gRNA_total: the total number of gRNAs in the experiment
normalize: if set to "true", the gRNA abundances (absolute frequencies) are converted into relative frequencies
visualize: if set to "true", a histogram of all gRNA abundances is plotted

"""
function gRNA_frequency_distribution(m, sd, l, u, n_gRNA_total; normalize = true, visualize=false)
    d_gRNA_freq = truncated(Normal(m, sd), l, u)  # gRNA frequency distribution
    p_gRNA_freq = collect(rand(d_gRNA_freq, n_gRNA_total))  # sample gRNA frequencies from distribution
    if normalize # convert into relative frequencies
        p_gRNA_freq /= sum(p_gRNA_freq)
    end
    if visualize
        return histogram(p_gRNA_freq, label="", 
            xlabel="Number of reads per gRNA", 
            ylabel="absolute frequency", title="Read distribution")
    else
        return p_gRNA_freq
    end
end

"""
    gRNA_edit_distribution(f_act, ϵ_edit_act, ϵ_edit_inact, sd_act, n_gRNA_total; visualize=false)   

Generates vector with genome editing efficiencies for all the gRNAs in the experiment. 

f_act: fraction of all gRNAs that is active
ϵ_edit_act: Average genome editing efficiency for active gRNAs - mean of the genome editing efficiency distribution for active gRNAs
ϵ_edit_inact: Average genome editing efficiency for inactive gRNAs - mean of the genome editing efficiency distribution for inactive gRNAs
sd_act: standard deviation of the genome editing efficiency distributions for active and inactive gRNAs
n_gRNA_total: the total number of gRNAs in the experiment
visualize: if set to "true", a histogram of all genome editing efficiency is plotted
"""
function gRNA_edit_distribution(f_act, ϵ_edit_act, ϵ_edit_inact, sd_act, n_gRNA_total; visualize=false)   
    d_act = Binomial(1, f_act) # there is a probability f_act that a gRNA is active
    d_high_act = truncated(Normal(ϵ_edit_act, sd_act), 0.01, 1)  # average genome editing efficiency for active gRNAs is equal to ϵ_edit_act
    d_low_act = truncated(Normal(ϵ_edit_inact, sd_act), 0.01, 1) # average genome editing efficiency for inactive gRNAs is equal to ϵ_edit_inact
    p_gRNA_edit = zeros(n_gRNA_total) # initialize vector with genome editing efficiencies for gRNAs
    for i in 1:n_gRNA_total
        if rand(d_act, 1) == [1]  # gRNA is active
            p_gRNA_edit[i] = rand(d_high_act, 1)[1]
        else  # gRNA is inactive
            p_gRNA_edit[i] = rand(d_low_act, 1)[1]
        end
    end
    if visualize
        return histogram(p_gRNA_edit, label="", xlabel="gRNA activity", ylabel="absolute frequency", title="gRNA activity distribution")

    else
        return p_gRNA_edit
    end
end

"""
    simulate_Nₓ₁(x, 
                g, 
                r, 
                n_gRNA_total, 
                p_gRNA_freq, 
                p_gRNA_edit, 
                ϵ_KO; 
                iter = 500)

Simulation-based approach for calculating Nₓ₁ of a CRISPR/Cas experiment.
"""
function simulate_Nₓ₁(x, 
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO; iter=500)
    """ 
    INPUT
  
    
    OUTPUT
    E: expected minimum number of plants to gRNA_read_distributionsee each pairwise combination at least once
    sd: standard deviation on the minimum number of plants
    """
    @assert x * g == n_gRNA_total
    
    T_vec = [] #stores number of plants required for each experiment
        for i in 1:iter       
            genes_vec = [] # Initialize matrix to count pairwise interactions
            T = 0
            while genes_vec != collect(1:x) # check if all pairwise combinations are present
             
                T += 1 # count how many plants must be sampled to fill pairwise interaction matrix
                
                # sample combinatorial gRNA/Cas9 construct
                gRNA_indices_construct = findall((rand(Multinomial(r, p_gRNA_freq))) .!= 0)
                
                # execute mutations
                gRNA_indices_mutations = [gRNA for gRNA in gRNA_indices_construct if rand(Binomial(1, p_gRNA_edit[gRNA])) == 1]
            
                # effective gene knockout (loss-of-function) ?
                gRNA_indices_KO = [gRNA for gRNA in gRNA_indices_mutations if rand(Binomial(1, ϵ_KO)) == 1]
            
                # which genes are knocked out?
                genes_indices_KO = Int.(ceil.(gRNA_indices_KO / g))
                 append!(genes_vec, genes_indices_KO)
                genes_vec = Int.(sort(unique(genes_vec)))
            end
            push!(T_vec, T)   
              
        end
        E = mean(T_vec); sd = std(T_vec)
    return E, sd
end
   

"""
    BioCCP_Nₓ₁(x, 
                g, 
                r, 
                n_gRNA_total, 
                p_gRNA_freq, 
                p_gRNA_edit, 
                ϵ_KO)

Employs the BioCCP.jl package to calculate the Nₓ₁ of a CRISPR/Cas experiment.
"""
function BioCCP_Nₓ₁(x, 
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO)
    
    p_gRNAs = p_gRNA_freq .* p_gRNA_edit * ϵ_KO
    p_genes = [sum(p_gRNAs[i:i+g-1]) for i in 1:g:n_gRNA_total]
    return expectation_minsamplesize(x; p=p_genes, r=r, normalize=false), 
    std_minsamplesize(x; p=p_genes, r=r, normalize=false)
end


"""
    simulate_Nₓ₂(x, 
                    g, 
                    r, 
                    n_gRNA_total, 
                    p_gRNA_freq, 
                    p_gRNA_edit, 
                    ϵ_KO; 
                    iter=500)

Simulation-based approach for calculating Nₓ₂ of multiplex CRISPR/Cas experiment.
"""
function simulate_Nₓ₂(x, 
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO; iter=500)
    """ 
    INPUT
  
    
    OUTPUT
    E: expected minimum number of plants to gRNA_read_distributionsee each pairwise combination at least once
    sd: standard deviation on the minimum number of plants
    """
    @assert x * g == n_gRNA_total
#     @assert sum(p_gRNA_freq) == 1
    
    T_vec = [] #stores number of plants required for each experiment
        for i in 1:iter     
            X_interactions_count = zeros(x, x) # Initialize matrix to count pairwise interactions
            T = 0
            while X_interactions_count != ones(x, x) # check if all pairwise combinations are present
                T += 1 # count how many plants must be sampled to fill pairwise interaction matrix
                
                # sample combinatorial gRNA/Cas9 construct
                gRNA_indices_construct = findall((rand(Multinomial(r, p_gRNA_freq))) .!= 0)
                
                # execute mutations
                gRNA_indices_mutations = [gRNA for gRNA in gRNA_indices_construct if rand(Binomial(1, p_gRNA_edit[gRNA])) == 1]
            
                # effective gene knockout (loss-of-function) ?
                gRNA_indices_KO = [gRNA for gRNA in gRNA_indices_mutations if rand(Binomial(1, ϵ_KO)) == 1]
            
                # which genes are knocked out?
                genes_indices_KO = Int.(ceil.(gRNA_indices_KO / g))
            
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
    BioCCP_Nₓ₂(x, 
                g, 
                r, 
                n_gRNA_total, 
                p_gRNA_freq, 
                p_gRNA_edit, ϵ_KO)

Employs the BioCCP.jl package to calculate the Nₓ₂ of a multiplex CRISPR/Cas experiment.
"""
function BioCCP_Nₓ₂(x, 
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO)
    
    # how many pairwise combinations of gRNAs
    ind_combinations_gRNA = collect(combinations(1:n_gRNA_total, 2))
    n_combinations_gRNA = length(ind_combinations_gRNA)
    
    # calculate probability and activity of gRNA combinations
    p_combinations_gRNA_library = zeros(n_combinations_gRNA)
    p_combinations_gRNA_act = zeros(n_combinations_gRNA)
    for i in 1:n_combinations_gRNA
        p_combinations_gRNA_library[i] = p_gRNA_freq[ind_combinations_gRNA[i][1]] * p_gRNA_freq[ind_combinations_gRNA[i][2]]
        p_combinations_gRNA_act[i] = p_gRNA_edit[ind_combinations_gRNA[i][1]] * p_gRNA_edit[ind_combinations_gRNA[i][2]]
    end
    
    # normalize probability gRNA combinations
    p_combinations_gRNA_library /= sum(p_combinations_gRNA_library)

    # select pairwise gRNA combinations of which each component codes for different gene (goal is to study combinations of knockouts in different genes)
    p_combinations_gRNA_library_interest = []
    p_combinations_gRNA_act_interest = []
    ind_combinations_gRNA_interest = []
    for i in 1:n_combinations_gRNA
        if ceil(ind_combinations_gRNA[i][1]/g) != ceil(ind_combinations_gRNA[i][2]/g)
            push!(p_combinations_gRNA_library_interest, p_combinations_gRNA_library[i])
            push!(p_combinations_gRNA_act_interest, p_combinations_gRNA_act[i])
            push!(ind_combinations_gRNA_interest, ind_combinations_gRNA[i])
        end
    end
        
    n_combinations_gRNA_interest = length(p_combinations_gRNA_library_interest)
    p_combinations_gRNA = p_combinations_gRNA_library_interest .* p_combinations_gRNA_act_interest * ϵ_KO^2

    # sum up probabilities or gRNA combinations for corresponding gene knockout combinations
    p_genes_matrix = zeros(x, x)
    for i in 1:n_combinations_gRNA_interest
        gene1 = Int(ceil(ind_combinations_gRNA_interest[i][1]/g))
        gene2 = Int(ceil(ind_combinations_gRNA_interest[i][2]/g))
        p_genes_matrix[gene1, gene2] += p_combinations_gRNA[i]
    end
    p_genes = collect([p_genes_matrix[i, j] for j in 2:size(p_genes_matrix, 1) for i in 1:j-1])  
    n_combinations_genes = length(p_genes)
    combinations_pp = length(collect(combinations(1:r, 2)))
    
    return expectation_minsamplesize(n_combinations_genes; p=p_genes, r=combinations_pp, normalize=false), std_minsamplesize(n_combinations_genes; p=p_genes, r=combinations_pp, normalize=false)
end

function simulate_Nₓ₂_countKOs(x, 
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO; iter=500)
     
    @assert x * g == n_gRNA_total
    
            n_KOs = []
       
            for j in 1:100000
                               
                # sample combinatorial gRNA/Cas9 construct
                gRNA_indices_construct = findall((rand(Multinomial(r, p_gRNA_freq))) .!= 0)
                
                # execute mutations
                gRNA_indices_mutations = [gRNA for gRNA in gRNA_indices_construct if rand(Binomial(1, p_gRNA_edit[gRNA])) == 1]
            
                # effective gene knockout (loss-of-function) ?
                gRNA_indices_KO = [gRNA for gRNA in gRNA_indices_mutations if rand(Binomial(1, ϵ_KO)) == 1]
            
                # which genes are knocked out?
                genes_indices_KO = Int.(ceil.(gRNA_indices_KO / g))
            
               push!(n_KOs, length(unique((genes_indices_KO))))
            end  
 
    return n_KOs
end

"""
    simulate_Nₓ₃(x, 
                    g, 
                    r, 
                    n_gRNA_total, 
                    p_gRNA_freq, 
                    p_gRNA_edit, 
                    ϵ_KO; 
                    iter=500)

Simulation-based approach for calculating Nₓ₃ of multiplex CRISPR/Cas experiment.
"""
function simulate_Nₓ₃(x, 
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO; iter=500)
    """ 
    INPUT
  
    
    OUTPUT
    E: expected minimum number of plants to gRNA_read_distributionsee each pairwise combination at least once
    sd: standard deviation on the minimum number of plants
    """
    @assert x * g == n_gRNA_total
#     @assert sum(p_gRNA_freq) == 1
    
    T_vec = [] #stores number of plants required for each experiment
        for i in 1:iter       
            X_interactions_count = zeros(x, x, x) # Initialize matrix to count triple interactions
            T = 0
            while X_interactions_count != ones(x, x, x) # check if all triple combinations are present
                T += 1 # count how many plants must be sampled to fill triple interaction matrix
                
                # sample combinatorial gRNA/Cas9 construct
                gRNA_indices_construct = findall((rand(Multinomial(r, p_gRNA_freq))) .!= 0)
                
                # execute mutations
                gRNA_indices_mutations = [gRNA for gRNA in gRNA_indices_construct if rand(Binomial(1, p_gRNA_edit[gRNA])) == 1]
            
                # effective gene knockout (loss-of-function) ?
                gRNA_indices_KO = [gRNA for gRNA in gRNA_indices_mutations if rand(Binomial(1, ϵ_KO)) == 1]
            
                # which genes are knocked out?
                genes_indices_KO = Int.(ceil.(gRNA_indices_KO / g))
            
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
    BioCCP_Nₓ₃(x, 
            g, 
            r, 
            n_gRNA_total, 
            p_gRNA_freq, 
            p_gRNA_edit, 
            ϵ_KO)

Employs the BioCCP.jl package to calculate the RES3 of a multiplex CRISPR/Cas experiment.
"""
function BioCCP_Nₓ₃(x, 
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO)
    
    # how many pairwise combinations of gRNAs
    ind_combinations_gRNA = collect(combinations(1:n_gRNA_total, 3))
    n_combinations_gRNA = length(ind_combinations_gRNA)
    
    # calculate probability and activity of triple gRNA combinations
    p_combinations_gRNA_library = zeros(n_combinations_gRNA)
    p_combinations_gRNA_act = zeros(n_combinations_gRNA)
    for i in 1:n_combinations_gRNA
        p_combinations_gRNA_library[i] = p_gRNA_freq[ind_combinations_gRNA[i][1]] * p_gRNA_freq[ind_combinations_gRNA[i][2]] * p_gRNA_freq[ind_combinations_gRNA[i][3]]
        p_combinations_gRNA_act[i] = p_gRNA_edit[ind_combinations_gRNA[i][1]] * p_gRNA_edit[ind_combinations_gRNA[i][2]] * p_gRNA_edit[ind_combinations_gRNA[i][3]]
    end
    
    # normalize probability gRNA combinations
    p_combinations_gRNA_library /= sum(p_combinations_gRNA_library)

    # select triple gRNA combinations of which each component codes for different gene (goal is to study combinations of knockouts in different genes)
    p_combinations_gRNA_library_interest = []
    p_combinations_gRNA_act_interest = []
    ind_combinations_gRNA_interest = []
    for i in 1:n_combinations_gRNA
        if ceil(ind_combinations_gRNA[i][1]/g) != ceil(ind_combinations_gRNA[i][2]/g) && ceil(ind_combinations_gRNA[i][1]/g) != ceil(ind_combinations_gRNA[i][3]/g) && ceil(ind_combinations_gRNA[i][3]/g) != ceil(ind_combinations_gRNA[i][2]/g)
            push!(p_combinations_gRNA_library_interest, p_combinations_gRNA_library[i])
            push!(p_combinations_gRNA_act_interest, p_combinations_gRNA_act[i])
            push!(ind_combinations_gRNA_interest, ind_combinations_gRNA[i])
        end
    end
        
    n_combinations_gRNA_interest = length(p_combinations_gRNA_library_interest)
    p_combinations_gRNA = p_combinations_gRNA_library_interest .* p_combinations_gRNA_act_interest * ϵ_KO^3

    # sum up probabilities or gRNA combinations for corresponding gene knockout combinations
    p_genes_matrix = zeros(x, x, x)
    for i in 1:n_combinations_gRNA_interest
        gene1 = Int(ceil(ind_combinations_gRNA_interest[i][1]/g))
        gene2 = Int(ceil(ind_combinations_gRNA_interest[i][2]/g))
        gene3 = Int(ceil(ind_combinations_gRNA_interest[i][3]/g))
        p_genes_matrix[gene1, gene2, gene3] += p_combinations_gRNA[i]
    end
    
    combinations_genes = collect(combinations(1:x, 3))
    p_genes = []
        for combination in combinations_genes
            push!(p_genes, p_genes_matrix[combination[1], combination[2], combination[3]])
        end
        
    n_combinations_genes = length(p_genes)
    combinations_pp = length(collect(combinations(1:r, 3)))
    
    return expectation_minsamplesize(n_combinations_genes; p=p_genes, r=combinations_pp, normalize=false), std_minsamplesize(n_combinations_genes; p=p_genes, r=combinations_pp, normalize=false)
end

"""
    BioCCP_Pₓ₁(x, 
                N,
                g, 
                r, 
                n_gRNA_total, 
                p_gRNA_freq, 
                p_gRNA_edit, 
                ϵ_KO)
Employs the BioCCP.jl package to calculate the probability of full coverage of all individual gene knockouts 
in the gene knockout library of a CRISPR/Cas experiment 
with respect to a specified plant library size (number of plants analyzed).
"""
function BioCCP_Pₓ₁(x, N,
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO)
    
    p_gRNAs = p_gRNA_freq .* p_gRNA_edit * ϵ_KO
    p_genes = [sum(p_gRNAs[i:i+g-1]) for i in 1:g:n_gRNA_total]
    return success_probability(x, N; p=p_genes, r=r, normalize=false) 
end

"""
BioCCP_γₓ₁(x, 
            N,
            g, 
            r, 
            n_gRNA_total, 
            p_gRNA_freq, 
            p_gRNA_edit, 
            ϵ_KO)

Employs the BioCCP.jl package to calculate the expected represented fraction of the design space, 
here consisting of all individual gene knockouts in a gene knockout library of a CRISPR/Cas experiment, 
with respect to a specified experimental scale (number of plants analyzed; ES).
"""
function BioCCP_γₓ₁(x, N,
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO)
    
    p_gRNAs = p_gRNA_freq .* p_gRNA_edit * ϵ_KO
    p_genes = [sum(p_gRNAs[i:i+g-1]) for i in 1:g:n_gRNA_total]
    return expectation_fraction_collected(x, N; p=p_genes, r=r, normalize=false) 
end


"""
    BioCCP_Pₓ₂(x, 
                N,
                g, 
                r, 
                n_gRNA_total, 
                p_gRNA_freq, 
                p_gRNA_edit, 
                ϵ_KO)

Employs the BioCCP.jl package to calculate the probability of full coverage of all pairwise combinations of gene knockouts 
in the combinatorial gene knockout library of a multiplex CRISPR/Cas experiment 
with respect to a specified experimental scale (number of plants analyzed; ES).
"""
function BioCCP_Pₓ₂(x, N,
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO)
    
    # how many pairwise combinations of gRNAs
    ind_combinations_gRNA = collect(combinations(1:n_gRNA_total, 2))
    n_combinations_gRNA = length(ind_combinations_gRNA)
    
    # calculate probability and activity of gRNA combinations
    p_combinations_gRNA_library = zeros(n_combinations_gRNA)
    p_combinations_gRNA_act = zeros(n_combinations_gRNA)
    for i in 1:n_combinations_gRNA
        p_combinations_gRNA_library[i] = p_gRNA_freq[ind_combinations_gRNA[i][1]] * p_gRNA_freq[ind_combinations_gRNA[i][2]]
        p_combinations_gRNA_act[i] = p_gRNA_edit[ind_combinations_gRNA[i][1]] * p_gRNA_edit[ind_combinations_gRNA[i][2]]
    end
    
    # normalize probability gRNA combinations
    p_combinations_gRNA_library /= sum(p_combinations_gRNA_library)

    # select pairwise gRNA combinations of which each component codes for different gene (goal is to study combinations of knockouts in different genes)
    p_combinations_gRNA_library_interest = []
    p_combinations_gRNA_act_interest = []
    ind_combinations_gRNA_interest = []
    for i in 1:n_combinations_gRNA
        if ceil(ind_combinations_gRNA[i][1]/g) != ceil(ind_combinations_gRNA[i][2]/g)
            push!(p_combinations_gRNA_library_interest, p_combinations_gRNA_library[i])
            push!(p_combinations_gRNA_act_interest, p_combinations_gRNA_act[i])
            push!(ind_combinations_gRNA_interest, ind_combinations_gRNA[i])
        end
    end
        
    n_combinations_gRNA_interest = length(p_combinations_gRNA_library_interest)
    p_combinations_gRNA = p_combinations_gRNA_library_interest .* p_combinations_gRNA_act_interest * ϵ_KO^2

    # sum up probabilities or gRNA combinations for corresponding gene knockout combinations
    p_genes_matrix = zeros(x, x)
    for i in 1:n_combinations_gRNA_interest
        gene1 = Int(ceil(ind_combinations_gRNA_interest[i][1]/g))
        gene2 = Int(ceil(ind_combinations_gRNA_interest[i][2]/g))
        p_genes_matrix[gene1, gene2] += p_combinations_gRNA[i]
    end

    p_genes = collect([p_genes_matrix[i, j] for j in 2:size(p_genes_matrix, 1) for i in 1:j-1])  
    n_combinations_genes = length(p_genes)
    combinations_pp = length(collect(combinations(1:r, 2)))
    
    return success_probability(n_combinations_genes, N; p=p_genes, r=combinations_pp, normalize=false)
end

"""
    BioCCP_γₓ₂(x, 
                N,
                g, 
                r, 
                n_gRNA_total, 
                p_gRNA_freq, 
                p_gRNA_edit, 
                ϵ_KO)

Employs the BioCCP.jl package to calculate the expected represented fraction of the design space, here consisting of all pairwise combinations of gene knockouts 
in the combinatorial gene knockout library of a multiplex CRISPR/Cas experiment, 
with respect to a specified plant library size (number of plants analyzed).
"""
function BioCCP_γₓ₂(x, N,
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO)
    
    # how many pairwise combinations of gRNAs
    ind_combinations_gRNA = collect(combinations(1:n_gRNA_total, 2))
    n_combinations_gRNA = length(ind_combinations_gRNA)
    
    # calculate probability and activity of gRNA combinations
    p_combinations_gRNA_library = zeros(n_combinations_gRNA)
    p_combinations_gRNA_act = zeros(n_combinations_gRNA)
    for i in 1:n_combinations_gRNA
        p_combinations_gRNA_library[i] = p_gRNA_freq[ind_combinations_gRNA[i][1]] * p_gRNA_freq[ind_combinations_gRNA[i][2]]
        p_combinations_gRNA_act[i] = p_gRNA_edit[ind_combinations_gRNA[i][1]] * p_gRNA_edit[ind_combinations_gRNA[i][2]]
    end
    
    # normalize probability gRNA combinations
    p_combinations_gRNA_library /= sum(p_combinations_gRNA_library)

    # select pairwise gRNA combinations of which each component codes for different gene (goal is to study combinations of knockouts in different genes)
    p_combinations_gRNA_library_interest = []
    p_combinations_gRNA_act_interest = []
    ind_combinations_gRNA_interest = []
    for i in 1:n_combinations_gRNA
        if ceil(ind_combinations_gRNA[i][1]/g) != ceil(ind_combinations_gRNA[i][2]/g)
            push!(p_combinations_gRNA_library_interest, p_combinations_gRNA_library[i])
            push!(p_combinations_gRNA_act_interest, p_combinations_gRNA_act[i])
            push!(ind_combinations_gRNA_interest, ind_combinations_gRNA[i])
        end
    end
        
    n_combinations_gRNA_interest = length(p_combinations_gRNA_library_interest)
    p_combinations_gRNA = p_combinations_gRNA_library_interest .* p_combinations_gRNA_act_interest * ϵ_KO^2

    # sum up probabilities or gRNA combinations for corresponding gene knockout combinations
    p_genes_matrix = zeros(x, x)
    for i in 1:n_combinations_gRNA_interest
        gene1 = Int(ceil(ind_combinations_gRNA_interest[i][1]/g))
        gene2 = Int(ceil(ind_combinations_gRNA_interest[i][2]/g))
        p_genes_matrix[gene1, gene2] += p_combinations_gRNA[i]
    end

    p_genes = collect([p_genes_matrix[i, j] for j in 2:size(p_genes_matrix, 1) for i in 1:j-1])  
    n_combinations_genes = length(p_genes)
    combinations_pp = length(collect(combinations(1:r, 2)))
    
    return expectation_fraction_collected(n_combinations_genes, N; p=p_genes, r=combinations_pp, normalize=false)
end


"""
    BioCCP_γₓ₃(x, 
                N,
                g, 
                r, 
                n_gRNA_total, 
                p_gRNA_freq, 
                p_gRNA_edit, 
                ϵ_KO)

Employs the BioCCP.jl package to calculate the RES3 of a multiplex CRISPR/Cas experiment.
"""
function BioCCP_γₓ₃(x, 
                N,
                g, 
                r, 
                n_gRNA_total, 
                p_gRNA_freq, 
                p_gRNA_edit, 
                ϵ_KO)
    
    # how many pairwise combinations of gRNAs
    ind_combinations_gRNA = collect(combinations(1:n_gRNA_total, 3))
    n_combinations_gRNA = length(ind_combinations_gRNA)
    
    # calculate probability and activity of triple gRNA combinations
    p_combinations_gRNA_library = zeros(n_combinations_gRNA)
    p_combinations_gRNA_act = zeros(n_combinations_gRNA)
    for i in 1:n_combinations_gRNA
        p_combinations_gRNA_library[i] = p_gRNA_freq[ind_combinations_gRNA[i][1]] * p_gRNA_freq[ind_combinations_gRNA[i][2]] * p_gRNA_freq[ind_combinations_gRNA[i][3]]
        p_combinations_gRNA_act[i] = p_gRNA_edit[ind_combinations_gRNA[i][1]] * p_gRNA_edit[ind_combinations_gRNA[i][2]] * p_gRNA_edit[ind_combinations_gRNA[i][3]]
    end
    
    # normalize probability gRNA combinations
    p_combinations_gRNA_library /= sum(p_combinations_gRNA_library)

    # select triple gRNA combinations of which each component codes for different gene (goal is to study combinations of knockouts in different genes)
    p_combinations_gRNA_library_interest = []
    p_combinations_gRNA_act_interest = []
    ind_combinations_gRNA_interest = []
    for i in 1:n_combinations_gRNA
        if ceil(ind_combinations_gRNA[i][1]/g) != ceil(ind_combinations_gRNA[i][2]/g) && ceil(ind_combinations_gRNA[i][1]/g) != ceil(ind_combinations_gRNA[i][3]/g) && ceil(ind_combinations_gRNA[i][3]/g) != ceil(ind_combinations_gRNA[i][2]/g)
            push!(p_combinations_gRNA_library_interest, p_combinations_gRNA_library[i])
            push!(p_combinations_gRNA_act_interest, p_combinations_gRNA_act[i])
            push!(ind_combinations_gRNA_interest, ind_combinations_gRNA[i])
        end
    end
        
    n_combinations_gRNA_interest = length(p_combinations_gRNA_library_interest)
    p_combinations_gRNA = p_combinations_gRNA_library_interest .* p_combinations_gRNA_act_interest * ϵ_KO^3

    # sum up probabilities or gRNA combinations for corresponding gene knockout combinations
    p_genes_matrix = zeros(x, x, x)
    for i in 1:n_combinations_gRNA_interest
        gene1 = Int(ceil(ind_combinations_gRNA_interest[i][1]/g))
        gene2 = Int(ceil(ind_combinations_gRNA_interest[i][2]/g))
        gene3 = Int(ceil(ind_combinations_gRNA_interest[i][3]/g))
        p_genes_matrix[gene1, gene2, gene3] += p_combinations_gRNA[i]
    end
    
    combinations_genes = collect(combinations(1:x, 3))
    p_genes = []
        for combination in combinations_genes
            push!(p_genes, p_genes_matrix[combination[1], combination[2], combination[3]])
        end
        
    n_combinations_genes = length(p_genes)
    combinations_pp = length(collect(combinations(1:r, 3)))
    
    return expectation_fraction_collected(n_combinations_genes, N; p=p_genes, r=combinations_pp, normalize=false)
end

"""
    BioCCP_Pₓ₃(x, 
                N,
                g, 
                r, 
                n_gRNA_total, 
                p_gRNA_freq, 
                p_gRNA_edit, 
                ϵ_KO)

Employs the BioCCP.jl package to calculate the RES3 of a multiplex CRISPR/Cas experiment.
"""
function BioCCP_Pₓ₃(x, 
                N,
                g, 
                r, 
                n_gRNA_total, 
                p_gRNA_freq, 
                p_gRNA_edit, 
                ϵ_KO)
    
    # how many pairwise combinations of gRNAs
    ind_combinations_gRNA = collect(combinations(1:n_gRNA_total, 3))
    n_combinations_gRNA = length(ind_combinations_gRNA)
    
    # calculate probability and activity of triple gRNA combinations
    p_combinations_gRNA_library = zeros(n_combinations_gRNA)
    p_combinations_gRNA_act = zeros(n_combinations_gRNA)
    for i in 1:n_combinations_gRNA
        p_combinations_gRNA_library[i] = p_gRNA_freq[ind_combinations_gRNA[i][1]] * p_gRNA_freq[ind_combinations_gRNA[i][2]] * p_gRNA_freq[ind_combinations_gRNA[i][3]]
        p_combinations_gRNA_act[i] = p_gRNA_edit[ind_combinations_gRNA[i][1]] * p_gRNA_edit[ind_combinations_gRNA[i][2]] * p_gRNA_edit[ind_combinations_gRNA[i][3]]
    end
    
    # normalize probability gRNA combinations
    p_combinations_gRNA_library /= sum(p_combinations_gRNA_library)

    # select triple gRNA combinations of which each component codes for different gene (goal is to study combinations of knockouts in different genes)
    p_combinations_gRNA_library_interest = []
    p_combinations_gRNA_act_interest = []
    ind_combinations_gRNA_interest = []
    for i in 1:n_combinations_gRNA
        if ceil(ind_combinations_gRNA[i][1]/g) != ceil(ind_combinations_gRNA[i][2]/g) && ceil(ind_combinations_gRNA[i][1]/g) != ceil(ind_combinations_gRNA[i][3]/g) && ceil(ind_combinations_gRNA[i][3]/g) != ceil(ind_combinations_gRNA[i][2]/g)
            push!(p_combinations_gRNA_library_interest, p_combinations_gRNA_library[i])
            push!(p_combinations_gRNA_act_interest, p_combinations_gRNA_act[i])
            push!(ind_combinations_gRNA_interest, ind_combinations_gRNA[i])
        end
    end
        
    n_combinations_gRNA_interest = length(p_combinations_gRNA_library_interest)
    p_combinations_gRNA = p_combinations_gRNA_library_interest .* p_combinations_gRNA_act_interest * ϵ_KO^3

    # sum up probabilities or gRNA combinations for corresponding gene knockout combinations
    p_genes_matrix = zeros(x, x, x)
    for i in 1:n_combinations_gRNA_interest
        gene1 = Int(ceil(ind_combinations_gRNA_interest[i][1]/g))
        gene2 = Int(ceil(ind_combinations_gRNA_interest[i][2]/g))
        gene3 = Int(ceil(ind_combinations_gRNA_interest[i][3]/g))
        p_genes_matrix[gene1, gene2, gene3] += p_combinations_gRNA[i]
    end
    
    combinations_genes = collect(combinations(1:x, 3))
    p_genes = []
        for combination in combinations_genes
            push!(p_genes, p_genes_matrix[combination[1], combination[2], combination[3]])
        end
        
    n_combinations_genes = length(p_genes)
    combinations_pp = length(collect(combinations(1:r, 3)))
    
    return success_probability(n_combinations_genes, N; p=p_genes, r=combinations_pp, normalize=false)
end
      
      
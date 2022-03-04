"""
    gRNA_frequency_distribution(m, sd, l, u, n_gRNA_total; normalize = true, visualize=false)

Generates vector with frequencies in the combinatorial gRNA/Cas9 construct library for all gRNAs

***INPUT***
m: the average abundance of the gRNAs (in terms of absolute or relative frequency)
sd: the standard deviation on the gRNA abundances (in terms of absolute or relative frequency)
l: minimal gRNA abundance (in terms of absolute or relative frequency)
u: maximal gRNA abundance (in terms of absolute or relative frequency)
n_gRNA_total: the total number of gRNAs in the experiment
normalize: if set to "true", the gRNA abundances (absolute frequencies) are converted into relative frequencies
visualize: if set to "true", a histogram of all gRNA abundances is plotted

***OUTPUT***
p_gRNA_freq: vector with frequencies for all gRNAs in the construct library
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
            ylabel="Frequency", title="gRNA frequency distribution")
    else
        return p_gRNA_freq
    end
end

"""
    gRNA_edit_distribution(f_act, ϵ_edit_act, ϵ_edit_inact, sd_act, n_gRNA_total; visualize=false)   

Generates vector with genome editing efficiencies for all the gRNAs in the experiment. 

***INPUT***
f_act: fraction of all gRNAs that is active
ϵ_edit_act: Average genome editing efficiency for active gRNAs - mean of the genome editing efficiency distribution for active gRNAs
ϵ_edit_inact: Average genome editing efficiency for inactive gRNAs - mean of the genome editing efficiency distribution for inactive gRNAs
sd_act: standard deviation of the genome editing efficiency distributions for active and inactive gRNAs
n_gRNA_total: the total number of gRNAs in the experiment
visualize: if set to "true", a histogram of all genome editing efficiency is plotted

***OUTPUT***
p_gRNA_edit: vector with genome editing efficiencies for all gRNAs
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

Computes the expected value and the standard deviation of the minimal plant library size for full coverage of all single gene knockouts (E[Nx,1] and σ[Nx,1]) using simulation

***INPUT***
x: number of target genes in the experiment
g: number of gRNAs designed per target gene
r: number of gRNA sequences per combinatorial gRNA/Cas construct
n_gRNA_total: total number of gRNAs in the experiment
p_gRNA_freq: vector with relative frequencies for all gRNAs in the construct library (normalized!)
p_gRNA_edit: vector with genome editing efficiencies for all gRNAs 
ϵ_KO: global knockout efficiency; fraction of mutations leading to effective gene knockout
iter: number of CRISPR/Cas experiments that are simulated to obtain E[Nₓ₁] and σ[Nₓ₁]

***OUTPUT***
E_Nₓ₁: expected value of the plant library size for full coverage of all single gene knockouts
sd_Nₓ₁: standard deviation on the plant library size for full coverage of all single gene knockouts
"""
function simulate_Nₓ₁(x, 
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO; iter=500)

    @assert x * g == n_gRNA_total     
    Nₓ₁_vec = [] #stores number of plants to reach full coverage for each simulated experiment
        for i in 1:iter       
            genes_vec = [] # Initialize vector to store single gene knockouts that are observed in plants
            Nₓ₁ = 0
            while genes_vec != collect(1:x) # check if all possible single gene knockouts are present: if no full coverage, sample an additional plant           
                Nₓ₁ += 1 # count how many plants must be sampled to observe all single gene knockouts
                
                # sample combinatorial gRNA/Cas9 construct
                gRNA_indices_construct = findall((rand(Multinomial(r, p_gRNA_freq))) .!= 0)
                
                # execute mutations
                gRNA_indices_mutations = [gRNA for gRNA in gRNA_indices_construct if rand(Binomial(1, p_gRNA_edit[gRNA])) == 1]
            
                # effective gene knockout (loss-of-function) ?
                gRNA_indices_KO = [gRNA for gRNA in gRNA_indices_mutations if rand(Binomial(1, ϵ_KO)) == 1]
            
                # which genes are knocked out?
                genes_indices_KO = Int.(ceil.(gRNA_indices_KO / g)) 
                append!(genes_vec, genes_indices_KO)

                # update vector with observed gene knockouts
                genes_vec = Int.(sort(unique(genes_vec)))
            end
            push!(Nₓ₁_vec, Nₓ₁) # add plant library size for full coverage of current experiment to vector         
        end

    # Calculate expected value and standard deviation
    E_Nₓ₁ = mean(Nₓ₁_vec); 
    sd_Nₓ₁ = std(Nₓ₁_vec)

    return E_Nₓ₁, sd_Nₓ₁
end
   

"""
    BioCCP_Nₓ₁(x, 
                g, 
                r, 
                n_gRNA_total, 
                p_gRNA_freq, 
                p_gRNA_edit, 
                ϵ_KO)

Computes the expected value and the standard deviation of the minimal plant library size for 
full coverage of all single gene knockouts (E[Nx,1] and σ[Nx,1]) using BioCCP

***INPUT***
x: number of target genes in the experiment
g: number of gRNAs designed per target gene
r: number of gRNA sequences per combinatorial gRNA/Cas construct
n_gRNA_total: total number of gRNAs in the experiment
p_gRNA_freq: vector with relative frequencies for all gRNAs in the construct library (normalized!)
p_gRNA_edit: vector with genome editing efficiencies for all gRNAs 
ϵ_KO: global knockout efficiency; fraction of mutations leading to effective gene knockout

***OUTPUT***
E_Nₓ₁ : expected value of the plant library size for full coverage of all single gene knockouts
sd_Nₓ₁ : standard deviation on the plant library size for full coverage of all single gene knockouts
"""
function BioCCP_Nₓ₁(x, 
                    g, 
                    r, 
                    n_gRNA_total, 
                    p_gRNA_freq, 
                    p_gRNA_edit, 
                    ϵ_KO)
    
    # prepare input for BioCCP functions
    p_gRNAs = p_gRNA_freq .* p_gRNA_edit * ϵ_KO  # calculate probability for each gRNA to induce effective gene knockout
    p_genes = [sum(p_gRNAs[i:i+g-1]) for i in 1:g:n_gRNA_total]  # obtain probability of single gene knockout by summing up probability of all corresponding gRNAs to induce effective gene knockout
    
    # Apply BioCCP functions
    E_Nₓ₁ = expectation_minsamplesize(x; p=p_genes, r=r, normalize=false)
    sd_Nₓ₁ = std_minsamplesize(x; p=p_genes, r=r, normalize=false)

    return E_Nₓ₁, sd_Nₓ₁
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

Computes the expected value and the standard deviation of the minimal plant library size for full coverage of all pairwise combinations of gene knockouts 
in a multiplex CRISPR/Cas experiment (E[Nx,2] and σ[Nx,2]) using simulation

***INPUT***
x: number of target genes in the experiment
g: number of gRNAs designed per target gene
r: number of gRNA sequences per combinatorial gRNA/Cas construct
n_gRNA_total: total number of gRNAs in the experiment
p_gRNA_freq: vector with relative frequencies for all gRNAs in the construct library (normalized!)
p_gRNA_edit: vector with genome editing efficiencies for all gRNAs 
ϵ_KO: global knockout efficiency; fraction of mutations leading to effective gene knockout
iter: number of CRISPR/Cas experiments that are simulated to obtain E[Nₓ₂] and σ[Nₓ₂]

***OUTPUT***
E_Nₓ₂ : expected value of the plant library size for full coverage of all pairwise combinations of gene knockouts
sd_Nₓ₂ : standard deviation on the plant library size for full coverage of all pairwise combinations of gene knockouts
"""
function simulate_Nₓ₂(x, 
                    g, 
                    r, 
                    n_gRNA_total, 
                    p_gRNA_freq, 
                    p_gRNA_edit, 
                    ϵ_KO; 
                    iter=500)

    @assert x * g == n_gRNA_total    
    Nₓ₂_vec = [] #stores number of plants to reach full coverage of all pairwise combinations of gene knockouts for each simulated experiment
       
    for i in 1:iter     
        X_interactions_count = zeros(x, x) # Initialize matrix to count pairwise interactions
        Nₓ₂ = 0
        while X_interactions_count != ones(x, x) # check if all pairwise combinations are present
            Nₓ₂ += 1 # count how many plants must be sampled to fill pairwise interaction matrix
                
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
        push!(Nₓ₂_vec, Nₓ₂)               
        end

    # Calculate expected value and standard deviation
    E_Nₓ₂ = mean(Nₓ₂_vec)
    sd_Nₓ₂ = std(Nₓ₂_vec)
        
    return E_Nₓ₂, sd_Nₓ₂
end

"""
    BioCCP_Nₓ₂(x, 
                g, 
                r, 
                n_gRNA_total, 
                p_gRNA_freq, 
                p_gRNA_edit, ϵ_KO)

Computes the expected value and the standard deviation of the minimal plant library size for full coverage of all pairwise combinations of gene knockouts in a multiplex CRISPR/Cas experiment (E[Nx,2] and σ[Nx,2]) using BioCCP

    ***INPUT***
x: number of target genes in the experiment
g: number of gRNAs designed per target gene
r: number of gRNA sequences per combinatorial gRNA/Cas construct
n_gRNA_total: total number of gRNAs in the experiment
p_gRNA_freq: vector with relative frequencies for all gRNAs in the construct library (normalized!)
p_gRNA_edit: vector with genome editing efficiencies for all gRNAs 
ϵ_KO: global knockout efficiency; fraction of mutations leading to effective gene knockout

***OUTPUT***
E_Nₓ₂: expected value of the plant library size for full coverage of all pairwise combinations of gene knockouts
sd_Nₓ₂: standard deviation on the plant library size for full coverage of all pairwise combinations of gene knockouts
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
    
    # Apply BioCCP functions
    E_Nₓ₂ = expectation_minsamplesize(n_combinations_genes; p=p_genes, r=combinations_pp, normalize=false)
    sd_Nₓ₂ = std_minsamplesize(n_combinations_genes; p=p_genes, r=combinations_pp, normalize=false)

    return E_Nₓ₂, sd_Nₓ₂
end

"""
    simulate_Nₓ₂_countKOs(x, 
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO; iter=100000)

Counts the number of knockouts per plant in the experiment.

***INPUT***
x: number of target genes in the experiment
g: number of gRNAs designed per target gene
r: number of gRNA sequences per combinatorial gRNA/Cas construct
n_gRNA_total: total number of gRNAs in the experiment
p_gRNA_freq: vector with relative frequencies for all gRNAs in the construct library (normalized!)
p_gRNA_edit: vector with genome editing efficiencies for all gRNAs 
ϵ_KO: global knockout efficiency; fraction of mutations leading to effective gene knockout
iter: number of plants that are sampled 

***OUTPUT***
n_KOs_vec: vector with the number of knockouts for each plant 
"""
function simulate_Nₓ₂_countKOs(x, 
                                g, 
                                r, 
                                n_gRNA_total, 
                                p_gRNA_freq, 
                                p_gRNA_edit, 
                                ϵ_KO; 
                                iter = 100000)
     
    @assert x * g == n_gRNA_total 
    n_KOs_vec = []      
    for j in 1:iter                           
        # sample combinatorial gRNA/Cas9 construct
        gRNA_indices_construct = findall((rand(Multinomial(r, p_gRNA_freq))) .!= 0)
                
        # execute mutations
        gRNA_indices_mutations = [gRNA for gRNA in gRNA_indices_construct if rand(Binomial(1, p_gRNA_edit[gRNA])) == 1]
            
        # effective gene knockout (loss-of-function) ?
        gRNA_indices_KO = [gRNA for gRNA in gRNA_indices_mutations if rand(Binomial(1, ϵ_KO)) == 1]
            
        # which genes are knocked out?
        genes_indices_KO = Int.(ceil.(gRNA_indices_KO / g))
            
        push!(n_KOs_vec, length(unique((genes_indices_KO))))
    end  
    return n_KOs_vec
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

Computes the expected value and the standard deviation of the minimal plant library size for full coverage of all triple combinations of gene knockouts in 
a multiplex CRISPR/Cas experiment (E[Nx,3] and σ[Nx,3]) using simulation

***INPUT***
x: number of target genes in the experiment
g: number of gRNAs designed per target gene
r: number of gRNA sequences per combinatorial gRNA/Cas construct
n_gRNA_total: total number of gRNAs in the experiment
p_gRNA_freq: vector with relative frequencies for all gRNAs in the construct library (normalized!)
p_gRNA_edit: vector with genome editing efficiencies for all gRNAs 
ϵ_KO: global knockout efficiency; fraction of mutations leading to effective gene knockout
iter: number of CRISPR/Cas experiments that are simulated to obtain E[Nₓ₃] and σ[Nₓ₃]

***OUTPUT***
E_Nₓ₃: expecteded value of the plant library size for full coverage of all triple combinations of gene knockouts
sd_Nₓ₃: standard deviation on the plant library size for full coverage of all triple combinations of gene knockouts
"""
function simulate_Nₓ₃(x, 
                    g, 
                    r, 
                    n_gRNA_total, 
                    p_gRNA_freq, 
                    p_gRNA_edit, 
                    ϵ_KO; 
                    iter = 500)

    @assert x * g == n_gRNA_total
    Nₓ₃_vec = [] # stores number of plants required for each experiment

    for i in 1:iter       
        X_interactions_count = zeros(x, x, x) # Initialize matrix to count triple interactions
        Nₓ₃ = 0

        while X_interactions_count != ones(x, x, x) # check if all triple combinations are present
            Nₓ₃ += 1 # count how many plants must be sampled to fill triple interaction matrix
                
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
        push!(Nₓ₃_vec, Nₓ₃)       
    end

    # calculate expected value and standard deviation
    E_Nₓ₃ = mean(Nₓ₃_vec); sd_Nₓ₃ = std(Nₓ₃_vec)

    return E_Nₓ₃, sd_Nₓ₃
end
   

"""
    BioCCP_Nₓ₃(x, 
            g, 
            r, 
            n_gRNA_total, 
            p_gRNA_freq, 
            p_gRNA_edit, 
            ϵ_KO)

Computes the expected value and the standard deviation of the minimal plant library size 
for full coverage of all triple combinations of gene knockouts in a multiplex CRISPR/Cas experiment (E[Nx,3] and σ[Nx,3]) using BioCCP

***INPUT***
x: number of target genes in the experiment
g: number of gRNAs designed per target gene
r: number of gRNA sequences per combinatorial gRNA/Cas construct
n_gRNA_total: total number of gRNAs in the experiment
p_gRNA_freq: vector with relative frequencies for all gRNAs in the construct library (normalized!)
p_gRNA_edit: vector with genome editing efficiencies for all gRNAs 
ϵ_KO: global knockout efficiency; fraction of mutations leading to effective gene knockout

***OUTPUT***
E_Nₓ₃: expecteded value of the plant library size for full coverage of all triple combinations of gene knockouts
sd_Nₓ₃: standard deviation on the plant library size for full coverage of all triple combinations of gene knockouts
"""
function BioCCP_Nₓ₃(x, 
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO)
    
    # how many triple combinations of gRNAs
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
    
    # apply BioCCP functions
    E_Nₓ₃ = expectation_minsamplesize(n_combinations_genes; p=p_genes, r=combinations_pp, normalize=false)
    sd_Nₓ₃ = std_minsamplesize(n_combinations_genes; p=p_genes, r=combinations_pp, normalize=false)

    return E_Nₓ₃, sd_Nₓ₃
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

Computes the probability of full coverage of all single gene knockouts (Px,1) for an experiment with given plant library size using BioCCP

***INPUT***
x: number of target genes in the experiment
N: plant library size
g: number of gRNAs designed per target gene
r: number of gRNA sequences per combinatorial gRNA/Cas construct
n_gRNA_total: total number of gRNAs in the experiment
p_gRNA_freq: vector with relative frequencies for all gRNAs in the construct library (normalized!)
p_gRNA_edit: vector with genome editing efficiencies for all gRNAs 
ϵ_KO: global knockout efficiency; fraction of mutations leading to effective gene knockout

***OUTPUT***
Pₓ₁: probability of full coverage of all single gene knockouts

"""
function BioCCP_Pₓ₁(x, 
                    N,
                    g, 
                    r, 
                    n_gRNA_total, 
                    p_gRNA_freq, 
                    p_gRNA_edit, 
                    ϵ_KO)
    
    # prepare input
    p_gRNAs = p_gRNA_freq .* p_gRNA_edit * ϵ_KO
    p_genes = [sum(p_gRNAs[i:i+g-1]) for i in 1:g:n_gRNA_total]
    
    # apply BioCCP function
    Pₓ₁ = success_probability(x, N; p=p_genes, r=r, normalize=false) 

    return Pₓ₁
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

Computes the expected coverage of all single gene knockouts (Px,1) for an experiment with given plant library size using BioCCP

***INPUT***
x: number of target genes in the experiment
N: plant library size
g: number of gRNAs designed per target gene
r: number of gRNA sequences per combinatorial gRNA/Cas construct
n_gRNA_total: total number of gRNAs in the experiment
p_gRNA_freq: vector with relative frequencies for all gRNAs in the construct library (normalized!)
p_gRNA_edit: vector with genome editing efficiencies for all gRNAs 
ϵ_KO: global knockout efficiency; fraction of mutations leading to effective gene knockout

***OUTPUT***
γₓ₁: expected coverage of all single gene knockouts
"""
function BioCCP_γₓ₁(x, N,
                                         g, 
                                         r, 
                                         n_gRNA_total, 
                                         p_gRNA_freq, 
                                         p_gRNA_edit, ϵ_KO)
    
    p_gRNAs = p_gRNA_freq .* p_gRNA_edit * ϵ_KO
    p_genes = [sum(p_gRNAs[i:i+g-1]) for i in 1:g:n_gRNA_total]

    # Apply BioCCP function
    γₓ₁ = expectation_fraction_collected(x, N; p=p_genes, r=r, normalize=false) 

    return γₓ₁ 
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

Computes the probability of full coverage of all pairwise combinations of gene knockouts (Px,2) for an experiment with given plant library size using BioCCP

***INPUT***
x: number of target genes in the experiment
N: plant library size
g: number of gRNAs designed per target gene
r: number of gRNA sequences per combinatorial gRNA/Cas construct
n_gRNA_total: total number of gRNAs in the experiment
p_gRNA_freq: vector with relative frequencies for all gRNAs in the construct library (normalized!)
p_gRNA_edit: vector with genome editing efficiencies for all gRNAs 
ϵ_KO: global knockout efficiency; fraction of mutations leading to effective gene knockout
                
***OUTPUT***
Pₓ₂: probability of full coverage of all pairwise combinations of gene knockouts
"""
function BioCCP_Pₓ₂(x, 
                    N,
                    g, 
                    r, 
                    n_gRNA_total, 
                    p_gRNA_freq, 
                    p_gRNA_edit, 
                    ϵ_KO)
    
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
    
    # Apply BioCCP function
    Pₓ₂ = success_probability(n_combinations_genes, N; p=p_genes, r=combinations_pp, normalize=false)

    return Pₓ₂ 
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

Computes the expected coverage of all pairwise combinations of gene knockouts (γx,2) for an experiment with given plant library size using BioCCP

***INPUT***
x: number of target genes in the experiment
N: plant library size
g: number of gRNAs designed per target gene
r: number of gRNA sequences per combinatorial gRNA/Cas construct
n_gRNA_total: total number of gRNAs in the experiment
p_gRNA_freq: vector with relative frequencies for all gRNAs in the construct library (normalized!)
p_gRNA_edit: vector with genome editing efficiencies for all gRNAs 
ϵ_KO: global knockout efficiency; fraction of mutations leading to effective gene knockout
                
***OUTPUT***
γₓ₂: expected coverage of all pairwise combinations of gene knockouts
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
    
    # Apply BioCCP function
    γₓ₂ = expectation_fraction_collected(n_combinations_genes, N; p=p_genes, r=combinations_pp, normalize=false)

    return γₓ₂ 
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

Computes the expected coverage of all triple combinations of gene knockouts (γx,3) for an experiment with given plant library size using BioCCP

***INPUT***
x: number of target genes in the experiment
N: plant library size
g: number of gRNAs designed per target gene
r: number of gRNA sequences per combinatorial gRNA/Cas construct
n_gRNA_total: total number of gRNAs in the experiment
p_gRNA_freq: vector with relative frequencies for all gRNAs in the construct library (normalized!)
p_gRNA_edit: vector with genome editing efficiencies for all gRNAs 
ϵ_KO: global knockout efficiency; fraction of mutations leading to effective gene knockout
                
***OUTPUT***
γₓ₃: expected coverage of all triple combinations of gene knockouts
"""
function BioCCP_γₓ₃(x, 
                N,
                g, 
                r, 
                n_gRNA_total, 
                p_gRNA_freq, 
                p_gRNA_edit, 
                ϵ_KO)
    
    # how many triple combinations of gRNAs
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
    
    # Apply BioCCP function
    γₓ₃ = expectation_fraction_collected(n_combinations_genes, N; p=p_genes, r=combinations_pp, normalize=false)

    return γₓ₃ 
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

Computes the probability of full coverage of all triple combinations of gene knockouts (Px,3) for an experiment with given plant library size using BioCCP

***INPUT***
x: number of target genes in the experiment
N: plant library size
g: number of gRNAs designed per target gene
r: number of gRNA sequences per combinatorial gRNA/Cas construct
n_gRNA_total: total number of gRNAs in the experiment
p_gRNA_freq: vector with relative frequencies for all gRNAs in the construct library (normalized!)
p_gRNA_edit: vector with genome editing efficiencies for all gRNAs 
ϵ_KO: global knockout efficiency; fraction of mutations leading to effective gene knockout
                
***OUTPUT***
Pₓ₃: probability of full coverage of all triple combinations of gene knockouts
"""
function BioCCP_Pₓ₃(x, 
                N,
                g, 
                r, 
                n_gRNA_total, 
                p_gRNA_freq, 
                p_gRNA_edit, 
                ϵ_KO)
    
    # how many triple combinations of gRNAs
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
    
    # Apply BioCCP function
    Pₓ₃ = success_probability(n_combinations_genes, N; p=p_genes, r=combinations_pp, normalize=false)

    return Pₓ₃ 
end
      
      
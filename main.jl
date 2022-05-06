using Random 
using Plots  
using Distributions 
using LinearAlgebra
using Combinatorics
using BioCCP
using ArgParse
using XLSX
using DataFrames
using Weave
using DataStructures
using PrettyTables

include("MultiplexCrisprDOE.jl");

function main(args)

    s = ArgParseSettings("MultiplexCrisprDOE")

    @add_arg_table! s begin
        "gfd" #, "gRNA_freq_dist"
            action = :command        # adds a command which will be read from an argument
            help = "gRNA/Cas9 frequencies"
        "ged" #, "gRNA_edit_dist" 
            action = :command        # adds a command which will be read from an argument
            help = "gRNA/Cas9 editing efficiencies"
        "sim"  # simulation
            action = :command
            help = "simulation-based approaches for computing the minimal plant library size that guarantees full combinatorial coverage (and other relevant statistics)"
        "ccp" # bioccp
            action = :command
            help = "BioCCP-based approaches for computing the minimal plant library size that guarantees full combinatorial coverage (and other relevant statistics)"
    end

    @add_arg_table! s["gfd"] begin    # add command arg_table: same as usual, but invoked on s["grna"]
        "m"
            arg_type = Int
            help = "plant library size"
        "sd"
            arg_type = Int
            help = "the standard deviation on the gRNA abundances (in terms of absolute or relative frequency)"
        "l"
            arg_type = Int
            help = "minimal gRNA abundance (in terms of absolute or relative frequency)"
        "u"
            arg_type = Int
            help = "maximal gRNA abundance (in terms of absolute or relative frequency)"
        "n" #, "--n_gRNA_total"
            arg_type = Int
            help = "the total number of gRNAs in the experiment"
        "--normalize"
            action = :store_true
            # arg_type = Bool
            # default = true
            help = "if set to 'true', the gRNA abundances (absolute frequencies) are converted into relative frequencies (default: true)"
        "--visualize"
            action = :store_true
            # arg_type = Bool
            # default = false
            help = "if set to 'true', a histogram of all gRNA abundances is plotted (default: false)"
    end

    @add_arg_table! s["ged"] begin    # add command arg_table: same as usual, but invoked on s["grna"]
        "f_act"
            arg_type = Float16
            help = "fraction of all gRNAs that is active"
        "eps_edit_act"
            arg_type = Float16
            help = "Average genome editing efficiency for active gRNAs - mean of the genome editing efficiency distribution for active gRNAs"
        "eps_edit_inact"
            arg_type = Float16
            help = "Average genome editing efficiency for inactive gRNAs - mean of the genome editing efficiency distribution for inactive gRNAs"
        "sd_act"
            arg_type = Float16
            help = "standard deviation of the genome editing efficiency distributions for active and inactive gRNAs"
        "n_gRNA_total"
            arg_type = Int
            help = "the total number of gRNAs in the experiment"
        "--visualize"
            action = :store_true
            # arg_type = Bool
            # default = false
            help = "if set to 'true', a histogram of all genome editing efficiency is plotted (default: false)"
    end
    
    @add_arg_table! s["sim"] begin
        "M" #, "--mode"
            # action = :command 
            # dest_name = "M"
            arg_type = Int
            range_tester = x -> 1 <= x <= 4
            help = """Select simulation mode (1: simulate_Nₓ₁; 2: simulate_Nₓ₂; 3: simulate_Nₓ₃; 4: simulate_Nₓ₂_countKOs)"""
        "x"
            arg_type = Int
            help = "number of target genes in the experiment"
        "g"
            arg_type = Int
            help = "number of gRNAs designed per target gene"
        "r"
            arg_type = Int
            help = "number of gRNA sequences per combinatorial gRNA/Cas construct"
        "t"#, "--n_gRNA_total"
            arg_type = Int
            help = "total number of gRNAs in the experiment"
        "f"#, "--p_gRNA_freq"
            arg_type = String #Vector{Float64}
            help = "vector with relative frequencies for all gRNAs in the construct library (normalized!)"
        "e"#, "--p_gRNA_edit"
            arg_type = String #Vector{Float64}
            help = "vector with genome editing efficiencies for all gRNAs"
        "E"#, "--ϵ_KO"
            arg_type=Float16
            help = "global knockout efficiency; fraction of mutations leading to effective gene knockout"
        "--i", "--iter"
            arg_type = Int
            default = 500
            help = "number of CRISPR/Cas experiments that are simulated"
    end

    @add_arg_table! s["ccp"] begin
        "M"#, "--mode"
            arg_type = Int
            range_tester = x -> 1 <= x <= 9
            help = """Select BioCCP mode (1: BioCCP_Nₓ₁; 2: BioCCP_Nₓ₂; 3: BioCCP_Nₓ₃; 4: BioCCP_Pₓ₁; 5: BioCCP_Pₓ₂ ;
            6: BioCCP_Pₓ₃; 7: BioCCP_γₓ₁; 8: BioCCP_γₓ₂; 9: BioCCP_γₓ₃)"""
        "x"
            arg_type = Int
            help = "number of target genes in the experiment"
        "N"
            arg_type = Int
            help = "plant library size"
        "g"
            arg_type = Int
            help = "number of gRNAs designed per target gene"
        "r"
            arg_type = Int
            help = "number of gRNA sequences per combinatorial gRNA/Cas construct"
        "t"#, "--n_gRNA_total"
            arg_type = Int
            help = "total number of gRNAs in the experiment"
        "f"#, "--p_gRNA_freq"
            arg_type = String #Vector{Float64}
            help = "File containing vector with relative frequencies for all gRNAs in the construct library (normalized!)"
        "e"#, "--p_gRNA_edit"
            arg_type = String #Vector{Float64}
            help = "File containing vector with genome editing efficiencies for all gRNAs"
        "E"#, "--ϵ_KO"
            arg_type=Float16
            help = "global knockout efficiency; fraction of mutations leading to effective gene knockout"
    end

    parsed_args = parse_args(args, s)
    command_args = parsed_args[parsed_args["%COMMAND%"]]

    tool_info = OrderedDict()
    args_info = OrderedDict()
    grna_dict = Dict()
    out_dict = Dict()
    if parsed_args["%COMMAND%"] == "gfd"
        tool_info["method"] = "gRNA_frequency_distribution"
        tool_info["description"] = "Generates vector with frequencies in the combinatorial "* 
                                    "gRNA/Cas9 construct library for all gRNAs"
        args_info["Plant library size"] = command_args["m"]
        args_info["SD on the gRNA abundances"] = command_args["sd"]
        args_info["Minimal gRNA abundance"] = command_args["l"]
        args_info["Maximal gRNA abundance"] = command_args["u"]
        args_info["Total number of gRNAs"] = command_args["n"]
        args_info["Convert gRNA abundances to relative frequencies"] = string(command_args["normalize"])
        args_info["Plot gRNA abundances"] = string(command_args["visualize"])

        p_gRNA_reads = gRNA_frequency_distribution(m, sd, l, u, n_gRNA_total; normalize = true, visualize = false)
        grna_dict["p_gRNA_reads"] = p_gRNA_reads

    elseif parsed_args["%COMMAND%"] == "ged"
        tool_info["method"] = "gRNA_edit_distribution"
        tool_info["description"] = "Generates vector with genome editing efficiencies "*
                                    "for all the gRNAs in the experiment"
        args_info["Fraction of active gRNAs"] = command_args["f_act"]
        args_info["Average genome editing efficiency of active gRNAs"] = command_args["eps_edit_act"]
        args_info["Average genome editing efficiency of inactive gRNAs"] = command_args["eps_edit_inact"]
        args_info["Standard deviation"] = command_args["sd_act"]
        args_info["Total number of gRNAs"] = command_args["n_gRNA_total"]
        args_info["Plot genome editing efficiency"] = string(command_args["visualize"])

        p_gRNA_edit = gRNA_frequency_distribution(m, sd, l, u, n_gRNA_total; normalize = true, visualize = false)
        grna_dict["p_gRNA_edit"] = p_gRNA_edit

    elseif parsed_args["%COMMAND%"] == "sim" || parsed_args["%COMMAND%"] == "ccp"

        filename = command_args["f"]
        sheet = 1
        data = DataFrame(XLSX.readtable(filename, 1)...)
        p_gRNA_reads = data[!,"gRNA_read"]
        p_gRNA_reads_normalized = p_gRNA_reads/sum(p_gRNA_reads)  # normalize
        f = p_gRNA_reads_normalized
        grna_dict["p_gRNA_reads"] = f

        filename = command_args["e"]
        sheet = 1
        data = DataFrame(XLSX.readtable(filename, 1)...)
        p_gRNA_edit = data[!,"gRNA_edit_efficiency"]
        e = p_gRNA_edit
        grna_dict["p_gRNA_edit"] = e

        x = command_args["x"]
        g = command_args["g"]
        r = command_args["r"]
        t = command_args["t"] # n_gRNA_total
        E = command_args["E"] # ϵ_KO # iter = 500
        
        args_info["# of target genes in the experiment"] = command_args["x"]
        args_info["# of gRNAs designed per target gene"] = command_args["g"]
        args_info["# of gRNAs / combi gRNA/Cas construct"] = command_args["r"]
        args_info["Total number of gRNAs"] = command_args["t"]
        args_info["Relative frequencies for all gRNAs"] = command_args["f"]
        args_info["Genome editing efficiencies for all gRNAs"] = command_args["e"]
        args_info["Global knockout efficiency"] = command_args["E"]

        if parsed_args["%COMMAND%"] == "sim"
            tool_info["method"] = "simulation"
            tool_info["description"] = "simulation-based approaches for computing the minimal "* 
                                        "plant library size that guarantees full combinatorial "*
                                        "coverage (and other relevant statistics)"
            i = command_args["i"] # iter = 500
            args_info["# of simulated experiments"] = command_args["i"]
            
            if command_args["M"] == 1
                tool_info["mode"] = "simulate_Nx1"
                tool_info["mode_description"] = "Computes the expected value and the standard deviation "* 
                                                "of the minimal plant library size for full coverage of "*
                                                "all single gene knockouts (E[Nx,1] and σ[Nx,1]) using simulation"
                E_sim, sd_sim = simulate_Nₓ₁(x, g, r, t, f, e, E; iter=i)
                out_dict["E_sim"] = E_sim
                out_dict["sd_sim"] = sd_sim

            elseif command_args["M"] == 2
                tool_info["mode"] = "simulate_Nx2"
                tool_info["mode_description"] = "Computes the expected value and the standard deviation of "* 
                    "the minimal plant library size for full coverage of "*
                    "all pairwise combinations of gene knockouts in a "*
                    "multiplex CRISPR/Cas experiment (E[Nx,2] and σ[Nx,2]) using simulation"

                E_sim, sd_sim = simulate_Nₓ₂(x, g, r, t, f, e, E; iter=i)
                out_dict["E_sim"] = E_sim
                out_dict["sd_sim"] = sd_sim

            elseif command_args["M"] == 3
                tool_info["mode"] = "simulate_Nx3"
                tool_info["mode_description"] = "Computes the expected value and the standard deviation of "*
                    "the minimal plant library size for full coverage of "*
                    "all triple combinations of gene knockouts in a "*
                    "multiplex CRISPR/Cas experiment (E[Nx,3] and σ[Nx,3]) using simulation"

                E_sim, sd_sim = simulate_Nₓ₃(x, g, r, t, f, e, E; iter=i)
                out_dict["E_sim"] = E_sim
                out_dict["sd_sim"] = sd_sim

            elseif command_args["M"] == 4
                tool_info["mode"] = "simulate_Nx2_countKOs"
                tool_info["mode_description"] = "Counts of the number of knockouts per plant in the experiment"

                n_KOs_vec = simulate_Nₓ₂_countKOs(x, g, r, t, f, e, E; iter=i)
                out_dict["n_KOs_vec"] = n_KOs_vec
            
            end

        elseif parsed_args["%COMMAND%"] == "ccp"
            method = "gRNA_frequency_distribution"
            description = "Generates vector with frequencies in the combinatorial "*
                                "gRNA/Cas9 construct library for all gRNAs"
            N = command_args["N"]
            args_info["Plant library size"] = command_args["N"]

            if command_args["M"] == 1
                tool_info["mode"] = "BioCCP_Nx1"
                tool_info["mode_description"] = "Computes the expected value and the standard deviation of "*
                    "the minimal plant library size for full coverage of all "*
                    "single gene knockouts (E[Nx,1] and σ[Nx,1]) using BioCCP"

                E_sim, sd_sim = BioCCP_Nₓ₁(x, g, r, t, f, e, E)
                out_dict["E_sim"] = E_sim
                out_dict["sd_sim"] = sd_sim

            elseif command_args["M"] == 2
                tool_info["mode"] = "BioCCP_Nx2"
                tool_info["mode_description"] = "Computes the expected value and the standard deviation of "*
                    "the minimal plant library size for full coverage of all "*
                    "pairwise combinations of gene knockouts in a multiplex "*
                    "CRISPR/Cas experiment (E[Nx,2] and σ[Nx,2]) using BioCCP"

                E_sim, sd_sim = BioCCP_Nₓ₂(x, g, r, t, f, e, E)
                out_dict["E_sim"] = E_sim
                out_dict["sd_sim"] = sd_sim

            elseif command_args["M"] == 3
                tool_info["mode"] = "BioCCP_Nx3"
                tool_info["mode_description"] = "Computes the expected value and the standard deviation of "*
                    "the minimal plant library size for full coverage of all triple combinations of "*
                    "gene knockouts in a multiplex CRISPR/Cas experiment (E[Nx,3] and σ[Nx,3]) using BioCCP"

                E_sim, sd_sim = BioCCP_Nₓ₃(x, g, r, t, f, e, E)
                out_dict["E_sim"] = E_sim
                out_dict["sd_sim"] = sd_sim

            elseif command_args["M"] == 4
                tool_info["mode"] = "BioCCP_Px1"
                tool_info["mode_description"] = "Computes the probability of full coverage of "* 
                    "all single gene knockouts (Px,1) for an experiment with given "*
                    "plant library size using BioCCP"

                E_sim, sd_sim = BioCCP_Pₓ₁(x, N, g, r, t, f, e, E)
                out_dict["E_sim"] = E_sim
                out_dict["sd_sim"] = sd_sim

            elseif command_args["M"] == 5
                tool_info["mode"] = "BioCCP_Px2"
                tool_info["mode_description"] = "Computes the probability of full coverage of "*
                    "all pairwise combinations of gene knockouts (Px,2) "* 
                    "for an experiment with given plant library size using BioCCP"

                E_sim, sd_sim = BioCCP_Pₓ₂(x, N, g, r, t, f, e, E)
                out_dict["E_sim"] = E_sim
                out_dict["sd_sim"] = sd_sim

            elseif command_args["M"] == 6
                tool_info["mode"] = "BioCCP_Px3"
                tool_info["mode_description"] = "Computes the probability of full coverage of all "*
                    "triple combinations of gene knockouts (Px,3) for an experiment "*
                    "with given plant library size using BioCCP"

                E_sim, sd_sim = BioCCP_Pₓ₃(x, N, g, r, t, f, e, E)
                out_dict["E_sim"] = E_sim
                out_dict["sd_sim"] = sd_sim

            elseif command_args["M"] == 7
                tool_info["mode"] = "BioCCP_γx1"
                tool_info["mode_description"] = "Computes the expected coverage of all "*
                    "single gene knockouts (E[γx,1]) for an experiment "*
                    "with given plant library size using BioCCP"

                E_sim, sd_sim = BioCCP_γₓ₁(x, N, g, r, t, f, e, E)
                out_dict["E_sim"] = E_sim
                out_dict["sd_sim"] = sd_sim

            elseif command_args["M"] == 8
                tool_info["mode"] = "BioCCP_γx2"
                tool_info["mode_description"] = "Computes the expected coverage of all "*
                    "pairwise combinations of gene knockouts (E[γx,2]) for an experiment with "*
                    "given plant library size using BioCCP"

                E_sim, sd_sim = BioCCP_γₓ₂(x, N, g, r, t, f, e, E)
                out_dict["E_sim"] = E_sim
                out_dict["sd_sim"] = sd_sim

            elseif command_args["M"] == 9
                tool_info["mode"] = "BioCCP_γx3"
                tool_info["mode_description"] = "Computes the expected coverage of all "*
                    "triple combinations of gene knockouts (E[γx,3]) for an experiment with "*
                    "given plant library size using BioCCP"

                E_sim, sd_sim = BioCCP_γₓ₃(x, N, g, r, t, f, e, E)
                out_dict["E_sim"] = E_sim
                out_dict["sd_sim"] = sd_sim
            end
        end
    end

    println(parsed_args)
    println("Parsed args:")
    for (key,val) in parsed_args
        println("  $key  =>  $(repr(val))")
    end
    println()
    println("Command: ", parsed_args["%COMMAND%"])
    h1 = histogram(grna_dict["p_gRNA_reads"], label="", 
                    xlabel="Number of reads per gRNA", 
                    linecolor="white", 
                    normalize=:probability,
                    xtickfontsize=10,ytickfontsize=10,
                    color=:mediumturquoise, size=(600,350), bins = 25,
                    ylabel="Relative frequency", 
                    title="gRNA frequency distribution")
    
    h2 = histogram(grna_dict["p_gRNA_edit"], 
                    normalize = :probability,
                    linecolor = "white",
                    label="", 
                    color=:turquoise4,
                    xtickfontsize=10,ytickfontsize=10, xlim = (0, 1),
                    xticks=(0:0.1:1),
                    bins = 150,
                    xlabel="gRNA editing efficiency", 
                    ylabel="Relative frequency", 
                    title="gRNA genome editing effiency distribution")


    weave("report.jmd", 
            args = (parsed_args = parsed_args,
                    tool_info = tool_info,
                    args_info = args_info,
                    grna_dict = grna_dict,
                    h1 = h1, h2 = h2,
                    output = out_dict); 
            doctype = "md2html", out_path = :pwd)
end

main(ARGS)

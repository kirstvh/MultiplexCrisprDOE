# MultiplexCrisprDOE
This repository accompanies the paper "Covering the Combinatorial Design Space of Multiplex CRISPR/Cas Experiments in Plants" by Kirsten Van Huffel, Michiel Stock, Thomas Jacobs, Tom Ruttink and Bernard De Baets.

It provides simulation- and BioCCP-based approaches for calculating the plant library size for full combinatorial coverage of multiplex CRISPR/Cas experiments in plants.

## Content
All functions belonging to the simulation- and BioCCP-based approaches can be found in the file `MultiplexCrisprDOE.jl`. A short description is provided in the Table below.

Function name    | Short description
---------------- | -----------------
`gRNA_frequency_distribution`        | Generates vector with frequencies in the combinatorial gRNA/Cas9 construct library for all gRNAs 
`gRNA_edit_distribution`      | Generates vector with genome editing efficiencies for all the gRNAs in the experiment 
`simulate_Nₓ₁`         | Computes E[N<sub>x,1</sub>] and σ[N<sub>x,1</sub>] of a CRISPR/Cas experiment, studying the minimal plant library size for full coverage of all single gene knockouts, using the simulation-based approach
`BioCCP_Nₓ₁` | Computes E[N<sub>x,1</sub>] and σ[N<sub>x,1</sub>] of a CRISPR/Cas experiment, studying the minimal plant library size for full coverage of all single gene knockouts, using the BioCCP-based approach
`simulate_Nₓ₂`      | Computes  E[N<sub>x,2</sub>] and σ[N<sub>x,2</sub>] of multiplex CRISPR/Cas experiment, studying the minimal plant library size for full coverage of all pairwise combinations of gene knockouts, using the simulation-based approach
`BioCCP_Nₓ₂`         | Computes  E[N<sub>x,2</sub>] and σ[N<sub>x,2</sub>] of a multiplex CRISPR/Cas experiment, studying the minimal plant library size for full coverage of all pairwise combinations of gene knockouts, using the BioCCP-based approach
`simulate_Nₓ₂_countKOs` | Counts of the number of knockouts per plant in the experiment 
`simulate_Nₓ₃` | Computes  E[N<sub>x,3</sub>] and σ[N<sub>x,3</sub>] of a multiplex CRISPR/Cas experiment, studying the minimal plant library size for full coverage of all pairwise combinations of gene knockouts, using the simulation-based approach
`BioCCP_Nₓ₃` | Computes  E[N<sub>x,3</sub>] and σ[N<sub>x,3</sub>] of a multiplex CRISPR/Cas experiment, studying the minimal plant library size for full coverage of all pairwise combinations of gene knockouts, using the BioCCP-based approach

The default values for the experimental design parameters used in this work can be found under `DefaultParameters_k=1.jl`, `DefaultParameters_k=2.jl` and `DefaultParameters_k=3.jl`.

## Reproducibility
The graphs and results from our work can be reproduced by running the Jupyter notebooks. 

Values of design parameters can be adjusted so that researchers can do computations for their own multiplex CRISPR/Cas experiments.

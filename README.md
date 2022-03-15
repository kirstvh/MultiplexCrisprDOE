# Design of Multiplex CRISPR/Cas Experiments in Plants (MultiplexCrisprDOE)
This repository accompanies the paper "Covering the Combinatorial Design Space of Multiplex CRISPR/Cas Experiments in Plants" by Kirsten Van Huffel, Michiel Stock, Thomas Jacobs, Tom Ruttink and Bernard De Baets.

It provides simulation- and BioCCP-based approaches for computing the minimal plant library size that guarantees full combinatorial coverage (and other relevant statistics) for  multiplex CRISPR/Cas experiments in plants. 

## Content
All functions belonging to the simulation- and BioCCP-based approaches are documented in the file `MultiplexCrisprDOE.jl`. A short description of all functions is provided in the Table below.

Function name    | Short description
---------------- | -----------------
`gRNA_frequency_distribution`        | Generates vector with frequencies in the combinatorial gRNA/Cas9 construct library for all gRNAs 
`gRNA_edit_distribution`      | Generates vector with genome editing efficiencies for all the gRNAs in the experiment 
`simulate_Nₓ₁`         | Computes the expected value and the standard deviation of the minimal plant library size for full coverage of all single gene knockouts (E[N<sub>x,1</sub>] and σ[N<sub>x,1</sub>]) using simulation 
`BioCCP_Nₓ₁` | Computes the expected value and the standard deviation of the minimal plant library size for full coverage of all single gene knockouts (E[N<sub>x,1</sub>] and σ[N<sub>x,1</sub>]) using BioCCP 
`BioCCP_Pₓ₁` | Computes the probability of full coverage of all single gene knockouts (P<sub>x,1</sub>) for an experiment with given plant library size using BioCCP 
`BioCCP_γₓ₁` | Computes the expected coverage of all single gene knockouts (E[γ<sub>x,1</sub>]) for an experiment with given plant library size using BioCCP 
`simulate_Nₓ₂`      | Computes the expected value and the standard deviation of the minimal plant library size for full coverage of all pairwise combinations of gene knockouts in a multiplex CRISPR/Cas experiment (E[N<sub>x,2</sub>] and σ[N<sub>x,2</sub>]) using simulation 
`BioCCP_Nₓ₂`         | Computes the expected value and the standard deviation of the minimal plant library size for full coverage of all pairwise combinations of gene knockouts in a multiplex CRISPR/Cas experiment (E[N<sub>x,2</sub>] and σ[N<sub>x,2</sub>]) using BioCCP 
`simulate_Nₓ₂_countKOs` | Counts of the number of knockouts per plant in the experiment 
`BioCCP_Pₓ₂` | Computes the probability of full coverage of all pairwise combinations of gene knockouts (P<sub>x,2</sub>) for an experiment with given plant library size using BioCCP 
`BioCCP_γₓ₂` |  Computes the expected coverage of all pairwise combinations of gene knockouts (E[γ<sub>x,2</sub>]) for an experiment with given plant library size using BioCCP 
`simulate_Nₓ₃` | Computes the expected value and the standard deviation of the minimal plant library size for full coverage of all triple combinations of gene knockouts in a multiplex CRISPR/Cas experiment (E[N<sub>x,3</sub>] and σ[N<sub>x,3</sub>]) using simulation 
`BioCCP_Nₓ₃` | Computes the expected value and the standard deviation of the minimal plant library size for full coverage of all triple combinations of gene knockouts in a multiplex CRISPR/Cas experiment (E[N<sub>x,3</sub>] and σ[N<sub>x,3</sub>]) using BioCCP 
`BioCCP_Pₓ₃` | Computes the probability of full coverage of all triple combinations of gene knockouts (P<sub>x,3</sub>) for an experiment with given plant library size using BioCCP 
`BioCCP_γₓ₃` | Computes the expected coverage of all triple combinations of gene knockouts (E[γ<sub>x,3</sub>]) for an experiment with given plant library size using BioCCP 

The default values for the experimental design parameters used in this work can be found under `DefaultParameters_k=1.jl`, `DefaultParameters_k=2.jl` and `DefaultParameters_k=3.jl`.

## Reproducibility
The graphs and results from our work can be reproduced by running the Jupyter notebooks. 

Values of design parameters can be adjusted so that researchers can do computations for their own multiplex CRISPR/Cas experiments.

## Animation
Consider a multiplex CRISPR/Cas experiment targeting pairwise combinations of gene knockouts, characterized by the experimental design parameters listed in the file `DefaultParameters_k=2.jl`.
The animation below illustrates the occurence of pairwise combinations of gene knockouts when including an increasing number of plants in a plant library. 

<img src="https://github.com/kirstvh/MultiplexCrisprDOE/blob/main/sampling_plants.gif" />



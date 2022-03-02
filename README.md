# MultiplexCrisprDOE
This repository accompanies the paper "Covering the Combinatorial Design Space of Multiplex CRISPR/Cas Experiments in Plants" by Kirsten Van Huffel, Michiel Stock, Thomas Jacobs, Tom Ruttink and Bernard De Baets.

It provides simulation- and BioCCP-based approaches for calculating the plant library size for full combinatorial coverage of multiplex CRISPR/Cas experiments in plants.

## Content
All functions belonging to the simulation- and BioCCP-based approaches can be found in the file `MultiplexCrisprDOE.jl`. An short description is provided in the Table below.

Function name    | Short description
---------------- | -----------------
`gRNA_frequency_distribution`        | description  
`gRNA_activity_distribution`      | description 
`simulate_Nₓ₁`         | description
`BioCCP_Nₓ₁` | description
`simulate_Nₓ₂`      | description
`BioCCP_Nₓ₂`         | description
`simulate_Nₓ₂_countKOs` | description
`simulate_Nₓ₃` | description
`BioCCP_Nₓ₃` | description

## Reproducibility
The graphs and results from our work can be regenerated by running the Jupyter notebooks. Values of design parameters can be adjusted so that researchers can do computations for their own multiplex CRISPR/Cas experiments.

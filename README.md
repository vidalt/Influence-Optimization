# Influence-Optimization
Code and data from the paper "Influence Optimization in Networks: New Formulations and Valid Inequalities", authored by V. Ferreira, A. Pessoa, and T. Vidal. This paper is under peer-review but its preprint is available on https://arxiv.org/abs/2209.13065.

## Requirements
- Julia 1.0.5
- JuMP 0.18
- CPLEX 12.9

> Disclaimer: this project does not use the latest versions of JuMP (and therefore CPLEX)

## Installation

Set CPLEX path:
```bash
export CPLEX_STUDIO_BINARIES="/path/to/cplex/bin/x86-64_*"
julia
```

Install the dependencies:
```julia
julia> ]
(v1.1) pkg> activate .
(InfluenceOptimization) pkg> instantiate
```

## Run

```bash
$ julia --project=. src/main.jl
usage: main.jl [-o OUTPUT] [-v] [-b] [-l] [-u UPCUTOFF] [-a ALPHA]
               [-g GAMMA] [-h] filepath {icc|icc+|licc+|cf}

commands:
  icc                   ICC
  icc+                  ICC+
  licc+                 LICC+
  cf                    CF

positional arguments:
  filepath              Instance file path or path directory (default:
                                                "/data/socnet-instances-v2/SW-n50-k4-b0.1-d1-10-g0.7-i3")

optional arguments:
  -o, --output OUTPUT   Output filename
  -v, --verbose         Verbose output
  -l, --latex           Latex Output
  -u, --upcutoff UPCUTOFF
                        Up cutoff (type: Float64, default: 1.0e7)
  -a, --alpha ALPHA     α|V| nodes that need to be activated (type:
                        Float64, default: 0.1)
  -g, --gamma GAMMA     Activation function Γ: 0.9, 1.0 or 1.1 (type:
                        Float64, default: 1.0)
  -h, --help            show this help message and exit
```

Usage example:
```bash
julia --project=. main.jl data/socnet-instances-v2/SW-n50-k8-b0.3-d1-10-g0.7-i1 -a 0.1 -u 45 -g 1.0 cf
```

To run the application inside the Julia environment, type:
```julia
include("src/main.jl")
```

# Reported results
Detailed results per instance and for each configuration of $n, k, \beta, \alpha, \Gamma$ can be seen in `results/final_results.xlsx`.

# References

```
@misc{https://doi.org/10.48550/arxiv.2209.13065,
  doi = {10.48550/ARXIV.2209.13065},
  url = {https://arxiv.org/abs/2209.13065},
  author = {Ferreira, Vinicius and Pessoa, Artur and Vidal, Thibaut},
  keywords = {Optimization and Control (math.OC), FOS: Mathematics, FOS: Mathematics},
  title = {Influence Optimization in Networks: New Formulations and Valid Inequalities},
  publisher = {arXiv},
  year = {2022},
  copyright = {arXiv.org perpetual, non-exclusive license}
}
```
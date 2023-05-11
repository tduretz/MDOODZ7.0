# [Julia](https://julialang.org) visualisation of MDoodz outputs using [Makie](https://docs.makie.org/stable/)

## Installation
1. Clone/Download the repository
2. Run Julia from within the folder 
3. In Julia's REPL switch to package mode: type `]`
4. Activate the environement: type `activate .`
5. Install all necessary dependencies: type `instantiate`

## Main visualisation 

There main file [Main_Visualisation_Makie_MD7.jl](./Main_Visualisation_Makie_MD7.jl) that contained several visualisation options.
Set the `path` to your file, e. g.:
```julia
  path ="/Users/imac/REPO_GIT/MDOODZ7.0/MDLIB/"
```
Select the time steps to be visualised, for example to visualise from step 0 to 100 with a step of 20:
```julia
  file_start = 0
  file_step  = 20
  file_end   = 100
```
Define which field, for example the phases:
```julia
  field = :phases
```

... customize and contribute!

## Visual tests

Results obtained with new code versions can be compared to reference models. Here is an example for `PinchSwellGSE` which models pinch-and-swell formation with grain size evolution:

<img src="https://github.com/tduretz/MDOODZ7.0/blob/update_aniso/JuliaVisualisation/_VisualTests/PinchSwellGSE.png" alt="PinchSwellGSE.png" width="600">

Another example with the `Shrinking` model configuration which models the stress perturbation due to densification:
![./_VisualTests/Shrinking.png](https://github.com/tduretz/MDOODZ7.0/blob/update_aniso/JuliaVisualisation/_VisualTests/Shrinking.png)

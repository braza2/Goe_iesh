using DrWatson
@quickactivate "Goe_iesh"

cd("./src/")
#include("Goe_iesh.jl")



curr_vers = "v.0.95"
using DrWatson
using LinearAlgebra
using Random
using StaticArrays
#using Plots
#using HDF5
using Statistics
using Dates
using Base.Threads


include("read_au_data.jl"); #executes
include("setup_parameters.jl");  #executes
include("types.jl");  #executes
include("misc_func.jl");
include("constructor.jl");
include("get_neighbors.jl");
include("get_energies.jl");
include("get_forces.jl");
include("propagation.jl");
include("run_trajectory.jl");

s = simulation_init()
simulation_constructor_x_v!(s)
simulation_constructor_nn!(s)
#simulation_constructor_x_300K!(s)
s.v[3:,:] .= 0.0
simulation_constructor_energy(s)
simulation_constructor_force(s)
propagate_init!(s)
#run simulation
simulate!(s)
module Goe_iesh


curr_vers = "v.0.95"
using DrWatson
using LinearAlgebra
using Random
using StaticArrays
using Plots
using HDF5
using Statistics
using Dates
using Base.Threads
using FastGaussQuadrature
using PolyChaos

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
simulation_constructor_x_300K!(s)
#s.v[3:end,:] .= 0.0
simulation_constructor_energy(s)
simulation_constructor_force(s)
propagate_init!(s)
#run simulation
simulate!(s)


multiple_trajectory()

end


sigma = 1.0*ev_kjmol *kjmol_seunit
sigma2 = 1.0*ev_kjmol *kjmol_seunit
em = -delta_E/2.0
ep = delta_E/2.0
supp= (em, ep)
w(t) = abs(t) <= delta_E/2.0 ? (exp(-(t-em)^2.0 / (2.0*sigma2^2))+exp(-(t)^2.0 / (2.0*sigma^2))+exp(-(t-ep)^2.0 / (2.0*sigma2^2)))^2 : error("not allowed")
my_meas = Measure("my_meas", w, supp, false, Dict())
my_op = OrthoPoly("my_op", 40, my_meas; Nquad=1000);
nodes, weight = gauss(my_op)



sigma = 1.0*ev_kjmol *kjmol_seunit
sigma2 = 1.0*ev_kjmol *kjmol_seunit
em = -delta_E/2.0
ep = delta_E/2.0
supp= (em, ep)
w(t) = abs(t) <= delta_E/2.0 ? 1000.0^2*(exp(-(t-em)^2.0 / (2.0*sigma2^2))+exp(-(t)^2.0 / (2.0*sigma^2))+exp(-(t-ep)^2.0 / (2.0*sigma2^2)))^2 : error("not allowed")
my_meas = Measure("my_meas", w, supp, false, Dict())
my_op = OrthoPoly("my_op", 40, my_meas; Nquad=1000);
nodes2, weight2 = gauss(my_op)







plot(nodes, weight*1000.0^2, seriestype = :scatter )
plot!(nodes2, weight2, seriestype = :scatter )

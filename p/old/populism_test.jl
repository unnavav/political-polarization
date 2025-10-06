#                           
#              (            
#    (     (   )\ (      )  
#    )\   ))\ ((_))\  ( /(  
#   ((_) /((_) _ ((_) )(_)) 
#     ! (_))( | | (_)((_)_  
#    | || || || | | |/ _` | 
#   _/ | \_,_||_| |_|\__,_| 
#  |__/                     
# 					by vaas
#
#	JAN 2025
# This file is part of the Political Polarization project. Script to import steady states and then
# transition between them. 


# ----------------------------------------------------
# Importing required packages
# import Pkg; Pkg.add("MAT"); Pkg.add("CSV");
using MAT, CSV, Tables

# importing modules from the project
include("./dependencies/predict_Mod.jl")
using .predict

# Loading liberalist household data
liberalism_populism_dir = "c:/Users/unnav/Dropbox/Education/OSU/Ongoing_Research/2nd Year Paper/political-polarization/d/liberalism_populism/"

cd(liberalism_populism_dir)

cd("results_t0.0000_eta1.0000")
VP = CSV.File("V.csv") |> Tables.matrix
GP = CSV.File("G.csv") |> Tables.matrix
adistrP = CSV.File("adistr.csv") |> Tables.matrix
rP = CSV.File("r.csv") |> Tables.matrix	
lP = CSV.File("l.csv") |> Tables.matrix

cd("../results_t0.0000_eta1.1000")
VL = CSV.File("V.csv") |> Tables.matrix
GL = CSV.File("G.csv") |> Tables.matrix
adistrL = CSV.File("adistr.csv") |> Tables.matrix
rL = CSV.File("r.csv") |> Tables.matrix
lL = CSV.File("l.csv") |> Tables.matrix

# Now to do perfect foresight transitions --------------------------------------------
# Give a series of interest rates that switch partway through

# lagg will change; first period and then the growth rate changes,
# which means new prices for everything
lt = repeat(lL, 500)

r0 = rP
r1 = rL
println("r0: ", r0)
println("r1: ", r1)

alpha = 0.36
delta = 0.06

terms = Dict("alpha" => alpha, "delta" => delta)

lambda = 0.5
dTol = 1e-4

rt = predict.transition(r0, r1, lt, terms, dTol, lambda)

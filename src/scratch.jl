using Revise
Pkg.activate(".")

using PCAmix

using DataFrames, TypedTables, RCall, LinearAlgebra, Statistics
D = Table(x1 = randn(1000), f1 = categorical(rand(["a", "b", "c"], 1000)), f2 = categorical(rand(["TRUE", "FALSE"], 1000)))

@btime U, Δ, Vt = pcamix(D);






import Juno.@enter
using BenchmarkTools

X ≈ U*Δ*V'

R"library(PCAmixdata)"
DF = DataFrame(D)
@btime R"PCAmix($DF[, 1, drop = FALSE], $DF[, 2:3], graph=FALSE)$eig"

R"pca = prcomp($X, scale = FALSE, center = FALSE)"

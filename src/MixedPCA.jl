module PCAmix
export pcamix
    using LinearAlgebra, TypedTables, Statistics, CategoricalArrays
    import StatsBase: countmap
    include("pca.jl")
end

function gsvd(Z::AbstractMatrix, N, M)
    Z̃ = N^(1/2)*Z*M^(1/2)
    r = svd(Z̃)
    Δ = Diagonal(r.S)
    U = N^(-1/2)*r.U
    Vt = r.Vt*M^(-1/2)
    return U, Δ, Vt
end

function centered_model_matrix(c::AbstractCategoricalArray)
    m = Matrix{Float64}(undef, length(c), length(c.pool))
    for j in 1:size(m, 2)
        l = (j .== c.refs)
        m[:, j] .= l .- mean(l)
    end
    return m
end

function pcamix(D::Table)
    data = []
    m = Float64[]
    max_dim = 0
    n = size(D, 1)
    for col in columns(D)
        if col isa Vector{T} where T<:Number
            max_dim += 1
            col .= (col .- mean(col))./sqrt(var(col, corrected = false))
            append!(m, 1.)
            push!(data, col)
        elseif col isa AbstractCategoricalVector
            ns = countmap(col.refs)
            max_dim += length(ns) - 1
            ns = [length(col)/ns[i] for i in 1:length(ns)]
            append!(m, ns)
            push!(data, centered_model_matrix(col))
        else
            @error "Columns must be numeric or categorical"
        end
    end
    Z = reduce(hcat, data)
    N = Diagonal(fill(1/n, n))
    M = Diagonal(m)
    U, Δ, Vt = gsvd(Z, N, M)
    return U[:, 1:max_dim], Δ[1:max_dim, 1:max_dim], Vt[1:max_dim, :]
end

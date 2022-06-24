using Distributions
using Graphs
using Base.Iterators

# Convenience functions from G.5
## Dictionares and named tupes functionalities
Base.Dict{Symbol,V}(a::NamedTuple) where V =
    Dict{Symbol,V}(n=>v for (n,v) in zip(keys(a), values(a)))

Base.convert(::Type{Dict{Symbol,V}}, a::NamedTuple) where V =
    Dict{Symbol,V}(a)

Base.isequal(a::Dict{Symbol,<:Any}, nt::NamedTuple) =
    length(a) == length(nt) &&
    all(a[n] == v for (n,v) in zip(keys(nt), values(nt)))

struct SetCategorical{S}
    elements::Vector{S} # Set elements (could be repeated)
    distr::Categorical # Categorical distribution over set elements
    function SetCategorical(elements::AbstractVector{S}) where S
        weights = ones(length(elements))
        return new{S}(elements, Categorical(normalize(weights, 1)))
    end
    function SetCategorical(
            elements::AbstractVector{S},
            weights::AbstractVector{Float64}
        ) where S
        l1 = norm(weights,1)
        if l1 < 1e-6 || isinf(l1)
            return SetCategorical(elements)
        end
        distr = Categorical(normalize(weights, 1))
        return new{S}(elements, distr)
    end
end

# Fix due to change in Graphs.jl
function topological_sort(bn)
    topological_sort_by_dfs(bn)
end



# Pg 26
struct Variable
    name::Symbol
    r::Int # number of possible values
end

const Assignment = Dict{Symbol,Int}

const FactorTable = Dict{Assignment,Float64}

struct Factor
    vars::Vector{Variable}
    table::FactorTable
end

variablenames(φ::Factor) = [var.name for var in φ.vars]

select(a::Assignment, varnames::Vector{Symbol}) =
    Assignment(n=>a[n] for n in varnames)

function assignments(vars::AbstractVector{Variable})
    names = [var.name for var in vars]
    return vec([Assignment(n=>v for (n,v) in zip(names, values))
                for values in product((1:v.r for v in vars)...)])
end

function normalize!(φ::Factor)
    z = sum(p for (a,p) in φ.table)
    for (a,p) in φ.table
        φ.table[a] = p/z
    end
    return φ 
end


# Pg 33

struct BayesianNetwork
    vars::Vector{Variable}
    factors::Vector{Factor}
    graph::SimpleDiGraph{Int64} # Directed graph 
end

function probability(bn::BayesianNetwork, assignment)
    subassignment(φ) = select(assignment, variablenames(φ))
    probability(φ) = get(φ.table, subassignment(φ), 0.0) # get(collection, key, default)
    return prod(probability(φ) for φ in bn.factors)
end


# Pg 45

function Base.:*(φ::Factor, ψ::Factor)
    φnames = variablenames(φ)
    ψnames = variablenames(ψ)
    ψonly = setdiff(ψ.vars, φ.vars)
    table = FactorTable()
    for (φa,φp) in φ.table
        for a in assignments(ψonly)
            a = merge(φa, a)
            ψa = select(a, ψnames)
            table[a] = φp * get(ψ.table, ψa, 0.0)
        end 
    end
    vars = vcat(φ.vars, ψonly)
    return Factor(vars, table)
end

function marginalize(φ::Factor, name)
    table = FactorTable()
    for (a, p) in φ.table
        a′ = delete!(copy(a), name)
        table[a′] = get(table, a′, 0.0) + p
    end
    vars = filter(v -> v.name != name, φ.vars)
    return Factor(vars, table)
end

in_scope(name, φ) = any(name == v.name for v in φ.vars)

function condition(φ::Factor, name, value)
    if !in_scope(name, φ)
        return φ 
    end
    table = FactorTable()
    for (a, p) in φ.table
        if a[name] == value
            table[delete!(copy(a), name)] = p
        end 
    end
    vars = filter(v -> v.name != name, φ.vars)
    return Factor(vars, table)
end

function condition(φ::Factor, evidence)
    for (name, value) in pairs(evidence)
        φ = condition(φ, name, value)
    end
    return φ 
end

# Pg 48

struct ExactInference end

function infer(M::ExactInference, bn, query, evidence)
    φ = prod(bn.factors)
    φ = condition(φ, evidence)
    for name in setdiff(variablenames(φ), query)
        φ = marginalize(φ, name)
    end
    return normalize!(φ)
end

# Pg 52

struct VariableElimination
    ordering # array of variable indices
end

function infer(M::VariableElimination, bn, query, evidence)
    Φ = [condition(φ, evidence) for φ in bn.factors]
    for i in M.ordering
        name = bn.vars[i].name
        if name ∉ query
            inds = findall(φ->in_scope(name, φ), Φ)
            if !isempty(inds)
                φ = prod(Φ[inds])
                deleteat!(Φ, inds)
                φ = marginalize(φ, name)
                push!(Φ, φ)
            end         
        end
    end
    return normalize!(prod(Φ))
end


#Pg 55

function Base.rand(φ::Factor)
    tot, p, w = 0.0, rand(), sum(values(φ.table))
    for (a,v) in φ.table
        tot += v/w
        if tot >= p
            return a 
        end
    end
    return Assignment()
end

function Base.rand(bn::BayesianNetwork)
    a = Assignment()
    for i in topological_sort(bn.graph)
        name, φ = bn.vars[i].name, bn.factors[i]
        a[name] = rand(condition(φ, a))[name]
    end
    return a 
end


struct DirectSampling
    m # number of samples
end


# Pg 56

function infer(M::DirectSampling, bn, query, evidence)
    table = FactorTable()
    for i in 1:(M.m)
        a = rand(bn)
        if all(a[k] == v for (k,v) in pairs(evidence))
            b = select(a, query)
            table[b] = get(table, b, 0) + 1
        end
    end
    vars = filter(v->v.name ∈ query, bn.vars)
    return normalize!(Factor(vars, table))
end


# Pg 58

 struct LikelihoodWeightedSampling
    m # number of samples
end

function infer(M::LikelihoodWeightedSampling, bn, query, evidence)
    table = FactorTable()
    ordering = topological_sort(bn.graph)
    for i in 1:(M.m)
        a, w = Assignment(), 1.0
        for j in ordering
            name, φ = bn.vars[j].name, bn.factors[j]
            if haskey(evidence, name)
                a[name] = evidence[name]
                w *= φ.table[select(a, variablenames(φ))]
            else
                a[name] = rand(condition(φ, a))[name]
            end
        end
        b = select(a, query)
        table[b] = get(table, b, 0) + w
    end
    vars = filter(v->v.name ∈ query, bn.vars)
    return normalize!(Factor(vars, table))
end


# Pg 61

function blanket(bn, a, i)
    name = bn.vars[i].name
    val = a[name]
    a = delete!(copy(a), name)
    Φ = filter(φ -> in_scope(name, φ), bn.factors)
    φ = prod(condition(φ, a) for φ in Φ)
    return normalize!(φ)
end

function update_gibbs_sample!(a, bn, evidence, ordering)
    for i in ordering
        name = bn.vars[i].name
        if !haskey(evidence, name)
            b = blanket(bn, a, i)
            a[name] = rand(b)[name]
        end
    end 
end

function gibbs_sample!(a, bn, evidence, ordering, m)
    for j in 1:m
        update_gibbs_sample!(a, bn, evidence, ordering)
    end
end


# Pg 62
struct GibbsSampling
    m_samples # number of samples to use
    m_burnin  # number of samples to discard during burn-in
    m_skip    # number of samples to skip for thinning
    ordering  # array of variable indices
end

function infer(M::GibbsSampling, bn, query, evidence)
    table = FactorTable()
    a = merge(rand(bn), pairs(evidence)) # this had a bug: missing pairs()
    println(a) 
    gibbs_sample!(a, bn, evidence, M.ordering, M.m_burnin)
    for i in 1:(M.m_samples)
        gibbs_sample!(a, bn, evidence, M.ordering, M.m_skip)
        b = select(a, query)
        table[b] = get(table, b, 0) + 1
    end
    vars = filter(v->v.name ∈ query, bn.vars)
    return normalize!(Factor(vars, table))
end


#Pg 64

function infer(D::MvNormal, query, evidencevars, evidence)
    μ, Σ = D.μ, D.Σ.mat
    b, μa, μb = evidence, μ[query], μ[evidencevars]
    A = Σ[query,query]
    B = Σ[evidencevars,evidencevars]
    C = Σ[query,evidencevars]
    μ = μ[query] + C * (B\(b - μb))
    Σ = A - C * (B \ C')
    return MvNormal(μ, Σ)
end
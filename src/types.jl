#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
abstract type AbstractCliffordAlgebra end 

"""
GA{D, metric, sym}

Geometric algebra in `D` dimensions with metric `metric` and basis vectors 
denoted `sym`. 
"""
struct GA{D, metric, sym} <: AbstractCliffordAlgebra
    function GA(D::Int, sym::Symbol=:σ)
        metric = SMatrix{D,D}(diagm(0 => ones(D)))
        new{D, metric, sym}()
    end
    function GA(m::Int, n::Int, sym::Symbol=:γ)
        D = m + n
        d = vcat(ones(m), -ones(n))
        metric = SMatrix{D,D}(diagm(0 => d))
        if m == 1
            metric = OffsetArray(metric, 0:n, 0:n)
        end
        new{D, metric, sym}()
    end
end

# dim(::SMatrix{D,D}) where {D} = D
# dim(::GA{metric,sym}) where {metric} = dim(metric)
# dim(g::OffsetArray) where {T,U,metric} = dim(g.parent)
# metric(::GA{metric}) where {metric} = metric

function Base.show(io::IO, ::MIME"text/plain", ::GA{D, metric, sym}) where {D, metric, sym}
    ltx = REPL.REPLCompletions.latex_symbols;
    println(io, "Geometric Algebra in $(D) dimensions")
    println(io, " • Basis vectors")
    println(io, "    ", [string(sym)*ltx["\\_$i"]*" " for i in axes(metric)[1]]...)
    println(" • Metric")
    show(io, MIME"text/plain"(), metric)
end


#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
abstract type CliffordElement{𝒢} <: Number end

struct BasisVector{𝒢} <: CliffordElement{𝒢} # 𝒢 is typed as either \scrG<TAB> or \mscrG<TAB>
    index::Int
    BasisVector(𝒢::GA, i::Int) = new{𝒢}(i)
end

sym(::GA{D, metric, T}) where {D, metric, T} = T

function Base.show(io::IO, σ::BasisVector{𝒢}) where {𝒢}
    ltx = REPL.REPLCompletions.latex_symbols
    sub = ltx["\\_$(σ.index)"]
    print(io, string(sym(𝒢))*sub)
end

function Base.show(io::IO, ::MIME"text/plain", σ::BasisVector{𝒢}) where {𝒢}
    show(io, σ)
    println(io, "\nBasis vector $(σ.index) in $(𝒢)")
    # show(io, MIME"text/plain"(), 𝒢)
end

function vectors(𝒢::GA{D, metric, sym}) where {D,metric,sym}  # 𝒢 is typed as either \scrG<TAB> or \mscrG<TAB>
    broadcast(i -> BasisVector(𝒢, i), axes(metric)[1])
end


#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# struct Blade{𝒢, indices} <: CliffordElement{𝒢}
#     coeff::Number
# end

# function Blade(coeff, elems::Vector{CliffordElement})
#     for i in 1:(length(elems)-1)
#         if 
# end


#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
struct MultiVector{𝒢, Tuple} <: CliffordElement{𝒢}
    terms::SortedDict
end

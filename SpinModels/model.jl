using Random
using LinearAlgebra

abstract type Model end

mutable struct Ising <: Model
    L::Integer
    lattice :: AbstractArray
    rng :: Random.MersenneTwister
    J :: Real
    UpdateMethod::String
    initialM::Real
    initialE::Real
    samples::Integer
    name::String

    function Ising(dims::Integer,L::Integer, rng::Random.AbstractRNG,J::Real,UpdateMethod::String,samples::Integer)     #inner construction
        model = new()
        model.name="Ising"
        model.L=L
        if dims==2
            model.lattice = ones(L,L)
        elseif dims==3
            model.lattice=ones(L,L,L)
        end
        model.rng = rng
        model.J = J
        model.UpdateMethod=UpdateMethod
        model.initialE=-J * 3 * (L^3)
        model.initialM=L^3
        model.samples=samples
        return model
    end

end


Ising(L::Integer,UpdateMethod::String)=Ising(3,L,Random.seed!(Random.MersenneTwister(41)),1.0,UpdateMethod,10^5)     #outer construction, default


mutable struct Heisenberg <: Model
    L::Integer
    lattice :: AbstractArray
    rng :: Random.MersenneTwister
    J :: Real
    UpdateMethod::String
    samples::Integer
    name::String

    function Heisenberg(L::Integer,rng::Random.AbstractRNG,J::Real,UpdateMethod::String,samples::Integer)
        model=new()
        model.name="Heisenberg"
        model.L=L
        model.lattice = Array{Vector}(undef,(L,L,L))
        for  i in eachindex(model.lattice)
            r=randn(rng,3)
            r /= norm(r)
            model.lattice[i]=r
        end
        model.rng=rng
        model.J=J
        model.UpdateMethod=UpdateMethod
        model.samples=samples
        return model
    end

end

Heisenberg(L::Integer,UpdateMethod::String)=Heisenberg(L,Random.seed!(Random.MersenneTwister(41)),1.0,UpdateMethod,10^5) #outer construction, default

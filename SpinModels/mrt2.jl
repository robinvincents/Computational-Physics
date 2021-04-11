#TODO: declare function such that I do not have to use argument beta for ising type and argument expf for heisenberg type, while using same function name?

#Functions for Metropolis Update for Ising and Heisenberg model


using Statistics: mean, std, var
using Random
using LinearAlgebra
include("model.jl")


#FUNCTIONS USED BOTH FOR ISING MODEL & HEISENBERG MODEL

function pbcminus(L::Integer,i::Integer)     #periodic boundary conditions for subtracting 1

    if i>1
        return -1

    else
        return L-1
    end

end

function pbcplus(L::Integer,i::Integer)   #periodic boundary conditions for adding 1

    if i<L
        return 1
    else
        return -L+1
    end
end

function mrt2_nb(L::Integer,lattice::AbstractArray, site::CartesianIndex)   #compute sum of neighbour spins with periodic boundary conditions

    return lattice[site+CartesianIndex(0,0,pbcplus(L,site[3]))].+lattice[site+CartesianIndex(0,0,pbcminus(L,site[3]))].+lattice[site+CartesianIndex(0,pbcplus(L,site[2]),0)].+lattice[site+CartesianIndex(0,pbcminus(L,site[2]),0)].+
    lattice[site+CartesianIndex(pbcplus(L,site[1]),0,0)].+lattice[site+CartesianIndex(pbcminus(L,site[1]),0,0)]

end


function right_nb(L::Integer,lattice::AbstractArray,site::CartesianIndex)  #sum of neighbours,  faster than mod1 with pbcplus/pbcminus, here only necessary to get 3 neighbours for energy computation as we loop over all sites

    return lattice[site+CartesianIndex(0,0,pbcplus(L,site[3]))].+lattice[site+CartesianIndex(0,pbcplus(L,site[2]),0)].+lattice[site+CartesianIndex(pbcplus(L,site[1]),0,0)]

end



#########################################################
#FUNCTIONS USED FOR ISING MODEL


function mrt2_step!(model::Ising,L::Int64,lattice::AbstractArray, expf::AbstractArray{Float64,1},beta::Real,J::Real)           #model has to be of type Ising

    #choose random lattice site
    randomsite=rand(CartesianIndices(lattice))      #choose random lattice site
    nspins = mrt2_nb(L,lattice, randomsite)
    index = trunc(Int64, (nspins + 8) ÷ 2)                 #computes index depending on number neighbours for lookuptable of the exp. function
    R = (lattice[randomsite] == 1) ? expf[8-index] : expf[index]     #depending on sign of lattice site spin chooses according probability (index of lookup table)

    if R>rand()                                            #metropolis spin flip step

        lattice[randomsite] *= -1

    end

    return nothing

end


function computephysics(model::Ising,L::Integer,lattice::AbstractArray,J::Real) #computes energy and magnetization for given lattice

    E::Real=0.0
    M::Real=0.0

    for index in CartesianIndices(lattice)
        E += -lattice[index]*right_nb(L,lattice,index)
        M += lattice[index]
    end

    return E*J,abs(M)

end



##############################################################
#FUNCTIONS USED FOR HEISENBERG MODEL

function mrt2_step!(model::Heisenberg,L::Int64,lattice::AbstractArray,expf::AbstractArray{Float64,1},beta::Real,J::Real) #model has to be of type Heisenberg


    randomsite=rand(CartesianIndices(lattice))
    newspin = lattice[randomsite] .+ 0.5 .* (randn(model.rng,3))
    newspin ./= norm(newspin)

    ΔE=dot(lattice[randomsite] .- newspin,mrt2_nb(L,lattice,randomsite))*J

    p=exp(-ΔE*beta)

    if p>rand(model.rng)
        lattice[randomsite]=newspin
    end

    return nothing

end


function computephysics(model::Heisenberg,L::Integer,lattice::AbstractArray,J::Real)

    E::Real=0.0
    M = zeros(3)

    for index in CartesianIndices(lattice)
        E += -dot(lattice[index],right_nb(L,lattice,index))
        M[1] += lattice[index][1]
        M[2] += lattice[index][2]
        M[3] += lattice[index][3]
    end

    M_abs=sqrt(sum(abs2,M))

    return E*J,M_abs
end

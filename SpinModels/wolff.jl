#TODO: declare functions such that I do not have to use argument beta for ising type and argument p for heisenberg type, while using same function name?

#Functions for Wolff Update of Ising and Heisenberg model


using Random
using LinearAlgebra
include("model.jl")

#FUNCTIONS USED BOTH FOR ISING MODEL & HEISENBERG MODEL


function right_nb(L::Integer,lattice::AbstractArray,site::CartesianIndex)  #sum of neighbours,  faster than mod1 with pbcplus/pbcminus, here only necessary to get 3 neighbours for energy computation as we loop over all sites

    return lattice[site+CartesianIndex(0,0,pbcplus(L,site[3]))].+lattice[site+CartesianIndex(0,pbcplus(L,site[2]),0)].+lattice[site+CartesianIndex(pbcplus(L,site[1]),0,0)]

end

# function right_nb(L::Integer,lattice::AbstractArray,site::CartesianIndex)  #sum of neighbours, alternative way-->same runtime
#
#     return lattice[site[1],site[2],(site[3])&(L-1)+1]+lattice[site[1],(site[2])&(L-1)+1,site[3]]+lattice[(site[1])&(L-1)+1,site[2],site[3]]
# end

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

#########################################################
#FUNCTIONS USED FOR ISING MODEL


function computephysics(model::Ising,L::Integer,lattice::AbstractArray,J::Real) #computes energy and magnetization for given lattice

    E::Real=0.0
    M::Real=0.0

    for index in CartesianIndices(lattice)
        E += -lattice[index]*right_nb(L,lattice,index)
        M += lattice[index]
    end

    return E*J,abs(M)

end

function bonds_Ising!(L::Int,lattice::AbstractArray,site::CartesianIndex,state::Float64,p::Float64,J::Real) #grow cluster for ising model


        for index in [site+CartesianIndex(0,0,pbcminus(L,site[3])),site+CartesianIndex(0,0,pbcplus(L,site[3])),site+CartesianIndex(0,pbcminus(L,site[2]),0),site+CartesianIndex(0,pbcplus(L,site[2]),0),site+CartesianIndex(pbcplus(L,site[1]),0,0),site+CartesianIndex(pbcminus(L,site[1]),0,0)]

            if lattice[index]==state && rand()<=p
                lattice[index]*=-1
                bonds_Ising!(L,lattice,index,state,p,J)
            end

        end

        return nothing
end

function wolff_clusterupdate!(model::Ising,L::Int,lattice::AbstractArray,p::Float64,beta::Float64,J::Float64)    #model has to be of type Ising

    randomsite=rand(CartesianIndices(lattice))      #choose random lattice site
    state=lattice[randomsite]
    lattice[randomsite]*=-1
    bonds_Ising!(L,lattice,randomsite,state,p,J)

    return nothing
end

##########################################################
#Metropolis update functions including tracking of the number of flipped spins


function bonds_Ising_trackspin!(L::Int,lattice::AbstractArray,site::CartesianIndex,state::Float64,p::Float64,J::Float64,flippedspins::Int) #grow cluster for ising model


        for index in [site+CartesianIndex(0,0,pbcminus(L,site[3])),site+CartesianIndex(0,0,pbcplus(L,site[3])),site+CartesianIndex(0,pbcminus(L,site[2]),0),site+CartesianIndex(0,pbcplus(L,site[2]),0),site+CartesianIndex(pbcplus(L,site[1]),0,0),site+CartesianIndex(pbcminus(L,site[1]),0,0)]

            if lattice[index]==state && rand()<=p
                lattice[index]*=-1
                flippedspins+=1
                flippedspins=bonds_Ising_trackspin!(L,lattice,index,state,p,J,flippedspins)
            end

        end

        return flippedspins
end




function wolff_clusterupdate_trackspin!(model::Ising,L::Int,lattice::AbstractArray,p::Float64,beta::Float64,J::Float64)    #model has to be of type Ising

    randomsite=rand(CartesianIndices(lattice))      #choose random lattice site
    state=lattice[randomsite]
    lattice[randomsite]*=-1
    flippedspins=1
    flippedspins=bonds_Ising_trackspin!(L,lattice,randomsite,state,p,J,flippedspins)

    return flippedspins        #return total number of flipped spins for Wolff step
end

##############################################################
#FUNCTIONS USED FOR HEISENBERG MODEL

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

function bonds_Heisenberg!(L::Int,lattice::AbstractArray,beta::Real,site::CartesianIndex,randomplane::Vector,J::Float64) #grow cluster for heisenberg model

    dot_center=dot(lattice[site],randomplane) #precompute for loop

    for index in [site+CartesianIndex(0,0,pbcminus(L,site[3])),site+CartesianIndex(0,0,pbcplus(L,site[3])),site+CartesianIndex(0,pbcminus(L,site[2]),0),site+CartesianIndex(0,pbcplus(L,site[2]),0),site+CartesianIndex(pbcplus(L,site[1]),0,0),site+CartesianIndex(pbcminus(L,site[1]),0,0)]

        dot_nb=dot(lattice[index],randomplane) #to avoid double computation
        if rand()<(1-exp(2*beta*J*dot_center*dot_nb))

                lattice[index]=lattice[index]-2 .* randomplane .* dot_nb
                bonds_Heisenberg!(L,lattice,beta,index,randomplane,J)

        end

    end

    return nothing

end

function wolff_clusterupdate!(model::Heisenberg,L::Int,lattice::AbstractArray,p::Real,beta::Real,J::Float64)  #model has to be of type Heisenberg

    randomplane=randn(model.rng,3)
    randomplane ./= norm(randomplane)
    randomsite=rand(CartesianIndices(lattice))


    lattice[randomsite]=lattice[randomsite]-2 .* randomplane .* dot(lattice[randomsite],randomplane)

    bonds_Heisenberg!(L,lattice,beta,randomsite,randomplane,J)

    return nothing
end


##############################################
#Wolff update functions including tracking of the number of flipped spins

function bonds_Heisenberg_trackspin!(L::Int,lattice::AbstractArray,beta::Real,site::CartesianIndex,randomplane::Vector,J::Float64,flippedspins::Int) #grow cluster for heisenberg model

    dot_center=dot(lattice[site],randomplane) #precompute for loop

    for index in [site+CartesianIndex(0,0,pbcminus(L,site[3])),site+CartesianIndex(0,0,pbcplus(L,site[3])),site+CartesianIndex(0,pbcminus(L,site[2]),0),site+CartesianIndex(0,pbcplus(L,site[2]),0),site+CartesianIndex(pbcplus(L,site[1]),0,0),site+CartesianIndex(pbcminus(L,site[1]),0,0)]

        dot_nb=dot(lattice[index],randomplane) #to avoid double computation
        if rand()<(1-exp(2*beta*J*dot_center*dot_nb))

                lattice[index]=lattice[index]-2 .* randomplane .* dot_nb
                flippedspins+=1
                flippedspins=bonds_Heisenberg_trackspin!(L,lattice,beta,index,randomplane,J,flippedspins)

        end

    end

    return flippedspins

end


function wolff_clusterupdate_trackspin!(model::Heisenberg,L::Int,lattice::AbstractArray,p::Real,beta::Real,J::Float64)  #model has to be of type Heisenberg

    randomplane=randn(model.rng,3)
    randomplane ./= norm(randomplane)
    randomsite=rand(CartesianIndices(lattice))


    lattice[randomsite]=lattice[randomsite]-2 .* randomplane .* dot(lattice[randomsite],randomplane)
    flippedspins=1
    flippedspins=bonds_Heisenberg_trackspin!(L,lattice,beta,randomsite,randomplane,J,flippedspins)

    return flippedspins
end

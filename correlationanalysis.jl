#CODE FOR AUTOCORRELATION ANALYSIS and RUNTIME ANALYSIS


using Plots; pyplot()
using StatsBase
using LsqFit
using DataFrames
using BenchmarkTools
using JLD
using GLM
include("wolff.jl")
include("mrt2.jl")
include("MCS.jl")
include("model.jl")



function Metropolis_autocorrelation(model,T::Real,n_e::Integer,n_0::Integer,stepspersweep::Integer)

    beta=1/T
    J = model.J
    L = model.L
    lattice = model.lattice
    modeltype=model.name


    expf = exp.(2*beta*J*[-6:2:6...])            #get lookup table for given temperature


    for i in 1:n_e*stepspersweep
        mrt2_step!(model,L,lattice,expf,beta,J)
    end

    energy=zeros(n_0)

    energy[1],_=computephysics(model,L,lattice,J)

    for i in 2:n_0
        for j in 1:stepspersweep
            mrt2_step!(model,L,lattice,expf,beta,J)
        end
        energy[i],_=computephysics(model,L,lattice,J)
    end

    correlation=autocor(energy,0:n_0-2)
    mu = mean(correlation)
    finalcorr=zeros(200)

    #average correlation function up to 100 sweeps
    for i in 1:(size(correlation,1)-200)
        diff = correlation[i]-mu
        for j in 1:200
            finalcorr[j] += diff * (correlation[i+j]-mu)
        end
    end
    #normalize
    norm = 1 ./ finalcorr[1]
    for i in 1:size(finalcorr,1)
        finalcorr[i] *= norm
    end

    tau=0.5                                     #compute tau of the correlation series, with cutoff at 0
    for i in 2:size(finalcorr,1)
        if finalcorr[i]<0
            break
        end
        tau+=finalcorr[i]
    end

    @show tau


    save("Ex5/$modeltype-correlation-Metropolis-$L-$T.jld","linearcorrelation",finalcorr,"tau",tau)

end

#test1=Heisenberg(12,"Metropolis")
#Metropolis_autocorrelation(test1,1.44,500,50000,12^3)

function Wolff_autocorrelation(model,T,n_e,n_0,stepspersweep)

    beta=1/T
    J = model.J
    L = model.L
    lattice = model.lattice
    modeltype=model.name

    p=1-exp(-2*beta*J)

    for i in 1:n_e*stepspersweep           #n_e*stepspersweep=number of Wolff steps to equilibrate
        wolff_clusterupdate!(model,L,lattice,p,beta,J)
    end

    E,M=computephysics(model,L,lattice,J)

    energy=zeros(n_0)
    energy[1]=E
    for j in 2:n_0                             #n_0 number of system sweeps to be performed
        for k in 1:stepspersweep              #number of Wolff steps for one sweep
            wolff_clusterupdate!(model,L,lattice,p,beta,J)
        end
        E,M=computephysics(model,L,lattice,J)
        energy[j]=E
    end

    correlation=autocor(energy,0:n_0-2)
    mu = mean(correlation)
    finalcorr=zeros(100)
    #average correlation function up to 100 sweeps
    for i in 1:(size(correlation,1)-100)
        diff = correlation[i]-mu
        for j in 1:100
            finalcorr[j] += diff * (correlation[i+j]-mu)
        end
    end
    #normalize
    norm = 1 ./ finalcorr[1]
    for i in 1:size(finalcorr,1)
        finalcorr[i] *= norm
    end

    tau=0.5                                     #compute tau of the correlation series, with cutoff at 0
    for i in 2:size(finalcorr,1)
        if finalcorr[i]<0
            break
        end
        tau+=finalcorr[i]
    end

    @show tau

    save("Ex5/$modeltype-correlation-Wolff-$L-$T.jld","linearcorrelation",finalcorr)

end

#test2=Heisenberg(12,"Wolff")
#Wolff_autocorrelation(test2,4.51,500,50000,18)


function dynamicalexponent(L_list,tau_list)

    X=collect(L_list)
    data=DataFrame(X=log.(L_list),Y=log.(tau_list))
    ols = lm(@formula(Y ~ X), data)
    print(ols)

    scatter((L_list),(tau_list),xlabel="System size",ylabel="tau",size=(400,400),legend=false,smooth=true)
    savefig("Ex5/dynamicalexponent_Ising.png")

end

L_list_Heisenberg=[4,6,8,10]
tau_list_Heisenberg=[11.1,29.3,37.8,48.0]
L_list_Ising=[6,8,10,12]
tau_list_Ising=[14.4,20.0,30.3,47.2]

dynamicalexponent(L_list_Heisenberg,tau_list_Heisenberg)


function runtime(model,T,n_e,n_0,stepspersweep)

    beta=1/T
    J = model.J
    L = model.L
    lattice = model.lattice
    updatemethod=model.UpdateMethod

    if updatemethod=="Metropolis"

        expf = exp.(2*beta*J*[-6:2:6...])            #lookup table for given temperature

        for i in 1:n_e*stepspersweep           #n_e*stepspersweep=number of Wolff steps to equilibrate
            mrt2_step!(model,L,lattice,expf,beta,J)
        end

        @time begin

            for j in 1:n_0*stepspersweep
                mrt2_step!(model,L,lattice,expf,beta,J)
            end

        end

    elseif updatemethod=="Wolff"

        p=1-exp(-2*beta*J)

        for i in 1:n_e*stepspersweep           #n_e*stepspersweep=number of Wolff steps to equilibrate
            wolff_clusterupdate!(model,L,lattice,p,beta,J)
        end

        @time begin

            for j in 1:n_0*stepspersweep
                wolff_clusterupdate!(model,L,lattice,p,beta,J)
            end

        end

    end

end

#runtime1=Ising(8,"Wolff")
# runtime(runtime1,4.51,500,500,5)

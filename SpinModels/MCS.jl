using JLD
include("model.jl")
include("mrt2.jl")
include("wolff.jl")


function setdecorr(model,beta::Real,n_0_list)       #set equilibration,decorrelation parameters

    updatemethod=model.UpdateMethod
    L=model.L

    if updatemethod=="Metropolis"

        n_0 = (2*L^(3))
        n_e = (500*L^(3))

    elseif updatemethod=="Wolff"

        n_0=n_0_list[L]
        n_e=500

    end

    return n_0, n_e

end



function ensemble(model, T::Real, n_0::Integer, n_e::Integer)       #computes magnetization,energy and susceptibility for given temperature

    beta=1/T
    J = model.J
    L = model.L
    lattice = model.lattice
    updatemethod=model.UpdateMethod


    samples = model.samples

    expf = exp.(2*beta*J*[-6:2:6...])            #lookup table for given temperature
    p=1-exp(-2*beta*J)                           #probability for Wolff update

    M_array = zeros(samples)                                      #array to store magnetization,energy values
    E_array = zeros(samples)

    if updatemethod=="Metropolis"

        for i in 1:n_e         #n_e steps to equilibrate
            mrt2_step!(model,L,lattice,expf,beta,J)
        end


        for i in 1:samples
            for j in 1:n_0         #perform for each samples decorrelation of n_0 steps
                mrt2_step!(model,L,lattice,expf,beta,J)
            end
            E_array[i],M_array[i]=computephysics(model,L,lattice,J)
        end

    elseif updatemethod=="Wolff"

        for i in 1:n_e                          #n_e steps to equilibrate
            wolff_clusterupdate!(model,L,lattice,p,beta,J)
        end

        totalflipped=0

        for _ in 1:1000                      #determine average number of flipped spin per update step
            flipped=wolff_clusterupdate_trackspin!(model,L,lattice,p,beta,J)
            totalflipped+=flipped
        end

        stepspersweep=trunc(Int,(L^3)*1000/(totalflipped))
        @show T,stepspersweep


        for i in 1:samples
            for _ in 1:n_0*stepspersweep                        #n_0*stepspersweep=number of system sweeps to decorrelate between samples
                wolff_clusterupdate!(model,L,lattice,p,beta,J)
            end
            E_array[i],M_array[i]=computephysics(model,L,lattice,J)
        end

    end

    M_array ./= L^3
    E_array ./= L^3

    M_mean = mean(M_array)
    E_mean = mean(E_array)
    M_var = var(M_array)
    E_var = var(E_array)
    C = (L^3 * E_var) * (beta^2)   #compute heat capacity
    S = (L^3 * M_var) * beta      #compute susceptibility
    M_array_4 = (M_array).^4
    M_array_2= (M_array).^2
    bindercumulant=1.0-(mean(M_array_4)/(3*mean(M_array_2)^2)) #compute binder cumulant

    return M_mean,E_mean,C,S,bindercumulant
end



function varytemperature(model,T_array,n_0_list)     #perform the ensemble computation for temperatures given in T_array,n_0_list is list of decorrelation sweeps needed for Wolff algorithm

    updatemethod=model.UpdateMethod
    L=model.L
    modeltype=model.name

    temperature_range = size(T_array,1)

    magnetization = zeros(temperature_range)
    energy = zeros(temperature_range)
    heatcapacity = zeros(temperature_range)
    susceptibility = zeros(temperature_range)
    bindercumulant=zeros(temperature_range)

    for i in 1:temperature_range
        n_0, n_e = setdecorr(model,T_array[i],n_0_list)
        magnetization[i], energy[i], heatcapacity[i], susceptibility[i], bindercumulant[i]=ensemble(model, T_array[i], n_0, n_e)
    end

    save("Ex5/$modeltype-$updatemethod-$L.jld", "T", T_array,"magnetization",magnetization,"energy",energy,"heatcapacity",heatcapacity,"susceptibility",susceptibility,"bindercumulant",bindercumulant)

    return nothing
end



# for L in [8]
#     test3=Heisenberg(L,"Wolff")
#     varytemperature(test3,[1.43,1.44,1.45],[0,0,0,15,0,30,0,40])
# end

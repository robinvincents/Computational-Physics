using Plots; pyplot()
using JLD
using LaTeXStrings
using GLM
using DataFrames
using LsqFit
using HDF5




function readdata(filename)                      #load data from jld files

    T=load(filename)["T"]
    magnetization=load(filename)["magnetization"]
    energy=load(filename)["energy"]
    heatcapacity=load(filename)["heatcapacity"]
    susceptibility=load(filename)["susceptibility"]
    bindercumulant=load(filename)["bindercumulant"]

    return T,magnetization,energy,heatcapacity,susceptibility,bindercumulant
end


function plotsystemsizes(L_list)    #plot values of magnetization,energy,susceptibility,heat capacity,bindercumulant as function of temperature for each system and together for varios system sizes
                                    #NOTE: specify datafilename in the function!


    p_energy=plot(size=(400,400))
    p_magnetization=plot(size=(400,400))
    p_susceptibility=plot(size=(400,400))
    p_heatcapacity=plot(size=(400,400))
    p_bindercumulant=plot(size=(400,400))


    for L in L_list

        file=load("Ex5/Heisenberg-Wolff-$L.jld")
        T_L,e_L,m_L,s_L,c_L,bindercumulant_L=file["T"],file["energy"],file["magnetization"],file["susceptibility"],file["heatcapacity"],file["bindercumulant"]


        plot(T_L,e_L)
        scatter!(T_L,e_L,label="E",xlabel="Temperature T",ylabel="Energy E per Spin",size=(400,400),legend=false)
        savefig("energy_$L.png")
        plot(T_L,m_L)
        scatter!(T_L,m_L,label="M",xlabel="Temperature T",ylabel="Magnetization M per Spin",size=(400,400),legend=false)
        savefig("magnetization_$L.png")
        plot(T_L,s_L)
        scatter!(T_L,s_L,label="M",xlabel="Temperature T",ylabel="Susceptibility",size=(400,400),legend=false)
        savefig("susceptibility_$L.png")
        plot(T_L,c_L)
        scatter!(T_L,c_L,label="M",xlabel="Temperature T",ylabel="Heat Capacity",size=(400,400),legend=false)
        savefig("heatcapacity_$L.png")
        plot(T_L,bindercumulant_L)
        scatter!(T_L,bindercumulant_L,label="M",xlabel="Temperature T",ylabel="Binder Cumulant",size=(400,400),legend=false)
        savefig("bindercumulant_$L.png")

        plot!(p_energy,T_L, e_L,
        xlabel="Temperature T",
        ylabel="Energy E per Spin",
        size=(400, 400),label="L=$L",marker=2)

        plot!(p_magnetization,T_L, m_L,
         xlabel="Temperature T",
         ylabel="Magnetization M per Spin",
         size=(400, 400),label="L=$L",marker=2)

        plot!(p_susceptibility,T_L, s_L,
        xlabel="Temperature T",
        ylabel="Susceptibility",
        size=(400, 400),label="L=$L",marker=2)

         plot!(p_heatcapacity,T_L, c_L,
         xlabel="Temperature T",
         ylabel="Heat Capacity",
         size=(400, 400),label="L=$L",marker=2)

         plot!(p_bindercumulant,T_L,bindercumulant_L,xlabel="Temperature T",
                  ylabel="Binder Cumulant",
                  size=(400, 400),label="L=$L",marker=2,xticks=T_L[1]:0.01:T_L[end])


    end

    savefig(p_energy,"energy_systemsizes.png")
    savefig(p_magnetization,"magnetization_systemsizes.png")
    savefig(p_susceptibility,"susceptibility_systemsizes.png")
    savefig(p_heatcapacity,"heatcapacity_systemsizes.png")
    savefig(p_bindercumulant,"bindercumulant_systemsizes.png")
end

#plotsystemsizes([4,6])

function getexponents(L_list,betalist)           #define range of inverse temperatures simulated in betalist and list of system sizes in L_list
                                                #plots data collapse for the system sizes and maximal susceptibility as a function system size

    beta_c=1/4.51
    maxsusc=[]

    p_datacollapse = plot(xlabel=L"Inverse temperature $\beta$", ylabel=L"Susceptibility $\chi$", size=(400, 400))

    for L in L_list
        s_L=load("Ising_critical_$L.jld")["susceptibility"]
        scatter!(p_datacollapse,L^(1/0.65)*(1 .- betalist ./ beta_c), s_L*L^(-1.27/0.65),
         xlabel=L"L$^{1/\nu}(1-\beta/\beta_c)$",
         ylabel=L"L$^{\gamma/\nu}\chi$",
         size=(400, 400),label="L=$L")
         append!(maxsusc,maximum(s_L))
    end

    X=collect(L_list)                                               #estimation of critical exponents ratio gamma/nu (linear fit-->slope of log plot)
    data=DataFrame(X=log.(L_list),Y=log.(maxsusc))
    ols = lm(@formula(Y ~ X), data)
    print(ols)

    savefig("Ex5/susceptibility_datacollapse.png")

    scatter((L_list),(maxsusc),xlabel="System size",ylabel=L"$\chi_{max}$",size=(400,400),xscale=:log10,yscale=:log10,legend=false,smooth=true)
    savefig("Ex5/maxsusp_systemsize.png")

end


function plotcorrelation(model::String,algorithm::String,L_array,T_array,n_0)  #NOTE: specify name of algorithm, L is system size, T_array list of temperatures to plot, n_e=equilibration sweeps, n_0=number of samples, stepspwersweep=number of steps per sweep

    plot()

    for L in L_array

        for T in T_array

            linearcorr_T=load("$model-correlation-$algorithm-$L-$T.jld")["linearcorrelation"]
            plot!(0:n_0-1,linearcorr_T,label="L=$L,T=$T",xlabel="System Sweeps",ylabel="Correlation")
            savefig("Ex5/$model-correlation-$algorithm-$L.png")

        end

    end

end

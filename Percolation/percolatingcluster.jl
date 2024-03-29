using Random
using Plots

function percolating(L,p)

    

    lattice = zeros(Int64, L, L)
    fire=true
    start=true
    
    spanningcluster=false
    k=1

    while spanningcluster==false && k<=L                  #repeat for different starting fire starting points [1,k] while no spanning 								  cluster formed and [1,k] still in lattice
	fire=true					  #reset fire indicator to enable entering while loop
	lattice[1,k]=0					  #set unsuccessful fire starting point to not belonging to the cluster
    	k+=1	
	start=true					  #reset start indicator						

	    while fire==true				  #looping over the grid while fire is still on
	    fire=false					  #resetting fire indicator

		for i in 1:L 				  #looping over lattice 
		    for j in 1:L
		        random = rand()

		        if start==true && random<=p	  #at first iteration fill lattice according to probability p
		            lattice[i,j] = 1
		        end
		        
	     		if j>=k && lattice[1,j]==1 && fire==false && start==true    #set occupied site in first row to fire, only if 											    no fire in current try to obtain spanning cluster 											    and in first iteration
			    lattice[1,j]=2
			    fire=true 
			end        	
		        
		        if start==false && lattice[i,j]==2                          #after first iteration update fire and burnt trees
		            
		                if i<L && lattice[i+1,j]==1
		                    lattice[i+1,j]+=1 
		                end

		                if j<L && lattice[i,j+1]==1
		                    lattice[i,j+1]+=1
		                end

		                if i>1 && lattice[i-1,j]==1
		                    lattice[i-1,j]+=1
		                end

		                if j>1 && lattice[i,j-1]==1
		                    lattice[i,j-1]+=1
		                end

		                lattice[i,j]+=1
		                
		                fire=true
		        end
		        
		        if lattice[L,j]==2			#if lattice on last row on fire,obtained spanningcluster
		            spanningcluster=true
		        end
		    end
		end
		
		start=false                                     #after first iteration set start to false
	    end
	end
return spanningcluster, lattice
end  



function obtaincluster(L,p)
	
	spanningcluster=false
	lattice=Int64[]
	while spanningcluster == false                         #obtain new samples until spanning cluster obtained
		spanningcluster, lattice = percolating(L,p)
	end
	
	for i in 1:L					       #set burnt trees to 1, all others to 0
		for j in 1:L
			if lattice[i,j]==3
				lattice[i,j]=1
			else lattice[i,j]=0
			end
			

		end
	end

	return lattice
end

function sandbox(L,p,R)                             #computes number of occupied sites for given box size R
	
	lattice=obtaincluster(L,p)
	counter=0
	
	while lattice[L÷2,L÷2]!=1                       #makes sure that center of box R is occupied
		lattice=obtaincluster(L,p)
	end

	for i in 0:R                                    #counts number of occupied sites in R
		for j in 0:R
			counter+=lattice[L÷2-R÷2+i,L÷2-R÷2+j]
		end
	end
	
	return counter
end

function sandboxstat(L,p,R,stepsize)               #for a given stepsize of increasing box size R, makes array of boxsizes of varying R and of number of occupied sites, return in log scaling
	
	boxsize=Int64[]
	occupiedsites=Int64[]
	
	while R<L
		push!(boxsize,R)
		push!(occupiedsites,sandbox(L,p,R))
		R+=stepsize
	end

	return log.(boxsize), log.(occupiedsites)
end



plot(sandboxstat(1000,0.6,10,2)[1],sandboxstat(1000,0.6,10,2)[2])
savefig("sandbox.png")

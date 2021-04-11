using Random






function HoshenKopelman(L,p)


    k = 2
    
    lattice=zeros(Int64,2,L)
    M=Int64[]
    push!(M,0)
    
    for i in 1:L
        
        lattice[1,:]=lattice[2,:]

        
        for j in 1:L
            

            
            random=rand()
            
            if random<=p
                lattice[2,j]=1
            end
            

            if i>1 && j>1 && lattice[2,j]==1
                
                if lattice[2,j-1]==0 && lattice[1,j]==0
                    lattice[2,j]=k
                    push!(M,1)
                    k+=1
                end
                
                
                if lattice[2,j-1]!=0 && lattice[1,j]==0
                    

                    lattice[2,j]=lattice[2,j-1]
                    
                                        
                    if M[lattice[2,j-1]]>=0
                        M[lattice[2,j-1]]+=1
                    end
                    
                    if M[lattice[2,j-1]]<0
                        realMleft=M[lattice[2,j-1]]
                        index=0
                        while realMleft<0
                            index=-1*realMleft
                            realMleft=M[-1*realMleft]
                        end
                        M[index]+=1
                    end
                    
                end
                
                if lattice[2,j-1]==0 && lattice[1,j]!=0
                
            
                    lattice[2,j]=lattice[1,j]
                    
                    if M[lattice[1,j]]>=0
                        M[lattice[1,j]]+=1
                    end
                    
                    if M[lattice[1,j]]<0
                        realMabove=M[lattice[1,j]]
                        index=0
                        while realMabove<0
                            index=-1*realMabove
                            realMabove=M[-1*realMabove]
                        end
                        M[index]+=1
                    end
                    
                    
                end
                
                if lattice[2,j-1]!=0 && lattice[1,j]!=0
                    
                
                    
                    lattice[2,j]=lattice[1,j]

                    if M[lattice[2,j-1]]>=0 && M[lattice[1,j]]>=0
                        M[lattice[1,j]]+=1+M[lattice[2,j-1]]
                    end
                    
                    if M[lattice[2,j-1]]<0 && M[lattice[1,j]]<0
                        realMabove=M[lattice[1,j]]
                        realMleft=M[lattice[2,j-1]]
                        indexabove=0
                        indexleft=0
                        while realMabove<0
                            index=-1*realMabove
                            realMabove=M[-1*realMabove]
                        end
                        while realMleft<0
                            index=-1*realMleft
                            realMleft=M[-1*realMleft]
                        end
                        M[indexabove]+=1+M[indexleft]
                    end
                            
                    
                    if M[lattice[2,j-1]]<0
                        realMleft=M[lattice[2,j-1]]
                        index=0
                        while realMleft<0
                            index=-1*realMleft
                            realMleft=M[-1*realMleft]
                        end 
                        M[lattice[1,j]]+=1+M[index]
                    end                        
                    
                    if M[lattice[1,j]]<0
                        realMabove=M[lattice[1,j]]
                        index=0
                        while realMabove<0
                            index=-1*realMabove
                            realMabove=M[-1*realMabove]
                        end
                        M[index]+=1+M[index]
                    end
                        
                        
                    M[lattice[2,j-1]]=-lattice[1,j]
                    
                end
                
            end
            
        end
    end
    
    return M
end

HoshenKopelman(10,0.8)
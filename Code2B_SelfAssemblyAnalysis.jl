###General Remarks###
#The Code first reads in a trajectory and creates structures for every Chain and Coild-Coil which contain the relevant positions as well as a System structure.
#The Analysis Calculations are ususally done for the whole trajectory at once.
#To understand the Code it is best to understand the initial strucutre types and see how the data is saved/strucutred or play around if the short example trajectories. 
#This code can run on multiple cores and if this is not nativels supported from a used IDE, you can start julia with additonal cores with julia -p x where x is the number of processors
#The results are saved as JLD2 format at each step and can be read out after the code finishes or in case of a crash. An example readout script is given on github
#The Code contains a lot of legacy parts of old calculations and no longer used features that might still be usefull if one wants to look at the data in detail
#The chains mentioned in the code are mostly the alpha helices of the coiled-coil
#This is not intended for general usage but can be used to recalculate the resulsts of my paper. There might be weaird hard coded parts that could be solved dynamically. Since I only wotk with this system for the foreseeble future I did not feel it necessary to change this. So be carfuel if you want to adept the code
#I generally tried to have discriptive names for my functions and variables. Since I worked on this code over a long time at different parts there is not really one uniform style. 
#Most of the comments where done much later than the coding. 

using Chemfiles
using LinearAlgebra
using StaticArrays
using Combinatorics
using DelimitedFiles # Do I Still nedd this one?
using StatsBase
using LsqFit
using Statistics
using JLD2

flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Start")

const Number_CA=8  # Number of C_alpha from Leucine 
const Number_NH3=27 #Number of H in NH3+ from N-Terminus and Lysine sidechain 
const Number_COO=10 #Number O in COO- from C-Terminus and Glutamatic acid sidechain

const Number_Hydrogen_Bonds=(floor(Int64,Number_NH3/3),floor(Int64,Number_COO/2))
const hydrogen_bond_cutoff=3.5 # Hydrogen Bond cutoff. 

#Stucture that contains the Atom id number for the different used atoms in a Chains
struct Chain
    ID::Int64
    CA::SVector{Number_CA,Int64}
    NH3::SVector{Number_NH3,Int64}
    COO::SVector{Number_COO,Int64}
end

#Structure that contains the positions for every frame for the used system
struct Chain_Traj
    ID::Int64
    CA::Array{Float64,3}
    NH3::Array{Float64,3}
    COO::Array{Float64,3}
end

#Structure that constains the two chains or helices that make up a coiled-coil. It contains two
struct Coil_Traj
    ID::Int64
    ChainA::Chain_Traj
    ChainB::Chain_Traj
    Beads::Array{Float64,3}             #Mean pbc Position of the 8 pairs of Leucine C_alpha atoms making up the bead chain.  
    Bead_distance::Array{Float64,2}     #Distance between the two C_alpha of each bead. Can be used to check if the coiled-coils untangle at the end. Currently unused
    Directions::Array{Float64,2}        #Direction vector between the first and last bead of a chain bead. Currently unused for this application
end

#Structure that contains variables of the System
struct System_parameters
    name::String 
    frames::Int64                       #Total number of frames
    Frames::Array{Int64,1}              #List of frames
    Number_Atoms::Int64                 #Number of Atoms in the System
    Number_Coils::Int64                 #Number of Coiled-Coils in the System
    Number_Chains::Int64                #Number of Chains/ alpha helices in the System
    Box::Array{Float64,2}               #Array with all Box size Parameters 
    RBox::Array{Float64,2}              #Array with the reciprocal Box size parameters. Needed for faster pbc calculation
    Atom_Names::Array{String,1}         
    Residue_Names::Array{String,1}
end

exp_fit(L,P)=exp.(-L./P)                                    #Exponential decay funtion to allow the fitting of the persistence length
#Krohn(L,P)=2 .*P.*L.*(1 .-P./L.*(1 .-exp.(-L./P)))
Krohn2(L,P)=2 .* P .* L .- 2 .*P.^2 .*(1 .- exp.(-L./P))    #Krohns-function to allow the fitting of the persistence length    

vector_length(Positions)=sqrt.(sum((Positions).^2,dims=1))[1,:,:]       #Function that takes in a array of position/distance vectors and gives out a vector with the distances

#Repositions a point into the box if it outside the periodic boundary condition
function pbc_point(A,box)
    A .-= floor.(Int64, A ./box) .* box
    return A
end

#Calculates the distance between two points under consideration of the periodic boundary conditions. Returns a vector with distances
function pbc_dist(A,B,box,r_box)
    dx=A.-B
    dx .-= round.(Int64, dx .* r_box) .* box
    dist=[sqrt.(mapslices(sum,(dx.^2),dims=1))...]
    return dist
end

#Calculates the distance between two points under consideration of the periodic boundary conditions. Keeps the original from
function pbc_dist2(A,B,box,r_box)
    dx=A.-B
    dx .-= round.(Int64, dx .* r_box) .* box
    dist=sqrt.(mapslices(sum,(dx.^2),dims=1))
    return dist
end

#Calculates a direction vector between two points under consderation of the periodic boundary conditions. It only uses the frames that are specified
function pbc_vector(A,B,System,frames)
    box=System.Box[frames]
    r_box=System.RBox[frames]
    dx=A.-B
    dx .-= floor.(Int64, dx .* r_box .+ 0.5) .* box
    return dx
end

#Calculates a direction vector between two points under consderation of the periodic boundary conditions.
function pbc_vector2(A,B,box)
    r_box=1 ./box
    dx=A.-B
    dx .-= floor.(Int64, dx .* r_box .+ 0.5) .* box
    return dx
end

#Calculates a direction vector between two points under consderation of the periodic boundary conditions. Calculation is possible for the whole trajectory at once
function pbc_traj_vec(AA,BB,box)
    N=size(AA)[2]
    box=reshape(repeat(box,inner=(1,N)),(3,N,size(AA)[3]))
    r_box=1 ./box
    dx=AA.-BB
    dx .-= floor.(Int64, dx .* r_box .+ 0.5) .* box
    return dx
end

#Calxulates all possible distance pairs between elements of to lists. Returns a Matrix with distances. Considers PBC
function all_distances_frame(PositionsA,PositionsB,box,r_box)
    NA=size(PositionsA)[2]
    NB=size(PositionsB)[2]
    Positions_A=repeat(PositionsA,inner=(1,NB))
    Positions_B=repeat(PositionsB,outer=(1,NA))
    distances=pbc_dist(Positions_A,Positions_B,box,r_box)
    return [minimum(reshape(distances,NB,NA),dims=1)...]
end

#Reads the first frame of a trajectory and creates a topology from it. Needs gro file for it. From this it creates a list of atomnames and creates a mask for each of the three important atom types
#It detects the different chains or helices and assignes the correct atom index to each. The main return is the index Chains structure
function create_index_Chains(traj)
    frame=read_step(traj,0)
    Number_Atoms=length(frame)
    atom_index=collect(1:Number_Atoms)
    top=Topology(frame)

    Atom_Names=Array{String}(undef,Number_Atoms)
    Residues_Names=Array{String}(undef,Number_Atoms)
    for i in 1:Number_Atoms
        Atom_Names[i]=name(Chemfiles.Atom(frame,i-1))
        Residues_Names[i]=name(residue_for_atom(top,i-1))
    end
    mask_CA=("CA" .== Atom_Names).&("LEU" .== Residues_Names)
    mask_NH3=("H1" .== Atom_Names).|("H2" .== Atom_Names).|("H3" .== Atom_Names).|("HZ1" .== Atom_Names).|("HZ2" .== Atom_Names).|("HZ3" .== Atom_Names)
    mask_COO=("OE1" .== Atom_Names).|("OE2" .== Atom_Names).|("OC1" .== Atom_Names).|("OC2" .== Atom_Names)

    CA_index=atom_index[mask_CA]
    NH3_index=atom_index[mask_NH3]
    COO_index=atom_index[mask_COO]

    Number_Chains=Int64(length(CA_index)/Number_CA)
    Number_Coils=floor(Int64,Number_Chains/2)
    Chains=Array{Chain}(undef,Number_Chains)

    for i in 1: Number_Chains
        CA=(CA_index[(i-1)*Number_CA+1:(i-1)*Number_CA+Number_CA])
        NH3=(NH3_index[(i-1)*Number_NH3+1:(i-1)*Number_NH3+Number_NH3])
        COO=(COO_index[(i-1)*Number_COO+1:(i-1)*Number_COO+Number_COO])
        Chains[i]=Chain(i,CA,NH3,COO)
    end
    return Chains,Number_Chains,Number_Coils,Number_Atoms,Atom_Names,Residues_Names
end

#Main read in functions. Reads a trajectory and gives back a Matrix with all positions per frame, index Chain structures for each Chain and creates the System structure with the box parameters for the trajectroy.
function read_System(filename,fname)
    traj=Trajectory(filename)
    frames=Int64(length(traj))
    Chains,Number_Chains,Number_Coils,Number_Atoms,Atom_Names,Residue_Names=create_index_Chains(traj)
    i=0
    Positions=zeros(3,Number_Atoms,frames)
    Box=zeros(3,frames)
    frame=Frame()
    while i < frames
        read_step!(traj,i,frame)
        box=diag(matrix(UnitCell(frame)))
        Positions[:,:,i+1]=positions(frame)
        Box[:,i+1]=box
        i+=1
    end
    RBox=1 ./Box

    close(traj)
    System=System_parameters(fname,frames,collect(1:frames),Number_Atoms,Number_Coils,Number_Chains,Box,RBox,Atom_Names,Residue_Names)
    return Positions,Chains,System
end

#This functions finds out which helices form a coiled-coil and returns a list of tuples
function Chain_pairings(Positions,Chains,System::System_parameters)
    Distance_Matrix=fill(Inf,(System.Number_Chains,System.Number_Chains))
    for (ChainA,ChainB) in combinations(Chains,2)
        Distance_Matrix[ChainA.ID,ChainB.ID]=sum(all_distances_frame(Positions[:,ChainA.CA,1],Positions[:,ChainB.CA,1],System.Box[:,1],System.RBox[:,1]))
    end
    Chain_pairings=[]
    for ChainA in Chains
        #global Coils,Chains
        if !any([Chain_pairings...] .== ChainA.ID )
            value,index=findmin(Distance_Matrix[ChainA.ID,:])
            append!(
            Chain_pairings,(ChainA.ID,index))
        end
    end
    return Chain_pairings
end

#No longer used
#function pbc_traj_mean(A,B,System::System_parameters)
#    N=size(A)[2]
#    box=reshape(repeat(System.Box,inner=(1,N)),(3,N,System.frames))
#    r_box=1 ./box
#    dx=A.-B
#    dx .-= floor.(Int64, dx .* r_box .+ 0.5) .* box
#    A=A-0.5*dx
#    return A
#end

#This function creates a direction Vector between the second and second to last bead of a coiled coil.
function pbc_traj_direction(A,B,System::System_parameters)
    N=size(A)[2]
    box=reshape(repeat(System.Box,inner=(1,N)),(3,N,System.frames))
    r_box=1 ./box
    dx=A.-B
    dx .-= floor.(Int64, dx .* r_box .+ 0.5) .* box
    absd=sqrt.(mapslices(sum,(dx.^2),dims=1))
    A=A-0.5*dx
    A=pbc_point(A,box)
    direction=A[:,end-1,:].-A[:,2,:]
    direction .-= floor.(Int64, direction .* System.RBox .+ 0.5) .* System.Box
    direction=direction./mapslices(norm,direction,dims=1)
    return A, direction,reshape(absd,(Number_CA,System.frames))
end

#This function creates the Coil_traj structures for each coiled-coil and chain_traj for every chain.
function create_Coils(Positions,Chains,System::System_parameters)
    Chains_Traj=Array{Chain_Traj}(undef,System.Number_Chains)
    for (i,chain) in enumerate(Chains)
        #global Positions
        CA=Positions[:,chain.CA,:]
        NH3=Positions[:,chain.NH3,:]
        COO=Positions[:,chain.COO,:]
        Chains_Traj[i]=Chain_Traj(Chains[i].ID,CA,NH3,COO)
    end
    pairings=Chain_pairings(Positions,Chains,System::System_parameters)
    Chains_Traj=Chains_Traj[pairings]
    coil_id=1
    Coils_Traj=Array{Coil_Traj}(undef,System.Number_Coils)
    for (chainA,chainB) in Iterators.partition(Chains_Traj,2)
        #global Box, RBox,coil_id
        Mean,Direction,bead_distance=pbc_traj_direction(chainA.CA,chainB.CA,System)
        Coils_Traj[coil_id]=Coil_Traj(coil_id,chainA,chainB,Mean,bead_distance,Direction)
        coil_id +=1
    end
    return Coils_Traj,Chains_Traj
end


#Meta function that creates the Coil_traj,Chain_traj and System structures of a given filename. Returns those and also saves them in JLD2 format.
function setup_trajectory(filename::String)
    fname=split(filename,".")[1]
    Positions,Chains,System=read_System("Traj/$filename",fname)
    Coils,Chains=create_Coils(Positions,Chains,System)
    mkpath("Output/")
    jldopen("Output/$(System.name).jld2","w+") do file
        file["Coils"]=Coils
        file["Chains"]=Chains
        file["System"]=System
    end
    return Coils,Chains,System
end

#Another function to calculate all distances between two set of positions. 
#This should work on while trajectories at once.
function distances_A_to_B(PositionsA,PositionsB,System::System_parameters)
    NA=size(PositionsA)[2]
    NB=size(PositionsB)[2]

    PositionsA=repeat(PositionsA,inner=(1,NB,1))
    PositionsB=repeat(PositionsB,outer=(1,NA,1))
    Box=reshape(repeat(System.Box,inner=(1,NA*NB)),(3,NA*NB,System.frames))
    RBox=1 ./Box

    distances=reshape(pbc_dist2(PositionsA,PositionsB,Box,RBox),(NB,NA,System.frames))
    return distances
end

#Function that checks wheter there are Hydrogenbonds between two given chains. It checks for each NH3-COO pair wheter there are any hydrogen bonds between all eleigle partners. Results are for whole trajectory
function calculate_chain_chain_hydrogen_bonds(ChainA::Chain_Traj,ChainB::Chain_Traj,System::System_parameters)
    distances=distances_A_to_B(ChainA.COO,ChainB.NH3,System)
    hydrogen_bonds=distances .< hydrogen_bond_cutoff
    hydrogen_bonds=hydrogen_bonds[1:3:end,:,:].+hydrogen_bonds[2:3:end,:,:].+hydrogen_bonds[3:3:end,:,:]
    hydrogen_bonds=hydrogen_bonds[:,1:2:end,:].+hydrogen_bonds[:,2:2:end,:]
    hydrogen_bonds=hydrogen_bonds .> 0
    return hydrogen_bonds
end


function extend_Chain(X::Int64,Matrix)
    N=size(Matrix)[1]
    NN=collect(1:N)
    X=[X]
    while !iszero(NN[Matrix[X[end],:]])
        if NN[Matrix[X[end],:]]... in X
            append!(X,NN[Matrix[X[end],:]])
            break
        else
            append!(X,NN[Matrix[X[end],:]])
        end
    end
    return X
end

function extend_Chain_split(X::Int64,Matrix)
    N=size(Matrix)[1]
    NN=collect(1:N)
    X=["0",[X]]
    for x in X
        if typeof(x)==String
            continue
        end
        while !iszero(NN[Matrix[x[end],:]])
            y=NN[Matrix[x[end],:]]
            if length(y)>1
                z=copy(x)
                push!(X,[z...,y[2]])
                append!(x,y[1])
            elseif NN[Matrix[x[end],:]][1] in x
                append!(x,NN[Matrix[x[end],:]])
                break
            else
                append!(x,NN[Matrix[x[end],:]])
            end
        end
    end
    return X[2:end]
end

function backwards_extend_Chain(X::Int64,Matrix::BitArray{2})
    N=size(Matrix)[1]
    NN=collect(1:N)
    X=[X]
    while !iszero(NN[Matrix[:,X[end]]])
        if NN[Matrix[:,X[end]]]... in X
            append!(X,NN[Matrix[:,X[end]]])
            break
        else
            append!(X,NN[Matrix[:,X[end]]])
        end
    end
    return X
end

function evaluate_connections_frame(Matrix)
    N=size(Matrix)[1]
    NN=collect(1:N)
    nn=copy(NN)
    Endings=NN[[sum(Matrix,dims=1)...] .== 0]
    Chain=[]
    Coils=[]
    for x in Endings
        chain=extend_Chain_split(x,Matrix)
        for z in chain
            if length(z)>1
                push!(Chain,z)
            end
            append!(Coils,z)
        end
    end
    pot_Circle=setdiff(NN[[sum(Matrix,dims=1)...] .== 1],unique(Coils))
    while !isempty(pot_Circle)
        try
            chain=extend_Chain(pot_Circle[1],Matrix)
            push!(Chain,chain)
            setdiff!(pot_Circle,chain)
            println("Looped Chain: ", chain)
        catch
            chain=backwards_extend_Chain(pot_Circle[end],Matrix)
            push!(Chain,chain[end:-1:1])
            setdiff!(pot_Circle,chain)
            println("Looped Chain: ", chain[end:-1:1])
        end
    end
    return Chain
end

function evaluate_connections(Matrix,System::System_parameters)
    Fragments=[]
    U_Fragments=[]
    for i in System.Frames
        fragments=evaluate_connections_frame(Matrix[:,:,i])
        push!(Fragments,fragments)
        append!(U_Fragments,fragments)
    end
    return Fragments,unique(U_Fragments),[length(x) for x in U_Fragments]
end

function connection_durations(Matrix::BitArray{3},System::System_parameters)
    Matrix=reshape(Matrix,System.Number_Coils^2,System.frames)
    sMatrix=[sum(Matrix,dims=2)...]
    bMatrix=Matrix[sMatrix.>0,:]
    Values=[]
    for i in 1:size(bMatrix)[1]
        k=0
        for frame in System.Frames
            if bMatrix[i,frame]==1
                k+=1
            else
                if k>0
                    append!(Values,k)
                    k=0
                end
            end
            if frame==System.frames
                if k>0
                    append!(Values,k)
                    k=0
                end
            end
        end
    end
    try
        jldopen("Output/$(System.name).jld2", "a+") do file
            file["Connection_duration"]=Values
            file["Mean_Connection_duration"]=mean(Values)

        end
    catch
        println("The file already has connection duration")
    end
end

#Create the dot products of a Array of vectors. The dot product is calculatet between each tuple in the table list.
#
function calculate_dot_products(vectors,table,n)
    F(X)=dot(vectors[:,X[1]],vectors[:,X[2]])
    cos_α=zeros(n,2)
    for i in table
        cos_α[i[2]-i[1],:]+=[F(i),1]
    end
    return cos_α[:,1]./cos_α[:,2]
end

#Uses least square fit to obtain a persistence length from the squared end to end distance R2. The fitted function is the Krohns-function.
function calculate_persistence_length_Krohn(R2,distance_max,initial_guess)
    R2_mean=mean(R2)
    Fit=curve_fit(Krohn2,[distance_max],[R2_mean],[initial_guess])
    return Fit.param[1 ]
end

#Finds the length at which the correlation function is 1/ℯ. If this is not inside tries to fit the points and interpolate the persistence length.
#If fitting is necessary it also means your chain is too short compared to your persistence length to get any good results.
#The fitting is quite bad.
function calculate_persistence_length_fit(correlation,distances,guess,N)
    N=floor(Int64,N/2)
    re=1/ℯ
    smaller_than_euler=correlation.<= re
    if any(smaller_than_euler)
        N2=findfirst(smaller_than_euler.<= re)
        N1=findlast(distances.>=1/ℯ)
        persistence_length=(distances[N2]-distances[N1])/(correlation[N2]-correlation[N1])*(re-correlation[N1])+distances[N1]
        return persistence_length,"Interp."
    else
        Fit=curve_fit(exp_fit,[distances[1:N]...],[correlation[1:N]...],[guess])
        return Fit.param[1],"Fit"
    end
end

#Centeres a chain so that its starting position is at the origin. The chain stays complete and does not consider periodic boundary conditions. Used for long bead chains to calculate the persistence length without creating error.
function centering(Positions,box)
    Start=Positions[:,1]
    Positions=pbc_vector2(Positions[:,2:end],Positions[:,1:end-1],box)
    Positions=hcat(Start,Positions)
    Positions=cumsum(Positions,dims=2)
    return Positions
end


function get_chain_beads(Coils,X,Fragments,Sequence,System)
    frames=System.Frames[Fragments[X,:]]
    Sequence=Sequence[X]
    Beads=zeros(Float64,(3,1,length(frames)))
    for coil in Sequence
        Beads=hcat(Beads,Coils[coil].Beads[:,:,frames])
    end
    Beads=Beads[:,2:end,:]
    for (i,l) in enumerate(frames)
        Beads[:,:,i]=centering(Beads[:,:,i],System.Box[:,l])
    end
    return Beads
end

#Calculate all possible intra chain hydrogen bonds of a given chain fot the whole trajectory.
function calculate_intra_chain_bonds(Chains::Vector{Chain_Traj},System::System_parameters)
    Intra_chain_bonds=zeros(Int64,(Number_Hydrogen_Bonds[1],Number_Hydrogen_Bonds[2],System.frames))
    Threads.@threads for Chain in Chains
        Intra_chain_bonds .+= calculate_chain_chain_hydrogen_bonds(Chain,Chain,System)
    end
    Total_intra_chain_bonds=sum(Intra_chain_bonds,dims=3)[:,:,1]
    try
        jldopen("Output/$(System.name).jld2", "a+") do file
            file["Intra_chain_bonds"]=Intra_chain_bonds
            file["Total_intra_chain_bonds"]=Total_intra_chain_bonds
        end
    catch
        println("The file already has intra chain bonds")
    end
end

#Calculate all possible inter chain hydrogen bonds of a given coiled-coil fot the whole trajectory.
function calculate_intra_coil_bonds(Coils::Vector{Coil_Traj},System::System_parameters)
    Intra_coil_bonds=zeros(Int64,(Number_Hydrogen_Bonds[1],Number_Hydrogen_Bonds[2],System.frames))
    Threads.@threads for Coil in Coils
        Intra_coil_bonds .+= calculate_chain_chain_hydrogen_bonds(Coil.ChainA,Coil.ChainB,System)
        Intra_coil_bonds .+= calculate_chain_chain_hydrogen_bonds(Coil.ChainB,Coil.ChainA,System)
    end
    Total_intra_coil_bonds=sum(Intra_coil_bonds,dims=3)[:,:,1]
    try
        jldopen("Output/$(System.name).jld2", "a+") do file
            file["Intra_coil_bonds"]=Intra_coil_bonds
            file["Total_intra_coil_bonds"]=Total_intra_coil_bonds
        end
    catch
        println("The file already has intra coil bonds")
    end
end

#Calculate all possible hydrogen bonds between chains that belong to two different coiled-coils for the whole trajectory.
function calculate_inter_coil_bonds(Coils::Vector{Coil_Traj},System::System_parameters)
    Inter_coil_bonds=zeros(Int64,(System.Number_Coils,System.Number_Coils,Number_Hydrogen_Bonds[1],Number_Hydrogen_Bonds[2],System.frames))
    n=binomial(System.Number_Coils,2)
    k=0
    Threads.@threads for (CoilA,CoilB) in collect(combinations(Coils,2))
        Inter_coil_bonds[CoilA.ID,CoilB.ID,:,:,:] .+= calculate_chain_chain_hydrogen_bonds(CoilA.ChainA,CoilB.ChainA,System)
        Inter_coil_bonds[CoilA.ID,CoilB.ID,:,:,:] .+= calculate_chain_chain_hydrogen_bonds(CoilA.ChainA,CoilB.ChainB,System)
        Inter_coil_bonds[CoilA.ID,CoilB.ID,:,:,:] .+= calculate_chain_chain_hydrogen_bonds(CoilA.ChainB,CoilB.ChainA,System)
        Inter_coil_bonds[CoilA.ID,CoilB.ID,:,:,:] .+= calculate_chain_chain_hydrogen_bonds(CoilA.ChainB,CoilB.ChainB,System)

        Inter_coil_bonds[CoilB.ID,CoilA.ID,:,:,:] .+= calculate_chain_chain_hydrogen_bonds(CoilB.ChainA,CoilA.ChainA,System)
        Inter_coil_bonds[CoilB.ID,CoilA.ID,:,:,:] .+= calculate_chain_chain_hydrogen_bonds(CoilB.ChainA,CoilA.ChainB,System)
        Inter_coil_bonds[CoilB.ID,CoilA.ID,:,:,:] .+= calculate_chain_chain_hydrogen_bonds(CoilB.ChainB,CoilA.ChainA,System)
        Inter_coil_bonds[CoilB.ID,CoilA.ID,:,:,:] .+= calculate_chain_chain_hydrogen_bonds(CoilB.ChainB,CoilA.ChainB,System)
        k+=1
        if (k/n*100)%10==0
            flush(stdout)
            println(k/n*100," % Progress")
        end
    end
    Total_inter_coil_bonds=sum(Inter_coil_bonds,dims=(1,2,5))[1,1,:,:,1]
    Any_bonds=sum(Inter_coil_bonds,dims=(3,4))[:,:,1,1,:].>0
    Chain_bonds=sum(Inter_coil_bonds[:,:,1:3,4:5,:],dims=(3,4))[:,:,1,1,:].>0
    try
        jldopen("Output/$(System.name).jld2", "a+") do file
            file["Inter_coil_bonds"]=Inter_coil_bonds
            file["Total_inter_coil_bonds"]=Total_inter_coil_bonds
            file["Any_bonds"]=Any_bonds
            file["Chain_bonds"]=Chain_bonds
        end
    catch
        println("The file already has inter coil bonds")
    end
    return Chain_bonds
end

function calculate_bead_chains(Coils,Chain_bonds,System)
    Chain_fragments,Unique_chain_fragments,Length_Unique=evaluate_connections(Chain_bonds,System)

    Fragments_per_frame=zeros(Bool,(length(Unique_chain_fragments),System.frames))
    for (index,frame) in enumerate(Chain_fragments)
    Fragments_per_frame[:,index]=[x in frame for x in Unique_chain_fragments]
    end

    Bead_chains=[]
    Bead_frames=[]
    for len in sort(unique(Length_Unique))
        bead_chains=zeros(Float64,(3,len*8,1))
        bead_frames=[]
        for k in 1:length(Unique_chain_fragments)
            if length(Unique_chain_fragments[k])==len
                append!(bead_frames,System.Frames[Fragments_per_frame[k,:]])
                bead_chains=cat(bead_chains,get_chain_beads(Coils,k,Fragments_per_frame,Unique_chain_fragments,System),dims=3)
            end
        end
        push!(Bead_chains,bead_chains[:,:,2:end])
        push!(Bead_frames,bead_frames)
    end
    try
        jldopen("Output/$(System.name).jld2", "a+") do file
            file["Chain_fragments"]=Chain_fragments
            file["Unique_chain_fragments"]=Unique_chain_fragments
            file["Length_Unique"]=Length_Unique
            file["Fragments_per_frame"]=Fragments_per_frame
            file["Bead_chains"]=Bead_chains
            file["Bead_frames"]=Bead_frames
        end
    catch
        println("The file already has bead chains")
    end
    return Bead_chains,Bead_frames,sort(unique(Length_Unique))
end

function calculate_persitence_length_chain(pos,frames,System::System_parameters)
    nbeads=size(pos)[2]
    nframe=size(pos)[3]
    nvector=nbeads-1
    ndotvector=nbeads-2

    table_full=collect(combinations(collect(1:nvector),2))       ## All Combinations for averaging over the chain positions and time
    table_single=table_full[1:ndotvector]

    R2=sum((pos[:,1,:].-pos[:,end,:]).^2,dims=1)

    dir_vector=pbc_traj_vec(pos[:,2:end,:],pos[:,1:end-1,:],System.Box[:,frames])
    dir_distance=vector_length(dir_vector)
    normvector= dir_vector./reshape(repeat(dir_distance,inner=(3,1)),(3,nvector,nframe))

    correlations=zeros(Float64,(ndotvector))
    correlations_s=zeros(Float64,(ndotvector))
    for i in 1:nframe
        correlations=hcat(correlations,calculate_dot_products(normvector[:,:,i],table_full,ndotvector))
        correlations_s=hcat(correlations_s,calculate_dot_products(normvector[:,:,i],table_single,ndotvector))
    end

    correlations=correlations[:,2:end]
    mcorr=mean(correlations,dims=2)[:,1]
    Persistence_length_summed=sum(mean(correlations,dims=2)[:,1])*mean(dir_distance)

    correlations_s=correlations_s[:,2:end]
    mcorr_s=mean(correlations_s,dims=2)[:,1]
    LP_s=sum(mean(correlations_s,dims=2)[:,1])*mean(dir_distance)
    Mean_Distances=Mean_Distances=mean(dir_distance,dims=2)
    Cum_Distances=[cumsum(Mean_Distances,dims=1)...]
    Contour_Length=Cum_Distances[end]

    Persistence_length_full_fit,mode_full=calculate_persistence_length_fit(mcorr,Cum_Distances[1:end-1],Persistence_length_summed,ndotvector)
    Persistence_length_Krohn=calculate_persistence_length_Krohn(R2[2:end],Contour_Length,LP_s)
    return Persistence_length_summed, Persistence_length_full_fit, Persistence_length_Krohn,mcorr,Mean_Distances[:,1]
end


function evaluate_trajectory(name::String)
    flush(stdout)
    println("Start with $name")
    @time Coils,Chains,System=setup_trajectory(name)
    flush(stdout)
    println("Finished reading and setting up system")
    @time calculate_intra_chain_bonds(Chains,System)
    flush(stdout)
    println("Finished calculating intra chain hydrogen bonds")
    @time calculate_intra_coil_bonds(Coils,System)
    flush(stdout)
    println("Finished calculating intra coil hydrogen bonds")
    @time Chain_bonds=calculate_inter_coil_bonds(Coils,System)
    flush(stdout)
    println("Finished calculating inter coil hydrogen bonds")
    @time connection_durations(Chain_bonds,System)
    flush(stdout)
    println("Finished calculating inter coil connection durations")
    @time Bead_chains,Bead_frames,Chain_length=calculate_bead_chains(Coils,Chain_bonds,System)
    flush(stdout)
    println("Finished calculating bead chains")

    N=size(Bead_chains)[1]
    Results=zeros(Float64,(N,5))
    Mean_correlations=Array{Vector{Float64},1}(undef,N)
    Mean_distances=Array{Vector{Float64},1}(undef,N)
    @time Threads.@threads for I in collect(1:N)
        try
            persistence_length_summed, persistence_length_full_fit, persistence_length_Krohn,mean_correlations,mean_distances=calculate_persitence_length_chain(Bead_chains[I],Bead_frames[I],System)
            Results[I,:]=[Chain_length[I]...,persistence_length_summed..., persistence_length_full_fit..., persistence_length_Krohn..., size(Bead_chains[I])[3]]
            Mean_correlations[I]=mean_correlations
            Mean_distances[I]=mean_distances
        catch
            Results[I,end]=size(Bead_chains[I])[3]
            println("Chains $I has a problem. The reason may be that it only has $(size(Bead_chains[I])[3]) frames")
        end
    end
    try
        jldopen("Output/$(System.name).jld2", "a+") do file
            file["Results"]=Results
            file["Mean_correlations"]=Mean_correlations
            file["Mean_distances"]=Mean_distances
        end
    catch
        println("The file already has persistence lengths")
    end
    flush(stdout)
    println("Finished calculating persistence lengths.")
    println("Finished with $(System.name)")
end


#@time evaluate_trajectory("NO.gro")
#@time evaluate_trajectory("PAB.gro")
@time evaluate_trajectory("OAB.gro")


GC.gc()
flush(stdout)
println("Finish")


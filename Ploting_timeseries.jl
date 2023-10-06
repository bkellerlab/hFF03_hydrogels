using Chemfiles
using LinearAlgebra
using StaticArrays
using PyCall
using Printf
using StatsBase
using Statistics
using JLD2


flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Start")
GC.gc()
plt = pyimport("matplotlib.pyplot")
cm = pyimport("matplotlib.cm")
pltc = pyimport("matplotlib.colors")
pltl = pyimport("matplotlib.lines")
pygui()

const Number_CA = 8
const Number_NH3 = 27
const Number_COO = 10

const hydrogen_bond_cutoff = 3.5
const running_number = 10

struct Chain
    ID::Int64
    CA::SVector{Number_CA,Int64}
    NH3::SVector{Number_NH3,Int64}
    COO::SVector{Number_COO,Int64}
end

struct Chain_Traj
    ID::Int64
    CA::Array{Float64,3}
    NH3::Array{Float64,3}
    COO::Array{Float64,3}
end

struct Coil_Traj
    ID::Int64
    ChainA::Chain_Traj
    ChainB::Chain_Traj
    Beads::Array{Float64,3}
    Bead_distance::Array{Float64,2}
    Directions::Array{Float64,2}
end

struct System_parameters
    name::String
    frames::Int64
    Frames::Array{Int64,1}
    Number_Atoms::Int64
    Number_Coils::Int64
    Number_Chains::Int64
    Box::Array{Float64,2}
    RBox::Array{Float64,2}
    Atom_Names::Array{String,1}
    Residue_Names::Array{String,1}
end

function block_not_bonded(X,blocksize::Int64)
    any_bonds=X["Any_bonds"]
    nc=X["System"].Number_Coils
    nf=X["System"].frames
    bonded=(reshape(sum(any_bonds,dims=1),(nc,nf)).+reshape(sum(any_bonds,dims=2),(nc,nf)) ).>0
    blocks=floor(Int64,nf/blocksize)
    not=nc.-[sum(bonded,dims=1)...]
    block_not=zeros(blocks)
    for x in 1:blocks
        start=(x-1)*blocksize+1
        ending=x*blocksize
        block_not[x]=floor(Int64,mean(not[start:ending]))
    end
    return block_not
end

function block_mean_coils(X,blocksize)
    values=X["Fragments_per_frame"]
    frames=X["System"].frames
    uniques=[length(x) for x in X["Unique_chain_fragments"]]
    mean_x=zeros(frames)
    std_x=zeros(frames)
    for x in 1:frames
        value=values[:,x].*uniques
        value=value[value.>0]
        mean_x[x]=mean(value)
        std_x[x]=std(value)
    end
    blocks=floor(Int64,frames/blocksize)
    block_mean=zeros(blocks)
    for x in 1:blocks
        start=(x-1)*blocksize+1
        ending=x*blocksize
        block_mean[x]=mean(mean_x[start:ending])
    end
    return block_mean
end

function create_means(name1,name2,name3,blocksize)
    X1 = load("Output/$name1.jld2")
    X2 = load("Output/$name2.jld2")
    X3 = load("Output/$name3.jld2")

    X_not=hcat( block_not_bonded(X1,blocksize),block_not_bonded(X2,blocksize),block_not_bonded(X3,blocksize))
    X_not_mean=[mean(X_not,dims=2)...]
    X_not_std=[std(X_not,dims=2)...]

    X_coils=hcat( block_mean_coils(X1,blocksize),block_mean_coils(X2,blocksize),block_mean_coils(X3,blocksize))
    X_coils_mean=[mean(X_coils,dims=2)...]
    X_coils_std=[std(X_coils,dims=2)...]
    return X_not_mean, X_not_std, X_coils_mean,X_coils_std
end
N_not_mean,N_not_std,N_coils_mean,N_coils_std=create_means("NO_32_1","NO_32_2","NO_32_3",50)
O_not_mean,O_not_std,O_coils_mean,O_coils_std=create_means("OABA_32_1","OABA_32_2","OABA_32_3",50)
P_not_mean,P_not_std,P_coils_mean,P_coils_std=create_means("PABA_32_1","PABA_32_2","PABA_32_3",50)
GC.gc()
times=collect(1:150)



fig, ax = plt.subplots()
ax2=ax.twinx()
ax.plot(times,N_not_mean,color="blue",label="no-hFF03")
ax.fill_between(times,N_not_mean.-N_not_std,N_not_mean.+N_not_std,color="blue",alpha=0.3)
ax2.plot(times,N_coils_mean,color="blue",linestyle="dashed")
ax2.fill_between(times,N_coils_mean.-N_coils_std,N_coils_mean.+N_coils_std,color="blue",alpha=0.2)

ax.plot(times,O_not_mean,color="red",label="oaba-hFF03")
ax.fill_between(times,O_not_mean.-O_not_std,O_not_mean.+O_not_std,color="red",alpha=0.3)
ax2.plot(times,O_coils_mean,color="red",linestyle="dashed")
ax2.fill_between(times,O_coils_mean.-O_coils_std,O_coils_mean.+O_coils_std,color="red",alpha=0.2)

ax.plot(times,P_not_mean,color="green",label="paba-hFF03")
ax.fill_between(times,P_not_mean.-P_not_std,P_not_mean.+P_not_std,color="green",alpha=0.3)
ax2.plot(times,P_coils_mean,color="green",linestyle="dashed")
ax2.fill_between(times,P_coils_mean.-P_coils_std,P_coils_mean.+P_coils_std,color="green",alpha=0.2)

ax.set_xlim(0,150)
ax.set_ylim(0,26)
ax.set_ylabel("unbound coiled-coils",fontsize=15)
ax.set_xlabel("t in ns",fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=16)

ax2.set_xlim(0,150)
ax2.set_ylim(0,10)
ax2.set_ylabel("coiled-coils per fibril",fontsize=15)
ax2.tick_params(axis="both", which="major", labelsize=16)

ax2.plot(0,0,color="blue",label="no-hFF03")
ax2.plot(0,0,color="red",label="oaba-hFF03")
ax2.plot(0,0,color="green",label="paba-hFF03")
ax2.plot(0,0,color="grey",label="mean unbound")
ax2.plot(0,0,color="grey",label="mean per fibril",linestyle="dashed")
ax2.fill(0,0,color="grey",label="standard deviation",alpha=0.2)

ax2.legend(loc="upper right",fontsize=15,ncol=2, columnspacing=0.5)
plt.tight_layout()
plt.savefig("Output/timeseries.pdf")
plt.close()


GC.gc()
flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Finish")


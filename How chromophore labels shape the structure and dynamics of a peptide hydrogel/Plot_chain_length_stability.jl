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

N1 = load("JLD/NO_32_end_1.jld2")
N2 = load("JLD/NO_32_end_2.jld2")
N3 = load("JLD/NO_32_end_3.jld2")

O1 = load("JLD/OABA_32_end_1.jld2")
O2 = load("JLD/OABA_32_end_2.jld2")
O3 = load("JLD/OABA_32_end_3.jld2")

P1 = load("JLD/PABA_32_end_1.jld2")
P2 = load("JLD/PABA_32_end_2.jld2")
P3 = load("JLD/PABA_32_end_3.jld2")


function not_bonded(X)
    any_bonds = X["Any_bonds"]
    nc = X["System"].Number_Coils
    nf = X["System"].frames
    bonded = (reshape(sum(any_bonds, dims=1), (nc, nf)) .+ reshape(sum(any_bonds, dims=2), (nc, nf))) .> 0
    not = nc .- [sum(bonded, dims=1)...]
    return not
end

function fiber_chain_distribution(X)
    values = X["Fragments_per_frame"]
    uniques = [length(x) for x in X["Unique_chain_fragments"]]
    coils = values .* uniques
    sorted = [1, sort(unique(sort(uniques)))...]
    summed = zeros(size(sorted))
    for (i, x) in enumerate(sorted)
        summed[i] = sum(coils .== x)
    end
    summed[1] = sum(not_bonded(X))
    p = summed ./ sum(summed) * 100
    return sorted, p
end



N_length_1,N_P_1=fiber_chain_distribution(N1)
N_length_2,N_P_2=fiber_chain_distribution(N2)
N_length_3,N_P_3=fiber_chain_distribution(N3)

N_length=N_length_1
N_P_3=[N_P_3...,0.0,0.0,0.0,0.0]

N_P=hcat(N_P_1,N_P_2,N_P_3)
N_P=hcat(N_P_2,N_P_3)
N_P_mean=[mean(N_P,dims=2)...]
N_P_std=[std(N_P,dims=2)...]


O_length_1,O_P_1=fiber_chain_distribution(O1)
O_length_2,O_P_2=fiber_chain_distribution(O2)
O_length_3,O_P_3=fiber_chain_distribution(O3)

O_length=O_length_1[1:end-1]
O_P_1=O_P_1[1:end-1]
O_P_2=[O_P_2...,0.0]
O_P_3=[O_P_3...,0.0,0.0,0.0]

O_P=hcat(O_P_1,O_P_2,O_P_3)
O_P=hcat(O_P_2,O_P_3)
O_P_mean=[mean(O_P,dims=2)...]
O_P_std=[std(O_P,dims=2)...]



P_length_1,P_P_1=fiber_chain_distribution(P1)
P_length_2,P_P_2=fiber_chain_distribution(P2)
P_length_3,P_P_3=fiber_chain_distribution(P3)

P_length=[1,2,3,4,5,6,7,8,9,10]

P_P_1=[sum(P_P_1[1:2]),P_P_1[3:end]...,0.0,0.0]
P_P_2=[P_P_2...,0.0]
P_P_3=[sum(P_P_3[1:2]),P_P_3[3:end]...]

P_P=hcat(P_P_1,P_P_2,P_P_3)
P_P=hcat(P_P_2,P_P_3)
P_P_mean=[mean(P_P,dims=2)...]
P_P_std=[std(P_P,dims=2)...]


fig, ax = plt.subplots()
ax.bar(N_length .- 0.2, N_P_mean,color="blue", width=0.2, label="no-hFF03")
ax.errorbar(N_length .- 0.2, N_P_mean,N_P_std, color="black",elinewidth=0.6,capsize=2,fmt=" ",markeredgewidth=0.6)
ax.bar(O_length, O_P_mean, color="red", width=0.2, label="oaba-hFF03")
ax.errorbar(O_length, O_P_mean,O_P_std, color="black",elinewidth=0.6,capsize=2,fmt=" ",markeredgewidth=0.6)
ax.bar(P_length .+ 0.2, P_P_mean, color="green", width=0.2, label="paba-hFF03")
ax.errorbar(P_length .+ 0.2, P_P_mean,P_P_std, color="black",elinewidth=0.6,capsize=2,fmt=" ",markeredgewidth=0.6)
ax.legend(fontsize=18)
ax.set_xlim(0, 13)
ax.set_xticks([1,2,3,4,5,6,7,8,9,10,11,12])
ax.set_ylim(0, 70)
ax.set_ylabel("P in %", fontsize=18)
ax.set_xlabel("Number of coiled-coils in fiber chain", fontsize=18)
ax.set_xticks(collect(1:13))
ax.tick_params(labelsize=18)
plt.tight_layout()
plt.savefig("Fiber_distribution2.pdf")
plt.close()


###############################################################################################################
###############################################################################################################
###############################################################################################################

N_d_1=N1["Connection_duration"]
N_d_2=N2["Connection_duration"]
N_d_3=N3["Connection_duration"]
Nd=vcat(N_d_1,N_d_2,N_d_3)
Nd_mean=mean(Nd)*0.02
Nd_std=std(Nd)*0.02

Nd_mean2=mean([mean(N_d_1),mean(N_d_2),mean(N_d_3)])*0.02
Nd_std2=std([mean(N_d_1),mean(N_d_2),mean(N_d_3)])*0.02


O_d_1=O1["Connection_duration"]
O_d_2=O2["Connection_duration"]
O_d_3=O3["Connection_duration"]
Od=vcat(O_d_1,O_d_2,O_d_3)
Od_mean=mean(Od)*0.02
Od_std=std(Od)*0.02

Od_mean2=mean([mean(O_d_1),mean(O_d_2),mean(O_d_3)])*0.02
Od_std2=std([mean(O_d_1),mean(O_d_2),mean(N_d_3)])*0.02


P_d_1=P1["Connection_duration"]
P_d_2=P2["Connection_duration"]
P_d_3=P3["Connection_duration"]
Pd=vcat(P_d_1,P_d_2,P_d_3)
Pd_mean=mean(Pd)*0.02
Pd_std=std(Pd)*0.02

Pd_mean2=mean([mean(P_d_1),mean(P_d_2),mean(P_d_3)])*0.02
Pd_std2=std([mean(P_d_1),mean(P_d_2),mean(P_d_3)])*0.02


using FHist
RR=collect(1:1:2501)

nH1=Hist1D(N_d_1,RR,overflow=true)|> normalize
nH2=Hist1D(N_d_2,RR,overflow=true)|> normalize
nH3=Hist1D(N_d_3,RR,overflow=true)|> normalize
nH=hcat(nH1.hist.weights,nH2.hist.weights,nH3.hist.weights)*100
nH_mean=[mean(nH,dims=2)...]
nH_std=[std(nH,dims=2)...]

oH1=Hist1D(O_d_1,RR,overflow=true)|> normalize
oH2=Hist1D(O_d_2,RR,overflow=true)|> normalize
oH3=Hist1D(O_d_3,RR,overflow=true)|> normalize
oH=hcat(oH1.hist.weights,oH2.hist.weights,oH3.hist.weights)*100
oH_mean=[mean(oH,dims=2)...]
oH_std=[std(oH,dims=2)...]


pH1=Hist1D(P_d_1,RR,overflow=true)|> normalize
pH2=Hist1D(P_d_2,RR,overflow=true)|> normalize
pH3=Hist1D(P_d_3,RR,overflow=true)|> normalize
pH=hcat(pH1.hist.weights,pH2.hist.weights,pH3.hist.weights)*100
pH_mean=[mean(pH,dims=2)...]
pH_std=[std(pH,dims=2)...]

xval=RR[1:end-1].*0.02

fig, ax=plt.subplots()
ax.step(xval,nH_mean,color="blue",label="no-hFF03",where="post")
ax.step(xval,cumsum(nH_mean),color="blue",where="post",linestyle="--")
ax.fill_between(xval, nH_mean .- nH_std, nH_mean .+ nH_std, color="blue", alpha=0.2,step="post")
ax.step(xval,oH_mean,color="red",label="oaba-hFF03",where="post")
ax.step(xval,cumsum(oH_mean),color="red",where="post",linestyle="--")
ax.fill_between(xval, oH_mean .- oH_std, oH_mean .+ oH_std, color="red", alpha=0.2,step="post")
ax.step(xval,pH_mean,color="green",label="paba-hFF03",where="post")
ax.step(xval,cumsum(pH_mean),color="green",where="post",linestyle="--")
ax.fill_between(xval, pH_mean .- pH_std, pH_mean .+ pH_std, color="green", alpha=0.2,step="post")
ax.set_xlim(0.02,51)
ax.set_ylim(0,101)
ax.set_xscale("log")
ax.set_xticks(ticks=[0.02,0.1,1,10,50],labels=[0.02,0.1,1,10,50])
ax.set_xlabel("Connection Time in ns",fontsize=16)
ax.set_ylabel("P in %",fontsize=16)
ax.tick_params(labelsize=16) 
ax.legend(loc="center right", fontsize=16,)
plt.tight_layout()
plt.savefig("connection_time.pdf")
plt.close()

############################################################################################################
############################################################################################################
############################################################################################################



function calculate_Salt_bridges(X,N)
    Intra_chain = reshape(sum(X["Intra_chain_bonds"], dims=3), (9, 5))
    Intra_coil  = reshape(sum(X["Intra_coil_bonds"], dims=3), (9, 5))
    Inter_coil = reshape(sum(X["Inter_coil_bonds"], dims=(1, 2,5)), (9, 5))
    if Inter_coil[2,4]>2*N
        Inter_coil[2,4]=Inter_coil[2,4]-maximum(X["Inter_coil_bonds"])
    end
    p_Intra_chain=sum(Intra_chain[1, :]) / N * 100
    p_Intra_coil=sum(Intra_coil[1, :]) / N * 100
    p_Inter_coil=sum(Inter_coil[1, :]) / N * 100
    p_Water=100-p_Intra_chain-p_Intra_coil-p_Inter_coil
    return [p_Intra_chain,p_Intra_coil,p_Inter_coil,p_Water]
end

NX=64*2501
N_Salt=hcat(calculate_Salt_bridges(N1,NX),calculate_Salt_bridges(N2,NX),calculate_Salt_bridges(N3,NX))
N_Salt_mean=mean(N_Salt,dims=2)
N_Salt_std=std(N_Salt,dims=2)

O_Salt=hcat(calculate_Salt_bridges(O1,NX),calculate_Salt_bridges(O2,NX),calculate_Salt_bridges(O3,NX))
O_Salt_mean=mean(O_Salt,dims=2)
O_Salt_std=std(O_Salt,dims=2)

P_Salt=hcat(calculate_Salt_bridges(P1,NX),calculate_Salt_bridges(P2,NX),calculate_Salt_bridges(P3,NX))
P_Salt_mean=mean(P_Salt,dims=2)
P_Salt_std=std(P_Salt,dims=2)

GC.gc()
flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Finish")


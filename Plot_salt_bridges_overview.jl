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
plt=pyimport("matplotlib.pyplot")
cm=pyimport("matplotlib.cm")
pltc=pyimport("matplotlib.colors")
pltl=pyimport("matplotlib.lines")
pygui()

const Number_CA=8
const Number_NH3=27
const Number_COO=10

const hydrogen_bond_cutoff=3.5
const running_number=10

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

Y_names=["N-Ter.","K2","K3","K9","K10","K16","K17","K23","K24"]
X_names=["E4","E11","E18","E25","C-Ter."]

function make_salt_bridge_plot(filename,typename,plotname,limit,bool)
    file_=load("JLD/$filename.jld2")
    NFrames=file_["System"].frames
    if bool
        N_=64
    else
        N_=32
    end

    traj=file_[typename]
    bonds = round.(traj./(N_*NFrames) .* 100,digits=3)
    max= ceil(maximum(bonds))
    fig,ax=plt.subplots()
    for x in 1:9
        for y in 1:5
            if traj[x, y] !=0.0
                plt.text(y - 0.5, x - 0.5, traj[x, y],horizontalalignment="center",verticalalignment="center",color="w", fontsize=16)
            end
        end
    end
    norm_=plt.Normalize(0.01,max)
    cmap = pltc.LinearSegmentedColormap.from_list("", ["white","orange","orangered","red","darkviolet","mediumblue"])
    cmap.set_under("white")
    cmap.set_over("midnightblue")
    cax=plt.pcolormesh(bonds, cmap=cmap,edgecolors="k", linewidth=0.1,norm=norm_)
    plt.clim(0,limit)
    cbar=fig.colorbar(cax,ticks=collect(0:5:limit))
    cbar.set_label("Salt bridge chance per frame in %", rotation=90, fontsize=16,y=0.45)
    cbar.ax.tick_params(labelsize=16) 
    ax.set_xticks(collect(0.5:1:5))
    ax.set_xticklabels(X_names,rotation=90, fontsize=16)
    ax.set_yticks(collect(0.5:1:9))
    ax.set_yticklabels(Y_names, fontsize=16)
    ax.set_aspect(0.5)
    plt.tight_layout()
    plt.savefig("Salt/$plotname.pdf")
    plt.close()
end

# make_salt_bridge_plot("NO_32_end_1","Total_intra_chain_bonds","Intra_chain_NO_1",60,true)
# make_salt_bridge_plot("NO_32_end_2","Total_intra_chain_bonds","Intra_chain_NO_2",60,true)
# make_salt_bridge_plot("NO_32_end_3","Total_intra_chain_bonds","Intra_chain_NO_3",60,true)

# make_salt_bridge_plot("OABA_32_end_1","Total_intra_chain_bonds","Intra_chain_OABA_1",60,true)
# make_salt_bridge_plot("OABA_32_end_2","Total_intra_chain_bonds","Intra_chain_OABA_2",60,true)
# make_salt_bridge_plot("OABA_32_end_3","Total_intra_chain_bonds","Intra_chain_OABA_3",60,true)

# make_salt_bridge_plot("PABA_32_end_1","Total_intra_chain_bonds","Intra_chain_PABA_1",60,true)
# make_salt_bridge_plot("PABA_32_end_2","Total_intra_chain_bonds","Intra_chain_PABA_2",60,true)
# make_salt_bridge_plot("PABA_32_end_3","Total_intra_chain_bonds","Intra_coil_PABA_3",60,true)



# make_salt_bridge_plot("NO_32_end_1","Total_intra_coil_bonds","Intra_coil_NO_1",30,true)
# make_salt_bridge_plot("NO_32_end_2","Total_intra_coil_bonds","Intra_coil_NO_2",30,true)
# make_salt_bridge_plot("NO_32_end_3","Total_intra_coil_bonds","Intra_coil_NO_3",30,true)

# make_salt_bridge_plot("OABA_32_end_1","Total_intra_coil_bonds","Intra_coil_OABA_1",30,true)
# make_salt_bridge_plot("OABA_32_end_2","Total_intra_coil_bonds","Intra_coil_OABA_2",30,true)
# make_salt_bridge_plot("OABA_32_end_3","Total_intra_coil_bonds","Intra_coil_OABA_3",30,true)

# make_salt_bridge_plot("PABA_32_end_1","Total_intra_coil_bonds","Intra_coil_PABA_1",30,true)
# make_salt_bridge_plot("PABA_32_end_2","Total_intra_coil_bonds","Intra_coil_PABA_2",30,true)
# make_salt_bridge_plot("PABA_32_end_3","Total_intra_coil_bonds","Intra_coil_PABA_3",30,true)


# make_salt_bridge_plot("NO_32_end_1","Total_inter_coil_bonds","Inter_coil_NO_1",60,false)
# make_salt_bridge_plot("NO_32_end_2","Total_inter_coil_bonds","Inter_coil_NO_2",60,false)
# make_salt_bridge_plot("NO_32_end_3","Total_inter_coil_bonds","Inter_coil_NO_3",60,false)

# make_salt_bridge_plot("OABA_32_end_1","Total_inter_coil_bonds","Inter_coil_OABA_1",60,false)
# make_salt_bridge_plot("OABA_32_end_2","Total_inter_coil_bonds","Inter_coil_OABA_2",60,false)
# make_salt_bridge_plot("OABA_32_end_3","Total_inter_coil_bonds","Inter_coil_OABA_3",60,false)

# make_salt_bridge_plot("PABA_32_end_1","Total_inter_coil_bonds","Inter_coil_PABA_1",60,false)
# make_salt_bridge_plot("PABA_32_end_2","Total_inter_coil_bonds","Inter_coil_PABA_2",60,false)
# make_salt_bridge_plot("PABA_32_end_3","Total_inter_coil_bonds","Inter_coil_PABA_3",60,false)

function make_salt_bridge_plot_frames(filename,typename,plotname,limit,bool,start_frame)
    file_=load("JLD/$filename.jld2")
    NFrames=file_["System"].frames-start_frame
    if bool
        N_=64
    else
        N_=32
    end

    traj=file_[typename]
    traj=reshape(sum(traj[:,:,:,:,start_frame:end],dims=(1,2,5)),(9,5))
    bonds = round.(traj./(N_*NFrames) .* 100,digits=3)
    max= ceil(maximum(bonds))
    fig,ax=plt.subplots()
    for x in 1:9
        for y in 1:5
            if traj[x, y] !=0.0
                plt.text(y - 0.5, x - 0.5, traj[x, y],horizontalalignment="center",verticalalignment="center",color="w", fontsize=16)
            end
        end
    end
    norm_=plt.Normalize(0.01,max)
    cmap = pltc.LinearSegmentedColormap.from_list("", ["white","orange","orangered","red","darkviolet","mediumblue"])
    cmap.set_under("white")
    cmap.set_over("midnightblue")
    cax=plt.pcolormesh(bonds, cmap=cmap,edgecolors="k", linewidth=0.1,norm=norm_)
    plt.clim(0,limit)
    cbar=fig.colorbar(cax,ticks=collect(0:5:limit))
    cbar.set_label("Salt bridge chance per frame in %", rotation=90, fontsize=16,y=0.45)
    cbar.ax.tick_params(labelsize=16) 
    ax.set_xticks(collect(0.5:1:5))
    ax.set_xticklabels(X_names,rotation=90, fontsize=16)
    ax.set_yticks(collect(0.5:1:9))
    ax.set_yticklabels(Y_names, fontsize=16)
    ax.set_aspect(0.5)
    plt.tight_layout()
    plt.savefig("Salt/$plotname.pdf")
    plt.close()
end

make_salt_bridge_plot_frames("OABA_32_2","Inter_coil_bonds","Inter_coil_OABA_2",60,false,5000)

make_salt_bridge_plot_frames("NO_32_1","Inter_coil_bonds","Inter_coil_NO_1",60,false,5000)
make_salt_bridge_plot_frames("NO_32_3","Inter_coil_bonds","Inter_coil_NO_3",60,false,5000)

GC.gc()
flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Finish")
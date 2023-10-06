using Chemfiles
using LinearAlgebra
using StaticArrays
using Combinatorics
using BenchmarkTools
using PyCall
using RollingFunctions
using CoordinateTransformations
using FFTW
using StatsBase

flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Start")
GC.gc()
plt=pyimport("matplotlib.pyplot")
pygui()

const start_water=1804
const m_O=15.999
const m_H=1.00784
const m_H20=m_O+2m_H

function center_of_mass(values)
    com=(m_O*values[:,1:3:end] .+m_H*values[:,2:3:end].+m_H*values[:,2:3:end])/m_H20
    return com
end

function pbc_mean(A,B,box,r_box)
    dx=A.-B
    dx .-= floor.(Int64, dx .* r_box .+ 0.5) .* box
    A=A-0.5*dx
    return A
end





filename="W_KF.trr"
const traj=Trajectory(filename)
const frames=Int64(length(traj))
topologyname="W_KF.gro"

top_traj=Trajectory(topologyname)
top_frame=read_step(top_traj,0)
top=Topology(top_frame)

Number_Atoms=length(top_frame)
Residues_Names=Array{String}(undef,Number_Atoms)
Atom_Names=Array{String}(undef,Number_Atoms)
for i in 1:Number_Atoms
    global Atom_Names, Residue_Names,Mass
    Residues_Names[i]=name(residue_for_atom(top,i-1))
    Atom_Names[i]=name(Chemfiles.Atom(top_frame,i-1))
end

const mask_water= (Atom_Names .== "OW") .|  (Atom_Names .== "HW1") .| (Atom_Names .== "HW2")

const N_water=Int64(sum(mask_water)/3)

Velocities=zeros(3,N_water,frames)

i=0
@time while i < frames
    global i,Velocities
    V=velocities(read_step(traj,i))
    V=center_of_mass(V[:,mask_water])
    Velocities[:,:,i+1]=V
    if i%1000==0
        println(i)
    end
    i+=1
end


function VAC(Values)
    N=length(Values)
    X=zeros(2*N)
    X[1:N]=Values
    X=fft(X)
    X=X.*conj.(X)
    X=ifft(X)
    X=real.(X[1:N])
    return X
end


Vac=zeros(frames)
k=1
K=N_water
@time while k<=K
    global Vac,k
    Vac+=(VAC(Velocities[1,k,:]).+VAC(Velocities[2,k,:]).+VAC(Velocities[3,k,:]))./3
    if k%1000==0
        println(k)
    end
    k+=1
end
Vac./=K
D=sum(Vac)/(10*frames)
plt.plot(Vac)

plt.show()

flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Finished")

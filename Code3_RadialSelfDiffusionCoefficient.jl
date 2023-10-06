using Chemfiles
using LinearAlgebra
using PyCall

using FFTW
using StatsBase

flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Start")
GC.gc()
plt=pyimport("matplotlib.pyplot")
pygui()


batchsize=500  #Change this number if you want to change the length of each VACF calculation. Longer batches increased the accuracy for the VACF but reduced the position accuracy.


const start_water=1804
const m_O=15.999
const m_H=1.00784
const m_H20=m_O+2m_H

#calculates the center of mass positionn and center of mass velocietes for Water molecules
function center_of_mass(values)
    com=(m_O*values[:,1:3:end] .+m_H*values[:,2:3:end].+m_H*values[:,2:3:end])/m_H20
    return com
end

#Transforms two points into an vector that conforms to periodic boundary conditions
function pbc_dir_vector(A,B,box,r_box)
    dx=A.-B
    dx .-= floor.(Int64, dx .* r_box .+ 0.5) .* box
    return dx
end

#Mean position between two points that conforms to periodic boundary conditions
function pbc_mean(A,B,box,r_box)
    dx=A.-B
    dx .-= floor.(Int64, dx .* r_box .+ 0.5) .* box
    A=A-0.5*dx
    return A
end

#Calculates the velocity autocorraletion function
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


topologyname="Traj/CC_0M.gro"    #Change this to your topology file. It needs to be a gro file otherwise the reading of the topology might not work. The coiled-coil is best cenbtered at the origin without rotation/movement during the simulation. GROMACS can calculate out such movements
filename="Traj/CC_0M.trr" #Chanfe this to your trajectry file. Can be any gromacs format. The coiled-coil is best cenbtered at the origin without rotation/movement during the simulation. GROMACS can calculate out such movements

top_traj=Trajectory(topologyname)
const top_frame=read_step(top_traj,0)
const top=Topology(top_frame)


#Setuot and read in of the topology

Number_Atoms=length(top_frame)
Residues_Names=Array{String}(undef,Number_Atoms)
Atom_Names=Array{String}(undef,Number_Atoms)
@time for i in 1:Number_Atoms
    global Atom_Names, Residue_Names,Mass
    Residues_Names[i]=name(residue_for_atom(top,i-1))
    Atom_Names[i]=name(Chemfiles.Atom(top_frame,i-1))
end


const mask_water= (Atom_Names .== "OW") .|  (Atom_Names .== "HW1") .| (Atom_Names .== "HW2")
const mask_CA=("CA" .== Atom_Names).&("LEU" .== Residues_Names)

const N_water=Int64(sum(mask_water)/3)
const N_CA=sum(mask_water)


#Evaluates a trajectory frame and returns the center of mass corrected velocity of each water molecule, their distance to the center of the coiled-coil and the simulation box parameters.
#For the distance, it assumes the coiled-coil is linear and spans a vector between the start and end bead. Each bead is the mean position of C_alpha of leucine. 
function evaluate_frame(frame)
    Pos=positions(frame)
    Pos_water=Pos[:,mask_water]
    Pos_CA=Pos[:,mask_CA]
    Vel=velocities(frame)
    Vel_water=Vel[:,mask_water]
    com_p=center_of_mass(Pos_water)
    com_v=center_of_mass(Vel_water)
    box=diag(matrix(UnitCell(frame)))
    r_box=1 ./box
    beads=pbc_mean(Pos_CA[:,1:8],Pos_CA[:,9:16],box,r_box)

    Tail=beads[:,1]
    vect=pbc_dir_vector(beads[:,end],beads[:,1],box,r_box)
    vect=beads[:,end].-Tail
    com_p=pbc_dir_vector(com_p,Tail,box,r_box)

    distances(X)=norm(cross(X,vect))/norm(vect)
    dist=[mapslices(distances,com_p,dims=1)... ]

    return com_v,dist,box[3]
end



const traj=Trajectory(filename)
const frames=Int64(length(traj))


Velocities=zeros(3,N_water,frames)
Distances=zeros(N_water,frames)
Box_Height=zeros(frames)

#loop over all frames and print out the frame number very 1000 frames
#This is an old implementation and read_step can lead to memory problems.

i=0
@time while i < frames
    global i
    V,R,Box_Z=evaluate_frame(read_step(traj,i))
    Velocities[:,:,i+1]=V
    Distances[:,i+1]=R
    Box_Height[i+1]=Box_Z
    if i%1000==0
        println(i)
    end
    i+=1
end



#Calculating of the velocity autocorrelation function 
flush(stdout)
 
batches=floor(Int64,frames/batchsize)
Rmax=ceil(Int64,maximum(Distances))



Vac2=zeros(Rmax,batchsize)
R_N=zeros(Rmax)
l=0
i=0
@time while l < batches
    global l,R_N,Vac2,Distances,i
    k=1
    while k<=N_water
        global Vac2,Distances,R_N,l,i
        vac=(VAC(Velocities[1,k,l*batchsize+1:(l+1)*batchsize]).+VAC(Velocities[2,k,l*batchsize+1:(l+1)*batchsize]).+VAC(Velocities[3,k,l*batchsize+1:(l+1)*batchsize]))./3
        r_i=round.(Int64,mean(Distances[k,l*batchsize+1:(l+1)*batchsize]))
        if r_i ==0
            r_i=1
        end
        Vac2[r_i,:]+=vac
        R_N[r_i]+=1
        k+=1
        i+=1
    end
    l+=1
end

#Histogram creation
R=collect(1:Rmax)
volume=2*R*mean(Box_Height)*π
@time hist=StatsBase.fit(Histogram,[Distances...],(0.5:1:Rmax+0.5))

Vac2=Vac2./R_N
#the Diffusion coefficient is calculated from the mean VACF for a certain distance. 
D=mapslices(sum,Vac2,dims=2)./batchsize .*0.01


mask_1=R_N.>200 #Number of points needed to be considered a valid results. 200 is chosen very arbitrary.
Na=6.022E23
u=1.66E-24
M_water=18.02

#This part is for plotting the results. Theresults for bulkwater come from a seperate simulation that was evaluated with gromacs. 
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

a=ax1.scatter(R[mask_1]/10,D[mask_1].*10,label="D water at radial distance",color="g")
c=ax1.plot([0,Rmax]/10,[6.78,6.78],label="D bulk water",color="k")
d=ax2.plot(R[1:31]/10,hist.weights[1:31]./volume[1:31]/frames*18*1.66,label="water density",color="b")
e=ax2.plot([0,Rmax]/10,[33,33]*18*1.66/1000,linestyle="dashed",label="bulk water density", color="b" )
f=ax1.plot([1.4,1.4],[0,8],color="cyan",linestyle="dotted")
f=ax1.plot([2.4,2.4],[0,8],color="magenta",linestyle="dotted")
f=ax1.plot([3.0,3.0],[0,8],color="red",linestyle="dotted")
L=[a,c[1],d[1],e[1]]
LL=[a.get_label(),c[1].get_label(),d[1].get_label(),e[1].get_label()]
ax1.legend(L,LL,loc="lower right",fontsize=16)
ax1.set_xlabel("radial distance from coiled-Coil core in nm",fontsize=16)
ax1.set_ylabel("diffusion coefficient in 10^-5 cm²/s",color="g",fontsize=16)
ax2.set_ylabel("water density in kg/l",color="b",fontsize=16)
ax1.set_ylim(bottom=0,top=7.5)
ax2.set_xlim(left=0, right=4.6)
ax1.set_xlim(left=0, right=4.6)
ax1.tick_params(axis="both", which="major",labelsize=16)
ax2.tick_params(axis="both", which="major",labelsize=16)
plt.tight_layout()
plt.savefig("RadialDiffusion.pdf")
#plt.show()
plt.close()


flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Finished")

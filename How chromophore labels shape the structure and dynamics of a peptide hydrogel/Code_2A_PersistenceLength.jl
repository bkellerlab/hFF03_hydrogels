###Version 1.03 from 3.9.2021
# This programm calculates the persistence length from a trajectory.
# In its current form the trajectory only contains the positions of Leucine Cα atoms.To obtain these see the srun.sh file
# To install a Julia package, start julia REPL pres ]  and write add PackageName . Julia needs Internet connection to download packages.
# To properly use PyCall a  Pythin 3.x installation is needed.
using ArgParse
using Chemfiles
using PyCall
using LinearAlgebra
using DelimitedFiles
using Printf
using StatsBase
using LsqFit
using Combinatorics
using Statistics

flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Start")
plt=pyimport("matplotlib.pyplot")
pygui()

# Functions for fitting. There is no difference between Krohn und Krohn2
exp_fit(L,P)=exp.(-L./P)
Krohn(L,P)=2 .*P.*L.*(1 .-P./L.*(1 .-exp.(-L./P)))
Krohn2(L,P)=2 .* P .* L .- 2 .*P.^2 .*(1 .- exp.(-L./P))


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--begin" , "-b"
            help = "Frame where to start the claculation."
            arg_type = Int
            default = 1
        "--end" , "-e"
            help = "Frame where to end the claculation."
            arg_type = Int
        "--out", "-o"
            help = "Name of the output file"
            arg_type = String
            default = "out.txt"
        "--file" , "-f"
            help = "Filename of trajectory. All trajectory have to follow the format Chain_TypeX_xxx.gro for correct recognition of the Chain Type. Supported are TypeA, TypeB , TypeC, TypeD(not jet implemented) and Coil. All files sould only contain C_alpha of Leucine."
            arg_type=String
            default="Traj/Chain_TypeC_1.gro"
        "--folder"
            help="Folder with multiple Trajectories for Input"
            arg_type=String
            #default="Test"
        "--fout"
            help="Outputfolder for the trajectories and plots"
            arg_type=String
            default="Output"
        "--traj"
            help="Set flag if you want a gro file of the mean chain trajectory"
            action=:store_true
            #default=true
        "--plot"
            help = "Flag wether Plots of for correlation function is made"
            action = :store_true
            #default=true
        "--nocentering"
            help="Set flag to center the Chain. Not all calculations work correctly if the chain is not centered"
            action=:store_true
            #default=true
        "--append"
            help="Set flag if the results should only be appended in the output file and not overwritten. If the outputfile cannot be found a new one will be generated"
            action=:store_false
            #default=false
    end
    return parse_args(s)
end

#Setup the calculation and read the ARGS from the command line. If you want to calculate interactivly change the default values above and make the lines active.
parsed_args = parse_commandline()
out=parsed_args["out"]
file = parsed_args["file"]
folder = parsed_args["folder"]
fout = parsed_args["fout"]


mkpath(fout)
if parsed_args["traj"]
    mkpath("$fout/Trajectories")
end
if parsed_args["plot"]
    mkpath("$fout/Plots")
end

#Write header for outputfile
if parsed_args["append"] || !isfile("$fout/$out")
    io=open("$fout/$out","w")
    write(io,"### Resutls of Persistence Length calculation for different types. Results in nm \n")
    write(io,"Filename & Type & Contour_Length & Persistence_Length_single & Persistence_Length_full & Persistence_Length_single_Fit & Modus & how Persistence_Length_full_Fit & Modus & Persistence_Length_Krohn & Frames & Notes \\\\ \n")
    close(io)
end


#Decapricated. Calculates the distance between points in a PBC system
function pbc_dist(A,B,box,r_box)
    dx=A.-B
    dx .-= floor.(Int64, dx .* r_box .+ 0.5) .* box
    absd=sqrt.(mapslices(sum,(dx.^2),dims=1))
    return dx./absd,mean(absd)
end

# Calcuates the direction vector between 2 points under consideration of periodic boundary conditions.
# A and B can be arrays instead of vectors.
function pbc_vector(A,B,box,r_box)
    dx=A.-B
    dx .-= floor.(Int64, dx .* r_box .+ 0.5) .* box
    return dx
end

#Calculate the squared end to end distance between the start and end of a list.
function squared_end_to_end_distance(Positions)
    return sqrt.( sum((Positions[:,1]-Positions[:,end]).^2) )^2
end

#Setup to get the correct indices for the different types. The trajectory needs to be ordered in such a way that first all positions on the Helix A and then all positions on Helix B are written
function index_types(type_,traj,file)
    N=size(positions(read_step(traj,0)))[2]
    try
        N2=Int64(N/2)
        if type_ == "TypeA"
            index_A=collect(3:N2)
            index_B=collect(N2+1:N-2)
        elseif type_ == "TypeB" ## This is a simplistic stop gap measure. Best to write down the chains above each other to see what it exactly does
            index_A=collect(4:N2)
            index_B=collect(N2+1:N-3)
        elseif type_ == "Coil" ## For pure coils the end begin to open slightly, so the first and alst Heptadt of every helix is cut from calculation. It is less of a problem if the whole system is more flexible like in the Types A-D
            index_A=collect(3:N2-2)
            index_B=collect(N2+3:N-2)
        elseif type_== "TypeD"
            index_A=collect(3:N2)
            index_B=collect(N2+1:N-4)
            splice!(index_A,17)
            splice!(index_A,17)
        else ##TypeD not jet implemented
            index_A=collect(1:N2)
            index_B=collect(N2+1:N)
        end
        return index_A,index_B,size(index_A)[1]
    catch
        println("There is an uneven number of positions in file: $file")
    end
end

#Reads a single frame, creates the mean of both helices and centeres the helix to the first atom.
function read_frame(x,index_A,index_B,traj)
    frame=read_step(traj,x)
    positions_A=positions(frame)[:,index_A]
    positions_B=positions(frame)[:,index_B]
    box=diag(matrix(UnitCell(frame)))
    r_box= 1 ./box
    Positions=positions_A.-0.5 .*pbc_vector(positions_A,positions_B,box,r_box)
    if parsed_args["nocentering"] != true
        Positions=centering(box,Positions)
    end
    return box,r_box,Positions
end

# Output of the mean chain trajectory as a gro file.
function write_trajectory(name,Mean,Box,frames,traj_begin_frame)
    if parsed_args["nocentering"] != true
        name="Centered_$name"
    end
    io=open("$fout/Trajectories/$name.gro","w")
    N=size(Mean)[2]
    for x in 1:frames
        t=x+traj_begin_frame-1
        write(io,"Chain t= $t\n$N\n")
        for i in collect(1:N)
            write(io,@sprintf "%5i%5s%5s%5i%8.3f% 8.3f% 8.3f\n" 1 "CHAIN" "X" i Mean[1,i,x]/10 Mean[2,i,x]/10 Mean[3,i,x]/10)
        end
        write(io,@sprintf "%8.5f% 8.5f% 8.5f\n" Box[1,x]/10 Box[2,x]/10 Box[3,x]/10 )
    end
    close(io)
end

# Write the results of each file into the output file.
function write_results(name,type,contour_length,persistence_length,mode_,begin_frame,end_frame)
    io=open("$fout/$out","a")
    write(io,@sprintf " %.20s & %.10s & %2.2f & %2.2f & %2.2f & %2.2f & %.10s & %2.2f & %.10s & %2.2f & %.5i:%.5i &  \\\\ \n" name type contour_length persistence_length[1] persistence_length[2] persistence_length[3] mode_[1] persistence_length[4] mode_[2] persistence_length[5] begin_frame end_frame)
    close(io)
end

# Crerate a plot for a file
function make_plot(correlation1,correlation2,distances,persistend_length1,persistence_length2,Name)
    plt.plot(distances,correlation1,label="Correlation",color="k")
    plt.plot(distances,exp_fit(distances,persistend_length1),label="Fit",color="k",linestyle=":")
    plt.plot(distances,correlation2,label="Correlation averaged over chain",color="b")
    plt.plot(distances,exp_fit(distances,persistence_length2),label="Fit",color="b",linestyle=":")
    plt.title("Peristence Length for the file: \n $Name")
    plt.legend()
    plt.xlabel("Length in Å")
    plt.ylabel("Chain correlation")
    plt.savefig("$fout/Plots/$Name.pdf")
    plt.close()
end

# Centeres the Chain. All positions are now relativ to the first atom in the list. (No more PBC)
# Start can also be a different Positions like [0.0, 0.0, 0.0 ]
function centering(box,Positions)
    r_box=1 ./box
    Start=Positions[:,1]
    Positions=pbc_vector(Positions[:,2:end],Positions[:,1:end-1],box,r_box)
    Positions=hcat(Start,Positions)
    Positions=cumsum(Positions,dims=2)
    return Positions
end

#Solve the Krohn equation for Persistence Length. Curve fit is simply the easiest to implement.
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

# Setup parameters that stay constant for each trajectory file
function setup_trajectory(file)
    flush(stdout) ##Only necessary if you run the program in Atom IDE
    traj_begin_frame=parsed_args["begin"]
    traj_end_frame=parsed_args["end"]
    try
        traj=Trajectory(file)
        traj_length=Int64(length(traj))


        if isnothing(traj_end_frame)
            traj_end_frame=traj_length
        elseif traj_end_frame>traj_length
            traj_end_frame=parsed_args["end"]
            println("The choosen endframe of $traj_end_frame is not inside the trajectory for the file: $file \n Last frame chosen as endframe")
            traj_end_frame=traj_length
        else
            traj_end_frame=parsed_args["end"]
        end

        if traj_begin_frame>traj_end_frame
            println("The choosen beginning frames is after the end frame for the file: $file")
            throw()
        elseif traj_begin_frame>traj_length
            println("The choosen beginning frame is longer than the trajectory for the file: $file")
            throw()
        else
            frames=traj_end_frame-(traj_begin_frame-1)
        end
        return traj,traj_begin_frame,traj_end_frame,frames
    catch
    println("Something went wrong while trying to load the file: $file")
    end
end

## Main function
function evaluate_traj(file)
    filename=split(file,"/")[end]
    name=split(filename,".")[1]
    type=split(filename,"_")[2]

    traj,traj_begin_frame,traj_end_frame,frames=setup_trajectory(file)
    index_A,index_B,number_chain_positions=index_types(type,traj,file)

    number_vectors=number_chain_positions-1
    number_dot_vectors=number_chain_positions-2
    table_full=collect(combinations(collect(1:number_vectors),2))       ## All Combinations for averaging over the chain positions and time
    table_single=table_full[1:number_dot_vectors]                       ## For only averaging over time.
    ## Single here and all further instances means the calculation is only done from the start to the end of the chain and not all possible chain element combiunations.
    ## You should take results from the full calculation, but single is helpfull to find error or strange behaviour.
    Box=zeros(3)
    Frame_positions=zeros((3,number_chain_positions))
    Distances=zeros(number_vectors)
    Cos_α_single=zeros(number_dot_vectors)
    Cos_α_full=zeros(number_dot_vectors)
    R2=0

    for i in traj_begin_frame:traj_end_frame
        box,r_box,Positions=read_frame(i-1,index_A,index_B,traj)
        Box=hcat(Box,box)
        Frame_positions=cat(Frame_positions,Positions,dims=3)
        ##calculate values for Persistence Length
        R2=vcat(R2,squared_end_to_end_distance(Positions))
        vectors=Positions[:,2:end]-Positions[:,1:end-1]
        distances=sqrt.(mapslices(sum,(vectors.^2),dims=1))
        vectors ./= distances                                       ## Now they are unit vectors
        Distances=hcat(Distances,[distances...])
        cos_α_single=calculate_dot_products(vectors,table_single,number_dot_vectors)
        cos_α_full=calculate_dot_products(vectors,table_full,number_dot_vectors)
        Cos_α_single=hcat(Cos_α_single,cos_α_single)
        Cos_α_full=hcat(Cos_α_full,cos_α_full)
    end
    Box=Box[:,2:end]
    Frame_positions=Frame_positions[:,:,2:end]
    Distances=Distances[:,2:end]
    Cos_α_single=Cos_α_single[:,2:end]
    Cos_α_full=Cos_α_full[:,2:end]
    R2=R2[2:end]

    Mean_Distance=mean(Distances)
    Mean_Distances=mean(Distances,dims=2)
    Cum_Distances=[cumsum(Mean_Distances,dims=1)...]
    Contour_Length=Cum_Distances[end]
    Correlation_single=[mean(Cos_α_single,dims=2)...]
    Correlation_full=[mean(Cos_α_full,dims=2)...]

    Persistence_Length_single_direct=sum(Correlation_single)*Mean_Distance
    Persistence_Length_full_direct=sum(Correlation_full)*Mean_Distance

    Persistence_Length_single_fit,mode_single=calculate_persistence_length_fit(Correlation_single,Cum_Distances[1:end-1],Persistence_Length_single_direct,number_dot_vectors)
    Persistence_Length_full_fit,mode_full=calculate_persistence_length_fit(Correlation_full,Cum_Distances[1:end-1],Persistence_Length_full_direct,number_dot_vectors)

    Persistence_Length_Krohn=calculate_persistence_length_Krohn(R2,Contour_Length,Persistence_Length_single_direct)
    Persistence_Length=[Persistence_Length_single_direct,Persistence_Length_full_direct,Persistence_Length_single_fit...,Persistence_Length_full_fit...,Persistence_Length_Krohn...]

    if parsed_args["traj"]
        write_trajectory(name,Frame_positions,Box,frames,traj_begin_frame)
    end
    if parsed_args["plot"]
        make_plot(Correlation_single,Correlation_full,Cum_Distances[1:end-1],Persistence_Length[3],Persistence_Length[4],name)
    end
    try
    write_results(name,type,Contour_Length./10,Persistence_Length./10,[mode_single,mode_full],traj_begin_frame,traj_end_frame)
    catch
        println("Error during writing the results for file: $file")
    end
    println("Finished with file: $file")
    return Correlation_full, Mean_Distances, Frame_positions,Contour_Length,R2,Persistence_Length_Krohn,Persistence_Length_single_direct,Distances
end

cor,dist,pos,con,R2,cohn,LPS,distances=evaluate_traj("Traj/Chain_TypeC_1.gro")

flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Calculations Finished")

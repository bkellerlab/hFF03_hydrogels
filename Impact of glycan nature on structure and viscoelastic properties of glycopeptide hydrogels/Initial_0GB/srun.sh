#!/bin/sh

#SBATCH --mail-user=f.heinz@fu-berlin.de
#SBATCH --job-name=OMA-1
#SBATCH --mail-type=all
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --qos=standard
#SBATCH --mem=6g
#SBATCH --gres=gpu:2
#SBATCH --partition=gpu

module add GROMACS/2020-fosscuda-2019b
cd /scratch/heif89/0GB/1

gmx insert-molecules -ci protein.gro -o box.gro -box 10 10 10 -nmol 1					
gmx solvate -cp box.gro -o solvate.gro -cs tip4p.gro -p topol.top				
gmx grompp -f ions.mdp -c solvate.gro -p topol.top -o ions.tpr -maxwarn 1				
echo SOL | gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.1
				
mkdir Setup												
mv * Setup
mkdir NVT
mkdir NPT
mkdir RUN

cp Setup/*.top NVT/
cp Setup/*.top NPT/
cp Setup/*.top RUN/

cp Setup/*.itp NVT/
cp Setup/*.itp NPT/
cp Setup/*.itp RUN/

cp Setup/nvtrun.mdp NVT/
cp Setup/nptrun.mdp NPT/
cp Setup/mdrun.mdp RUN/

cp -r Setup/amber99sb-ildn.ff NVT/
cp -r Setup/amber99sb-ildn.ff NPT/
cp -r Setup/amber99sb-ildn.ff RUN/

cd Setup

gmx grompp -f minim.mdp -c ions.gro -p topol.top -o em.tpr					#Run Energy minimasation
gmx mdrun -v -deffnm em -nb gpu

cp em.gro ../NVT/.
cd ../NVT

gmx grompp -f nvtrun.mdp -c em.gro -p topol.top -o nvt.tpr					#Run NVT relaxation
gmx mdrun -v -deffnm nvt -nb gpu

cp nvt.gro ../NPT/.
cd ../NPT

gmx grompp -f nptrun.mdp -c nvt.gro -p topol.top -o npt.tpr					#Run NPT relaxation
gmx mdrun -v -deffnm npt -nb gpu

cp npt.gro ../RUN/.
cd ../RUN
			
gmx grompp -f mdrun.mdp -c npt.gro -p topol.top -o run.tpr					#Run System
gmx mdrun -v -deffnm run -nb gpu

18 | gmx trjconv -f run.trr -o 0MA.gro -pbc mol -s run.tpr


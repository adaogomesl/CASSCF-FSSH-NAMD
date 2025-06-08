# 1. Single-point HF - Generate initial orbitals
## References
GATEWAY: https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/gateway.html
SEWARD: https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/seward.html    
SCF: https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/scf.html


## Files necessary
```.inp```,```.xyz``` (from dft opt) and ```.sh``` 

#### .inp example
```
>> EXPORT MOLCAS_MAXITER=10000
&GATEWAY
 coord=dbhme-T-sp-HF.xyz
 basis=ano-s-vdzp
 Group=c1
 RICD

&SEWARD
doanalytic

&SCF
```

#### .xyz example
```
18
dbh
 C                  1.49605900   -0.55258700    0.78740200
 C                  1.49540400   -0.56520300   -0.77950200
 C                  0.05621200   -0.05471400   -1.11121600
 C                  0.05686800   -0.03605400    1.11192000
 C                 -0.79550600   -0.71081500    0.00565100
 N                  0.09210500    1.40950500   -0.64606700
 N                  0.09271600    1.41816200    0.62448100
 H                 -0.63772200   -1.79563800    0.01546000
 H                 -0.26299900   -0.09421300    2.14897300
 H                 -0.26448300   -0.12705800   -2.14708400
 H                  2.25600200    0.09028200   -1.21102000
 H                  2.25642900    0.11074400    1.20732500
 H                  1.64624000   -1.55223700    1.20659400
 H                  1.64467900   -1.57196100   -1.18189700
 C                 -2.29233100   -0.39597800    0.00323900
```

#### .sh example
```
#!/bin/sh
## script for OpenMalCas
## $INPUT and $WORKDIR do not belong to OpenMolCas
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=23:00:00
#SBATCH --job-name=dbh-me-T-sp-HF
#SBATCH --partition=short
#SBATCH --mem=6Gb
#SBATCH --output=%j.o.slurm
#SBATCH --error=%j.e.slurm

export INPUT=dbhme-T-sp-HF
export WORKDIR=$PWD

export MOLCAS_NPROCS=$SLURM_NTASKS
export MOLCAS_MEM=3200
export MOLCAS_PRINT=3
export MOLCAS_PROJECT=$INPUT
export MOLCAS_WORKDIR=$WORKDIR
export OMP_NUM_THREADS=1
export MOLCAS=/work/lopez/Molcas

module load python/3.7.1

cd $WORKDIR
$MOLCAS/bin/pymolcas -f $INPUT.inp -b 1
rm -r $MOLCAS_PROJECT
```
## Analysis
Check .log and .status should have ```Happy landing!```

--- 
# 2. Single-point CASSCF - Generate CAS orbitals

## References
RASSCF: [https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/gateway.html](https://molcas.gitlab.io/OpenMolcas/sphinx/users.guide/programs/rassi.html)

## Files necessary
```.inp```,```.xyz```,```.sh``` and ```.StrOrb``` (ScfOrb  from previous step)

#### .inp example
```
&GATEWAY
 coord=$Project.xyz
 basis=ano-s-vdzp
 group=C1
 RICD

&SEWARD
 doanalytic

&RASSCF
 FileOrb=$Project.StrOrb
 Spin=1
 Charge=0
 NActel=8
 Ras2=9
 CIRoot=8 8 1
```

## Analysis
1. Check .log and .status should have ```Happy landing!```
2. Open .rasscf.molden on ChemCraft to vizualize orbitals in the active space (Energy = 0)
3. If active space is not as planned, a new cas-sp will be submitted:
   
   a. Use the orbitals generated from previous cas-sp (RasOrb) as initial orbitals (StrOrb)
   
   b. Use the alter command to get orbitals in the correct order
   
   c. Repeat the steps a and b as needed


#### .inp example with alter
```
&GATEWAY
 coord=$Project.xyz
 basis=ano-s-vdzp
 group=C1
 RICD

&SEWARD
 doanalytic

&RASSCF
 FileOrb=$Project.StrOrb
 Spin=1
 alter
 1
 1 58 33
 Charge=0
 NActel=8
 Ras2=9
 CIRoot=8 8 1
```


--- 
# 3. Optimization and Frequency 
## References
SLAPAF: https://www.molcas.org/documentation/manual/node115.html

# 4. CASSCF VEE
## Analysis
1. Check .log and .status should have ```Happy landing!```
   
2. Check if the orbitals in the active space are still the same

3. Extract the following data by going to the end of ```.log``` by pressing ```SHIFT G``` and searching for the following keywords:
   
   a. Energies - convert to eV (multiply by 27.2114) and subtract the S0 Energy to get the VEE
   ```
   ?RASSI
   ```
   b. Oscillator strength
   ```
   ?Dipole transition strengths
   ```
   c. Transitions
   ```
   ?Wave function printout
   ```






# 5. CASPT2 VEE
## Analysis
1. Check .log and .status should have ```Happy landing!```
   
2. Check if the orbitals in the active space are still the same

3. Extract the following data by going to the end of ```.log``` by pressing ```SHIFT G``` and searching for the following keywords:
   
   a. Energies - convert to eV (multiply by 27.2114) and subtract the S0 Energy to get the VEE
   ```
   ?XMS-CASPT2
   ```
   

# 6. Absorption Spectrum
## Files necessary
```.inp```,```.xyz```,```.sh```, ```.StrOrb```, ```control```,Gen-FSSH.py```. Input is the same as XMS-CASPT2

Scripts
```
/work/lopez/share_from_Leticia/model-scripts/spectrum
```

## Submission

1. Request resources
```
srun -N 1 --exclusive --partition=short --time=23:59:59 --pty /bin/bash
```
2. Load python
```
module load python/3.7.1
```

4. Generate input file and wigner samples - It will generate the wigner.xyz
```
python Gen-FSSH.py -x control
```
5. To run dynamics, add in the first line of the runall.sh:  
```
 #!/bin/sh
```
6. Submit runall

## Analysis
1. Use ```collector.sh``` to extract energies
   ```
   sbatch collector.sh
   ```
2. Check generated ```data.txt```, and replace blank oscillator to 0
3. Plot spectrum
   ```
   python3 plot-abs.py
   ````


# 7. MEP
## References
SLAPAF: https://www.molcas.org/documentation/manual/node115.html

## Files necessary
```.inp```,```.xyz```,```.sh```, ```.StrOrb```

#### .inp example
```
&GATEWAY
 coord=$Project.xyz
 basis=ano-s-vdzp
 group=C1
 RICD

>> do while

&SEWARD
 doanalytic

&RASSCF
 FileOrb=$Project.StrOrb
 Spin=1
 Symmetry=1
 Charge=0
 NActel=8
 Ras2=9
 CIRoot=4 4 1
 Rlxroot=2

&SLAPAF
 MEP-search
 Iterations=200
 THRShld= 0.0 5.0D-4

>> end do
```
If small amount of geometries are generated, change MEPStep




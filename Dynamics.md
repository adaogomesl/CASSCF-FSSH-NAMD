## Submission
### Model scripts
All files mentioned on this documentation are uploaded on https://github.com/adaogomesl/Leticia-LopezLab/tree/main/Molcas/Dynamics

### Necessary files
*From Optimization:*
1. .freq.molden
2. .xyz
3. StrOrb

*Script/File*
1. Gen-FSH.py
2. Control (number of molecules = nodes x cores)

*Input file*
```
&GATEWAY
 coord=$MOLCAS_PROJECT.xyz
 basis=ano-s-vdzp 
 Group=c1
 RICD

>> FOREACH ITER in (1 .. 2000)

&SEWARD
doanalytic

&RASSCF
 Fileorb=$MOLCAS_PROJECT.StrOrb
 Spin=1
 Nactel=10 0 0 
 Charge=0
 Ras1=0
 Ras2=9
 Ras3=0
 ITERATIONS=200,100
 CIRoot=8 8 1
 MDRLXR=2

&Surfacehop
 tully
 decoherence = 0.1

&ALASKA

&Dynamix
 velver
 dt = 20
 velo = 3
 ther= 2
 temp = 300

>> End do
```

### Running dynamics
1. Request resources
```
srun -N 1 --exclusive --partition=short --time=23:59:59 --pty /bin/bash
```
2. Load python
```
module load python/3.7.1
```
3. Edit the name on the control file
```
      input      dbh-unsub.freq.molden
      temp       300
      method     wigner
      partition  lopez
      time       30-00:00:00
      memory     2000
      nodes      20
      cores      25
      jobs       25
      index      1
```
5. Generate input file and wigner samples - It will generate the wigner.xyz
```
python Gen-FSSH.py -x control
```
5. To run dynamics, add in the first line of the runall.sh:  
```
 #!/bin/sh
```
6. Submit runall

## Analysis
1. Request resources
```
srun -N 1 --exclusive --partition=short --time=23:59:59 --pty /bin/bash
```
2. Load python
```
module load python/3.7.1
```
### Diagnostic - Check final timestep
```
python3 HOP-FSSH.py -x 01-diagnostic
```
1. Create txt file with timesteps printed on the screen
2. Make a histogram to get the distribution of the timestep
   
   <img width="391" alt="Screenshot 2023-09-04 at 7 34 36 AM" src="https://github.com/adaogomesl/Leticia-LopezLab/assets/100699955/58159589-e995-42ec-8052-15941f4d4d8e">

   
### Determine Final State
```
python3 HOP-FSSH.py -x 02-finalstates
```
xyz files will be generated with the final structure for each state and also for the CI

1. Classify final products using classify_product.py
```
python3 classify_product.py Prod.S0.xyz
```
2. Use results on ```.csv``` to classify products


### State Population Analysis
1. Create file ```index-S0``` that contains the path to the folders with trajectories that landed on S0
   
   a. Create file ```range-S0.txt``` that contains the S0 range. Example file:
   ```
   1 3-9 11-16 18-26 28-35 37-45 47-82 84-88 90-106 108-117 120 121 123 124 126-134 136-138 140-152 154-179 181-184 186-193 195-214 216-231 233-235 237-245 247 248 250-266 268-328 330-338 340-347 349-355 357-365 367-369 371-376 378-409 411-453 455-465 467-472 474-476 480-493 495-500
   ```
   
   b. Convert range to index using ```range2index.py```
   ```
   python3 range2index.py
   ```

   c. Copy the indexes to a file named ```states.txt```, it will be used on later steps
   
   e. Copy the indexes list to ```generatepathtoS0.py``` and update the path

   f. Generate the paths:
   ```
   python3 generatepathtoS0.py
   ```

   g. Copy paths to a file named ```index-S0```

2. Extract population data, ```.dat``` file generated
   ```
   python3 HOP-FSSH.py -x 03-pop
   ```

3. Edit the ```plot-state-pop-2d.py``` to the correct number of states. Example 4 states
   ```
   #Format .dat file timestep, Ekin, Epot0-3, Etot0-3,Pop0-3
   #Population info will start at index  9
   p1=data[:,9]
   p2=data[:,10]
   p3=data[:,11]
   p4=data[:,12]
   ```

3. Create a plot showing the population with ```plot-state-pop-2d.py```
```
export python="$PATH:/work/lopez/Python-3.7.4/lib"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/lopez/Python-3.7.4/lib
/work/lopez/Python-3.7.4/bin/python3 plot-state-pop-2d.py average-dbh-unsub.dat
```
![dbh-unsub-pop](https://github.com/adaogomesl/Leticia-LopezLab/assets/100699955/ac8a2797-5bb4-4e5f-aa79-b7e274ef02b2)

The timestep of the half-life of S0 and S1 will be printed, to compute the time constant, divide by 2.


### Spaghetti Plots
1. Extract geometrical parameters, ```.json```file created 
```
python3 HOP-FSSH.py -x 04-geometricalpar
```
2. Generate plots
   ```
   python3 chart_v2.py
   ```
   
<img width="894" alt="Screenshot 2023-09-04 at 10 16 40 AM" src="https://github.com/adaogomesl/Leticia-LopezLab/assets/100699955/818f25ec-afae-476d-9349-b336a53115f0">
  

3. Generate 3D-plot
   a. Edit the ```chart_v2_3D-byproduct.py```, to have the correct index for each product. To do that, I copied the index from the Excel sheet with classified products and pasted it on a browse to get rid of the formation. Using TextMate I replaced the spaced by ,.
   
   b. Also, change the energy by the energy obtained on the optimization ```RASSCF state energy```

   c. Generate plot
   ```
   python3 chart_v2_3D-byproduct.py
   ``` 









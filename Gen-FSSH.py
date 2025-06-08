## Molcas TSH Generator now support NewtonX/Bagel
## Jingbai Li Dec 16 2019
## minor fix, load module openblas/0.3.6, copy velocity.xyz, Dec 19 2019 Jingbai Li
## add new option for lopez_new partition, Dec 30 2019 Jingbai Li
## add restart function, Jan 14 2020 Jingbai Li
## adapted for new Discovery partitions Feb 2 2020 Jingbai Li
## add initial excited state selection Feb 6 2020 Jingbai Li
## modification on tmp folder detection Feb 12 2020 Jingbai Li
## add 5s delay in submission loop Feb 17 2020 Jingbai Li
## add function to sample initial conditions using molden,bagel,newton,xyz format May 2 2020 Jingbai Li
## add option to use NewtonX/Bagel May 2 2020 Jingbai Li
## small fix on individual submission script Aug 7 2020 Jingbai Li
## support PyRAI2MD Sep 8 2020 Jingbai Li
## change BAGEL defaut mpi to mvapich2-2.3.4 Sep 24 Jingbai Li
## add convenient function to print the running set name in stdout Oct 25 Jingbai Li
## add optional therm.inp for Bagel  Nov 18 2020 Jingbai Li
## improve the PyRAI2MD-Molcas/NNs generation
## support PyRAI2MD/Molcas hybrid trajectory files
## support Molcas customized basis set Feb 22 2021 Jingbai Li
## add ld_library enviroment variable for python lib Mar 03 2021 Jingbai
## support Molcas/Tinker initcondition Oct 10 2021 Jingbai Li
## update the pyrai2md functions Oct 14 2021 Jingbai Li

import sys,os,shutil
import numpy as np
from optparse import OptionParser

def main():
    ##  This is the main function
    ##  It read all options from command line or control file and 
    ##  pass them to other module for sampling and generating structures
    ##  The parser code has clear description for all used variables
    ##  The last part writes sbatch file for calculations running on each nodes

    usage="""

    FSSH Generator

    Usage:
      python3 Gen-FSSH.py -i input.file_name.molden
      python3 Gen-FSSH.py -h for help

    Or prepare a text file that contains control parameters

      input file_name.freq.molden
      seed       -1
      temp       298.15  
      method     boltzmann
      partition  long
      time       4-23:59:59
      memory     3000
      nodes      1
      cores      28
      jobs       28
      index      1
      restart    0
      initex     0
      prog       molcas
      pyqd        /work/lopez/PyQD_1.0
      molcas      /work/lopez/Molcas

      Options below are only avaiable via control file

      newton        /work/lopez/NX-2.2-B08
      bagel         /work/lopez/Bagel-mvapich
      lib_blas      /work/lopez/BLAS
      lib_lapack    /work/lopez/BLAS
      lib_scalapack /work/lopez/BLAS
      lib_boost     /work/lopez/Boost
      mkl           /work/lopez/intel/mkl
      mpi           /work/lopez/mvapich2-2.3.4
      
    For more information, please see Gen-FSSH-readme.txt 
    
    """
    description=''
    #sampling control
    parser = OptionParser(usage=usage, description=description)
    parser.add_option('-i', dest='input',   type=str,   nargs=1, help='Frequency file name.')
    parser.add_option('-s', dest='iseed',   type=int,   nargs=1, help='Random number seed (0 - +inf). Default is random.',default=-1)
    parser.add_option('-t', dest='temp',    type=float, nargs=1, help='Temperature in K. Default is 298.15K.',default=298.15)
    parser.add_option('-d', dest='dist',    type=str,   nargs=1, help='Sampling method: wigner or boltzmann (low-cases). Default is boltzmann.',default='boltzmann')
    parser.add_option('-v', dest='prog',    type=str,   nargs=1, help='Choose program for FSSH. Default is molcas.',default='molcas')
    
    #slurm control
    parser.add_option('-p', dest='slpt',    type=str,   nargs=1, help='Slurm partionInput. Default is fullnode.',default='fullnode')
    parser.add_option('-l', dest='sltm',    type=str,   nargs=1, help='Slurm time limit. Default is 2-23:59:59.',default='2-23:59:59') 
    parser.add_option('-m', dest='slmm',    type=int,   nargs=1, help='Slurm memory limit per node in MB. Default is 3000.',default=3000)
    parser.add_option('-n', dest='slnd',    type=int,   nargs=1, help='Slurm nodes number. Default is 1.',default=1)
    parser.add_option('-c', dest='slcr',    type=int,   nargs=1, help='Slurm cores per node. Default is 28.',default=28)
    parser.add_option('-j', dest='sljb',    type=int,   nargs=1, help='Slurm jobs per node. Default is 28.',default=28)
    parser.add_option('-e', dest='slin',    type=int,   nargs=1, help='Slurm jobs initial index. Default is 1.',default=1)
    #additional input
    parser.add_option('-r', dest='restart', type=int,   nargs=1, help='Number of restarted calculation. Default is 0.',default=0)
    parser.add_option('-y', dest='initex',  type=int,   nargs=1, help='Select initial excited state. Default is 0.',default=0)
    parser.add_option('-f', dest='tomlcs',  type=str,   nargs=1, help='Path to Molcas directory. Defualt is /work/lopez/Molcas.',default='/work/lopez/Molcas')
    parser.add_option('-x', dest='extnl',   type=str,   nargs=1, help='A text file that contains control parameters')

    (options, args) = parser.parse_args()

    input=options.input
    iseed=options.iseed
    temp=options.temp
    dist=options.dist
    prog=options.prog

    slpt=options.slpt
    sltm=options.sltm
    slmm=options.slmm
    slnd=options.slnd
    slcr=options.slcr
    sljb=options.sljb
    slin=options.slin

    restart=options.restart
    initex=options.restart
    tomlcs=options.tomlcs
    extnl=options.extnl

    ## defaults for newton/begal'
    pyqd='/work/lopez/PyQD_1.0/'
    tontx='/work/lopez/NX-2.2-B08'
    tobgl='/work/lopez/Bagel-mvapich'
    lbbls='/work/lopez/BLAS'
    lblpk='/work/lopez/BLAS'
    lbslp='/work/lopez/BLAS'
    lbbst='/work/lopez/Boost'
    tomkl='/work/lopez/intel/mkl'
    tompi='/work/lopez/mvapich2-2.3.4'

    ## defaults for PyRAI2MD
    topyd='/work/lopez/PyRAI2MD/'
    topython='/work/lopez/Python-3.7.4'

    if   extnl != None:   #read controls if provided
        with open(extnl,'r') as ext:
            for line in ext:
                if   'input' in line:
                    input=line.split()[1]                
                elif 'seed'  in line:
                    iseed=int(line.split()[1])
                elif 'temp'  in line:
       	       	    temp=float(line.split()[1])
                elif 'method' in line:
       	       	    dist=line.split()[1].lower()
                elif 'prog' in line:
       	       	    prog=line.split()[1].lower()
                elif 'partition' in line:
       	       	    slpt=line.split()[1]
                    if slpt=='lopez_new':
                        slpt='lopez\n#SBATCH --constraint="ib"' # This is temporary solution to lopez_new partition
                elif 'time' in line:
       	       	    sltm=line.split()[1]
       	       	elif 'memory'  in line:
       	       	    slmm=int(line.split()[1])
       	       	elif 'nodes' in line:
       	       	    slnd=int(line.split()[1])
       	       	elif 'cores' in line:
       	       	    slcr=int(line.split()[1])
       	       	elif 'jobs' in line:
       	       	    sljb=int(line.split()[1])
                elif 'index' in line:
       	       	    slin=int(line.split()[1])
                elif 'restart' in line:
                    restart=int(line.split()[1])
                elif 'initex' in line:
                    initex=int(line.split()[1])
                elif 'pyqd'   in line:
                    pyqd=line.split()[1]
                elif 'molcas' in line:
       	       	    tomlcs=line.split()[1]
                elif 'newton' in line:
                    tontx=line.split()[1]
       	       	elif 'bagel' in line:
                    tobgl=line.split()[1]
       	       	elif 'lib_blas' in line:
                    lbbls=line.split()[1]
       	       	elif 'lib_lapack' in line:
                    lblpk=line.split()[1]
       	       	elif 'lib_scalapack' in line:
                    lbslp=line.split()[1]
       	       	elif 'lib_boost' in line:
                    lbbst=line.split()[1]
       	       	elif 'mkl' in line:
                    tomkl=line.split()[1]
       	       	elif 'mpi' in line:
                    tompi=line.split()[1]
                elif 'pyrai2md' in line:
                    topyd=line.split()[1]
                elif 'python'  in line:
                    topython=line.split()[1]

    if   input != None and os.path.exists(input) == True:
        print ('\n>>> %s' % (input))
    else:
        print ('\n!!! File %s not found !!!' % (input))
        print (usage)
        print ('!!! File %s not found !!!' % (input))
        exit()

    format=input.split('.')[-1]
    input=input.split('.')[0]

    if   prog == 'molcas':
        if os.path.exists('%s.inp' % (input)) == False:
            print ('\n!!! Molcas template input %s.inp not found !!!' % (input))
            print (usage)
       	    print ('!!! Molcas template input %s.inp not found !!!\n' % (input))
            exit()
        if os.path.exists('%s.StrOrb' % (input)) == False and os.path.exists('%s.JobIph' % (input)) == False:
       	    print ('\n!!! Molcas orbital file %s.StrOrb or JobIph not found !!!' % (input))
            print (usage)
       	    print ('!!! Molcas orbital file %s.StrOrb or JobIph not found !!!\n' % (input))
       	    exit()
    elif prog == 'nxbagel':
        if os.path.exists('control.dyn') == False:
            print ('\n!!! NewtonX: control.dyn not found !!!')
            print (usage)
            print ('!!! NewtonX: control.dyn not found !!!')
            exit()
        if os.path.exists('bagelinput.basis.inp') == False:
            print ('\n!!! Bagel: bagelinput.basis.inp not found !!!')
            print (usage)
            print ('!!! Bagel: bagelinput.basis.inp not found !!!')
            exit()
        if os.path.exists('bagelinput.part1.inp') == False or os.path.exists('bagelinput.part2.inp') == False or os.path.exists('bagelinput.part3.inp') == False:
            print ('\n!!! Bagel: bagelinput.part1-3.inp not found !!!')
            print (usage)
            print ('!!! Bagel: bagelinput.part1-3.inp not found !!!')
            exit()
    elif prog == 'pyrai2mdnn':
        if os.path.exists('input') == False:
            print ('\n!!! PyRAI2MD input not found !!!')
            print (usage)
            print ('!!! PyRAI2MD input not found !!!')
            exit()
    elif prog == 'pyrai2mdmolcas' or prog == 'pyrai2mdhybrid':
        if os.path.exists('input') == False:
            print ('\n!!! PyRAI2MD input not found !!!')
            print (usage)
            print ('!!! PyRAI2MD input not found !!!')
            exit()
        if os.path.exists('%s.molcas' % (input)) == False:
            print ('\n!!! PyRAI2MD molcas template not found !!!')
            print (usage)
            print ('!!! PyRAI2MD molcas template not found !!!')
            exit()
        if os.path.exists('%s.StrOrb' % (input)) == False and os.path.exists('%s.JobIph' % (input)) == False:
            print ('\n!!! PyRAI2MD molcas orbital file StrOrb or JobIph not found !!!')
            print (usage)
            print ('!!! PyRAI2MD molcas orbital file StrOrb or JobIph not found !!!')
            exit()

    else:
        print ('\n!!! Program %s not found !!!' % (prog))
        print (usage)
       	print ('!!! Program %s not found !!!' % (prog))
        exit()

    ## import initial condition sampling module
    sys.path.append('%s/bin/' % (pyqd))
    from dynamixsampling import Sampling

    callsample=['molden','bagel','g16','orca'] 
    skipsample=['newtonx','xyz']
    nesmb=int(sljb*slnd)

    if   format in callsample:
        print ("""
        Start Sampling
  ------------------------------------------------
    Seed         %10d
    Method       %10s
    Trajectories %10d = %3d Nodes X %3d Jobs
    Temperature  %10.2f    
    """ % (iseed,dist,nesmb,slnd,sljb,temp))

    elif format in skipsample:
        print ("""
        Read Sampled Initial Conditions

    """)


    ensemble=Sampling(input,nesmb,iseed,temp,dist,format) #generate initial conditions

    print ("""

    Additional Info
  ------------------------------------------------
    Restart after the first run:  %s (Only for Molcas)
    Select initial excited state: %s (Only for Molcas)
    """ % (restart,initex))

    # These functions take too many arguments, which can be optimized by making a dictionary. I'm just too lazy to do it:)
    if prog   == 'molcas':
        gen_molcas(ensemble,input,slpt,sltm,slmm,slnd,slcr,sljb,slin,restart,initex,tomlcs,format)
    elif prog == 'nxbagel':
        gen_nxbagel(ensemble,input,slpt,sltm,slmm,slnd,slcr,sljb,slin,tontx,tobgl,lbbls,lblpk,lbslp,lbbst,tomkl,tompi)
    elif prog == 'pyrai2mdnn':
        gen_pyrai2md(ensemble,input,slpt,sltm,slmm,slnd,slcr,sljb,slin,restart,initex,topyd,topython,'nn')
    elif prog == 'pyrai2mdmolcas':
        gen_pyrai2md(ensemble,input,slpt,sltm,slmm,slnd,slcr,sljb,slin,restart,initex,topyd,topython,'molcas',format)
    elif prog == 'pyrai2mdhybrid':
        gen_pyrai2md(ensemble,input,slpt,sltm,slmm,slnd,slcr,sljb,slin,restart,initex,topyd,topython,'hybrid',format)

def gen_molcas(ensemble,input,slpt,sltm,slmm,slnd,slcr,sljb,slin,restart,initex,tomlcs,format):
    ## This function will group Mocas calculations to individule runset
    ## this function will call molcas_batch and molcas to prepare files

    marks=[]
    if os.path.exists('%s.basis' % (input)) == True:
        with open('%s.basis' % (input)) as atommarks:
            marks=atommarks.read().splitlines()
            natom=int(marks[0])
            marks=marks[2:2+natom]

    in_temp=open('%s.inp' % (input),'r').read()
    if os.path.exists('%s.StrOrb' % (input)):
        in_orb=open('%s.StrOrb' % (input),'r').read()
    else:
        in_orb=None
    if os.path.exists('%s.JobIph' % (input)):
        jobiph='%s.JobIph' % (input)
    else:
        jobiph=None
    if os.path.exists('%s.key' % (input)):
        qmmmkey='%s.key' % (input)
    else:
        qmmmkey=None

    in_path=os.getcwd()

    runall=''
    for j in range(slnd):
        start=slin+j*sljb
        end=start+sljb-1
        for i in range(sljb):
            if format != 'xz':
                in_xyz,in_velo=Unpack(ensemble[i+j*sljb],'molcas')  # unpack initial condition to xyz and velocity
                if len(marks)>0:
                    in_xyz=Markatom(in_xyz,marks)
            else:
                in_xyz,in_velo=UnpackXZ(ensemble[i+j*sljb])
            inputname='%s-%s' % (input,i+start)
            inputpath='%s/%s' % (in_path,inputname)
            molcas(input,inputname,inputpath,slmm,in_temp,in_path,in_orb,jobiph,qmmmkey,in_xyz,in_velo,tomlcs)  #prepare calculations
            sys.stdout.write('Setup Calculation: %.2f%%\r' % ((i+j*sljb+1)*100/(sljb*slnd)))
        batch=molcas_batch(input,j,start,end,in_path,slcr,sltm,slpt,slmm,tomlcs,initex,restart)
        run=open('./runset-%d.sh' % (j+1),'w')
        run.write(batch)
        run.close()
        os.system("chmod 777 runset-%d.sh" % (j+1))
        runall+='sbatch runset-%d.sh\n' % (j+1)
    all=open('./runall.sh','w')
    all.write(runall)
    all.close()
    os.system("chmod 777 runall.sh")
    print('\n\n Done\n')

def molcas_batch(input,j,start,end,in_path,slcr,sltm,slpt,slmm,tomlcs,initex,restart):
    ## This function will be called by gen_molcas function
    ## This function generates runset for MolCas calculation


    ## copy neccessary module for to Molcas Tools
    if restart != 0 and os.path.exists('%s/Tools/TSHrestart.py' % (tomlcs)) == False:
        shutil.copy2('%s/bin/TSHrestart.py' % (pyqd),'%s/Tools/TSHrestart.py' % (tomlcs))
    if initex !=0 and os.path.exists('%s/Tools/InitEx.py' % (tomlcs)) == False:
        shutil.copy2('%s/bin/InitEx.py' % (pyqd),'%s/Tools/InitEx.py' % (tomlcs))

    if initex == 0:
        pri=''
    else:
        pri="""
  if [ "$STEP" == "0" ]
  then
    python3 $MOLCAS/Tools/InitEx.py prep $INPUT.inp
    cd INITEX
    $MOLCAS/bin/pymolcas -f $INPUT.inp -b 1
    cd ../  
    python3 $MOLCAS/Tools/InitEx.py read $INPUT.inp 
  fi
"""

    if restart == 0:
        add=' '
        addend=' '
    else:
        add="""
  python3 $MOLCAS/Tools/TSHrestart.py PROG
""" 
        addend="""
STEP=`tail -n1 PROG|awk '{print $1}'`
if [ "$STEP" -lt "$MAX" ]
then
cd ../
sbatch runset-%d.sh
fi
""" % (j+1)

    batch="""#!/bin/sh
## script for OpenMalCas
#SBATCH --nodes=1
#SBATCH --ntasks=%d
#SBATCH --time=%s
#SBATCH --job-name=%s-%d
#SBATCH --partition=%s
#SBATCH --mem=%dmb
#SBATCH --output=%%j.o.slurm
#SBATCH --error=%%j.e.slurm
#SBATCH --constraint="[broadwell|haswell|cascadelake|ib]"

export MOLCAS_NPROCS=1
export MOLCAS_MEM=%d
export MOLCAS_PRINT=2
export OMP_NUM_THREADS=1
export MOLCAS=%s
export PATH=$MOLCAS/bin:$PATH

echo $SLURM_JOB_NAME

if [ -d "/srv/tmp" ]
then
 export LOCAL_TMP=/srv/tmp
else
 export LOCAL_TMP=/tmp
fi

module load python/3.7.1

RunMolcas(){
  if [ -a "PROG" ]
  then
    STEP=`tail -n1 PROG|awk '{print $1}'`
  else
    STEP=0
    echo "$INPUT $MAX" > PROG
    echo "$STEP" >> PROG
  fi
%s
  $MOLCAS/bin/pymolcas -f $INPUT.inp -b 1
  STEP=`expr $STEP + 1`
  echo "$STEP" >> PROG
%s
  rm -r $MOLCAS_WORKDIR/$MOLCAS_PROJECT
}

MAX=%s

for ((i=%d;i<=%d;i++))
do
  export INPUT="%s-$i"
  export WORKDIR="%s/%s-$i"
  export MOLCAS_PROJECT=$INPUT
  export MOLCAS_WORKDIR=$LOCAL_TMP/$USER/$SLURM_JOB_ID
  mkdir -p $MOLCAS_WORKDIR/$MOLCAS_PROJECT
  cd $WORKDIR
  RunMolcas &
  sleep 5
done
wait
rm -r $MOLCAS_WORKDIR

%s
""" % (slcr,sltm,input,j+1,slpt,int(slmm*slcr*1.333),slmm,tomlcs,pri,add,restart+1,start,end,input,in_path,input,addend)

    return batch

def molcas(input,inputname,inputpath,slmm,in_temp,in_path,in_orb,jobiph,qmmmkey,in_xyz,in_velo,tomlcs):    
    ## This function prepares MolCas calculation
    ## It generates .inp .StrOrb .xyz .velocity.xyz
    ## This function generates a backup slurm batch file for each calculation

    if os.path.exists('%s' % (inputpath)) == False:
        os.makedirs('%s' % (inputpath))

    runscript="""#!/bin/sh
## backup script for OpenMalCas
## $INPUT and $WORKDIR do not belong to OpenMolCas
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=23:50:00
#SBATCH --job-name=%s
#SBATCH --partition=short
#SBATCH --mem=%dmb
#SBATCH --output=%%j.o.slurm
#SBATCH --error=%%j.e.slurm
#SBATCH --constraint="[broadwell|haswell|cascadelake|ib]"

if [ -d "/srv/tmp" ]
then
 export LOCAL_TMP=/srv/tmp
else
 export LOCAL_TMP=/tmp
fi

export INPUT=%s
export WORKDIR=%s
export MOLCAS_NPROCS=1
export MOLCAS_MEM=%d
export MOLCAS_PRINT=2
export MOLCAS_PROJECT=$INPUT
export OMP_NUM_THREADS=1
export MOLCAS=%s
export TINKER=$MOLCAS/tinker-6.3.3/bin
export MOLCAS_WORKDIR=$LOCAL_TMP/$USER/$SLURM_JOB_ID
export PATH=$MOLCAS/bin:$PATH

module load python/3.7.1

mkdir -p $MOLCAS_WORKDIR/$MOLCAS_PROJECT
cd $WORKDIR
$MOLCAS/bin/pymolcas -f $INPUT.inp -b 1
rm -r $MOLCAS_WORKDIR
""" % (inputname,int(slmm*1.333),inputname,inputpath,slmm,tomlcs)

    out=open('%s/%s.inp'           % (inputpath,inputname),'w')
    out.write(in_temp)
    out.close()
    out=open('%s/%s.sh'            % (inputpath,inputname),'w')
    out.write(runscript)
    out.close()

    if in_orb != None:
        out=open('%s/%s.StrOrb'        % (inputpath,inputname),'w')
        out.write(in_orb)
        out.close()

    if jobiph != None:
        shutil.copy2(jobiph,'%s/%s.JobIph' % (inputpath, inputname))

    if qmmmkey != None:
        shutil.copy(qmmmkey, '%s/%s.key' % (inputpath, inputname))

    out=open('%s/%s.xyz'           % (inputpath,inputname),'w')
    out.write(in_xyz)
    out.close()
    out=open('%s/%s.velocity.xyz'  % (inputpath,inputname),'w')
    out.write(in_velo)
    out.close()
    os.system("chmod 777 %s/%s.sh" % (inputpath,inputname))

def Markatom(xyz,marks):
    ## This function marks atoms for different basis set specification of Molcas
    ## here xyz, e, x, y, z are all strings

    xyz=xyz.splitlines()
    new_xyz='%s\n%s\n' % (xyz[0],xyz[1])
    for n,line in enumerate(xyz[2:]):
        e,x,y,z=line.split()[0:4]
        e = marks[n].split()[0]
        new_xyz+='%-5s%30s%30s%30s\n' % (e,x,y,z)

    return new_xyz

def gen_nxbagel(ensemble,input,slpt,sltm,slmm,slnd,slcr,sljb,slin,tontx,tobgl,lbbls,lblpk,lbslp,lbbst,tomkl,tompi):
    ## This function groups nxbagel calculations in individule runset
    ## This function will call nxbagel_batch and nxbagel to prepare files
    ## This function need control.dyn for NxBagel calculation, therm.inp, sh.inp and jiri.inp are optional

    in_jiri=''
    in_sh=''
    in_orb=''
    in_bagel0=open('bagelinput.basis.inp','r').read()
    in_bagel1=open('bagelinput.part1.inp','r').read()
    in_bagel2=open('bagelinput.part2.inp','r').read()
    in_bagel3=open('bagelinput.part3.inp','r').read()
    in_control=open('control.dyn','r').read()

    if os.path.exists('%s.archive' % (input)) == True:
        in_orb='%s.archive' % (input)
    if os.path.exists('jiri.inp') == True:
        in_jiri=open('jiri.inp','r').read()
    if os.path.exists('sh.inp') == True:
        in_sh=open('sh.inp','r').read()
    if os.path.exists('therm.inp') == True:
        in_therm=open('therm.inp','r').read()

    in_temp={ 
    'jiri'   : in_jiri,
    'sh'     : in_sh,
    'therm'  : in_therm,
    'orb'    : in_orb,
    'control': in_control,
    'basis'  : in_bagel0,
    'part1'  : in_bagel1,
    'part2'  : in_bagel2,
    'part3'  : in_bagel3
    }

    in_path=os.getcwd()
    runall=''

    for j in range(slnd):
        start=slin+j*sljb
        end=start+sljb-1
        for i in range(sljb):
            in_xyz,in_velo=Unpack(ensemble[i+j*sljb],'newton')  # unpack initial condition to xyz and velocity
            inputname='%s-%s' % (input,i+start)
            inputpath='%s/%s' % (in_path,inputname)
            nxbagel(input,inputname,inputpath,slcr,sljb,sltm,slpt,slmm,in_temp,in_xyz,in_velo,tontx,tobgl,lbbls,lblpk,lbslp,lbbst,tomkl,tompi) #prepare calculations
            sys.stdout.write('Setup Calculation: %.2f%%\r' % ((i+j*sljb+1)*100/(sljb*slnd)))
        batch=nxbagel_batch(input,j,start,end,in_path,slcr,sljb,sltm,slpt,slmm,tontx,tobgl,lbbls,lblpk,lbslp,lbbst,tomkl,tompi)
        run=open('./runset-%d.sh' % (j+1),'w')
        run.write(batch)
        run.close()
        os.system("chmod 777 runset-%d.sh" % (j+1))
        runall+='sbatch runset-%d.sh\n' % (j+1)
    all=open('./runall.sh','w')
    all.write(runall)
    all.close()
    os.system("chmod 777 runall.sh")
    print('\n\n Done\n')

def nxbagel_batch(input,j,start,end,in_path,slcr,sljb,sltm,slpt,slmm,tontx,tobgl,lbbls,lblpk,lbslp,lbbst,tomkl,tompi):
    ## This function will be called by gen_nxbagel
    ## This function generates runset for NxBagel calculation

    bagelpal=int(slcr/sljb) ### Note, this line doesn't check if slcr and sljb are appropriate for parallelization, be careful!!!

    batch="""#!/bin/sh
## script for NX/BAGEL
#SBATCH --nodes=1
#SBATCH --ntasks=%d
#SBATCH --time=%s
#SBATCH --job-name=%s-%d
#SBATCH --partition=%s
#SBATCH --mem=%dmb
#SBATCH --output=%%j.o.slurm
#SBATCH --error=%%j.e.slurm

export MKL_DEBUG_CPU_TYPE=5

export BAGEL_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BAGELPAL=%s

export NX=%s/bin
export BAGEL=%s/bin/BAGEL
export BAGEL_LIB=%s/lib
export BLAS_LIB=%s
export LAPACK_LIB=%s
export SCALAPACK_LIB=%s
export BOOST_LIB=%s/lib
source %s/bin/mklvars.sh intel64
export MPI=%s
export LD_LIBRARY_PATH=$MPI/lib:$BAGEL_LIB:$BLAS_LIB:$LAPACK_LIB:$SCALAPACK_LIB:$BOOST_LIB:$LD_LIBRARY_PATH
export PATH=$MPI/bin:$PATH

echo $SLURM_JOB_NAME

for ((i=%d;i<=%d;i++))
do
  export WORKDIR=%s/%s-$i
  cd $WORKDIR
  $NX/moldyn.pl > moldyn.log &
  sleep 5
done
wait
""" % (slcr,sltm,input,j+1,slpt,int(slmm*slcr*1.1),bagelpal,tontx,tobgl,tobgl,lbbls,lblpk,lbslp,lbbst,tomkl,tompi,start,end,in_path,input)

    return batch
def nxbagel(input,inputname,inputpath,slcr,sljb,sltm,slpt,slmm,in_temp,in_xyz,in_velo,tontx,tobgl,lbbls,lblpk,lbslp,lbbst,tomkl,tompi):
    ## This function prepares NxBagel calculation
    ## It generates TRAJECTORIES/TRAJ#/JOB_NAD, geom, veloc,controd.dyn for NxBagel calculations
    ## This function generates a backup slurm batch file for each calculation

    bagelpal=int(slcr/sljb) ### Note, this line doesn't check if slcr and sljb are appropriate for parallelization, be careful!!!

    if os.path.exists('%s/JOB_NAD' % (inputpath)) == False:
        os.makedirs('%s/JOB_NAD' % (inputpath))

    runscript="""#!/bin/sh
## script for NX/BAGEL
#SBATCH --nodes=1
#SBATCH --ntasks=%s
#SBATCH --time=%s
#SBATCH --job-name=%s
#SBATCH --partition=%s
#SBATCH --mem=%dmb
#SBATCH --output=%%j.o.slurm
#SBATCH --error=%%j.e.slurm
#SBATCH --constraint="[broadwell|haswell|cascadelake|ib]"

export BAGEL_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BAGELPAL=$SLURM_NTASKS

export NX=%s/bin
export BAGEL=%s/bin/BAGEL
export BAGEL_LIB=%s/lib
export BLAS_LIB=%s
export LAPACK_LIB=%s
export SCALAPACK_LIB=%s
export BOOST_LIB=%s/lib
source %s/bin/mklvars.sh intel64
export MPI=%s
export LD_LIBRARY_PATH=$MPI/lib:$BAGEL_LIB:$BLAS_LIB:$LAPACK_LIB:$SCALAPACK_LIB:$BOOST_LIB:$LD_LIBRARY_PATH
export PATH=$MPI/bin:$PATH

export WORKDIR=%s
cd $WORKDIR
$NX/moldyn.pl > moldyn.log
""" % (bagelpal,sltm,inputname,slpt,int(slmm*1.1),tontx,tobgl,tobgl,lbbls,lblpk,lbslp,lbbst,tomkl,tompi,inputpath)

    out=open('%s/%s.sh'            % (inputpath,inputname),'w')
    out.write(runscript)
    out.close()
    if in_temp['orb'] != '':
        shutil.copy2('%s' % (in_temp['orb']), '%s/%s' % (inputpath,in_temp['orb']))
    if in_temp['jiri'] != '':
        out=open('%s/jiri.inp'         % (inputpath),'w')
        out.write(in_temp['jiri'])
        out.close()
    if in_temp['sh'] != '':
        out=open('%s/sh.inp'           % (inputpath),'w')
        out.write(in_temp['sh'])
        out.close()
    if in_temp['therm'] != '':
        out=open('%s/therm.inp'           % (inputpath),'w')
        out.write(in_temp['therm'])
        out.close()

    out=open('%s/control.dyn'      % (inputpath),'w')
    out.write(in_temp['control'])
    out.close()
    out=open('%s/JOB_NAD/bagelinput.basis.inp' % (inputpath),'w')
    out.write(in_temp['basis'])
    out.close()
    out=open('%s/JOB_NAD/bagelinput.part1.inp' % (inputpath),'w')
    out.write(in_temp['part1'])
    out.close()
    out=open('%s/JOB_NAD/bagelinput.part2.inp' % (inputpath),'w')
    out.write(in_temp['part2'])
    out.close()
    out=open('%s/JOB_NAD/bagelinput.part3.inp' % (inputpath),'w')
    out.write(in_temp['part3'])
    out.close()
    out=open('%s/geom'  % (inputpath),'w')
    out.write(in_xyz)
    out.close()
    out=open('%s/veloc' % (inputpath),'w')
    out.write(in_velo)
    out.close()
    os.system("chmod 777 %s/%s.sh" % (inputpath,inputname))


def Unpack(ensemble,prog):  
    ## This function unpacks initial condition to xyz and velocity

    xyz=''
    velo=''
    if   prog == 'molcas':
        natom=int(len(ensemble))
        xyz='%d\n\n' % (natom)
        for i in ensemble:
            xyz+='%-5s%30.16f%30.16f%30.16f\n' % (i[0],float(i[1]),float(i[2]),float(i[3]))
            velo+='%30.16f%30.16f%30.16f\n'    % (float(i[4]),float(i[5]),float(i[6]))
    elif prog == 'newton':
        for i in ensemble:
            xyz+='%-5s%6.1f%30.16f%30.16f%30.16f%30.16f\n' % (i[0],float(i[8]),float(i[1])*1.88973,float(i[2])*1.88973,float(i[3])*1.88973,float(i[7]))
            velo+='%30.16f%30.16f%30.16f\n'    % (float(i[4]),float(i[5]),float(i[6]))

    return xyz,velo

def gen_pyrai2md(ensemble,input,slpt,sltm,slmm,slnd,slcr,sljb,slin,restart,initex,topyd,topython,qm,format):
    ## This function will group PyRAI2MD calculations to individule runset
    ## this function will call pyrai2md_batch and pyrai2md to prepare files

    in_temp=open('input','r').read()
    in_path=os.getcwd()
    runall=''
    runall2=''
    for j in range(slnd):
        start=slin+j*sljb
        end=start+sljb-1
        for i in range(sljb):
            if format != 'xz':
                in_xyz,in_velo=Unpack(ensemble[i+j*sljb],'molcas')  # unpack initial condition to xyz and velocity
            else:
                in_xyz,in_velo=UnpackXZ(ensemble[i+j*sljb])
            inputname='%s-%s' % (input,i+start)
            inputpath='%s/%s' % (in_path,inputname)
            pyrai2md(input,inputname,inputpath,slcr,sljb,sltm,slpt,slmm,in_temp,in_xyz,in_velo,topyd,topython,qm)   #prepare calculations
            sys.stdout.write('Setup Calculation: %.2f%%\r' % ((i+j*sljb+1)*100/(sljb*slnd)))
            runall2+='cd %s\nsbatch run_PyRAI2MD.sh\n' % (inputpath)
        batch=pyrai2md_batch(input,j,start,end,in_path,slcr,sljb,sltm,slpt,slmm,topyd,topython)
        run=open('./runset-%d.sh' % (j+1),'w')
        run.write(batch)
        run.close()
        os.system("chmod 777 runset-%d.sh" % (j+1))
        runall+='sbatch runset-%d.sh\n' % (j+1)
    all=open('./runall.sh','w')
    all.write(runall)
    all.close()
    os.system("chmod 777 runall.sh")
    all=open('./runall2.sh','w')
    all.write(runall2)
    all.close()
    os.system("chmod 777 runall2.sh")

    print('\n\n Done\n')

def pyrai2md_batch(input,j,start,end,in_path,slcr,sljb,sltm,slpt,slmm,topyd,topython):
    ## This function will be called by gen_pyrai2md
    ## This function generates runset for PyRAI2MD calculation

    bagelpal=int(slcr/sljb) ### Note, this line doesn't check if slcr and sljb are appropriate for parallelization, be careful!!!

    batch="""#!/bin/sh
## script for PyRAI2MD
#SBATCH --nodes=1
#SBATCH --cpus-per-task=%d
#SBATCH --time=%s
#SBATCH --job-name=%s-%d
#SBATCH --partition=%s
#SBATCH --mem=%dmb
#SBATCH --output=%%j.o.slurm
#SBATCH --error=%%j.e.slurm
#SBATCH --constraint="[broadwell|haswell|cascadelake|ib]"

export INPUT=input
export PYRAI2MD=%s
export PATH=%s/bin:$PATH
export LD_LIBRARY_PATH=%s/lib:$LD_LIBRARY_PATH

for ((i=%d;i<=%d;i++))
do
  export WORKDIR=%s/%s-$i
  cd $WORKDIR
  python3 $PYRAI2MD/pyrai2md.py $INPUT &
  sleep 5
done
wait
""" % (slcr,sltm,input,j+1,slpt,int(slmm*slcr*1.1),topyd,topython,topython,start,end,in_path,input)

    return batch

def pyrai2md(input,inputname,inputpath,slcr,sljb,sltm,slpt,slmm,in_temp,in_xyz,in_velo,topyd,topython,qm):
    ## This function prepares PyRAI2MD calculation
    ## This function generates a backup slurm batch file for each calculation

    ncpu=int(slcr/sljb) ### Note, this line doesn't check if slcr and sljb are appropriate for parallelization, be careful!!!

    if os.path.exists('%s' % (inputpath)) == False:
        os.makedirs('%s' % (inputpath))

    in_temp=update_pyrai2md_input(in_temp,inputname)

    if os.path.exists('%s/%s.molcas' % (inputpath,input)) == False and (qm == 'molcas' or qm == 'hybrid'):
        shutil.copy2('%s.molcas' % (input),'%s/%s.molcas' % (inputpath,inputname))

    if  os.path.exists('%s.StrOrb' % (input)) == True and\
        os.path.exists('%s/%s.StrOrb' % (inputpath,input)) == False and\
        (qm == 'molcas' or qm == 'hybrid'):
        shutil.copy2('%s.StrOrb' % (input),'%s/%s.StrOrb' % (inputpath,inputname))

    if  os.path.exists('%s.JobIph' % (input)) == True and\
        os.path.exists('%s/%s.JobIph' % (inputpath,input)) == False and\
        (qm == 'molcas' or qm == 'hybrid'):
        shutil.copy2('%s.JobIph' % (input),'%s/%s.JobIph' % (inputpath,inputname))

    runscript="""#!/bin/sh
## script for PyRAI2MD
#SBATCH --nodes=1
#SBATCH --cpus-per-task=%s
#SBATCH --time=%s
#SBATCH --job-name=%s
#SBATCH --partition=%s
#SBATCH --mem=%dmb
#SBATCH --output=%%j.o.slurm
#SBATCH --error=%%j.e.slurm
#SBATCH --constraint="[broadwell|haswell|cascadelake|ib]"

export INPUT=input
export WORKDIR=%s
export PYRAI2MD=%s
export PATH=%s/bin:$PATH
export LD_LIBRARY_PATH=%s/lib:$LD_LIBRARY_PATH

cd $WORKDIR
python3 $PYRAI2MD/pyrai2md.py $INPUT

""" % (ncpu,sltm,inputname,slpt,int(slmm*1.1),inputpath,topyd,topython,topython)

    out=open('%s/run_PyRAI2MD.sh'            % (inputpath),'w')
    out.write(runscript)
    out.close()

    out=open('%s/input'	   % (inputpath),'w')
    out.write(in_temp)
    out.close()

    out=open('%s/%s.xyz'           % (inputpath,inputname),'w')
    out.write(in_xyz)
    out.close()

    out=open('%s/%s.velo'           % (inputpath,inputname),'w')
    out.write(in_velo)
    out.close()

    os.system("chmod 777 %s/run_PyRAI2MD.sh" % (inputpath))


def update_pyrai2md_input(in_temp,inputname):
    ## first edit input - change the title
    ## second copy neural network
    input=''
    for line in in_temp.splitlines():
        if 'title' in line:
            input+='title %s\n' % (inputname)
        else:
            input+='%s\n' % (line)

    return input

def UnpackXZ(ensemble):
    xyz = ensemble['txyz']
    velo = ensemble['velo']

    xyz = '\n'.join(xyz) + '\n'
    velo = '\n'.join(velo) + '\n'

    return xyz, velo

if __name__ == '__main__':
    main()


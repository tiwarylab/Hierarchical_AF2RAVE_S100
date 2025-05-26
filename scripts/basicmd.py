import numpy as np
from openmmplumed import PlumedForce
from openmm.app import *
from openmm import *
from openmm.unit import *
import pdbfixer
import openmm_utils as op            # openmm functions
import S100CVs as scv
import sys 
import os
import mdtraj as md

oss=os.system


if __name__=='__main__':
    '''
    Python file to perform Molecular dynamics using openmm - basic example to use openmm_utils package
    
        (A better script will be updated)
    
    Args:
    
    Input:
     -file            : PDB file to start
     -nsteps          : Number simulation steps (default=25000000)
     -temp            : Temperature of simulation (default=300K)
     -pressure        : Pressure of simulation (default=1bar)
     -plumed          : Plumed file to bias or calculate CVs on the fly
    
    Output:
      fixed_1.pdb     : PDB after ficing the missing residues and terminal atoms
      fixedH_1.pdb    : PDB after adding hydrogen
      solvated_1.pdb  : After solvating adding ions
      minim_1.pdb     : After running a small energy minimization step
      mdout_1.pdb     : final structure file after runnung 'run_NVT()' once
      mdlog_1.txt     : Log file after runnung 'run_NVT()' once
      mdout_1.dcd     : trajectory file after runnung 'run_NVT()' once
    '''
    
    if '-file' in sys.argv:
        pdbf = sys.argv[sys.argv.index('-file')+1]
        flag=1
    else:
        print("An initial coordinate file is required")
        flag=0

    if '-nsteps' in sys.argv:
        nsteps = int(sys.argv[sys.argv.index('-nsteps')+1])
    elif 'tns' in sys.argv:
        nsteps = int(sys.argv[sys.argv.index('-tns')+1])*500000
    else:
        nsteps=5000000

    if '-temp' in sys.argv:
        temp = int(sys.argv[sys.argv.index('-temp')+1])   #Integer value
    else:
        temp = 300              # 300K temperature
    
    if '-press' in sys.argv:
        pressure = float(sys.argv[sys.argv.index('-press')+1])   #float value
    else:
        pressure = 1.0              # 1.0 bar pressure
    
       
    if flag:
            
            oss("cp %s init.pdb"%pdbf)
            
            positions,topology = op.fix_pdb("init.pdb")                             # Fixes pdb with misiing atoms and terminal atoms
            
            #ff = ForceField('charmm36.xml', 'charmm36/water.xml')         # Creates forcefield
            #ff = ForceField('amber99sbildn.xml', 'tip3p.xml')              # Creates forcefield - for DDR1 kinase
            ff = ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', "amber/tip3p_IOD_multivalent.xml")

            positions,topology = op.add_hydrogen(positions,topology,ff)                          # Adds hydrogen according to the forcefield
            positions,topology = op.solvate_me(positions,topology,ff,True,1,'tip3p')     # (Default) neutralizes the charge with Na+ and Cl- ions

            #Create simulation
            simulation=op.get_LangevinM_system(topology,ff,temp,0.002)


            # Perform equilibration

            # Adding restrain
            simulation = op.add_pos_res(positions,topology,simulation,10)    # 10 KJ/(mol.A^2) contrain

            # Minimize
            positions,simulation = op.ener_Minimize(positions,simulation)

            # NVT - restrain
            positions,velocities,simulation=op.run_MD(positions,simulation,150000,ens='NVT',run_type='restrainedNVT')   #No '_' in run_type
            
            positions,velocities,simulation=op.run_MD(positions,simulation,150000,ens='NPT',run_type='restrainedNPT',\
                                                   velocities=velocities,cont=True)                            

            # Remove restrain
            simulation.context.getSystem().removeForce(simulation.context.getSystem().getNumForces()-1)
            simulation.context.reinitialize(preserveState=True)

            # inputs for ions restrain
            traj=md.load("restrainedNPT_0.pdb")
            chelates = scv.ion_oxy_pairs(traj.top)

            # NVT - no restrain
            positions,velocities,simulation=op.run_MD(positions,simulation,250000,ens='NVT',run_type='equilNVT',\
                                                   velocities=velocities,cont=True,chelates=chelates)                            
            
            # NPT - no restrain
            positions,velocities,simulation=op.run_MD(positions,simulation,500000,ens='NPT',run_type='equilNPT',\
                                                   velocities=velocities,cont=True,chelates=chelates)                            


            # NPT - no restrain prodrun
            positions,velocities,simulation=op.run_MD(positions,simulation,nsteps,ens='NPT',run_type='prodrun_NPT',\
                                                    velocities=velocities,cont=True,chelates=chelates)



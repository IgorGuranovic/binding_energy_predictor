from rosetta import *
from pyrosetta import *
import glob

files = glob.glob("autodock_automation/PDB_Files/*")

for file in files:
    init()
    p = Pose()
    p = pose_from_pdb(file)

    scorefxn = get_fa_scorefxn()

    ncycles = 50
    kT = 1.0
    mc = MonteCarlo(p, scorefxn, kT)

    movemap = MoveMap()
    movemap.set_bb(True)

    small_mover = pyrosetta.rosetta.protocols.simple_moves.SmallMover(movemap, kT, 5)

    for i in range(1, ncycles):
        small_mover.apply(p)
        mc.boltzmann(p)

    mc.recover_low(p)
    file_temp = file.replace('PDB_Files','PDB_Files_Minimized')
    
    dump_pdb(p, file_temp)



import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj

file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")] #Select all PDB files in this folder
file_names=sorted(file_names)
#rc("open protein.gro")
i=0
ind=0
SEL=":36.A,42.A,43.A,44.A,45.A,47.A,52.A,55.A,56.A,59.A,60.A,76.A,83.A,84.A,87.A,0.B,8.B,11.B,12.B"
while ind<len(file_names):
        fn=file_names[ind]
        print("Loading:",fn)
        fnnew=fn.replace(".pdb","_charged.mol2")
        replyobj.status("Processing " + fn)  #show what file we're working on
        rc("open " + fn)
        rc("addcharge std chargeModel ff14SB")
        rc("write format mol2 #0 "+fnnew)
        rc("close all")
        ind+=1
	# save image to a file that ends in .png rather than .pdb
	#png_name = fn[:-3] + "png"
	#rc("copy file " + png_name + " supersample 3")
	#rc("close all")

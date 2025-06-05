import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj

file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")] #Select all PDB files in this folder
file_names=sorted(file_names)
#rc("open protein.gro")
i=0
ind=0
#SEL=":37.A,43.A,44.A,45.A,46.A,48.A,53.A,56.A,57.A,60.A,61.A,77.A,80.A,81.A,84.A,85.A,88.A,89.A,1.B,5.B,6.B,8.B,9.B,12.B,13.B"
#SEL=":36.A,42.A,43.A,44.A,45.A,47.A,48.A,49.A,52.A,55.A,56.A,60.A,76.A,79.A,80.A,83.A,84.A,87.A,0.B,11.B,12.B,13.B,5.B,7.B,8.B,9.B" # 16.B?
#SEL=":36.A,42.A,43.A,44.A,45.A,46.A,47.A,52.A,55.A,56.A,59.A,60.A,76.A,79.A,80.A,83.A,84.A,85.A,87.A,0.B,11.B,12.B,8.B"
#SEL=":37.A,43.A,44.A,45.A,46.A,47.A,48.A,53.A,56.A,57.A,60.A,61.A,77.A,80.A,81.A,84.A,85.A,86.A,88.A,1.B,12.B,13.B,9.B"
SEL=":37.A,41.A,43.A,44.A,45.A,46.A,47.A,48.A,53.A,56.A,57.A,58.A,59.A,60.A,61.A,74.A,77.A,80.A,81.A,84.A,85.A,86.A,88.A,89.A,1.B,12.B,13.B,5.B,8.B,9.B"
#SEL=":30.A,31.A,32.A,33.A,34.A,35.A,36.A,37.A,38.A,39.A,40.A,41.A,42.A,43.A,44.A,45.A,46.A,47.A,48.A,49.A,50.A,51.A,52.A,53.A,54.A,55.A,56.A,57.A,58.A,59.A,60.A,61.A,71.A,72.A,73.A,74.A,75.A,76.A,77.A,78.A,79.A,80.A,81.A,82.A,83.A,84.A,85.A,86.A,87.A,88.A,1.B,2.B,11.B,12.B,13.B,14.B,3.B,4.B,5.B,6.B,7.B,8.B,9.B,10.B"
while ind<len(file_names):
        fn=file_names[ind]
        print("Loading:",fn)
        fnnew=fn.replace(".pdb","_site.mol2")
        replyobj.status("Processing " + fn)  #show what file we're working on
        rc("open " + fn)
        rc("sel "+SEL)
        rc("write format mol2 relative #0 selected #0 "+fnnew)
        rc("close all")
        ind+=1
	# save image to a file that ends in .png rather than .pdb
	#png_name = fn[:-3] + "png"
	#rc("copy file " + png_name + " supersample 3")
	#rc("close all")

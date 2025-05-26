import sys
import numpy as np
import pandas as pd
import os
import glob

prot = sys.argv[1]
bid = int(sys.argv[2])
if len(sys.argv) == 4:
    op_id = [int(x) for x in np.loadtxt(sys.argv[3])]
    folder = f'{sys.argv[3][:2]}_{prot}_b{bid}'
else:
    op_id = [int(x) for x in range(49)]
    folder = f'{prot}_b{bid}'
print(op_id, 'cv#:', len(op_id))


listindices = []
df_traj = pd.read_csv('basin_traj.csv', header=0,sep=',')
df_tmp = df_traj[(df_traj.name.str.split('_').str[0] == prot) & (df_traj.basin == bid)]
print(df_tmp)
for i in range(len(df_tmp)):
    listindices.append(f"/home/xg23/scratch/S100/rave/{df_tmp.iloc[i]['name']}/pk_Dist{df_tmp.iloc[i]['chain']}.txt")
print(listindices, 'example data shape:', np.loadtxt(listindices[0]).shape)

split = 5
num_state = len(listindices)*split

os.mkdir(folder)
os.chdir(folder)

input_folder = f'{prot}_b{bid}_input'
os.mkdir(input_folder)


cvs = []
for i in range(len(listindices)):
    cvs.append(np.loadtxt(listindices[i])[:,op_id])
    lentraj=len(cvs[i])
    state_indices = np.hstack([np.full(int(lentraj / split), ii, dtype=np.int8) for ii in range(split)])
    if len(state_indices) < lentraj:
        state_indices = np.hstack([state_indices, np.full(lentraj - len(state_indices), split - 1, dtype=np.int8)])
    initlabels = np.eye(num_state)[state_indices + int(i * split)]
    np.save(f"{input_folder}/labels_{i}_unb.npy",initlabels)

# min-max scaling
op_max = np.max(np.concatenate(cvs), axis=0)
op_min = np.min(np.concatenate(cvs), axis=0)
np.save(f"{input_folder}/max_unb.npy",op_max)
np.save(f"{input_folder}/min_unb.npy",op_min)
for i in range(len(listindices)):
    np.save(f"{input_folder}/colvar_{i}_unb.npy", (cvs[i]-op_min)/(op_max-op_min))


curr = os.getcwd()
f_nodt=open(f"/home/xg23/scratch/S100/rave/pocket_spib/sample_config.ini")
f=open(f"config.ini","w")
lines=f_nodt.readlines()
f.writelines(lines)
f_nodt.close()
unbpath = f"{curr}/{input_folder}"
f.write("\n traj_data = [%s]\n"%",".join(["%s/colvar_%i_unb.npy"%(unbpath,i) for i in range(len(cvs))]))
f.write("\n initial_labels = [%s]\n"%",".join(["%s/labels_%i_unb.npy"%(unbpath,i) for i in range(len(cvs))]))
f.write("\n traj_weights \n")
f.close()

f=open(f"{unbpath}/listindices.txt","w")
f.write(f"{listindices}")
f.close()

for i in range(len(listindices)):
    #os.system(f'cp {listindices[i]} {unbpath}/CVs_{i}.txt')
    np.savetxt(f'{unbpath}/CVs_{i}.txt', cvs[i], fmt='%.5f')

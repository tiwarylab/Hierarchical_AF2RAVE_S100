#!/bin/bash
#

GRIDFOL="grids"
SPHFOL="spheres"
TEMPFOL="temp_gridgen"
PREFIX="S100B_b1_"
SUFFIX="_ca.pdb"

mkdir -p $GRIDFOL


rm $TEMPFOL/*.pdb -f
for f in `ls *.pdb`
do
        echo "Processing $f"
        cp $f $TEMPFOL/rec.pdb
	cp `echo $f|sed s/".pdb"/"_charged.mol2"/` $TEMPFOL/rec_charged.mol2
        fid=`echo $f|sed s/$PREFIX//g|sed s/$SUFFIX//g`
        ls $SPHFOL/$fid"_sel.sph"
	cp $SPHFOL/$fid"_sel.sph" $TEMPFOL/region_spheres.sph
        cd $TEMPFOL
        ./script_demo

        cd ..
        cp $TEMPFOL/grid.nrg $GRIDFOL/"$fid""_grid.nrg"
        cp $TEMPFOL/grid.bmp $GRIDFOL/"$fid""_grid.bmp"
        echo "Wrote grid to $GRIDFOL/"$fid"_grid"
        cd $TEMPFOL
        ./script_clean
        rm *.pdb region_spheres.sph rec_charged.mol2 -f
        cd ..
done

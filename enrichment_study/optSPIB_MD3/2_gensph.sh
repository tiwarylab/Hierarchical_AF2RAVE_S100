#!/bin/bash
#
SPHFOL="spheres"
MSFOL="MS_spheres"
TEMPFOL="temp_sphgen"
PREFIX="S100B_b1_"
SUFFIX="_ca.pdb"

mkdir -p $SPHFOL
mkdir -p "$TEMPFOL"

rm $TEMPFOL/*.pdb
for f in `ls *.pdb`
do
	echo "Processing $f"
	cp $f $TEMPFOL/rec.pdb
	fid=`echo $f|sed s/$PREFIX//g|sed s/$SUFFIX//g`
	echo $fid
	cp $MSFOL/"$fid".ms $TEMPFOL/rec.ms
	sitefile=`echo $f|sed s/".pdb"//g`
	echo $sitefile
	cp $sitefile"_site.mol2" $TEMPFOL/region.mol2
	cd $TEMPFOL
	./script_demo
	
	cd ..
	cp $TEMPFOL/selected_spheres.sph $SPHFOL/"$fid""_sel.sph"
	echo "Wrote spheres to $SPHFOL/"$fid"_sel.sph"
	cd $TEMPFOL
	./script_clean
	rm *.pdb rec.ms region.mol2 OUTSPH -f
	cd ..
done

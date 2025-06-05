#!/bin/bash
#
# Dock stuff

LIGFILE=$1
if test -z "$LIGFILE"
then
	echo "Please pick a mol2 file for ligands"
	exit 1
fi

OUTFOL="dockout"
GRIDFOL="grids"
SPHFOL="spheres"
TEMPFOL="temp_dock"
PREFIX="S100B_b1_"
SUFFIX=".pdb"

mkdir -p $OUTFOL
echo $LIGFILE

for f in `ls *.pdb`
do
	fid=`echo $f|sed s/$PREFIX//g|sed s/$SUFFIX//g`
	echo $fid
	outf=$fid"_results.mol2"
	if test -f "$OUTFOL/$outf"
	then
		echo " Exists"
		continue
	else
		echo " Starting"
	fi

	ls $SPHFOL/$fid"_sel.sph"
        cp $SPHFOL/$fid"_sel.sph" $TEMPFOL/spheres.sph
        cp $GRIDFOL/"$fid""_grid.nrg" $TEMPFOL/grid.nrg
        cp $GRIDFOL/"$fid""_grid.bmp" $TEMPFOL/grid.bmp
	
	cd $TEMPFOL
	ln -s ../$LIGFILE ./todock.mol2
	../../../dock6/bin/dock6 -i anchor_and_grow.in
	cd ..
	cp $TEMPFOL/anc_results_ranked.mol2 $OUTFOL/$outf
	cd $TEMPFOL
	./script_clean
	cd ..
done

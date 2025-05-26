for lag in {5..200..5}
do
	mkdir  dt-$lag
	cd dt-$lag 
	cp ../train.slurm .
	cp ../config.ini .
	sed -i 's/Time/'$lag'/g' config.ini
	sbatch train.slurm
	cd ..
done

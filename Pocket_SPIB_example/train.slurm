#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --time=4:00:00
#SBATCH --job-name=spib

module purge
module load cuda
python /home/xg23/scratch/State-Predictive-Information-Bottleneck/test_model_advanced.py -config config.ini

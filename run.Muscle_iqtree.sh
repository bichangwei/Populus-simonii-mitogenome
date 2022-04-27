#!/bin/bash
#########################################################################
# FileName: run.Muscle_iqtree.sh
# Version: dc211e42-373c-43f3-9970-520d098bd1bc
# Author: Changwei <bichwei@163.cn>
# CreatedTime: Thu Sep  9 09:24:22 2021
#########################################################################
muscle -in Populus_simonii.phylogeny.fa -out Populus_simonii.muscle.alignment.afa -maxiters 2
iqtree -s Populus_simonii.muscle.alignment.afa -m MFP -B 1000 --bnni -T AUTO

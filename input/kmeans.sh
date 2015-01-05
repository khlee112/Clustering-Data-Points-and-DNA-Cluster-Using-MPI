#!/bin/sh
#Number of Points
b_point=200
#Number of Cluster
k_point=100
echo ********GENERATING $b INPUT POINTS EACH IN $k CLUSTERS
python ./randomclustergen/generatepoint.py -c $k_point  -p $b_point -o data/points.csv -v 10
#Number of Points
b_dna=400
#Number of Cluster
k_dna=50
echo ********GENERATING $b INPUT DNA STRAND EACH IN $k CLUSTERS
python ./randomclustergen/generatednastrand.py -c $k_dna  -p $b_dna -o data/DNA_strand.txt -v 100
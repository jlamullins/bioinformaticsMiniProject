#!/bin/bash
for i in `seq 1 70`;
	do
		echo ./gibbs ../bioData/dataset$i/sequences.fa ../bioData/dataset$i/motiflength.txt ../results/dataset$i/predictedsites.txt ../results/dataset$i/predictedmotif.txt $i
		./gibbs ../bioData/dataset$i/sequences.fa ../bioData/dataset$i/motiflength.txt ../results/dataset$i/predictedsites.txt ../results/dataset$i/predictedmotif.txt $i
	done

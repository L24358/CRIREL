#!/bin/bash

##################################
for d in ./*/ ; do
	
	cd $d
	echo Motif: $d
	
	for line in ./*/ ; do
		cd $line
		./flysim.out -pro network.pro -conf network.conf -t 4 -s moderate -nmodel LIF
		rm ConfsInfo.log
		rm network.log
		rm flysim.out
		cd ..

	done
	cd ..

done
##################################
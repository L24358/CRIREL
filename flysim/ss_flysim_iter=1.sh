#!/bin/bash

##################################
for d in ./*/ ; do
	
	cp flysim.out $d
	cd $d
	echo File: $d

	./flysim.out -pro network.pro -conf network.conf -t 4 -s moderate -nmodel LIF
	rm ConfsInfo.log
	rm network.log
	rm flysim.out
	cd ..

done
##################################

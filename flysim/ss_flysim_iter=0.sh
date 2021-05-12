#!/bin/bash

##################################
#run flysim
./flysim.out -pro network.pro -conf network.conf -t 4 -s moderate -nmodel LIF
rm ConfsInfo.log
rm network.log
##################################

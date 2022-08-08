#!/bin/bash

#allow to run parallel computation using openMP
#export OMP_NUM_THREADS=8
g++ KEbbReweight.c $(root-config --cflags --libs) -fopenmp -o KEbbReweight
./KEbbReweight

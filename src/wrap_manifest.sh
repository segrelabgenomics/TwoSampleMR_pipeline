#!/bin/bash

# Authors: Skanda Rajasundaram, Puja Mehta
# Segre lab, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
# Date: November 16, 2023

#############################
### Default UGER Requests ###
#############################

# This section specifies uger requests.
# This is good for jobs you need to run multiple times so you don't forget what it needs.

#$ -cwd
#$ -N MR

# Memory request for 15G
#$ -l h_vmem=125G

######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts

# Use your dotkit
reuse UGER
reuse R-4.1

##################
### Run script ###
##################

Rscript --vanilla  MR/src/MR_Script.R $1 $2 $3 $4 $5 $6 $7 $8 $9

#Arguments:
#Trait name, 
#Gene ID, 
#Gene Symbol, 
#Tissue, 
#QTL type, 
#p-Value cutoff, 
#File name of the trait, 
#path to the QTL file and 
#QTL file extension

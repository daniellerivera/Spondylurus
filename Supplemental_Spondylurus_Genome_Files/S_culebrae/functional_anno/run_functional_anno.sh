#!/bin/bash

source /home/drivera/Spondylurus_culebrae/anno/functional_anno/functional_anno_routines.sh

#############################################################################################################################
################################  Created by functional_anno_setup.sh  ######################################################
#############################################################################################################################

# vertebrate or arthropoda vars can have any value to trigger orthodb runs

vertebrate=yes
#arthropoda=yes

# use the more-sensitive diamond blastp for orthodb. comment this out to use the much slower blastp
diamond_orthodb=yes

threads=48

run_anno


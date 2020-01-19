#!/bin/bash
# Transfers behavioral and functional and structural data to mac for local analysis
#  usage: transferRawDataForFmri.sh PPID
#
# PPID as numeral, no project id
#

# sync behavioral and phys data
rsync -tuvaz colombo09:DataSPIFI/SPIFI_000$1/logs \
../Results/RawFmri/SPIFI_000$1

# sync anatomical
rsync -tuvaz colombo09:DataSPIFI/SPIFI_000$1/scandata/* \
../Results/RawFmri/SPIFI_000$1/anat

# sync functional data, is reconstructed
rsync -tuvaz --exclude="**phase_*" colombo09:ResultsSPIFI/RawFmri/SPIFI_000$1/func \
../Results/RawFmri/SPIFI_000$1


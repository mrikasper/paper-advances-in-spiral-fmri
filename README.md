Spiral Fmri 7T Paper
====================

Summary
-------

This repository contains the analysis code and [manuscript](Manuscript.md) 
(including representation code for all figures) for the paper 
"[Advances in Spiral fMRI - a High-resolution Study with Single-shot Acquisition](https://www.biorxiv.org/content/10.1101/842179v1)".

This work is part of the Spiral Functional Imaging (SPIFI) project at the
Institute for Biomedical Engineering conducted at the ETH Zurich and University
of Zurich.

This paper gives the first account of spiral fMRI at ultra-high field (7 Tesla)
and with sub-millimeter (0.8 mm) spatial resolution. This progress was enabled
by expanded signal modeling of the MRI signal, including static and dynamic
magnetic fields, their accurate concurrent measurement utilizing NMR field
probes, and the corresponding iterative model inversion using conjugate-gradient
SENSE.

The employed task is a simple visual quarter-field stimulation paradigm.

Installation
------------

The preferred installation method is to clone this repository with the
*recursive* option, in order to also checkout the required third-party tools as
submodules. On Unix/Mac command line, this

`git clone --recursive https://tnrurepository.ethz.ch/spiral_fmri_7t_paper.git`

If you download the `.zip` file instead, please make sure to also download the
packages mentioned in the Requirements section below, and place them in the
`Toolboxes` subfolder of this repository.

Getting Started
---------------

To reproduce the fMRI analysis, run `Code/main.m` in Matlab. This should perform

1.  The preprocessing, including slice timing correction, realignment and
    smoothing

2.  The first level GLM analysis of the visual quarterfield task, including
    physiological noise modeling (RETROICOR).

3.  The generation of all figure content for the accompanying 
    [manuscript](Manuscript.md). This can also be run separately in
    `Code/Representation/main_create_figures.m`.

Timeline
--------

-   Project started in September 2016

-   First successful imaging results November 2016

-   Started Paper in October 2017

    -   General structure repo: October 2017

    -   Introduction: February 2018

-   Paper Data

    -   Started Acquisition in October 2017

    -   Started Analysis (maps) in December 2017

    -   Started pipeline recon on Euler in February 2018

    -   Started Analysis Code in February 2018

    -   Improved B0 map processing for LAYMM and SPIFI: October - December 2018

    -   Rerun Analysis of all subjects: March 2019

        -   with LAYMM map improvements

    -   2nd Rerun Analysis of all subjects: June 2019

        -   with LAYMM recon phase 2 map improvements

-   First Full Story Bullets: September 06, 2019

-   First Complete Figure Drafts: October 04, 2019

-   First Full Draft: November 10, 2019

-   First Preprint ([BiorXiv](https://www.biorxiv.org/content/10.1101/842179v1)): November 15, 2019 

-   First Submission (Neuroimage): January 18, 2020

-   Starting Revision:

-   Submitted Revision:

-   Accepted:

-   Proofs Accepted:

First Author: Lars Kasper

Requirements
------------

This code is written in Matlab, R2018b. It further relies on the following
open-source Matlab toolboxes

-   Statistical Parametric Mapping [SPM12.4](https://github.com/spm-central/spm12)

-   Unified NeuroImaging Quality Control Toolbox [UniQC](https://gitlab.ethz.ch/uniqc/uniqc-code)

-   PhysIO Toolbox for Physiological Noise Modeling, part of [TAPAS](https://translationalneuromodeling.github.io/tapas)

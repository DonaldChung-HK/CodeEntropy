---
title: 'CodeEntropy: A Python package for Multiscale Entropy and Structure Quantification from Molecular Dynamics Simulation'
tags:
  - Python
  - entropy
  - molecular dynamics
  - molecular simulations
authors:
  - name: Arghya Chakravorty
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Jas Kalayan
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Donald Chung
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 3
affiliations:
 - name: University of Michigan, Ann Arbor, USA
   index: 1
 - name: University of Manchester, United Kingdom
   index: 2
 - name: Science and Technology Facilities Council, United Kingdom 
   index: 3
date: 21 July 2022
bibliography: paper.bib
---
# Summary

Entropy is a fundamental property of any system that quantifies the structural flexibility of the system and together with energy governs system stability. It is particularly important in biomolecular systems because of their highly flexible and complex nature. Many properties are a delicate balance of entropy versus energy, necessitating the determination of entropy in order to understand stability. Moreover, entropy provides a way to quantify the structural flexibility of a system over all its degrees of freedom. ``CodeEntropy`` is a code based on the Multiscale Cell Correlation (MCC) method which is a novel solution to the problems encountered by other methods by providing a single, scalable and general framework applicable to all molecules in the system.

``CodeEntropy`` is a code based on a combination of ``CodeEntropy`` [@argoRepo] by Dr Arghya Chakravorty and ``POSEIDON`` [@jasRepo] by Jas Kalayan. The code written by Arghya Chakravorty  accounts for the vibrational and conformational entropy in a multiscale formulation using an Applications Programming Interface that makes it highly customisable. To make ``CodeEntropy`` fully applicable to biomolecular systems, the ``POSEIDON`` code will calculate the topographical entropy terms for solvents and mixtures. The topology/trajectory parser, atom selector and distance calculation is performed by the ``MDAnalysis`` [@mda1; @mda2] package.

# Statement of needs
There are a range of existing methods to calculate entropy from molecular dynamics simulations but they suffer from a number of limitations: they may only work for particular kinds of system or degrees of freedom, they may require additional calculations, they can be difficult to interpret, or they do not scale well to large and complex systems. Some methods only work for water, such as Inhomogeneous Solvation Theory, for liquids such as 2-Phase Thermodynamics, for only some degrees of freedom such as dihedral binning, or for single molecules such as Quasiharmonic Analysis, Normal Mode Analysis or non-parametric methods such as Minimal Spanning Tree or K-Nearest-Neighbours.

Given the widespread use of free-energy calculations and molecular dynamics simulations, there is a large user-community for software to calculate entropy and quantify full structural flexibility of biomolecular systems. Multiscale Cell Correlation (MCC) provides a novel solution to the problems encountered by other methods by providing a single, scalable and general framework applicable to all molecules in the system. It utilises a judicial synthesis of mean-field cell theory and covariance matrices over a range of length scales: 
- Correlations are considered between groups of locally connected atoms as in a mean-field cell, and longer-range correlations are accounted for using a coarser representation of the groups, a framework that is scaled to higher length scales. 
- At each length scale, the potential energy surface is discretised into energy wells for translational and rotational motion. These are represented as an average energy well and an energy-well distribution, denoted as vibrational and topographical, respectively.
- The decomposition over molecules, length scales, type of motion and energy-well size and distribution provides an exquisite level of detail in explaining the entropies obtained

MCC has been applied by the group of RH to a wide range of systems, namely liquids, aqueous and octanol solutions, host-guest complexes, chemical reactions and large biomolecules such as proteins, DNA and membrane bilayers in aqueous electrolytes. 

# References
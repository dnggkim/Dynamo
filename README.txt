DYNAMO - Quantum Dynamic Optimization Package

Version 1.4.0 alpha1
Released 2014-??-??


Introduction
============

DYNAMO is a flexible framework for quantum optimal control algorithms.
It requires MATLAB R2011a or newer (matlab.mixin.Copyable).

The latest development version can be downloaded from our Git repository at https://github.com/smite/Dynamo
For the latest version of this software, guides and information, visit http://www.qlib.info

The user manual can be found here: docs/dynamo_manual.tex

If you use DYNAMO in your research, please add an attribution in the form of the following reference:
S. Machnes et-al, arXiv 1011.4874


Getting started
===============

Initialize the package by running init.m

The best way to learn how to use DYNAMO is to review the demos in the
examples/ directory.


The design is modularized and easily extendable.

DYNAMO attempts to minimize calculations by way of Delayed Calculations.
Whenever a control field is modified, DYNAMO only marks that value as changed,
and the values calculated from it (single-slice exponents, start-to-T and 
T-to-end propagators, etc) as stale. Only when a specific value is needed
(such as the current fidelity), are the calculations performed. At this point
DYNAMO attempts to perform the minimal number of matrix exponentiations and 
multiplications to arrive at the desired result.

See cache.m for the nitty-gritty details.


License
=======

Released under the terms of the Lesser GNU Public License and
Creative-Commons Attribution Share-Alike (see "LICENSE.txt" for details).


Authors
=======

Shai Machnes         2010-2012
Ville Bergholm       2011-2014

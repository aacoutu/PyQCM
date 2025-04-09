Included in the "scripts" folder is a set of Python code for running an acoustic transfer-matrix analysis of films measured using a quartz crystal microbalance (QCM).  This analysis (in "Modeling_Film_Properties.py") makes use of a continuum mechanics model of shear wave propagation, where the film deposited on the QCM sensor can be modeled using an arbitrary number of layers in order to identify local gradients in viscoelastic properties or density.  See our study using this approach: "Analyzing QCM Data Using a New Transfer-Matrix Model: Long-Ranged Asymmetric Gradient in Shear Modulus Identified Across Immiscible Glassy-Rubbery Polymer Interface," Alexander A. Couturier, Justin C. Burton, Connie B. Roth, Macromolecules 2025, 58, 3520-3536 (https://doi.org/10.1021/acs.macromol.4c02847).  In this study, we use the code included here to determine a large gradient in shear modulus across a glassy/rubbery polymer interface of polystyrene (PS) and polybutadiene (PB).  Detailed information on the mathematical reasoning underlying the code can be found in our publication. 

The input data one needs to provide "Modeling_Film_Properties.py" are the shifts in resonance frequency and dissipation of the QCM sensor as a result of the deposition of the studied film.  Such data can be extracted from raw resonance traces using "Circuit_Analysis.py" (also included in the "scripts" folder), where the specific circuit equation used corresponds to our experimental setup.  This code analyzes a given set of QCM resonance traces at different harmonics and outputs the associated resonance frequency and dissipation data.  For the details of our circuit equation and how it corresponds to our experimental setup, see the publication: "Physically Intuitive Continuum Mechanics Model for QCM:  Viscoelasticity of Rubbery Polymers at MHz Frequencies,"  Yannic J. Gagnon, Justin C. Burton, and Connie B. Roth, Journal of Polymer Science 2022, 60, 244-257 (https://doi.org/10.1002/pol.20210763).  

-------

Make sure that all imported csv files are the basic comma delimited file format, and not another csv format.

lmfit library installation in command prompt: if using pip, run- "pip install lmfit" ; if using conda, run one of the lines found here- https://anaconda.org/conda-forge/lmfit

This should be the only library installation needed. In case of issues with version dependence, or if a virtual environment is desired, below are instructions for a stable installation of libraries to a virtual environment.

Install Anaconda and run these lines in order from the Anaconda command prompt:


conda create --name QCM_analysis_env  
conda activate QCM_analysis_env

conda install python=3.8  
conda install spyder=5.5.1  
conda install numpy=1.24.3  
conda install matplotlib=3.7.2  
conda install pandas=2.0.3  
conda install tk=8.6.14  
conda install sympy=1.12  
conda install conda-forge::lmfit=1.2.2

<br>
Then run Spyder from the Anaconda command prompt while in the virtual environment.

-------

lmfit copyright information:

"BSD-3

Copyright 2022 Matthew Newville, The University of Chicago
               Renee Otten, Brandeis University
               Till Stensitzki, Freie Universitat Berlin
               A. R. J. Nelson, Australian Nuclear Science and Technology Organisation
               Antonino Ingargiola, University of California, Los Angeles
               Daniel B. Allen, Johns Hopkins University
               Michal Rawlik, Eidgenossische Technische Hochschule, Zurich

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  3. Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.


Some code has been taken from the scipy library whose licence is below.

Copyright (c) 2001, 2002 Enthought, Inc.
All rights reserved.

Copyright (c) 2003-2019 SciPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of Enthought nor the names of the SciPy Developers
     may be used to endorse or promote products derived from this software
     without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

Some code has been taken from the AMPGO library of Andrea Gavana, which was
released under a MIT license."

Time-Harmonic-Elastography
=========================

Overview
--------
***Time-Harmonic-Elastography*** is a Matlab toolbox for time-harmonic optical flow.
It implements reconstruction algorithms and simulation of synthetic data for time-harmonic optical flow.

__Reconstruction Methods__
* Least squares solver for L_2^2 data fidelity and L_2^2 regularizer
* Iteratively reweighted least squares algorithm for L_1 data fidelity and L_1 regularizer
* Iteratively reweighted least squares algorithm for L_1 data fidelity and L_2^2 regularizer 
* Adoption of Matlab Horn-Schunk implementation and primal-dual hybrid gradient method for time-harmonic optical flow

References
----------
When you are using this code, please cite the paper

[1] Oleh Melnyk, Michael Quellmalz, Gabriele Steidl, Noah Jaitner, Jakob Jordan, Ingolf Sack:
__Time-Harmonic Optical Flow with Applications in Elastography__.
*to appear*

This paper also explains the algorithms in more detail.

The code for primal-dual hybrid gradient method is described in the paper

[2] Frank Balle, Tilmann Beck, Dietmar Eifler, Jan Henrik Fitschen, Sebastian Schuff, Gabriele Steidl:
__Strain analysis by a total generalized variation regularized optical flow model__.
*Inverse Problems in Science and Engineering* 27.4 (2019), pp. 540–564. [doi: 10.1080/17415977.2018.1475479](https://doi.org/10.1080/17415977.2018.1475479). 

Dependencies
------------
This software uses [export_fig](https://de.mathworks.com/matlabcentral/fileexchange/23629-export_fig) for all source files in folder "figures" 

Directory structure
-------------------

File/Folder                | Purpose
--------------------------:| ---------------------------------------------------------------------
OpticalFlow (dir)          | Contains implementation of algortihms for time-harmonic optical flow 
Util (dir)                 | Contains utility functions such as computation of discrete derivatives and simulation of the displaced images
figures (dir) 	           | Source code for the numerical experiments in [1] and gel dataset
flow_programs_matlab (dir) | Implementation of primal-dual hybrid gradient method from [2]
COPYING                    | License information
README.md                  | This file

Feedback
--------
Your comments are welcome! This is the first version of the library and may
not be as robust or well documented as it should be. Please keep track of bugs
or missing/confusing instructions in Guthub Issues.
Alternatively, you might contact
[Oleh Melnyk](mailto:oleh.melnyk@tu-berlin.de)
or
[Michael Quellmalz](mailto:quellmalz@math.tu-berlin.de).

Legal Information & Credits
---------------------------
Copyright (c) 2024 Oleh Melnyk and Michael Quellmalz

This software was written by [Oleh Melnyk](https://olehmelnyk.xyz/) and [Michael Quellmalz](https://page.math.tu-berlin.de/~quellm/).
It was developed at the Institute of Mathematics, TU Berlin.

This research was funded in whole or in part by the German Research Foundation (DFG): STE 571/19-1, project number 495365311, within the Austrian Science Fund (FWF) [SFB Tomography Across the Scales](https://tomography.univie.ac.at/). 
The authors acknowledge funding from the German Research Foundation (DFG) within the project BIOQIC (GRK2260).

Time-Harmonic-Elastography is free software. You can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. If not stated otherwise, this applies to all files contained in this
package and its sub-directories.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

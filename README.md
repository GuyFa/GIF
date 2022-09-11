# Globally Injective Flattening via a Reduced Harmonic Subspace
--------------------------------------------------------------------------------------------------------------------------------------------------------

This code includes the implementation of the SIGGRAPH ASIA 2022 paper "Globally Injective Flattening via a Reduced Harmonic Subspace" authored by Guy Fargion and Ofir Weber. 

The use of this application is limited to academic use only!

The code is provided as-is and without any guarantees.

This code is partially based on https://github.com/eden-fed/HGP#harmonic-global-parametrizations.

For questions or comments about the code please contact:
Guy Fargion (guy.fargion@gmail.com)

----------------------------------------------------------------------------
A Visual Studio 2017 (MSVC 14.1) project is provided for easy compilation on Windows machines.
The code should be platform independent though we never tested it on other than Windows OS.

The following prerequisites are necessary for building and running the code:

1) Matlab R2019b

2) Boost 1.80.0 - I downloaded boost from here https://www.boost.org/. My code only requires the headers of the Boost libraries. Hence, there is no need to build boost as well.

3) CGAL 5.5 - CGAL has to be compiled with 64 bits dynamic (also called shared) - in order to do that we used cmake.

4) GMM C++ template library version 4.2 (http://getfem.org/download.html).

5) PARDISO 6.0

6) Eigen 3.7.7


Other versions of the above listed tools might be compatible but weren't tested!

----------------------------------------------------------------------------

How to use:

1) Install and build the above mentioned prerequisites.

2) Add the following environment variables to your system (see some possible paths):

GMM_INCLUDE_DIR		  (%your GMM folder path%)\gmm-4.2\include

MATLAB_64_DIR		    C:\Program Files\MATLAB\R2019b

CGAL_64_DIR		      C:\Program Files\CGAL-5.5

BOOST_64_DIR		    C:\Program Files\boost\boost_1_80_0

EIGEN_DIR  		      (%your Eigen folder path%)

PARDISO_BIN         (%path to the folder where PARDISO dll and lib are%)

PARDISO_LIC_PATH    (%path to the folder with PARDISO licence%)

OMP_NUM_THREADS  	  number of cores in your CPU (for PARDISO)

3) Add the the folder "MatlabScripts" to your Matlab path.

4) Make sure all the required dlls can be loacted by including the relevant paths into the system PATH variable.
For example:
PATH=
%MATLAB_64_DIR%\bin\win64;
%MATLAB_64_DIR%\extern\include;
%MATLAB_64_DIR%\extern\lib\win64\microsoft;
%BOOST_64_DIR%\libs;
%CGAL_64_DIR%\lib;
%CGAL_64_DIR%\auxiliary\gmp\lib;
%CGAL_64_DIR%\include\CGAL;
%PARDISO_BIN%;
%MATLAB_64_DIR%

5) Execute GIF.exe.
This should invoke Matlab and open a dialog box.
You will be asked to choose your .obj file and your settings (the default settings for the isometric and conformal options are already given there).
	

# 2D-Thermal-heat-transfer-in-electrically-heated-pavement

Provided here is a Fortran90 code to calculate temperature levels within a 2D domain representing an electrically-heated pavement.

The Fortran90 code takes as input the files: (i) "Airtemp.dat"; (ii) "Windspeed.dat"; and (iii) "Longwave.dat". These correspond to measured and calculated weather effects. The file "para.dat" is also required.

The Fortran90 code produces as output Tdepth.XX.dat files, which can be treated with the Python code "Data_analysis.py". This code takes as input all TdepthXX.dat files and generates a "Summary.dat" file. "Summary.dat" contains the maximal and minimal temperature levels across the width of the pavement.
To launch the Fortran90 code, place all the input files and "temp.f90" within the same folder. Open a terminal, navigate to the folder, and type "gfortran -o exe temp_updated_structure.f90", wait for the compilation (less than a second) and then type "exe".

To run the Python code, make sure you change appropriately the path folder.

This code is linked to a manuscript which will be published. Once published, I will share the citation information if you decide to utilize the code.

Fortran90 compiler : gfortran (GNU Fortran GCC version 12.2.0) obtained on http://www.equation.com/servlet/equation.cmd?fa=fortran Python version 3.11.0

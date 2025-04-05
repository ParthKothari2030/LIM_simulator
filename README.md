# LIM_simulator
This tool simulates galaxy-line intensity maps from a given halo catalog.

## Prerequisites
"g++" compiler is required to compile the code. Any other compiler can be used, requiring the makefile to be altered accordingly.

## Format of halo catalogs
The halo catalogs need to be in the following format:
There should be a directory inside which seperate files hold data for the halo masses and halo positions. For example, "Halos_z7.0" can be directory. Inside this there need to be four binary files, "Pos_x.bin", "Pos_y.bin", "Pos_z.bin", and "Mass.bin", with no header. The positions are those of the halos in grid units. Mass can be in any unit. Additionally, there should be a file named "num_halos.txt", which should contain a single integer number which is the number of halos in the files.

## Intensity mapping models
Inside the "src" directory, the source file "sfr_buffered.cc" contains the default model of star-formation rate for a given halo mass (Eq. (26) of [Silva et al. 2013](https://doi.org/10.1088/0004-637X/763/2/132)). A different model of star-formation rate can be included and used in the source file, "process_buffered.cc". The source file, "lum_buffered.cc", converts SFR to line-luminosities using the pre-defined functions. New models can be added and used in "process_buffered.cc".

**Revised code** -  Added the HI model prescribing HI to individual halo based on its mass (Scheme 3 of [Bagla et al. 2010](https://academic.oup.com/mnras/article/407/1/567/986082), Eq. (6). Also check [Sarkar et al. 2016](https://academic.oup.com/mnras/article/460/4/4310/2609038) Eq. (5) )

## How to use
Inside the directory "LIM_simulator", type "make" and hit "Enter", which will compile the code and generate an executable "LIM_simulator". This executable needs two text files as inputs to function. The first one (e.g. "param.txt") contains the necessary paramters for the simulation. An example file named "param.txt" is given here, with the meaning of the parameters explained inside. The second file (e.g. "paths_CO.txt") should contain the paths for the halo catalogs, redshifts and names of the output files aranged in seperate columns. For example, this might look like:
```
Halos_z7.000 7.000 CO_map_z7.000.bin
Halos_z6.000 6.000 CO_map_z6.000.bin
Halos_z5.000 5.000 CO_map_z5.000.bin
```
Here, the first column represents the directory that holds the halo catalogs. The second column represents the redshift and the third one represents the output binary files that contain the intensity mapping cube. To run the executable, use command "./LIM_simulator param.txt paths_CO.txt". 


**Revised executable** - The code now externally takes the parameters of CO and HI. To pass the parameter for CO and HI, look at the respective executables "./LIM_simulator param.txt paths_CO.txt <alpha_param> <beta_param>" or "./LIM_simulator param.txt paths_HI.txt <Omega_HI_param>".


Note that the path file name is changed for clarification purposes.


The format of the output files are such that it can it can be read with the following python script:
```
import numpy as np
f = open ('CO_map_z6.000.bin', 'rb')
N = np.fromfile (f, count=3, dtype=np.uint64)
N1,N2,N3 = N
grid_size = np.fromfile (f, count=1, dtype=np.float32)
data = np.fromfile (f, count=N1*N2*N3, dtype=np.float32)
f.close ()
data = data.reshape ((N1,N2,N3), order='C')
```
"N1","N2", and "N3" are the number of grids along each axis of the intensity mapping cube, while "grid_size" is the size of the smallest cell in Mpc.

## Acknowledgements
If you are using this tool, please acknowledge the following reference:[https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2500M/abstract](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2500M/abstract)

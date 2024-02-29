1. Run scf calculatin with ORCA (get water.gbw file)
2. Run orca with the option keepdens xyzfile
3. python mep.py water 40 32  # where 40 is the grid size and 32 is the number of cores per node
4. This will generate water_mep.cube
5. Run orca_plot water.gbw -i
   1 - Enter type of plot
   2 -   (scf) electron density      ......  (scfp  )  - available
   The default name of the density would be: water.scfp
   Is this the one you want (y/n)? y
   4 - Enter number of grid intervals
   Enter a number: 4
   Enter NGRID: 120
   5 - Select output file format
        7 -  3D   Gaussian cube
   10 - Generate the plot
   11 - exit this program
6. This will create water.eldens.cube
7. Now plot the density and potential cube with jmol
8. Have a nice day!  

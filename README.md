1. Run scf calculatin with ORCA (get water.gbw file)
2. Run orca with the option keepdens xyzfile
3. python mep.py water 40 32  # where 40 is the grid size and 32 is the number of cores per node
4. This will generate water_mep.cube
5. Run orca_plot water.gbw -i
6. 1 - Enter type of plot
7. 2 -   (scf) electron density      ......  (scfp  )  - available
8. The default name of the density would be: water.scfp
9. Is this the one you want (y/n)? y
10. 4 - Enter number of grid intervals
11. Enter a number: 4
12. Enter NGRID: 120
13. 5 - Select output file format
14. 7 -  3D   Gaussian cube
15. 10 - Generate the plot
16. 11 - exit this program
17. This will create water.eldens.cube
18. Now plot the density and potential cube with jmol
19. Have a nice day!  

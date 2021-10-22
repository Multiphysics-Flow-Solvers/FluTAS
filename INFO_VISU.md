### Output files
The code performs output through the library 2DECOMP, which dumps one binary file for each field. The results can be then visualized through the files available in the folder `src/data/`. To do so:

  1. copy in the folder `src/data/` the following files: `param.h90`, `genview.sh` and `gen_xdmf.f90`.
     If you don't have them, you can take those files from one of the examples. E.g.:
         cp ../../examples/abc_triperiodic/visu_ex/* .

  2. change `param.h90` consistently with the input files (the .in files) and the number of fields you want to visualize;
 
  3. on the command line, type: `./genview.sh`
     This command generates a file (viewfld.xdmf) in the same directory;

  4. To open it, type on the command line `paraview viewfld.xdmf`.

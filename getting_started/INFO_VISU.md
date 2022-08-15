# how to visualize the output binary files

### the easy way

In addition to the binary files for visualization, the code now generates a log file that contains information about the saved data (see `out2d.h90` and `out3d.h90` for more details); this new approach uses that log file to generate the `Xdmf` visualization file.

The steps are as follows:

1. after the simulation has run, copy the contents of `utils/visualize_fields/gen_xdmf_easy/write_xdmf.py` to the simulation `data` folder;
2. run the file with `python write_xdmf.py` in the `data` folder. If successful, this operation generates `viewfld_DNS.xmf` or `viewfld_DNS_2d.xmf` (see below) in the same folder;
3. load the generated Xdmf (`*.xmf`) file using paraview/visit or other visualization software, i.e., on the command line `paraview viewfld_DNS.xmf`. If requested, choose XDMF Reader in Paraview!

## example: how to visualize the default binary output

### 3D fields

when running the script, `write_xdmf.py` we get the following prompts:

~~~
 $ python write_xdmf.py
 Name of the log file written by FluTAS [log_visu_3d.out]:
 Name to be appended to the grid files to prevent overwriting []:
 Name of the output file [viewfld_DNS.xmf]:
~~~

* the first value is the name of the file that logged the saved data;
* the second is a name to append to the grid files that are generated, which should change for different log files to prevent conflicts;
* the third is the name of the visualization file.

By pressing <kbd>enter</kbd> three times, the default values in the square brackets are assumed by the script; these correspond to the default steps required for visualizing 3D field data.

### 2D fields

the procedure for visualizing 2D field data that is saved in `out2d.h90` is exactly the same; it is just that the correct log file should be selected. The code saves by default field data in a plane of constant `y=ly/2`, and logs the saves to a file named `log_visu_2d_slice_1.out`. If more planes are saved, the user should make sure that one log file per plane is saved by the code (e.g. if another plane is saved, the log file written in `out2d.h90` could be named `log_visu_2d_slice_2.out`); see `out2d.h90` for more details. The corresponding steps to generate the Xdmf file would be, for instance:

~~~
 $ python write_xdmf.py
 Name of the log file written by FluTAS [log_visu_3d.out]: log_visu_2d_slice_1.out
 Name to be appended to the grid files to prevent overwriting []: 2d
 Name of the output file [viewfld_DNS.xmf]: viewfld_DNS_2d.xmf
~~~

## Alternative Approach
The code performs output through the library 2DECOMP, which dumps one binary file for each field. The results can be then visualized through the files available in the folder `src/data/`. To do so:

1. Copy in the folder `src/data/` the following files: `param.h90`, `genview.sh` and `gen_xdmf.f90`. If you do not have them, you can take those files from one of the examples, e.g.: cp ../../examples/two_phase_inc_isot/rising_bubble_3d/visu_ex/* .
2. Change `param.h90` consistently with the input files (the .in files) and the number of fields you want to visualize;
3. On the command line, type: `./genview.sh`. This command generates a file (viewfld.xdmf) in the same directory;
4. To open it, type on the command line `paraview viewfld.xdmf`.

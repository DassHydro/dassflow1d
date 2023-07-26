###########################################################
# HOW TO RUN THIS EXAMPLE
###########################################################

1) Compile the executable
  - setup the compilation parameters in Makefile.inc, at least:
    CASEDIR=cases/forward/dam_break
  - Get to the DassFlow_SQ directory and run the following command:
    $ make
    
2) Run the executable in forward mode
  - run the following command:
    $ make run
  
3) Create the results plots
  - run the following command:
    $ python scripts/plot_results.py -dir cases/forward/dam_break
  
Remark: In order to run again the case from start, run the following commands:
    $ make cleanres
    $ make runexe


###########################################################
# HOW TO CREATE MORE MESHES
###########################################################
1) Generate a new mesh
  - modify parameter N in file create_mesh.py
  - run the script create_mesh.py
  
2) Modify the input.txt file
  - replace mesh_name "mesh.geo" with the name of the created file (channel_[N].geo)

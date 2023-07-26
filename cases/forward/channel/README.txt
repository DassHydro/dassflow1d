###########################################################
# HOW TO RUN THIS EXAMPLE
###########################################################

1) Compile the executable
  - setup the compilation parameters in Makefile.inc, at least:
    CASEDIR=cases/forward/channel
  - Get to the DassFlow_SQ directory and runthe following command:
    $ make
    
2) Run the executable in forward mode
  - run the following command:
    $ make run
  
3) Create the results plots
  - run the following command:
    $ python scripts/plot_results.py -dir cases/forward/channel
  
Remark: In order to run again the case from start, run the following commands:
    $ make cleanres
    $ make runexe

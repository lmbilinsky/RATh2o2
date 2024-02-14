This script contains the in vitro computational model of endogenous H2O2 metabolism in rat hepatocytes (rat hepatocytes
existing in an in vitro preparation). It is identical to the script RATh2o2_dma3.m, minus everything having to do with
arsenic. 

To run this script, all files in this repository must be in the same directory (folder). The script will load the .mat file
yinitialcondition.mat, which contains the variable yinit. yinit dictates the initial concentrations of all chemical species
in both hepatocytes and medium, and the initial number of cells (element 65, set here to 20,000). 

RATh2o2 describes the situation of 20,000 rat hepatocytes existing per 100 microliters of medium. If you wish to model a 
different in vitro setup, then the constants cellnum and/or vmed in RATh2o2.m ,and possibly element 65 of yinit, must be 
changed. 

IMPORTANT: the initial concentrations of cystine and cysteine in medium (elements 10 and 11 of yinit, respectively) are 
set to values intended to describe the composition of Dulbecco's Modified Eagle medium. Only the total cysteine content
is specified on its list of ingredients, and I had to estimate how this is divvied up between cysteine and cystine. Please
see the Supplementary Material of my paper for how this was done. It may be that these values are not accurate, but this
isn't important for model behavior. HOWEVER: if you wish to use RATh2o2.m to simulate a different kind of medium, it is very
important that you look at its list of ingredients and alter elements 10 and 11 of yinit so that 2*[cystine]+[cysteine] equals
the total listed cysteine content.

The initial concentrations of chemicals species in hepatocytes stored in yinit reflect the normal (in vivo) values reported in the
physiology literature. The model is constructed to yield a steady-state solution which agrees with these values. HOWEVER: solutions
will only approach this steady state if equation 11 of RATh2o2.m is commented out (effectively holding medium cysteine constant). 
If this is done, the file SteadyState.txt, generated by RATh2o2.m, will print out the steady-state concentrations and also, the steady
state reaction and transport velocities. If equation 11 is not commented out, solutions will not approach a steady state. Instead,
GSH in hepatocytes will decline with time, H2O2 in hepatocytes will rise, cysteine/cystine in medium will decline, and GSH in medium 
will rise. This is because cellular medium lacks the enzymes necessary to break GSH (exported from hepatocytes) down to its component
amino acids, one of which is cysteine. Given how large the volume of medium is compared to the total volume of the 20,000 
hepatocytes, this effect takes a while to see.

To simulate the effects of a total halt in GSH synthesis, set ``V" in VGCLholorat.m to zero. 


 
Dulbecco's Modified Eagle Medium. It is identical to 

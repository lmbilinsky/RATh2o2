The script RATh2o2.m contains the in vitro computational model of endogenous H2O2 metabolism in rat hepatocytes (rat hepatocytes
existing in an in vitro preparation). It is identical to the script RATh2o2_dma3.m, minus everything having to do with
arsenic. 

To run this script, all files in this repository must be in the same directory (folder). The script will load the .mat file
yinitialcondition.mat, which contains the variable yinit. yinit dictates the initial concentrations of all chemical species
in both hepatocytes and medium, and the initial number of cells (element 65, set here to 20,000). 

RATh2o2.m assumes that there are 20,000 rat hepatocytes for every 100 microliters of medium, and that the medium formula is
Dulbecco's Modified Eagle Medium (DMEM). If you wish to model a different in vitro setup, then the constants cellnum (and element 
65 of yinit) and/or vmed in RATh2o2.m  must be changed. Regarding the initial concentrations of cystine and cysteine 
in medium (elements 10 and 11 of yinit, respectively), they are set to values intended to describe the composition of DMEM. 
Only the total cysteine/cystine content is specified on its list of ingredients, which I deduced was equivalent to 260.5 micromolar
cystine equivalents (0.5 [cysteine]+[cystine]). I had to estimate how this is divvied up between cysteine and cystine. I assumed that 
the redox balance bewteen the two forms is such that cysteine is present in medium at its normal concentration in rat blood plasma, about
20 micromolar. If you wish to use RATh2o2.m to simulate a different kind of medium, it is very important that you look at its list of 
ingredients and alter elements 10 and 11 of yinit so that 0.5 [cysteine]+[cystine] equals the total listed cystine content.

The concentrations of biochemical species in rat hepatocytes stored in yinit reflect the normal (in vivo) values reported in the
physiology literature. The model is constructed to yield a steady-state solution consistent with these values. HOWEVER: solutions
will only approach this steady state if equation 11 of RATh2o2.m is commented out (effectively holding medium cysteine constant). 
If this is done, the file SteadyState.txt, generated by RATh2o2.m, will print out the steady-state concentrations and also, the steady
state reaction and transport velocities. If equation 11 is not commented out, solutions will not approach a steady state. Instead,
GSH in hepatocytes will decline with time and GSH in medium will rise. This is because cellular medium lacks the enzymes necessary to 
break down exported GSH and regenerate its component amino acids, recycling cysteine. Given how large the volume of medium is compared to 
the total volume of the 20,000 hepatocytes, this effect takes a while to become significant. After 72 hours, intracellular GSH is 
significantly depleted. The larger is the volume of medium compared to the volume (equivalently, number) of hepatocytes, the longer it 
takes for this effect to become significant. To use an analogy, if I am fishing in an aquarium, I will quickly deplete the number of fish; 
if I am fishing in a lake, this effect is trivial on relevant time scales.

To simulate the effects of a total halt in GSH synthesis, set ``V" in VGCLholorat.m to zero. This generates Fig. 6 in the
paper. Note that there is a typo in the text of section 3.3: where it reads ``N(24)=20,000" it should instead read ``N(0)=20,000."


 


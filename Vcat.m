%catalase removal of H2O2

function a = Vcat(h2o2);



Vmax=1e6;



km = 32300;   %for mouse liver 


a = Vmax*h2o2./(km + h2o2);  %Chen 2015, for mouse liver




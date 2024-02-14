%cysteine dioxygenase, the major route by which cysteine is converted to taurine (Weinstein88)
%according to Stipanuk92, at low levels of cysteine 90% is used for glutathione synthesis. At high levels of
%cysteine, 20% is used for GSH synthesis. According to Weinstein88 the catabolism of cysteine at high levels
%is about half to taurine and about half to pyruvate.

function output = VCYSDIOXYGENASEm(cys);


a=0.006;
b=(194-104^2*a)/104;



output=a*cys.*cys+b*cys;




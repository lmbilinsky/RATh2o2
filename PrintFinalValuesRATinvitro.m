function PrintFinalValues(vcell,vmed,M);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRINT OUT CONCENTRATIONS AND VELOCITIES AT FINAL TIME
t=datestr(now); % today's date
fid=fopen('FinalState.txt','w');  %this opens a file and then the vector components are put in  it.
fprintf(fid, '%s \n\n', t);



fprintf(fid, 'CONCENTRATION (cell) \t\t\t VELOCITY (cell) \t\t\t  CONCENTRATION (medium) \t  VELOCITY (medium) \n\n');

  

fprintf(fid, 'H2O2cyt= %4.3f \t\t\t\t\t VGPXcyt=%4.2f \n',M(99),M(30));
fprintf(fid,' \t\t\t\t\t\t\t\t VCAT=%4.2f \n\n',M(113));
H2O2out=M(30)+M(113);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTout=%4.2f \n\n',-H2O2out);
fprintf(fid, ' \t\t\t\t\t\t\t\t vh2o2other=%4.2f     \n',M(103));
fprintf(fid, ' \t\t\t\t\t\t\t\t diffh2o2mito_to_cyt=%4.2f     \n',  (0.15/0.85)*M(98));
H2O2in=(0.15/0.85)*M(98)+M(103);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTin=%4.2f \n\n',H2O2in);                                   
                           
                              
fprintf(fid, 'H2O2mito= %4.3f \t\t\t\t VGPXmito=%4.2f \n',M(91),M(53));
fprintf(fid, ' \t\t\t\t\t\t\t\t diffh2o2mito_to_cyt=%4.2f \n',M(98));
H2O2out=M(53)+M(98);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTout=%4.2f \n\n',-H2O2out);
fprintf(fid, ' \t\t\t\t\t\t\t\t VH2O2prod=%4.2f     \n',  M(89));
H2O2in=M(89);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTin=%4.2f \n\n',H2O2in);                                   




fprintf(fid, 'lCYS= %4.2f \t\t\t\t\t VCYSDIOXY=%4.2f \t\t\t  medCYS=%4.2f \t\t\t\t VmedCYSloss=%4.2f \n',M(4),M(17),M(11),M(37));
fprintf(fid, ' \t\t\t\t\t\t\t\t VGCL=%4.2f \t\t\t\t\t\t\t\t\t\t\t vmedCYSliv=%4.2f                     \n',M(19),M(27));
CYSoutliver=M(17)+M(19);
CYSoutblood=M(37)+M(27);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTout=%4.2f  \t\t\t\t\t\t\t\t\t\t TOTout=%4.2f \n\n',-CYSoutliver,-CYSoutblood);
fprintf(fid, ' \t\t\t\t\t\t\t\t VmedCYSliv=%4.2f  \t\t\t\t\t\t\t\t\t VCYSinput=%4.2f            \n',(vmed/vcell)*M(27),0);
fprintf(fid, ' \t\t\t\t\t\t\t\t VCYSmeth=%4.2f \t\t\t\t\t\t\t\t\t\t VmedGSHdeg=%4.2f                     \n',100,M(33));
CYSinliver=(vmed/vcell)*M(27)+100;
CYSinblood=M(33);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTin=%4.2f \t\t\t\t\t\t\t\t\t\t\t TOTin=%4.2f \n\n\n',CYSinliver,CYSinblood);


fprintf(fid, 'GluCys=%4.2f \t\t\t\t\t VGSS=%4.2f \n',M(6),M(20));
GluCysout=M(20);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTout=%4.2f  \n\n',-GluCysout);
fprintf(fid, '\t\t\t\t\t\t\t     VGCL=%4.2f  \n',M(19)/0.85);
GluCysin=M(19)/0.85;
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTin=%4.2f   \n\n\n',GluCysin);


fprintf(fid, 'cytGSH= %4.2f \t\t\t\t gshcyttomito=%4.2f  \t\t  medGSH=%6.4f \t\t\t VmedGSHdeg=%4.2f  \n',M(5),M(51),M(8),M(33));
fprintf(fid, ' \t\t\t\t\t\t\t\t VlivGSHb=%4.2f \t\t\t\t\t\t\t\t\t  \n',M(28));
fprintf(fid, ' \t\t\t\t\t\t\t\t cytVGPX=%4.2f  \n',2*M(30));
GSHoutliver=M(28)+2*M(30)+M(51);
GSHoutblood=M(33);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTout=%4.2f \t\t\t\t\t\t\t\t\t\t TOTout=%4.2f \n\n',-GSHoutliver,-GSHoutblood);
fprintf(fid, ' \t\t\t\t\t\t\t\t VGSS=%4.2f \t\t\t\t\t\t\t\t\t\t\t VlivGSHb=%4.2f \n',  M(20),(vcell/vmed)*M(28));
fprintf(fid, ' \t\t\t\t\t\t\t\t \t\t\t\t\t\t\t\t\t\t\t\t\t\t VlivGSSGb=%4.2f \n',  2*(vcell/vmed)*M(29));
fprintf(fid, ' \t\t\t\t\t\t\t\t cytVGR= %4.2f   \n',2*M(31));
fprintf(fid, ' \t\t\t\t\t\t\t\t gshmitotocyt= %4.2f   \n',M(55)*(0.15/0.85));
GSHinliver=M(20)+2*M(31)+M(55)*(0.15/0.85);
GSHinblood=(vcell/vmed)*M(28)+(vcell/vmed)*2*M(29);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTin=%4.2f \t\t\t\t\t\t\t\t\t\t TOTin=%4.2f \n\n\n',GSHinliver,GSHinblood);

fprintf(fid, 'mitoGSH= %4.2f  \t\t\t\t mitoVGPX=%4.2f \n',M(57),2*M(53));
fprintf(fid, ' \t\t\t\t\t\t\t\t gshmitotocyt= %4.2f   \n',M(55));
mitoGSHoutliver=2*M(53)+M(55);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTout=%4.2f \n\n',-mitoGSHoutliver);
fprintf(fid, ' \t\t\t\t\t\t\t\t mitoVGR= %4.2f    \n',2*M(54));
fprintf(fid, ' \t\t\t\t\t\t\t\t gshcyttomito= %4.2f   \n',M(51)*(0.85/0.15));
mitoGSHinliver=2*M(54)+M(51)*(0.85/0.15);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTin=%4.2f  \n\n\n',mitoGSHinliver);

fprintf(fid, 'cytGSSG=%4.2f \t\t\t\t\t cytVGR=%4.2f \n',M(12),M(31));
fprintf(fid, '\t\t\t\t\t\t\t     VlivGSSGb=%4.2f  \n',M(29));
GSSGoutliver=M(31)+M(29);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTout=%4.2f  \n\n',-GSSGoutliver);
fprintf(fid, '\t\t\t\t\t\t\t\t cytVGPX=%4.2f\n', M(30));
GSSGinliver=M(30);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTin=%4.2f   \n\n\n',GSSGinliver);

fprintf(fid, 'mitoGSSG=%4.2f \t\t\t\t\t mitoVGR=%4.2f \n',M(58),M(54));
mitoGSSGoutliver=M(54);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTout=%4.2f  \n\n',-mitoGSSGoutliver);
fprintf(fid, '\t\t\t\t\t\t\t\t mitoVGPX=%4.2f\n', M(53));
mitoGSSGinliver=M(53);
fprintf(fid, ' \t\t\t\t\t\t\t\t TOTin=%4.2f   \n\n\n',mitoGSSGinliver);





fclose(fid);
open FinalState.txt;   %this opens a window with the info in it
  

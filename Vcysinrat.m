function a = Vcysinrat(b);




Vmax=273;


km = 2100; %ASC transporter (Kilberg80)


a = Vmax*b./(km + b);




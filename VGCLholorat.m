% Glutamate cysteine ligase (GCL) 
% for rat kidney; could not find info in liver; the Huang93 source is cited in the Shelly Lu review

% Note that inhibition here is competitive; in Chen05 mouse paper it is noncompetitive.

function a = VGCLholorat(cys, gsh_cyt);

V=13702;         %normal conditions: V=13702

atp=3000;
glu=2000;

k1 = 200;       %Km cys, from Huang93
k2 = 1400;       %Km glut, from Huang93
k3= 200;        %km ATP, from Meister76
K2is=8200;     %ki for competitive inhibition of Glu by GSH, from Huang93

     

a = V*(cys./(k1+cys))*(atp./(k3+atp))*(glu./(k2+k2.*gsh_cyt./K2is+glu));


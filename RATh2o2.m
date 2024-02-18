function [T gsh gssg h2o2] =RATh2o2

close all



%MODEL PARAMETERS
%Volumes of compartments
vmed=100e-6; %volume in liters
cellnum=20000;
cellvol=7000e-9; %B10NUMB3R5 website states that the volume of a rat hepatocyte is 7000 cubic micrometers. I have converted this to microL.
vcell=cellnum*cellvol*1e-6; %starting volume of seeded cells, in liters

%Baseline hepatocyte parameters
NADPH = 50;
diffh2o2=2284581;
vh2o2other=3026;
vh2o2prod=28567;

%Baseline medium parameters
gsh_deg=0;







load yinitialcondition
[T,Y] = ode15s(@rath2o2, [0 3*24], yinit, [], vmed, vcell, diffh2o2, vh2o2other, vh2o2prod, gsh_deg, NADPH);





gsh=0.85*Y(:,5)+0.15*Y(:,3);
gssg=0.85*Y(:,12)+0.15*Y(:,30);
h2o2=0.85*Y(:,32)+0.15*Y(:,21);



figure
subplot(4,1,1)
plot(T,gsh*100/9435)
xlim([8 16])
title('to simulate a total halting of GSH synthesis, set V=0 in VGCLholorat.m') 
xlabel('Hours')
ylabel('GSH (% normal)')

subplot(4,1,2)
plot(T,gsh./gssg)
xlim([8 16])
ylabel('GSH:GSSG ratio')
xlabel('Hours')


subplot(4,1,3)
plot(T,h2o2)
xlim([8 16])
xlabel('Hours')
ylabel('H_{2}O_{2} (\mu M)')

subplot(4,1,4)
plot(T,Y(:,65)*100/cellnum)
xlim([8 16])
%xlabel('Hours')
ylabel('% hepatocyte survival')


format long


%CONCENTRATIONS
cystine_med=Y(length(T),10);
cys_med = Y(length(T),11);
gsh_med = Y(length(T),8);

cys_mz  = Y(length(T),4);
glutamylcys_mz = Y(length(T),6);
gsh_cyt_mz  = Y(length(T),5);
gsh_mito_mz=Y(length(T),3);
gssg_cyt_mz = Y(length(T),12);
gssg_mito_mz=Y(length(T),30);
h2o2_cyt_mz=Y(length(T),32);
h2o2_mito_mz=Y(length(T),21);


%VELOCITIES
vmedCYSl=Vcysinrat(cys_med);
vmedGSHdeg=gsh_deg*gsh_med;

vcysdioxy_mz=Vcysdioxygenase(cys_mz);
vgcl_mz=  VGCLholorat(cys_mz,gsh_cyt_mz);
vgss_mz = VGSSrat(glutamylcys_mz);
vgshcyt_to_mito_mz=Vgshcyttomito(gsh_cyt_mz);
vgshmito_to_cyt_mz=Vgshmitotocyt(gsh_mito_mz);
vdiffh2o2mito_to_cyt_mz=diffh2o2*h2o2_mito_mz;
vmzGSHb=vGSHout_l(gsh_cyt_mz) + vGSHout_h(gsh_cyt_mz); 
vcat_mzcyt=Vcat(h2o2_cyt_mz);
vgpx_mzcyt=VGPXcyt(gsh_cyt_mz,h2o2_cyt_mz);
vgr_mzcyt=VGRcyt(gssg_cyt_mz,NADPH);
vgpx_mzmito=VGPXmito(gsh_mito_mz,h2o2_mito_mz);
vgr_mzmito=VGRmito(gssg_mito_mz,NADPH);
vmzGSSGbile= vGSSGout(gssg_cyt_mz);


M=zeros(116,1);
M(4)=cys_mz;
M(5) = gsh_cyt_mz;
M(6) = glutamylcys_mz;
M(8) = gsh_med;
M(10)=cystine_med;
M(11)=cys_med;
M(12)=gssg_cyt_mz;
M(17)=vcysdioxy_mz;
M(19)=vgcl_mz;
M(20)=vgss_mz;
M(27)=vmedCYSl;
M(28)=vmzGSHb;
M(29)=vmzGSSGbile;
M(30)=vgpx_mzcyt;
M(31)=vgr_mzcyt;
M(33)=vmedGSHdeg;
M(51)=vgshcyt_to_mito_mz;
M(53)=vgpx_mzmito;
M(54)=vgr_mzmito;
M(55)=vgshmito_to_cyt_mz;
M(57)=gsh_mito_mz;
M(58)=gssg_mito_mz;
M(89)=vh2o2prod;
M(91)=h2o2_mito_mz;
M(98)=vdiffh2o2mito_to_cyt_mz;
M(99)=h2o2_cyt_mz;
M(103)=vh2o2other;
M(113)=vcat_mzcyt;


PrintFinalValuesRATinvitro(vcell, vmed, M);






end


%------------------------------------------------------------------------
function dy=rath2o2(t, y, vmed, vcell, diffh2o2, vh2o2other, vh2o2prod, gsh_deg, NADPH);

cellnum=20000;
fracliver=y(65)/cellnum;







death=0.0009;

h2o2=0.85*y(32)+0.15*y(21);


limit=150;


if h2o2 < 1
     deathrate=0;
elseif ((h2o2 >= 1) & (h2o2 < limit))
     deathrate=death*h2o2;
else
     deathrate=death*limit;
end


%if h2o2<limit
%    deathrate=death*h2o2;
%else
%    deathrate=death*limit;
%end




dy=zeros(107,1);

dy(3)=Vgshcyttomito(y(5))*(0.85/0.15)-Vgshmitotocyt(y(3)) - 2*VGPXmito(y(3), y(21)) + 2*VGRmito(y(30), NADPH); %liver mito GSH

dy(30)= VGPXmito(y(3),y(21)) - VGRmito(y(30), NADPH); % liver mito GSSG 

dy(4)=(vmed/vcell)*Vcysinrat(y(11)) - VGCLholorat(y(4),y(5)) + 100 - Vcysdioxygenase(y(4)); % liver cysteine. Extra 100 is input from MET cycle

dy(5)=VGSSrat(y(6)) - 2*VGPXcyt(y(5), y(32)) + 2*VGRcyt(y(12), NADPH) - vGSHout_h(y(5)) - vGSHout_l(y(5)) - Vgshcyttomito(y(5)) + Vgshmitotocyt(y(3))*(0.15/0.85); %liver cytosolic GSH. 

dy(6)=VGCLholorat(y(4),y(5))/0.85 - VGSSrat(y(6)); %liver gamma-GC. 

dy(8)=(0.85*vcell*fracliver/vmed)*vGSHout_h(y(5)) + (0.85*vcell*fracliver/vmed)*vGSHout_l(y(5)) - gsh_deg*y(8) + 2*(0.85*vcell*fracliver/vmed)*vGSSGout(y(12)); % med gsh

dy(10) = 100*y(11) - 100*(20/250.5)*y(10); %med cystine

dy(11)= gsh_deg*y(8) - fracliver*Vcysinrat(y(11)) - 2*100*y(11) + 2*100*(20/250.5)*y(10); % med cysteine

dy(12)=VGPXcyt(y(5),y(32)) - VGRcyt(y(12), NADPH) - vGSSGout(y(12)); % liver cytosolic GSSG

dy(21)= vh2o2prod - diffh2o2*y(21) - VGPXmito(y(3),y(21)); %H2O2 in liver mitochondria. 

dy(32)= (0.15/0.85)*diffh2o2*y(21) - VGPXcyt(y(5), y(32)) + vh2o2other - Vcat(y(32)); %H2O2 in liver cytosol



dy(57) = sin(100*t);  %this is just to ensure numerical accuracy

dy(65) = - deathrate*y(65); % living hepatocytes (neglect regeneration)














end






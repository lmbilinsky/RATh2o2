% we assume VGR doesn't depend on GSH

function a = VGRmito(gssg_mito, cNADPH);

km1 = 107;        % for gssg
km2 = 10.4;       % for NADPH





vmax=44035;






a = vmax*((gssg_mito/km1).*(cNADPH/km2))./(1 + gssg_mito/km1 + cNADPH/km2 + (gssg_mito/km1).*(cNADPH/km2));


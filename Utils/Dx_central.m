% Central derivative for x axis
function ux = Dx_central(u)
h = [-1,-2,-1; 0,0,0; 1,2,1] / 8;
ux = imfilter(u,h,"circular"); % circular
ux(1,:,:) = 0;
ux(end,:,:) = 0;
end
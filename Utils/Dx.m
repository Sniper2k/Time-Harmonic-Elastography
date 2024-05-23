% Forward difference derivative for x axis
function ux = Dx(u)
h = [1; -1; 0]; % forward derivative; flipped because we use conv2 instead of imfilter
ux = conv2(u,h,"same");
ux(end,:,:) = 0;
end
% Adjoint operator to forward difference derivative for x axis
function ux = Dxt(u)
h = [0; -1; 1]; % forward derivative
u(end,:,:) = 0;
ux = conv2(u,h,"same");
end
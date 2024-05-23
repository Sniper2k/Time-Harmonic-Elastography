% Adjoint operator to forward difference derivative for y axis
function uy = Dyt(u)
h = [0;-1; 1]; % forward derivative
u(:,end,:) = 0;
uy = conv2(u,h',"same");
end
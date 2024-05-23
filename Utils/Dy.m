% Forward difference derivative for y axis
function uy = Dy(u)
h = [1; -1; 0]; % forward derivative
uy = conv2(u,h',"same");
uy(:,end,:) = 0;
end
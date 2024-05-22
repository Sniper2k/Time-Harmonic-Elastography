function uy = Dy(u)
% uy = circshift(u,-1,2) - u;
% uy(:,end,:) = 0; 

% Central difference
% uy = (circshift(u,-1,2) - circshift(u,1,2))/2;

h = [0,0,0; -1,-2,-1; 1,2,1] / 4;
h = [1; -1; 0]; % forward derivative
% uy = imfilter(u,h',"circular");
uy = conv2(u,h',"same");
% uy(1,:,:) = uy(2,:,:);
% uy(end,:,:) = uy(end-1,:,:);
% uy(:,1,:,:) = uy(:,2,:,:);
% uy(:,end,:,:) = uy(:,end-1,:,:);

% uy(:,1,:) = 0;
uy(:,end,:) = 0;
end
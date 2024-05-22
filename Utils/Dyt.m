function uy = Dyt(u)
% u(:,end,:) = 0; 
% uy = circshift(u,1,2) - u;

% Central difference
% uy = (circshift(u,1,2) - circshift(u,-1,2))/2;

h = [-1,-2,-1; 1,2,1; 0,0,0] / 4;
h = [0;-1; 1]; % forward derivative
% u(:,1,:) = 0;
u(:,end,:) = 0;

% u(2,:,:) = u(2,:,:)+u(1,:,:);
% u(1,:,:) = 0;
% u(end-1,:,:) = u(end-1,:,:)+u(end,:,:);
% u(end,:,:) = 0;
% u(:,2,:,:) = u(:,2,:,:)+u(:,1,:,:);
% u(:,1,:,:) = 0;
% u(:,end-1,:,:) = u(:,end-1,:,:)+u(:,end,:,:);
% u(:,end,:,:) = 0;
% uy = imfilter(u,-h',"circular");
uy = conv2(u,h',"same");
end
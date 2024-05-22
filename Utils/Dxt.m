function ux = Dxt(u)
% u(end,:,:) = 0; 
% ux = circshift(u,1,1) - u;

% Central difference
% ux = (circshift(u,1,1) - circshift(u,-1,1))/2;

% h = [-1,-2,-1; 0,0,0; 1,2,1] / 8;
h = [-1,-2,-1; 1,2,1; 0,0,0] / 4;
h = [0; -1; 1]; % forward derivative

% u(1,:,:) = 0;
u(end,:,:) = 0;

% u(2,:,:) = u(2,:,:)+u(1,:,:);
% u(1,:,:) = 0;
% u(end-1,:,:) = u(end-1,:,:)+u(end,:,:);
% u(end,:,:) = 0;
% u(:,2,:,:) = u(:,2,:,:)+u(:,1,:,:);
% u(:,1,:,:) = 0;
% u(:,end-1,:,:) = u(:,end-1,:,:)+u(:,end,:,:);
% u(:,end,:,:) = 0;
% ux = imfilter(u,-h,"circular");
ux = conv2(u,h,"same");
end
function ux = Dx(u)
% Forward differences
% ux = circshift(u,-1,1) - u;
% ux(end,:,:) = 0; 

% Central difference (seems to produce some artifacts
% ux = (circshift(u,-1,1) - circshift(u,1,1))/2;

% h = [-1,-2,-1; 0,0,0; 1,2,1] / 8;
h = [0,0,0; -1,-2,-1; 1,2,1] / 4; % forward derivative
h = [1; -1; 0]; % forward derivative; flipped because we use conv2 instead of imfilter
% ux = imfilter(u,h,"circular"); %,"circular"
ux = conv2(u,h,"same");
% ux(1,:,:) = ux(2,:,:);
% ux(end,:,:) = ux(end-1,:,:);
% ux(:,1,:,:) = ux(:,2,:,:);
% ux(:,end,:,:) = ux(:,end-1,:,:);

% ux(1,:,:) = 0;
ux(end,:,:) = 0;
end
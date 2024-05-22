function [a_x, a_y] = optical_flow_HS_matlab(I,bin,lambda, weights, itmax, tol, lsqr_maxit, a_0)
d = size(I,[1 2]);
nt = size(I,3);

v_x = zeros([d nt]);
v_y = zeros([d nt]);

flowModel = opticalFlowHS;%Farneback;
flowModel.MaxIteration = itmax;
flowModel.Smoothness = lambda;

% flowModel = opticalFlowFarneback;
% flowModel.NumIterations = itmax;
% flowModel.NumPyramidLevels = 5;
% flowModel.PyramidScale = .8;

% Initialize with first image as we assume periodic movement
estimateFlow(flowModel, I(:,:,1));

for k = 1:nt
  flow = estimateFlow(flowModel, I(:,:,mod(k,nt)+1));
  % Swap directions to be consistent with other code
  v_x(:,:,k) = flow.Vy;
  v_y(:,:,k) = flow.Vx;
%   v(:,:,k) = flow.Magnitude;
end

tmp = fft(v_x, [], 3);
a_x = 2 * tmp(:,:, (bin-1)+1);
tmp = fft(v_y, [], 3);
a_y = 2 * tmp(:,:, (bin-1)+1);


%% plot the images
% figure(102)
% for k = 1:size(v_x,3)
%   imagesc(abs(v_x(:,:,k)))
% %   colorbar
%   title('abs v_x')
%   colorbar()
%   clim([0,max(abs(v_x(:)))])
%   drawnow
%   pause(0.02)
% end
end
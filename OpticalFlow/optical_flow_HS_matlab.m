%% Adoption of Matlab optical flow implementation for time-harmonic optical flow
% Inputs:
%
% I: d_1 x d_2 x nt volume containing nt images of size d_1 x d_2  
%
% bin: in {1, ..., nt}, frequency of the time-harmonic oscilation.
% Equivalently, 1 + the number of period repetitions in the observed images
% 
% lambda: >= 0, regularization parameter for Horn-Schunck functional between two consequtive images 
%
% weights: unused
%
% itmax: integer>= 0, the number of Horn-Schunk iterations.
% 
% tol: >=0, unused
%
% lsqr_maxit: unused
%
% a_0: unused 

function [a_x, a_y] = optical_flow_HS_matlab(I,bin,lambda, weights, itmax, tol, lsqr_maxit, a_0)
d = size(I,[1 2]);
nt = size(I,3);

v_x = zeros([d nt]);
v_y = zeros([d nt]);

flowModel = opticalFlowHS;
flowModel.MaxIteration = itmax;
flowModel.Smoothness = lambda;

% Solve optical flow for each pair of images

% Initialize with first image as we assume periodic movement
estimateFlow(flowModel, I(:,:,1));

for k = 1:nt
  flow = estimateFlow(flowModel, I(:,:,mod(k,nt)+1));
  % Swap directions to be consistent with other code
  v_x(:,:,k) = flow.Vy;
  v_y(:,:,k) = flow.Vx;
end

% Extract single frequency component

tmp = fft(v_x, [], 3);
a_x = 2 * tmp(:,:, (bin-1)+1);
tmp = fft(v_y, [], 3);
a_y = 2 * tmp(:,:, (bin-1)+1);

end
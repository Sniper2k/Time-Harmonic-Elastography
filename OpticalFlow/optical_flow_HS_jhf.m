%% Adoption of PDHGM optical flow implementation for time-harmonic optical flow
% Inputs:
%
% I: d_1 x d_2 x nt volume containing nt images of size d_1 x d_2  
%
% bin: in {1, ..., nt}, frequency of the time-harmonic oscilation.
% Equivalently, 1 + the number of period repetitions in the observed images
% 
% lambda: >= 0, regularization parameter for regularized determined by "method"
%
% weights: unused
%
% itmax: integer>= 0, the number of Horn-Schunk iterations.
% 
% tol: >=0, unused
%
% omega: >0, the standard deviation of the Gaussian filter applied to the
% images when constructing the pyramid. Recommended omega = 1/sqrt(2*scale)
%
% scale: 0<scale<1, compression rate of the coarse-to-fine pyramid.
%
% level: number of levels in the pyramid
%
% method: the choice of regularizer, see switch below. 
%
% par: [0<theta<1, tau>0, sigma >0], PDHGM parameters. 

function [a_x, a_y] = optical_flow_HS_jhf(I,bin,lambda,weights, itmax, tol, omega, scale, level, method, par)
% Use Matlab's own optical flow to compute the velocity between consecutive
% images and take the Fourier transform afterwards
d = size(I,[1 2]);
nt = size(I,3);

v_x = zeros([d nt]);
v_y = zeros([d nt]);

if isempty(par)  
    theta      = 1; 
    tau        = 1/8; 
    sigma      = 1/8;
else
    theta      = par(1);
    sigma      = par(2);
    tau        = par(3);
end

% set parameters (standard values)
% level      = 15;
% iterations = 600;
% scale      = 0.95;
% omega      = 0.7;

% method = 'TV'; % TV, H1, TGV, TV-TV2 or IC

% switch method
% 	case 'TGV'
% 		lambda     = [0.1 5];
% 	case 'TV'
% 		lambda     = [0.1];
% 	case 'TV-TV2'
% 		lambda     = [0.1 0.02];
% 	case 'IC'
% 		lambda     = [0.1 5 0.00001];
% 	case 'H1'
% 		lambda     = [50];
% 	otherwise
% 		error('Wrong method.')
% end

% Solve optical flow for each pair of images

for k = 1:nt
  im1 = I(:,:,k);
  im2 = I(:,:,mod(k,nt)+1);
  [u,v] = coarsetf_all(im1,im2,omega,scale,level,tau,sigma,lambda,theta,itmax,method);

  v_x(:,:,k) = - 2*u;
  v_y(:,:,k) = - 2*v;
end

% Extract single frequency component

tmp = fft(v_x, [], 3);
a_x = tmp(:,:, (bin-1)+1);
tmp = fft(v_y, [], 3);
a_y = tmp(:,:, (bin-1)+1);
end

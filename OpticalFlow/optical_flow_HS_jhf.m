function [a_x, a_y] = optical_flow_HS_jhf(I,bin,lambda,weights, itmax, tol, omega, scale, level, method, par)
% Use Matlab's own optical flow to compute the velocity between consecutive
% images and take the Fourier transform afterwards
d = size(I,[1 2]);
nt = size(I,3);

v_x = zeros([d nt]);
v_y = zeros([d nt]);

% add path
% p1 = genpath('flow_programs_matlab');
% p2 = genpath('flow_programs_matlab/disp_helper');
% p2 = genpath('flow_programs_matlab/disp_helper/cpp');
% addpath(p1,p2);

if isempty(par)  
    theta      = 1;%0.25; % this parameter is never used by JHF
    tau        = 1/8;%10^-4; % parameters of PDHG
    sigma      = 1/8;%1/tau;%10^0;%1/8;
else
    theta      = par(1);%0.25; % this parameter is never used by JHF
    sigma      = par(2);%10^-4; % parameters of PDHG
    tau        = par(3);%1/tau;%10^0;%1/8;
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


for k = 1:nt
%   fprintf('Images: %d %d', k,k+1)
  im1 = I(:,:,k);
  im2 = I(:,:,mod(k,nt)+1);
  [u,v] = coarsetf_all(im1,im2,omega,scale,level,tau,sigma,lambda,theta,itmax,method);

  % Minus in order to be consistent with other coder
  % I don't know why we have to multiply with 2
  v_x(:,:,k) = - 2*u;
  v_y(:,:,k) = - 2*v;
end

tmp = fft(v_x, [], 3);
a_x = tmp(:,:, (bin-1)+1);
tmp = fft(v_y, [], 3);
a_y = tmp(:,:, (bin-1)+1);

% Remove path to avoid problems with afun
% rmpath(p2);
end

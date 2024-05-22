%% Implementation of least squares solver for time-harmonic optical flow with 
% L_2^2 data fidelity and L_2^2 norm of the gradient as regularizer (Model I).
%
% Authors: Oleh Melnyk, Michael Quellmalz 
% Based on 
% [1] Melnyk, Quellmalz, Steidl, Jaitner, Jordan, Sack, Time-Harmonic
% Optical Flow with Applications in Elastography, to appear
%
% Inputs:
%
% I: d_1 x d_2 x nt volume containing nt images of size d_1 x d_2  
%
% bin: in {1, ..., nt}, frequency of the time-harmonic oscilation.
% Equivalently, 1 + the number of period repetitions in the observed images
% 
% lambda: >= 0, regularization parameter for the L_2^2 penalty 
%
% weights: d_1 x d_2 matrix with non-negatice entries. For a ceration pixels, 
% sets the importance to fit errors or, if zero, ignore the error completely. 
% Recommended default w = ones(d), however, can be useful when some regions on 
% the image should be ignored. For instance, in gel_dataset, one could ignore 
% points beyond the circular rim.
%
% itmax: integer>= 0, the number of CG iterations.
% 
% tol: >=0, tolarance parameter for preconditioned conjugate gradient
% method, see details in "pcg"
%
% a_0: d_1 x d_2 x 4, initialization for amplitude to be recovered.
% Contains 4 real components, 1 - x axis real part, 2 - x axis imaginary
% part, 3 - y axis real part, 4 - y axis imaginary part. 

function [a_x, a_y] = optical_flow_HS_harmonic_lsqr(I,bin,lambda, weights, itmax, tol, a_0)
% Set up dimensional parameters
d = size(I,[1 2]);
nt = size(I,3);

% Compute warped derivatives as described in Section 4 of [1]
[Ix,Iy,It] = derivatives_warped(I,a_0,bin);

% Precompute 0, bin and 2*bin-1 Fourier coefficients appearing in equation (11) of [1]

exp1 = repmat(reshape(exp(-2j*pi*(bin-1)*(0:(nt-1))/nt),1,1,[]), d(1), d(2),1);
exp2 = repmat(reshape(exp(-2j*pi*2*(bin-1)*(0:(nt-1))/nt),1,1,[]), d(1), d(2),1);

IxIxf0 = sum(Ix.*Ix,3);
IxIyf0 = sum(Ix.*Iy,3);
IyIyf0 = sum(Iy.*Iy,3);

IxIxf2w = sum(Ix.*Ix.* exp2,3);
IxIyf2w = sum(Ix.*Iy.* exp2,3);
IyIyf2w = sum(Iy.*Iy.* exp2,3);

IxItf = sum(Ix.*It .* exp1,3);
IyItf = sum(Iy.*It .* exp1,3);

% Define a function performing  multiplication with a system matrix C as in (24) and appendix C of [1] 
C = @(x) afun(x,flag,d,nt,lambda,IxIxf0,IxIxf2w,IxIyf0,IxIyf2w,IyIyf0,IyIyf2w, weights);
% Compute b from (24) and appendix C of [1] 
b = [-nt*reshape(real(IxItf),[],1); -nt*reshape(imag(IxItf),[],1); -nt*reshape(real(IyItf),[],1); -nt*reshape(imag(IyItf),[],1)];
dp = prod(d);

% Run CG to solve (24) in [1]
a_long_0 = reshape(a_0,[],1);
[a_long, ~] = pcg(C,b,tol,itmax,[],[],a_long_0);

% Reshape the obtained solution
a_x_re = reshape(a_long(1:dp),d);
a_x_im = reshape(a_long((dp+1):2*dp),d);
a_y_re = reshape(a_long((2*dp+1):3*dp),d);
a_y_im = reshape(a_long((3*dp+1):4*dp),d);

a_x = a_x_re + 1i*a_x_im;
a_y = a_y_re + 1i*a_y_im;
end

% Computes the product with matrix C given the Fourier coefficients
function v = afun(x,flag,d,nt,lambda,IxIxf0,IxIxf2w,IxIyf0,IxIyf2w,IyIyf0,IyIyf2w, weights)
% Reshape the vector
dp = prod(d);
a_x_re = reshape(x(1:dp),d);
a_x_im = reshape(x((dp+1):2*dp),d);
a_y_re = reshape(x((2*dp+1):3*dp),d);
a_y_im = reshape(x((3*dp+1):4*dp),d);

% C_11 block
y_1 = 0.5 * (real(IxIxf0) + real(IxIxf2w)) .* a_x_re;
y_1 = y_1 + 0.5* imag(IxIxf2w) .* a_x_im;
y_1 = y_1 + 0.5* (real(IxIyf0) + real(IxIyf2w)) .* a_y_re;
y_1 = y_1 + 0.5* imag(IxIyf2w) .* a_y_im;
y_1 = y_1 + nt*lambda * Dxt(weights.* Dx(a_x_re));
y_1 = y_1 + nt*lambda * Dyt(weights.*Dy(a_x_re));

% C_12 block
y_2 = 0.5 * imag(IxIxf2w) .* a_x_re;
y_2 = y_2 + 0.5* (real(IxIxf0) - real(IxIxf2w)) .* a_x_im;
y_2 = y_2 + 0.5* imag(IxIyf2w) .* a_y_re;
y_2 = y_2 + 0.5* (real(IxIyf0) - real(IxIyf2w)) .* a_y_im;
y_2 = y_2 + nt*lambda * Dxt(weights.*Dx(a_x_im));
y_2 = y_2 + nt*lambda * Dyt(weights.*Dy(a_x_im));

% C_21 block
y_3 = 0.5 * (real(IxIyf0) + real(IxIyf2w)) .* a_x_re;
y_3 = y_3 + 0.5* imag(IxIyf2w) .* a_x_im;
y_3 = y_3 + 0.5* (real(IyIyf0) + real(IyIyf2w)) .* a_y_re;
y_3 = y_3 + 0.5* imag(IyIyf2w) .* a_y_im;
y_3 = y_3 + nt*lambda * Dxt(weights.*Dx(a_y_re));
y_3 = y_3 + nt*lambda * Dyt(weights.*Dy(a_y_re));

% C_22 block
y_4 = 0.5 * imag(IxIyf2w) .* a_x_re;
y_4 = y_4 + 0.5* (real(IxIyf0) - real(IxIyf2w)) .* a_x_im;
y_4 = y_4 + 0.5* imag(IyIyf2w) .* a_y_re;
y_4 = y_4 + 0.5* (real(IyIyf0) - real(IyIyf2w)) .* a_y_im;
y_4 = y_4 + nt*lambda * Dxt(weights.*Dx(a_y_im));
y_4 = y_4 + nt*lambda * Dyt(weights.*Dy(a_y_im));

% Reshape back to vector
v = [reshape(y_1,dp,1);reshape(y_2,dp,1);reshape(y_3,dp,1);reshape(y_4,dp,1)];
end

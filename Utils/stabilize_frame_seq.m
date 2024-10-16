%% Stilify the frame sequence by subtracting the estimated motion 
% Authors: Oleh Melnyk, Michael Quellmalz 
% Based on constancy equation
% I(t+1, phi(t+1,x)) = I(t, phi(t,x)),
% see Equation 1 from
% [1] Melnyk, Quellmalz, Steidl, Jaitner, Jordan, Sack, Time-Harmonic
% Optical Flow with Applications in Elastography, to appear
%
% Inputs:
%
% I: d_1 x d_2 x nt frame sequence.    
%
% v_x, v_y: d_1 x d_2 x nt estimated velocity fields in both x and y axis.  
%
% varargin: >0, if given, denotes the index of the dimension corresponding to the time component.

function I_stable = stabilize_frame_seq(I, v_x, v_y, varargin)
nt = size(I,3);
% Approximate phi(t+1, x) \approx v(t,phi(t,x)) + phi(t,x) 
%                         \approx v(t,x) + phi(t,x)
%                         = ... = sum_{s=0}^{t} v(s,x)

w = zeros([size(v_x, [1, 2]) 2 size(v_x, 3)], 'double');
if nargin>3 && int64(varargin{1})>0
    k = varargin{1};
    w(:, :, 1, 1:end) = cumsum(v_y, k);
    w(:, :, 2, 1:end) = cumsum(v_x, k);
else
    w(:, :, 1, 1:end) = v_y;
    w(:, :, 2, 1:end) = v_x;
end

% Use image wrapping to compute I(t+1, phi(t+1,x))
I_stable = zeros(size(I));
I_stable(:,:,1) = I(:,:,1);
for i = 2:size(v_x, 3)
    I_stable(:, :, i) = imwarp(double(I(:, :, i)), w(:, :, :, i-1)/nt);
end
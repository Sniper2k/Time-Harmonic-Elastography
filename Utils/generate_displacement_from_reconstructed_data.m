%% Simulation of image displacement
% Authors: Oleh Melnyk, Michael Quellmalz 
% Based on Section 5.1 of
% [1] Melnyk, Quellmalz, Steidl, Jaitner, Jordan, Sack, Time-Harmonic
% Optical Flow with Applications in Elastography, to appear
%
% Inputs:
%
% I_0: d_1 x d_2 image of the initial frame.    
%
% a_x, a_y: d_1 x d_2 amplitudes of the velocity fields in both x and y axis.  
%
% nt: integer>0, the number of frames to generate
%
% bin: in {1, ..., nt}, frequency of the time-harmonic oscilation.
% Equivalently, 1 + the number of period repetitions in the observed images
% 
% varargin: >0, if given denotes the standard deviation of Gauss filter
% applied to derivatives. If not given, defaul is 100.

function [Irec, varargout] = generate_displacement_from_reconstructed_data(I_0, a_x, a_y, nt, bin, varargin)
d = size(I_0, [1 2]);

% Initialize variables

[X, Y] = meshgrid(1:d(1), 1:d(2));
X = X';
Y = Y';

v_x = zeros(d(1),d(2),nt);
v_y = zeros(d(1),d(2),nt);

% Initial condition (25) in [1] for psi as in 

psi_x = zeros(d(1),d(2),nt+1);
psi_x(:,:,1) = X;
psi_y = zeros(d(1),d(2),nt+1);
psi_y(:,:,1) = Y;


if nargin > 5 && varargin{1} > 0
    scale = varargin{1};
else
    scale = 100;
end

for k=2:(nt+1)
    % Compute velocity at time k-1
    v_x(:,:,k-1) = real( a_x .* exp(2.0j*pi *(bin-1)*(k-2)/nt ));
    v_y(:,:,k-1) = real( a_y .* exp(2.0j*pi *(bin-1)*(k-2)/nt ));

    % Perform a step in the Euler scheme for D_t psi = - Grad psi v
    psi_x(:,:,k) = psi_x(:,:,k-1) - imgaussfilt(Dx_central(psi_x(:,:,k-1)),scale).*v_x(:,:,k-1)/nt - imgaussfilt(Dy_central(psi_x(:,:,k-1)),scale).*v_y(:,:,k-1)/nt;
    psi_y(:,:,k) = psi_y(:,:,k-1) - imgaussfilt(Dx_central(psi_y(:,:,k-1)),scale).*v_x(:,:,k-1)/nt - imgaussfilt(Dy_central(psi_y(:,:,k-1)),scale).*v_y(:,:,k-1)/nt;
end

% Compute displaced images via (26) in [1]
Irec = zeros(d(1),d(2),nt);
Irec(:,:,1) = I_0(:,:);
I2 = griddedInterpolant(X,Y,squeeze(Irec(:,:,1)), 'spline', 'nearest');
for k = 2:nt
  Irec(:,:,k) = I2(squeeze(psi_x(:,:,k)),squeeze(psi_y(:,:,k)));
end

% Store computed velocities for further use
varargout = cell(1,2);
varargout{1} = v_x;
varargout{2} = v_y;
end


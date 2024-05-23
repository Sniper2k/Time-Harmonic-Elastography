% This function computes computes the warped derivatives for given images sequence  
% and approximation of the amplitudes. For details, see Section 4 in
% [1] Melnyk, Quellmalz, Steidl, Jaitner, Jordan, Sack, Time-Harmonic
% Optical Flow with Applications in Elastography, to appear
function [Ix,Iy,It] = derivatives_warped(I,a_0,bin)
Ix = zeros(size(I));
Iy = zeros(size(I));
It = zeros(size(I));
nt = size(I,3);

% For each pair of images:
for k=1:nt 
    % Compute velocity corresponding to a
    u = -(a_0(:,:,1) * cos(2.0*pi *(bin-1)*(k-1)/nt) - a_0(:,:,2) * sin(2.0*pi *(bin-1)*(k-1)/nt)) / nt / 2;
    v = -(a_0(:,:,3) * cos(2.0*pi *(bin-1)*(k-1)/nt) - a_0(:,:,4) * sin(2.0*pi *(bin-1)*(k-1)/nt)) / nt / 2;
    I1 = I(:,:,k);
    I2 = I(:,:,mod(k,nt)+1);
    
    % Compute discrete dervatives 
    Ix(:,:,k) = Dx_central(I2);
    Iy(:,:,k) = Dy_central(I2);

    % Compute warped spatial derivatives
    Ix(:,:,k) = warp_image(Ix(:,:,k), -v, -u);
    Iy(:,:,k) = warp_image(Iy(:,:,k), -v, -u);
    I2 = warp_image(I2, -v, -u);

    % Compute warped time derivative
    It(:,:,k) = -(I1 - I2 - Iy(:,:,k).*v - Ix(:,:,k).*u);
end
end


function [Ix,Iy,It] = derivatives_warped(I,a_0,bin)
  % Adapted from JHF, seems to improve things here
  Ix = zeros(size(I));
  Iy = zeros(size(I));
  It = zeros(size(I));
  nt = size(I,3);
  for k=1:nt 
    u = -(a_0(:,:,1) * cos(2.0*pi *(bin-1)*(k-1)/nt) - a_0(:,:,2) * sin(2.0*pi *(bin-1)*(k-1)/nt)) / nt / 2;
    v = -(a_0(:,:,3) * cos(2.0*pi *(bin-1)*(k-1)/nt) - a_0(:,:,4) * sin(2.0*pi *(bin-1)*(k-1)/nt)) / nt / 2;
    I1 = I(:,:,k);
    I2 = I(:,:,mod(k,nt)+1);
    Ix(:,:,k) = Dx_central(I2);
    Iy(:,:,k) = Dy_central(I2);
    Ix(:,:,k) = warp_image(Ix(:,:,k), -v, -u);
    Iy(:,:,k) = warp_image(Iy(:,:,k), -v, -u);
    I2 = warp_image(I2, -v, -u);
    It(:,:,k) = -(I1 - I2 - Iy(:,:,k).*v - Ix(:,:,k).*u);
  %   It(:,:,k) = -(I1 - I2);
  end
end

% Use central derivatives here to improve accuracy 
% (they seem to cause problems in when used inside LSQR though)
% function ux = Dx(u)
%   h = [-1,-2,-1; 0,0,0; 1,2,1] / 8;
%   ux = imfilter(u,h,"circular"); %,"circular"
%   % ux(1,:,:) = ux(2,:,:);
%   % ux(end,:,:) = ux(end-1,:,:);
%   % ux(:,1,:,:) = ux(:,2,:,:);
%   % ux(:,end,:,:) = ux(:,end-1,:,:);
% end
% 
% function uy = Dy(u)
%   h = [-1,-2,-1; 0,0,0; 1,2,1] / 8;
%   uy = imfilter(u,h',"circular");
%   % uy(1,:,:) = uy(2,:,:);
%   % uy(end,:,:) = uy(end-1,:,:);
%   % uy(:,1,:,:) = uy(:,2,:,:);
%   % uy(:,end,:,:) = uy(:,end-1,:,:);
% end
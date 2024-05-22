
function [Irec, varargout] = generate_displacement_from_reconstructed_data(I, a_x, a_y, nt, bin, varargin)
d = size(I, [1 2]);

[X, Y] = meshgrid(1:d(1), 1:d(2));
X = X';
Y = Y';

v_x = zeros(d(1),d(2),nt);
v_y = zeros(d(1),d(2),nt);
phi_x = zeros(d(1),d(2),nt+1);
phi_x(:,:,1) = X;
phi_y = zeros(d(1),d(2),nt+1);
phi_y(:,:,1) = Y;

if nargin>5 && varargin{1}<0 % Use old code for negative smoothing parameter

a2_x_re = griddedInterpolant(X,Y,squeeze(real(a_x)), 'spline', 'nearest');
a2_x_im = griddedInterpolant(X,Y,squeeze(imag(a_x)), 'spline', 'nearest');
a2_y_re = griddedInterpolant(X,Y,squeeze(real(a_y)), 'spline', 'nearest');
a2_y_im = griddedInterpolant(X,Y,squeeze(imag(a_y)), 'spline', 'nearest');

for k=2:nt
    a_phi_x = a2_x_re(squeeze(phi_x(:,:,k-1)), squeeze(phi_y(:,:,k-1))) + 1j*a2_x_im(squeeze(phi_x(:,:,k-1)), squeeze(phi_y(:,:,k-1)));
    a_phi_y = a2_y_re(squeeze(phi_x(:,:,k-1)), squeeze(phi_y(:,:,k-1))) + 1j*a2_y_im(squeeze(phi_x(:,:,k-1)), squeeze(phi_y(:,:,k-1)));

    v_x(:,:,k-1) = real( a_phi_x .* exp(2.0j*pi *(bin-1)*(k-2)/nt ) );
    v_y(:,:,k-1) = real( a_phi_y .* exp(2.0j*pi *(bin-1)*(k-2)/nt ) );

    phi_x(:,:,k) = phi_x(:,:,k-1) - v_x(:,:,k-1)/nt;
    phi_y(:,:,k) = phi_y(:,:,k-1) - v_y(:,:,k-1)/nt;
end
a_phi_x = a2_x_re(squeeze(phi_x(:,:,nt)), squeeze(phi_y(:,:,nt))) + 1j*a2_x_im(squeeze(phi_x(:,:,nt)), squeeze(phi_y(:,:,nt)));
a_phi_y = a2_y_re(squeeze(phi_x(:,:,nt)), squeeze(phi_y(:,:,nt))) + 1j*a2_y_im(squeeze(phi_x(:,:,nt)), squeeze(phi_y(:,:,nt)));
v_x(:,:,nt) = real( a_phi_x .* exp(2.0j*pi *(bin-1)*(k-2)/nt) );
v_y(:,:,nt) = real( a_phi_y .* exp(2.0j*pi *(bin-1)*(k-2)/nt) );

else

if nargin > 5
    scale = varargin{1};
else
    scale = 100;
end

for k=2:(nt+1)
    v_x(:,:,k-1) = real( a_x .* exp(2.0j*pi *(bin-1)*(k-2)/nt ));
    v_y(:,:,k-1) = real( a_y .* exp(2.0j*pi *(bin-1)*(k-2)/nt ));

    phi_x(:,:,k) = phi_x(:,:,k-1) - imgaussfilt(Dx_central(phi_x(:,:,k-1)),scale).*v_x(:,:,k-1)/nt - imgaussfilt(Dy_central(phi_x(:,:,k-1)),scale).*v_y(:,:,k-1)/nt;
    phi_y(:,:,k) = phi_y(:,:,k-1) - imgaussfilt(Dx_central(phi_y(:,:,k-1)),scale).*v_x(:,:,k-1)/nt - imgaussfilt(Dy_central(phi_y(:,:,k-1)),scale).*v_y(:,:,k-1)/nt;
end

% for k=2:nt  % Approximates Dphi by one
% %     a_phi_x = a2_x_re(squeeze(phi_x(:,:,k-1)), squeeze(phi_y(:,:,k-1))) + 1j*a2_x_im(squeeze(phi_x(:,:,k-1)), squeeze(phi_y(:,:,k-1)));
% %     a_phi_y = a2_y_re(squeeze(phi_x(:,:,k-1)), squeeze(phi_y(:,:,k-1))) + 1j*a2_y_im(squeeze(phi_x(:,:,k-1)), squeeze(phi_y(:,:,k-1)));
% 
%     v_x(:,:,k-1) = real( a_x .* exp(2.0j*pi *(bin-1)*(k-2)/nt ));
%     v_y(:,:,k-1) = real( a_y .* exp(2.0j*pi *(bin-1)*(k-2)/nt ));
% 
%     phi_x(:,:,k) = phi_x(:,:,k-1) - v_x(:,:,k-1)/nt;
%     phi_y(:,:,k) = phi_y(:,:,k-1) - v_y(:,:,k-1)/nt;
% end
end

Irec = zeros(d(1),d(2),nt);
Irec(:,:,1) = I(:,:);
I2 = griddedInterpolant(X,Y,squeeze(Irec(:,:,1)), 'spline', 'nearest');
for k = 2:nt
  Irec(:,:,k) = I2(squeeze(phi_x(:,:,k)),squeeze(phi_y(:,:,k)));
end

varargout = cell(1,2);
varargout{1} = v_x;
varargout{2} = v_y;
end


function [ u, v ] = pdhg_h1( I1, I2, u, v, sigma, tau, lambda, theta, iterations )
% Computes optical flow on one scale of the ctf scheme with a variational model and H1 regularizer.
% For minimization the PDHG algorithm is used.
%
% The data term is based on gray-value constancy assumption.
%
% Parameters:
% I1, I2	two images
% u,v		initial flow (from last scale of ctf scheme)
% sigma, tau	parameters of PDHG
% lambda	regularization parameter
%
% written by
% 	Jan Henrik Fitschen
% 	04/04/2017

% Initialization
[nx,ny] = size(I1);
Dx = spdiags([-ones(nx,1), ones(nx,1)],[0,1],nx,nx);
Dx(end,end)=0;
Dy = spdiags([-ones(ny,1), ones(ny,1)],[0,1],ny,ny);
Dy(end,end)=0;
Dxt = Dx';
Dyt = Dy';

Ix = Dx *   I2      ;
Iy =        I2 * Dyt;
Ix = warp_image(Ix, -v, -u);
Iy = warp_image(Iy, -v, -u);
I2 = warp_image(I2, -v, -u);
It = (I1 - I2 - Iy.*v - Ix.*u);
% It = I1-I2;

p = zeros(4*nx,ny);
p_bar = zeros(4*nx,ny);

% fprintf('-----------------------------------------------------\n')
% 
% data_fid = sum(abs(Ix .* u + Iy .* v + It),'all');
% reg = sum(abs([Dx * u; u * Dyt; Dx * v; v * Dyt]).^2,'all');
% total = data_fid + lambda*reg;
% fprintf('Data Fid.: %f, \t Reg: %f. \t Total: %f\n', data_fid, reg, total)

for k=1:iterations    
    % Step 1: Update u, v

    [u, v] = cMatrixShrink(Ix, Iy, -It, ...
               u - tau * (Dxt * p_bar(1:nx,:) + p_bar((nx+1):2*nx,:) * Dy), ...
               v - tau * (Dxt * p_bar(2*nx+1:3*nx,:) + p_bar(3*nx+1:end,:) * Dy), ...
               1/tau); 
    
    % Step 2: Update p1, p2
    p_old = p;
    %p = (p + sigma * [Dx * u; u * Dyt; Dx * v; v * Dyt])/(1+lambda);  

    % MQ: exclude case with I1 ~ I2
    if abs(sum(sum((p + sigma * [Dx * u; u * Dyt; Dx * v; v * Dyt]).^2))) > 1e-6
    p = (p + sigma * [Dx * u; u * Dyt; Dx * v; v * Dyt])*lambda/(sum(sum((p + sigma * [Dx * u; u * Dyt; Dx * v; v * Dyt]).^2)));
    end

    % Step 3: Extrapolation
    p_bar = p + theta*(p - p_old);
    
%     if mod(k,100) == 1
%         data_fid = sum(abs(Ix .* u + Iy .* v + It),'all');
%         reg = sum(abs([Dx * u; u * Dyt; Dx * v; v * Dyt]).^2,'all');
%         total = data_fid + lambda*reg;
%         step_p = theta.^2*sum(abs(p-p_old).^2,'all')/sum(abs(p_old).^2,'all');
%         fprintf('Data Fid.: %f, \t Reg: %f. \t Total: %f\n', data_fid, reg, total)
%     end
end

end

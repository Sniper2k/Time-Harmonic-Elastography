function [a_x, a_y] = optical_flow_HS_harmonic_l1_l2_irls(I,bin,lambda, weights, itmax, tol, lsqr_maxit, a_0)
d = size(I,[1 2]);
nt = size(I,3);

% I = real(I);
%bin = 11;
% disp(['Bin is ' num2str(bin)])
% internal_medfilt = OME_parameters.parameters.harmonic_hs.internal_medfilt; %[5 5]; is standard
% lambda = OME_parameters.parameters.harmonic_hs.lambda; % smoothing strenght
% internal_medfilt = OME_parameters.parameters.harmonic_hs.internal_medfilt;
% loopNumber = OME_parameters.parameters.harmonic_hs.loopNumber;

%OME_parameters.parameters.harmonic_hs.internal_medfilt(1) = 8;
% res_adapted_filter = round(OME_parameters.parameters.harmonic_hs.internal_medfilt(1)*2.0000e-06/OME_parameters.values.resolution(1),TieBreaker="odd");
% res_adapted_filter = 8;
% internal_medfilt = [res_adapted_filter, res_adapted_filter];
% % kernel_laplacian = OME_parameters.parameters.harmonic_hs.kernel_laplacian;
% % kernel_laplacian = kernel_laplacian/sum(kernel_laplacian(:)); % normalizing kernel
% 
% OME_parameters.parameters.harmonic_hs.kernel_laplacian = ones([res_adapted_filter, res_adapted_filter]);
% OME_parameters.parameters.harmonic_hs.kernel_laplacian(ceil(res_adapted_filter/2), ceil(res_adapted_filter/2)) = 0;
% kernel_laplacian = OME_parameters.parameters.harmonic_hs.kernel_laplacian;
% kernel_laplacian = kernel_laplacian / sum(kernel_laplacian(:));
% try
%     if OME_parameters.parameters.harmonic_hs.Harcode_filter.flag
%         internal_medfilt = OME_parameters.parameters.harmonic_hs.Harcode_filter.values;
% %         filter_size = internal_medfilt(1);
%     end
% catch
%     disp('no hardcoded parameters')
% end
% 
% disp(['Internal mdefilter ist: ' num2str(internal_medfilt)])
% % disp(kernel_laplacian)

% Gauss_filter = gauss2D(internal_medfilt(1),0.8,[0 0],2); % ingolfs version
% Gauss_filter = gauss2D(internal_medfilt(1),2,[0 0],2);
% Gauss_filter = Gauss_filter / sum(Gauss_filter(:));

% I = I / max(I(:));


[Ix,Iy,It] = derivatives_warped(I,a_0,bin);
% Ix = Dx(I);
% Iy = Dy(I);
% It = Dt(I);
% [Iy, Ix, It] = imgradientxyz(I, 'central');

a_long = reshape(a_0,[],1);

% expected_sparsity = round(0.05 * numel(I));

exp1 = repmat(reshape(exp(-2j*pi*(bin-1)*(0:(nt-1))/nt),1,1,[]), d(1), d(2),1);
exp2 = repmat(reshape(exp(-2j*pi*2*(bin-1)*(0:(nt-1))/nt),1,1,[]), d(1), d(2),1);
% exp2_1d = exp(-2j*pi*2*(bin-1)*(0:(nt-1))/nt).';
dp = prod(d);

epsilon = realmax;
epsilon2 = realmax;
% epsilon = 10^3;


% weights for IRLS smoothing parameters
% Dn = sqrt(nt)* max(abs([Ix Iy]),[],"all");
% Rn = 2*prod(d);

% epsilon = Dn*Rn;

for irls_it = 1:itmax
    res = residual(a_long, Ix, Iy, It, bin);
    res2 = residual2(a_long, d, bin, nt);
    
%     eps_new = sum(res,"all") + lambda * sum(res2,'all');
%     eps_new = eps_new / lambda/ (nt*dp+nt);

%     eps_new = sum(res,"all") + lambda * sum(res2,'all');
%     eps_new = eps_new / (nt*dp+lambda*nt * Rn);

%     eps_new = sum(res,"all") / (Dn*nt*dp);
%     eps_new = eps_new + lambda * sum(res2,'all') / Rn / nt;

      eps_new = 0.1*mean(res,"all") ;
      epsilon = min(epsilon, eps_new/sqrt(irls_it)); 
      epsilon = max(epsilon, 10^-8 / sqrt(irls_it));

      eps_new2 = 0.1*mean(res2,'all');
      epsilon2 = min(epsilon2, eps_new2/sqrt(irls_it)); 
      epsilon2 = max(epsilon2, 10^-8 / sqrt(irls_it));

%     eps_new = eps_new / (nt*dp+lambda*nt * Rn);

%     epsilon = min(epsilon, eps_new);
%     epsilon = min(epsilon, eps_new/sqrt(irls_it)); 
%     epsilon = epsilon * 0.1;

%     sigma = max(Dn *lambda*epsilon, res);
%     sigma2 = max(Rn * epsilon,res2);

    sigma = max(epsilon, res);
    sigma2 = max(epsilon2,res2);

%     sigma = ones(size(sigma));
%     sigma2 = ones(size(sigma2));

    IxIxf0 = sum(Ix.*Ix ./ sigma,3);
    IxIyf0 = sum(Ix.*Iy ./ sigma,3);
    IyIyf0 = sum(Iy.*Iy ./ sigma,3);
    
    IxIxf2w = sum(Ix.*Ix ./ sigma .* exp2,3);
    IxIyf2w = sum(Ix.*Iy ./ sigma .* exp2,3);
    IyIyf2w = sum(Iy.*Iy ./ sigma .* exp2,3);
    
    IxItf = sum(Ix.*It ./ sigma .* exp1,3);
    IyItf = sum(Iy.*It ./ sigma .* exp1,3);

    fwR0 = sum( 1.0 ./ sigma2,3);
    fwR2w = sum( exp2 ./ sigma2,3); 

    b = [-nt*reshape(real(IxItf),[],1); -nt*reshape(imag(IxItf),[],1); -nt*reshape(real(IyItf),[],1); -nt*reshape(imag(IyItf),[],1)];

    M = @(x) afun(x,flag,d,nt,lambda,IxIxf0,IxIxf2w,IxIyf0,IxIyf2w,IyIyf0,IyIyf2w, weights,fwR0, fwR2w);
%     warning off
    [a_long,~] = pcg(M,b,tol,lsqr_maxit,[],[],a_long);
%     warning on
end

% a_long_0 = randn(dp*4,1);

% u = u_long(1:d) + 1i*u_long((d+1):2*d);

a_x_re = reshape(a_long(1:dp),d);
a_x_im = reshape(a_long((dp+1):2*dp),d);
a_y_re = reshape(a_long((2*dp+1):3*dp),d);
a_y_im = reshape(a_long((3*dp+1):4*dp),d);

% M(:,:,1,1) = 0.5 * real(IxIxf0) + 0.5* real(IxIxf2w) + 0.5*lambda;
% M(:,:,1,2) = 0.5* imag(IxIxf2w); %-0.5 * imag(IxIxf0) + 0.5* imag(IxIxf2w);
% M(:,:,1,3) = 0.5 * real(IxIyf0) + 0.5* real(IxIyf2w);
% M(:,:,1,4) = 0.5* imag(IxIyf2w); %-0.5 * imag(IxIyf0) + 0.5* imag(IxIyf2w);

% M(:,:,2,1) = 0.5* imag(IxIxf2w); %0.5 * imag(IxIxf0) + 0.5* imag(IxIxf2w);
% M(:,:,2,2) = 0.5 * real(IxIxf0) - 0.5* real(IxIxf2w) + 0.5*lambda;
% M(:,:,2,3) = 0.5* imag(IxIyf2w); % 0.5 * imag(IxIyf0) + 0.5* imag(IxIyf2w);
% M(:,:,2,4) = 0.5 * real(IxIyf0) - 0.5* real(IxIyf2w);

% M(:,:,3,1) = 0.5 * real(IxIyf0) + 0.5* real(IxIyf2w);
% M(:,:,3,2) = 0.5* imag(IxIyf2w);%-0.5 * imag(IxIyf0) + 0.5* imag(IxIyf2w);
% M(:,:,3,3) = 0.5 * real(IyIyf0) + 0.5* real(IyIyf2w) + 0.5*lambda;
% M(:,:,3,4) = 0.5* imag(IyIyf2w);%-0.5 * imag(IyIyf0) + 0.5* imag(IyIyf2w);

% M(:,:,4,1) = 0.5* imag(IxIyf2w);%0.5 * imag(IxIyf0) + 0.5* imag(IxIyf2w);
% M(:,:,4,2) = 0.5 * real(IxIyf0) - 0.5* real(IxIyf2w);
% M(:,:,4,3) = 0.5* imag(IyIyf2w);%0.5 * imag(IyIyf0) + 0.5* imag(IyIyf2w);
% M(:,:,4,4) = 0.5 * real(IyIyf0) - 0.5* real(IyIyf2w) + 0.5*lambda;

% u_re = zeros(si(1),si(2));
% u_im = zeros(si(1),si(2));
% v_re = zeros(si(1),si(2));
% v_im = zeros(si(1),si(2));

% a_x_re = rand(si(1),si(2)) + 1i *rand(si(1),si(2));
% a_x_im = rand(si(1),si(2)) + 1i *rand(si(1),si(2));
% a_y_re = rand(si(1),si(2)) + 1i *rand(si(1),si(2));
% a_y_im = rand(si(1),si(2)) + 1i *rand(si(1),si(2));

a_x = a_x_re + 1i*a_x_im;
a_y = a_y_re + 1i*a_y_im;

end

function res = residual(a_long, Ix, Iy, It, bin)
d = size(Ix, [1 2]);
nt = size(Ix,3);
dp = prod(d);
res = zeros(size(Ix));
a_x = reshape(a_long(1:dp),d) + 1j*reshape(a_long((dp+1):2*dp),d);
a_y = reshape(a_long((2*dp+1):3*dp),d) + 1j*reshape(a_long((3*dp+1):4*dp),d);

for t=1:nt
    res(:,:,t) = abs( Ix(:,:,t) .* real(a_x * exp(2j*pi*(bin-1)*(t-1)/nt)) + Iy(:,:,t) .* real(a_y * exp(2j*pi*(bin-1)*(t-1)/nt)) + nt* It(:,:,t));
end
end

function res = residual2(a_long, d,bin,nt)
dp = prod(d);
a_x = reshape(a_long(1:dp),d) + 1j * reshape(a_long((dp+1):2*dp),d);
a_y = reshape(a_long((2*dp+1):3*dp),d) + 1j * reshape(a_long((3*dp+1):4*dp),d);
res = zeros([d nt]);

for t=1:nt
    v_x = real(a_x * exp(2j*pi*(bin-1)*(t-1)/nt));
    v_y = real(a_y * exp(2j*pi*(bin-1)*(t-1)/nt));
%     res(t) = norm([Dx(v_x) Dx(v_y) Dy(v_x) Dy(v_y)],'fro');
    res(:,:,t) = sqrt(abs(Dx(v_x)).^2 + abs(Dx(v_y)).^2 + abs(Dy(v_x)).^2 + abs(Dy(v_y)).^2);
end
end


function v = afun(x,flag,d,nt,lambda,IxIxf0,IxIxf2w,IxIyf0,IxIyf2w,IyIyf0,IyIyf2w, weights, fwR0, fwR2w)
dp = prod(d);
a_x_re = reshape(x(1:dp),d);
a_x_im = reshape(x((dp+1):2*dp),d);
a_y_re = reshape(x((2*dp+1):3*dp),d);
a_y_im = reshape(x((3*dp+1):4*dp),d);

y_1 = 0.5 * (real(IxIxf0) + real(IxIxf2w)) .* a_x_re;
y_1 = y_1 + 0.5* imag(IxIxf2w) .* a_x_im;
y_1 = y_1 + 0.5* (real(IxIyf0) + real(IxIyf2w)) .* a_y_re;
y_1 = y_1 + 0.5* imag(IxIyf2w) .* a_y_im;
y_1 = y_1 + nt*lambda * Dxt( (fwR0 + real(fwR2w)) .* weights.* Dx(a_x_re));
y_1 = y_1 + nt*lambda * Dyt( (fwR0 + real(fwR2w)) .* weights.* Dy(a_x_re));
y_1 = y_1 + nt*lambda * Dxt( imag(fwR2w) .* weights .* Dx(a_x_im));
y_1 = y_1 + nt*lambda * Dyt( imag(fwR2w) .* weights .* Dy(a_x_im));

y_2 = 0.5 * imag(IxIxf2w) .* a_x_re;
y_2 = y_2 + 0.5* (real(IxIxf0) - real(IxIxf2w)) .* a_x_im;
y_2 = y_2 + 0.5* imag(IxIyf2w) .* a_y_re;
y_2 = y_2 + 0.5* (real(IxIyf0) - real(IxIyf2w)) .* a_y_im;
y_2 = y_2 + nt*lambda * Dxt( (fwR0 - real(fwR2w)) .* weights.*Dx(a_x_im));
y_2 = y_2 + nt*lambda * Dyt( (fwR0 - real(fwR2w)) .* weights.*Dy(a_x_im));
y_2 = y_2 + nt*lambda * Dxt( imag(fwR2w) .* weights.*Dx(a_x_re));
y_2 = y_2 + nt*lambda * Dyt( imag(fwR2w) .* weights.*Dy(a_x_re));


y_3 = 0.5 * (real(IxIyf0) + real(IxIyf2w)) .* a_x_re;
y_3 = y_3 + 0.5* imag(IxIyf2w) .* a_x_im;
y_3 = y_3 + 0.5* (real(IyIyf0) + real(IyIyf2w)) .* a_y_re;
y_3 = y_3 + 0.5* imag(IyIyf2w) .* a_y_im;
y_3 = y_3 + nt*lambda * Dxt((fwR0 + real(fwR2w)) .* weights.*Dx(a_y_re));
y_3 = y_3 + nt*lambda * Dyt((fwR0 + real(fwR2w)) .* weights.*Dy(a_y_re));
y_3 = y_3 + nt*lambda * Dxt(imag(fwR2w) .* weights.*Dx(a_y_im));
y_3 = y_3 + nt*lambda * Dyt(imag(fwR2w) .* weights.*Dy(a_y_im));

y_4 = 0.5 * imag(IxIyf2w) .* a_x_re;
y_4 = y_4 + 0.5* (real(IxIyf0) - real(IxIyf2w)) .* a_x_im;
y_4 = y_4 + 0.5* imag(IyIyf2w) .* a_y_re;
y_4 = y_4 + 0.5* (real(IyIyf0) - real(IyIyf2w)) .* a_y_im;
y_4 = y_4 + nt*lambda * Dxt( (fwR0 - real(fwR2w)) .* weights.*Dx(a_y_im));
y_4 = y_4 + nt*lambda * Dyt( (fwR0 - real(fwR2w)) .* weights.*Dy(a_y_im));
y_4 = y_4 + nt*lambda * Dxt( imag(fwR2w) .* weights.*Dx(a_y_re));
y_4 = y_4 + nt*lambda * Dyt( imag(fwR2w) .* weights.*Dy(a_y_re));

v = [reshape(y_1,dp,1);reshape(y_2,dp,1);reshape(y_3,dp,1);reshape(y_4,dp,1)];
end


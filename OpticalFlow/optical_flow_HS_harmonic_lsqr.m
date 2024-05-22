function [a_x, a_y] = optical_flow_HS_harmonic_lsqr(I,bin,lambda, weights, itmax, tol, a_0)
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


% Ix = Dx(I);
% Iy = Dy(I);
% It = Dt(I);
[Ix,Iy,It] = derivatives_warped(I,a_0,bin);

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

% % FFT der gradienten produkte
% tmp = fft(Ix.*Ix,[],3);
% IxIxf0 = tmp(:,:,1);
% IxIxf2w = tmp(:,:,mod(2*(bin-1),nt)+1);
% tmp = fft(Ix.*Iy,[],3);
% IxIyf0 = tmp(:,:,1);
% IxIyf2w = tmp(:,:,mod(2*(bin-1),nt)+1);
% tmp = fft(Iy.*Iy,[],3);
% IyIyf0 = tmp(:,:,1);
% IyIyf2w = tmp(:,:,mod(2*(bin-1),nt)+1);
% 
% tmp = fft(Ix.*It,[],3);
% IxItf = tmp(:,:,bin);
% tmp = fft(Iy.*It,[],3);
% IyItf = tmp(:,:,bin);

% Compute the system matrix for each pixel
% d=gpuArray(d);
% nt=gpuArray(nt);
% lambda=gpuArray(lambda);
% IxIxf0=gpuArray(IxIxf0);
% IxIxf2w=gpuArray(IxIxf2w);
% IxIyf0=gpuArray(IxIyf0);
% IxIyf2w=gpuArray(IxIyf2w);
% IyIyf0=gpuArray(IyIyf0);
% IyIyf2w=gpuArray(IyIyf2w);
% weights=gpuArray(weights);
M = @(x) afun(x,flag,d,nt,lambda,IxIxf0,IxIxf2w,IxIyf0,IxIyf2w,IyIyf0,IyIyf2w, weights);
b = [-nt*reshape(real(IxItf),[],1); -nt*reshape(imag(IxItf),[],1); -nt*reshape(real(IyItf),[],1); -nt*reshape(imag(IyItf),[],1)];
dp = prod(d);

% a_long_0 = randn(dp*4,1);
a_long_0 = reshape(a_0,[],1);
% a_long_0=gpuArray(a_long_0);
% b=gpuArray(b);

% [a_long, ~] = lsqr(M,b,tol,itmax,[],[],a_long_0);
[a_long, ~] = pcg(M,b,tol,itmax,[],[],a_long_0);
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

function v = afun(x,flag,d,nt,lambda,IxIxf0,IxIxf2w,IxIyf0,IxIyf2w,IyIyf0,IyIyf2w, weights)
dp = prod(d);
a_x_re = reshape(x(1:dp),d);
a_x_im = reshape(x((dp+1):2*dp),d);
a_y_re = reshape(x((2*dp+1):3*dp),d);
a_y_im = reshape(x((3*dp+1):4*dp),d);

y_1 = 0.5 * (real(IxIxf0) + real(IxIxf2w)) .* a_x_re;
y_1 = y_1 + 0.5* imag(IxIxf2w) .* a_x_im;
y_1 = y_1 + 0.5* (real(IxIyf0) + real(IxIyf2w)) .* a_y_re;
y_1 = y_1 + 0.5* imag(IxIyf2w) .* a_y_im;
y_1 = y_1 + nt*lambda * Dxt(weights.* Dx(a_x_re));
y_1 = y_1 + nt*lambda * Dyt(weights.*Dy(a_x_re));

y_2 = 0.5 * imag(IxIxf2w) .* a_x_re;
y_2 = y_2 + 0.5* (real(IxIxf0) - real(IxIxf2w)) .* a_x_im;
y_2 = y_2 + 0.5* imag(IxIyf2w) .* a_y_re;
y_2 = y_2 + 0.5* (real(IxIyf0) - real(IxIyf2w)) .* a_y_im;
y_2 = y_2 + nt*lambda * Dxt(weights.*Dx(a_x_im));
y_2 = y_2 + nt*lambda * Dyt(weights.*Dy(a_x_im));

y_3 = 0.5 * (real(IxIyf0) + real(IxIyf2w)) .* a_x_re;
y_3 = y_3 + 0.5* imag(IxIyf2w) .* a_x_im;
y_3 = y_3 + 0.5* (real(IyIyf0) + real(IyIyf2w)) .* a_y_re;
y_3 = y_3 + 0.5* imag(IyIyf2w) .* a_y_im;
y_3 = y_3 + nt*lambda * Dxt(weights.*Dx(a_y_re));
y_3 = y_3 + nt*lambda * Dyt(weights.*Dy(a_y_re));

y_4 = 0.5 * imag(IxIyf2w) .* a_x_re;
y_4 = y_4 + 0.5* (real(IxIyf0) - real(IxIyf2w)) .* a_x_im;
y_4 = y_4 + 0.5* imag(IyIyf2w) .* a_y_re;
y_4 = y_4 + 0.5* (real(IyIyf0) - real(IyIyf2w)) .* a_y_im;
y_4 = y_4 + nt*lambda * Dxt(weights.*Dx(a_y_im));
y_4 = y_4 + nt*lambda * Dyt(weights.*Dy(a_y_im));

v = [reshape(y_1,dp,1);reshape(y_2,dp,1);reshape(y_3,dp,1);reshape(y_4,dp,1)];
end

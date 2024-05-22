clear
addpath ../../OpticalFlow/ ../../Utils/ 

diary noiseless.log

rng(1234)

%% Load I_0

I = imread('../gel_dataset/Medimex_Sq20_1_7V_5x_surface_C001H001S0001000001.bmp');
I = double(I);
I = imresize(I,0.201);

figure(3)
imagesc(rot90(I))
title('Original function ')
colorbar()

%% Parameters Setup

bin = 1 + 3;
nsteps = 100;
repetitions = 3;
nt = nsteps*repetitions;

d = size(I);

%% Construct amplitude a

[X, Y] = meshgrid(1:d(1), 1:d(2));
X = X';
Y = Y';

sigma = 1650;
h = round(d./2);
H = zeros(d(1),d(2),2); 
H(:,:,1) = h(1);
H(:,:,2) = h(2);

a_x = @(x) + sin(sum((x-H).^2,3)/2000 ) .* 4*0.2.* (x(:,:,1)-H(:,:,1)).*exp( - sum((x-H).^2,3) /2/sigma); 
a_y = @(x) 10.*exp( -sum((x-0.5*H).^2,3)/sigma) ;

coor = zeros(d(1),d(2),2);
coor(:,:,1) = X;
coor(:,:,2) = Y;

a_x_grid = a_x(coor);
a_y_grid = a_y(coor);

cx = [0 21];
cy = [0 12];

% Store ground truth

f = figure(1);
imagesc(rot90(abs(a_x_grid)))
colorbar()
clim(cx)
axis('off')

export_fig a_x_true -png -transparent


figure(2)
imagesc(rot90(abs(a_y_grid)))
colorbar()
clim(cy)
axis('off')

export_fig a_y_true -png -transparent

%% Generate displacement

Ipg = generate_displacement_from_reconstructed_data(I,a_x_grid,a_y_grid, nt, bin, 1);
Ipg(Ipg < 0) = 0;

Ip = Ipg;
Ipn = Ip;

for k = 1:nt
  Ip(:,:,k) = imgaussfilt(Ip(:,:,k),0.65);
end
fprintf('Rel. Noise %.04f \n',sum(abs(Ipg - Ip).^2,"all")/sum(abs(Ipg).^2,"all"))

%% Show images 

figure(4)
for k = 1:nt
  imagesc(squeeze(Ip(:,:,k)))
  title('Displaced function')
  colorbar()
  clim([0,255])
  drawnow
end

%% Reconstruction parameters

weights = ones(d);

%% Model I
tic
[a_x_rec, a_y_rec] = optical_flow_HS_harmonic_lsqr_pyramid(Ip,bin,1000,weights, 50, 10^-6, 1/sqrt(1.6), 0.8, 2,'lsqr');
time = toc;

fprintf('Runtime %.03f (Model I)\n',time)
l2e = sum(abs(a_x_rec(:)-a_x_grid(:)).^2 + abs(a_y_rec(:)-a_y_grid(:)).^2);
l2e = l2e / sum(abs(a_x_grid(:)).^2 + abs(a_y_grid(:)).^2);
ss = density_ssim(a_x_rec,a_y_rec,a_x_grid,a_y_grid);
Irec = generate_displacement_from_reconstructed_data(Ipn(:,:,1),a_x_rec,a_y_rec, nt, bin,1);
m_err = rel_measurement_error(Irec, double(Ipn));
m_ssim = images_ssim(Irec, Ipn);

fprintf('Error: %.03f Rel.meas.err %.07f  SSIM %.03f meas SSIM %0.07f \n', l2e, m_err, ss, m_ssim)

figure(5)
imagesc(rot90(abs(a_x_rec)))
% title('Magnitudes of spatial density a_x reconstructed LSQR')
colorbar()
clim(cx)
axis('off')

export_fig a_x_lsqr -png -transparent

figure(6)
imagesc(rot90(abs(a_y_rec)))
% title('Magnitudes of spatial density a_y reconstructed LSQR')
colorbar()
clim(cy)
axis('off')

export_fig a_y_lsqr -png -transparent

%% Model II
tic
[a_x_rec, a_y_rec] = optical_flow_HS_harmonic_lsqr_pyramid(Ip,bin,0.0008,weights, [500, 100], 10^-6, 1/sqrt(1.6), 0.8, 4,'l1+l2');
time = toc;

fprintf('Runtime %.03f (Model II) \n',time)
l2e = sum(abs(a_x_rec(:)-a_x_grid(:)).^2 + abs(a_y_rec(:)-a_y_grid(:)).^2);
l2e = l2e / sum(abs(a_x_grid(:)).^2 + abs(a_y_grid(:)).^2);
ss = density_ssim(a_x_rec,a_y_rec,a_x_grid,a_y_grid);
Irec = generate_displacement_from_reconstructed_data(Ipn(:,:,1),a_x_rec,a_y_rec, nt, bin,1);
m_err = rel_measurement_error(Irec, double(Ipn));
m_ssim = images_ssim(Irec, Ipn);

fprintf('Error: %.03f Rel.meas.err %.07f  SSIM %.03f meas SSIM %0.07f \n', l2e, m_err, ss, m_ssim)

figure(9)
imagesc(rot90(abs(a_x_rec)))
% title('Magnitudes of spatial density a_x reconstructed Pyramid')
colorbar()
clim(cx)
axis('off')

export_fig a_x_l1l2 -png -transparent

figure(10)
imagesc(rot90(abs(a_y_rec)))
% title('Magnitudes of spatial density a_y reconstructed Pyramid')
colorbar()
clim(cy)
axis('off')

export_fig a_y_l1l2 -png -transparent

%% Model III
tic
[a_x_rec, a_y_rec] = optical_flow_HS_harmonic_lsqr_pyramid(Ip,bin,20,weights, [100, 25], 10^-6, 1/sqrt(1.6), 0.8, 4,'l1');
time = toc;

fprintf('Runtime %.03f (Model III) \n',time)
l2e = sum(abs(a_x_rec(:)-a_x_grid(:)).^2 + abs(a_y_rec(:)-a_y_grid(:)).^2);
l2e = l2e / sum(abs(a_x_grid(:)).^2 + abs(a_y_grid(:)).^2);
ss = density_ssim(a_x_rec,a_y_rec,a_x_grid,a_y_grid);

Irec = generate_displacement_from_reconstructed_data(Ipn(:,:,1),a_x_rec,a_y_rec, nt, bin, 1);
m_err = rel_measurement_error(Irec, double(Ipn));
m_ssim = images_ssim(Irec, Ipn);

fprintf('Error: %.03f Rel.meas.err %.07f  SSIM %.03f meas SSIM %0.07f \n', l2e, m_err, ss, m_ssim)

figure(17)
imagesc(rot90(abs(a_x_rec)))
% title('Magnitudes of spatial density a_x reconstructed Pyramid l1')
colorbar()
clim(cx)
axis('off')

export_fig a_x_l1 -png -transparent

figure(18)
imagesc(rot90(abs(a_y_rec)))
% title('Magnitudes of spatial density a_y reconstructed Pyramid l1')
colorbar()
clim(cy)
axis('off')

export_fig a_y_l1 -png -transparent

%% Matlab
tic
[a_x_rec, a_y_rec] = optical_flow_HS_matlab(Ip,bin,15000, [], 1000);
a_x_rec = medfilt2(real(a_x_rec),[5,5],'symmetric') + 1j*medfilt2(imag(a_x_rec),[5,5],'symmetric');
a_y_rec = medfilt2(real(a_y_rec),[5,5],'symmetric') + 1j*medfilt2(imag(a_y_rec),[5,5],'symmetric');
fprintf('Runtime %.03f (Matlab Flow)\n',toc)
l2e = sum(abs(a_x_rec(:)-a_x_grid(:)).^2 + abs(a_y_rec(:)-a_y_grid(:)).^2);
l2e = l2e / sum(abs(a_x_grid(:)).^2 + abs(a_y_grid(:)).^2);
ss = density_ssim(a_x_rec,a_y_rec,a_x_grid,a_y_grid);
Irec = generate_displacement_from_reconstructed_data(Ipn(:,:,1),a_x_rec,a_y_rec, nt, bin,1);
m_err = rel_measurement_error(Irec, double(Ipn));
m_ssim = images_ssim(Irec, Ipn);

fprintf('Error: %.03f Rel.meas.err %.07f  SSIM %.03f meas SSIM %0.07f \n', l2e, m_err, ss, m_ssim)

figure(13)
imagesc(rot90(abs(a_x_rec)))
% title('Magnitudes of spatial density a_x reconstructed Matlab')
% legend({'Groundtruth', 'reconstructed'})
colorbar()
clim(cx)
axis('off')

export_fig a_x_matlab -png -transparent

figure(14)
imagesc(rot90(abs(a_y_rec)))
% title('Magnitudes of spatial density a_y reconstructed Matlab')
% legend({'Groundtruth', 'reconstructed'})
colorbar()
clim(cy)
axis('off')

export_fig a_y_matlab -png -transparent


%% PDHGM

addpath ../../flow_programs_matlab/ ../../flow_programs_matlab/disp_helper/ ../../flow_programs_matlab/disp_helper/cpp

tic
[a_x_rec, a_y_rec] = optical_flow_HS_jhf(Ip,bin, 1000, [], 50, [], 0.7, .95, 15, 'TV', [0.25 10^-3 10^-3]);
fprintf('Runtime %.03f (JHF)\n',toc)

l2e = sum(abs(a_x_rec(:)-a_x_grid(:)).^2 + abs(a_y_rec(:)-a_y_grid(:)).^2);
l2e = l2e / sum(abs(a_x_grid(:)).^2 + abs(a_y_grid(:)).^2);
ss = density_ssim(a_x_rec,a_y_rec,a_x_grid,a_y_grid);
Irec = generate_displacement_from_reconstructed_data(Ipn(:,:,1),a_x_rec,a_y_rec, nt, bin,1);
m_err = rel_measurement_error(Irec, double(Ipn));
m_ssim = images_ssim(Irec, Ipn);

fprintf('Error: %.03f Rel.meas.err %.07f  SSIM %.03f meas SSIM %0.07f \n', l2e, m_err, ss, m_ssim)

figure(15)
imagesc(rot90(abs(a_x_rec)))
% title('Magnitudes of spatial density a_x reconstructed JHF')
colorbar()
clim(cx)
axis('off')

export_fig a_x_jhf -png -transparent

figure(16)
imagesc(rot90(abs(a_y_rec)))
% title('Magnitudes of spatial density a_y reconstructed JHF')
colorbar()
clim(cy)
axis('off')

export_fig a_y_jhf -png -transparent


addpath ../../OpticalFlow
addpath ../../Utils
%addpath ../../export_fig-3.44/

diary diary.log
rng(1234)

%% Load I

dirname = '../gel_dataset/';

files = dir([dirname '*.bmp']);
I = imread([dirname files(1).name]);
for k = 2:length(files)
  I(:,:,k) = imread([dirname files(k).name]);
end
time_steps = length(files);
imsize = [size(I,1), size(I,2)];

%% Parameters Setup
nsteps = 8;
repetitions = 3;
nt = nsteps*repetitions; 
bin = 1 + repetitions;

%% Show images
figure(1)
for k = 1:nt
  imagesc(I(:,:,k))
  colorbar
  drawnow
end

I = I(:,:,1:nt);

% reduce dimension for faster tests
scale = 1; %0.25;
d = round(size(I, [1,2]) * scale);

% Preprocess with Gaussian filter

Ip = zeros(d(1), d(2), nt);
for k = 1:nt
  Ip(:,:,k) = imresize(imgaussfilt(double(I(:,:,k)),3),scale);
end

% d = size(Ip, [1 2]);

%% Reconstruction parameters
weights = ones(d);

%% Model I

tic
[a_x_rec, a_y_rec] = optical_flow_HS_harmonic_lsqr_pyramid(Ip,bin,100,weights, 500, 10^-9, 1, 0.5, 1,'lsqr');
time = toc;

save a_rec_lsqr a_x_rec a_y_rec

%% Model I evaluation
load a_rec_lsqr.mat

fprintf('Runtime %.03f (Model I)\n',time)
Irec = generate_displacement_from_reconstructed_data(I(:,:,1),a_x_rec,a_y_rec, nt, bin,10);
m_err = rel_measurement_error(Irec, double(I));
m_ss = ssim(Irec, double(I), DynamicRange=max(I(:)));
fprintf('Rel.meas.err %.07f  SSIM: %.03f \n', m_err, m_ss)

figure(5)
imagesc(abs(a_x_rec))
% title('Magnitudes of spatial density a_x reconstructed LSQR')
colorbar()
axis('off')

export_fig a_x_lsqr -png -transparent

figure(6)
imagesc(abs(a_y_rec))
% title('Magnitudes of spatial density a_y reconstructed LSQR')
colorbar()
axis('off')

export_fig a_y_lsqr -png -transparent

a_rec_lsqr_x = a_x_rec;
a_rec_lsqr_y = a_y_rec;

%% Model II

tic
[a_x_rec, a_y_rec] = optical_flow_HS_harmonic_lsqr_pyramid(Ip,bin,0.2,weights, [1000 200], 10^-6, 1, 0.5, 5,'l1+l2');
time = toc;

save a_rec_l1l2 a_x_rec a_y_rec

%% Model II evaluation
load a_rec_l1l2.mat

fprintf('Runtime %.03f (Model II)\n',time)
Irec = generate_displacement_from_reconstructed_data(I(:,:,1),a_x_rec,a_y_rec, nt, bin,10);
m_err = rel_measurement_error(Irec, double(I));
m_ss = ssim(Irec, double(I), DynamicRange=max(I(:)));
fprintf('Rel.meas.err %.07f  SSIM: %.03f \n', m_err, m_ss)

figure(9)
imagesc(abs(a_x_rec))
% title('Magnitudes of spatial density a_x reconstructed LSQR')
colorbar()
axis('off')

export_fig a_x_l1l2 -png -transparent

figure(10)
imagesc(abs(a_y_rec))
% title('Magnitudes of spatial density a_y reconstructed LSQR')
colorbar()
axis('off')

export_fig a_y_l1l2 -png -transparent

a_rec_l1l2_x = a_x_rec;
a_rec_l1l2_y = a_y_rec;

%% Model III
tic
[a_x_rec, a_y_rec] = optical_flow_HS_harmonic_lsqr_pyramid(Ip,bin,20,weights, [300 50], 10^-6, 1, 0.5, 5,'l1');
time = toc;

save a_rec_l1 a_x_rec a_y_rec

%% Model III evaluation
load a_rec_l1.mat

fprintf('Runtime %.03f (Model III)\n',time)
Irec = generate_displacement_from_reconstructed_data(I(:,:,1),a_x_rec,a_y_rec, nt, bin,10);
m_err = rel_measurement_error(Irec, double(I));
m_ss = ssim(Irec, double(I), DynamicRange=max(I(:)));
fprintf('Rel.meas.err %.07f  SSIM: %.03f \n', m_err, m_ss)

figure(7)
imagesc(abs(a_x_rec))
% title('Magnitudes of spatial density a_x reconstructed LSQR')
colorbar()
axis('off')

export_fig a_x_l1 -png -transparent

figure(8)
imagesc(abs(a_y_rec))
% title('Magnitudes of spatial density a_y reconstructed LSQR')
colorbar()
axis('off')

export_fig a_y_l1 -png -transparent

a_rec_l1_x = a_x_rec;
a_rec_l1_y = a_y_rec;

%% Matlab
tic
[a_x_rec, a_y_rec] = optical_flow_HS_matlab(Ip,bin,1000, [], 1000);
time = toc;

save a_rec_matlab a_x_rec a_y_rec

%% Matlab evaluation
load a_rec_matlab.mat

fprintf('Runtime %.03f (Matlab)\n',time)
Irec = generate_displacement_from_reconstructed_data(I(:,:,1),a_x_rec,a_y_rec, nt, bin,10);
m_err = rel_measurement_error(Irec, double(I));
m_ss = ssim(Irec, double(I), DynamicRange=max(I(:)));
fprintf('Rel.meas.err %.07f  SSIM: %.03f \n', m_err, m_ss)

figure(11)
imagesc(abs(a_x_rec))
% title('Magnitudes of spatial density a_x reconstructed LSQR')
colorbar()
axis('off')

export_fig a_x_matlab -png -transparent

figure(12)
imagesc(abs(a_y_rec))
% title('Magnitudes of spatial density a_y reconstructed LSQR')
colorbar()
axis('off')

export_fig a_y_matlab -png -transparent

a_rec_matlab_x = a_x_rec;
a_rec_matlab_y = a_y_rec;

%% PDHGM

addpath ../../flow_programs_matlab/ ../../flow_programs_matlab/disp_helper/ ../../flow_programs_matlab/disp_helper/cpp

tic
[a_x_rec, a_y_rec] = optical_flow_HS_jhf(Ip,bin, 1000, [], 50, [], .7, .95, 15, 'TV',[0.25 10^-3 10^-3]);
time = toc;

save a_rec_jhf a_x_rec a_y_rec

%% PDHGM evaluation
load a_rec_jhf.mat

fprintf('Runtime %.03f (JHF)\n',time)

Irec = generate_displacement_from_reconstructed_data(I(:,:,1),a_x_rec,a_y_rec, nt, bin,10);
m_err = rel_measurement_error(Irec, double(I));
m_ss = ssim(Irec, double(I), DynamicRange=max(I(:)));
fprintf('Rel.meas.err %.07f  SSIM: %.03f \n', m_err, m_ss)

figure(13)
imagesc(abs(a_x_rec))
% title('Magnitudes of spatial density a_x reconstructed LSQR')
colorbar()
axis('off')

export_fig a_x_jhf -png -transparent

figure(14)
imagesc(abs(a_y_rec))
% title('Magnitudes of spatial density a_y reconstructed LSQR')
colorbar()
axis('off')

export_fig a_y_jhf -png -transparent

space = zeros(d(1),5,nt);
save_video(repmat([double(I) space Irec],1,1,10), strcat('gel_jhf'))

a_rec_jhf_x = a_x_rec;
a_rec_jhf_y = a_y_rec;

diary off

%% Save video files

[Irec_lsqr,vx_lsqr,vy_lsqr] = generate_displacement_from_reconstructed_data(I(:,:,1),a_rec_lsqr_x,a_rec_lsqr_y, nt, bin,10);
Irec_lsqr = rescale(Irec_lsqr, 0.25);
v_lsqr = rescale(sqrt(vx_lsqr.^2+vy_lsqr.^2), 0.25);

[Irec_l1,vx_l1,vy_l1] = generate_displacement_from_reconstructed_data(I(:,:,1),a_rec_l1_x,a_rec_l1_y, nt, bin,10);
Irec_l1 = rescale(Irec_l1, 0.25);
v_l1 = rescale(sqrt(vx_l1.^2+vy_l1.^2), 0.25);

[Irec_l1l2,vx_l1l2,vy_l1l2] = generate_displacement_from_reconstructed_data(I(:,:,1),a_rec_l1l2_x,a_rec_l1l2_y, nt, bin,10);
Irec_l1l2 = rescale(Irec_l1l2, 0.25);
v_l1l2 = rescale(sqrt(vx_l1l2.^2+vy_l1l2.^2), 0.25);

[Irec_matlab,vx_matlab,vy_matlab] = generate_displacement_from_reconstructed_data(I(:,:,1),a_rec_matlab_x,a_rec_matlab_y, nt, bin,10);
Irec_matlab = rescale(Irec_matlab, 0.25);
v_matlab = rescale(sqrt(vx_matlab.^2+vy_matlab.^2), 0.25);

[Irec_jhf,vx_jhf,vy_jhf] = generate_displacement_from_reconstructed_data(I(:,:,1),a_rec_jhf_x,a_rec_jhf_y, nt, bin,10);
Irec_jhf = rescale(Irec_jhf, 0.25);
v_jhf = rescale(sqrt(vx_jhf.^2+vy_jhf.^2), 0.25);

I_small = rescale(double(I), 0.25);

vs = zeros(size(Irec_lsqr,1),5,nt);
hs = zeros(5, 3*size(Irec_lsqr,2) + 10, nt);

frames = [I_small vs Irec_matlab vs Irec_jhf ;
          hs;
          Irec_lsqr vs Irec_l1l2 vs Irec_l1];

frames_v = [I_small/max(I_small,[],'all')*max([v_lsqr,v_l1,v_l1l2],[],'all') vs v_matlab vs v_jhf;
          hs;
          v_lsqr vs v_l1l2 vs v_l1];

orig_label = {'Original', 'Matlab' 'PDHGM', 'Model I', 'Model II', 'Model III'};
h = size(I_small, [2 1]) / 2;
s = size(I_small, [2 1]);
s(1) = s(1) + size(hs,1);
s(2) = s(2) + size(vs,2);
orig_pos = [h(1) 0; h(1) + s(1) 0; h(1) + 2*s(1) 0; h(1) s(2); h(1) + s(1) s(2); h(1) + 2*s(1) s(2)];

frames_lab = video_add_label(frames, orig_label, orig_pos);
save_video(repmat(frames_lab,1,1,5), 'gel_video', 24)
[~,cmdout] = system('ffmpeg -y -i gel_video.avi -c:v libx264 -preset slow -crf 16 gel_video.mp4');

frames_lab = video_add_label(frames_v, orig_label, orig_pos);
save_video(repmat(frames_lab,1,1,5), 'gel_velocity', 6)
[~,cmdout] = system('ffmpeg -y -i gel_velocity.avi -c:v libx264 -preset slow -crf 16 gel_velocity.mp4');


function Ip = rescale(I,scale) 
d =round(size(I, [1,2]) * scale);
nt = size(I,3);

Ip = zeros(d(1), d(2), nt);
for k = 1:nt
  Ip(:,:,k) = imresize(I(:,:,k),scale);
end
end

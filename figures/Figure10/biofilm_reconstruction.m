addpath ../../OpticalFlow
addpath ../../Utils

dataset = 'biofilm';
diary diary.log

dirname = '../agar_dataset/';

rng(1234)

%% Load I

files = dir([dirname '*.tif']);
I = imread([dirname files(1).name]);
for k = 2:length(files)
  I(:,:,k) = imread([dirname files(k).name]);
end
time_steps = size(I,3);
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

ci = [min(I(:)), max(I(:))];

% reduce dimension for faster tests
scale = 1; %0.25;
plot_scale = .5;
d =round(size(I, [1,2]) * scale);

% for k = 1:nt
%   save_image(rescale(I(:,:,k),plot_scale), sprintf('I_%s_t=%d',dataset,k));
% end

I = I(:,:,1:nt);

% Preprocess with Gaussian filter

Ip = zeros(d(1), d(2), nt);
for k = 1:nt
  Ip(:,:,k) = imresize(imgaussfilt(double(I(:,:,k)),3),scale);
end

%d = size(Ip, [1 2]);

%% Reconstruction parameters
weights = ones(d);


%% Baseline
% Colorbars
cx = [0 scale*100];
cy = [0 scale*100];

m_err = rel_measurement_error(repmat(double(I(:,:,1)),[1,1,nt]), double(I));
m_ss = ssim(repmat(double(I(:,:,1)),[1,1,nt]), double(I), DynamicRange=max(I(:)));
fprintf('No movement: Rel.meas.err %.07f  SSIM: %.03f \n', m_err, m_ss)

%% Model 1

tic
[a_x_rec, a_y_rec] = optical_flow_HS_harmonic_lsqr_pyramid(Ip,bin,2000,weights, 100, 10^-9, 1, 0.8, 1,'lsqr');
time = toc;

save a_rec_lsqr a_x_rec a_y_rec

%% Model I evaluation
load a_rec_lsqr.mat

name = sprintf('lambda=%g_its=%d', 2000, 100);
fprintf('%s Runtime %.03f (LSQR)\n', name,time)
[Irec,vx,vy] = generate_displacement_from_reconstructed_data(Ip(:,:,1),a_x_rec,a_y_rec, nt, bin,10);
m_err = rel_measurement_error(Irec, double(Ip));
m_ss = ssim(Irec, double(Ip), DynamicRange=max(Ip(:)));
fprintf('Rel.meas.err %.06f  SSIM: %.04f \n', m_err, m_ss)
    
figure(5)
imagesc(abs(a_x_rec))
% title('Magnitudes of spatial density a_x reconstructed LSQR')
colorbar()
clim(cx)
axis('off')
export_fig(['a_x_lsqr_' name], '-png', '-transparent')

figure(6)
imagesc(abs(a_y_rec))
% title('Magnitudes of spatial density a_y reconstructed LSQR')
colorbar()
clim(cy)
axis('off')
export_fig(['a_y_lsqr_' name], '-png', '-transparent')

a_rec_lsqr_x = a_x_rec;
a_rec_lsqr_y = a_y_rec;

%% Model II
tic
[a_x_rec, a_y_rec] = optical_flow_HS_harmonic_lsqr_pyramid(Ip,bin,0.5,weights, [2000 200], 10^-6, 1, 0.5, 5,'l1+l2');
time = toc;

save a_rec_l1l2 a_x_rec a_y_rec

%% Model II evaluation
load a_rec_l1l2.mat

name = sprintf('lambda=%g_its=%d_iti=%d', 0.5, 2000, 200);
fprintf('%s Runtime %.03f (L1L2)\n',name,time)
[Irec,vx,vy] = generate_displacement_from_reconstructed_data(Ip(:,:,1),a_x_rec,a_y_rec, nt, bin,10);
m_err = rel_measurement_error(Irec, double(Ip));
m_ss = ssim(Irec, double(Ip), DynamicRange=max(Ip(:)));
fprintf('Rel.meas.err %.06f  SSIM: %.03f \n', m_err, m_ss)

figure(9)
imagesc(abs(a_x_rec))
colorbar()
clim(cx)
axis('off')
export_fig(['a_x_l1l2_' name], '-png', '-transparent')

figure(10)
imagesc(abs(a_y_rec))
colorbar()
clim(cy)
axis('off')
export_fig(['a_y_l1l2_' name], '-png', '-transparent')

a_rec_l1l2_x = a_x_rec;
a_rec_l1l2_y = a_y_rec;

%% Model III

tic
[a_x_rec, a_y_rec] = optical_flow_HS_harmonic_lsqr_pyramid(Ip,bin,20,weights, [300 50], 10^-6, 1, 0.5, 5,'l1');
time = toc;

save a_rec_l1 a_x_rec a_y_rec

%% Model III evaluation
load a_rec_l1.mat
    
name = sprintf('lambda=%g_its=%d_iti=%d', 20, 300, 50);
fprintf('%s Runtime %.03f (L1)\n',name,time)
[Irec,vx,vy] = generate_displacement_from_reconstructed_data(Ip(:,:,1),a_x_rec,a_y_rec, nt, bin,10);
m_err = rel_measurement_error(Irec, double(Ip));
m_ss = ssim(Irec, double(Ip), DynamicRange=max(Ip(:)));
fprintf('Rel.meas.err %.06f  SSIM: %.03f \n', m_err, m_ss)
    
figure(7)
imagesc(abs(a_x_rec))
colorbar()
clim(cx)
axis('off')
export_fig(['a_x_l1_' name], '-png', '-transparent')

figure(8)
imagesc(abs(a_y_rec))
colorbar()
clim(cy)
axis('off')
export_fig(['a_y_l1_' name], '-png', '-transparent')

a_rec_l1_x = a_x_rec;
a_rec_l1_y = a_y_rec;

%% Matlab
tic
[a_x_rec, a_y_rec] = optical_flow_HS_matlab(Ip,bin,10000, [], 100); % 100000 for scale 0.5
time = toc;

save a_rec_matlab a_x_rec a_y_rec

%% Matlab evaluation
load a_rec_matlab.mat

fprintf('Runtime %.03f (Matlab)\n',time)
[Irec,vx,vy] = generate_displacement_from_reconstructed_data(Ip(:,:,1),a_x_rec,a_y_rec, nt, bin,10);
m_err = rel_measurement_error(Irec, double(Ip));
m_ss = ssim(Irec, double(Ip), DynamicRange=max(Ip(:)));
fprintf('Rel.meas.err %.06f  SSIM: %.03f \n', m_err, m_ss)

figure(11)
imagesc(abs(a_x_rec))
colorbar()
clim(cx)
axis('off')

export_fig a_x_matlab -png -transparent

figure(12)
imagesc(abs(a_y_rec))
colorbar()
clim(cy)
axis('off')

export_fig a_y_matlab -png -transparent

a_rec_matlab_x = a_x_rec;
a_rec_matlab_y = a_y_rec;

%% JHF
addpath ../../flow_programs_matlab/ ../../flow_programs_matlab/disp_helper/ ../../flow_programs_matlab/disp_helper/cpp

tic
[a_x_rec, a_y_rec] = optical_flow_HS_jhf(Ip,bin, 1000, [], 50, [], .7, .95, 15, 'TV',[0.25 10^-3 10^-3]);
time = toc;

save a_rec_jhf a_x_rec a_y_rec

%% PDHGM evaluation
load a_rec_jhf.mat

fprintf('Runtime %.03f (JHF)\n',time)
[Irec,vx,vy] = generate_displacement_from_reconstructed_data(Ip(:,:,1),a_x_rec,a_y_rec, nt, bin,10);
m_err = rel_measurement_error(Irec, double(Ip));
m_ss = ssim(Irec, double(Ip), DynamicRange=max(Ip(:)));
fprintf('Rel.meas.err %.06f  SSIM: %.03f \n', m_err, m_ss)

figure(13)
imagesc(abs(a_x_rec))
colorbar()
clim(cx)
axis('off')
 
export_fig a_x_jhf -png -transparent
 
figure(14)
imagesc(abs(a_y_rec))
colorbar()
clim(cy)
axis('off')

export_fig a_y_jhf -png -transparent

a_rec_jhf_x = a_x_rec;
a_rec_jhf_y = a_y_rec;

diary off

%% Plot velocity

[~, v_x,v_y] = generate_displacement_from_reconstructed_data(Ip(:,:,1), a_rec_lsqr_x, a_rec_lsqr_y, nt, bin);

norm_v = sqrt(v_x.^2+v_y.^2);
for k = 1:nt
  figure(2)
  imagesc(norm_v(:,:,k));colorbar
%   surf(sqrt(v_x(:,:,k).^2+v_y(:,:,k).^2),EdgeColor='none');colorbar
  figure(3)
  U=imresize(v_x(:,:,k), .03);
  V=imresize(v_y(:,:,k), .03);
  quiver(U, V, 4); axis tight

  figure(4)
  plot(norm_v(:,100,k))
  pause(.2)
end

%% Save video files

plot_scale = .5;

[Irec_lsqr,vx_lsqr,vy_lsqr] = generate_displacement_from_reconstructed_data(Ip(:,:,1),a_rec_lsqr_x,a_rec_lsqr_y, nt, bin,10);
Irec_lsqr = rescale(Irec_lsqr, plot_scale);
v_lsqr = rescale(sqrt(vx_lsqr.^2+vy_lsqr.^2), plot_scale);

[Irec_l1,vx_l1,vy_l1] = generate_displacement_from_reconstructed_data(Ip(:,:,1),a_rec_l1_x,a_rec_l1_y, nt, bin,10);
Irec_l1 = rescale(Irec_l1, plot_scale);
v_l1 = rescale(sqrt(vx_l1.^2+vy_l1.^2), plot_scale);

[Irec_l1l2,vx_l1l2,vy_l1l2] = generate_displacement_from_reconstructed_data(Ip(:,:,1),a_rec_l1l2_x,a_rec_l1l2_y, nt, bin,10);
Irec_l1l2 = rescale(Irec_l1l2, plot_scale);
v_l1l2 = rescale(sqrt(vx_l1l2.^2+vy_l1l2.^2), plot_scale);

[Irec_matlab,vx_mat,vy_mat] = generate_displacement_from_reconstructed_data(Ip(:,:,1),a_rec_matlab_x,a_rec_matlab_y, nt, bin,10);
Irec_matlab = rescale(Irec_matlab, plot_scale);
v_matlab = rescale(sqrt(vx_mat.^2+vy_mat.^2), plot_scale);

[Irec_jhf,vx_jhf,vy_jhf] = generate_displacement_from_reconstructed_data(Ip(:,:,1),a_rec_jhf_x,a_rec_jhf_y, nt, bin,10);
Irec_jhf = rescale(Irec_jhf, plot_scale);
v_jhf = rescale(sqrt(vx_jhf.^2+vy_jhf.^2), plot_scale);

I_small = rescale(double(Ip), plot_scale);

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
save_video(repmat(frames_lab,1,1,5), [dataset '_video'], 24)

frames_lab = video_add_label(frames_v, orig_label, orig_pos);
save_video(repmat(frames_lab,1,1,5), [dataset '_velocity'], 6)

%%
cv = .5*max([v_lsqr,v_l1,v_l1l2],[],'all');
for t=1:6
  figure(31)
  imagesc(v_lsqr(:,:,t))
  axis('off')
  clim([0,cv])
  export_fig(sprintf('agar_v_lsqr_t=%d',t),'-png','-transparent')

  figure(32)
  imagesc(v_l1(:,:,t))
  axis('off')
  clim([0,cv])
  export_fig(sprintf('agar_v_l1_t=%d',t),'-png','-transparent')

  figure(33)
  imagesc(v_l1l2(:,:,t))
  axis('off')
  clim([0,cv])
  export_fig(sprintf('agar_v_l1l2_t=%d',t),'-png','-transparent')

  figure(34)
  imagesc(v_matlab(:,:,t))
  axis('off')
  clim([0,cv])
  export_fig(sprintf('agar_v_mat_t=%d',t),'-png','-transparent')

  figure(35)
  imagesc(v_jhf(:,:,t))
  axis('off')
  clim([0,cv])
  export_fig(sprintf('agar_v_jhf_t=%d',t),'-png','-transparent')
end

%% Reverting the distortion

I_stable_lsqr = rescale(stabilize_frame_seq(Ip, vx_lsqr, vy_lsqr,3),plot_scale);
I_stable_mat  = rescale(stabilize_frame_seq(Ip, vx_mat, vy_mat,3),plot_scale);
I_stable_jhf  = rescale(stabilize_frame_seq(Ip, vx_jhf, vy_jhf,3),plot_scale);
I_stable_l1   = rescale(stabilize_frame_seq(Ip, vx_l1, vy_l1,3),plot_scale);
I_stable_l1l2 = rescale(stabilize_frame_seq(Ip, vx_l1l2, vy_l1l2,3),plot_scale);

vs = zeros(size(I_stable_lsqr,1),5,nt);
hs = zeros(5, 3*size(I_stable_lsqr,2) + 10, nt);

frames_still2 = [I_small vs I_stable_mat vs I_stable_jhf;...
  hs; I_stable_lsqr vs I_stable_l1l2 vs I_stable_l1]/16;
frames_still2 = video_add_label(frames_still2, orig_label, orig_pos);

save_video(repmat(frames_still2,1,1,1), [dataset '_still_video'], 6)

function Ip = rescale(I,scale) 
d =round(size(I, [1,2]) * scale);
nt = size(I,3);

Ip = zeros(d(1), d(2), nt);
for k = 1:nt
  Ip(:,:,k) = imresize(I(:,:,k),scale);
end
end

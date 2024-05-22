clear all
addpath ../../OpticalFlow/ ../../Utils/ 

%% Parameters Setup

bin = 1 + 1;

nsteps = 2000;
repetitions = 1;
nt = nsteps*repetitions;

d = [200 200];

%% Construct amplitude a
[X, Y] = meshgrid(1:d(1), 1:d(2));
X = X';
Y = Y';
sigma = 0.01*sum(d.^2);
h = round(d./2);
H = zeros(d(1),d(2),2); 
H(:,:,1) = h(1);
H(:,:,2) = h(2);

a_x = @(x) -240*0.3.*exp( - sum((x-reshape([0.7 1], [1 1 2]).*H).^2,3) /2/sigma); %- 240*0.3j*exp( -(x-3*h/2).^2/2/sigma);
a_y = @(x) abs(-1000*0.3.*exp( - sum((x-H).^2,3) /2/sigma) + 1000*0.3.*exp( - sum((x-H).^2,3) /3/sigma)); 


coor = zeros(d(1),d(2),2);
coor(:,:,1) = X;
coor(:,:,2) = Y;

a_x_grid = a_x(coor);
a_y_grid = a_y(coor);

figure(1)
imagesc(abs(a_x_grid))
title('Magnitude of spatial density $a_x$')
colorbar()

figure(2)
imagesc(abs(a_y_grid))
title('Magnitude of spatial density $a_y$')
colorbar()

%% Construct I_0 and its displaced sequence I

I = 2 + sin(2 * pi * coor(:,:,1) ./ d(1))  + cos(2 * pi * coor(:,:,2) ./ d(2)) + 4*sqrt( ((coor(:,:,1) - h(1)) /d(1)).^2 + ((coor(:,:,2)-h(2))/d(2)).^2);
nvals = 10;
vals = linspace(min(I,[],'all'),max(I,[],'all'),nvals);
[~, idx] = min(abs(repmat(I,[1 1 nvals]) - repmat(reshape(vals,[1 1 nvals]),[d(1) d(2) 1])),[],3);
I = vals(idx);


[Ip, v_x, v_y] = generate_displacement_from_reconstructed_data(I, a_x_grid, a_y_grid, nt, bin, 10);


%% Generate velocity v directly and the function phi via Euler scheme

[X, Y] = meshgrid(1:d(1), 1:d(2));
X = X';
Y = Y';

v_x_func = @(x,t) real(a_x(x) * exp(2.0j*pi *(bin-1)*(t-1)/nt) );
v_y_func = @(x,t) real(a_y(x) * exp(2.0j*pi *(bin-1)*(t-1)/nt) );

phi_x = zeros([d nt+1]);
phi_y = zeros([d nt+1]);
phi_x(:,:,1) = X;
phi_y(:,:,1) = Y;

for k=2:(nt+1)
    coor = zeros(d(1),d(2),2);
    coor(:,:,1) = phi_x(:,:,k-1);
    coor(:,:,2) = phi_y(:,:,k-1);

    phi_x(:,:,k) = phi_x(:,:,k-1) + v_x_func(coor,k-1)/nt;
    phi_y(:,:,k) = phi_y(:,:,k-1) + v_y_func(coor,k-1)/nt;
end

figure(3)
imagesc(I)
title('Original function ')
colorbar()

clim_vals = [min(I,[],'all'),max(I,[],'all')];

%% Show images

figure(4)
for k = 1:10:nt
  imagesc(squeeze(Ip(:,:,k)))
%   colorbar
  clim(clim_vals)
  title('Displaced function')
  colorbar()
  drawnow
%   pause(0.33)
end

%% Select one point and box around it for pictures later

pos = [120 83];
x_vals = squeeze(phi_x(pos(1),pos(2),:));
y_vals = squeeze(phi_y(pos(1),pos(2),:));

band = [8 8];

figure(20)
scatter(x_vals,y_vals)
daspect([2*band(2) + 1 2*band(1) + 1 1])

%% Extract the box around the selected point and rescale it. 

I_small = Ip((pos(1)-band(1)):(pos(1)+band(1)),(pos(2)-band(2):(pos(2)+band(2))),:);
figure(21)
imagesc(I_small(:,:,1))
clim([min(I,[],'all'),max(I,[],'all')])

factor = 10;
I_scaled = zeros(size(I_small,1)*factor, size(I_small,2)*factor, size(I_small,3)); 
for k = 1:nt
    I_scaled(:,:,k) = imresize(I_small(:,:,k),factor,'nearest');
end

%% Generate images:
% 1) Displaced image
% 2) Box of the displace image depicting the trajectory of a point and its
% current position 

timestamps = [1 500 1000 1500 2000];

for k = 1:length(timestamps)
    t = timestamps(k);

    fig = figure(101);
    imagesc(Ip(:,:,t))
    clim(clim_vals)
    daspect([d(2) d(1) 1])
    hold on
    plot([pos(2)-band(2) pos(2)+band(2) pos(2)+band(2) pos(2)-band(2) pos(2)-band(2)]', ...
        [pos(1)-band(1) pos(1)-band(1) pos(1)+band(1) pos(1)+band(1) pos(1)-band(1)]', '-r','LineWidth',3)
    hold off

    name = strcat('image_',int2str(t),'.png');
    img = frame2im(getframe(gca));
    imwrite(img, name);
    
    x_vals_scaled = factor*(x_vals - pos(1) + band(1) + 0.5);
    y_vals_scaled = factor*(y_vals - pos(2) + band(2) + 0.5);

    fig = figure(105);
    imagesc(I_scaled(:,:,t))
    clim(clim_vals)
    daspect([2*band(2) + 1 2*band(1) + 1 1])
    hold on
    plot(y_vals_scaled(1:t), x_vals_scaled(1:t),'-r.','LineWidth',5, 'MarkerIndices', [t], 'MarkerSize', 80, 'MarkerEdgeColor', 'r')
    hold off

    name = strcat('trajectory_',int2str(t),'.png');
    img = frame2im(getframe(gca));
    imwrite(img, name);
end
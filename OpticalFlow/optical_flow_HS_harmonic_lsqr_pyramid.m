function [a_x, a_y] = optical_flow_HS_harmonic_lsqr_pyramid(I,bin,lambda,weights, itmax, tol, omega, scale, level, method)
% d = size(I, [1 2]);
% nt = size(I, 3);

pyr = pyramidize_mult(I, omega, scale, level);
pyr_weights = pyramidize_mult(weights, omega, scale, level);
d_min = size(pyr{level},[1 2]);
a_0 = randn(d_min(1),d_min(2),4);

itmax_scale = itmax(1);
if numel(itmax)>1
  it_lsqr = itmax(2);
else
  it_lsqr = 50;
end


for i = level:-1:1
    temp1 = pyr{i};
    lambda_lvl = lambda * scale^(i-1);
    d_prev = size(a_0, [1 2]);
    d_cur = size(temp1, [1,2]);
    fac = d_cur./d_prev;
    a_0_old = a_0;
    a_0 = zeros(d_cur(1),d_cur(2),4);
    for k=1:4
        a_0(:,:,k) = imresize(a_0_old(:,:,k), d_cur)*fac(1);
    end

    it_outer = round(itmax_scale/it_lsqr);
    if it_outer<1; it_outer=1; end

    switch method
        case 'lsqr'
            [a_x_rec, a_y_rec] = optical_flow_HS_harmonic_lsqr(temp1,bin,lambda_lvl,pyr_weights{i}, itmax_scale, tol, a_0);
        case 'l2'
            [a_x_rec, a_y_rec] = optical_flow_HS_harmonic_l2_irls(temp1,bin,lambda_lvl,pyr_weights{i}, it_outer, tol, it_lsqr, a_0);
        case 'l1'
            [a_x_rec, a_y_rec] = optical_flow_HS_harmonic_l1_irls(temp1,bin,lambda_lvl,pyr_weights{i}, it_outer, tol, it_lsqr, a_0);
        case 'l1+l2'
            [a_x_rec, a_y_rec] = optical_flow_HS_harmonic_l1_l2_irls(temp1,bin,lambda_lvl,pyr_weights{i}, it_outer, tol, it_lsqr, a_0);
        case 'l1+l2a'
            [a_x_rec, a_y_rec] = optical_flow_HS_harmonic_l1_l2a_irls(temp1,bin,lambda_lvl,pyr_weights{i}, it_outer, tol, it_lsqr, a_0);
    end
        
    
    
    a_0(:,:,1) = medfilt2(real(a_x_rec),[5,5],'symmetric');
    a_0(:,:,2) = medfilt2(imag(a_x_rec),[5,5],'symmetric');
    a_0(:,:,3) = medfilt2(real(a_y_rec),[5,5],'symmetric');
    a_0(:,:,4) = medfilt2(imag(a_y_rec),[5,5],'symmetric');

%     itmax_scale = round(itmax_scale * scale);
% 
%     a_0(:,:,1) = real(a_x_rec);
%     a_0(:,:,2) = imag(a_x_rec);
%     a_0(:,:,3) = real(a_y_rec);
%     a_0(:,:,4) = imag(a_y_rec);

    figure(100)
    imagesc(abs(a_x_rec))
    title(sprintf('Magnitudes of spatial density a_x level %d',i))
    % legend({'Groundtruth', 'reconstructed'})
    colorbar()
    drawnow
    
    figure(101)
    imagesc(abs(a_y_rec))
    title(sprintf('Magnitudes of spatial density a_y level %d',i))
    % legend({'Groundtruth', 'reconstructed'})
    colorbar()
    drawnow
% 
%     pause(0.5)
end

a_x = squeeze(a_0(:,:,1) + 1j*a_0(:,:,2));
a_y = squeeze(a_0(:,:,3) + 1j*a_0(:,:,4));
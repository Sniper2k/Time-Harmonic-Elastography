%% Implementation of the coarse-to-fine approach for time-harmonic optical flow
% Authors: Oleh Melnyk, Michael Quellmalz 
% Based on 
% [1] Melnyk, Quellmalz, Steidl, Jaitner, Jordan, Sack, Time-Harmonic
% Optical Flow with Applications in Elastography, to appear
%
% Inputs:
%
% I: d_1 x d_2 x nt volume containing nt images of size d_1 x d_2  
%
% bin: in {1, ..., nt}, frequency of the time-harmonic oscilation.
% Equivalently, 1 + the number of period repetitions in the observed images
% 
% lambda: >= 0, regularization parameter for the penalty determined by "method"
%
% weights: d_1 x d_2 matrix with non-negatice entries. For a ceration pixels, 
% sets the importance to fit errors or, if zero, ignore the error completely. 
% Recommended default w = ones(d), however, can be useful when some regions on 
% the image should be ignored. For instance, in gel_dataset, one could ignore 
% points beyond the circular rim.
%
% itmax: integer>= 0 or [integer>= 0 integer>= 0], 
% In the first case, it is transformed into [itmax 50].
% The first value denotes the total number of CG iterations per pyramid level. 
% The second value is the number of CG iteration per IRLS iteration. 
% For Model I, second parameter is ignored.
% 
% tol: >=0, tolarance parameter for preconditioned conjugate gradient
% method, see details in "pcg"
%
% omega: >0, the standard deviation of the Gaussian filter applied to the
% images when constructing the pyramid. Recommended omega = 1/sqrt(2*scale)
%
% scale: 0<scale<1, compression rate of the coarse-to-fine pyramid.
%
% level: number of levels in the pyramid
%
% method: data fidelity and regularization terms to be optimized. Possible
% choices:
% 'lsqr' (Model I in [1]) L_2^2 data fidelity and L_2^2 regularization, 
% 'l1+l2' (Model II in [1]) L_1 data fidelity and L_1 regularization, 
% 'l1' (Model III in [1]) L_1 data fidelity and L_2^2 regularization. 

function [a_x, a_y] = optical_flow_HS_harmonic_lsqr_pyramid(I,bin,lambda,weights, itmax, tol, omega, scale, level, method)
% Construcct the pyramid 
pyr = pyramidize_mult(I, omega, scale, level);
pyr_weights = pyramidize_mult(weights, omega, scale, level);

% Random initialization
d_min = size(pyr{level},[1 2]);
a_0 = randn(d_min(1),d_min(2),4);

% If itmax does not contain "CG iterations per IRLS iteration", set to 50
itmax_scale = itmax(1);
if numel(itmax)>1
  it_lsqr = itmax(2);
else
  it_lsqr = 50;
end



for i = level:-1:1
    % Prepare level parameters

    temp1 = pyr{i};
    lambda_lvl = lambda * scale^(i-1);
    
    % Rescale the outcome from the previous level to the size of the
    % current one

    d_prev = size(a_0, [1 2]);
    d_cur = size(temp1, [1,2]);
    fac = d_cur./d_prev;
    a_0_old = a_0;
    a_0 = zeros(d_cur(1),d_cur(2),4);
    for k=1:4
        a_0(:,:,k) = imresize(a_0_old(:,:,k), d_cur)*fac(1);
    end

    % Compute number of IRLS iterations

    it_outer = round(itmax_scale/it_lsqr);
    if it_outer<1; it_outer=1; end

    % Reconstruct with algorithm of choice
    switch method
        case 'lsqr' % Model I
            [a_x_rec, a_y_rec] = optical_flow_HS_harmonic_lsqr(temp1,bin,lambda_lvl,pyr_weights{i}, itmax_scale, tol, a_0);
        case 'l1+l2' % Model II
            [a_x_rec, a_y_rec] = optical_flow_HS_harmonic_l1_l2_irls(temp1,bin,lambda_lvl,pyr_weights{i}, it_outer, tol, it_lsqr, a_0);
        case 'l1' % Model III
            [a_x_rec, a_y_rec] = optical_flow_HS_harmonic_l1_irls(temp1,bin,lambda_lvl,pyr_weights{i}, it_outer, tol, it_lsqr, a_0);
    end
    
    % Postprocess the outcomes with medfilters

    a_0(:,:,1) = medfilt2(real(a_x_rec),[5,5],'symmetric');
    a_0(:,:,2) = medfilt2(imag(a_x_rec),[5,5],'symmetric');
    a_0(:,:,3) = medfilt2(real(a_y_rec),[5,5],'symmetric');
    a_0(:,:,4) = medfilt2(imag(a_y_rec),[5,5],'symmetric');

    % Plot to see the between-level progression 
    figure(100)
    imagesc(abs(a_x_rec))
    title(sprintf('Magnitudes of spatial density a_x level %d',i))
    colorbar()
    drawnow
    
    figure(101)
    imagesc(abs(a_y_rec))
    title(sprintf('Magnitudes of spatial density a_y level %d',i))
    colorbar()
    drawnow
end

% Combine into complex amplitudes

a_x = squeeze(a_0(:,:,1) + 1j*a_0(:,:,2));
a_y = squeeze(a_0(:,:,3) + 1j*a_0(:,:,4));
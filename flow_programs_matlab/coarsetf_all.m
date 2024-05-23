function [ u, v ] = coarsetf_all( im1, im2, omega, scale, level, sigma, tau, lambda, theta, iterations,method)
% Computes optical flow using a variational model with various regularizers.
% For minimization a coarse-to-fine (ctf) scheme is used in combination with the PDHG algorithm
% for the minimization on each scale.
%
% The data term is based on gray-value constancy assumption.
%
% Parameters:
% im1, im2	two images
% omega		blur parameter for creation of the image pyramid
% scale		down-scaling factor in ctf scheme
% level		number of scales in ctf scheme
% sigma, tau	parameters of PDHG
% lambda	array containing regularization parameters (1 parameter except TV-TV2 2 params, TGV 2 params and IC 3 params)
% method	string containing the name of the regularization term ('TV', 'H1', 'TV-TV2', 'TGV', 'IC')
%
% written by
% 	Jan Henrik Fitschen
% 	04/04/2017

pyr1 = pyramidize(im1, omega, scale, level);
pyr2 = pyramidize(im2, omega, scale, level);

u = zeros(size(pyr1{end}));
v = u;

for i = level:-1:1
    temp1 = pyr1{i};
    temp2 = pyr2{i};    
    switch method
	case 'TV'
		fac = size(temp1)./size(u);
		u = imresize(u,size(temp1))*fac(1);
		v = imresize(v,size(temp1))*fac(2);
                [u,v] = pdhg_tv	( temp1, temp2, u, v, 		sigma, tau, lambda, theta, iterations );
	    	v=medfilt2(v,[5,5],'symmetric');
	    	u=medfilt2(u,[5,5],'symmetric');
	case 'H1'
		fac = size(temp1)./size(u);
		u = imresize(u,size(temp1))*fac(1);
		v = imresize(v,size(temp1))*fac(2);
                [u,v] = pdhg_h1		( temp1, temp2, u, v, 		sigma, tau, lambda, theta,iterations );
	    	v=medfilt2(v,[5,5],'symmetric');
	    	u=medfilt2(u,[5,5],'symmetric');
	otherwise 
	    warning('Unknown method.')
    end    
%     disp(['Level ' num2str(level-i+1) ' of ' num2str(level)]);
end

end

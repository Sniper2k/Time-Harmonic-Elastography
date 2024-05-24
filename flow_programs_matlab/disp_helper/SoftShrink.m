function q = SoftShrink(p,lambda)
% Computes componentwise soft-shrinkage of p with parameter lambda.
%
% Authors: 		Jan Henrik Fitschen, Markus Löckel, 04/09/2013

q = (abs(p) > lambda) .* (p - lambda .* sign(p));

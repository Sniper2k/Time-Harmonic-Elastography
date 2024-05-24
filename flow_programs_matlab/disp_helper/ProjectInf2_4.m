function [p] = ProjectInf2_4(q, lambda)
% Project onto inf,2 ball where vectors with 4 entries are grouped 
% (standard tv would be vectors with 2 entries instead of 4)
% e.g. sqrt(q_i^2+q_(i+nx/4)^2+q_(i+nx/2)^2+q_(i+3*nx/4)^2) <= lambda
% written by:           Jan Henrik Fitschen             11/04/2014

[nx, ~] = size(q);
q_abs = sqrt(q(1:nx/4,:).^2 + q(nx/4+1:nx/2,:).^2 + q(nx/2+1:3*nx/4,:).^2 + q(3*nx/4+1:nx,:).^2);
p = q .* repmat(min(1,lambda./q_abs),[4,1]); 
end
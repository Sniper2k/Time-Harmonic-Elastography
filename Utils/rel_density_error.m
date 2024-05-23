% This function computes RE for given pair of 
% reconstructed amplitudes a_rec = (a_rec_x, a_rec_y) 
% and ground truth amplitude a = (a_x, a_y)
function err = rel_density_error(a_x_rec, a_y_rec, a_x, a_y)
err = norm([a_x_rec a_y_rec] - [a_x a_y], 'fro').^2/norm([a_x a_y], 'fro').^2;
end
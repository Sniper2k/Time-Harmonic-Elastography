% This function computes RIE for given pair of 
% simulated images I_rec and original images I
function err = rel_measurement_error(Irec, I)
err = sum(abs(Irec - I).^2, 'all')/sum(abs(I).^2, 'all');
end
function err = rel_measurement_error(Irec, I)
err = sum(abs(Irec - I).^2, 'all')/sum(abs(I).^2, 'all');
end
function err = rel_density_error(a_x_rec, a_y_rec, a_x, a_y)
err = norm([a_x_rec a_y_rec] - [a_x a_y], 'fro')/norm([a_x a_y], 'fro');
end
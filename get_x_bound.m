function [x_bound_all, x_bound_mean, filtered_mean] = get_x_bound(ydata, xdata, n_pixels, n_loci)
    % Convert inputs to arrays
    ydata = double(ydata);
    xdata = double(xdata);
    n_pixels = double(n_pixels);
    n_loci = double(n_loci);

    % % Filter out values where ydata < xdata
    valid_indices = ydata >= xdata;
    ydata(~valid_indices) = 0;
    xdata(~valid_indices) = 0;
    n_pixels(~valid_indices) = 0;
    % ydata = ydata(valid_indices);
    % xdata = xdata(valid_indices);
    % n_pixels = n_pixels(valid_indices);

    if length(n_loci) > 1
       n_loci(~valid_indices) = 0;
    end

    % Calculate x_bound_all and x_bound_xdata
    x_bound_all = sqrt((ydata - xdata) .* n_pixels .* n_loci);
    x_bound_mean = sqrt((ydata - xdata) .* mean(n_pixels) .* n_loci);

    % Calculate the filtered mean
    filtered_mean = xdata;
end
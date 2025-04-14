function variance = variance_model(x, n_pixels, K, O, slope, exponent, n_loci)
% Convert inputs to double precision (equivalent to `np.array(...)`)
    n_pixels = double(n_pixels);
    x = double(x);
    n_loci = double(n_loci);
    
    % Ensure x is non-negative (equivalent to `x_mean = np.abs(x)`)
    x_mean = abs(x);  
   
    % Error checking (equivalent to `if np.any(x_mean < 0) or K < 0: raise ValueError(...)`)
    if any(x_mean < 0) || K < 0
        error('Concentration values (x_mean) and K must be non-negative.');
    end
    
    % Small value to prevent division by zero in case n_loci=0 (equivalent to `epsilon = 1e-9`)
    epsilon = 1e-12;
    
    % Hill equation to calculate fraction of bound proteins
    % Equivalent to:
    % `x_bound = (np.power(x_mean, exponent) * O) / (np.power(K, exponent) + np.power(x_mean, exponent))`
    x_bound = (x_mean.^exponent .* O) ./ (K^exponent + x_mean.^exponent);
    
    % Calculate the variance based on the model
    % Equivalent to:
    % `variance = (x_mean + (x_bound**2)/(n_pixels*(n_loci+epsilon)))`
    variance = x_mean + (x_bound.^2) ./ (n_pixels .* (n_loci + epsilon));
    
    % Scale variance by the slope (equivalent to `return variance * slope`)
    % In most cases, slope = 1 since we convert data from units of
    % intensity/pxl to counts/pxl before fitting
    variance = variance * slope;
end
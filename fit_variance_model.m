function fitReport = fit_variance_model(dataFile, meanField, varField, n_pixels, slope, exponent, n_loci, p0, lb, ub)
% fit_variance_model performs nonlinear curve fitting using lsqcurvefit
% on a model for the variance. It loads the specified .mat file, extracts the data using
% the provided field names, fits the model, and then displays a detailed fit report.
%
% USAGE:
%   fitReport = fit_variance_model(dataFile, meanField, varField, ...
%         n_pixels, slope, exponent, n_loci, p0, lb, ub)
%
% If any argument is omitted, the function will interactively ask for it.
%
% INPUTS:
%   dataFile  - String with the name of the .mat file containing the data.
%   meanField - String with the field name for the independent variable (means).
%   varField  - String with the field name for the dependent variable (variance).
%   n_pixels  - Fixed number of pixels (e.g., 1035).
%   slope     - Fixed slope value (e.g., 1).
%   exponent  - Fixed exponent (e.g., 2).
%   n_loci    - Fixed number of loci (e.g., 1).
%   p0        - Initial guesses for parameters to be fitted [K, O].
%   lb        - Lower bounds for the parameters (e.g., [0, 0]).
%   ub        - Upper bounds for the parameters (e.g., [Inf, Inf]).
%
% OUTPUT:
%   fitReport - Structure containing the fit report details.

    %% Check for missing inputs and interactively prompt if necessary
    if nargin < 1 || isempty(dataFile)
        dataFile = input('Enter the name of the data file (e.g., simulated_binding_data.mat): ', 's');
    end
    if nargin < 2 || isempty(meanField)
        meanField = input('Enter the field name for the mean values (e.g., binding.means): ', 's');
    end
    if nargin < 3 || isempty(varField)
        varField = input('Enter the field name for the variance values (e.g., binding.vars): ', 's');
    end
    if nargin < 4 || isempty(n_pixels)
        n_pixels = input('Enter the value for n_pixels (e.g., 1035): ');
    end
    if nargin < 5 || isempty(slope)
        slope = input('Enter the slope (e.g., 1): ');
    end
    if nargin < 6 || isempty(exponent)
        exponent = input('Enter the exponent (e.g., 2): ');
    end
    if nargin < 7 || isempty(n_loci)
        n_loci = input('Enter the number of loci (e.g., 1): ');
    end
    if nargin < 8 || isempty(p0)
        p0 = input('Enter the initial guesses for [K, O] (e.g., [0.1, 250]): ');
    end
    if nargin < 9 || isempty(lb)
        lb = input('Enter the lower bounds for [K, O] (e.g., [0, 0]): ');
    end
    if nargin < 10 || isempty(ub)
        ub = input('Enter the upper bounds for [K, O] (e.g., [Inf, Inf]): ');
    end

    %% Load Data
    fprintf('\nLoading data from %s...\n', dataFile);
    dataStruct = load(dataFile);
    
    % Use dynamic field referencing to get xdata and ydata
    % For example, if meanField is 'binding.means', we split the string:
    xdata = getFieldByPath(dataStruct, meanField);
    ydata = getFieldByPath(dataStruct, varField);

    %% Define the Variance Model Function
    % Create an anonymous function handle for the model.
    % p(1) corresponds to K, p(2) corresponds to O.
    model = @(p, x) variance_model(x, n_pixels, p(1), p(2), slope, exponent, n_loci);

    %% Set Optimization Options
    opts = optimoptions('lsqcurvefit', ...
        'Display', 'iter', ...
        'Algorithm', 'levenberg-marquardt', ...
        'StepTolerance', 1e-100, ...
        'FunctionTolerance', 1e-100, ...
        'MaxIterations', 100000, ...
        'MaxFunctionEvaluations', 500000);

    %% Perform the Nonlinear Curve Fit
    fprintf('\nStarting lsqcurvefit...\n');
    [p_fit, resnorm, residual, exitflag, output, lambda, jacobian] = ...
        lsqcurvefit(model, p0, xdata, ydata, lb, ub, opts);

    %% Compute Fit Statistics
    dof = length(ydata) - length(p_fit);       % degrees of freedom
    mse = resnorm / dof;                       % mean squared error

    % Estimate covariance matrix and parameter standard errors
    if ~isempty(jacobian)
        covar = inv(jacobian' * jacobian)* mse;
        % Convert standard errors to full array to avoid sparse issues
        param_se = full(sqrt(diag(covar)));
    else
        covar = [];
        param_se = [];
    end

    % R-squared calculation
    ss_tot = sum((ydata - mean(ydata)).^2);
    ss_res = resnorm;
    R_squared = 1 - ss_res / ss_tot;

    %% Create Fit Report Structure
    fitReport.exitFlag = exitflag;
    fitReport.funcCount = output.funcCount;
    fitReport.iterations = output.iterations;
    fitReport.resnorm = resnorm;
    fitReport.degreesOfFreedom = dof;
    fitReport.mse = mse;
    fitReport.R_squared = R_squared;
    fitReport.parameters = struct('K', p_fit(1), 'O', p_fit(2));
    if ~isempty(param_se)
        fitReport.paramSE = struct('K', param_se(1), 'O', param_se(2));
    else
        fitReport.paramSE = [];
    end
    fitReport.covariance = covar;

    %% Display Detailed Fit Report (Dialogue)
    fprintf('\n----- Fit Report -----\n');
    fprintf('Exit Flag: %d\n', fitReport.exitFlag);
    fprintf('Function Evaluations: %d\n', fitReport.funcCount);
    fprintf('Iterations: %d\n', fitReport.iterations);
    fprintf('Residual Norm (Sum of Squared Residuals): %.4e\n', fitReport.resnorm);
    fprintf('Degrees of Freedom: %d\n', fitReport.degreesOfFreedom);
    fprintf('Mean Squared Error: %.4e\n', fitReport.mse);
    fprintf('R-squared: %.4f\n\n', fitReport.R_squared);

    fprintf('Fitted Parameters:\n');
    fprintf('K = %.4e', fitReport.parameters.K);
    if ~isempty(fitReport.paramSE)
        % Convert to full for printing
        fprintf(' ± %.4e\n', full(fitReport.paramSE.K));
    else
        fprintf('\n');
    end
        fprintf('O = %.4e', fitReport.parameters.O);
    if ~isempty(fitReport.paramSE)
        fprintf(' ± %.4e\n', full(fitReport.paramSE.O));
    else
        fprintf('\n');
    end

    if ~isempty(covar)
        fprintf('\nCovariance Matrix:\n');
        disp(covar);
    end
    %% Optional: Plot the Fit Result
    
    % Calculate binned data using quantile-based binning
    num_bins = 100; % Number of bins
    edges = quantile(xdata, linspace(0, 1, num_bins + 1)); % Define bin edges based on quantiles
    bin_centers = (edges(1:end-1) + edges(2:end)) / 2; % Calculate bin centers
    binned_means = zeros(1, num_bins); % Initialize array for binned x means
    binned_stds = zeros(1, num_bins); % Initialize array for binned x standard deviations
    binned_y_means = zeros(1, num_bins); % Initialize array for binned y means
    binned_y_stds = zeros(1, num_bins); % Initialize array for binned y standard deviations
    
    % Loop through each bin to calculate statistics
    for i = 1:num_bins
        bin_indices = xdata >= edges(i) & xdata < edges(i+1); % Find indices of data in the current bin
        if i == num_bins % Include the upper edge for the last bin
            bin_indices = xdata >= edges(i) & xdata <= edges(i+1);
        end
        binned_means(i) = mean(xdata(bin_indices)); % Mean of xdata in the bin
        binned_stds(i) = std(xdata(bin_indices)); % Standard deviation of xdata in the bin
        binned_y_means(i) = mean(ydata(bin_indices)); % Mean of ydata in the bin
        binned_y_stds(i) = std(ydata(bin_indices)); % Standard deviation of ydata in the bin
    end
    % Plot the data
    figure;
    
    % Scatter plot of original data with larger markers
    scatter(xdata, ydata, 100, 'MarkerEdgeColor', '#edb232', 'MarkerFaceColor', '#edb232'); 
    % '50' specifies the marker size
    hold on;
    
    % Plot binned data with error bars
    % Use black color for both binned data points and error bars
    errorbar(binned_means, binned_y_means, binned_y_stds, binned_y_stds, binned_stds, binned_stds, ...
         'o', 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 10);
    % 'LineWidth' controls the thickness of the error bars
    % 'MarkerSize' controls the size of the binned data points

    % Plot the fit curve
    x_fit = linspace(min(xdata), max(xdata), 1000); % Generate x values for the fit curve
    y_fit = model(p_fit, x_fit); % Evaluate the model at the generated x values
    plot(x_fit, y_fit, 'r-', 'LineWidth', 5); % Plot the fit curve as a solid black line

    % Add labels and title
    xlabel('Mean (#/pxl)','FontSize',18); % Label for the x-axis
    ylabel('Variance (#/pxl)','FontSize',18); % Label for the y-axis
    title('Variance Model Fit for TetR-mYFP Full Binding Curve','FontSize',24); % Title of the plot

    % Add legend to distinguish data, binned data, and fit
    legend('Data', 'Binned Data', 'Fit','fontsize',18);

    %Make the x and y ticks bigger
    ax = gca; % Get current axes
    ax.FontSize = 16; % Set the font size for tick labels

    % Activate grid
    grid on
    ax.GridLineWidth = 3;
    hold off; % Release the hold on the current figure

end

%% Helper Function: getFieldByPath
function value = getFieldByPath(dataStruct, fieldPath)
% getFieldByPath extracts a nested field value from a structure.
%   fieldPath is a dot-separated string, e.g., 'binding.means'.
    fields = strsplit(fieldPath, '.');
    value = dataStruct;
    for i = 1:length(fields)
        value = value.(fields{i});
    end
end

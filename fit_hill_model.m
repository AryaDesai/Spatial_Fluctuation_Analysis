function fitReport = fit_hill_model(dataFile, meanField, varField, p0, lb, ub)
% fit_hill_model performs nonlinear curve fitting using lsqcurvefit
% on a Hill function model.
%
% USAGE:
%   fitReport = fit_hill_model(dataFile, meanField, varField, p0, lb, ub)
%
% INPUTS:
%   dataFile  - Name of the .mat file containing data.
%   meanField - Field name for independent variable (x_mean).
%   varField  - Field name for the dependent variable (response).
%   p0        - Initial guesses for [Kd, n, A].
%   lb        - Lower bounds for [Kd, n, A].
%   ub        - Upper bounds for [Kd, n, A].
%
% OUTPUT:
%   fitReport - Structure containing the fit report details.

    %% Check for missing inputs and interactively prompt if necessary
    if nargin < 1 || isempty(dataFile)
        dataFile = input('Enter the data file name: ', 's');
    end
    if nargin < 2 || isempty(meanField)
        meanField = input('Enter the field name for mean (x) values: ', 's');
    end
    if nargin < 3 || isempty(varField)
        varField = input('Enter the field name for dependent values: ', 's');
    end
    if nargin < 4 || isempty(p0)
        p0 = input('Enter initial guesses for [Kd, n, A]: ');
    end
    if nargin < 5 || isempty(lb)
        lb = input('Enter lower bounds for [Kd, n, A]: ');
    end
    if nargin < 6 || isempty(ub)
        ub = input('Enter upper bounds for [Kd, n, A]: ');
    end

    %% Load Data
    fprintf('\nLoading data from %s...\n', dataFile);
    dataStruct = load(dataFile);
    
    % Try to get the independent variable from the saved structure.
    % If accessing the nested fields (using dot-separated string) fails,
    % default to retrieving the entire field with that name.
    try
        xdata = getFieldByPath(dataStruct, meanField);
    catch
        if isfield(dataStruct, meanField)
            xdata = dataStruct.(meanField);
        elseif evalin('base', ['exist(''', meanField, ''', ''var'')'])
            xdata = evalin('base', meanField);
        else
            error('Field or variable "%s" not found.', meanField);
        end
    end

    try
        ydata = getFieldByPath(dataStruct, varField);
    catch
        if isfield(dataStruct, varField)
            ydata = dataStruct.(varField);
        elseif evalin('base', ['exist(''', varField, ''', ''var'')'])
            ydata = evalin('base', varField);
        else
            error('Field or variable "%s" not found.', varField);
        end
    end
    
    %% Ensure that the data are column vectors    
    xdata = xdata(:);
    ydata = ydata(:);
    validIdx = (xdata >= 0) & ~isnan(xdata) & ~isnan(ydata);
    xdata = xdata(validIdx);
    ydata = ydata(validIdx);
    xdata = rmmissing(xdata);
    ydata = rmmissing(ydata);
    
    %% Define the Hill Function Model (without rounding n)
    % The Hill function is defined as:
    %   hill_function(x, Kd, n, A) = A * (x.^n) ./ (Kd.^n + x.^n)
    % Here n can take decimal values.
    %
    % p = [Kd, n, A]
    model = @(p, x) hill_function(x, p(1), p(2), p(3));

    %% Set Optimization Options
    opts = optimoptions('lsqcurvefit', ...
        'Display', 'iter', ...
        'Algorithm', 'levenberg-marquardt', ...
        'StepTolerance', 1e-100, ...
        'FunctionTolerance', 1e-100, ...
        'MaxIterations', 10000000, ...
        'MaxFunctionEvaluations', 50000000);


    %% Perform the Nonlinear Curve Fit
    fprintf('\nStarting lsqcurvefit for Hill function...\n');
    [p_fit, resnorm, residual, exitflag, output, lambda, jacobian] = ...
        lsqcurvefit(model, p0, xdata, ydata, lb, ub, opts);

    %% Compute Fit Statistics
    dof = length(ydata) - length(p_fit);
    mse = resnorm / dof;
    if ~isempty(jacobian)
        covar = inv(jacobian' * jacobian) * mse;
        param_se = full(sqrt(diag(covar)));
    else
        covar = [];
        param_se = [];
    end
    ss_tot = sum((ydata - mean(ydata)).^2);
    R_squared = 1 - resnorm / ss_tot;

    %% Create Fit Report Structure
    fitReport.exitFlag         = exitflag;
    fitReport.funcCount        = output.funcCount;
    fitReport.iterations       = output.iterations;
    fitReport.resnorm          = resnorm;
    fitReport.degreesOfFreedom = dof;
    fitReport.mse              = mse;
    fitReport.R_squared        = R_squared;
    fitReport.parameters       = struct('Kd', p_fit(1), 'n', p_fit(2), 'A', p_fit(3));
    if ~isempty(param_se)
        fitReport.paramSE = struct('Kd', param_se(1), 'n', param_se(2), 'A', param_se(3));
    else
        fitReport.paramSE = [];
    end
    fitReport.covariance = covar;

    %% Display Detailed Fit Report
    fprintf('\n----- Fit Report -----\n');
    fprintf('Exit Flag: %d\n', fitReport.exitFlag);
    fprintf('Function Evaluations: %d\n', fitReport.funcCount);
    fprintf('Iterations: %d\n', fitReport.iterations);
    fprintf('Residual Norm: %.4e\n', fitReport.resnorm);
    fprintf('Degrees of Freedom: %d\n', fitReport.degreesOfFreedom);
    fprintf('MSE: %.4e\n', fitReport.mse);
    fprintf('R-squared: %.4f\n\n', fitReport.R_squared);
    fprintf('Fitted Parameters:\n');
    fprintf('Kd = %.4e', fitReport.parameters.Kd);
    if ~isempty(fitReport.paramSE)
        fprintf(' ± %.4e\n', full(fitReport.paramSE.Kd));
    else
        fprintf('\n');
    end
    fprintf('n = %.4f', fitReport.parameters.n);
    if ~isempty(fitReport.paramSE)
        fprintf(' ± %.4f\n', full(fitReport.paramSE.n));
    else
        fprintf('\n');
    end
    fprintf('O = %.4f', fitReport.parameters.A);
    if ~isempty(fitReport.paramSE)
        fprintf(' ± %.4f\n', full(fitReport.paramSE.A));
    else
        fprintf('\n');
    end
    if ~isempty(covar)
        fprintf('\nCovariance Matrix:\n');
        disp(covar);
    end

    %% Plot the Fit Result
    % Bin the data (for plotting error bars)
    num_bins = 50;
    edges = quantile(xdata, linspace(0, 1, num_bins+1));
    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
    binned_means = zeros(1, num_bins);
    binned_stds = zeros(1, num_bins); % Now will hold standard error for xdata
    binned_y_means = zeros(1, num_bins);
    binned_y_stds = zeros(1, num_bins); % Now will hold standard error for ydata
    
    for i = 1:num_bins
    % For non-final bins, include left endpoint but exclude right endpoint.
    % For the final bin, include both endpoints.
    if i < num_bins
    bin_indices = xdata >= edges(i) & xdata < edges(i+1);
    else
    bin_indices = xdata >= edges(i) & xdata <= edges(i+1);
    end
    
    if any(bin_indices)
        % Calculate the number of data points in the current bin
        N = sum(bin_indices);
        
        % Compute the mean for xdata and ydata in the bin
        binned_means(i)   = mean(xdata(bin_indices));
        binned_y_means(i) = mean(ydata(bin_indices));
        
        % Compute standard error =
        % (standard deviation)/sqrt(n) if there is more than one point,
        % else set the error to zero.
        if N > 1
            binned_stds(i)   = std(xdata(bin_indices)) / sqrt(N);
            binned_y_stds(i) = std(ydata(bin_indices)) / sqrt(N);
        end
        else
            binned_stds(i)   = 0;
            binned_y_stds(i) = 0;
        end
    end
    


    figure;
    scatter(xdata* 1000 / (6.022e23 * (44e-9)^2 * sqrt(2 * 3.1415)), ydata, 100, 'MarkerEdgeColor', '#808080', 'MarkerFaceColor', '#808080');
    alpha(.5);
    hold on;
    errorbar(binned_means* 1000 / (6.022e23 * (44e-9)^2 * sqrt(2 * 3.1415)), binned_y_means, binned_y_stds, binned_y_stds, ...
             binned_stds* 1000 / (6.022e23 * (44e-9)^2 * sqrt(2 * 3.1415)), binned_stds* 1000 / (6.022e23 * (44e-9)^2 * sqrt(2 * 3.1415)), 'o', 'Color', 'k',...
             'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k',...
             'LineWidth', 3, 'MarkerSize', 15);
    x_fit = linspace(min(xdata), max(xdata), 1000);
    y_fit = model(p_fit, x_fit);
    plot(x_fit* 1000 / (6.022e23 * (44e-9)^2 * sqrt(2 * 3.1415)), y_fit, 'r-', 'LineWidth', 5);
    xlabel('Mean Concentration (M)', 'FontSize', 22);
    ylabel('Bound CI-mVenus monomers (#)', 'FontSize', 22);
    title('Binding curve for CI-mVenus', 'FontSize', 28);
    % Create a dynamic legend string for the fitted parameters
    fitString = sprintf(['Model fit $\\frac{O\\cdot\\langle x \\rangle^n}{K_d^n+\\langle x \\rangle^n}$:\n' ...
    '$K = %.3e \\pm %.3e M, \\; n = %.3f \\pm %.3f, \\; O = %.3f \\pm %.3f$'], ...
    fitReport.parameters.Kd * 1000 / (6.022e23 * (44e-9)^2 * sqrt(2 * 3.1415)), ...
    fitReport.paramSE.Kd * 1000 / (6.022e23 * (44e-9)^2 * sqrt(2 * 3.1415)), ...
    fitReport.parameters.n, fitReport.paramSE.n, ...
    fitReport.parameters.A, fitReport.paramSE.A);

    % Then, use the dynamic string in the legend call:
    legend('Data', 'Binned Data', fitString, 'FontSize', 22,'Interpreter','latex');
    ax = gca;
    ax.FontSize = 20;
    grid on;
    ax.GridLineWidth = 3;
    hold off;

end

%% Helper Function: hill_function
function f = hill_function(x, Kd, n, A)
% hill_function evaluates the Hill function:
%   f = A * (x.^n) ./ (Kd.^n + x.^n)
%
% Unlike previous versions, we do not round n so that decimal values are retained.
fval = A .* (x.^n) ./ (Kd.^n + x.^n);
% If tiny imaginary values occur due to complex-step derivative evaluations,
% discard them if they are negligible.
if max(abs(imag(fval))) < 1e-10
    f = real(fval);
else
    f = fval;
end
end

function value = getFieldByPath(dataStruct, fieldPath)
% getFieldByPath attempts to access nested fields in a structure
% using a dot-separated string. If an error occurs, it defaults to using
% the field directly.
try
    fields = strsplit(fieldPath, '.');
    value = dataStruct;
    for i = 1:length(fields)
        value = value.(fields{i});
    end
catch
    if isfield(dataStruct, fieldPath)
        value = dataStruct.(fieldPath);
    else
        error('Field "%s" not found in the data structure.', fieldPath);
    end
end
end

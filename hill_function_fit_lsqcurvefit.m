%% Modified MATLAB Code

function fitReport = fit_hill_model(dataFile, meanField, varField, p0, lb, ub)
% fit_hill_model performs nonlinear curve fitting using lsqcurvefit
% on a Hill function model.
%
% USAGE:
% fitReport = fit_hill_model(dataFile, meanField, varField, p0, lb, ub)
%
% INPUTS:
% dataFile - Name of the .mat file containing data.
% meanField - Field name for independent variable (x_mean).
% varField - Field name for the dependent variable (response).
% p0 - Initial guesses for [Kd, n, A].
% lb - Lower bounds for [Kd, n, A].
% ub - Upper bounds for [Kd, n, A].
%
% OUTPUT:
% fitReport - Structure containing the fit report details.

text
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
% default to retrieving the entire field with that name (from an already
% loaded vector).
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

%% Define the Hill Function Model
% The Hill function is defined as:
%   hill_function(x, Kd, n, A) = A * (x.^n) ./ (Kd.^n + x.^n)
% with n rounded to the nearest integer.
%
% p = [Kd, n, A]
model = @(p, x) hill_function(x, p(1), p(2), p(3));

%% Set Optimization Options
opts = optimoptions('lsqcurvefit', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt');

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
fprintf('n = %.4e', fitReport.parameters.n);
if ~isempty(fitReport.paramSE)
    fprintf(' ± %.4e\n', full(fitReport.paramSE.n));
else
    fprintf('\n');
end
fprintf('A = %.4e', fitReport.parameters.A);
if ~isempty(fitReport.paramSE)
    fprintf(' ± %.4e\n', full(fitReport.paramSE.A));
else
    fprintf('\n');
end
if ~isempty(covar)
    fprintf('\nCovariance Matrix:\n');
    disp(covar);
end

%% Plot the Fit Result
% Bin the data (for plotting error bars)
num_bins = 100;
edges = quantile(xdata, linspace(0, 1, num_bins+1));
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
binned_means = zeros(1, num_bins);
binned_stds = zeros(1, num_bins);
binned_y_means = zeros(1, num_bins);
binned_y_stds = zeros(1, num_bins);
for i = 1:num_bins
    if i < num_bins
        bin_indices = xdata >= edges(i) & xdata < edges(i+1);
    else
        bin_indices = xdata >= edges(i) & xdata <= edges(i+1);
    end
    if any(bin_indices)
        binned_means(i)   = mean(xdata(bin_indices));
        binned_stds(i)    = std(xdata(bin_indices));
        binned_y_means(i) = mean(ydata(bin_indices));
        binned_y_stds(i)  = std(ydata(bin_indices));
    end
end

figure;
scatter(xdata, ydata, 100, 'MarkerEdgeColor', '#edb232', 'MarkerFaceColor', '#edb232');
hold on;
errorbar(binned_means, binned_y_means, binned_y_stds, binned_y_stds, ...
         binned_stds, binned_stds, 'o', 'Color', 'k',...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k',...
         'LineWidth', 2, 'MarkerSize', 10);
x_fit = linspace(min(xdata), max(xdata), 1000);
y_fit = model(p_fit, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 5);
xlabel('Mean (x)', 'FontSize', 18);
ylabel('Response', 'FontSize', 18);
title('Hill Function Model Fit', 'FontSize', 24);
legend('Data', 'Binned Data', 'Fit', 'FontSize', 18);
ax = gca;
ax.FontSize = 16;
grid on;
ax.GridLineWidth = 3;
hold off;

end

%% Helper Function: hill_function
function f = hill_function(x, Kd, n, A)
% hill_function evaluates the Hill function:
% f = A * (x.^n) ./ (Kd.^n + x.^n)
f = A .* (x.^n) ./ (Kd.^n + x.^n);
end

%% Helper Function: getFieldByPath
function value = getFieldByPath(dataStruct, fieldPath)
% getFieldByPath attempts to access nested fields in a structure
% using a dot-separated string. If an error occurs (for example if one
% of the fields does not exist), it defaults to trying to access the field
% directly (i.e. without dot separation).
%
% Example:
% For fieldPath = 'group.subgroup.value' the function first tries:
% value = dataStruct.group.subgroup.value
% If an error is raised, it will then check if dataStruct has a field
% 'group.subgroup.value'.

text
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

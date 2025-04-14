% animate_simulation
% This script animates the binding simulation results, showing the cell grid and variance vs mean plot.

% % Run the simulation
% results = simple_binding2();

% Extract PSF-applied data
means = results.means;
vars = results.vars;
cell_grids = results.cell_grids_psf; % Dimensions: [n_N, n_pixels_x, n_pixels_y]

% Define Ns based on default parameters (0:1000)
n_N = length(means);
Ns = 0:(n_N-1); % Assumes N starts at 0 and increments by 1

% Determine color limits for consistent scaling
min_val = min(cell_grids(:));
max_val = max(cell_grids(:));

% Create figure with two subplots
figure;

% Top subplot: Simulated cell grid
subplot(2,1,1);
h_img = imagesc(squeeze(cell_grids(1,:,:))'); % Transpose for correct orientation
axis image;
clim([min_val, max_val]);
colorbar;
title(sprintf('Simulated Cell (N = %d)', Ns(1)));

% Bottom subplot: Variance vs Mean
subplot(2,1,2);
h_plot = plot(means(1), vars(1), 'bo', 'MarkerFaceColor', 'b');
xlabel('Mean Intensity');
ylabel('Variance');
title('Variance vs Mean Intensity');
grid on;
hold on;

% Set axis limits based on all data
xlim([min(means), max(means)]);
ylim([min(vars), max(vars)]);

% Animate each time step
for i = 1:n_N
    % Update cell image
    subplot(2,1,1);
    imagesc(squeeze(cell_grids(i,:,:))');
    title(sprintf('Simulated Cell (N = %d)', Ns(i)));
    axis image;
    caxis([min_val, max_val]);
    colormap("jet")
    % Update variance vs mean plot
    subplot(2,1,2);
    set(h_plot, 'XData', means(1:i), 'YData', vars(1:i));
    
    % Pause for animation
    drawnow;
    pause(0.01); % Adjust speed as needed
end
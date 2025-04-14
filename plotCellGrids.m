function plotCellGrids(data, numRows, numCols, cellStep, figTitle)
% plotSimulatedCellGrids Visualizes simulated cell grids in a grid layout.
% This version avoids an empty dummy subplot for the colorbar by creating a
% dedicated colorbar axes adjacent to the last image in each row.
%
% PARAMETERS:
%   data     - 3D array of size (numCells x height x width), where each slice
%              data(i,:,:) represents a single cell grid.
%   numRows  - (Optional) Number of rows in the subplot grid (default = 5).
%   numCols  - (Optional) Number of columns in the subplot grid (default = 5).
%   cellStep - (Optional) Step size between cells to plot (default = 40). For example,
%              if cellStep is 40, the function plots cell 1, then cell 41, etc.
%   figTitle - (Optional) Overall title for the figure (default = 'Simulated Cell Grids').
%
% EXAMPLE USAGE:
%   % For a 1000x45x23 dataset, to plot every 40th cell:
%   plotSimulatedCellGrids(data, 5, 5, 40, 'Cell Grid Visualization');

    % Set default parameters if not provided.
    if nargin < 5 || isempty(figTitle)
        figTitle = 'Simulated Cell Grids';
    end
    if nargin < 4 || isempty(cellStep)
        cellStep = 40;
    end
    if nargin < 3 || isempty(numCols)
        numCols = 5;
    end
    if nargin < 2 || isempty(numRows)
        numRows = 5;
    end

    % Calculate total number of images to plot.
    totalImages = numRows * numCols;
    if 1 + (totalImages - 1) * cellStep > size(data, 1)
        warning('Not enough cells in data for the requested grid and step. Adjusting total images.');
        totalImages = floor((size(data, 1) - 1) / cellStep) + 1;
        numCols = min(numCols, totalImages);
        numRows = ceil(totalImages / numCols);
    end

    % Create a new figure and set an overall title.
    figure;
    if exist('sgtitle', 'file')
        sgtitle(figTitle);
    end

    % Initialize a waitbar to show progress.
    hWait = waitbar(0, 'Plotting cell grids...');
    plotCounter = 0;

    % Loop over each row.
    for r = 1:numRows
        % -----------------------------------------------
        % First pass: determine common color limits for this row.
        % -----------------------------------------------
        rowMin = inf;
        rowMax = -inf;
        for c = 1:numCols
            idx = (r - 1) * numCols + c;     % Overall subplot index
            cellIndex = 1 + (idx - 1) * cellStep;  % Cell index in the data array
            if cellIndex > size(data, 1)
                break;
            end
            currImage = squeeze(data(cellIndex, :, :));
            rowMin = min(rowMin, min(currImage(:)));
            rowMax = max(rowMax, max(currImage(:)));
        end

        lastAx = [];
        % -----------------------------------------------
        % Second pass: create subplots for each cell in the row.
        % -----------------------------------------------
        for c = 1:numCols
            idx = (r - 1) * numCols + c;
            cellIndex = 1 + (idx - 1) * cellStep;
            if cellIndex > size(data, 1)
                break;
            end

            % Create subplot for the current cell.
            ax = subplot(numRows, numCols, (r - 1) * numCols + c);
            imagesc(squeeze(data(cellIndex, :, :)));
            axis image off;  % Clean up the axes display
            caxis([rowMin rowMax]);
            title(sprintf('Num of Proteins = %d', cellIndex));
            lastAx = ax;  % Store reference to the last subplot in this row

            % Update progress.
            plotCounter = plotCounter + 1;
            waitbar(plotCounter / totalImages, hWait, ...
                sprintf('Plotting image %d of %d...', plotCounter, totalImages));
        end

        % -----------------------------------------------
        % Create a dedicated colorbar for the row without a dummy plot.
        % -----------------------------------------------
        if ~isempty(lastAx)
            % Get the position of the last subplot in the current row.
            pos = get(lastAx, 'Position');
            % Define new axes for the colorbar next to the row.
            % Adjust 'cbGap' and 'cbWidth' as needed for aesthetics.
            cbGap = 0.01;
            cbWidth = 0.03;
            pos_cb = [pos(1) + pos(3) + cbGap, pos(2), cbWidth, pos(4)];
            ax_cb = axes('Position', pos_cb);
            set(ax_cb, 'Visible', 'off');  % Hide the axis ticks and box
            
            
            caxis([rowMin rowMax]);
            colorbar('peer', ax_cb, 'Location', 'eastoutside');
            ylabel(colorbar, 'Number', 'FontSize', 12);
        end
    end

    % Close the waitbar once plotting is complete.
    close(hWait);
end

CellArea = ci_yfp_fixed.binding.Area;

s = 'CI-YFP Fixed Binding';
sorted_areas = sort(CellArea);
n = length(sorted_areas);
percentages = zeros(n, 1);

% Precompute the maximum area to avoid repeated calls
max_area = sorted_areas(end);

% Calculate percentages for each A0
for i = 1:n
    A0 = sorted_areas(i);
    upper_bound = 2 * A0;
    % Use binary search to find the last index <= upper_bound
    low = i;
    high = n;
    j = i; % Initialize j to the smallest possible value
    while low <= high
        mid = floor((low + high) / 2);
        if sorted_areas(mid) <= upper_bound
            j = mid;
            low = mid + 1;
        else
            high = mid - 1;
        end
    end
    count = j - i + 1;
    percentages(i) = (count / n) * 100;
end

% Find the optimal A0
[max_percentage, idx] = max(percentages);
optimal_A0 = sorted_areas(idx);

% Plotting
figure;
subplot(2,1,1);
plot(sorted_areas, percentages, 'b-','LineWidth',4);
xlabel('A0','FontSize',24);
ylabel('Percentage (%)','FontSize',24);
title('Percentage of Cells Between A0 and 2A0','FontSize',26);
hold on;
xline(optimal_A0, 'r--', 'LineWidth', 3);
grid on;

subplot(2,1,2);
histogram(sorted_areas, 'BinMethod', 'auto','DisplayStyle','bar');
xlabel('Cell Area','FontSize',24);
ylabel('Count','FontSize',24);
title(['Cell Size Distribution for',varName],'FontSize',26);
hold on;
xline(optimal_A0, 'r--', 'LineWidth', 3);
xline(2.*double(optimal_A0),'r--','LineWidth',3)
grid on;
legend(['Optimal A0: ', num2str(optimal_A0),'um^2'],['Cells between A0 and 2A0: ', num2str(max_percentage), '%'],'FontSize',18)

% Display optimal A0
disp(['Optimal A0: ', num2str(optimal_A0)]);
disp(['Maximum Percentage: ', num2str(max_percentage), '%']);
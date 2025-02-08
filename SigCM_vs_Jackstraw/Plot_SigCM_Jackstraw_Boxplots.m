% Define rowNames and the data matrix
rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf', ...
    'Pt', 'De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce'};

%% Partitions_random_pvs
load('Partitions_random_pvs.mat');
maxLength_random_pvs = max(cellfun(@length, pvs_random));
% Initialize a new matrix to store the padded data
dataMatrix_random_pvs = NaN(maxLength_random_pvs, numel(pvs_random));
% Fill the data matrix, each column corresponds to a dataset
for i = 1:numel(pvs_random)
    data = pvs_random{i};
    dataMatrix_random_pvs(1:length(data), i) = data; % Fill each dataset into the corresponding column
end
%% Partitions_GT_pvs
load('Partitions_GT_pvs.mat');
maxLength_GT_pvs = max(cellfun(@length, pvs_GT));
% Initialize a new matrix to store the padded data
dataMatrix_GT_pvs = NaN(maxLength_GT_pvs, numel(pvs_GT));
% Fill the data matrix, each column corresponds to a dataset
for i = 1:numel(pvs_GT)
    data = pvs_GT{i};
    dataMatrix_GT_pvs(1:length(data), i) = data; % Fill each dataset into the corresponding column
end
%% Jackstraw_random_pvs
Jackstraw_pvs_random = cell(18,1);
for I = 1:18
    Jackstraw_pvs_random{I,1} = load(['Random_Jackstraw_all_pvs_', rowNames{I}, '.txt']);
end
maxLength_Jackstraw_random_pvs = max(cellfun(@length, Jackstraw_pvs_random));
% Initialize a new matrix to store the padded data
dataMatrix_Jackstraw_random_pvs = NaN(maxLength_Jackstraw_random_pvs, numel(Jackstraw_pvs_random));
% Fill the data matrix, each column corresponds to a dataset
for i = 1:numel(Jackstraw_pvs_random)
    data = Jackstraw_pvs_random{i};
    dataMatrix_Jackstraw_random_pvs(1:length(data), i) = data; % Fill each dataset into the corresponding column
end
%% Jackstraw_GT_pvs
Jackstraw_pvs_GT = cell(18,1);
for I = 1:18
    Jackstraw_pvs_GT{I,1} = load(['GT_Jackstraw_all_pvs_', rowNames{I}, '.txt']);
end
maxLength_Jackstraw_GT_pvs = max(cellfun(@length, Jackstraw_pvs_GT));
% Initialize a new matrix to store the padded data
dataMatrix_Jackstraw_GT_pvs = NaN(maxLength_Jackstraw_GT_pvs, numel(Jackstraw_pvs_GT));
% Fill the data matrix, each column corresponds to a dataset
for i = 1:numel(Jackstraw_pvs_GT)
    data = Jackstraw_pvs_GT{i};
    dataMatrix_Jackstraw_GT_pvs(1:length(data), i) = data; % Fill each dataset into the corresponding column
end

% Get screen resolution to automatically center the figure
screenSize = get(0, 'ScreenSize');  % Get screen size [x, y, width, height]
figWidth = 1600;  % Set the figure window width
figHeight = 1200;  % Set the figure window height

% Calculate the position to center the figure
xPos = (screenSize(3) - figWidth) / 2;
yPos = (screenSize(4) - figHeight) / 2;

% Create the figure window, centered on the screen
figure('Position', [xPos, yPos, figWidth, figHeight]);  % Use the calculated centered position

% Size of each subplot: 80mm x 123.5mm, convert to pixels (300 dpi)
width = 80 * 300 / 25.4;   % Width 80mm -> pixels
height = 123.5 * 300 / 25.4; % Height 123.5mm -> pixels
 
% Plot four boxplots
subplot(2, 2, 1);  % First plot
h1 = boxplot(dataMatrix_random_pvs, 'Labels', rowNames, 'Whisker', 1.5, 'Colors', 'k', ...
             'BoxStyle', 'outline', 'MedianStyle', 'target', 'Symbol', ' ', 'Widths', 0.5);
title('SigCM (Random partition)');
ylabel('p-value');
hold on;
yline(0.05, '--r', 'LineWidth', 1.5);  % Draw a horizontal dashed line at 0.05
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
set(gca, 'XTickLabelRotation', 0);  % Ensure x-axis labels are horizontal
set(gca, 'TickLength', [0.02, 0.02]);
% Adjust Position for elongating the width (increase the width value)
set(gca, 'Position', [0.05, 0.55, 0.4, 0.4]);

subplot(2, 2, 2);  % Second plot
h2 = boxplot(dataMatrix_GT_pvs, 'Labels', rowNames, 'Whisker', 1.5, 'Colors', 'k', ...
             'BoxStyle', 'outline', 'MedianStyle', 'target', 'Symbol', ' ', 'Widths', 0.5);
title('SigCM (Ground-truh partition)');
ylabel('p-value');
hold on;
yline(0.05, '--r', 'LineWidth', 1.5);
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
set(gca, 'XTickLabelRotation', 0);
set(gca, 'TickLength', [0.02, 0.02]);
% Adjust Position for elongating the width (increase the width value)
set(gca, 'Position', [0.55, 0.55, 0.4, 0.4]);

subplot(2, 2, 3);  % Third plot
h3 = boxplot(dataMatrix_Jackstraw_random_pvs, 'Labels', rowNames, 'Whisker', 1.5, 'Colors', 'k', ...
             'BoxStyle', 'outline', 'MedianStyle', 'target', 'Symbol', ' ', 'Widths', 0.5);
title('Jackstraw (Random partition)');
ylabel('p-value');
hold on;
yline(0.05, '--r', 'LineWidth', 1.5);
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
set(gca, 'XTickLabelRotation', 0);
set(gca, 'TickLength', [0.02, 0.02]);
% Adjust Position for elongating the width (increase the width value)
set(gca, 'Position', [0.05, 0.05, 0.4, 0.4]);

subplot(2, 2, 4);  % Fourth plot
h4 = boxplot(dataMatrix_Jackstraw_GT_pvs, 'Labels', rowNames, 'Whisker', 1.5, 'Colors', 'k', ...
             'BoxStyle', 'outline', 'MedianStyle', 'target', 'Symbol', ' ', 'Widths', 0.5);
title('Jackstraw (Ground-truh partition)');
ylabel('p-value');
hold on;
yline(0.05, '--r', 'LineWidth', 1.5);
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
set(gca, 'XTickLabelRotation', 0);
set(gca, 'TickLength', [0.02, 0.02]);
% Adjust Position for elongating the width (increase the width value)
set(gca, 'Position', [0.55, 0.05, 0.4, 0.4]);

% Adjust the layout of the figure to avoid overlap
set(gcf, 'Position', [xPos, yPos, figWidth, figHeight]);  % Set the figure position to keep it centered

% Global adjustment of outlier transparency and color
outliers = findobj(gcf, 'Tag', 'Outliers');  % Find all outliers in the current figure
set(outliers, 'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor', [0.7, 0.7, 0.7], 'MarkerSize', 3);
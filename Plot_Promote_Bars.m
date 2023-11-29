function Plot_Promote_Bars(threeBars)
% Define horizontal axis labels
rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf', 'Pt', 'De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce'};

% Draw stacked bar chart and get handles
bh = bar(threeBars, 'stacked');

% Set color for each stack
set(bh(1), 'FaceColor', [0.7, 0.1, 0.1]); % Set the first stack to dark red
set(bh(2), 'FaceColor', [0.8, 0.8, 0.8]); % Set the second stack to light gray
set(bh(3), 'FaceColor', [0.1, 0.1, 0.7]); % Set the third stack to dark blue

% Customize the appearance of axes and labels
set(gca, 'XTick', 1:length(rowNames), 'XTickLabel', rowNames, 'FontSize', 10, 'FontName', 'Arial');
ylabel('%', 'FontSize', 12, 'FontName', 'Arial');
legend('Better', 'Same', 'Worse', 'Location', 'best');
% title('Cluster enhancement via outlier removal', 'FontSize', 14, 'FontName', 'Arial');
% title('Cluster refinement via reassignment', 'FontSize', 14, 'FontName', 'Arial');

% Add values to each sub-bar
for i = 1:size(threeBars, 1) % Iterate over each bar
    cumulativeHeight = 0; % Cumulative height
    for j = 1:size(threeBars, 2) % Iterate over each sub-bar
        height = threeBars(i, j);
        if height ~= 0 % Display values only for non-zero heights
            cumulativeHeight = cumulativeHeight + height; % Update cumulative height
            x = bh(j).XEndPoints(i);
            y = cumulativeHeight - (height / 2); % Calculate the Y coordinate of the sub-bar midpoint
            % Add integer value and percentage sign to the middle of the sub-bar, text rotated 90 degrees and bold
            text(x, y, [num2str(height, '%d'), '%'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'white', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 8);
        end
    end
end
end
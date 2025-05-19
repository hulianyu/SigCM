% Select data based on condition 
x = com_FSC_Het2Hom(com_CMI_Het2Hom > 0.5);
y = com_FSC_Het2Hom(com_CMI_Het2Hom < 0.5);

% Prepare data and group labels
data = [x; y];
groups = [ones(size(x)); 2 * ones(size(y))];

% Seaborn color palette (e.g., 'muted' or 'deep' are popular choices)
color_palette = [0.8, 0.3, 0.5; 0.3, 0.5, 0.8]; % Using a 'muted' color palette with blue and red tones

% Draw boxplot with customizations
figure;
boxplot(data, groups, 'BoxStyle', 'outline', 'Colors', color_palette, 'MedianStyle', 'target', 'Symbol', 'o', 'Widths', 0.5);

% Adjust plot appearance
set(gca, 'XTickLabel', {'CMI > 0.5', 'CMI < 0.5'}, 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
ylabel('FSC (Het2Hom)', 'FontSize', 14, 'FontWeight', 'bold');

% Tight layout
set(gca, 'LooseInset', get(gca, 'TightInset'));
set(gca, 'Box', 'on'); % Add box around the plot for better readability

% Test for significance using Mann-Whitney U test
[p,h] = ranksum(x, y, 'tail', 'right');

% Add significance text with exact p-value (no rounding)
ylim_vals = get(gca, 'YLim');
text(1.5, ylim_vals(2)*0.95, ...
    ['Mann-Whitney test (one-tailed): p = ', num2str(p)], ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');

% Improve boxplot spacing
set(gca, 'XTickLabelRotation', 0); % Ensure x labels are horizontal
set(gca, 'TickLength', [0.02, 0.02]); % Adjust tick lengths for cleaner look

% Make figure more compact (tighten layout) and adjust height
set(gcf, 'Position', [100, 100, 500, 600]); % Adjust figure size to make it taller

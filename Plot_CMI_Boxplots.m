% Select data based on condition
x = com_ACC_kmodes(com_CMI_kmodes > 0.5);
y = com_ACC_kmodes(com_CMI_kmodes < 0.5);
% Prepare data and group labels
data = [x; y];
groups = [ones(size(x)); 2 * ones(size(y))];
% Convert color codes
color1 = [194, 255, 184]/255; % #C2FFB8
color2 = [184, 238, 255]/255; % #B8EEFF
% Draw boxplot
boxplot(data, groups, 'BoxStyle', 'outline', 'Colors', [color1; color2], 'MedianStyle', 'target', 'Symbol', 'o');
set(gca, 'XTickLabel', {'CMI > 0.5', 'CMI < 0.5'}, 'FontSize', 12, 'FontName', 'Arial');
% Set axis labels and title
ylabel('ACC');
% Set margins
set(gca, 'LooseInset', get(gca, 'TightInset'));
% Test for significance
[p,h] = ranksum(x, y, 'tail','right');
ylim = get(gca, 'YLim');
text(1.5, ylim(2)*0.95, ['Mann-Whitney test (one-tailed): p = ' num2str(p)], 'HorizontalAlignment', 'center');

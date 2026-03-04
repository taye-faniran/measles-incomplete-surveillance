function Figure2_plot_unknown_detection
% Plot Figure 2 from the saved non-identifiability results

clear; clc;

baseFolder = fileparts(mfilename('fullpath'));
load(fullfile(baseFolder, 'Figure2_unknown_detection_results.mat'), 'XX');

pdfFile = fullfile(folderPath, 'Figure2_unknown_detection.pdf');

p_vals = [0.25 0.5 0.75 1.0];

Rhat = round(XX(:,:,1), 2);
khat = round(XX(:,:,2), 2);


Rlow  = flip(round(XX(:,:,3), 3));
Rhigh = flip(round(XX(:,:,4), 3));
klow  = flip(round(XX(:,:,5), 3));
khigh = flip(round(XX(:,:,6), 3));

[rows, cols] = size(Rlow);

topGray = [0.92 0.92 0.92];
bottomTeal = [0.00 0.55 0.50];
highlightYellow = [0.93 0.90 0.35];

fig = figure('Color','w', 'Position', [100 100 1100 800]);
tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

% Top-left: Estimated R_eff
nexttile
draw_text_panel(Rhat, rows, cols, 'Estimated R_{eff}', p_vals, topGray, 'black');

% Top-right: Estimated k
nexttile
draw_text_panel(khat, rows, cols, 'Estimated k', p_vals, topGray, 'black');

% Bottom-left: 95% CI for R_eff
nexttile
draw_ci_panel(Rlow, Rhigh, rows, cols, ...
    '95% Conf. Interval for R_{eff}', p_vals, bottomTeal, highlightYellow);

% Bottom-right: 95% CI for k
nexttile
draw_ci_panel(klow, khigh, rows, cols, ...
    '95% Conf. Interval for k', p_vals, bottomTeal, highlightYellow);

exportgraphics(fig, pdfFile, 'ContentType', 'vector');
fprintf('Saved combined Figure 2 PDF to:\n%s\n', pdfFile);

end


% -------------------------------------------------------------------------
function draw_text_panel(mat, rows, cols, panelTitle, tickVals, bgColor, textColor)

imagesc(ones(rows, cols));
ax = gca;
ax.Colormap = bgColor;
caxis([0 1]);

axis equal
axis([0.5 cols+0.5 0.5 rows+0.5]);
set(gca, 'YDir', 'reverse');

for r = 1:rows
    for c = 1:cols
        text(c, r, sprintf('%.2f', mat(r,c)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Color', textColor, ...
            'FontSize', 14);
    end
end

hold on
for i = 1:rows-1
    yline(i + 0.5, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
end
for i = 1:cols-1
    xline(i + 0.5, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
end
hold off

title(panelTitle, 'FontSize', 16, 'FontWeight', 'bold');
xlabel('p_p', 'FontSize', 16);
ylabel('p_a', 'FontSize', 16);

xticks(1:cols);
yticks(1:rows);
xticklabels(tickVals);
yticklabels(flip(tickVals));

set(gca, 'FontSize', 13, 'LineWidth', 1);
end


% -------------------------------------------------------------------------
function draw_ci_panel(lowerMat, upperMat, rows, cols, panelTitle, tickVals, bgColor, hiColor)

imagesc(ones(rows, cols));
ax = gca;
ax.Colormap = bgColor;
caxis([0 1]);

axis equal
axis([0.5 cols+0.5 0.5 rows+0.5]);
set(gca, 'YDir', 'normal');

hold on

% Highlight only the extreme cells (the yellow cells in the paper figure)
for r = 1:rows
    for c = 1:cols
        if upperMat(r,c) >= 999
            rectangle('Position', [c-0.5, r-0.5, 1, 1], ...
                      'FaceColor', hiColor, ...
                      'EdgeColor', 'none');
        end
    end
end

for r = 1:rows
    for c = 1:cols
        label_string = sprintf('%.3f\n%.3f', lowerMat(r,c), upperMat(r,c));
        text(c, r, label_string, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Color', 'white', ...
            'FontSize', 13);
    end
end

for i = 1:rows-1
    yline(i + 0.5, 'w-', 'LineWidth', 2);
end
for i = 1:cols-1
    xline(i + 0.5, 'w-', 'LineWidth', 2);
end
hold off

title(panelTitle, 'FontSize', 16, 'FontWeight', 'bold');
xlabel('p_p', 'FontSize', 16);
ylabel('p_a', 'FontSize', 16);

xticks(1:cols);
yticks(1:rows);
xticklabels(tickVals);
yticklabels(flip(tickVals));

set(gca, 'FontSize', 13, 'LineWidth', 1);
end
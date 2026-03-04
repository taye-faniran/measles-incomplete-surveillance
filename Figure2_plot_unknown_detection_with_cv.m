function Figure2_plot_unknown_detection_with_cv
% Plot Figure 2 results plus the extra CV panel

clear; clc;

baseFolder = fileparts(mfilename('fullpath'));
load(fullfile(baseFolder, 'Figure2_unknown_detection_results.mat'), 'XX');

% Cap very large CI values for display only
XX = min(XX, 1000);

pdfFile = fullfile(folderPath, 'Figure2_unknown_detection_with_cv.pdf');

p_vals = [0.25 0.5 0.75 1.0];

Rhat  = flip(round(XX(:,:,1), 2));
khat  = flip(round(XX(:,:,2), 2));
Rlow  = flip(round(XX(:,:,3), 3));
Rhigh = flip(round(XX(:,:,4), 3));
klow  = flip(round(XX(:,:,5), 3));
khigh = flip(round(XX(:,:,6), 3));
CVhat = flip(round(XX(:,:,7), 2));

[rows, cols] = size(Rlow);

topGray = [0.92 0.92 0.92];
bottomTeal = [0.00 0.55 0.50];
highlightYellow = [0.93 0.90 0.35];

fig = figure('Color','w', 'Position', [100 100 1200 1100]);
tiledlayout(3,2, 'TileSpacing','compact', 'Padding','compact');

% 1. Estimated R_eff
nexttile
draw_text_panel(Rhat, rows, cols, 'Estimated R_{eff}', p_vals, topGray, 'black');

% 2. 95% CI for R_eff
nexttile
draw_ci_panel(Rlow, Rhigh, rows, cols, ...
    '95% Conf. Interval for R_{eff}', p_vals, bottomTeal, highlightYellow);

% 3. Estimated k
nexttile
draw_text_panel(khat, rows, cols, 'Estimated k', p_vals, topGray, 'black');

% 4. 95% CI for k
nexttile
draw_ci_panel(klow, khigh, rows, cols, ...
    '95% Conf. Interval for k', p_vals, bottomTeal, highlightYellow);

% 5. Estimated CV
nexttile
draw_text_panel(CVhat, rows, cols, 'Estimated CV', p_vals, topGray, 'black');

% 6. Leave blank or use for note
nexttile
axis off
text(0.5, 0.5, 'Supplementary diagnostic panel set', ...
    'HorizontalAlignment', 'center', 'FontSize', 15);

exportgraphics(fig, pdfFile, 'ContentType', 'vector');
fprintf('Saved supplementary Figure 2 PDF to:\n%s\n', pdfFile);

end


% -------------------------------------------------------------------------
function draw_text_panel(mat, rows, cols, panelTitle, tickVals, bgColor, textColor)

imagesc(ones(rows, cols));
ax = gca;
ax.Colormap = bgColor;
caxis([0 1]);

axis equal
axis([0.5 cols+0.5 0.5 rows+0.5]);
set(gca, 'YDir', 'normal');

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
xlabel('p_p', 'FontSize', 15);
ylabel('p_a', 'FontSize', 15);

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
xlabel('p_p', 'FontSize', 15);
ylabel('p_a', 'FontSize', 15);

xticks(1:cols);
yticks(1:rows);
xticklabels(tickVals);
yticklabels(flip(tickVals));

set(gca, 'FontSize', 13, 'LineWidth', 1);
end
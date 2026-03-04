clear; clc;

% -------------------------------------------------------------------------
% Load saved Figure 1 results
% -------------------------------------------------------------------------
outFile = fullfile(fileparts(mfilename('fullpath')), 'Figure1_known_detection_results.mat');

% Detection grid labels
p_vals = [0.25 0.5 0.75 1.0];
xvalues = p_vals;
ylabels = flip(p_vals);

% Extract pages from XX
Rhat   = flip(round(XX(:,:,1), 2));
khat   = flip(round(XX(:,:,2), 2));
Rlow   = flip(round(XX(:,:,3), 3));
Rhigh  = flip(round(XX(:,:,4), 3));
klow   = flip(round(XX(:,:,5), 3));
khigh  = flip(round(XX(:,:,6), 3));

[rows, cols] = size(Rlow);

topGray = [0.92 0.92 0.92];
bottomTeal = [0.00 0.55 0.50];


fig = figure('Color','w', 'Position', [100 100 1100 800]);
t = tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

% -------------------------------------------------------------------------
% Panel 1: Estimated R_eff
% -------------------------------------------------------------------------
nexttile
draw_text_panel(Rhat, rows, cols, ...
    'Estimated R_{eff}', p_vals, topGray, 'black');

% -------------------------------------------------------------------------
% Panel 2: Estimated k
% -------------------------------------------------------------------------
nexttile
draw_text_panel(khat, rows, cols, ...
    'Estimated k', p_vals, topGray, 'black');

% -------------------------------------------------------------------------
% Panel 3: 95% CI for R_eff
% -------------------------------------------------------------------------
nexttile
draw_ci_panel(Rlow, Rhigh, rows, cols, ...
    '95% Conf. Interval for R_{eff}', p_vals, bottomTeal);

% -------------------------------------------------------------------------
% Panel 4: 95% CI for k
% -------------------------------------------------------------------------
nexttile
draw_ci_panel(klow, khigh, rows, cols, ...
    '95% Conf. Interval for k', p_vals, bottomTeal);

% -------------------------------------------------------------------------
% Save combined figure automatically as PDF
% -------------------------------------------------------------------------
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
fprintf('Saved combined Figure 1 PDF to:\n%s\n', pdfFile);


% =========================================================================
% Helper: top panels with light gray background
% =========================================================================
function draw_text_panel(mat, rows, cols, panelTitle, tickVals, bgColor, textColor)

    bg = ones(rows, cols);
    imagesc(bg);

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


% =========================================================================
% Helper: bottom panels with teal background
% =========================================================================
function draw_ci_panel(lowerMat, upperMat, rows, cols, panelTitle, tickVals, bgColor)

    bg = ones(rows, cols);
    imagesc(bg);

    ax = gca;
    ax.Colormap = bgColor;
    caxis([0 1]);

    axis equal
    axis([0.5 cols+0.5 0.5 rows+0.5]);
    set(gca, 'YDir', 'normal');

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

    hold on
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
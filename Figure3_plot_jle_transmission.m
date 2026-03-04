function Figure3_plot_jle_transmission
% Plot Figure 3 from saved JLE transmission results.
% This script loads the saved JLE results matrix XX and creates a 2x2 figure:
%   (1) Estimated R_0
%   (2) Estimated k
%   (3) 95% Conf. Interval for R_0
%   (4) 95% Conf. Interval for k
% The figure is saved automatically as a PDF.

clear; clc; close all;

% -------------------------------------------------------------------------
% Paths
% -------------------------------------------------------------------------
baseFolder = fileparts(mfilename('fullpath'));
dataFile   = fullfile(baseFolder, 'Figure3_4_jle_results.mat');
outPdf     = fullfile(baseFolder, 'Figure3.pdf');% -------------------------------------------------------------------------
% Load saved results
% -------------------------------------------------------------------------
S = load(dataFile, 'XX');
XX = S.XX;

% Display order: p_a from 1.00 at top to 0.25 at bottom
xLabels = {'0.25','0.50','0.75','1.00'};
yLabels = {'1.00','0.75','0.50','0.25'};

% Figure 3 uses pages 1-6
R0_est  = flipud(round(XX(:,:,1), 2));
R0_low  = flipud(round(XX(:,:,2), 2));
R0_high = flipud(round(XX(:,:,3), 2));

k_est   = flipud(round(XX(:,:,4), 2));
k_low   = flipud(round(XX(:,:,5), 2));
k_high  = flipud(round(XX(:,:,6), 2));

topGray    = [0.90 0.90 0.90];
bottomTeal = [0.18 0.63 0.63];
gridDark   = [0.35 0.35 0.35];
txtDark    = [0.10 0.10 0.10];
txtWhite   = [1.00 1.00 1.00];

fig = figure('Color','w','Position',[80 80 1050 820]);

% ===================== Panel 1: Estimated R_0 ============================
ax1 = subplot(2,2,1);
draw_value_panel(ax1, R0_est, xLabels, yLabels, topGray, gridDark, txtDark, 'Estimated R_0', '%.2f');

% ===================== Panel 2: Estimated k ==============================
ax2 = subplot(2,2,2);
draw_value_panel(ax2, k_est, xLabels, yLabels, topGray, gridDark, txtDark, 'Estimated k', '%.2f');

% ===================== Panel 3: 95% CI for R_0 ===========================
ax3 = subplot(2,2,3);
draw_ci_panel(ax3, R0_low, R0_high, xLabels, yLabels, bottomTeal, txtWhite, '95% Conf. Interval for R_0');

% ===================== Panel 4: 95% CI for k =============================
ax4 = subplot(2,2,4);
draw_ci_panel(ax4, k_low, k_high, xLabels, yLabels, bottomTeal, txtWhite, '95% Conf. Interval for k');

% Improve spacing
set(fig, 'PaperPositionMode', 'auto');

% Save automatically to PDF
print(fig, outPdf, '-dpdf', '-bestfit');

fprintf('Figure saved to:\n%s\n', outPdf);

end



function draw_value_panel(ax, M, xLabels, yLabels, bgColor, gridColor, txtColor, ttl, fmt)

axes(ax); %#ok<LAXES>
hold(ax, 'on');

imagesc(ax, ones(size(M)));
colormap(ax, repmat(bgColor, 2, 1));
caxis(ax, [0 1]);

axis(ax, 'square');
set(ax, 'XTick', 1:4, 'XTickLabel', xLabels, ...
        'YTick', 1:4, 'YTickLabel', yLabels, ...
        'FontSize', 11, 'LineWidth', 1.0, ...
        'TickLength', [0 0], ...
        'XColor', txtColor, 'YColor', txtColor, ...
        'Box', 'on');

xlabel(ax, 'p_p', 'FontSize', 11, 'Color', txtColor);
ylabel(ax, 'p_a', 'FontSize', 11, 'Color', txtColor);
title(ax, ttl, 'FontSize', 13, 'FontWeight', 'normal', 'Color', txtColor);

% Grid lines
for i = 0.5:1:4.5
    plot(ax, [0.5 4.5], [i i], '-', 'Color', gridColor, 'LineWidth', 1.0);
    plot(ax, [i i], [0.5 4.5], '-', 'Color', gridColor, 'LineWidth', 1.0);
end

% Cell text
for r = 1:4
    for c = 1:4
        text(ax, c, r, sprintf(fmt, M(r,c)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 11, ...
            'FontWeight', 'normal', ...
            'Color', txtColor);
    end
end

xlim(ax, [0.5 4.5]);
ylim(ax, [0.5 4.5]);
set(ax, 'YDir', 'normal');
hold(ax, 'off');

end


% =========================================================================
% Draw bottom-row panels: darker teal background with white text
% =========================================================================
function draw_ci_panel(ax, L, U, xLabels, yLabels, bgColor, txtColor, ttl)

axes(ax); %#ok<LAXES>
hold(ax, 'on');

imagesc(ax, ones(size(L)));
colormap(ax, repmat(bgColor, 2, 1));
caxis(ax, [0 1]);

axis(ax, 'square');
set(ax, 'XTick', 1:4, 'XTickLabel', xLabels, ...
        'YTick', 1:4, 'YTickLabel', yLabels, ...
        'FontSize', 11, 'LineWidth', 1.0, ...
        'TickLength', [0 0], ...
        'XColor', txtColor, 'YColor', txtColor, ...
        'Box', 'on');

xlabel(ax, 'p_p', 'FontSize', 11, 'Color', txtColor);
ylabel(ax, 'p_a', 'FontSize', 11, 'Color', txtColor);
title(ax, ttl, 'FontSize', 13, 'FontWeight', 'normal', 'Color', txtColor);

% White grid lines
for i = 0.5:1:4.5
    plot(ax, [0.5 4.5], [i i], 'w-', 'LineWidth', 1.2);
    plot(ax, [i i], [0.5 4.5], 'w-', 'LineWidth', 1.2);
end

% CI text
for r = 1:4
    for c = 1:4
        label = sprintf('%.2f\n%.2f', L(r,c), U(r,c));
        text(ax, c, r, label, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 10.5, ...
            'FontWeight', 'normal', ...
            'Color', txtColor);
    end
end

xlim(ax, [0.5 4.5]);
ylim(ax, [0.5 4.5]);
set(ax, 'YDir', 'normal');
hold(ax, 'off');

end
function Figure4_plot_jle_detection
% Plot Figure 4 from saved JLE detection results.
% This script loads the saved JLE results matrix XX and creates a 2x2 figure:
%   (1) Estimated p_p
%   (2) Estimated p_a
%   (3) 95% Conf. Interval for p_p
%   (4) 95% Conf. Interval for p_a
% The figure is saved automatically as a PDF.

clear; clc; close all;

% -------------------------------------------------------------------------
% Paths
% -------------------------------------------------------------------------
baseFolder = fileparts(mfilename('fullpath'));
dataFile   = fullfile(baseFolder, 'Figure3_4_jle_results.mat');
outPdf     = fullfile(baseFolder, 'Figure4.pdf');
% -------------------------------------------------------------------------
% Load saved results
% -------------------------------------------------------------------------
S = load(dataFile, 'XX');
XX = S.XX;

% Display order: p_a from 1.00 at top to 0.25 at bottom
xLabels = {'0.25','0.50','0.75','1.00'};
yLabels = {'1.00','0.75','0.50','0.25'};

% Figure 4 uses pages 7-12
pp_est  = flipud(round(XX(:,:,7), 2));
pp_low  = flipud(round(XX(:,:,8), 2));
pp_high = flipud(round(XX(:,:,9), 2));

pa_est  = flipud(round(XX(:,:,10), 2));
pa_low  = flipud(round(XX(:,:,11), 2));
pa_high = flipud(round(XX(:,:,12), 2));

topGray    = [0.90 0.90 0.90];
bottomTeal = [0.18 0.63 0.63];
gridDark   = [0.35 0.35 0.35];
txtDark    = [0.10 0.10 0.10];
txtWhite   = [1.00 1.00 1.00];

fig = figure('Color','w','Position',[80 80 1050 820]);

% ===================== Panel 1: Estimated p_p ============================
ax1 = subplot(2,2,1);
draw_value_panel(ax1, pp_est, xLabels, yLabels, topGray, gridDark, txtDark, 'Estimated p_p', '%.2f');

% ===================== Panel 2: Estimated p_a ============================
ax2 = subplot(2,2,2);
draw_value_panel(ax2, pa_est, xLabels, yLabels, topGray, gridDark, txtDark, 'Estimated p_a', '%.2f');

% ===================== Panel 3: 95% CI for p_p ===========================
ax3 = subplot(2,2,3);
draw_ci_panel(ax3, pp_low, pp_high, xLabels, yLabels, bottomTeal, txtWhite, '95% Conf. Interval for p_p');

% ===================== Panel 4: 95% CI for p_a ===========================
ax4 = subplot(2,2,4);
draw_ci_panel(ax4, pa_low, pa_high, xLabels, yLabels, bottomTeal, txtWhite, '95% Conf. Interval for p_a');

set(fig, 'PaperPositionMode', 'auto');
print(fig, outPdf, '-dpdf', '-bestfit');

fprintf('Figure saved to:\n%s\n', outPdf);

end


% =========================================================================
% Draw top-row panels: light gray background with dark text/grid
% =========================================================================
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

for i = 0.5:1:4.5
    plot(ax, [0.5 4.5], [i i], '-', 'Color', gridColor, 'LineWidth', 1.0);
    plot(ax, [i i], [0.5 4.5], '-', 'Color', gridColor, 'LineWidth', 1.0);
end

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

for i = 0.5:1:4.5
    plot(ax, [0.5 4.5], [i i], 'w-', 'LineWidth', 1.2);
    plot(ax, [i i], [0.5 4.5], 'w-', 'LineWidth', 1.2);
end

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
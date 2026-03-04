function Figure5_plot_map_transmission
% Plot Figure 5 from saved MAP results.
% Figure 5 contains:
%   (1) Estimated R_eff
%   (2) Estimated k
%   (3) 95% Conf. Interval for R_eff
%   (4) 95% Conf. Interval for k

clear; clc; close all;

baseFolder = fileparts(mfilename('fullpath'));
dataFile   = fullfile(baseFolder, 'Figure5_6_map_results.mat');
outPdf     = fullfile(baseFolder, 'Figure5.pdf');
S = load(dataFile, 'XX');
XX = S.XX;

xLabels = {'0.25','0.50','0.75','1.00'};
yLabels = {'1.00','0.75','0.50','0.25'};

% Figure 5 uses pages 1-6 only
R_eff_est  = flipud(round(min(XX(:,:,1),1000), 2));
R_eff_low  = flipud(round(min(XX(:,:,2),1000), 2));
R_eff_high = flipud(round(min(XX(:,:,3),1000), 2));

k_est      = flipud(round(min(XX(:,:,4),1000), 2));
k_low      = flipud(round(min(XX(:,:,5),1000), 2));
k_high     = flipud(round(min(XX(:,:,6),1000), 2));

topGray    = [0.90 0.90 0.90];
bottomTeal = [0.18 0.63 0.63];
gridDark   = [0.35 0.35 0.35];
txtDark    = [0.10 0.10 0.10];
txtWhite   = [1.00 1.00 1.00];

fig = figure('Color','w','Position',[80 80 1050 820]);

ax1 = subplot(2,2,1);
draw_value_panel(ax1, R_eff_est, xLabels, yLabels, topGray, gridDark, txtDark, 'Estimated R_{eff}', '%.2f');

ax2 = subplot(2,2,2);
draw_value_panel(ax2, k_est, xLabels, yLabels, topGray, gridDark, txtDark, 'Estimated k', '%.2f');

ax3 = subplot(2,2,3);
draw_ci_panel(ax3, R_eff_low, R_eff_high, xLabels, yLabels, bottomTeal, txtWhite, '95% Conf. Interval for R_{eff}', '%.2f');

ax4 = subplot(2,2,4);
draw_ci_panel(ax4, k_low, k_high, xLabels, yLabels, bottomTeal, txtWhite, '95% Conf. Interval for k', '%.2f');

set(fig, 'PaperPositionMode', 'auto');
print(fig, outPdf, '-dpdf', '-bestfit');

fprintf('Figure saved to:\n%s\n', outPdf);

end

function draw_value_panel(ax, M, xLabels, yLabels, bgColor, gridColor, txtColor, ttl, fmt)
axes(ax);
hold(ax, 'on');

imagesc(ax, ones(size(M)));
colormap(ax, repmat(bgColor, 2, 1));
caxis(ax, [0 1]);

axis(ax, 'square');
set(ax, 'XTick', 1:4, 'XTickLabel', xLabels, ...
        'YTick', 1:4, 'YTickLabel', yLabels, ...
        'FontSize', 11, ...
        'LineWidth', 1.0, ...
        'TickLength', [0 0], ...
        'XColor', txtColor, ...
        'YColor', txtColor, ...
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
            'Color', txtColor);
    end
end

xlim(ax, [0.5 4.5]);
ylim(ax, [0.5 4.5]);
set(ax, 'YDir', 'normal');
hold(ax, 'off');
end

function draw_ci_panel(ax, L, U, xLabels, yLabels, bgColor, txtColor, ttl, fmt)
axes(ax);
hold(ax, 'on');

imagesc(ax, ones(size(L)));
colormap(ax, repmat(bgColor, 2, 1));
caxis(ax, [0 1]);

axis(ax, 'square');
set(ax, 'XTick', 1:4, 'XTickLabel', xLabels, ...
        'YTick', 1:4, 'YTickLabel', yLabels, ...
        'FontSize', 11, ...
        'LineWidth', 1.0, ...
        'TickLength', [0 0], ...
        'XColor', txtColor, ...
        'YColor', txtColor, ...
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
        text(ax, c, r, sprintf([fmt '\n' fmt], L(r,c), U(r,c)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 10.5, ...
            'Color', txtColor);
    end
end

xlim(ax, [0.5 4.5]);
ylim(ax, [0.5 4.5]);
set(ax, 'YDir', 'normal');
hold(ax, 'off');
end
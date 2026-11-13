function Figure3_JLE_16Panel_Correct
% FIGURE3_JLE_16PANEL_CORRECT
%
% SCIENTIFIC PURPOSE
% ------------------
% This program tests whether joint likelihood estimation (JLE) can recover
% transmission when surveillance is unknown but an auxiliary measurement
% of the combined observation probability is available.
%
% For each unique surveillance setting:
%   1. Generate 5,000 final observed clusters using
%          R_eff = 0.50, k = 0.15,
%          p_p in {0.25,0.50,0.75,1.00},
%          p_a in {1.00,0.75,0.50,0.25}.
%   2. Generate an auxiliary binomial surveillance measurement
%          Y ~ Binomial(N_aux,p_obs),
%          p_obs = p_p + (1-p_p)*p_a.
%   3. Estimate R_eff, k, p_p, and p_a jointly.
%   4. Profile over p_p and p_a when drawing the joint R_eff-k contour.
%
% The program does not repeat simulations until a favorable dataset is
% obtained. The random seed is fixed before one honest simulation run.
%
% FIGURE SYMBOLS
% --------------
%   black curve   : joint 2-D 95% profile-likelihood contour for (k,R_eff)
%   blue circle   : JLE maximum-likelihood estimate
%   red diamond   : generating transmission pair (0.15,0.50)
%   horizontal CI : separate 1-D 95% profile CI for k
%   vertical CI   : separate 1-D 95% profile CI for R_eff

clc;
close all;

%% 1. Main folder and automatically created output folder

programFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(programFolder);
baseFolder = fullfile(projectRoot,'outputs');

if ~isfolder(baseFolder)
    mkdir(baseFolder);
end

if ~isfolder(baseFolder)
    error('The main folder was not found:\n%s',baseFolder);
end

outputName = 'Figure3_JLE_16Panel_PublicationReady_v2';
outputFolder = fullfile(baseFolder,outputName);

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

logFile = fullfile(outputFolder,[outputName '_Command_Window_Log.txt']);
if isfile(logFile)
    delete(logFile);
end
diary(logFile);
diaryCleanup = onCleanup(@() safe_diary_off()); %#ok<NASGU>

startTimer = tic;

%% 2. Reproducible scientific settings

rng(20260817,'twister');

R_true = 0.50;
k_true = 0.15;

p_p_grid = [0.25 0.50 0.75 1.00];
p_a_grid = [1.00 0.75 0.50 0.25];

N_observed = 5000;
N_auxiliary = 5000;

% Numerical support.
Jmax = 2500;
Omax = 500;

maximumAllowedExpectedOmittedClusters = 0.01;

% Parameter bounds.
% p_p must remain strictly positive because the likelihood conditions on
% a cluster being observed. At p_p=0, that conditional likelihood is
% undefined. p_a is permitted to equal 0 exactly.
bounds.lower = [0.02, log(0.002), 0.0001, 0.0];
bounds.upper = [0.95, log(5.000), 1.0, 1.0];

% Likelihood-ratio thresholds expressed as increases in NLL.
threshold1D = 0.5*3.841458820694124;
threshold2D = 0.5*5.991464547107979;
maximumAllowedContourThresholdError = 0.02;

% Direct contour settings.
numberContourAngles = 64;
contourBisectionIterations = 14;

% Matching cache makes subsequent runs fast.
cacheVersion = 'Figure3-JLE-publication-v3-2026-08-17';
cacheFile = fullfile(outputFolder,[outputName '_Analysis_Cache.mat']);
forceRecompute = false;

fprintf('============================================================\n');
fprintf('STARTING CORRECT FIGURE 3: JOINT LIKELIHOOD ESTIMATION\n');
fprintf('True transmission: R_eff = %.4f, k = %.4f\n',R_true,k_true);
fprintf('Observed clusters per unique surveillance setting = %d\n',N_observed);
fprintf('Auxiliary binomial trials per unique setting = %d\n',N_auxiliary);
fprintf('All four parameters R_eff, k, p_p, and p_a are estimated.\n');
fprintf('The R_eff-k contours profile over p_p and p_a.\n');
fprintf('The auxiliary Binomial measurement is independent of the cluster data.\n');
fprintf('Exact p_a=0 and probability-one likelihood boundaries are handled safely.\n');
fprintf('p_p is kept strictly positive because inference conditions on observation.\n');
fprintf('No favorable-dataset selection and no spline interpolation are used.\n\n');
fprintf('ALL FIGURES AND FILES WILL BE SAVED AUTOMATICALLY IN:\n');
fprintf('%s\n',outputFolder);
fprintf('============================================================\n\n');

%% 3. Load a matching cache or perform the complete calculation

loadedCache = false;

if ~forceRecompute && isfile(cacheFile)
    cacheInfo = load(cacheFile,'cacheVersion');
    if isfield(cacheInfo,'cacheVersion') && ...
            strcmp(cacheInfo.cacheVersion,cacheVersion)
        load(cacheFile);
        loadedCache = true;
        fprintf('A matching saved JLE analysis cache was found.\n');
        fprintf('The simulated data, estimates, CIs, and contours were loaded.\n\n');
    end
end

if ~loadedCache

    numberRows = numel(p_a_grid);
    numberColumns = numel(p_p_grid);
    numberPanels = numberRows*numberColumns;

    panelLabel = cell(numberPanels,1);
    truePP = zeros(numberPanels,1);
    truePA = zeros(numberPanels,1);
    truePobs = zeros(numberPanels,1);
    auxiliarySuccesses = zeros(numberPanels,1);

    R_hat = NaN(numberPanels,1);
    R_lower = NaN(numberPanels,1);
    R_upper = NaN(numberPanels,1);
    k_hat = NaN(numberPanels,1);
    k_lower = NaN(numberPanels,1);
    k_upper = NaN(numberPanels,1);
    pp_hat = NaN(numberPanels,1);
    pa_hat = NaN(numberPanels,1);
    pobs_hat = NaN(numberPanels,1);
    minimumNLL = NaN(numberPanels,1);

    R_lower_at_bound = false(numberPanels,1);
    R_upper_at_bound = false(numberPanels,1);
    k_lower_at_bound = false(numberPanels,1);
    k_upper_at_bound = false(numberPanels,1);
    contourTouchesParameterBoundary = false(numberPanels,1);
    maximumContourThresholdError = NaN(numberPanels,1);

    omittedTrueMass = NaN(numberPanels,1);
    omittedObservedPositiveMass = NaN(numberPanels,1);
    panelRuntimeSeconds = NaN(numberPanels,1);

    observedSizesCell = cell(numberPanels,1);
    observedFrequenciesCell = cell(numberPanels,1);
    dataCell = cell(numberPanels,1);
    thetaHatCell = cell(numberPanels,1);
    contourKCell = cell(numberPanels,1);
    contourRCell = cell(numberPanels,1);

    reusedPP1Panel = [];
    panel = 0;

    fprintf('Generating data and fitting the 16 displayed JLE panels...\n\n');

    for row = 1:numberRows
        for column = 1:numberColumns
            panel = panel+1;
            panelTimer = tic;

            p_p_true = p_p_grid(column);
            p_a_true = p_a_grid(row);

            panelLabel{panel} = sprintf('(%c)',char(96+panel));
            truePP(panel) = p_p_true;
            truePA(panel) = p_a_true;
            truePobs(panel) = p_p_true+(1-p_p_true)*p_a_true;

            % When p_p=1, p_a has no effect.
            if abs(p_p_true-1)<1e-12 && ~isempty(reusedPP1Panel)
                source = reusedPP1Panel;

                auxiliarySuccesses(panel) = auxiliarySuccesses(source);
                R_hat(panel) = R_hat(source);
                R_lower(panel) = R_lower(source);
                R_upper(panel) = R_upper(source);
                k_hat(panel) = k_hat(source);
                k_lower(panel) = k_lower(source);
                k_upper(panel) = k_upper(source);
                pp_hat(panel) = pp_hat(source);
                pa_hat(panel) = pa_hat(source);
                pobs_hat(panel) = pobs_hat(source);
                minimumNLL(panel) = minimumNLL(source);
                R_lower_at_bound(panel) = R_lower_at_bound(source);
                R_upper_at_bound(panel) = R_upper_at_bound(source);
                k_lower_at_bound(panel) = k_lower_at_bound(source);
                k_upper_at_bound(panel) = k_upper_at_bound(source);
                contourTouchesParameterBoundary(panel) = ...
                    contourTouchesParameterBoundary(source);
                maximumContourThresholdError(panel) = ...
                    maximumContourThresholdError(source);
                omittedTrueMass(panel) = omittedTrueMass(source);
                omittedObservedPositiveMass(panel) = ...
                    omittedObservedPositiveMass(source);
                observedSizesCell{panel} = observedSizesCell{source};
                observedFrequenciesCell{panel} = ...
                    observedFrequenciesCell{source};
                dataCell{panel} = dataCell{source};
                thetaHatCell{panel} = thetaHatCell{source};
                contourKCell{panel} = contourKCell{source};
                contourRCell{panel} = contourRCell{source};
                panelRuntimeSeconds(panel) = 0;

                fprintf('%s true p_p=%.2f, p_a=%.2f: ', ...
                    panelLabel{panel},p_p_true,p_a_true);
                fprintf('reused the identical p_p=1 data and JLE analysis.\n');
                continue;
            end

            [data,tailInfo] = generate_observed_JLE_data( ...
                R_true,k_true,p_p_true,p_a_true,N_observed, ...
                N_auxiliary,Jmax,Omax);

            auxiliarySuccesses(panel) = data.auxSuccesses;
            omittedTrueMass(panel) = tailInfo.omittedTrueMass;
            omittedObservedPositiveMass(panel) = ...
                tailInfo.omittedObservedPositiveMass;
            observedSizesCell{panel} = data.observedSizes;
            observedFrequenciesCell{panel} = data.counts;
            dataCell{panel} = data;

            [thetaHat,minNLL] = fit_JLE_multistart(data,bounds);
            thetaHatCell{panel} = thetaHat;

            R_hat(panel) = thetaHat(1);
            k_hat(panel) = exp(thetaHat(2));
            pp_hat(panel) = thetaHat(3);
            pa_hat(panel) = thetaHat(4);
            pobs_hat(panel) = pp_hat(panel)+ ...
                (1-pp_hat(panel))*pa_hat(panel);
            minimumNLL(panel) = minNLL;

            [Rci,Rflags] = profile_CI_one_parameter( ...
                thetaHat,minNLL,1,data,bounds,threshold1D);
            [logKci,Kflags] = profile_CI_one_parameter( ...
                thetaHat,minNLL,2,data,bounds,threshold1D);

            R_lower(panel) = Rci(1);
            R_upper(panel) = Rci(2);
            k_lower(panel) = exp(logKci(1));
            k_upper(panel) = exp(logKci(2));
            R_lower_at_bound(panel) = Rflags.lower;
            R_upper_at_bound(panel) = Rflags.upper;
            k_lower_at_bound(panel) = Kflags.lower;
            k_upper_at_bound(panel) = Kflags.upper;

            [contourKCell{panel},contourRCell{panel}, ...
                contourTouchesParameterBoundary(panel), ...
                maximumContourThresholdError(panel)] = ...
                profile_contour_R_k(thetaHat,minNLL,data,bounds, ...
                Rci,logKci,threshold2D,numberContourAngles, ...
                contourBisectionIterations);

            if abs(p_p_true-1)<1e-12
                reusedPP1Panel = panel;
            end

            panelRuntimeSeconds(panel) = toc(panelTimer);

            fprintf(['%s true p_p=%.2f, p_a=%.2f, Y=%d/%d: ' ...
                'R_eff=%.4f [%s,%s], k=%.5g [%s,%s], ' ...
                'p_p_hat=%.5f, p_a_hat=%.5f, %.1f s\n'], ...
                panelLabel{panel},p_p_true,p_a_true, ...
                data.auxSuccesses,data.Naux, ...
                R_hat(panel), ...
                interval_text(R_lower(panel), ...
                    R_lower_at_bound(panel),'lower'), ...
                interval_text(R_upper(panel), ...
                    R_upper_at_bound(panel),'upper'), ...
                k_hat(panel), ...
                interval_text(k_lower(panel), ...
                    k_lower_at_bound(panel),'lower'), ...
                interval_text(k_upper(panel), ...
                    k_upper_at_bound(panel),'upper'), ...
                pp_hat(panel),pa_hat(panel),panelRuntimeSeconds(panel));
        end
    end

    %% Scientific and numerical checks before saving cache

    assert(all(cellfun(@(x) sum(x), ...
        observedFrequenciesCell)==N_observed), ...
        'At least one displayed panel does not contain 5,000 clusters.');

    pp1Panels = find(abs(truePP-1)<1e-12);

    if max(abs(R_hat(pp1Panels)-R_hat(pp1Panels(1))))>1e-10 || ...
            max(abs(k_hat(pp1Panels)-k_hat(pp1Panels(1))))>1e-10
        error('The p_p=1 panels are not identical, but they must be.');
    end

    expectedOmittedTrueClusters = N_observed*omittedTrueMass;
    expectedOmittedPositiveClusters = ...
        N_observed*omittedObservedPositiveMass;

    if max(expectedOmittedTrueClusters) > ...
            maximumAllowedExpectedOmittedClusters || ...
            max(expectedOmittedPositiveClusters) > ...
            maximumAllowedExpectedOmittedClusters
        error(['Numerical support is too short: the expected omitted ' ...
            'cluster count exceeds %.4g. Increase Jmax or Omax.'], ...
            maximumAllowedExpectedOmittedClusters);
    end

    fprintf(['\nNumerical-tail check passed: expected omitted count ' ...
        '< %.4g cluster per 5,000-cluster dataset.\n'], ...
        maximumAllowedExpectedOmittedClusters);

    if max(maximumContourThresholdError) > ...
            maximumAllowedContourThresholdError
        error(['A joint contour missed its requested likelihood ' ...
            'threshold by more than %.4g NLL units. Increase ' ...
            'contourBisectionIterations before saving.'], ...
            maximumAllowedContourThresholdError);
    end

    fprintf(['Contour check passed: maximum threshold error ' ...
        '<= %.4g NLL units.\n'], ...
        maximumAllowedContourThresholdError);

    save(cacheFile,'cacheVersion','R_true','k_true', ...
        'p_p_grid','p_a_grid','N_observed','N_auxiliary', ...
        'Jmax','Omax','maximumAllowedExpectedOmittedClusters', ...
        'bounds','threshold1D','threshold2D', ...
        'maximumAllowedContourThresholdError', ...
        'numberContourAngles','contourBisectionIterations', ...
        'panelLabel','truePP','truePA','truePobs', ...
        'auxiliarySuccesses','R_hat','R_lower','R_upper', ...
        'k_hat','k_lower','k_upper','pp_hat','pa_hat','pobs_hat', ...
        'minimumNLL','R_lower_at_bound','R_upper_at_bound', ...
        'k_lower_at_bound','k_upper_at_bound', ...
        'contourTouchesParameterBoundary', ...
        'maximumContourThresholdError','omittedTrueMass', ...
        'omittedObservedPositiveMass','panelRuntimeSeconds', ...
        'observedSizesCell','observedFrequenciesCell','dataCell', ...
        'thetaHatCell','contourKCell','contourRCell','-v7.3');

    fprintf('The verified JLE calculation was saved in the reusable cache.\n\n');
end

%% 4. Recheck loaded or newly calculated results

numberPanels = numel(panelLabel);

assert(numberPanels==16, ...
    'The JLE figure must contain 16 panels.');

assert(all(isfinite(R_hat)) && all(isfinite(k_hat)), ...
    'At least one JLE transmission estimate is not finite.');

assert(all(R_lower<=R_hat & R_hat<=R_upper), ...
    'At least one R_eff estimate lies outside its own profile CI.');

assert(all(k_lower<=k_hat & k_hat<=k_upper), ...
    'At least one k estimate lies outside its own profile CI.');

assert(max(maximumContourThresholdError) <= ...
    maximumAllowedContourThresholdError, ...
    ['A joint contour is not sufficiently close to the requested ' ...
    'likelihood threshold. Increase contourBisectionIterations.']);

finiteKUpperExceeds1 = isfinite(k_upper) & k_upper>1;

truthCoveredR = R_lower<=R_true & R_true<=R_upper;
truthCoveredK = k_lower<=k_true & k_true<=k_upper;
truthCoveredBoth = truthCoveredR & truthCoveredK;

truthCoveredJoint = false(numberPanels,1);

for panel = 1:numberPanels
    truthCoveredJoint(panel) = inpolygon(k_true,R_true, ...
        contourKCell{panel},contourRCell{panel});
end

expectedOmittedTrueClusters = N_observed*omittedTrueMass;
expectedOmittedPositiveClusters = ...
    N_observed*omittedObservedPositiveMass;

assert(max(expectedOmittedTrueClusters) <= ...
    maximumAllowedExpectedOmittedClusters, ...
    'The expected omitted true-cluster count exceeds its tolerance.');

assert(max(expectedOmittedPositiveClusters) <= ...
    maximumAllowedExpectedOmittedClusters, ...
    'The expected omitted observed-cluster count exceeds its tolerance.');

%% 5. Choose common k-axis

finitePositiveK = [k_true;k_hat; ...
    k_lower(isfinite(k_lower) & k_lower>0); ...
    k_upper(isfinite(k_upper) & k_upper>0)];

smallestReportedK = min(finitePositiveK);
largestReportedK = max(finitePositiveK);
reportedKRangeFactor = largestReportedK/smallestReportedK;

useLogK = reportedKRangeFactor>=100 || any(k_upper_at_bound);

[commonKLimits,commonRLimits] = common_axis_limits( ...
    k_true,R_true,k_hat,k_lower,k_upper,R_hat,R_lower,R_upper, ...
    contourKCell,contourRCell,bounds,useLogK);

%% 6. Save numerical CSV files

resultsTable = table(panelLabel,truePP,truePA,truePobs, ...
    repmat(N_observed,numberPanels,1), ...
    repmat(N_auxiliary,numberPanels,1),auxiliarySuccesses, ...
    R_hat,R_lower,R_upper,R_lower_at_bound,R_upper_at_bound, ...
    k_hat,k_lower,k_upper,k_lower_at_bound,k_upper_at_bound, ...
    pp_hat,pa_hat,pobs_hat,minimumNLL,truthCoveredR,truthCoveredK, ...
    truthCoveredJoint,contourTouchesParameterBoundary, ...
    maximumContourThresholdError,omittedTrueMass, ...
    omittedObservedPositiveMass,expectedOmittedTrueClusters, ...
    expectedOmittedPositiveClusters,panelRuntimeSeconds, ...
    'VariableNames',{'Panel','True_p_p','True_p_a','True_p_obs', ...
    'ObservedClusters','AuxiliaryTrials','AuxiliarySuccesses', ...
    'R_hat','R_lower_1D95','R_upper_1D95', ...
    'R_lower_reached_parameter_bound', ...
    'R_upper_reached_parameter_bound', ...
    'k_hat','k_lower_1D95','k_upper_1D95', ...
    'k_lower_reached_parameter_bound', ...
    'k_upper_reached_parameter_bound', ...
    'p_p_hat','p_a_hat','p_obs_hat','MinimumJointNLL', ...
    'True_R_covered','True_k_covered', ...
    'True_pair_covered_by_joint_2D_region', ...
    'JointContourTouchesParameterBoundary', ...
    'MaximumContourThresholdError','OmittedTrueMass', ...
    'OmittedObservedPositiveMass','ExpectedOmittedTrueClusters', ...
    'ExpectedOmittedPositiveClusters','RuntimeSeconds'});

resultsCSV = fullfile(outputFolder,[outputName '_results.csv']);
writetable(resultsTable,resultsCSV);

panelColumn = cell(0,1);
truePPColumn = zeros(0,1);
truePAColumn = zeros(0,1);
sizeColumn = zeros(0,1);
frequencyColumn = zeros(0,1);

for panel = 1:numberPanels
    q = numel(observedSizesCell{panel});

    panelColumn = [panelColumn; ...
        repmat(panelLabel(panel),q,1)]; %#ok<AGROW>

    truePPColumn = [truePPColumn; ...
        repmat(truePP(panel),q,1)]; %#ok<AGROW>

    truePAColumn = [truePAColumn; ...
        repmat(truePA(panel),q,1)]; %#ok<AGROW>

    sizeColumn = [sizeColumn; ...
        observedSizesCell{panel}(:)]; %#ok<AGROW>

    frequencyColumn = [frequencyColumn; ...
        observedFrequenciesCell{panel}(:)]; %#ok<AGROW>
end

observedDataTable = table(panelColumn,truePPColumn,truePAColumn, ...
    sizeColumn,frequencyColumn, ...
    'VariableNames',{'Panel','True_p_p','True_p_a', ...
    'ObservedClusterSize','Frequency'});

observedDataCSV = fullfile(outputFolder, ...
    [outputName '_observed_cluster_counts.csv']);

writetable(observedDataTable,observedDataCSV);

contourPanelColumn = cell(0,1);
contourPointColumn = zeros(0,1);
contourKColumn = zeros(0,1);
contourRColumn = zeros(0,1);

for panel = 1:numberPanels
    q = numel(contourKCell{panel});

    contourPanelColumn = [contourPanelColumn; ...
        repmat(panelLabel(panel),q,1)]; %#ok<AGROW>

    contourPointColumn = [contourPointColumn;(1:q)']; %#ok<AGROW>

    contourKColumn = [contourKColumn; ...
        contourKCell{panel}(:)]; %#ok<AGROW>

    contourRColumn = [contourRColumn; ...
        contourRCell{panel}(:)]; %#ok<AGROW>
end

contourTable = table(contourPanelColumn,contourPointColumn, ...
    contourKColumn,contourRColumn, ...
    'VariableNames',{'Panel','ContourPoint','k','R_eff'});

contourCSV = fullfile(outputFolder, ...
    [outputName '_profile_contours.csv']);

writetable(contourTable,contourCSV);

%% 7. Draw the readable 16-panel Figure 3

fig = figure('Color','w','Name',outputName, ...
    'Units','inches','Position',[0.4 0.4 16 9.5], ...
    'Renderer','painters','InvertHardcopy','on');

tl = tiledlayout(fig,4,4, ...
    'TileSpacing','compact','Padding','compact');

try
    tl.OuterPosition = [0.055 0.075 0.93 0.80];
catch
end

blue = [0.02 0.32 0.95];
red = [0.95 0.05 0.05];

for panel = 1:numberPanels
    ax = nexttile(tl,panel);
    hold(ax,'on');

    plot(ax,contourKCell{panel},contourRCell{panel}, ...
        'k-','LineWidth',2.0);

    draw_CI_cross(ax,k_hat(panel),R_hat(panel), ...
        k_lower(panel),k_upper(panel), ...
        R_lower(panel),R_upper(panel), ...
        commonKLimits,commonRLimits,blue);

    plot(ax,k_hat(panel),R_hat(panel),'o', ...
        'Color',blue,'MarkerFaceColor',blue, ...
        'MarkerSize',6.8,'LineWidth',1.0);

    plot(ax,k_true,R_true,'d', ...
        'Color',red,'MarkerFaceColor',[1.00 0.08 0.08], ...
        'MarkerSize',7.2,'LineWidth',1.0);

    xlim(ax,commonKLimits);
    ylim(ax,commonRLimits);

    if useLogK
        set(ax,'XScale','log');
    end

    title(ax,sprintf('%s  p_p=%.2f, p_a=%.2f', ...
        panelLabel{panel},truePP(panel),truePA(panel)), ...
        'FontSize',10,'FontWeight','bold','Interpreter','tex');

    box(ax,'on');
    grid(ax,'on');
    ax.GridAlpha = 0.12;
    ax.FontSize = 9.5;
    ax.LineWidth = 0.9;

    if mod(panel-1,4)~=0
        ax.YTickLabel = [];
    end

    if panel<=12
        ax.XTickLabel = [];
    end

    hold(ax,'off');
end

xlabel(tl,axis_label_k(useLogK), ...
    'FontSize',17,'FontWeight','bold');

ylabel(tl,'Effective reproduction number, R_{eff}', ...
    'FontSize',17,'FontWeight','bold','Interpreter','tex');

annotation(fig,'textbox',[0.05 0.955 0.90 0.035], ...
    'String',['JLE with auxiliary surveillance: transmission recovery ' ...
    'under 16 surveillance settings'], ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',18,'FontWeight','bold', ...
    'EdgeColor','none','Interpreter','tex');

annotation(fig,'textbox',[0.05 0.920 0.90 0.030], ...
    'String',sprintf(['Generating transmission: ' ...
    'R_{eff}=%.2f, k=%.2f; %d observed clusters and ' ...
    '%d independent auxiliary trials per unique setting'], ...
    R_true,k_true,N_observed,N_auxiliary), ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',11.5,'FontWeight','normal', ...
    'EdgeColor','none','Interpreter','tex');

legendAxis = axes('Parent',fig, ...
    'Position',[0.08 0.875 0.84 0.035], ...
    'Visible','off');

hold(legendAxis,'on');

h1 = plot(legendAxis,nan,nan,'k-','LineWidth',2.0);

h2 = plot(legendAxis,nan,nan,'o', ...
    'Color',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',6.8);

h3 = plot(legendAxis,nan,nan,'d', ...
    'Color',red,'MarkerFaceColor',[1 0.08 0.08], ...
    'MarkerSize',7.2);

h4 = plot(legendAxis,nan,nan,'--', ...
    'Color',blue,'LineWidth',1.7);

legend(legendAxis,[h1 h2 h3 h4], ...
    {'Joint 2-D 95% profile contour', ...
    'JLE estimate', ...
    'True transmission pair', ...
    'Separate 1-D 95% profile CIs'}, ...
    'Orientation','horizontal', ...
    'Location','north', ...
    'Box','off','FontSize',10.5);

hold(legendAxis,'off');

%% 8. Export PNG, PDF and MATLAB FIG

pngFile = fullfile(outputFolder,[outputName '.png']);
pdfFile = fullfile(outputFolder,[outputName '.pdf']);
figFile = fullfile(outputFolder,[outputName '.fig']);

fprintf('\nSaving the MATLAB FIG file...\n');
savefig(fig,figFile);

set(fig,'PaperUnits','inches');
set(fig,'PaperPosition',[0 0 16 9.5]);
set(fig,'PaperSize',[16 9.5]);
set(fig,'Color','w','InvertHardcopy','on');

fprintf('Saving the full-size PNG figure...\n');
print(fig,pngFile,'-dpng','-r300');

fprintf('Saving the vector PDF figure...\n');
print(fig,pdfFile,'-dpdf','-painters');

fprintf('PNG, PDF, and FIG files were saved successfully.\n');

%% 9. Save code copy and README

sourceFile = [mfilename('fullpath') '.m'];
savedCodeFile = fullfile(outputFolder,[mfilename '.m']);

if isfile(sourceFile) && ~strcmpi(sourceFile,savedCodeFile)
    copyfile(sourceFile,savedCodeFile,'f');
end

readmeFile = fullfile(outputFolder,'README_Figure3_JLE.txt');
fid = fopen(readmeFile,'w');

if fid~=-1
    fprintf(fid, ...
        'CORRECT FIGURE 3: JOINT LIKELIHOOD ESTIMATION\r\n\r\n');

    fprintf(fid, ...
        'Generating R_eff = %.4f and k = %.4f.\r\n', ...
        R_true,k_true);

    fprintf(fid, ...
        ['Each unique setting has %d final observed clusters and ' ...
        '%d independent auxiliary binomial trials.\r\n'], ...
        N_observed,N_auxiliary);

    fprintf(fid, ...
        ['The auxiliary model is Y~Binomial(N_aux,p_obs), where ' ...
        'p_obs=p_p+(1-p_p)*p_a.\r\n']);

    fprintf(fid, ...
        ['JLE estimates R_eff, k, p_p, and p_a jointly. ' ...
        'The black R_eff-k contour profiles over p_p and p_a.\r\n']);

    fprintf(fid, ...
        ['Thirteen independent datasets are used: one p_p=1 ' ...
        'dataset is reused across the four equivalent p_a panels.\r\n']);

    fprintf(fid, ...
        'Total independently simulated observed clusters = %d.\r\n', ...
        13*N_observed);

    fprintf(fid,'k-axis: %s.\r\n',axis_label_k(useLogK));

    fprintf(fid, ...
        ['No simulation was repeated to force favorable estimates, ' ...
        'and no spline interpolation was used.\r\n']);

    fprintf(fid, ...
        ['The exact p_a=0 and probability-one Binomial boundaries ' ...
        'are handled explicitly. p_p is constrained to be strictly ' ...
        'positive because the cluster likelihood is conditional on ' ...
        'a cluster being observed.\r\n']);

    fprintf(fid, ...
        ['The 1-D and 2-D 95%% profile regions use asymptotic ' ...
        'chi-square likelihood-ratio cutoffs with 1 and 2 degrees ' ...
        'of freedom, respectively. They are not ' ...
        'bootstrap-calibrated.\r\n']);

    fprintf(fid, ...
        ['The numerical support passed an expected-omitted-count ' ...
        'tolerance of %.4g cluster per dataset.\r\n'], ...
        maximumAllowedExpectedOmittedClusters);

    fprintf(fid, ...
        ['Every joint contour passed a maximum likelihood-threshold ' ...
        'error tolerance of %.4g NLL units.\r\n'], ...
        maximumAllowedContourThresholdError);

    fclose(fid);
end

%% 10. Final command-window report

fprintf('\n============================================================\n');
fprintf('DONE: THE CORRECT JLE FIGURE 3 WAS CREATED.\n');
fprintf('Elapsed time: %.1f seconds.\n',toc(startTimer));
fprintf('Unique independently simulated datasets = 13.\n');
fprintf('Observed clusters per unique setting = %d.\n',N_observed);

fprintf('Total independently simulated observed clusters = %d.\n', ...
    13*N_observed);

fprintf('Auxiliary trials per unique setting = %d.\n', ...
    N_auxiliary);

fprintf('R_eff estimates range from %.4f to %.4f.\n', ...
    min(R_hat),max(R_hat));

fprintf('k estimates range from %.5g to %.5g.\n', ...
    min(k_hat),max(k_hat));

fprintf('Panels whose R_eff CI contains the truth = %d of 16.\n', ...
    sum(truthCoveredR));

fprintf('Panels whose k CI contains the truth = %d of 16.\n', ...
    sum(truthCoveredK));

fprintf(['Panels whose two separate CIs both contain the truth ' ...
    '= %d of 16.\n'],sum(truthCoveredBoth));

fprintf(['Joint 2-D regions containing the true transmission pair ' ...
    '= %d of 16.\n'],sum(truthCoveredJoint));

fprintf('Finite k upper limits exceeding 1 = %d.\n', ...
    sum(finiteKUpperExceeds1));

fprintf('k upper intervals reaching the parameter bound = %d.\n', ...
    sum(k_upper_at_bound));

fprintf('Reported k range factor = %.3f.\n', ...
    reportedKRangeFactor);

fprintf('The selected k-axis is: %s.\n', ...
    axis_label_k(useLogK));

fprintf('Joint contours touching a parameter boundary = %d.\n', ...
    sum(contourTouchesParameterBoundary));

fprintf('Largest contour threshold error = %.5g NLL units.\n', ...
    max(maximumContourThresholdError));

fprintf('Allowed contour threshold error = %.5g NLL units.\n', ...
    maximumAllowedContourThresholdError);

fprintf('Largest omitted true-cluster tail mass = %.3g.\n', ...
    max(omittedTrueMass));

fprintf('Largest omitted positive-observation tail mass = %.3g.\n', ...
    max(omittedObservedPositiveMass));

fprintf('Largest expected omitted true-cluster count = %.3g.\n', ...
    max(expectedOmittedTrueClusters));

fprintf(['Largest expected omitted positive-observation count ' ...
    '= %.3g.\n'],max(expectedOmittedPositiveClusters));

fprintf('Allowed expected omitted count per dataset = %.3g.\n', ...
    maximumAllowedExpectedOmittedClusters);

fprintf('\nALL FIGURES AND FILES WERE SAVED SUCCESSFULLY IN:\n');
fprintf('%s\n',outputFolder);

fprintf('\nPNG FIGURE:\n%s\n',pngFile);
fprintf('\nPDF FIGURE:\n%s\n',pdfFile);
fprintf('\nMATLAB FIG FILE:\n%s\n',figFile);
fprintf('\nRESULTS CSV:\n%s\n',resultsCSV);
fprintf('\nOBSERVED-COUNTS CSV:\n%s\n',observedDataCSV);
fprintf('\nPROFILE-CONTOUR CSV:\n%s\n',contourCSV);
fprintf('\nREUSABLE MAT CACHE:\n%s\n',cacheFile);
fprintf('\nSAVED CODE COPY:\n%s\n',savedCodeFile);
fprintf('\nREADME:\n%s\n',readmeFile);
fprintf('\nCOMMAND-WINDOW LOG:\n%s\n',logFile);
fprintf('============================================================\n');

end

%% =====================================================================
function [data,tailInfo] = generate_observed_JLE_data( ...
    R,k,p_p,p_a,Nobserved,Naux,Jmax,Omax)

pS = true_cluster_pmf(R,k,Jmax);
tailInfo.omittedTrueMass = max(0,1-sum(pS));

K = observation_kernel_vectorized(p_p,p_a,Jmax,0:Omax);
pO = pS*K;
probabilityObserved = sum(pS)-pO(1);

if probabilityObserved<=0
    error('The probability of observing a cluster is not positive.');
end

pPositive = pO(2:end)/probabilityObserved;
tailInfo.omittedObservedPositiveMass = ...
    max(0,1-sum(pPositive));

pSimulation = max(pPositive,0);
pSimulation = pSimulation/sum(pSimulation);

countsFull = mnrnd(Nobserved,pSimulation);
lastObservedSize = find(countsFull>0,1,'last');

if isempty(lastObservedSize)
    error('No positive observed clusters were generated.');
end

countsFull = countsFull(1:lastObservedSize);
observedSizes = find(countsFull>0);
counts = countsFull(observedSizes);

p_obs = p_p+(1-p_p)*p_a;
auxSuccesses = binornd(Naux,p_obs);

data.observedSizes = observedSizes(:)';
data.counts = counts(:)';
data.Naux = Naux;
data.auxSuccesses = auxSuccesses;
data.Jmax = Jmax;
data.s = (1:Jmax)';
data.logChoose = precompute_log_choose( ...
    Jmax,data.observedSizes);

end

%% =====================================================================
function logChoose = precompute_log_choose(Jmax,observedSizes)

s = (1:Jmax)';
o = observedSizes(:)';
remaining = s-o;
valid = remaining>=0;

safeRemaining = max(remaining,0);

logChoose = gammaln(s+1)-gammaln(o+1) ...
    -gammaln(safeRemaining+1);

logChoose(~valid) = -Inf;

end

%% =====================================================================
function [thetaHat,minNLL] = fit_JLE_multistart(data,bounds)

RkStarts = [ ...
    0.25 0.04; ...
    0.35 0.15; ...
    0.50 0.15; ...
    0.65 0.08; ...
    0.75 0.50];

detectionStarts = [ ...
    0.15 0.15; ...
    0.15 0.75; ...
    0.50 0.50; ...
    0.75 0.15; ...
    0.75 0.75];

numberStarts = ...
    size(RkStarts,1)*size(detectionStarts,1);

starts = zeros(numberStarts,4);
q = 0;

for i = 1:size(RkStarts,1)
    for j = 1:size(detectionStarts,1)
        q = q+1;

        starts(q,:) = [ ...
            RkStarts(i,1), ...
            log(RkStarts(i,2)), ...
            detectionStarts(j,1), ...
            detectionStarts(j,2)];
    end
end

options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','off', ...
    'MaxFunctionEvaluations',1800, ...
    'MaxIterations',350, ...
    'OptimalityTolerance',1e-7, ...
    'StepTolerance',1e-9);

objective = @(theta) JLE_nll(theta,data,bounds);

minNLL = Inf;
thetaHat = NaN(1,4);

for i = 1:numberStarts
    x0 = min(max(starts(i,:),bounds.lower),bounds.upper);

    try
        [candidate,value] = fmincon( ...
            objective,x0,[],[],[],[], ...
            bounds.lower,bounds.upper,[],options);

        if isfinite(value) && value<minNLL
            minNLL = value;
            thetaHat = candidate(:)';
        end
    catch
    end
end

if any(~isfinite(thetaHat))
    error(['The four-parameter JLE optimization failed ' ...
        'from every starting value.']);
end

tightOptions = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','off', ...
    'MaxFunctionEvaluations',3500, ...
    'MaxIterations',600, ...
    'OptimalityTolerance',1e-9, ...
    'StepTolerance',1e-10);

[thetaHat,minNLL] = fmincon( ...
    objective,thetaHat,[],[],[],[], ...
    bounds.lower,bounds.upper,[],tightOptions);

thetaHat = thetaHat(:)';

end

%% =====================================================================
function value = JLE_nll(theta,data,bounds)

theta = theta(:)';

if numel(theta)~=4 || any(~isfinite(theta)) || ...
        any(theta<bounds.lower) || ...
        any(theta>bounds.upper)
    value = 1e100;
    return;
end

R = theta(1);
k = exp(theta(2));
p_p = theta(3);
p_a = theta(4);

pS = true_cluster_pmf(R,k,data.Jmax);

[positiveSelected,pZero] = ...
    selected_observation_probabilities( ...
    pS,p_p,p_a,data);

denominator = sum(pS)-pZero;

if denominator<=0 || ...
        any(positiveSelected<=0) || ...
        any(~isfinite(positiveSelected))
    value = 1e100;
    return;
end

pConditional = positiveSelected/denominator;

nllClusters = ...
    -sum(data.counts.*log(pConditional));

p_obs = p_p+(1-p_p)*p_a;
Y = data.auxSuccesses;
N = data.Naux;

logAuxiliary = ...
    exact_binomial_log_likelihood(Y,N,p_obs);

if ~isfinite(logAuxiliary)
    value = 1e100;
    return;
end

nllAuxiliary = -logAuxiliary;
value = nllClusters+nllAuxiliary;

if ~isfinite(value) || ~isreal(value)
    value = 1e100;
end

end

%% =====================================================================
function logLikelihood = ...
    exact_binomial_log_likelihood(Y,N,p)

if p<0 || p>1 || Y<0 || Y>N
    logLikelihood = -Inf;
    return;
end

if p==0
    if Y==0
        logLikelihood = 0;
    else
        logLikelihood = -Inf;
    end
    return;
end

if p==1
    if Y==N
        logLikelihood = 0;
    else
        logLikelihood = -Inf;
    end
    return;
end

logBinomialCoefficient = ...
    gammaln(N+1)-gammaln(Y+1)-gammaln(N-Y+1);

logLikelihood = logBinomialCoefficient;

if Y>0
    logLikelihood = logLikelihood+Y*log(p);
end

if N-Y>0
    logLikelihood = ...
        logLikelihood+(N-Y)*log1p(-p);
end

end

%% =====================================================================
function [positiveSelected,pZero] = ...
    selected_observation_probabilities( ...
    pS,p_p,p_a,data)

s = data.s;
o = data.observedSizes;
logChoose = data.logChoose;
remaining = s-o;

z = p_p+(1-p_p)*p_a;

if z<=0
    term1 = zeros(size(logChoose));
elseif z>=1
    term1 = double(s==o);
else
    logTerm1 = logChoose ...
        +o*log(z) ...
        +remaining*log1p(-z);

    term1 = exp(logTerm1);
end

if p_p>=1 || p_a<=0
    term2 = zeros(size(logChoose));
elseif p_a>=1
    term2 = (1-p_p).^s.*double(s==o);
else
    logTerm2 = ...
        s*log1p(-p_p) ...
        +logChoose ...
        +o*log(p_a) ...
        +remaining*log1p(-p_a);

    term2 = exp(logTerm2);
end

kernelPositive = max(term1-term2,0);

positiveSelected = pS*kernelPositive;

pZero = sum(pS(:).*((1-p_p).^s));

end

%% =====================================================================
function [CI,flags] = profile_CI_one_parameter( ...
    thetaHat,minNLL,targetIndex,data,bounds,threshold)

target = minNLL+threshold;
estimate = thetaHat(targetIndex);

[lower,flags.lower] = find_profile_crossing( ...
    thetaHat,targetIndex,estimate, ...
    bounds.lower(targetIndex), ...
    target,data,bounds);

[upper,flags.upper] = find_profile_crossing( ...
    thetaHat,targetIndex,estimate, ...
    bounds.upper(targetIndex), ...
    target,data,bounds);

CI = [lower upper];

end

%% =====================================================================
function [crossing,reachedBoundary] = ...
    find_profile_crossing( ...
    thetaHat,targetIndex,estimate,boundary, ...
    target,data,bounds)

screen = linspace(estimate,boundary,14);

previousX = screen(1);

previousDifference = ...
    profile_nll_fixed_parameter( ...
    previousX,targetIndex,thetaHat,data,bounds)-target;

bracketFound = false;

for i = 2:numel(screen)
    currentX = screen(i);

    currentDifference = ...
        profile_nll_fixed_parameter( ...
        currentX,targetIndex,thetaHat,data,bounds)-target;

    if currentDifference>=0 && previousDifference<0
        lowerBracket = min(previousX,currentX);
        upperBracket = max(previousX,currentX);
        bracketFound = true;
        break;
    end

    previousX = currentX;
    previousDifference = currentDifference;
end

if ~bracketFound
    crossing = boundary;
    reachedBoundary = true;
    return;
end

options = optimset( ...
    'Display','off','TolX',2e-5);

crossing = fzero( ...
    @(x) profile_nll_fixed_parameter( ...
    x,targetIndex,thetaHat,data,bounds)-target, ...
    [lowerBracket upperBracket],options);

reachedBoundary = false;

end

%% =====================================================================
function bestValue = ...
    profile_nll_fixed_parameter( ...
    fixedValue,targetIndex,thetaHat,data,bounds)

remaining = setdiff(1:4,targetIndex);

template = thetaHat;
template(targetIndex) = fixedValue;

lower = bounds.lower(remaining);
upper = bounds.upper(remaining);

start1 = min(max( ...
    thetaHat(remaining),lower),upper);

midpoint = 0.5*(lower+upper);

logKposition = find(remaining==2,1);

if ~isempty(logKposition)
    midpoint(logKposition) = log(0.15);
end

starts = [start1;midpoint];

options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','off', ...
    'MaxFunctionEvaluations',900, ...
    'MaxIterations',220, ...
    'OptimalityTolerance',2e-7, ...
    'StepTolerance',2e-8);

bestValue = Inf;

for i = 1:size(starts,1)
    objective = @(x) assemble_and_evaluate( ...
        x,remaining,template,data,bounds);

    try
        [~,value] = fmincon( ...
            objective,starts(i,:),[],[],[],[], ...
            lower,upper,[],options);

        bestValue = min(bestValue,value);
    catch
    end
end

if ~isfinite(bestValue)
    bestValue = 1e100;
end

end

%% =====================================================================
function value = assemble_and_evaluate( ...
    freeValues,freeIndices,template,data,bounds)

theta = template;
theta(freeIndices) = freeValues;

value = JLE_nll(theta,data,bounds);

end

%% =====================================================================
function [kCurve,RCurve, ...
    touchesBoundary,maxThresholdError] = ...
    profile_contour_R_k( ...
    thetaHat,minNLL,data,bounds, ...
    Rci,logKci,threshold, ...
    numberAngles,bisectionIterations)

Rhat = thetaHat(1);
logKhat = thetaHat(2);

scaleRlower = max(Rhat-Rci(1),0.005);
scaleRupper = max(Rci(2)-Rhat,0.005);
scaleKlower = max(logKhat-logKci(1),0.03);
scaleKupper = max(logKci(2)-logKhat,0.03);

angles = linspace(0,2*pi,numberAngles+1);

kCurve = NaN(size(angles));
RCurve = NaN(size(angles));
thresholdErrors = NaN(size(angles));

touchesBoundary = false;

for q = 1:numel(angles)
    directionK = cos(angles(q));
    directionR = sin(angles(q));

    if directionK>=0
        scaleK = scaleKupper;
    else
        scaleK = scaleKlower;
    end

    if directionR>=0
        scaleR = scaleRupper;
    else
        scaleR = scaleRlower;
    end

    maximumRadius = Inf;

    if directionK>1e-12
        maximumRadius = min( ...
            maximumRadius, ...
            (bounds.upper(2)-logKhat) ...
            /(directionK*scaleK));

    elseif directionK<-1e-12
        maximumRadius = min( ...
            maximumRadius, ...
            (bounds.lower(2)-logKhat) ...
            /(directionK*scaleK));
    end

    if directionR>1e-12
        maximumRadius = min( ...
            maximumRadius, ...
            (bounds.upper(1)-Rhat) ...
            /(directionR*scaleR));

    elseif directionR<-1e-12
        maximumRadius = min( ...
            maximumRadius, ...
            (bounds.lower(1)-Rhat) ...
            /(directionR*scaleR));
    end

    maximumRadius = ...
        max(0,0.999*maximumRadius);

    radiusLow = 0;
    radiusHigh = min(1.45,maximumRadius);
    deltaLow = 0;

    deltaHigh = profiled_R_k_delta( ...
        radiusHigh,directionK,directionR, ...
        scaleK,scaleR,thetaHat,minNLL,data,bounds);

    while deltaHigh<threshold && ...
            radiusHigh<maximumRadius

        newRadius = min( ...
            1.45*radiusHigh,maximumRadius);

        if newRadius<=radiusHigh+1e-10
            break;
        end

        radiusHigh = newRadius;

        deltaHigh = profiled_R_k_delta( ...
            radiusHigh,directionK,directionR, ...
            scaleK,scaleR,thetaHat,minNLL,data,bounds);
    end

    if deltaHigh<threshold
        radiusBoundary = radiusHigh;

        bestThresholdError = ...
            abs(deltaHigh-threshold);

        touchesBoundary = true;
    else
        bestRadius = radiusHigh;

        bestThresholdError = ...
            abs(deltaHigh-threshold);

        for iteration = 1:bisectionIterations
            radiusMiddle = ...
                0.5*(radiusLow+radiusHigh);

            deltaMiddle = profiled_R_k_delta( ...
                radiusMiddle,directionK,directionR, ...
                scaleK,scaleR,thetaHat,minNLL,data,bounds);

            currentError = ...
                abs(deltaMiddle-threshold);

            if currentError<bestThresholdError
                bestThresholdError = currentError;
                bestRadius = radiusMiddle;
            end

            if deltaMiddle<threshold
                radiusLow = radiusMiddle;
                deltaLow = deltaMiddle;
            else
                radiusHigh = radiusMiddle;
                deltaHigh = deltaMiddle;
            end
        end

        if abs(deltaLow-threshold)<bestThresholdError
            bestThresholdError = ...
                abs(deltaLow-threshold);

            bestRadius = radiusLow;
        end

        if abs(deltaHigh-threshold)<bestThresholdError
            bestThresholdError = ...
                abs(deltaHigh-threshold);

            bestRadius = radiusHigh;
        end

        radiusBoundary = bestRadius;
    end

    RCurve(q) = ...
        Rhat+radiusBoundary*directionR*scaleR;

    logKvalue = ...
        logKhat+radiusBoundary*directionK*scaleK;

    kCurve(q) = exp(logKvalue);

    thresholdErrors(q) = bestThresholdError;
end

maxThresholdError = max(thresholdErrors);

end

%% =====================================================================
function delta = profiled_R_k_delta( ...
    radius,directionK,directionR, ...
    scaleK,scaleR,thetaHat,minNLL,data,bounds)

Rvalue = ...
    thetaHat(1)+radius*directionR*scaleR;

logKvalue = ...
    thetaHat(2)+radius*directionK*scaleK;

if Rvalue<bounds.lower(1) || ...
        Rvalue>bounds.upper(1) || ...
        logKvalue<bounds.lower(2) || ...
        logKvalue>bounds.upper(2)

    delta = 1e100;
    return;
end

delta = profile_R_k_at_values( ...
    Rvalue,logKvalue,thetaHat,data,bounds)-minNLL;

end

%% =====================================================================
function bestValue = profile_R_k_at_values( ...
    Rvalue,logKvalue,thetaHat,data,bounds)

template = thetaHat;
template(1) = Rvalue;
template(2) = logKvalue;

lower = bounds.lower(3:4);
upper = bounds.upper(3:4);

starts = [ ...
    min(max(thetaHat(3:4),lower),upper); ...
    0.50 0.50];

options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','off', ...
    'MaxFunctionEvaluations',600, ...
    'MaxIterations',160, ...
    'OptimalityTolerance',3e-7, ...
    'StepTolerance',3e-8);

bestValue = Inf;

for i = 1:size(starts,1)
    objective = @(detection) ...
        evaluate_fixed_R_k( ...
        detection,template,data,bounds);

    try
        [~,value] = fmincon( ...
            objective,starts(i,:),[],[],[],[], ...
            lower,upper,[],options);

        bestValue = min(bestValue,value);
    catch
    end
end

if ~isfinite(bestValue)
    bestValue = 1e100;
end

end

%% =====================================================================
function value = evaluate_fixed_R_k( ...
    detection,template,data,bounds)

theta = template;
theta(3:4) = detection;

value = JLE_nll(theta,data,bounds);

end

%% =====================================================================
function pS = true_cluster_pmf(R,k,Jmax)

j = 1:Jmax;

logP = ...
    -log(j) ...
    +gammaln((k+1).*j-1) ...
    -gammaln(k.*j) ...
    -gammaln(j) ...
    +(j-1).*log(R./k) ...
    -((k+1).*j-1).*log1p(R./k);

pS = exp(logP);

pS(~isfinite(pS)) = 0;

end

%% =====================================================================
function K = observation_kernel_vectorized( ...
    p_p,p_a,Jmax,observedSizes)

observedSizes = observedSizes(:)';

K = zeros(Jmax,numel(observedSizes));

s = (1:Jmax)';

z = p_p+(1-p_p)*p_a;

for q = 1:numel(observedSizes)
    o = observedSizes(q);

    if o==0
        K(:,q) = (1-p_p).^s;
        continue;
    end

    valid = s>=o;
    sv = s(valid);

    first = ...
        binomial_probability_fixed_o(sv,o,z);

    second = ...
        (1-p_p).^sv.* ...
        binomial_probability_fixed_o(sv,o,p_a);

    K(valid,q) = max(first-second,0);
end

end

%% =====================================================================
function values = ...
    binomial_probability_fixed_o(s,o,p)

values = zeros(size(s));

if p<=0
    if o==0
        values(:) = 1;
    end

elseif p>=1
    values(s==o) = 1;

else
    logChoose = ...
        gammaln(s+1) ...
        -gammaln(o+1) ...
        -gammaln(s-o+1);

    values = exp( ...
        logChoose ...
        +o.*log(p) ...
        +(s-o).*log1p(-p));
end

end

%% =====================================================================
function [kLimits,rLimits] = common_axis_limits( ...
    trueK,trueR,khat,klo,khi,Rhat,Rlo,Rhi, ...
    contourKCell,contourRCell,bounds,useLogK)

allContourK = cell2mat( ...
    cellfun(@(x) x(:),contourKCell, ...
    'UniformOutput',false));

allContourR = cell2mat( ...
    cellfun(@(x) x(:),contourRCell, ...
    'UniformOutput',false));

finiteK = [ ...
    trueK;khat(:);klo(:);khi(:);allContourK(:)];

finiteK = finiteK( ...
    isfinite(finiteK) & finiteK>0);

if useLogK
    logMin = log10(min(finiteK));
    logMax = log10(max(finiteK));

    padding = max( ...
        0.08,0.08*(logMax-logMin));

    kLimits = [ ...
        10^(logMin-padding), ...
        10^(logMax+padding)];
else
    kMin = min(finiteK);
    kMax = max(finiteK);

    span = max( ...
        kMax-kMin,0.18*kMax);

    kLimits = [ ...
        max(exp(bounds.lower(2)), ...
        kMin-0.35*span), ...
        min(exp(bounds.upper(2)), ...
        kMax+0.35*span)];
end

finiteR = [ ...
    trueR;Rhat(:);Rlo(:);Rhi(:);allContourR(:)];

finiteR = finiteR(isfinite(finiteR));

rMin = min(finiteR);
rMax = max(finiteR);

rSpan = max(rMax-rMin,0.10);

rLimits = [ ...
    max(bounds.lower(1), ...
    rMin-0.25*rSpan), ...
    min(bounds.upper(1), ...
    rMax+0.25*rSpan)];

end

%% =====================================================================
function draw_CI_cross( ...
    ax,khat,Rhat,klo,khi,Rlo,Rhi, ...
    kLimits,rLimits,color)

kloPlot = max(klo,kLimits(1));
khiPlot = min(khi,kLimits(2));

RloPlot = max(Rlo,rLimits(1));
RhiPlot = min(Rhi,rLimits(2));

plot(ax,[kloPlot khiPlot],[Rhat Rhat],'--', ...
    'Color',color,'LineWidth',1.7);

plot(ax,[khat khat],[RloPlot RhiPlot],'--', ...
    'Color',color,'LineWidth',1.7);

capR = 0.018*(rLimits(2)-rLimits(1));
capK = 0.012*(kLimits(2)-kLimits(1));

plot(ax,[kloPlot kloPlot], ...
    [Rhat-capR Rhat+capR],'-', ...
    'Color',color,'LineWidth',1.3);

plot(ax,[khiPlot khiPlot], ...
    [Rhat-capR Rhat+capR],'-', ...
    'Color',color,'LineWidth',1.3);

plot(ax,[khat-capK khat+capK], ...
    [RloPlot RloPlot],'-', ...
    'Color',color,'LineWidth',1.3);

plot(ax,[khat-capK khat+capK], ...
    [RhiPlot RhiPlot],'-', ...
    'Color',color,'LineWidth',1.3);

end

%% =====================================================================
function textValue = interval_text(value,atBound,side)

if atBound
    if strcmpi(side,'lower')
        textValue = sprintf('<=%.5g',value);
    else
        textValue = sprintf('>=%.5g',value);
    end
else
    textValue = sprintf('%.5g',value);
end

end

%% =====================================================================
function label = axis_label_k(useLogK)

if useLogK
    label = ...
        'Dispersion parameter, k (logarithmic scale)';
else
    label = ...
        'Dispersion parameter, k (linear scale)';
end

end

%% =====================================================================
function safe_diary_off

try
    diary off;
catch
end

end
function Figure4_JLE_Detection_Profiles_PublicationReady
% FIGURE4_JLE_DETECTION_PROFILES_PUBLICATIONREADY
%
% PURPOSE
% -------
% This program implements the requested revision of Figure 4. It displays
% the JLE surveillance estimates in the same 16-panel format used for the
% transmission results, but with passive detection p_p on the horizontal
% axis and active detection p_a on the vertical axis.
%
% The program uses the exact approved Figure 3 analysis cache. It does not
% simulate new clusters and it does not change any saved JLE point estimate.
% It performs the new calculations required for the revised figure:
%   1. separate 95% profile-likelihood intervals for p_p and p_a; and
%   2. a joint 95% profile-likelihood region for (p_p,p_a), profiling over
%      R_eff and k at every directly evaluated grid point.
%
% FIGURE SYMBOLS
% --------------
%   black boundary : joint 2-D 95% profile-likelihood region for (p_p,p_a)
%   blue circle    : JLE estimate of (p_p,p_a)
%   red diamond    : generating detection pair
%   blue dashed    : separate 1-D 95% profile-likelihood intervals
%
% When the generating p_p=1, active detection is structurally irrelevant.
% In those four panels, the code shows a blue vertical maximum-likelihood
% set at p_p=1 and labels p_a as N/I instead of displaying an arbitrary
% optimizer-dependent active-detection point estimate.

clc;
close all;

%% 1. Folders and file names

programFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(programFolder);
baseFolder = fullfile(projectRoot,'outputs');

if ~isfolder(baseFolder)
    mkdir(baseFolder);
end

if ~isfolder(baseFolder)
    error('The main manuscript folder was not found:\n%s',baseFolder);
end

sourceName = 'Figure3_JLE_16Panel_PublicationReady_v2';
sourceFolder = fullfile(baseFolder,sourceName);
sourceCacheFile = fullfile(sourceFolder,[sourceName '_Analysis_Cache.mat']);
sourceResultsCSV = fullfile(sourceFolder,[sourceName '_results.csv']);

if ~isfile(sourceCacheFile)
    error(['The approved Figure 3 analysis cache was not found.\n\n' ...
        'Expected file:\n%s\n\n' ...
        'Run Figure3_JLE_16Panel_Correct once before running this program.'], ...
        sourceCacheFile);
end

if ~isfile(sourceResultsCSV)
    error(['The approved Figure 3 results CSV was not found.\n\n' ...
        'Expected file:\n%s'],sourceResultsCSV);
end

outputName = 'Figure4_JLE_Detection_Profiles_PublicationReady';
outputFolder = fullfile(baseFolder,outputName);

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

figureBase = fullfile(outputFolder,'Figure 4');
pngFile = [figureBase '.png'];
pdfFile = [figureBase '.pdf'];
figFile = [figureBase '.fig'];
resultsCSV = fullfile(outputFolder,[outputName '_results.csv']);
surfaceCSV = fullfile(outputFolder,[outputName '_profile_surface.csv']);
contourCSV = fullfile(outputFolder,[outputName '_profile_contours.csv']);
profileCacheFile = fullfile(outputFolder,[outputName '_Profile_Cache.mat']);
readmeFile = fullfile(outputFolder,'README_Figure4_JLE_Detection_Profiles.txt');
logFile = fullfile(outputFolder,[outputName '_Command_Window_Log.txt']);

if isfile(logFile)
    delete(logFile);
end

diary(logFile);
diaryCleanup = onCleanup(@() safe_diary_off()); %#ok<NASGU>

startTimer = tic;

fprintf('============================================================\n');
fprintf('STARTING THE REVISED JLE FIGURE 4\n');
fprintf('Horizontal axis: passive detection probability, p_p.\n');
fprintf('Vertical axis: active detection probability, p_a.\n');
fprintf('No new cluster or auxiliary data will be simulated.\n');
fprintf('No approved JLE point estimate will be changed.\n');
fprintf('New detection profile-likelihood uncertainty will be calculated.\n');
fprintf('The first run is checkpointed panel by panel.\n');
fprintf('Later matching runs will load the saved profile cache.\n\n');
fprintf('APPROVED FIGURE 3 CACHE:\n%s\n\n',sourceCacheFile);
fprintf('ALL REVISED FIGURE 4 FILES WILL BE SAVED IN:\n%s\n',outputFolder);
fprintf('============================================================\n\n');

%% 2. Load and verify the exact approved Figure 3 analysis

requiredSourceVariables = { ...
    'cacheVersion','R_true','k_true','p_p_grid','p_a_grid', ...
    'N_observed','N_auxiliary','bounds','threshold1D','threshold2D', ...
    'panelLabel','truePP','truePA','truePobs','auxiliarySuccesses', ...
    'R_hat','k_hat','pp_hat','pa_hat','pobs_hat','minimumNLL', ...
    'dataCell','thetaHatCell'};

sourceInfo = whos('-file',sourceCacheFile);
sourceNames = {sourceInfo.name};

missingSourceVariables = requiredSourceVariables( ...
    ~ismember(requiredSourceVariables,sourceNames));

if ~isempty(missingSourceVariables)
    error('The Figure 3 cache is missing the variable: %s', ...
        missingSourceVariables{1});
end

S = load(sourceCacheFile,requiredSourceVariables{:});
Tsource = readtable(sourceResultsCSV);

numberPanels = 16;

if numel(S.panelLabel)~=numberPanels || ...
        numel(S.dataCell)~=numberPanels || ...
        numel(S.thetaHatCell)~=numberPanels
    error('The Figure 3 cache must contain exactly 16 displayed panels.');
end

if S.N_observed~=5000 || S.N_auxiliary~=5000
    error('The approved Figure 3 sample sizes are not 5,000 and 5,000.');
end

if height(Tsource)~=numberPanels
    error('The approved Figure 3 results CSV must contain 16 rows.');
end

for panel = 1:numberPanels
    if sum(S.dataCell{panel}.counts)~=S.N_observed
        error('Panel %d does not contain 5,000 observed clusters.',panel);
    end

    if S.dataCell{panel}.Naux~=S.N_auxiliary
        error('Panel %d does not contain 5,000 auxiliary trials.',panel);
    end
end

if max(abs(Tsource.p_p_hat-S.pp_hat))>1e-10 || ...
        max(abs(Tsource.p_a_hat-S.pa_hat))>1e-10 || ...
        max(abs(Tsource.MinimumJointNLL-S.minimumNLL))>1e-7
    error('The Figure 3 CSV and analysis cache do not agree.');
end

if max(abs(S.truePobs-(S.truePP+(1-S.truePP).*S.truePA)))>1e-12
    error('The generating combined probabilities failed their identity check.');
end

fprintf('Approved Figure 3 cache checks passed.\n');
fprintf('Displayed panels = 16.\n');
fprintf('Unique independently simulated analyses = 13.\n');
fprintf('Observed clusters per unique analysis = 5000.\n');
fprintf('Independent auxiliary trials per unique analysis = 5000.\n\n');

%% 3. Detection-profile settings and optional parallel pool

% The directly evaluated surface includes probability values from 0 to 1.
% The p_p=0 plotting point is represented numerically by the exact positive
% lower bound used in Figure 3 because the conditional cluster likelihood
% is undefined at p_p=0.

gridStep = 0.025;
ppGrid = unique([S.bounds.lower(3),gridStep:gridStep:1.0]);
paGrid = 0:gridStep:1.0;

profileCacheVersion = ...
    'Figure4-JLE-detection-profiles-v1-2026-08-25';
sourceCacheVersion = S.cacheVersion;

useParallel = false;

try
    if license('test','Distrib_Computing_Toolbox')
        pool = gcp('nocreate');

        if isempty(pool)
            pool = parpool('local');
        end

        useParallel = ~isempty(pool);
    end
catch poolError
    fprintf('Parallel pool was not started: %s\n',poolError.message);
    useParallel = false;
end

if useParallel
    fprintf('Parallel processing is available for the profile surface.\n\n');
else
    fprintf(['Parallel processing is unavailable; ' ...
        'a serial calculation will be used.\n\n']);
end

%% 4. Initialize or resume the checkpointed detection profiles

ppLower = NaN(numberPanels,1);
ppUpper = NaN(numberPanels,1);
paLower = NaN(numberPanels,1);
paUpper = NaN(numberPanels,1);

ppLowerAtBound = false(numberPanels,1);
ppUpperAtBound = false(numberPanels,1);
paLowerAtBound = false(numberPanels,1);
paUpperAtBound = false(numberPanels,1);

truePairDeltaNLL = NaN(numberPanels,1);
truthCoveredJoint = false(numberPanels,1);
jointRegionTouchesBoundary = false(numberPanels,1);
profileDeltaCell = cell(numberPanels,1);
panelProfileRuntimeSeconds = NaN(numberPanels,1);
completedProfiles = false(numberPanels,1);

loadedProfileCache = false;

if isfile(profileCacheFile)
    P = load(profileCacheFile);

    cacheMatches = isfield(P,'profileCacheVersion') && ...
        isfield(P,'sourceCacheVersion') && ...
        isfield(P,'gridStep') && ...
        strcmp(P.profileCacheVersion,profileCacheVersion) && ...
        strcmp(P.sourceCacheVersion,sourceCacheVersion) && ...
        abs(P.gridStep-gridStep)<1e-14;

    if cacheMatches
        ppLower = P.ppLower;
        ppUpper = P.ppUpper;
        paLower = P.paLower;
        paUpper = P.paUpper;

        ppLowerAtBound = P.ppLowerAtBound;
        ppUpperAtBound = P.ppUpperAtBound;
        paLowerAtBound = P.paLowerAtBound;
        paUpperAtBound = P.paUpperAtBound;

        truePairDeltaNLL = P.truePairDeltaNLL;
        truthCoveredJoint = P.truthCoveredJoint;
        jointRegionTouchesBoundary = P.jointRegionTouchesBoundary;
        profileDeltaCell = P.profileDeltaCell;
        panelProfileRuntimeSeconds = P.panelProfileRuntimeSeconds;
        completedProfiles = P.completedProfiles;

        loadedProfileCache = true;

        fprintf('A matching Figure 4 detection-profile cache was found.\n');
        fprintf('Completed panels in the cache = %d of 16.\n\n', ...
            sum(completedProfiles));
    end
end

if ~loadedProfileCache
    fprintf('No matching detection-profile cache was found.\n');
    fprintf('The required profiles will be calculated and checkpointed.\n\n');
end

% Panels 4, 8, 12, and 16 use the same p_p=1 data and likelihood.

canonicalPanel = (1:numberPanels)';
canonicalPanel([8 12 16]) = 4;

uniquePanels = [1 2 3 4 5 6 7 9 10 11 13 14 15];

for index = 1:numel(uniquePanels)

    panel = uniquePanels(index);

    if completedProfiles(panel)
        fprintf('%s detection profile already complete; loaded from checkpoint.\n', ...
            S.panelLabel{panel});
        continue;
    end

    panelTimer = tic;

    fprintf(['%s profiling detection for true p_p=%.2f, p_a=%.2f ' ...
        '(%d of %d unique analyses)...\n'], ...
        S.panelLabel{panel},S.truePP(panel),S.truePA(panel), ...
        index,numel(uniquePanels));

    thetaHat = S.thetaHatCell{panel};
    minNLL = S.minimumNLL(panel);
    data = S.dataCell{panel};

    [ppCI,ppFlags] = profile_CI_one_parameter( ...
        thetaHat,minNLL,3,data,S.bounds,S.threshold1D);

    [paCI,paFlags] = profile_CI_one_parameter( ...
        thetaHat,minNLL,4,data,S.bounds,S.threshold1D);

    ppLower(panel) = ppCI(1);
    ppUpper(panel) = ppCI(2);
    paLower(panel) = paCI(1);
    paUpper(panel) = paCI(2);

    ppLowerAtBound(panel) = ppFlags.lower;
    ppUpperAtBound(panel) = ppFlags.upper;
    paLowerAtBound(panel) = paFlags.lower;
    paUpperAtBound(panel) = paFlags.upper;

    profileNLL = profile_detection_surface( ...
        thetaHat,data,S.bounds,ppGrid,paGrid,useParallel);

    rawDelta = profileNLL-minNLL;
    minimumRawDelta = min(rawDelta(:));

    if minimumRawDelta < -0.03
        error(['The fixed-detection profile found an NLL %.5g below the ' ...
            'approved global minimum in panel %s. ' ...
            'Recheck the Figure 3 fit.'], ...
            minimumRawDelta,S.panelLabel{panel});
    end

    profileDeltaCell{panel} = max(rawDelta,0);

    truePairNLL = profile_nll_fixed_detection( ...
        S.truePP(panel),S.truePA(panel),thetaHat,data,S.bounds);

    truePairDeltaNLL(panel) = max(0,truePairNLL-minNLL);

    truthCoveredJoint(panel) = ...
        truePairDeltaNLL(panel)<=S.threshold2D+0.02;

    accepted = profileDeltaCell{panel}<=S.threshold2D;

    jointRegionTouchesBoundary(panel) = ...
        any(accepted(1,:)) || any(accepted(end,:)) || ...
        any(accepted(:,1)) || any(accepted(:,end));

    panelProfileRuntimeSeconds(panel) = toc(panelTimer);
    completedProfiles(panel) = true;

    fprintf('    p_p 95%% profile CI = [%s, %s].\n', ...
        interval_text(ppLower(panel),ppLowerAtBound(panel),'lower'), ...
        interval_text(ppUpper(panel),ppUpperAtBound(panel),'upper'));

    fprintf('    p_a 95%% profile CI = [%s, %s].\n', ...
        interval_text(paLower(panel),paLowerAtBound(panel),'lower'), ...
        interval_text(paUpper(panel),paUpperAtBound(panel),'upper'));

    fprintf('    True detection-pair profile delta NLL = %.5f.\n', ...
        truePairDeltaNLL(panel));

    fprintf('    Panel profile runtime = %.1f seconds.\n', ...
        panelProfileRuntimeSeconds(panel));

    save(profileCacheFile, ...
        'profileCacheVersion','sourceCacheVersion','gridStep', ...
        'ppGrid','paGrid','ppLower','ppUpper','paLower','paUpper', ...
        'ppLowerAtBound','ppUpperAtBound', ...
        'paLowerAtBound','paUpperAtBound', ...
        'truePairDeltaNLL','truthCoveredJoint', ...
        'jointRegionTouchesBoundary','profileDeltaCell', ...
        'panelProfileRuntimeSeconds','completedProfiles','-v7.3');

    fprintf('    Checkpoint saved.\n\n');
end

% Copy the single p_p=1 likelihood analysis to its equivalent panels.

for panel = [8 12 16]

    sourcePanel = canonicalPanel(panel);

    ppLower(panel) = ppLower(sourcePanel);
    ppUpper(panel) = ppUpper(sourcePanel);
    paLower(panel) = paLower(sourcePanel);
    paUpper(panel) = paUpper(sourcePanel);

    ppLowerAtBound(panel) = ppLowerAtBound(sourcePanel);
    ppUpperAtBound(panel) = ppUpperAtBound(sourcePanel);
    paLowerAtBound(panel) = paLowerAtBound(sourcePanel);
    paUpperAtBound(panel) = paUpperAtBound(sourcePanel);

    profileDeltaCell{panel} = profileDeltaCell{sourcePanel};
    panelProfileRuntimeSeconds(panel) = 0;

    truePairNLL = profile_nll_fixed_detection( ...
        S.truePP(panel),S.truePA(panel), ...
        S.thetaHatCell{sourcePanel}, ...
        S.dataCell{sourcePanel},S.bounds);

    truePairDeltaNLL(panel) = ...
        max(0,truePairNLL-S.minimumNLL(sourcePanel));

    truthCoveredJoint(panel) = ...
        truePairDeltaNLL(panel)<=S.threshold2D+0.02;

    jointRegionTouchesBoundary(panel) = ...
        jointRegionTouchesBoundary(sourcePanel);

    completedProfiles(panel) = true;
end

if ~all(completedProfiles)
    error('Not all 16 displayed detection profiles were completed.');
end

save(profileCacheFile, ...
    'profileCacheVersion','sourceCacheVersion','gridStep', ...
    'ppGrid','paGrid','ppLower','ppUpper','paLower','paUpper', ...
    'ppLowerAtBound','ppUpperAtBound', ...
    'paLowerAtBound','paUpperAtBound', ...
    'truePairDeltaNLL','truthCoveredJoint', ...
    'jointRegionTouchesBoundary','profileDeltaCell', ...
    'panelProfileRuntimeSeconds','completedProfiles','-v7.3');

fprintf('All 16 displayed detection profiles are complete.\n\n');

%% 5. Final scientific checks and summary measures

if any(S.pp_hat<ppLower-1e-6 | S.pp_hat>ppUpper+1e-6)
    error('At least one p_p estimate lies outside its own 95%% profile interval.');
end

identifiablePA = S.truePP<1-1e-12;

if any(S.pa_hat(identifiablePA)<paLower(identifiablePA)-1e-6 | ...
        S.pa_hat(identifiablePA)>paUpper(identifiablePA)+1e-6)
    error('At least one identifiable p_a estimate lies outside its own interval.');
end

if any(~truthCoveredJoint)
    fprintf(['CAUTION: the generating detection pair is outside the ' ...
        'joint 95%% profile region in %d displayed panel(s).\n'], ...
        sum(~truthCoveredJoint));
end

ppMAE = mean(abs(S.pp_hat-S.truePP));

paMAE = mean(abs( ...
    S.pa_hat(identifiablePA)-S.truePA(identifiablePA)));

pobsMAE = mean(abs(S.pobs_hat-S.truePobs));

fprintf('Mean absolute error for p_p = %.6f.\n',ppMAE);
fprintf('Mean absolute error for identifiable p_a = %.6f.\n',paMAE);
fprintf('Mean absolute error for combined detection = %.6f.\n',pobsMAE);

fprintf('Joint detection regions containing truth = %d of 16.\n', ...
    sum(truthCoveredJoint));

fprintf('Joint detection regions touching a probability boundary = %d.\n\n', ...
    sum(jointRegionTouchesBoundary));

%% 6. Save the numerical result table

displayedPAHat = S.pa_hat;
displayedPAHat(~identifiablePA) = NaN;

resultTable = table( ...
    S.panelLabel(:),S.truePP(:),S.truePA(:),S.truePobs(:), ...
    S.pp_hat(:),ppLower,ppUpper,ppLowerAtBound,ppUpperAtBound, ...
    displayedPAHat(:),S.pa_hat(:),paLower,paUpper, ...
    paLowerAtBound,paUpperAtBound,identifiablePA(:), ...
    S.pobs_hat(:),truePairDeltaNLL,truthCoveredJoint, ...
    jointRegionTouchesBoundary,panelProfileRuntimeSeconds, ...
    'VariableNames',{ ...
    'Panel','True_p_p','True_p_a','True_combined_detection', ...
    'JLE_p_p_hat','p_p_lower_1D95','p_p_upper_1D95', ...
    'p_p_lower_reached_bound','p_p_upper_reached_bound', ...
    'Displayed_JLE_p_a_hat','Raw_optimizer_p_a_hat', ...
    'p_a_lower_1D95','p_a_upper_1D95', ...
    'p_a_lower_reached_bound','p_a_upper_reached_bound', ...
    'p_a_structurally_identifiable','JLE_combined_detection_hat', ...
    'True_pair_profile_delta_NLL', ...
    'True_pair_covered_by_joint_2D_region', ...
    'Joint_region_touches_probability_boundary', ...
    'Profile_runtime_seconds'});

writetable(resultTable,resultsCSV);

%% 7. Save the profile surface and contour coordinates

numberGridRows = numberPanels*numel(paGrid)*numel(ppGrid);

surfacePanel = cell(numberGridRows,1);
surfacePP = zeros(numberGridRows,1);
surfacePA = zeros(numberGridRows,1);
surfaceDelta = zeros(numberGridRows,1);
surfaceInside = false(numberGridRows,1);

rowCounter = 0;

for panel = 1:numberPanels

    [PPmesh,PAmesh] = meshgrid(ppGrid,paGrid);
    numberPoints = numel(PPmesh);
    rows = rowCounter+(1:numberPoints);

    surfacePanel(rows) = ...
        repmat(S.panelLabel(panel),numberPoints,1);

    surfacePP(rows) = PPmesh(:);
    surfacePA(rows) = PAmesh(:);
    surfaceDelta(rows) = profileDeltaCell{panel}(:);

    surfaceInside(rows) = ...
        profileDeltaCell{panel}(:)<=S.threshold2D;

    rowCounter = rowCounter+numberPoints;
end

surfaceTable = table( ...
    surfacePanel,surfacePP,surfacePA,surfaceDelta,surfaceInside, ...
    'VariableNames',{'Panel','p_p','p_a','ProfileDeltaNLL', ...
    'InsideJoint95Region'});

writetable(surfaceTable,surfaceCSV);

contourPanel = cell(0,1);
contourSegment = zeros(0,1);
contourPoint = zeros(0,1);
contourPP = zeros(0,1);
contourPA = zeros(0,1);

for panel = 1:numberPanels

    segments = extract_contour_segments( ...
        ppGrid,paGrid,profileDeltaCell{panel},S.threshold2D);

    for segment = 1:numel(segments)

        numberPoints = numel(segments{segment}.x);

        contourPanel = [contourPanel; ...
            repmat(S.panelLabel(panel),numberPoints,1)]; %#ok<AGROW>

        contourSegment = [contourSegment; ...
            repmat(segment,numberPoints,1)]; %#ok<AGROW>

        contourPoint = [contourPoint;(1:numberPoints)']; %#ok<AGROW>
        contourPP = [contourPP;segments{segment}.x(:)]; %#ok<AGROW>
        contourPA = [contourPA;segments{segment}.y(:)]; %#ok<AGROW>
    end
end

contourTable = table( ...
    contourPanel,contourSegment,contourPoint,contourPP,contourPA, ...
    'VariableNames',{'Panel','ContourSegment','ContourPoint','p_p','p_a'});

writetable(contourTable,contourCSV);

%% 8. Draw the revised clean 16-panel Figure 4

fig = figure( ...
    'Color','w', ...
    'Name',outputName, ...
    'Units','inches', ...
    'Position',[0.25 0.25 16 10], ...
    'Renderer','painters', ...
    'InvertHardcopy','on');

leftMargin = 0.070;
rightMargin = 0.025;
bottomMargin = 0.085;
topOfPanels = 0.850;
horizontalGap = 0.028;
verticalGap = 0.055;

panelWidth = ...
    (1-leftMargin-rightMargin-3*horizontalGap)/4;

panelHeight = ...
    (topOfPanels-bottomMargin-3*verticalGap)/4;

legendHandles = gobjects(1,4);

for panel = 1:numberPanels

    row = ceil(panel/4);
    column = mod(panel-1,4)+1;

    xPosition = ...
        leftMargin+(column-1)*(panelWidth+horizontalGap);

    yPosition = ...
        topOfPanels-row*panelHeight-(row-1)*verticalGap;

    ax = axes( ...
        'Parent',fig, ...
        'Position',[xPosition yPosition panelWidth panelHeight]);

    hold(ax,'on');

    [PPmesh,PAmesh] = meshgrid(ppGrid,paGrid);

    [~,hContour] = contour( ...
        ax,PPmesh,PAmesh,profileDeltaCell{panel}, ...
        [S.threshold2D S.threshold2D], ...
        'k-','LineWidth',1.4);

    hTruth = plot( ...
        ax,S.truePP(panel),S.truePA(panel),'d', ...
        'MarkerSize',6.5, ...
        'MarkerFaceColor',[1 0.10 0.08], ...
        'MarkerEdgeColor',[1 0.10 0.08], ...
        'LineWidth',1.0);

    if identifiablePA(panel)

        hHorizontal = plot( ...
            ax,[ppLower(panel) ppUpper(panel)], ...
            [S.pa_hat(panel) S.pa_hat(panel)],'--', ...
            'Color',[0.05 0.35 0.95], ...
            'LineWidth',1.15);

        plot( ...
            ax,[S.pp_hat(panel) S.pp_hat(panel)], ...
            [paLower(panel) paUpper(panel)],'--', ...
            'Color',[0.05 0.35 0.95], ...
            'LineWidth',1.15);

        hEstimate = plot( ...
            ax,S.pp_hat(panel),S.pa_hat(panel),'o', ...
            'MarkerSize',5.8, ...
            'MarkerFaceColor',[0.02 0.36 0.95], ...
            'MarkerEdgeColor',[0.02 0.30 0.80], ...
            'LineWidth',0.9);

    else

        % At p_p=1 the complete vertical segment is a likelihood-equivalent
        % maximum set for p_a. No arbitrary p_a estimate is displayed.

        hEstimate = plot( ...
            ax,[1 1],[0 1],'-', ...
            'Color',[0.02 0.36 0.95], ...
            'LineWidth',2.0);

        hHorizontal = plot( ...
            ax,[ppLower(panel) ppUpper(panel)], ...
            [0.06 0.06],'--', ...
            'Color',[0.05 0.35 0.95], ...
            'LineWidth',1.15);

        text( ...
            ax,0.96,0.13,'p_a N/I', ...
            'HorizontalAlignment','right', ...
            'VerticalAlignment','middle', ...
            'FontSize',8.3, ...
            'FontWeight','bold', ...
            'Interpreter','tex');
    end

    if panel==1
        legendHandles = ...
            [hContour hEstimate hTruth hHorizontal];
    end

    % Small margins prevent boundary markers from being clipped.
    xlim(ax,[-0.025 1.025]);
    ylim(ax,[-0.025 1.025]);

    xticks(ax,[0 0.25 0.50 0.75 1.00]);
    yticks(ax,[0 0.25 0.50 0.75 1.00]);

    xticklabels(ax,{'0','0.25','0.50','0.75','1'});
    yticklabels(ax,{'0','0.25','0.50','0.75','1'});

    grid(ax,'on');
    ax.GridAlpha = 0.16;
    ax.LineWidth = 0.75;
    ax.FontSize = 8.2;
    ax.Layer = 'top';

    box(ax,'on');

    if row<4
        ax.XTickLabel = [];
    end

    if column>1
        ax.YTickLabel = [];
    end

    title( ...
        ax, ...
        sprintf('(%c) p_p=%.2f, p_a=%.2f', ...
        char(96+panel),S.truePP(panel),S.truePA(panel)), ...
        'FontSize',9.3, ...
        'FontWeight','bold', ...
        'Interpreter','tex');
end

annotation( ...
    fig,'textbox',[0.04 0.943 0.92 0.045], ...
    'String','JLE detection estimates', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',18, ...
    'FontWeight','bold', ...
    'EdgeColor','none', ...
    'Interpreter','tex');

legendObject = legend( ...
    legendHandles, ...
    {'95% joint region','JLE estimate','Truth','95% 1-D intervals'}, ...
    'Orientation','horizontal', ...
    'Box','off', ...
    'FontSize',10.2);

legendObject.Position = [0.245 0.892 0.51 0.035];

annotation( ...
    fig,'textbox',[0.29 0.010 0.42 0.045], ...
    'String','Passive detection probability, p_p', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',12.5, ...
    'FontWeight','bold', ...
    'EdgeColor','none', ...
    'Interpreter','tex');

labelAxes = axes( ...
    'Parent',fig, ...
    'Position',[0 0 1 1], ...
    'Visible','off', ...
    'HitTest','off');

text( ...
    labelAxes,0.016,0.47, ...
    'Active detection probability, p_a', ...
    'Rotation',90, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',12.5, ...
    'FontWeight','bold', ...
    'Interpreter','tex', ...
    'Visible','on');

drawnow;

%% 9. Save PDF, PNG, FIG, code, README, and log

fprintf('Saving the revised Figure 4...\n');

savefig(fig,figFile);

set(fig,'PaperUnits','inches');
set(fig,'PaperPosition',[0 0 16 10]);
set(fig,'PaperSize',[16 10]);
set(fig,'Color','w','InvertHardcopy','on');

print(fig,pngFile,'-dpng','-r300');
print(fig,pdfFile,'-dpdf','-painters');

close(fig);

sourceCodeFile = [mfilename('fullpath') '.m'];
savedCodeFile = fullfile(outputFolder,[mfilename '.m']);

if isfile(sourceCodeFile) && ~strcmpi(sourceCodeFile,savedCodeFile)
    copyfile(sourceCodeFile,savedCodeFile,'f');
end

fid = fopen(readmeFile,'w');

if fid~=-1

    fprintf(fid,'REVISED FIGURE 4: JLE DETECTION PROFILES\r\n\r\n');

    fprintf(fid,'Scientific purpose\r\n');

    fprintf(fid,['This figure displays passive detection p_p on the ' ...
        'horizontal axis and active detection p_a on the vertical axis.\r\n']);

    fprintf(fid,['It uses the exact approved Figure 3 data and JLE point ' ...
        'estimates. No new data were simulated.\r\n']);

    fprintf(fid,['The separate 95%% intervals and joint 95%% region are new ' ...
        'profile-likelihood calculations requested for this figure.\r\n']);

    fprintf(fid,['At every fixed p_p-p_a grid point, R_eff and k were ' ...
        're-optimized.\r\n']);

    fprintf(fid,['The 1-D and 2-D regions use the chi-square ' ...
        'likelihood-ratio thresholds for one and two parameters, ' ...
        'respectively.\r\n']);

    fprintf(fid,['When generating p_p=1, p_a is structurally not ' ...
        'identifiable; the blue vertical line represents the ' ...
        'likelihood-equivalent maximum set.\r\n\r\n']);

    fprintf(fid,'Source Figure 3 cache: %s\r\n',sourceCacheFile);
    fprintf(fid,'Observed clusters per unique analysis: %d\r\n',S.N_observed);
    fprintf(fid,'Auxiliary trials per unique analysis: %d\r\n',S.N_auxiliary);
    fprintf(fid,'Direct profile grid step: %.3f\r\n',gridStep);
    fprintf(fid,'Mean absolute error for p_p: %.6f\r\n',ppMAE);

    fprintf(fid,'Mean absolute error for identifiable p_a: %.6f\r\n', ...
        paMAE);

    fprintf(fid,'Mean absolute error for combined detection: %.6f\r\n', ...
        pobsMAE);

    fprintf(fid,'Joint detection regions containing truth: %d of 16\r\n', ...
        sum(truthCoveredJoint));

    fclose(fid);
end

fprintf('\n============================================================\n');
fprintf('DONE: THE REVISED JLE FIGURE 4 WAS CREATED.\n');
fprintf('Elapsed time for this run: %.1f seconds.\n',toc(startTimer));
fprintf('No new cluster or auxiliary data were simulated.\n');
fprintf('No approved JLE point estimate was changed.\n');
fprintf('New p_p and p_a profile-likelihood uncertainty was calculated.\n');
fprintf('Displayed surveillance settings = 16.\n');
fprintf('Unique independently simulated analyses = 13.\n');

fprintf('Joint detection regions containing truth = %d of 16.\n', ...
    sum(truthCoveredJoint));

fprintf('Structurally non-identifiable p_a panels = 4.\n');
fprintf('Mean absolute error for p_p = %.6f.\n',ppMAE);

fprintf('Mean absolute error for identifiable p_a = %.6f.\n', ...
    paMAE);

fprintf('Mean absolute error for combined detection = %.6f.\n', ...
    pobsMAE);

fprintf('\nALL REVISED FIGURE 4 FILES WERE SAVED IN:\n%s\n', ...
    outputFolder);

fprintf('\nPDF FIGURE:\n%s\n',pdfFile);
fprintf('\nPNG FIGURE:\n%s\n',pngFile);
fprintf('\nMATLAB FIG FILE:\n%s\n',figFile);
fprintf('\nRESULTS CSV:\n%s\n',resultsCSV);
fprintf('\nPROFILE-SURFACE CSV:\n%s\n',surfaceCSV);
fprintf('\nPROFILE-CONTOUR CSV:\n%s\n',contourCSV);
fprintf('\nREUSABLE PROFILE CACHE:\n%s\n',profileCacheFile);
fprintf('\nSAVED CODE COPY:\n%s\n',savedCodeFile);
fprintf('\nREADME:\n%s\n',readmeFile);
fprintf('\nCOMMAND-WINDOW LOG:\n%s\n',logFile);
fprintf('============================================================\n');

end

%% =====================================================================
function value = JLE_nll(theta,data,bounds)

theta = theta(:)';

if numel(theta)~=4 || any(~isfinite(theta)) || ...
        any(theta<bounds.lower) || any(theta>bounds.upper)

    value = 1e100;
    return;
end

R = theta(1);
k = exp(theta(2));
p_p = theta(3);
p_a = theta(4);

pS = true_cluster_pmf(R,k,data.Jmax);

[positiveSelected,pZero] = ...
    selected_observation_probabilities(pS,p_p,p_a,data);

denominator = sum(pS)-pZero;

if denominator<=0 || any(positiveSelected<=0) || ...
        any(~isfinite(positiveSelected))

    value = 1e100;
    return;
end

pConditional = positiveSelected/denominator;
nllClusters = -sum(data.counts.*log(pConditional));

pObserved = p_p+(1-p_p)*p_a;

logAuxiliary = exact_binomial_log_likelihood( ...
    data.auxSuccesses,data.Naux,pObserved);

if ~isfinite(logAuxiliary)
    value = 1e100;
    return;
end

value = nllClusters-logAuxiliary;

if ~isfinite(value) || ~isreal(value)
    value = 1e100;
end

end

%% =====================================================================
function logLikelihood = exact_binomial_log_likelihood(Y,N,p)

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

logLikelihood = ...
    gammaln(N+1)-gammaln(Y+1)-gammaln(N-Y+1);

if Y>0
    logLikelihood = logLikelihood+Y*log(p);
end

if N-Y>0
    logLikelihood = logLikelihood+(N-Y)*log1p(-p);
end

end

%% =====================================================================
function [positiveSelected,pZero] = ...
    selected_observation_probabilities(pS,p_p,p_a,data)

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
    logTerm1 = ...
        logChoose+o*log(z)+remaining*log1p(-z);

    term1 = exp(logTerm1);
end

if p_p>=1 || p_a<=0
    term2 = zeros(size(logChoose));
elseif p_a>=1
    term2 = (1-p_p).^s.*double(s==o);
else
    logTerm2 = ...
        s*log1p(-p_p)+logChoose+ ...
        o*log(p_a)+remaining*log1p(-p_a);

    term2 = exp(logTerm2);
end

kernelPositive = max(term1-term2,0);
positiveSelected = pS*kernelPositive;
pZero = sum(pS(:).*((1-p_p).^s));

end

%% =====================================================================
function pS = true_cluster_pmf(R,k,Jmax)

j = 1:Jmax;

logProbability = ...
    -log(j) ...
    +gammaln(j*k+j-1) ...
    -gammaln(j*k) ...
    -gammaln(j) ...
    +j*k*log(k/(k+R)) ...
    +(j-1)*log(R/(k+R));

pS = exp(logProbability);
pS(~isfinite(pS)) = 0;
pS = pS(:)';

end

%% =====================================================================
function [CI,flags] = profile_CI_one_parameter( ...
    thetaHat,minNLL,targetIndex,data,bounds,threshold)

target = minNLL+threshold;
estimate = thetaHat(targetIndex);

[lower,flags.lower] = find_profile_crossing( ...
    thetaHat,targetIndex,estimate,bounds.lower(targetIndex), ...
    target,data,bounds);

[upper,flags.upper] = find_profile_crossing( ...
    thetaHat,targetIndex,estimate,bounds.upper(targetIndex), ...
    target,data,bounds);

CI = [lower upper];

end

%% =====================================================================
function [crossing,reachedBoundary] = find_profile_crossing( ...
    thetaHat,targetIndex,estimate,boundary,target,data,bounds)

if abs(estimate-boundary)<1e-12

    boundaryDifference = profile_nll_fixed_parameter( ...
        boundary,targetIndex,thetaHat,data,bounds)-target;

    if boundaryDifference<=0
        crossing = boundary;
        reachedBoundary = true;
        return;
    end
end

screen = linspace(estimate,boundary,16);

previousX = screen(1);

previousDifference = profile_nll_fixed_parameter( ...
    previousX,targetIndex,thetaHat,data,bounds)-target;

bracketFound = false;

for i = 2:numel(screen)

    currentX = screen(i);

    currentDifference = profile_nll_fixed_parameter( ...
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

options = optimset('Display','off','TolX',2e-5);

crossing = fzero( ...
    @(x) profile_nll_fixed_parameter( ...
    x,targetIndex,thetaHat,data,bounds)-target, ...
    [lowerBracket upperBracket],options);

reachedBoundary = false;

end

%% =====================================================================
function bestValue = profile_nll_fixed_parameter( ...
    fixedValue,targetIndex,thetaHat,data,bounds)

remaining = setdiff(1:4,targetIndex);

template = thetaHat;
template(targetIndex) = fixedValue;

lower = bounds.lower(remaining);
upper = bounds.upper(remaining);

start1 = min(max(thetaHat(remaining),lower),upper);
start2 = 0.5*(lower+upper);

logKposition = find(remaining==2,1);

if ~isempty(logKposition)
    start2(logKposition) = log(0.15);
end

starts = [start1;start2];

options = optimoptions( ...
    'fmincon', ...
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
function profileNLL = profile_detection_surface( ...
    thetaHat,data,bounds,ppGrid,paGrid,useParallel)

[PPmesh,PAmesh] = meshgrid(ppGrid,paGrid);

ppValues = PPmesh(:);
paValues = PAmesh(:);

numberPoints = numel(ppValues);
values = NaN(numberPoints,1);

if useParallel

    parfor q = 1:numberPoints
        values(q) = profile_nll_fixed_detection( ...
            ppValues(q),paValues(q),thetaHat,data,bounds);
    end

else

    for q = 1:numberPoints
        values(q) = profile_nll_fixed_detection( ...
            ppValues(q),paValues(q),thetaHat,data,bounds);
    end

end

profileNLL = reshape(values,size(PPmesh));

end

%% =====================================================================
function bestValue = profile_nll_fixed_detection( ...
    p_p,p_a,thetaHat,data,bounds)

persistent fixedDetectionOptions

if isempty(fixedDetectionOptions)

    fixedDetectionOptions = optimoptions( ...
        'fmincon', ...
        'Algorithm','sqp', ...
        'Display','off', ...
        'MaxFunctionEvaluations',650, ...
        'MaxIterations',180, ...
        'OptimalityTolerance',3e-7, ...
        'StepTolerance',3e-8);
end

if p_p<bounds.lower(3) || p_p>bounds.upper(3) || ...
        p_a<bounds.lower(4) || p_a>bounds.upper(4)

    bestValue = 1e100;
    return;
end

template = thetaHat;
template(3:4) = [p_p p_a];

lower = bounds.lower(1:2);
upper = bounds.upper(1:2);

starts = [ ...
    min(max(thetaHat(1:2),lower),upper); ...
    0.50 log(0.15)];

bestValue = Inf;

for i = 1:size(starts,1)

    objective = @(transmission) evaluate_fixed_detection( ...
        transmission,template,data,bounds);

    try
        [~,value] = fmincon( ...
            objective,starts(i,:),[],[],[],[], ...
            lower,upper,[],fixedDetectionOptions);

        bestValue = min(bestValue,value);
    catch
    end
end

if ~isfinite(bestValue)
    bestValue = 1e100;
end

end

%% =====================================================================
function value = evaluate_fixed_detection( ...
    transmission,template,data,bounds)

theta = template;
theta(1:2) = transmission;
value = JLE_nll(theta,data,bounds);

end

%% =====================================================================
function segments = extract_contour_segments(x,y,z,level)

C = contourc(x,y,z,[level level]);

segments = cell(0,1);
column = 1;

while column<size(C,2)

    numberPoints = round(C(2,column));

    if numberPoints<=0 || column+numberPoints>size(C,2)
        break;
    end

    segment.x = C(1,column+(1:numberPoints));
    segment.y = C(2,column+(1:numberPoints));

    segments{end+1,1} = segment; %#ok<AGROW>

    column = column+numberPoints+1;
end

end

%% =====================================================================
function textValue = interval_text(value,atBound,side)

if atBound

    if strcmpi(side,'lower')
        textValue = sprintf('<=%.4f',value);
    else
        textValue = sprintf('>=%.4f',value);
    end

else
    textValue = sprintf('%.4f',value);
end

end

%% =====================================================================
function safe_diary_off

try
    diary off;
catch
end

end
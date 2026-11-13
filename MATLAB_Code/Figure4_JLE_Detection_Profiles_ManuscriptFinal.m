function Figure4_JLE_Detection_Profiles_ManuscriptFinal
% FIGURE4_JLE_DETECTION_PROFILES_MANUSCRIPTFINAL
%
% PURPOSE
% -------
% This program creates the manuscript-sized detection-parameter version of
% Figure 4.
% It uses passive detection p_p on the horizontal axis and active detection
% p_a on the vertical axis, matching the inferential layout of Figure 3.
%
% IMPORTANT SCIENTIFIC RULES
% --------------------------
% 1. No cluster or auxiliary observation is simulated here.
% 2. No approved JLE point estimate or one-dimensional interval is changed.
% 3. The program loads the corrected Figure 3 likelihood and the completed
%    adaptively refined Figure 4 profile surfaces.
% 4. Panel (f) is checked one additional time on a directly evaluated
%    half-step grid (spacing 0.003125) around the thin confidence ridge.
% 5. The two pieces are joined only if direct likelihood calculations show
%    that they are connected. No component is deleted or hidden.
% 6. All probability-boundary portions of a confidence region are drawn.
% 7. When p_p=1, p_a is structurally not identifiable. The blue vertical
%    line then represents the full maximum-likelihood set; an arbitrary
%    optimizer value of p_a is not presented as an estimate.
% 8. The physical canvas is 7.0 by 8.25 inches so that the text remains
%    readable when the figure is inserted at manuscript page width.

clc;
close all;

%% 1. User folder and required saved analyses

programFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(programFolder);
baseFolder = fullfile(projectRoot,'outputs');

if ~isfolder(baseFolder)
    mkdir(baseFolder);
end

if ~isfolder(baseFolder)
    error('The main manuscript folder was not found:\n%s',baseFolder);
end

figure3Name = 'Figure3_JLE_16Panel_PublicationReady_v2';
figure3Folder = fullfile(baseFolder,figure3Name);
figure3Cache = fullfile(figure3Folder, ...
    [figure3Name '_Analysis_Cache.mat']);
figure3Results = fullfile(figure3Folder,[figure3Name '_results.csv']);

sourceFigure4Name = 'Figure4_JLE_Detection_Profiles_Final';
sourceFigure4Folder = fullfile(baseFolder,sourceFigure4Name);
sourceFigure4Cache = fullfile(sourceFigure4Folder, ...
    [sourceFigure4Name '_Refined_Profile_Cache.mat']);
sourceFigure4Results = fullfile(sourceFigure4Folder, ...
    [sourceFigure4Name '_results.csv']);

outputName = 'Figure4_JLE_Detection_Profiles_ManuscriptFinal';
outputFolder = fullfile(baseFolder,outputName);

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

figureBase = fullfile(outputFolder,'Figure 4');
pdfFile = [figureBase '.pdf'];
pngFile = [figureBase '.png'];
figFile = [figureBase '.fig'];
resultsCSV = fullfile(outputFolder,[outputName '_results.csv']);
surfaceCSV = fullfile(outputFolder,[outputName '_refined_profile_surface.csv']);
contourCSV = fullfile(outputFolder,[outputName '_profile_contours.csv']);
refinedCache = fullfile(outputFolder,[outputName '_Profile_Cache.mat']);
readmeFile = fullfile(outputFolder, ...
    'README_Figure4_JLE_Detection_Profiles_ManuscriptFinal.txt');
logFile = fullfile(outputFolder,[outputName '_Command_Window_Log.txt']);

requiredFiles = {figure3Cache,figure3Results, ...
    sourceFigure4Cache,sourceFigure4Results};
for q = 1:numel(requiredFiles)
    if ~isfile(requiredFiles{q})
        error(['A required saved file was not found:\n%s\n\n' ...
            'Do not rerun a simulation. Confirm that the corrected Figure 3 ' ...
            'and completed final Figure 4 folders are inside the main ' ...
            'manuscript folder.'], ...
            requiredFiles{q});
    end
end

if isfile(logFile)
    delete(logFile);
end
diary(logFile);
diaryGuard = onCleanup(@() safe_diary_off()); %#ok<NASGU>
startTimer = tic;

fprintf('============================================================\n');
fprintf('STARTING THE MANUSCRIPT-FINAL JLE FIGURE 4\n');
fprintf('Horizontal axis: passive detection probability, p_p.\n');
fprintf('Vertical axis: active detection probability, p_a.\n');
fprintf('No cluster or auxiliary observation will be simulated.\n');
fprintf('No approved JLE point estimate will be changed.\n');
fprintf('No approved one-dimensional interval will be changed.\n');
fprintf('No spline interpolation or cosmetic smoothing will be used.\n');
fprintf(['Panel (f) will receive one localized direct half-step ' ...
    'likelihood check.\n']);
fprintf('Every other completed profile surface will be retained.\n');
fprintf('The figure will be redrawn at its intended manuscript size.\n\n');
fprintf('CORRECTED FIGURE 3 CACHE:\n%s\n\n',figure3Cache);
fprintf('COMPLETED REFINED FIGURE 4 CACHE:\n%s\n\n',sourceFigure4Cache);
fprintf('ALL MANUSCRIPT-FINAL FIGURE 4 MATERIALS WILL BE SAVED IN:\n%s\n', ...
    outputFolder);
fprintf('============================================================\n\n');

%% 2. Load and verify the corrected Figure 3 calculation

requiredFigure3Variables = { ...
    'cacheVersion','R_true','k_true','N_observed','N_auxiliary', ...
    'bounds','threshold2D','panelLabel','truePP','truePA','truePobs', ...
    'R_hat','k_hat','pp_hat','pa_hat','pobs_hat','minimumNLL', ...
    'dataCell','thetaHatCell'};

sourceInfo = whos('-file',figure3Cache);
sourceNames = {sourceInfo.name};
missing = requiredFigure3Variables( ...
    ~ismember(requiredFigure3Variables,sourceNames));
if ~isempty(missing)
    error('The corrected Figure 3 cache is missing %s.',missing{1});
end

S = load(figure3Cache,requiredFigure3Variables{:});
T3 = readtable(figure3Results);

numberPanels = 16;
if numel(S.panelLabel)~=numberPanels || ...
        numel(S.dataCell)~=numberPanels || ...
        numel(S.thetaHatCell)~=numberPanels
    error('The Figure 3 cache must contain 16 displayed panels.');
end
if S.N_observed~=5000 || S.N_auxiliary~=5000
    error('The Figure 3 sample sizes are not 5,000 and 5,000.');
end
if height(T3)~=numberPanels
    error('The Figure 3 result table does not contain 16 rows.');
end
if max(abs(T3.p_p_hat-S.pp_hat))>1e-10 || ...
        max(abs(T3.p_a_hat-S.pa_hat))>1e-10 || ...
        max(abs(T3.MinimumJointNLL-S.minimumNLL))>1e-7
    error('The corrected Figure 3 MAT and CSV files do not agree.');
end
if max(abs(S.truePobs-(S.truePP+(1-S.truePP).*S.truePA)))>1e-12
    error('The combined-detection identity failed.');
end
for panel = 1:numberPanels
    if sum(S.dataCell{panel}.counts)~=5000
        error('Panel %s does not contain 5,000 observed clusters.', ...
            S.panelLabel{panel});
    end
    if S.dataCell{panel}.Naux~=5000
        error('Panel %s does not contain 5,000 auxiliary trials.', ...
            S.panelLabel{panel});
    end
end

% This check prevents accidental use of the old local solution for panel c.
if S.minimumNLL(3)>5684.6852 || ...
        abs(S.pp_hat(3)-0.6489231)>2e-5 || ...
        abs(S.pa_hat(3)-1)>2e-8
    error(['The Figure 3 cache does not contain the verified boundary-aware ' ...
        'correction for panel (c). Run the correction program first.']);
end

fprintf('Corrected Figure 3 checks passed.\n');
fprintf('Displayed panels = 16.\n');
fprintf('Unique independent analyses = 13.\n');
fprintf('Observed clusters per unique analysis = 5000.\n');
fprintf('Independent auxiliary trials per unique analysis = 5000.\n');
fprintf('Corrected panel (c) NLL = %.10f.\n\n',S.minimumNLL(3));

%% 3. Load and verify the completed refined detection profiles

F = load(sourceFigure4Cache);
Told = readtable(sourceFigure4Results);

requiredSourceVariables = { ...
    'refinedCacheVersion','sourceCacheVersion', ...
    'ppFineCell','paFineCell','refinedDeltaCell', ...
    'directEvaluationMaskCell','refinementFactor', ...
    'connectedComponents','componentSizesCell', ...
    'newLikelihoodEvaluations','refinementRuntimeSeconds', ...
    'completedRefinement'};
for q = 1:numel(requiredSourceVariables)
    if ~isfield(F,requiredSourceVariables{q})
        error('The completed Figure 4 cache is missing %s.', ...
            requiredSourceVariables{q});
    end
end
if ~all(F.completedRefinement) || numel(F.refinedDeltaCell)~=16
    error('The completed Figure 4 profile calculation is incomplete.');
end
if height(Told)~=16
    error('The completed Figure 4 result table must contain 16 rows.');
end
if max(abs(Told.JLE_p_p_hat-S.pp_hat))>1e-10 || ...
        max(abs(Told.Raw_optimizer_p_a_hat-S.pa_hat))>1e-10
    error('The completed Figure 4 results do not match corrected Figure 3.');
end
for panel = 1:numberPanels
    expectedSize = [numel(F.paFineCell{panel}), ...
        numel(F.ppFineCell{panel})];
    if ~isequal(size(F.refinedDeltaCell{panel}),expectedSize) || ...
            ~isequal(size(F.directEvaluationMaskCell{panel}),expectedSize)
        error('The completed profile surface for panel %s has the wrong size.', ...
            S.panelLabel{panel});
    end
end

% Keep the already approved one-dimensional intervals in one small
% compatibility structure. They are read, checked and plotted; they are
% never recalculated by this program.
P.ppLower = Told.p_p_lower_1D95;
P.ppUpper = Told.p_p_upper_1D95;
P.paLower = Told.p_a_lower_1D95;
P.paUpper = Told.p_a_upper_1D95;

fprintf('Completed refined Figure 4 profile checks passed.\n');
fprintf('Saved panel (f) grid spacing = %.5f.\n', ...
    median(diff(F.paFineCell{6})));
fprintf('Saved panel (f) accepted components = %d.\n', ...
    F.connectedComponents(6));
fprintf(['Only the narrow neighborhood responsible for this component ' ...
    'count will be checked at the next half step.\n\n']);

%% 4. Start or resume direct adaptive refinement

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
end

if useParallel
    fprintf('Parallel processing is available.\n\n');
else
    fprintf('Parallel processing is unavailable; serial refinement will be used.\n\n');
end

manuscriptCacheVersion = ...
    'Figure4-JLE-manuscript-final-panel-f-half-step-v1-2026-08-25';
sourceCacheVersion = S.cacheVersion;
sourceFigure4CacheVersion = F.refinedCacheVersion;

% Begin with the completed, directly evaluated surfaces. Only panel (f)
% can change resolution in this program.
ppFineCell = F.ppFineCell;
paFineCell = F.paFineCell;
refinedDeltaCell = F.refinedDeltaCell;
directEvaluationMaskCell = F.directEvaluationMaskCell;
refinementFactor = F.refinementFactor;
connectedComponents = F.connectedComponents;
componentSizesCell = F.componentSizesCell;
newLikelihoodEvaluations = zeros(numberPanels,1);
refinementRuntimeSeconds = zeros(numberPanels,1);

% Mark every genuinely calculated point in the existing surfaces. The old
% cache stores new adaptive evaluations separately; its original 0.025
% grid points occur every refinementFactor entries.
evaluatedPointMaskCell = cell(numberPanels,1);
for panel = 1:numberPanels
    evaluatedPointMaskCell{panel} = ...
        directEvaluationMaskCell{panel};
    oldCoarseRows = 1:refinementFactor(panel): ...
        numel(paFineCell{panel});
    oldCoarseColumns = 1:refinementFactor(panel): ...
        numel(ppFineCell{panel});
    evaluatedPointMaskCell{panel}(oldCoarseRows,oldCoarseColumns) = true;
end

loadedManuscriptCache = false;
if isfile(refinedCache)
    C = load(refinedCache);
    cacheMatches = isfield(C,'manuscriptCacheVersion') && ...
        isfield(C,'sourceCacheVersion') && ...
        isfield(C,'sourceFigure4CacheVersion') && ...
        strcmp(C.manuscriptCacheVersion,manuscriptCacheVersion) && ...
        strcmp(C.sourceCacheVersion,sourceCacheVersion) && ...
        strcmp(C.sourceFigure4CacheVersion,sourceFigure4CacheVersion);
    if cacheMatches
        requiredCachedFields = {'ppFineCell','paFineCell', ...
            'refinedDeltaCell','directEvaluationMaskCell', ...
            'evaluatedPointMaskCell','refinementFactor', ...
            'connectedComponents','componentSizesCell', ...
            'newLikelihoodEvaluations','refinementRuntimeSeconds'};
        cacheComplete = true;
        for q = 1:numel(requiredCachedFields)
            cacheComplete = cacheComplete && ...
                isfield(C,requiredCachedFields{q});
        end
        if cacheComplete
            ppFineCell = C.ppFineCell;
            paFineCell = C.paFineCell;
            refinedDeltaCell = C.refinedDeltaCell;
            directEvaluationMaskCell = C.directEvaluationMaskCell;
            evaluatedPointMaskCell = C.evaluatedPointMaskCell;
            refinementFactor = C.refinementFactor;
            connectedComponents = C.connectedComponents;
            componentSizesCell = C.componentSizesCell;
            newLikelihoodEvaluations = C.newLikelihoodEvaluations;
            refinementRuntimeSeconds = C.refinementRuntimeSeconds;
            loadedManuscriptCache = true;
            fprintf(['A matching manuscript-final cache was found.\n' ...
                'No likelihood point needs to be recalculated.\n\n']);
        end
    end
end

if ~loadedManuscriptCache
    panel = 6;
    panelTimer = tic;
    fprintf(['(f) directly checking the narrow confidence ridge at ' ...
        'half the saved grid spacing...\n']);

    sourceSpacing = median(diff(paFineCell{panel}));
    targetSpacing = sourceSpacing/2;
    searchMargin = 0.20;

    [ppNew,paNew,deltaNew,newDirectMask,newEvaluatedMask, ...
        numberNew] = refine_existing_surface_one_level( ...
        ppFineCell{panel},paFineCell{panel}, ...
        refinedDeltaCell{panel},evaluatedPointMaskCell{panel}, ...
        S.threshold2D,searchMargin,S.thetaHatCell{panel}, ...
        S.minimumNLL(panel),S.dataCell{panel},S.bounds,useParallel);

    [numberComponents,componentSizes] = ...
        count_grid_components(deltaNew<=S.threshold2D);

    ppFineCell{panel} = ppNew;
    paFineCell{panel} = paNew;
    refinedDeltaCell{panel} = deltaNew;
    directEvaluationMaskCell{panel} = newDirectMask;
    evaluatedPointMaskCell{panel} = newEvaluatedMask;
    refinementFactor(panel) = 2*F.refinementFactor(panel);
    connectedComponents(panel) = numberComponents;
    componentSizesCell{panel} = componentSizes;
    newLikelihoodEvaluations(panel) = numberNew;
    refinementRuntimeSeconds(panel) = toc(panelTimer);

    if abs(median(diff(paNew))-targetSpacing)>1e-12
        error('The panel (f) half-step grid spacing is incorrect.');
    end

    fprintf('    Previous spacing = %.6f.\n',sourceSpacing);
    fprintf('    Checked spacing  = %.6f.\n',targetSpacing);
    fprintf('    New direct likelihood evaluations = %d.\n',numberNew);
    fprintf('    Accepted components after checking = %d.\n', ...
        numberComponents);
    if numberComponents==1
        fprintf(['    RESULT: the tiny second piece was a grid-resolution ' ...
            'artifact and the directly checked region is connected.\n']);
    else
        fprintf(['    RESULT: %d pieces remain after direct half-step ' ...
            'checking; every piece will be retained and displayed.\n'], ...
            numberComponents);
    end
    fprintf('    Runtime = %.1f seconds.\n\n', ...
        refinementRuntimeSeconds(panel));

    save_manuscript_checkpoint(refinedCache, ...
        manuscriptCacheVersion,sourceCacheVersion, ...
        sourceFigure4CacheVersion,ppFineCell,paFineCell, ...
        refinedDeltaCell,directEvaluationMaskCell, ...
        evaluatedPointMaskCell,refinementFactor, ...
        connectedComponents,componentSizesCell, ...
        newLikelihoodEvaluations,refinementRuntimeSeconds);
end

fprintf('All manuscript-final detection-profile surfaces are ready.\n');
fprintf('New likelihood points evaluated by this program = %d.\n', ...
    sum(newLikelihoodEvaluations));
fprintf('Panel (f) accepted components = %d.\n', ...
    connectedComponents(6));
fprintf('Largest accepted-component count = %d.\n\n', ...
    max(connectedComponents));

%% 5. Verify the statistical interpretation retained from the saved analysis

identifiablePA = S.truePP<1-1e-12;
truePairDelta = Told.True_pair_profile_delta_NLL;
truthCovered = truePairDelta<=S.threshold2D+0.02;

if any(~truthCovered)
    error('At least one generating detection pair is outside its joint region.');
end
if any(S.pp_hat<P.ppLower-1e-6 | S.pp_hat>P.ppUpper+1e-6)
    error('At least one p_p estimate lies outside its saved 1-D interval.');
end
if any(S.pa_hat(identifiablePA)<P.paLower(identifiablePA)-1e-6 | ...
        S.pa_hat(identifiablePA)>P.paUpper(identifiablePA)+1e-6)
    error('At least one identifiable p_a estimate lies outside its interval.');
end

ppMAE = mean(abs(S.pp_hat-S.truePP));
paMAE = mean(abs(S.pa_hat(identifiablePA)-S.truePA(identifiablePA)));
pobsMAE = mean(abs(S.pobs_hat-S.truePobs));

fprintf('FINAL SCIENTIFIC CHECKS\n');
fprintf('Mean absolute error for p_p = %.6f.\n',ppMAE);
fprintf('Mean absolute error for identifiable p_a = %.6f.\n',paMAE);
fprintf('Mean absolute error for combined detection = %.6f.\n',pobsMAE);
fprintf('Joint detection regions containing truth = %d of 16.\n', ...
    sum(truthCovered));
fprintf('Structurally non-identifiable p_a panels = 4.\n\n');

%% 6. Save final numerical tables

Tfinal = Told;
Tfinal.RefinementFactor = refinementFactor;
Tfinal.AcceptedGridComponents = connectedComponents;
Tfinal.AdditionalDirectLikelihoodEvaluations = ...
    newLikelihoodEvaluations;
Tfinal.AdditionalRefinementRuntimeSeconds = ...
    refinementRuntimeSeconds;
Tfinal.FinalGridSpacing = zeros(numberPanels,1);
for panel = 1:numberPanels
    Tfinal.FinalGridSpacing(panel) = ...
        median(diff(ppFineCell{panel}));
end
writetable(Tfinal,resultsCSV);

surfacePanel = cell(0,1);
surfacePP = zeros(0,1);
surfacePA = zeros(0,1);
surfaceDelta = zeros(0,1);
surfaceInside = false(0,1);
surfaceDirect = false(0,1);
surfacePreviouslyEvaluated = false(0,1);

contourPanel = cell(0,1);
contourType = cell(0,1);
contourSegment = zeros(0,1);
contourPoint = zeros(0,1);
contourPP = zeros(0,1);
contourPA = zeros(0,1);

for panel = 1:numberPanels
    ppFine = ppFineCell{panel};
    paFine = paFineCell{panel};
    deltaFine = refinedDeltaCell{panel};
    if panel==6
        directMask = directEvaluationMaskCell{panel};
    else
        directMask = false(size(deltaFine));
    end
    evaluatedMask = evaluatedPointMaskCell{panel};
    [PPmesh,PAmesh] = meshgrid(ppFine,paFine);

    % Save only genuinely evaluated points. The high finite values placed
    % outside the directly checked corridor are deliberately omitted.
    keep = evaluatedMask;
    q = sum(keep(:));
    surfacePanel = [surfacePanel; ...
        repmat(S.panelLabel(panel),q,1)]; %#ok<AGROW>
    surfacePP = [surfacePP;PPmesh(keep)]; %#ok<AGROW>
    surfacePA = [surfacePA;PAmesh(keep)]; %#ok<AGROW>
    surfaceDelta = [surfaceDelta;deltaFine(keep)]; %#ok<AGROW>
    surfaceInside = [surfaceInside; ...
        deltaFine(keep)<=S.threshold2D]; %#ok<AGROW>
    surfaceDirect = [surfaceDirect;directMask(keep)]; %#ok<AGROW>
    surfacePreviouslyEvaluated = [surfacePreviouslyEvaluated; ...
        evaluatedMask(keep) & ~directMask(keep)]; %#ok<AGROW>

    segments = extract_all_region_boundaries( ...
        ppFine,paFine,deltaFine,S.threshold2D);
    for s = 1:numel(segments)
        q = numel(segments{s}.x);
        contourPanel = [contourPanel; ...
            repmat(S.panelLabel(panel),q,1)]; %#ok<AGROW>
        contourType = [contourType; ...
            repmat({segments{s}.type},q,1)]; %#ok<AGROW>
        contourSegment = [contourSegment;repmat(s,q,1)]; %#ok<AGROW>
        contourPoint = [contourPoint;(1:q)']; %#ok<AGROW>
        contourPP = [contourPP;segments{s}.x(:)]; %#ok<AGROW>
        contourPA = [contourPA;segments{s}.y(:)]; %#ok<AGROW>
    end
end

surfaceTable = table(surfacePanel,surfacePP,surfacePA,surfaceDelta, ...
    surfaceInside,surfacePreviouslyEvaluated,surfaceDirect, ...
    'VariableNames',{'Panel','p_p','p_a','ProfileDeltaNLL', ...
    'InsideJoint95Region','PreviouslyEvaluatedPoint', ...
    'NewPanelFHalfStepEvaluation'});
writetable(surfaceTable,surfaceCSV);

contourTable = table(contourPanel,contourType,contourSegment, ...
    contourPoint,contourPP,contourPA, ...
    'VariableNames',{'Panel','BoundaryType','ContourSegment', ...
    'ContourPoint','p_p','p_a'});
writetable(contourTable,contourCSV);

%% 7. Draw the clean manuscript-sized Figure-3-style 16-panel figure

fprintf('Drawing the final Figure 4...\n');

canvasWidth = 7.00;
canvasHeight = 8.25;
fig = figure('Color','w','Name',outputName, ...
    'Units','inches','Position',[0.25 0.25 canvasWidth canvasHeight], ...
    'Renderer','painters','InvertHardcopy','on');

leftMargin = 0.105;
rightMargin = 0.025;
bottomMargin = 0.070;
topOfPanels = 0.840;
horizontalGap = 0.025;
verticalGap = 0.035;
panelWidth = (1-leftMargin-rightMargin-3*horizontalGap)/4;
panelHeight = (topOfPanels-bottomMargin-3*verticalGap)/4;

legendHandles = gobjects(1,4);

for panel = 1:numberPanels
    row = ceil(panel/4);
    column = mod(panel-1,4)+1;
    xPosition = leftMargin+(column-1)*(panelWidth+horizontalGap);
    yPosition = topOfPanels-row*panelHeight-(row-1)*verticalGap;

    ax = axes('Parent',fig,'Position', ...
        [xPosition yPosition panelWidth panelHeight]);
    hold(ax,'on');

    ppFine = ppFineCell{panel};
    paFine = paFineCell{panel};
    deltaFine = refinedDeltaCell{panel};
    [PPmesh,PAmesh] = meshgrid(ppFine,paFine);

    [~,hJoint] = contour(ax,PPmesh,PAmesh,deltaFine, ...
        [S.threshold2D S.threshold2D],'k-','LineWidth',1.45);
    draw_probability_boundary_segments(ax,ppFine,paFine, ...
        deltaFine<=S.threshold2D,[0 0 0],1.45);

    if identifiablePA(panel)
        hIntervals = plot(ax,[P.ppLower(panel) P.ppUpper(panel)], ...
            [S.pa_hat(panel) S.pa_hat(panel)],'--', ...
            'Color',[0.02 0.34 0.96],'LineWidth',1.20);
        plot(ax,[S.pp_hat(panel) S.pp_hat(panel)], ...
            [P.paLower(panel) P.paUpper(panel)],'--', ...
            'Color',[0.02 0.34 0.96],'LineWidth',1.20);
        hEstimate = plot(ax,S.pp_hat(panel),S.pa_hat(panel),'o', ...
            'MarkerSize',5.4,'MarkerFaceColor',[0.02 0.36 0.95], ...
            'MarkerEdgeColor',[0.01 0.26 0.78],'LineWidth',0.9);
    else
        % p_a is not identifiable when p_p=1. Show the complete MLE set.
        hEstimate = plot(ax,[1 1],[0 1],'-', ...
            'Color',[0.02 0.36 0.95],'LineWidth',1.8);
        intervalY = 0.10;
        hIntervals = plot(ax,[P.ppLower(panel) P.ppUpper(panel)], ...
            [intervalY intervalY],'--', ...
            'Color',[0.02 0.34 0.96],'LineWidth',1.20);
        plot(ax,[P.ppLower(panel) P.ppLower(panel)], ...
            [intervalY-0.025 intervalY+0.025],'-', ...
            'Color',[0.02 0.34 0.96],'LineWidth',1.0);
        text(ax,0.50,0.50,'p_a not identifiable', ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize',8.0,'FontWeight','bold','Interpreter','tex', ...
            'BackgroundColor','w');
        text(ax,0.69,intervalY+0.065,'p_p interval', ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',7.6,'Interpreter','tex');
    end

    hTruth = plot(ax,S.truePP(panel),S.truePA(panel),'d', ...
        'MarkerSize',5.8,'MarkerFaceColor',[1.00 0.08 0.06], ...
        'MarkerEdgeColor',[0.90 0.03 0.03],'LineWidth',0.9);

    if panel==1
        legendHandles = [hJoint hEstimate hTruth hIntervals];
    end

    xlim(ax,[-0.035 1.035]);
    ylim(ax,[-0.035 1.035]);
    xticks(ax,[0 0.25 0.50 0.75 1.00]);
    yticks(ax,[0 0.25 0.50 0.75 1.00]);
    xticklabels(ax,{'0','0.25','0.50','0.75','1'});
    yticklabels(ax,{'0','0.25','0.50','0.75','1'});
    pbaspect(ax,[1 1 1]);
    grid(ax,'on');
    ax.GridAlpha = 0.15;
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

    title(ax,sprintf('(%c)  p_p=%.2f, p_a=%.2f', ...
        char(96+panel),S.truePP(panel),S.truePA(panel)), ...
        'FontSize',8.2,'FontWeight','bold','Interpreter','tex');
end

annotation(fig,'textbox',[0.04 0.955 0.92 0.030], ...
    'String','JLE detection estimates', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',14.5,'FontWeight','bold','EdgeColor','none', ...
    'Interpreter','tex');

legendObject = legend(legendHandles, ...
    {'95% joint region','JLE estimate or MLE set', ...
    'Truth','Separate 95% 1-D intervals'}, ...
    'Orientation','horizontal','Box','off','FontSize',8.5);
legendObject.Position = [0.145 0.900 0.71 0.031];

annotation(fig,'textbox',[0.29 0.012 0.42 0.030], ...
    'String','Passive detection probability, p_p', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',11.5,'FontWeight','bold','EdgeColor','none', ...
    'Interpreter','tex');

labelAxes = axes('Parent',fig,'Position',[0 0 1 1], ...
    'Visible','off','HitTest','off');
text(labelAxes,0.025,0.46,'Active detection probability, p_a', ...
    'Rotation',90,'HorizontalAlignment','center', ...
    'VerticalAlignment','middle','FontSize',11.5, ...
    'FontWeight','bold','Interpreter','tex','Visible','on');

drawnow;

%% 8. Save every final output

fprintf('Saving the final Figure 4...\n');

savefig(fig,figFile);
set(fig,'PaperUnits','inches');
set(fig,'PaperPosition',[0 0 canvasWidth canvasHeight]);
set(fig,'PaperSize',[canvasWidth canvasHeight]);
set(fig,'Color','w','InvertHardcopy','on');
print(fig,pngFile,'-dpng','-r360');
print(fig,pdfFile,'-dpdf','-painters');
close(fig);

sourceCodeFile = [mfilename('fullpath') '.m'];
savedCodeFile = fullfile(outputFolder,[mfilename '.m']);
if isfile(sourceCodeFile) && ~strcmpi(sourceCodeFile,savedCodeFile)
    copyfile(sourceCodeFile,savedCodeFile,'f');
end

fid = fopen(readmeFile,'w');
if fid~=-1
    fprintf(fid,'MANUSCRIPT-FINAL JLE DETECTION-PROFILE FIGURE 4\r\n\r\n');
    fprintf(fid,'No cluster or auxiliary observation was simulated.\r\n');
    fprintf(fid,'No approved point estimate or 1-D interval was changed.\r\n');
    fprintf(fid,['The completed direct profile surfaces were retained. ' ...
        'Panel (f) alone was checked at half its previous spacing.\r\n']);
    fprintf(fid,['The small pieces were joined only if the direct ' ...
        'likelihood evaluations supported a connection.\r\n']);
    fprintf(fid,'No spline interpolation or cosmetic smoothing was used.\r\n');
    fprintf(fid,['When p_p=1, p_a is structurally not identifiable and the ' ...
        'blue vertical line is the maximum-likelihood set.\r\n\r\n']);
    fprintf(fid,'Corrected Figure 3 cache: %s\r\n',figure3Cache);
    fprintf(fid,'Source Figure 4 cache: %s\r\n',sourceFigure4Cache);
    fprintf(fid,'Joint regions containing truth: %d of 16\r\n', ...
        sum(truthCovered));
    fprintf(fid,'Mean absolute error for p_p: %.6f\r\n',ppMAE);
    fprintf(fid,'Mean absolute error for identifiable p_a: %.6f\r\n',paMAE);
    fprintf(fid,'Mean absolute error for combined detection: %.6f\r\n', ...
        pobsMAE);
    fprintf(fid,'New likelihood evaluations: %d\r\n', ...
        sum(newLikelihoodEvaluations));
    fprintf(fid,'Panel (f) accepted components: %d\r\n', ...
        connectedComponents(6));
    fprintf(fid,'Physical canvas: %.2f by %.2f inches\r\n', ...
        canvasWidth,canvasHeight);
    fclose(fid);
end

fprintf('\n============================================================\n');
fprintf('DONE: THE MANUSCRIPT-FINAL JLE FIGURE 4 WAS CREATED.\n');
fprintf('Elapsed time for this run: %.1f seconds.\n',toc(startTimer));
fprintf('No new cluster or auxiliary data were simulated.\n');
fprintf('No approved JLE point estimate was changed.\n');
fprintf('No approved 1-D interval was changed.\n');
fprintf('No spline interpolation or cosmetic smoothing was used.\n');
fprintf('Direct adaptive likelihood refinement was used: YES.\n');
fprintf('Panel (f) accepted components = %d.\n', ...
    connectedComponents(6));
fprintf('Physical figure size = %.2f by %.2f inches.\n', ...
    canvasWidth,canvasHeight);
fprintf('Joint detection regions containing truth = %d of 16.\n', ...
    sum(truthCovered));
fprintf('Mean absolute error for p_p = %.6f.\n',ppMAE);
fprintf('Mean absolute error for identifiable p_a = %.6f.\n',paMAE);
fprintf('Mean absolute error for combined detection = %.6f.\n',pobsMAE);
fprintf('Largest accepted-component count = %d.\n', ...
    max(connectedComponents));

fprintf('\nALL MANUSCRIPT-FINAL FIGURE 4 MATERIALS WERE SAVED IN:\n%s\n', ...
    outputFolder);
fprintf('\nPDF FIGURE:\n%s\n',pdfFile);
fprintf('\nPNG FIGURE:\n%s\n',pngFile);
fprintf('\nMATLAB FIG FILE:\n%s\n',figFile);
fprintf('\nFINAL RESULTS CSV:\n%s\n',resultsCSV);
fprintf('\nREFINED PROFILE-SURFACE CSV:\n%s\n',surfaceCSV);
fprintf('\nPROFILE-CONTOUR CSV:\n%s\n',contourCSV);
fprintf('\nREUSABLE REFINED PROFILE CACHE:\n%s\n',refinedCache);
fprintf('\nSAVED CODE COPY:\n%s\n',savedCodeFile);
fprintf('\nREADME:\n%s\n',readmeFile);
fprintf('\nCOMMAND-WINDOW LOG:\n%s\n',logFile);
fprintf('============================================================\n');

end

%% =====================================================================
function [xFine,yFine,deltaFine,newDirectMask,evaluatedMask, ...
    totalNewEvaluations] = refine_existing_surface_one_level( ...
    xSource,ySource,deltaSource,sourceEvaluatedMask,threshold, ...
    searchMargin,thetaHat,minNLL,data,bounds,useParallel)

% Split every saved interval into two. Values are copied only at points
% that were genuinely evaluated in the preceding analysis. Missing points
% close to the confidence region are then evaluated from the likelihood.
factor = 2;
xFine = subdivide_grid(xSource,factor);
yFine = subdivide_grid(ySource,factor);
numberRows = numel(yFine);
numberColumns = numel(xFine);

deltaFine = NaN(numberRows,numberColumns);
newDirectMask = false(numberRows,numberColumns);
inheritedMask = false(numberRows,numberColumns);

[sourceRows,sourceColumns] = find(sourceEvaluatedMask);
fineRows = factor*(sourceRows-1)+1;
fineColumns = factor*(sourceColumns-1)+1;
sourceLinear = sub2ind(size(deltaSource),sourceRows,sourceColumns);
fineLinear = sub2ind(size(deltaFine),fineRows,fineColumns);
deltaFine(fineLinear) = deltaSource(sourceLinear);
inheritedMask(fineLinear) = true;

% Direct calculation is requested around every saved point whose profile
% value is at or moderately above the 95% threshold. The extra dilation
% covers a ridge that passes between adjacent saved grid points.
sourceNear = deltaSource<=threshold+searchMargin;
sourceNear = sourceNear & sourceEvaluatedMask;
sourceNear = dilate_logical_mask(sourceNear,1);
refineMask = false(numberRows,numberColumns);
for r = 1:size(sourceNear,1)
    for c = 1:size(sourceNear,2)
        if sourceNear(r,c)
            rf = factor*(r-1)+1;
            cf = factor*(c-1)+1;
            rowRange = max(1,rf-factor):min(numberRows,rf+factor);
            columnRange = max(1,cf-factor):min(numberColumns,cf+factor);
            refineMask(rowRange,columnRange) = true;
        end
    end
end

totalNewEvaluations = 0;
maximumExpansionRounds = 12;

for roundNumber = 1:maximumExpansionRounds
    evaluatedMask = inheritedMask | newDirectMask;
    toEvaluate = refineMask & ~evaluatedMask;
    numberNew = sum(toEvaluate(:));

    if numberNew>0
        [rowIndex,columnIndex] = find(toEvaluate);
        ppValues = xFine(columnIndex);
        paValues = yFine(rowIndex);
        newValues = NaN(numberNew,1);

        if useParallel
            parfor q = 1:numberNew
                profileNLL = profile_nll_fixed_detection( ...
                    ppValues(q),paValues(q),thetaHat,data,bounds);
                newValues(q) = profileNLL-minNLL;
            end
        else
            for q = 1:numberNew
                profileNLL = profile_nll_fixed_detection( ...
                    ppValues(q),paValues(q),thetaHat,data,bounds);
                newValues(q) = profileNLL-minNLL;
            end
        end

        if min(newValues)<-0.03
            error(['The panel (f) half-step check found an NLL more than ' ...
                '0.03 below the corrected Figure 3 minimum.']);
        end
        newValues = max(newValues,0);
        newLinear = sub2ind(size(deltaFine),rowIndex,columnIndex);
        deltaFine(newLinear) = newValues;
        newDirectMask(newLinear) = true;
        totalNewEvaluations = totalNewEvaluations+numberNew;
    end

    evaluatedMask = inheritedMask | newDirectMask;
    accepted = evaluatedMask & deltaFine<=threshold;
    unevaluatedNeighbor = dilate_logical_mask(accepted,1) & ...
        ~evaluatedMask;
    if ~any(unevaluatedNeighbor(:))
        break;
    end
    refineMask = refineMask | ...
        dilate_logical_mask(unevaluatedNeighbor,2);
end

evaluatedMask = inheritedMask | newDirectMask;
accepted = evaluatedMask & deltaFine<=threshold;
remainingNeighbor = dilate_logical_mask(accepted,1) & ~evaluatedMask;
if roundNumber==maximumExpansionRounds && any(remainingNeighbor(:))
    error('The panel (f) direct refinement did not enclose the region.');
end

% Unevaluated points lie outside the enclosed threshold corridor. A high
% finite value lets MATLAB trace the boundary without creating a new low
% likelihood region. Such points are excluded from the numerical CSV.
deltaFine(~evaluatedMask) = threshold+searchMargin+20;

% Every accepted point must now be surrounded by evaluated points. This
% prevents a displayed boundary from being caused by the numerical fill.
accepted = deltaFine<=threshold;
remainingNeighbor = dilate_logical_mask(accepted,1) & ~evaluatedMask;
if any(remainingNeighbor(:))
    error('An accepted panel (f) point touches an unevaluated point.');
end
end

%% =====================================================================
function save_manuscript_checkpoint(cacheFile, ...
    manuscriptCacheVersion,sourceCacheVersion, ...
    sourceFigure4CacheVersion,ppFineCell,paFineCell, ...
    refinedDeltaCell,directEvaluationMaskCell,evaluatedPointMaskCell, ...
    refinementFactor,connectedComponents,componentSizesCell, ...
    newLikelihoodEvaluations,refinementRuntimeSeconds)

temporaryFile = [cacheFile '.temporary.mat'];
if isfile(temporaryFile)
    delete(temporaryFile);
end

save(temporaryFile, ...
    'manuscriptCacheVersion','sourceCacheVersion', ...
    'sourceFigure4CacheVersion','ppFineCell','paFineCell', ...
    'refinedDeltaCell','directEvaluationMaskCell', ...
    'evaluatedPointMaskCell','refinementFactor', ...
    'connectedComponents','componentSizesCell', ...
    'newLikelihoodEvaluations','refinementRuntimeSeconds','-v7.3');

test = load(temporaryFile,'manuscriptCacheVersion', ...
    'sourceFigure4CacheVersion','connectedComponents');
if ~strcmp(test.manuscriptCacheVersion,manuscriptCacheVersion) || ...
        ~strcmp(test.sourceFigure4CacheVersion, ...
        sourceFigure4CacheVersion) || ...
        numel(test.connectedComponents)~=16
    delete(temporaryFile);
    error('The temporary manuscript-final checkpoint failed verification.');
end
movefile(temporaryFile,cacheFile,'f');
end

function fineGrid = subdivide_grid(coarseGrid,factor)

fineGrid = zeros(1,(numel(coarseGrid)-1)*factor+1);
position = 1;
for q = 1:numel(coarseGrid)-1
    segment = linspace(coarseGrid(q),coarseGrid(q+1),factor+1);
    fineGrid(position:position+factor-1) = segment(1:end-1);
    position = position+factor;
end
fineGrid(end) = coarseGrid(end);
end

%% =====================================================================
function expanded = dilate_logical_mask(mask,radius)

expanded = mask;
for q = 1:radius
    padded = false(size(expanded)+2);
    padded(2:end-1,2:end-1) = expanded;
    next = false(size(expanded));
    for dr = 0:2
        for dc = 0:2
            next = next | padded(1+dr:end-2+dr,1+dc:end-2+dc);
        end
    end
    expanded = next;
end
end

%% =====================================================================
function [numberComponents,componentSizes] = count_grid_components(mask)

visited = false(size(mask));
componentSizes = zeros(0,1);
numberComponents = 0;
[numberRows,numberColumns] = size(mask);

for startRow = 1:numberRows
    for startColumn = 1:numberColumns
        if ~mask(startRow,startColumn) || visited(startRow,startColumn)
            continue;
        end

        numberComponents = numberComponents+1;
        queueRows = zeros(nnz(mask),1);
        queueColumns = zeros(nnz(mask),1);
        head = 1;
        tail = 1;
        queueRows(1) = startRow;
        queueColumns(1) = startColumn;
        visited(startRow,startColumn) = true;
        sizeOfComponent = 0;

        while head<=tail
            r = queueRows(head);
            c = queueColumns(head);
            head = head+1;
            sizeOfComponent = sizeOfComponent+1;

            for dr = -1:1
                for dc = -1:1
                    if dr==0 && dc==0
                        continue;
                    end
                    rr = r+dr;
                    cc = c+dc;
                    if rr>=1 && rr<=numberRows && ...
                            cc>=1 && cc<=numberColumns && ...
                            mask(rr,cc) && ~visited(rr,cc)
                        tail = tail+1;
                        queueRows(tail) = rr;
                        queueColumns(tail) = cc;
                        visited(rr,cc) = true;
                    end
                end
            end
        end
        componentSizes(numberComponents,1) = sizeOfComponent;
    end
end
end

%% =====================================================================
function draw_probability_boundary_segments(ax,x,y,inside,colorValue,lineWidth)

% A contour routine draws threshold crossings inside a grid but does not
% automatically draw the part of a region that lies on p=0 or p=1.
draw_runs(ax,x,repmat(y(1),size(x)),inside(1,:),colorValue,lineWidth);
draw_runs(ax,x,repmat(y(end),size(x)),inside(end,:),colorValue,lineWidth);
draw_runs(ax,repmat(x(1),size(y)),y,inside(:,1)',colorValue,lineWidth);
draw_runs(ax,repmat(x(end),size(y)),y,inside(:,end)',colorValue,lineWidth);
end

%% =====================================================================
function draw_runs(ax,x,y,mask,colorValue,lineWidth)

mask = logical(mask(:)');
edges = diff([false mask false]);
starts = find(edges==1);
stops = find(edges==-1)-1;
for q = 1:numel(starts)
    if stops(q)>starts(q)
        plot(ax,x(starts(q):stops(q)),y(starts(q):stops(q)),'-', ...
            'Color',colorValue,'LineWidth',lineWidth);
    end
end
end

%% =====================================================================
function segments = extract_all_region_boundaries(x,y,z,level)

segments = extract_contour_segments(x,y,z,level);
inside = z<=level;
segments = append_boundary_runs(segments,x,repmat(y(1),size(x)), ...
    inside(1,:),'p_a_lower_boundary');
segments = append_boundary_runs(segments,x,repmat(y(end),size(x)), ...
    inside(end,:),'p_a_upper_boundary');
segments = append_boundary_runs(segments,repmat(x(1),size(y)),y, ...
    inside(:,1)','p_p_lower_boundary');
segments = append_boundary_runs(segments,repmat(x(end),size(y)),y, ...
    inside(:,end)','p_p_upper_boundary');
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
    segment.type = 'likelihood_threshold';
    segments{end+1,1} = segment; %#ok<AGROW>
    column = column+numberPoints+1;
end
end

%% =====================================================================
function segments = append_boundary_runs(segments,x,y,mask,typeName)

mask = logical(mask(:)');
edges = diff([false mask false]);
starts = find(edges==1);
stops = find(edges==-1)-1;
for q = 1:numel(starts)
    if stops(q)>starts(q)
        segment.x = x(starts(q):stops(q));
        segment.y = y(starts(q):stops(q));
        segment.type = typeName;
        segments{end+1,1} = segment; %#ok<AGROW>
    end
end
end

%% =====================================================================
function bestValue = profile_nll_fixed_detection( ...
    p_p,p_a,thetaHat,data,bounds)

persistent optimizationOptions
if isempty(optimizationOptions)
    optimizationOptions = optimoptions( ...
        'fmincon','Algorithm','sqp','Display','off', ...
        'MaxFunctionEvaluations',700,'MaxIterations',190, ...
        'OptimalityTolerance',2e-7,'StepTolerance',2e-8);
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
for q = 1:size(starts,1)
    objective = @(transmission) evaluate_fixed_detection( ...
        transmission,template,data,bounds);
    try
        [~,value] = fmincon(objective,starts(q,:),[],[],[],[], ...
            lower,upper,[],optimizationOptions);
        bestValue = min(bestValue,value);
    catch
    end
end
if ~isfinite(bestValue)
    bestValue = 1e100;
end
end

%% =====================================================================
function value = evaluate_fixed_detection(transmission,template,data,bounds)

theta = template;
theta(1:2) = transmission;
value = JLE_nll(theta,data,bounds);
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
[positiveSelected,pZero] = selected_observation_probabilities( ...
    pS,p_p,p_a,data);
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

logLikelihood = gammaln(N+1)-gammaln(Y+1)-gammaln(N-Y+1);
if Y>0
    logLikelihood = logLikelihood+Y*log(p);
end
if N-Y>0
    logLikelihood = logLikelihood+(N-Y)*log1p(-p);
end
end

%% =====================================================================
function [positiveSelected,pZero] = selected_observation_probabilities( ...
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
    logTerm1 = logChoose+o*log(z)+remaining*log1p(-z);
    term1 = exp(logTerm1);
end

if p_p>=1 || p_a<=0
    term2 = zeros(size(logChoose));
elseif p_a>=1
    term2 = (1-p_p).^s.*double(s==o);
else
    logTerm2 = s*log1p(-p_p)+logChoose+ ...
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
logProbability = -log(j)+gammaln(j*k+j-1)-gammaln(j*k)- ...
    gammaln(j)+j*k*log(k/(k+R))+(j-1)*log(R/(k+R));
pS = exp(logProbability);
pS = pS(:)';
end

function safe_diary_off

try
    diary off;
catch
end
end

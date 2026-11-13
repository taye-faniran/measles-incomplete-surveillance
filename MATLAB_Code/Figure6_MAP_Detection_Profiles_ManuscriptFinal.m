function Figure6_MAP_Detection_Profiles_ManuscriptFinal
% FIGURE6_MAP_DETECTION_PROFILES_MANUSCRIPTFINAL
% -------------------------------------------------------------------------
% Creates the manuscript-final MAP detection figure with passive detection
% probability p_p on the horizontal axis and active detection probability
% p_a on the vertical axis.
%
% SCIENTIFIC RULES USED BY THIS PROGRAM
%   1. The exact saved Figure 5/6 simulated datasets are reused.
%   2. No cluster is simulated.
%   3. The saved four-parameter MAP point estimates are not changed.
%   4. Separate 95% profile-posterior intervals are calculated for p_p and
%      p_a by profiling over the other three parameters.
%   5. The joint 95% p_p-p_a region is calculated by profiling over
%      R_eff and k.
%   6. The same likelihood-ratio thresholds used for Figure 5 are used:
%         2 Delta NLP <= chi2inv(0.95,1) for a separate 1-D interval;
%         2 Delta NLP <= chi2inv(0.95,2) for a joint 2-D region.
%   7. Completed panels are checkpointed. A later matching run loads the
%      cache and redraws the figure quickly.
%   8. If profiling finds an objective below a saved MAP minimum by more
%      than the numerical audit tolerance, the program stops instead of
%      silently changing the saved analysis.
%
% REQUIRED TOOLBOXES
%   Optimization Toolbox
%   Statistics and Machine Learning Toolbox
%
% SAVE THIS FILE AS
%   Figure6_MAP_Detection_Profiles_ManuscriptFinal.m
% -------------------------------------------------------------------------

clc;
close all;

programTimer = tic;

programName = 'Figure6_MAP_Detection_Profiles_ManuscriptFinal';
outputFolderName = 'Figure6_MAP_Detection_Profiles_ManuscriptFinal';

programFolder = fileparts(mfilename('fullpath'));

if isempty(programFolder)
    programFolder = pwd;
end

repositoryRoot = fileparts(programFolder);
projectRoot = repositoryRoot;
outputFolder = fullfile(repositoryRoot,'outputs',outputFolderName);

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

logFile = fullfile(outputFolder, ...
    'Figure6_MAP_Detection_Profiles_ManuscriptFinal_Command_Window_Log.txt');

diary(logFile);
diaryCleanup = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('\n============================================================\n');
fprintf('STARTING THE MANUSCRIPT-FINAL MAP FIGURE 6\n');
fprintf('Horizontal axis: passive detection probability, p_p.\n');
fprintf('Vertical axis: active detection probability, p_a.\n');
fprintf('No cluster will be simulated.\n');
fprintf('No saved MAP point estimate will be changed.\n');
fprintf('New detection profile-posterior uncertainty will be calculated.\n');
fprintf('The first run is checkpointed panel by panel.\n');
fprintf('Later matching runs load the completed cache and redraw quickly.\n');
fprintf('\nPROJECT FOLDER:\n%s\n',projectRoot);
fprintf('\nALL FIGURE 6 MATERIALS WILL BE SAVED IN:\n%s\n', ...
    outputFolder);
fprintf('============================================================\n\n');

verify_required_toolboxes();

[sourceMat,observedCountsCSV,baselineCSV] = ...
    locate_map_sources(projectRoot,outputFolder);

fprintf('SAVED MAP ANALYSIS OBJECT:\n%s\n\n',sourceMat);
fprintf('EXACT OBSERVED-CLUSTER COUNT FILE:\n%s\n\n',observedCountsCSV);

if ~isempty(baselineCSV)
    fprintf('MAP SUMMARY FILE USED FOR AN INDEPENDENT CROSS-CHECK:\n%s\n\n', ...
        baselineCSV);
end

saved = load(sourceMat,'cfg','baseline');

if ~isfield(saved,'cfg') || ~isfield(saved,'baseline')
    error(['The saved MAT file must contain variables named cfg and ' ...
        'baseline: %s'],sourceMat);
end

cfg = saved.cfg;
baseline = saved.baseline(:);

cfg = add_figure6_configuration(cfg);

countsByPanel = read_exact_observed_counts( ...
    observedCountsCSV,baseline);

verify_saved_map_inputs( ...
    baseline,countsByPanel,cfg,baselineCSV);

fprintf('Saved MAP input checks passed.\n');
fprintf('Displayed surveillance settings = 16.\n');
fprintf('One fixed population of true clusters = %d.\n', ...
    cfg.nTrueClusters);
fprintf('Observed clusters range from %d to %d.\n', ...
    min([baseline.nObserved]),max([baseline.nObserved]));
fprintf('The MAP estimates and original Figure 5 intervals are unchanged.\n\n');

sourceSignature = make_source_signature( ...
    sourceMat,observedCountsCSV,baseline,cfg);

cacheFile = fullfile(outputFolder, ...
    'Figure6_MAP_Detection_Profiles_ManuscriptFinal_Profile_Cache.mat');

profileResults = repmat(empty_profile_result(),16,1);

if exist(cacheFile,'file')
    cached = load(cacheFile);

    if isfield(cached,'sourceSignature') && ...
            strcmp(cached.sourceSignature,sourceSignature) && ...
            isfield(cached,'profileResults') && ...
            numel(cached.profileResults)==16

        profileResults = cached.profileResults(:);

        fprintf('A matching Figure 6 profile cache was found.\n');
        fprintf('Completed panels in the cache = %d of 16.\n\n', ...
            sum([profileResults.complete]));
    else
        fprintf(['An older or nonmatching Figure 6 cache was found. ' ...
            'It will not be used.\n\n']);
    end
end

if any(~[profileResults.complete])
    maybe_start_parallel_pool(cfg);
else
    fprintf(['All 16 detection profiles are already complete. ' ...
        'No parallel pool is needed for this redraw.\n\n']);
end

for panelIndex = 1:16

    if profileResults(panelIndex).complete
        fprintf('(%s) detection profile already complete; loaded from cache.\n', ...
            baseline(panelIndex).panel);
        continue;
    end

    panelTimer = tic;
    result = baseline(panelIndex);
    counts = countsByPanel{panelIndex};
    kernel = make_likelihood_kernel(counts,cfg);

    fprintf(['(%s) profiling detection for true p_p=%.2f, p_a=%.2f ' ...
        '(%d of 16)...\n'], ...
        result.panel,result.truePP,result.truePA,panelIndex);

    recomputedMinimum = negative_log_posterior( ...
        result.mapTheta,counts,kernel,result.prior,cfg);

    minimumTolerance = objective_match_tolerance(result.minimumObjective);

    if abs(recomputedMinimum-result.minimumObjective)>minimumTolerance
        error(['Panel (%s) does not reproduce its saved minimum.\n' ...
            'Saved negative log posterior = %.12f.\n' ...
            'Recomputed value = %.12f.'], ...
            result.panel,result.minimumObjective,recomputedMinimum);
    end

    [ppCI,ppLowerAtBound,ppUpperAtBound,ppThresholdError] = ...
        direct_profile_interval( ...
        3,result.mapTheta,result.minimumObjective, ...
        counts,kernel,result.prior,cfg);

    [paCI,paLowerAtBound,paUpperAtBound,paThresholdError] = ...
        direct_profile_interval( ...
        4,result.mapTheta,result.minimumObjective, ...
        counts,kernel,result.prior,cfg);

    contour = direct_joint_detection_contour( ...
        result.mapTheta,result.minimumObjective,ppCI,paCI, ...
        counts,kernel,result.prior,cfg);

    evaluatedTruePP = min(max(result.truePP,cfg.lowerBounds(3)), ...
        cfg.upperBounds(3));
    evaluatedTruePA = min(max(result.truePA,cfg.lowerBounds(4)), ...
        cfg.upperBounds(4));

    truePairProfileValue = profile_at_detection_pair( ...
        evaluatedTruePP,evaluatedTruePA,result.mapTheta, ...
        result.minimumObjective,counts,kernel,result.prior,cfg);

    truePairDelta = max(0,truePairProfileValue-result.minimumObjective);
    jointThresholdIncrease = 0.5*chi2inv(0.95,2);

    truePairInsideJoint95 = ...
        truePairDelta <= jointThresholdIncrease+cfg.figure6AuditTolerance;

    truePPInsideSeparate95 = ...
        evaluatedTruePP>=ppCI(1)-cfg.figure6BoundaryTolerance && ...
        evaluatedTruePP<=ppCI(2)+cfg.figure6BoundaryTolerance;

    truePAInsideSeparate95 = ...
        evaluatedTruePA>=paCI(1)-cfg.figure6BoundaryTolerance && ...
        evaluatedTruePA<=paCI(2)+cfg.figure6BoundaryTolerance;

    panelProfile = empty_profile_result();
    panelProfile.complete = true;
    panelProfile.panel = result.panel;
    panelProfile.truePP = result.truePP;
    panelProfile.truePA = result.truePA;
    panelProfile.evaluatedTruePP = evaluatedTruePP;
    panelProfile.evaluatedTruePA = evaluatedTruePA;
    panelProfile.nObserved = sum(counts);
    panelProfile.mapTheta = result.mapTheta;
    panelProfile.minimumObjective = result.minimumObjective;
    panelProfile.ppCI = ppCI;
    panelProfile.paCI = paCI;
    panelProfile.ppLowerAtBound = ppLowerAtBound;
    panelProfile.ppUpperAtBound = ppUpperAtBound;
    panelProfile.paLowerAtBound = paLowerAtBound;
    panelProfile.paUpperAtBound = paUpperAtBound;
    panelProfile.profileThresholdError = ...
        max(ppThresholdError,paThresholdError);
    panelProfile.contourPP = contour.pp;
    panelProfile.contourPA = contour.pa;
    panelProfile.contourThresholdError = contour.maxThresholdError;
    panelProfile.contourTouchesBoundary = contour.touchesBoundary;
    panelProfile.contourBoundaryShortfall = contour.maxBoundaryShortfall;
    panelProfile.truePairProfileDelta = truePairDelta;
    panelProfile.truePairInsideJoint95 = truePairInsideJoint95;
    panelProfile.truePPInsideSeparate95 = truePPInsideSeparate95;
    panelProfile.truePAInsideSeparate95 = truePAInsideSeparate95;
    panelProfile.elapsedSeconds = toc(panelTimer);

    profileResults(panelIndex) = panelProfile;

    save_profile_cache_atomic( ...
        cacheFile,profileResults,sourceSignature,sourceMat, ...
        observedCountsCSV,cfg);

    fprintf('    MAP p_p=%.4f; 95%% profile interval [%.4f, %.4f].\n', ...
        result.mapTheta(3),ppCI(1),ppCI(2));
    fprintf('    MAP p_a=%.4f; 95%% profile interval [%.4f, %.4f].\n', ...
        result.mapTheta(4),paCI(1),paCI(2));
    fprintf('    True-pair profile delta NLP = %.5f.\n',truePairDelta);
    fprintf('    Joint 95%% region contains truth: %s.\n', ...
        yes_no(truePairInsideJoint95));
    fprintf('    Joint region touches a probability bound: %s.\n', ...
        yes_no(contour.touchesBoundary));
    fprintf('    Panel runtime = %.1f seconds.\n', ...
        panelProfile.elapsedSeconds);
    fprintf('    Checkpoint saved.\n\n');
end

if ~all([profileResults.complete])
    error('At least one Figure 6 detection profile is incomplete.');
end

quality = summarize_profile_quality(profileResults,baseline,cfg);

fprintf('\nFINAL SCIENTIFIC CHECKS\n');
fprintf('Mean absolute error for p_p = %.6f.\n',quality.ppMAE);
fprintf(['Mean absolute error for p_a, excluding the four generating ' ...
    'p_p=1 settings = %.6f.\n'],quality.paMAEIdentifiable);
fprintf('Separate p_p intervals containing truth = %d of 16.\n', ...
    quality.ppSeparateContainsTruth);
fprintf(['Separate p_a intervals containing truth among data-identifiable ' ...
    'settings = %d of 12.\n'],quality.paSeparateContainsTruthIdentifiable);
fprintf('Joint p_p-p_a regions containing truth = %d of 16.\n', ...
    quality.jointContainsTruth);
fprintf('Joint regions touching a probability bound = %d of 16.\n', ...
    quality.jointTouchesBoundary);
fprintf('Largest non-boundary contour threshold error = %.8g NLP units.\n', ...
    quality.maxContourThresholdError);
fprintf('Largest 1-D threshold error = %.8g NLP units.\n\n', ...
    quality.maxProfileThresholdError);

write_results_csv(profileResults,outputFolder);
write_contour_csv(profileResults,outputFolder);

allResultsFile = fullfile(outputFolder, ...
    'Figure6_MAP_Detection_Profiles_ManuscriptFinal_All_Results.mat');

save(allResultsFile,'cfg','baseline','profileResults','quality', ...
    'sourceSignature','sourceMat','observedCountsCSV','-v7.3');

fprintf('Drawing the manuscript-final Figure 6...\n');

fig = draw_figure6(profileResults);

fprintf('Saving the manuscript-final Figure 6...\n');

save_publication_figure(fig,outputFolder,'Figure 6');
close(fig);

copy_running_code(outputFolder,programName);
write_readme(outputFolder,sourceMat,observedCountsCSV,quality,cfg);

fprintf('\n============================================================\n');
fprintf('DONE: THE MANUSCRIPT-FINAL MAP FIGURE 6 WAS CREATED.\n');
fprintf('Elapsed time for this run: %.1f seconds.\n',toc(programTimer));
fprintf('No cluster was simulated.\n');
fprintf('No saved MAP point estimate was changed.\n');
fprintf('Detection profile-posterior uncertainty was calculated: YES.\n');
fprintf('Displayed surveillance settings = 16.\n');
fprintf('Joint p_p-p_a regions containing truth = %d of 16.\n', ...
    quality.jointContainsTruth);
fprintf('\nALL FIGURE 6 MATERIALS WERE SAVED IN:\n%s\n',outputFolder);
fprintf('\nPDF FIGURE:\n%s\n',fullfile(outputFolder,'Figure 6.pdf'));
fprintf('\nPNG FIGURE:\n%s\n',fullfile(outputFolder,'Figure 6.png'));
fprintf('\nMATLAB FIG FILE:\n%s\n',fullfile(outputFolder,'Figure 6.fig'));
fprintf('\nRESULTS CSV:\n%s\n',fullfile(outputFolder, ...
    'Figure6_MAP_Detection_Profiles_ManuscriptFinal_results.csv'));
fprintf('\nPROFILE-CONTOUR CSV:\n%s\n',fullfile(outputFolder, ...
    'Figure6_MAP_Detection_Profiles_ManuscriptFinal_profile_contours.csv'));
fprintf('\nREUSABLE PROFILE CACHE:\n%s\n',cacheFile);
fprintf('\nCOMPLETE RESULTS MAT FILE:\n%s\n',allResultsFile);
fprintf('\nSAVED CODE COPY:\n%s\n',fullfile(outputFolder, ...
    [programName '.m']));
fprintf('\nREADME:\n%s\n',fullfile(outputFolder, ...
    'README_Figure6_MAP_Detection_Profiles_ManuscriptFinal.txt'));
fprintf('\nCOMMAND-WINDOW LOG:\n%s\n',logFile);
fprintf('============================================================\n');

end

%% ========================================================================
% CONFIGURATION AND SOURCE LOCATION
% =========================================================================
function cfg = add_figure6_configuration(cfg)

requiredFields = {'lowerBounds','upperBounds','trueClusterTruncation', ...
    'nTrueClusters'};

for i = 1:numel(requiredFields)
    if ~isfield(cfg,requiredFields{i})
        error('The saved configuration is missing %s.',requiredFields{i});
    end
end

cfg.figure6ProfileBracketPoints = 13;
cfg.figure6ProfileBisectionIterations = 20;
cfg.figure6ContourAngles = 48;
cfg.figure6ContourBisectionIterations = 14;
cfg.figure6MaximumContourThresholdError = 0.02;
cfg.figure6AuditTolerance = 1e-4;
cfg.figure6BoundaryTolerance = 2e-5;

if ~isfield(cfg,'useParallel')
    cfg.useParallel = true;
end

end

function verify_required_toolboxes()

if isempty(ver('optim'))
    error(['Optimization Toolbox is required because the program uses ' ...
        'fmincon for direct profiling.']);
end

if isempty(ver('stats'))
    error(['Statistics and Machine Learning Toolbox is required because ' ...
        'the program uses chi-square quantiles.']);
end

end

function [sourceMat,observedCountsCSV,baselineCSV] = ...
    locate_map_sources(projectRoot,outputFolder)

matNames = { ...
    'MAP_Figures_5_7_Analysis_Results.mat', ...
    'Figure5_7_MAP_methodologicalComplete_All_Results.mat'};

sourceMat = locate_best_named_file( ...
    projectRoot,matNames,outputFolder,'saved MAP analysis MAT file');

releaseRoot = '';
matParent = fileparts(sourceMat);
[parentOfMatParent,matParentName] = fileparts(matParent);

if strcmpi(matParentName,'Saved_Analysis_Objects')
    releaseRoot = parentOfMatParent;
end

observedNames = { ...
    'Figure5_7_MAP_Observed_Cluster_Counts.csv', ...
    'Figure5_7_observed_cluster_counts.csv'};

if ~isempty(releaseRoot)
    candidate = fullfile(releaseRoot,'Data_and_Numerical_Results', ...
        observedNames{1});
else
    candidate = '';
end

if ~isempty(candidate) && exist(candidate,'file')
    observedCountsCSV = candidate;
else
    observedCountsCSV = locate_best_named_file( ...
        projectRoot,observedNames,outputFolder, ...
        'Figure 5/7 observed-cluster count CSV');
end

baselineCSV = '';
baselineNames = { ...
    'Figure5_6_MAP_Baseline_Estimates.csv', ...
    'Figure5_6_MAP_baseline_results.csv'};

if ~isempty(releaseRoot)
    candidate = fullfile(releaseRoot,'Data_and_Numerical_Results', ...
        baselineNames{1});
else
    candidate = '';
end

if ~isempty(candidate) && exist(candidate,'file')
    baselineCSV = candidate;
else
    baselineCSV = locate_optional_named_file( ...
        projectRoot,baselineNames,outputFolder);
end

end

function selected = locate_best_named_file( ...
    rootFolder,names,excludedFolder,description)

selected = locate_optional_named_file(rootFolder,names,excludedFolder);

if isempty(selected)
    nameList = strjoin(names,sprintf('\n  '));
    error(['The required %s was not found beneath:\n%s\n' ...
        'Expected one of:\n  %s'],description,rootFolder,nameList);
end

end

function selected = locate_optional_named_file( ...
    rootFolder,names,excludedFolder)

allCandidates = [];

for nameIndex = 1:numel(names)
    matches = dir(fullfile(rootFolder,'**',names{nameIndex}));
    matches = matches(~[matches.isdir]);

    for i = 1:numel(matches)
        fullName = fullfile(matches(i).folder,matches(i).name);

        if starts_with_path(fullName,excludedFolder)
            continue;
        end

        candidate = matches(i);
        candidate.fullName = fullName;
        candidate.namePriority = numel(names)-nameIndex+1;

        if isempty(allCandidates)
            allCandidates = candidate;
        else
            allCandidates(end+1) = candidate; %#ok<AGROW>
        end
    end
end

if isempty(allCandidates)
    selected = '';
    return;
end

scores = zeros(numel(allCandidates),1);

for i = 1:numel(allCandidates)
    pathText = lower(allCandidates(i).fullName);
    scores(i) = 1e6*allCandidates(i).namePriority;

    if contains(pathText,'final_release_')
        scores(i) = scores(i)+1e5;
    elseif contains(pathText,'release_')
        scores(i) = scores(i)+5e4;
    end

    if contains(pathText,'saved_analysis_objects') || ...
            contains(pathText,'data_and_numerical_results')
        scores(i) = scores(i)+1e4;
    end

    if contains(pathText,'before_') || contains(pathText,'backup')
        scores(i) = scores(i)-1e5;
    end

    scores(i) = scores(i)+allCandidates(i).datenum;
end

[~,bestIndex] = max(scores);
selected = allCandidates(bestIndex).fullName;

end

function answer = starts_with_path(candidate,parentFolder)

if isempty(parentFolder)
    answer = false;
    return;
end

candidate = char(candidate);
parentFolder = char(parentFolder);

answer = strncmpi(candidate,parentFolder,numel(parentFolder));

end

%% ========================================================================
% EXACT SAVED-DATA RECONSTRUCTION AND CHECKS
% =========================================================================
function countsByPanel = read_exact_observed_counts(csvFile,baseline)

T = readtable(csvFile,'TextType','string');

requiredNames = { ...
    'Panel','True_pp','True_pa','ObservedClusterSize','Count'};

for i = 1:numel(requiredNames)
    if ~ismember(requiredNames{i},T.Properties.VariableNames)
        error('The observed-count CSV is missing column %s.', ...
            requiredNames{i});
    end
end

countsByPanel = cell(16,1);
panelColumn = string(T.Panel);

for panelIndex = 1:16
    panelText = string(baseline(panelIndex).panel);
    rows = strcmpi(panelColumn,panelText);

    if ~any(rows)
        error('The observed-count CSV has no rows for panel %s.',panelText);
    end

    sizes = double(T.ObservedClusterSize(rows));
    frequencies = double(T.Count(rows));

    if any(~isfinite(sizes)) || any(~isfinite(frequencies)) || ...
            any(sizes<1) || any(mod(sizes,1)~=0) || ...
            any(frequencies<0) || any(mod(frequencies,1)~=0)
        error('Panel %s contains an invalid size or count.',panelText);
    end

    if numel(unique(sizes))~=numel(sizes)
        error('Panel %s contains a repeated cluster-size row.',panelText);
    end

    counts = zeros(max(sizes),1);
    counts(sizes) = frequencies;

    truePPValues = unique(double(T.True_pp(rows)));
    truePAValues = unique(double(T.True_pa(rows)));

    if numel(truePPValues)~=1 || numel(truePAValues)~=1 || ...
            abs(truePPValues-baseline(panelIndex).truePP)>1e-12 || ...
            abs(truePAValues-baseline(panelIndex).truePA)>1e-12
        error('The generating detection values disagree in panel %s.', ...
            panelText);
    end

    countsByPanel{panelIndex} = counts;
end

end

function verify_saved_map_inputs( ...
    baseline,countsByPanel,cfg,baselineCSV)

if numel(baseline)~=16
    error('The saved baseline analysis must contain exactly 16 panels.');
end

for i = 1:16
    result = baseline(i);

    requiredValues = [ ...
        result.truePP;result.truePA;result.nObserved; ...
        result.mapTheta(:);result.minimumObjective; ...
        result.rCI(:);result.kCI(:)];

    if any(~isfinite(requiredValues))
        error('Saved panel %d contains a missing or non-finite value.',i);
    end

    if numel(result.mapTheta)~=4
        error('Saved panel %d does not contain four MAP parameters.',i);
    end

    if sum(countsByPanel{i})~=result.nObserved
        error(['Panel (%s) contains %d clusters in the CSV but %d in ' ...
            'the saved MAP result.'],result.panel,sum(countsByPanel{i}), ...
            result.nObserved);
    end

    if any(result.mapTheta<cfg.lowerBounds) || ...
            any(result.mapTheta>cfg.upperBounds)
        error('The saved MAP estimate for panel (%s) is outside its bounds.', ...
            result.panel);
    end

    requiredPriorFields = {'ppAlpha','ppBeta','paAlpha','paBeta', ...
        'ppMean','paMean'};

    for f = 1:numel(requiredPriorFields)
        if ~isfield(result.prior,requiredPriorFields{f})
            error('Panel (%s) prior is missing %s.', ...
                result.panel,requiredPriorFields{f});
        end
    end
end

if ~isempty(baselineCSV)
    verify_baseline_csv(baseline,baselineCSV);
end

end

function verify_baseline_csv(baseline,csvFile)

T = readtable(csvFile,'TextType','string');

required = {'Panel','MAP_R_eff','MAP_k','MAP_pp','MAP_pa', ...
    'MinimumNegativeLogPosterior'};

for i = 1:numel(required)
    if ~ismember(required{i},T.Properties.VariableNames)
        error('The MAP summary CSV is missing column %s.',required{i});
    end
end

if height(T)~=16
    error('The MAP summary CSV must contain 16 panel rows.');
end

for i = 1:16
    row = find(strcmpi(string(T.Panel),string(baseline(i).panel)),1);

    if isempty(row)
        error('The MAP summary CSV is missing panel %s.',baseline(i).panel);
    end

    csvValues = [double(T.MAP_R_eff(row)) double(T.MAP_k(row)) ...
        double(T.MAP_pp(row)) double(T.MAP_pa(row))];

    if max(abs(csvValues-baseline(i).mapTheta))>1e-10
        error('The MAP summary CSV disagrees with the MAT file in panel %s.', ...
            baseline(i).panel);
    end

    if abs(double(T.MinimumNegativeLogPosterior(row))- ...
            baseline(i).minimumObjective)>1e-8
        error('The saved minimum disagrees in panel %s.',baseline(i).panel);
    end
end

fprintf('CSV-versus-MAT MAP estimate check passed exactly.\n');

end

function signature = make_source_signature( ...
    sourceMat,countsCSV,baseline,cfg)

matInfo = dir(sourceMat);
csvInfo = dir(countsCSV);
mapMatrix = vertcat(baseline.mapTheta);
minimumVector = [baseline.minimumObjective]';
weightedMapChecksum = sum(mapMatrix(:).*(1:numel(mapMatrix))');
weightedMinimumChecksum = sum(minimumVector.*(1:numel(minimumVector))');

signature = sprintf([ ...
    'Figure6-MAP-detection-profile-v1|' ...
    'mat=%d|csv=%d|map=%.15g|min=%.15g|J=%d|angles=%d|bisect=%d'], ...
    matInfo.bytes,csvInfo.bytes,weightedMapChecksum, ...
    weightedMinimumChecksum,cfg.trueClusterTruncation, ...
    cfg.figure6ContourAngles,cfg.figure6ContourBisectionIterations);

end

%% ========================================================================
% LIKELIHOOD AND POSTERIOR
% =========================================================================
function kernel = make_likelihood_kernel(counts,cfg)

counts = counts(:);
positiveX = find(counts>0);
Jmax = cfg.trueClusterTruncation;
j = 1:Jmax;
x = positiveX(:);
logCombination = -inf(numel(x),Jmax);

for row = 1:numel(x)
    valid = j>=x(row);
    logCombination(row,valid) = ...
        gammaln(j(valid)+1)-gammaln(x(row)+1) ...
        -gammaln(j(valid)-x(row)+1);
end

kernel = struct();
kernel.counts = counts;
kernel.positiveX = positiveX;
kernel.positiveCounts = counts(positiveX);
kernel.j = j;
kernel.x = x;
kernel.logCombination = logCombination;

end

function value = negative_log_posterior( ...
    theta,counts,kernel,prior,cfg)

if any(~isfinite(theta)) || any(theta<cfg.lowerBounds) || ...
        any(theta>cfg.upperBounds)
    value = 1e100;
    return;
end

probability = observed_size_probabilities(theta,kernel);

if any(~isfinite(probability)) || any(probability<=0)
    value = 1e100;
    return;
end

nll = -sum(kernel.positiveCounts.*log(probability));
pp = theta(3);
pa = theta(4);

logPriorPP = (prior.ppAlpha-1)*log(pp) ...
    +(prior.ppBeta-1)*log(1-pp) ...
    -betaln(prior.ppAlpha,prior.ppBeta);

logPriorPA = (prior.paAlpha-1)*log(pa) ...
    +(prior.paBeta-1)*log(1-pa) ...
    -betaln(prior.paAlpha,prior.paBeta);

value = nll-logPriorPP-logPriorPA;

if ~isfinite(value) || ~isreal(value)
    value = 1e100;
end

end

function probability = observed_size_probabilities(theta,kernel)

R = theta(1);
k = theta(2);
pp = theta(3);
pa = theta(4);
j = kernel.j;
x = kernel.x;

pTrue = true_cluster_pmf_untruncated(R,k,j);
z = pp+(1-pp)*pa;

logFirst = kernel.logCombination+x*log(z)+(j-x)*log(1-z);

logSecond = kernel.logCombination+x*log(pa)+(j-x)*log(1-pa) ...
    +ones(numel(x),1)*(j*log(1-pp));

conditionalMatrix = exp(logFirst)-exp(logSecond);

conditionalMatrix(conditionalMatrix<0 & ...
    conditionalMatrix>-1e-13) = 0;

if any(conditionalMatrix(:)<0)
    probability = nan(numel(x),1);
    return;
end

unconditionalPositive = conditionalMatrix*pTrue(:);
zeroProbability = sum(((1-pp).^j).*pTrue);
positiveProbability = 1-zeroProbability;

if positiveProbability<=0
    probability = nan(numel(x),1);
else
    probability = unconditionalPositive/positiveProbability;
end

end

function pTrue = true_cluster_pmf_untruncated(R,k,j)

logP = -log(j)+gammaln(k*j+j-1)-gammaln(k*j)-gammaln(j) ...
    +(j-1)*log(R/k)-(k*j+j-1)*log(1+R/k);

pTrue = exp(logP);

end

%% ========================================================================
% ONE-DIMENSIONAL DETECTION PROFILE INTERVALS
% =========================================================================
function [interval,lowerAtBound,upperAtBound,thresholdError] = ...
    direct_profile_interval( ...
    parameterIndex,mapTheta,minObjective,counts,kernel,prior,cfg)

threshold = minObjective+0.5*chi2inv(0.95,1);
mapValue = mapTheta(parameterIndex);
lowerBound = cfg.lowerBounds(parameterIndex);
upperBound = cfg.upperBounds(parameterIndex);

lowerCandidates = linspace(mapValue,lowerBound, ...
    cfg.figure6ProfileBracketPoints);
upperCandidates = linspace(mapValue,upperBound, ...
    cfg.figure6ProfileBracketPoints);

[lowerValue,lowerAtBound,lowerError] = find_profile_crossing( ...
    parameterIndex,mapTheta,minObjective,threshold,counts,kernel, ...
    prior,cfg,lowerCandidates);

[upperValue,upperAtBound,upperError] = find_profile_crossing( ...
    parameterIndex,mapTheta,minObjective,threshold,counts,kernel, ...
    prior,cfg,upperCandidates);

interval = [lowerValue upperValue];
finiteErrors = [lowerError upperError];
finiteErrors = finiteErrors(isfinite(finiteErrors));

if isempty(finiteErrors)
    thresholdError = 0;
else
    thresholdError = max(finiteErrors);
end

end

function [crossing,atBound,errorAtCrossing] = find_profile_crossing( ...
    parameterIndex,mapTheta,minObjective,threshold,counts,kernel, ...
    prior,cfg,candidates)

previousValue = candidates(1);
previousDifference = minObjective-threshold;
crossing = candidates(end);
atBound = true;
errorAtCrossing = nan;

for i = 2:numel(candidates)
    currentValue = candidates(i);
    currentObjective = profile_with_one_fixed( ...
        parameterIndex,currentValue,mapTheta,minObjective, ...
        counts,kernel,prior,cfg);
    currentDifference = currentObjective-threshold;

    if currentDifference>=0
        left = previousValue;
        right = currentValue;
        leftDifference = previousDifference;
        rightDifference = currentDifference;

        for iteration = 1:cfg.figure6ProfileBisectionIterations
            middle = 0.5*(left+right);
            middleObjective = profile_with_one_fixed( ...
                parameterIndex,middle,mapTheta,minObjective, ...
                counts,kernel,prior,cfg);
            middleDifference = middleObjective-threshold;

            if sign(middleDifference)==sign(leftDifference)
                left = middle;
                leftDifference = middleDifference;
            else
                right = middle;
                rightDifference = middleDifference;
            end
        end

        if abs(leftDifference)<abs(rightDifference)
            crossing = left;
            errorAtCrossing = abs(leftDifference);
        else
            crossing = right;
            errorAtCrossing = abs(rightDifference);
        end

        atBound = false;
        return;
    end

    previousValue = currentValue;
    previousDifference = currentDifference;
end

end

function value = profile_with_one_fixed( ...
    fixedIndex,fixedValue,mapTheta,minObjective,counts,kernel,prior,cfg)

remaining = setdiff(1:4,fixedIndex);
start = mapTheta(remaining);
lower = cfg.lowerBounds(remaining);
upper = cfg.upperBounds(remaining);

options = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
    'MaxFunctionEvaluations',2500,'MaxIterations',700, ...
    'OptimalityTolerance',1e-7,'StepTolerance',1e-7);

objective = @(free) fixed_parameter_objective( ...
    free,remaining,fixedIndex,fixedValue,mapTheta, ...
    counts,kernel,prior,cfg);

[~,value] = fmincon(objective,start,[],[],[],[], ...
    lower,upper,[],options);

value = audit_profile_value(value,minObjective,cfg, ...
    sprintf('parameter %d fixed at %.8g',fixedIndex,fixedValue));

end

function value = fixed_parameter_objective( ...
    free,remaining,fixedIndex,fixedValue,mapTheta, ...
    counts,kernel,prior,cfg)

theta = mapTheta;
theta(fixedIndex) = fixedValue;
theta(remaining) = free;
value = negative_log_posterior(theta,counts,kernel,prior,cfg);

end

%% ========================================================================
% JOINT p_p-p_a PROFILE REGION
% =========================================================================
function contour = direct_joint_detection_contour( ...
    mapTheta,minObjective,ppCI,paCI,counts,kernel,prior,cfg)

target = minObjective+0.5*chi2inv(0.95,2);
scalePP = max(0.5*(ppCI(2)-ppCI(1)),0.04);
scalePA = max(0.5*(paCI(2)-paCI(1)),0.04);

angles = linspace(0,2*pi,cfg.figure6ContourAngles+1);
angles(end) = [];

contourPP = nan(size(angles));
contourPA = nan(size(angles));
errors = nan(size(angles));
boundaryShortfalls = zeros(size(angles));
touches = false(size(angles));

canParallel = cfg.useParallel && ~isempty(ver('parallel')) && ...
    ~isempty(gcp('nocreate'));

if canParallel
    parfor angleIndex = 1:numel(angles)
        [contourPP(angleIndex),contourPA(angleIndex), ...
            errors(angleIndex),touches(angleIndex), ...
            boundaryShortfalls(angleIndex)] = contour_point_for_angle( ...
            angles(angleIndex),mapTheta,minObjective,target, ...
            scalePP,scalePA,counts,kernel,prior,cfg);
    end
else
    for angleIndex = 1:numel(angles)
        [contourPP(angleIndex),contourPA(angleIndex), ...
            errors(angleIndex),touches(angleIndex), ...
            boundaryShortfalls(angleIndex)] = contour_point_for_angle( ...
            angles(angleIndex),mapTheta,minObjective,target, ...
            scalePP,scalePA,counts,kernel,prior,cfg);
    end
end

finiteCrossingErrors = errors(~touches & isfinite(errors));

if isempty(finiteCrossingErrors)
    maxThresholdError = 0;
else
    maxThresholdError = max(finiteCrossingErrors);
end

if maxThresholdError>cfg.figure6MaximumContourThresholdError
    error(['A detection contour threshold error (%.6g) exceeds the ' ...
        'prespecified tolerance (%.6g).'],maxThresholdError, ...
        cfg.figure6MaximumContourThresholdError);
end

contour = struct();
contour.pp = [contourPP contourPP(1)];
contour.pa = [contourPA contourPA(1)];
contour.maxThresholdError = maxThresholdError;
contour.touchesBoundary = any(touches);
contour.maxBoundaryShortfall = max(boundaryShortfalls);

end

function [pointPP,pointPA,thresholdError,touchesBoundary, ...
    boundaryShortfall] = contour_point_for_angle( ...
    angle,mapTheta,minObjective,target,scalePP,scalePA, ...
    counts,kernel,prior,cfg)

directionPP = cos(angle);
directionPA = sin(angle);
maximumRadius = maximum_detection_radius( ...
    mapTheta,scalePP,scalePA,directionPP,directionPA,cfg);

lowRadius = 0;
lowDifference = minObjective-target;
highRadius = min(0.5,maximumRadius);

highDifference = profiled_detection_difference( ...
    highRadius,mapTheta,minObjective,target,scalePP,scalePA, ...
    directionPP,directionPA,counts,kernel,prior,cfg);

while highDifference<0 && highRadius<maximumRadius
    lowRadius = highRadius;
    lowDifference = highDifference;
    highRadius = min(maximumRadius,2*highRadius);
    highDifference = profiled_detection_difference( ...
        highRadius,mapTheta,minObjective,target,scalePP,scalePA, ...
        directionPP,directionPA,counts,kernel,prior,cfg);
end

if highDifference<0
    radius = maximumRadius;
    touchesBoundary = true;
    thresholdError = nan;
    boundaryShortfall = abs(highDifference);
else
    touchesBoundary = false;
    boundaryShortfall = 0;

    for iteration = 1:cfg.figure6ContourBisectionIterations
        middleRadius = 0.5*(lowRadius+highRadius);
        middleDifference = profiled_detection_difference( ...
            middleRadius,mapTheta,minObjective,target,scalePP,scalePA, ...
            directionPP,directionPA,counts,kernel,prior,cfg);

        if middleDifference>=0
            highRadius = middleRadius;
            highDifference = middleDifference;
        else
            lowRadius = middleRadius;
            lowDifference = middleDifference;
        end
    end

    if abs(lowDifference)<abs(highDifference)
        radius = lowRadius;
        thresholdError = abs(lowDifference);
    else
        radius = highRadius;
        thresholdError = abs(highDifference);
    end
end

pointPP = mapTheta(3)+radius*scalePP*directionPP;
pointPA = mapTheta(4)+radius*scalePA*directionPA;

pointPP = min(max(pointPP,cfg.lowerBounds(3)),cfg.upperBounds(3));
pointPA = min(max(pointPA,cfg.lowerBounds(4)),cfg.upperBounds(4));

end

function maximumRadius = maximum_detection_radius( ...
    mapTheta,scalePP,scalePA,directionPP,directionPA,cfg)

limits = inf(1,2);

if directionPP>0
    limits(1) = (cfg.upperBounds(3)-mapTheta(3)) ...
        /(scalePP*directionPP);
elseif directionPP<0
    limits(1) = (cfg.lowerBounds(3)-mapTheta(3)) ...
        /(scalePP*directionPP);
end

if directionPA>0
    limits(2) = (cfg.upperBounds(4)-mapTheta(4)) ...
        /(scalePA*directionPA);
elseif directionPA<0
    limits(2) = (cfg.lowerBounds(4)-mapTheta(4)) ...
        /(scalePA*directionPA);
end

maximumRadius = min(limits);

if ~isfinite(maximumRadius) || maximumRadius<=0
    error('A valid maximum contour radius could not be calculated.');
end

end

function difference = profiled_detection_difference( ...
    radius,mapTheta,minObjective,target,scalePP,scalePA, ...
    directionPP,directionPA,counts,kernel,prior,cfg)

fixedPP = mapTheta(3)+radius*scalePP*directionPP;
fixedPA = mapTheta(4)+radius*scalePA*directionPA;

fixedPP = min(max(fixedPP,cfg.lowerBounds(3)),cfg.upperBounds(3));
fixedPA = min(max(fixedPA,cfg.lowerBounds(4)),cfg.upperBounds(4));

profileValue = profile_at_detection_pair( ...
    fixedPP,fixedPA,mapTheta,minObjective,counts,kernel,prior,cfg);

difference = profileValue-target;

end

function profileValue = profile_at_detection_pair( ...
    fixedPP,fixedPA,mapTheta,minObjective,counts,kernel,prior,cfg)

start = mapTheta(1:2);
lower = cfg.lowerBounds(1:2);
upper = cfg.upperBounds(1:2);

options = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
    'MaxFunctionEvaluations',1800,'MaxIterations',500, ...
    'OptimalityTolerance',1e-7,'StepTolerance',1e-7);

objective = @(transmission) detection_fixed_objective( ...
    transmission,fixedPP,fixedPA,mapTheta,counts,kernel,prior,cfg);

[~,profileValue] = fmincon(objective,start,[],[],[],[], ...
    lower,upper,[],options);

profileValue = audit_profile_value(profileValue,minObjective,cfg, ...
    sprintf('p_p=%.8g and p_a=%.8g',fixedPP,fixedPA));

end

function value = detection_fixed_objective( ...
    transmission,fixedPP,fixedPA,mapTheta,counts,kernel,prior,cfg)

theta = mapTheta;
theta(1:2) = transmission;
theta(3) = fixedPP;
theta(4) = fixedPA;
value = negative_log_posterior(theta,counts,kernel,prior,cfg);

end

function value = audit_profile_value(value,minObjective,cfg,description)

if value<minObjective-cfg.figure6AuditTolerance
    error(['A profile calculation found a negative log posterior below ' ...
        'the saved MAP minimum.\nProfile condition: %s.\n' ...
        'Saved minimum = %.12f.\nProfiled value = %.12f.\n' ...
        'The saved MAP analysis must be audited before Figure 6 is used.'], ...
        description,minObjective,value);
end

if value<minObjective
    value = minObjective;
end

end

function tolerance = objective_match_tolerance(objective)

tolerance = max(1e-6,1e-10*abs(objective));

end

%% ========================================================================
% CHECKPOINTING AND QUALITY SUMMARIES
% =========================================================================
function result = empty_profile_result()

result = struct( ...
    'complete',false, ...
    'panel','', ...
    'truePP',nan, ...
    'truePA',nan, ...
    'evaluatedTruePP',nan, ...
    'evaluatedTruePA',nan, ...
    'nObserved',nan, ...
    'mapTheta',nan(1,4), ...
    'minimumObjective',nan, ...
    'ppCI',nan(1,2), ...
    'paCI',nan(1,2), ...
    'ppLowerAtBound',false, ...
    'ppUpperAtBound',false, ...
    'paLowerAtBound',false, ...
    'paUpperAtBound',false, ...
    'profileThresholdError',nan, ...
    'contourPP',nan, ...
    'contourPA',nan, ...
    'contourThresholdError',nan, ...
    'contourTouchesBoundary',false, ...
    'contourBoundaryShortfall',nan, ...
    'truePairProfileDelta',nan, ...
    'truePairInsideJoint95',false, ...
    'truePPInsideSeparate95',false, ...
    'truePAInsideSeparate95',false, ...
    'elapsedSeconds',nan);

end

function save_profile_cache_atomic( ...
    cacheFile,profileResults,sourceSignature,sourceMat, ...
    observedCountsCSV,cfg)

temporaryFile = [cacheFile '.temporary.mat'];

if exist(temporaryFile,'file')
    delete(temporaryFile);
end

save(temporaryFile,'profileResults','sourceSignature','sourceMat', ...
    'observedCountsCSV','cfg','-v7.3');

test = load(temporaryFile,'profileResults','sourceSignature');

if numel(test.profileResults)~=16 || ...
        ~strcmp(test.sourceSignature,sourceSignature)
    error('The temporary Figure 6 cache did not pass its read-back check.');
end

[moveSucceeded,moveMessage] = movefile(temporaryFile,cacheFile,'f');

if ~moveSucceeded
    error('The Figure 6 cache could not be finalized: %s',moveMessage);
end

end

function quality = summarize_profile_quality(profileResults,baseline,cfg)

mapMatrix = vertcat(baseline.mapTheta);
truePP = [baseline.truePP]';
truePA = [baseline.truePA]';
identifiablePA = truePP<1-cfg.figure6BoundaryTolerance;

quality = struct();
quality.ppMAE = mean(abs(mapMatrix(:,3)-truePP));
quality.paMAEIdentifiable = ...
    mean(abs(mapMatrix(identifiablePA,4)-truePA(identifiablePA)));
quality.ppSeparateContainsTruth = ...
    sum([profileResults.truePPInsideSeparate95]);
quality.paSeparateContainsTruthIdentifiable = ...
    sum([profileResults(identifiablePA).truePAInsideSeparate95]);
quality.jointContainsTruth = ...
    sum([profileResults.truePairInsideJoint95]);
quality.jointTouchesBoundary = ...
    sum([profileResults.contourTouchesBoundary]);
quality.maxContourThresholdError = ...
    max([profileResults.contourThresholdError]);
quality.maxProfileThresholdError = ...
    max([profileResults.profileThresholdError]);

end

function answer = yes_no(value)

if value
    answer = 'YES';
else
    answer = 'NO';
end

end

%% ========================================================================
% FIGURE
% =========================================================================
function fig = draw_figure6(results)

fig = figure('Color','w','Units','pixels', ...
    'Position',[30 30 2400 1500],'Renderer','painters', ...
    'Name','Figure 6 - MAP detection estimates');

fontName = 'Helvetica';
left = 0.065;
right = 0.985;
bottom = 0.105;
top = 0.815;
horizontalGap = 0.025;
verticalGap = 0.040;
panelWidth = (right-left-3*horizontalGap)/4;
panelHeight = (top-bottom-3*verticalGap)/4;

for i = 1:16
    row = ceil(i/4);
    col = mod(i-1,4)+1;
    xPosition = left+(col-1)*(panelWidth+horizontalGap);
    yPosition = top-row*panelHeight-(row-1)*verticalGap;

    ax = axes('Parent',fig,'Position', ...
        [xPosition yPosition panelWidth panelHeight]);
    hold(ax,'on');

    plot(ax,results(i).contourPP,results(i).contourPA, ...
        'k-','LineWidth',2.0);

    plot(ax,results(i).ppCI, ...
        [results(i).mapTheta(4) results(i).mapTheta(4)], ...
        'b--','LineWidth',1.6);

    plot(ax,[results(i).mapTheta(3) results(i).mapTheta(3)], ...
        results(i).paCI,'b--','LineWidth',1.6);

    plot(ax,results(i).mapTheta(3),results(i).mapTheta(4),'o', ...
        'MarkerSize',7,'MarkerFaceColor',[0.05 0.35 0.95], ...
        'MarkerEdgeColor',[0 0.15 0.65]);

    plot(ax,results(i).truePP,results(i).truePA,'d', ...
        'MarkerSize',8,'MarkerFaceColor',[1 0.12 0.12], ...
        'MarkerEdgeColor',[0.75 0 0]);

    xlim(ax,[-0.03 1.03]);
    ylim(ax,[-0.03 1.03]);
    set(ax,'XTick',[0 0.25 0.50 0.75 1.00], ...
        'YTick',[0 0.25 0.50 0.75 1.00]);
    grid(ax,'on');
    box(ax,'on');

    set(ax,'FontName',fontName,'FontSize',9.5,'LineWidth',0.9, ...
        'Layer','top','TickDir','out');

    title(ax,sprintf('(%s)  p_p=%.2f, p_a=%.2f', ...
        results(i).panel,results(i).truePP,results(i).truePA), ...
        'Interpreter','tex','FontName',fontName, ...
        'FontWeight','bold','FontSize',10.5);

    if col>1
        set(ax,'YTickLabel',[]);
    end

    if row<4
        set(ax,'XTickLabel',[]);
    end
end

add_canvas_text(fig,[0.05 0.945 0.90 0.04], ...
    'MAP detection estimates',20,'bold',0,fontName,'none');

draw_figure_key(fig,fontName);

add_canvas_text(fig,[0.33 0.018 0.38 0.035], ...
    'Passive detection probability, p_p',15,'bold',0,fontName,'tex');

add_canvas_text(fig,[0.004 0.31 0.035 0.40], ...
    'Active detection probability, p_a',15,'bold',90,fontName,'tex');

add_canvas_text(fig,[0.13 0.055 0.74 0.028], ...
    ['For generating p_p=1, p_a is not identified by the cluster data; ' ...
    'its displayed MAP result reflects prior regularization.'], ...
    10,'normal',0,fontName,'tex');

end

function draw_figure_key(fig,fontName)

ax = axes('Parent',fig,'Position',[0.12 0.855 0.76 0.055], ...
    'Visible','off','XLim',[0 1],'YLim',[0 1]);
hold(ax,'on');

h1 = plot(ax,nan,nan,'k-','LineWidth',2.0);
h2 = plot(ax,nan,nan,'o','MarkerSize',7, ...
    'MarkerFaceColor',[0.05 0.35 0.95], ...
    'MarkerEdgeColor',[0 0.15 0.65]);
h3 = plot(ax,nan,nan,'d','MarkerSize',8, ...
    'MarkerFaceColor',[1 0.12 0.12], ...
    'MarkerEdgeColor',[0.75 0 0]);
h4 = plot(ax,nan,nan,'b--','LineWidth',1.6);

legend(ax,[h1 h2 h3 h4], ...
    {'95% joint region','MAP','Truth','95% 1-D intervals'}, ...
    'Orientation','horizontal','Location','north', ...
    'Box','off','FontName',fontName,'FontSize',11);

end

function add_canvas_text( ...
    fig,position,textValue,fontSize,fontWeight,rotation,fontName,interpreter)

% Use an invisible axes plus a text object instead of an annotation
% textbox.  Some MATLAB releases do not support the Rotation property for
% annotation textboxes, whereas ordinary text objects support Rotation.
% The axes position is expressed in normalized figure coordinates, so the
% title, shared axis labels, and footnote retain their intended locations.
textAxes = axes('Parent',fig, ...
    'Units','normalized', ...
    'Position',position, ...
    'XLim',[0 1], ...
    'YLim',[0 1], ...
    'Visible','off');

text(textAxes,0.5,0.5,textValue, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontName',fontName, ...
    'FontSize',fontSize, ...
    'FontWeight',fontWeight, ...
    'Rotation',rotation, ...
    'Interpreter',interpreter, ...
    'Clipping','off');

end

function save_publication_figure(fig,outputFolder,baseName)

pngFile = fullfile(outputFolder,[baseName '.png']);
pdfFile = fullfile(outputFolder,[baseName '.pdf']);
figFile = fullfile(outputFolder,[baseName '.fig']);

savefig(fig,figFile);
set(fig,'PaperPositionMode','auto');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 16 10], ...
    'PaperSize',[16 10]);

print(fig,pngFile,'-dpng','-r300');
set(fig,'Renderer','painters');
print(fig,pdfFile,'-dpdf','-painters');

if ~exist(pngFile,'file') || dir(pngFile).bytes<10000
    error('The PNG export failed or is unexpectedly small: %s',pngFile);
end

if ~exist(pdfFile,'file') || dir(pdfFile).bytes<5000
    error('The PDF export failed or is unexpectedly small: %s',pdfFile);
end

pngInfo = imfinfo(pngFile);
fprintf('    PNG canvas = %d by %d pixels.\n', ...
    pngInfo.Width,pngInfo.Height);

end

%% ========================================================================
% CSV, CODE, AND README OUTPUTS
% =========================================================================
function write_results_csv(results,outputFolder)

n = numel(results);
Panel = strings(n,1);
True_pp = nan(n,1);
True_pa = nan(n,1);
ObservedClusters = nan(n,1);
MAP_R_eff = nan(n,1);
MAP_k = nan(n,1);
MAP_pp = nan(n,1);
pp_Lower95 = nan(n,1);
pp_Upper95 = nan(n,1);
MAP_pa = nan(n,1);
pa_Lower95 = nan(n,1);
pa_Upper95 = nan(n,1);
pp_LowerAtBound = false(n,1);
pp_UpperAtBound = false(n,1);
pa_LowerAtBound = false(n,1);
pa_UpperAtBound = false(n,1);
True_pp_In_Separate95 = false(n,1);
True_pa_In_Separate95 = false(n,1);
TruePair_ProfileDeltaNLP = nan(n,1);
TruePair_In_Joint95 = false(n,1);
JointRegionTouchesBoundary = false(n,1);
MaximumContourThresholdError = nan(n,1);
MaximumProfileThresholdError = nan(n,1);
ElapsedSeconds = nan(n,1);

for i = 1:n
    Panel(i) = string(results(i).panel);
    True_pp(i) = results(i).truePP;
    True_pa(i) = results(i).truePA;
    ObservedClusters(i) = results(i).nObserved;
    MAP_R_eff(i) = results(i).mapTheta(1);
    MAP_k(i) = results(i).mapTheta(2);
    MAP_pp(i) = results(i).mapTheta(3);
    pp_Lower95(i) = results(i).ppCI(1);
    pp_Upper95(i) = results(i).ppCI(2);
    MAP_pa(i) = results(i).mapTheta(4);
    pa_Lower95(i) = results(i).paCI(1);
    pa_Upper95(i) = results(i).paCI(2);
    pp_LowerAtBound(i) = results(i).ppLowerAtBound;
    pp_UpperAtBound(i) = results(i).ppUpperAtBound;
    pa_LowerAtBound(i) = results(i).paLowerAtBound;
    pa_UpperAtBound(i) = results(i).paUpperAtBound;
    True_pp_In_Separate95(i) = results(i).truePPInsideSeparate95;
    True_pa_In_Separate95(i) = results(i).truePAInsideSeparate95;
    TruePair_ProfileDeltaNLP(i) = results(i).truePairProfileDelta;
    TruePair_In_Joint95(i) = results(i).truePairInsideJoint95;
    JointRegionTouchesBoundary(i) = results(i).contourTouchesBoundary;
    MaximumContourThresholdError(i) = results(i).contourThresholdError;
    MaximumProfileThresholdError(i) = results(i).profileThresholdError;
    ElapsedSeconds(i) = results(i).elapsedSeconds;
end

T = table(Panel,True_pp,True_pa,ObservedClusters,MAP_R_eff,MAP_k, ...
    MAP_pp,pp_Lower95,pp_Upper95,MAP_pa,pa_Lower95,pa_Upper95, ...
    pp_LowerAtBound,pp_UpperAtBound,pa_LowerAtBound,pa_UpperAtBound, ...
    True_pp_In_Separate95,True_pa_In_Separate95, ...
    TruePair_ProfileDeltaNLP,TruePair_In_Joint95, ...
    JointRegionTouchesBoundary,MaximumContourThresholdError, ...
    MaximumProfileThresholdError,ElapsedSeconds);

writetable(T,fullfile(outputFolder, ...
    'Figure6_MAP_Detection_Profiles_ManuscriptFinal_results.csv'));

end

function write_contour_csv(results,outputFolder)

Panel = strings(0,1);
Vertex = zeros(0,1);
Profile_pp = zeros(0,1);
Profile_pa = zeros(0,1);

for i = 1:numel(results)
    nVertices = numel(results(i).contourPP);
    Panel = [Panel; repmat(string(results(i).panel),nVertices,1)]; %#ok<AGROW>
    Vertex = [Vertex; (1:nVertices)']; %#ok<AGROW>
    Profile_pp = [Profile_pp; results(i).contourPP(:)]; %#ok<AGROW>
    Profile_pa = [Profile_pa; results(i).contourPA(:)]; %#ok<AGROW>
end

T = table(Panel,Vertex,Profile_pp,Profile_pa);

writetable(T,fullfile(outputFolder, ...
    'Figure6_MAP_Detection_Profiles_ManuscriptFinal_profile_contours.csv'));

end

function copy_running_code(outputFolder,programName)

sourceCode = [mfilename('fullpath') '.m'];
destinationCode = fullfile(outputFolder,[programName '.m']);

if ~exist(sourceCode,'file')
    warning('The running MATLAB code could not be copied: %s',sourceCode);
    return;
end

if ~strcmpi(sourceCode,destinationCode)
    copyfile(sourceCode,destinationCode,'f');
end

end

function write_readme( ...
    outputFolder,sourceMat,observedCountsCSV,quality,cfg)

readmeFile = fullfile(outputFolder, ...
    'README_Figure6_MAP_Detection_Profiles_ManuscriptFinal.txt');

fid = fopen(readmeFile,'w');

if fid<0
    error('The Figure 6 README could not be created.');
end

cleanupFile = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid,'MANUSCRIPT-FINAL MAP FIGURE 6\n');
fprintf(fid,'================================\n\n');
fprintf(fid,'Purpose\n');
fprintf(fid,['Display the MAP detection estimates in the same inferential ' ...
    'format as Figure 5, with p_p on the horizontal axis and p_a on the ' ...
    'vertical axis.\n\n']);
fprintf(fid,'Scientific provenance\n');
fprintf(fid,'Saved MAP object: %s\n',sourceMat);
fprintf(fid,'Observed-count CSV: %s\n\n',observedCountsCSV);
fprintf(fid,'No cluster was simulated by this program.\n');
fprintf(fid,'No saved MAP point estimate was changed.\n');
fprintf(fid,['Separate 95%% intervals profile over the other three ' ...
    'parameters.\n']);
fprintf(fid,['Joint 95%% p_p-p_a regions profile over R_eff and k.\n']);
fprintf(fid,'One-dimensional threshold: 2 Delta NLP <= %.8f.\n', ...
    chi2inv(0.95,1));
fprintf(fid,'Two-dimensional threshold: 2 Delta NLP <= %.8f.\n\n', ...
    chi2inv(0.95,2));
fprintf(fid,'Quality summary\n');
fprintf(fid,'Mean absolute error for p_p: %.8f\n',quality.ppMAE);
fprintf(fid,['Mean absolute error for data-identifiable p_a settings: ' ...
    '%.8f\n'],quality.paMAEIdentifiable);
fprintf(fid,'Separate p_p intervals containing truth: %d/16\n', ...
    quality.ppSeparateContainsTruth);
fprintf(fid,['Separate p_a intervals containing truth among generating ' ...
    'p_p<1 settings: %d/12\n'], ...
    quality.paSeparateContainsTruthIdentifiable);
fprintf(fid,'Joint regions containing truth: %d/16\n', ...
    quality.jointContainsTruth);
fprintf(fid,'Joint regions touching a bound: %d/16\n', ...
    quality.jointTouchesBoundary);
fprintf(fid,'Largest contour threshold error: %.10g\n', ...
    quality.maxContourThresholdError);
fprintf(fid,'Required contour threshold-error tolerance: %.10g\n\n', ...
    cfg.figure6MaximumContourThresholdError);
fprintf(fid,'Interpretation of generating p_p=1\n');
fprintf(fid,['When every case is detected passively, no passively missed ' ...
    'case remains for active surveillance to find. The cluster data then ' ...
    'do not identify p_a. Any displayed MAP location or posterior-profile ' ...
    'restriction for p_a in that column arises from prior regularization ' ...
    'and its interaction with p_p, not from direct p_a information in the ' ...
    'cluster data.\n']);

end

%% ========================================================================
% PARALLEL POOL
% =========================================================================
function maybe_start_parallel_pool(cfg)

if ~cfg.useParallel || isempty(ver('parallel'))
    fprintf(['Parallel processing is unavailable or disabled. ' ...
        'The profiles will be calculated serially.\n\n']);
    return;
end

try
    pool = gcp('nocreate');

    if isempty(pool)
        parpool('local');
    end

    fprintf('Parallel processing is available for the joint contours.\n\n');
catch err
    fprintf('Parallel pool could not be started: %s\n',err.message);
    fprintf('The workflow will continue serially.\n\n');
end

end

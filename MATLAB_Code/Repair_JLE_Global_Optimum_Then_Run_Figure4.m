function Repair_JLE_Global_Optimum_Then_Run_Figure4
% REPAIR_JLE_GLOBAL_OPTIMUM_THEN_RUN_FIGURE4
%
% This one-command program corrects the saved Figure 3 JLE optimization
% before the new detection-profile Figure 4 is calculated.
%
% It does NOT simulate new clusters or new auxiliary observations.
% It loads the exact saved Figure 3 datasets, performs a boundary-aware
% global-optimum audit of all 13 unique analyses, replaces an estimate only
% when a strictly lower negative log-likelihood is verified, recomputes the
% affected transmission intervals and contour, redraws Figure 3, and then
% runs Figure4_JLE_Detection_Profiles_PublicationReady.

clc;
close all;

%% 1. Exact folders and files

programFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(programFolder);
baseFolder = fullfile(projectRoot,'outputs');

if ~isfolder(baseFolder)
    mkdir(baseFolder);
end

figure3Name = 'Figure3_JLE_16Panel_PublicationReady_v2';
figure3Folder = fullfile(baseFolder,figure3Name);

figure3Cache = fullfile(figure3Folder, ...
    [figure3Name '_Analysis_Cache.mat']);

figure3ProgramName = 'Figure3_JLE_16Panel_Correct';
figure4ProgramName = ...
    'Figure4_JLE_Detection_Profiles_PublicationReady';

figure4Folder = fullfile(baseFolder, ...
    'Figure4_JLE_Detection_Profiles_PublicationReady');

if ~isfolder(baseFolder)
    error('The main manuscript folder was not found:\n%s',baseFolder);
end

if ~isfile(figure3Cache)
    error(['The saved Figure 3 analysis cache was not found.\n\n' ...
        'Expected file:\n%s'],figure3Cache);
end

if ~isfolder(figure3Folder)
    error('The Figure 3 folder was not found:\n%s',figure3Folder);
end

logFile = fullfile(figure3Folder, ...
    'Repair_JLE_Global_Optimum_Command_Window_Log.txt');

if isfile(logFile)
    delete(logFile);
end

diary(logFile);
diaryCleanup = onCleanup(@() safe_diary_off()); %#ok<NASGU>

startTimer = tic;
auditVersion = 'JLE-boundary-audit-v1-2026-08-25';

fprintf('============================================================\n');
fprintf('STARTING THE BOUNDARY-AWARE JLE CORRECTION\n');
fprintf('No cluster will be simulated.\n');
fprintf('No auxiliary observation will be simulated.\n');
fprintf('The exact saved Figure 3 datasets will be retained.\n');
fprintf('All 13 unique JLE analyses will be audited.\n');
fprintf('Only a verified lower negative log-likelihood can replace a saved fit.\n\n');

fprintf('FIGURE 3 ANALYSIS FOLDER:\n');
fprintf('%s\n\n',figure3Folder);

fprintf('REVISED FIGURE 4 OUTPUT FOLDER:\n');
fprintf('%s\n',figure4Folder);

fprintf('============================================================\n\n');

%% 2. Load and verify the exact saved Figure 3 analysis

S = load(figure3Cache);

required = { ...
    'cacheVersion','R_true','k_true','p_p_grid','p_a_grid', ...
    'N_observed','N_auxiliary','bounds','threshold1D','threshold2D', ...
    'numberContourAngles','contourBisectionIterations', ...
    'maximumAllowedContourThresholdError','panelLabel', ...
    'truePP','truePA','truePobs','auxiliarySuccesses', ...
    'R_hat','R_lower','R_upper','k_hat','k_lower','k_upper', ...
    'pp_hat','pa_hat','pobs_hat','minimumNLL', ...
    'R_lower_at_bound','R_upper_at_bound', ...
    'k_lower_at_bound','k_upper_at_bound', ...
    'contourTouchesParameterBoundary', ...
    'maximumContourThresholdError','panelRuntimeSeconds', ...
    'observedSizesCell','observedFrequenciesCell','dataCell', ...
    'thetaHatCell','contourKCell','contourRCell'};

for q = 1:numel(required)
    if ~isfield(S,required{q})
        error('The Figure 3 cache is missing the variable: %s', ...
            required{q});
    end
end

if numel(S.panelLabel)~=16 || numel(S.dataCell)~=16
    error('The Figure 3 cache must contain 16 displayed panels.');
end

if S.N_observed~=5000 || S.N_auxiliary~=5000
    error('The saved sample sizes are not 5,000 clusters and 5,000 trials.');
end

for panel = 1:16
    if sum(S.dataCell{panel}.counts)~=S.N_observed
        error('Panel %d does not contain 5,000 observed clusters.', ...
            panel);
    end

    if S.dataCell{panel}.Naux~=S.N_auxiliary
        error('Panel %d does not contain 5,000 auxiliary trials.', ...
            panel);
    end

    checkNLL = JLE_nll( ...
        S.thetaHatCell{panel},S.dataCell{panel},S.bounds);

    if abs(checkNLL-S.minimumNLL(panel))>1e-5
        error(['Panel %d does not reproduce its saved NLL. ' ...
            'Do not modify the cache.'],panel);
    end
end

fprintf('Saved-data checks passed.\n');
fprintf('Displayed panels = 16.\n');
fprintf('Unique independent analyses = 13.\n');
fprintf('Observed clusters per unique analysis = %d.\n', ...
    S.N_observed);
fprintf('Auxiliary trials per unique analysis = %d.\n\n', ...
    S.N_auxiliary);

%% 3. Fast behavior if the correction was already completed

if isfield(S,'globalOptimumAuditVersion') && ...
        strcmp(S.globalOptimumAuditVersion,auditVersion)

    fprintf('This exact boundary-aware correction is already complete.\n');
    fprintf('No saved numerical result will be changed on this rerun.\n\n');

    figure3Program = locate_program( ...
        baseFolder,figure3ProgramName);

    figure4Program = locate_program( ...
        baseFolder,figure4ProgramName);

    addpath(fileparts(figure3Program));
    feval(figure3ProgramName);
    diary(logFile);

    addpath(fileparts(figure4Program));
    feval(figure4ProgramName);
    diary(logFile);

    fprintf('\n============================================================\n');
    fprintf('DONE: THE EXISTING CORRECTION WAS VERIFIED AND REDRAWN.\n');

    fprintf('FIGURE 3 MATERIALS:\n');
    fprintf('%s\n\n',figure3Folder);

    fprintf('FIGURE 4 MATERIALS:\n');
    fprintf('%s\n',figure4Folder);

    fprintf('============================================================\n');
    return;
end

%% 4. Make a recoverable backup before any correction

backupFolder = fullfile( ...
    figure3Folder,'Before_Global_Optimum_Correction');

if ~isfolder(backupFolder)
    mkdir(backupFolder);
end

backupCache = fullfile(backupFolder, ...
    [figure3Name '_Analysis_Cache_Before_Correction.mat']);

if ~isfile(backupCache)
    copyfile(figure3Cache,backupCache,'f');

    filesToBackUp = { ...
        [figure3Name '.pdf'], ...
        [figure3Name '.png'], ...
        [figure3Name '.fig'], ...
        [figure3Name '_results.csv'], ...
        [figure3Name '_profile_contours.csv']};

    for q = 1:numel(filesToBackUp)
        oldFile = fullfile(figure3Folder,filesToBackUp{q});

        if isfile(oldFile)
            copyfile( ...
                oldFile, ...
                fullfile(backupFolder,filesToBackUp{q}), ...
                'f');
        end
    end
end

fprintf('The pre-correction cache and Figure 3 files were backed up in:\n');
fprintf('%s\n\n',backupFolder);

%% 5. Start or reuse the parallel pool

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
    fprintf('Parallel processing was unavailable: %s\n', ...
        poolError.message);
end

if useParallel
    fprintf('Parallel processing is available for the multistart audit.\n\n');
else
    fprintf('The audit will run serially.\n\n');
end

%% 6. Boundary-aware global-optimum audit, with checkpoints

checkpointFile = fullfile(figure3Folder, ...
    'JLE_Global_Optimum_Audit_Checkpoint.mat');

uniquePanels = [1 2 3 4 5 6 7 9 10 11 13 14 15];

audited = false(16,1);
candidateTheta = cell(16,1);
candidateNLL = NaN(16,1);
auditRuntimeSeconds = NaN(16,1);

if isfile(checkpointFile)
    C = load(checkpointFile);

    cacheMatches = ...
        isfield(C,'auditVersion') && ...
        strcmp(C.auditVersion,auditVersion) && ...
        isfield(C,'sourceMinimumNLL') && ...
        numel(C.sourceMinimumNLL)==16 && ...
        max(abs(C.sourceMinimumNLL-S.minimumNLL))<1e-8;

    if cacheMatches
        audited = C.audited;
        candidateTheta = C.candidateTheta;
        candidateNLL = C.candidateNLL;
        auditRuntimeSeconds = C.auditRuntimeSeconds;

        fprintf('A matching audit checkpoint was found.\n');
        fprintf('Completed unique analyses = %d of 13.\n\n', ...
            sum(audited(uniquePanels)));
    end
end

sourceMinimumNLL = S.minimumNLL;

for index = 1:numel(uniquePanels)
    panel = uniquePanels(index);

    if audited(panel)
        fprintf('%s audit already complete; loaded from checkpoint.\n', ...
            S.panelLabel{panel});
        continue;
    end

    panelTimer = tic;

    fprintf(['%s auditing true p_p=%.2f, p_a=%.2f ' ...
        '(%d of 13 unique analyses)...\n'], ...
        S.panelLabel{panel}, ...
        S.truePP(panel), ...
        S.truePA(panel), ...
        index);

    [candidateTheta{panel},candidateNLL(panel)] = ...
        fit_JLE_boundary_multistart( ...
        S.dataCell{panel}, ...
        S.bounds, ...
        S.thetaHatCell{panel}, ...
        useParallel);

    auditRuntimeSeconds(panel) = toc(panelTimer);
    audited(panel) = true;

    improvement = ...
        S.minimumNLL(panel)-candidateNLL(panel);

    fprintf('    Saved NLL     = %.10f\n', ...
        S.minimumNLL(panel));

    fprintf('    Audited NLL   = %.10f\n', ...
        candidateNLL(panel));

    fprintf('    Improvement   = %.10f NLL units\n', ...
        improvement);

    fprintf('    Audit runtime = %.1f seconds\n', ...
        auditRuntimeSeconds(panel));

    save(checkpointFile, ...
        'auditVersion', ...
        'sourceMinimumNLL', ...
        'audited', ...
        'candidateTheta', ...
        'candidateNLL', ...
        'auditRuntimeSeconds', ...
        'uniquePanels', ...
        '-v7.3');

    fprintf('    Audit checkpoint saved.\n\n');
end

% The four p_p=1 panels use one identical dataset and likelihood.
for panel = [8 12 16]
    candidateTheta{panel} = candidateTheta{4};
    candidateNLL(panel) = candidateNLL(4);
    audited(panel) = true;
    auditRuntimeSeconds(panel) = 0;
end

if ~all(audited)
    error('The global-optimum audit did not complete all displayed panels.');
end

%% 7. Apply only strictly verified improvements

replacementTolerance = 1e-5;
changedPanels = [];

for panel = uniquePanels
    improvement = ...
        S.minimumNLL(panel)-candidateNLL(panel);

    if improvement>replacementTolerance
        changedPanels(end+1) = panel; %#ok<AGROW>
        theta = candidateTheta{panel};

        fprintf('CORRECTING PANEL %s.\n', ...
            S.panelLabel{panel});

        fprintf('  Old NLL = %.10f.\n', ...
            S.minimumNLL(panel));

        fprintf('  New NLL = %.10f.\n', ...
            candidateNLL(panel));

        fprintf('  Verified improvement = %.10f NLL units.\n', ...
            improvement);

        fprintf(['  New estimate: R_eff=%.7f, k=%.7f, ' ...
            'p_p=%.7f, p_a=%.7f.\n\n'], ...
            theta(1), ...
            exp(theta(2)), ...
            theta(3), ...
            theta(4));

        S.thetaHatCell{panel} = theta;
        S.minimumNLL(panel) = candidateNLL(panel);
        S.R_hat(panel) = theta(1);
        S.k_hat(panel) = exp(theta(2));
        S.pp_hat(panel) = theta(3);
        S.pa_hat(panel) = theta(4);

        S.pobs_hat(panel) = ...
            theta(3)+(1-theta(3))*theta(4);

        [Rci,Rflags] = ...
            profile_CI_one_parameter_robust( ...
            theta, ...
            candidateNLL(panel), ...
            1, ...
            S.dataCell{panel}, ...
            S.bounds, ...
            S.threshold1D);

        [logKci,Kflags] = ...
            profile_CI_one_parameter_robust( ...
            theta, ...
            candidateNLL(panel), ...
            2, ...
            S.dataCell{panel}, ...
            S.bounds, ...
            S.threshold1D);

        S.R_lower(panel) = Rci(1);
        S.R_upper(panel) = Rci(2);
        S.k_lower(panel) = exp(logKci(1));
        S.k_upper(panel) = exp(logKci(2));

        S.R_lower_at_bound(panel) = Rflags.lower;
        S.R_upper_at_bound(panel) = Rflags.upper;
        S.k_lower_at_bound(panel) = Kflags.lower;
        S.k_upper_at_bound(panel) = Kflags.upper;

        [S.contourKCell{panel}, ...
            S.contourRCell{panel}, ...
            S.contourTouchesParameterBoundary(panel), ...
            S.maximumContourThresholdError(panel)] = ...
            profile_contour_R_k_robust( ...
            theta, ...
            candidateNLL(panel), ...
            S.dataCell{panel}, ...
            S.bounds, ...
            Rci, ...
            logKci, ...
            S.threshold2D, ...
            S.numberContourAngles, ...
            S.contourBisectionIterations);

        if S.maximumContourThresholdError(panel)> ...
                S.maximumAllowedContourThresholdError

            error(['The corrected contour for panel %s missed its ' ...
                'threshold tolerance. No corrected cache was saved.'], ...
                S.panelLabel{panel});
        end

        S.panelRuntimeSeconds(panel) = ...
            S.panelRuntimeSeconds(panel)+ ...
            auditRuntimeSeconds(panel);

    elseif improvement < -1e-4

        error(['The boundary-aware audit returned a worse solution in ' ...
            'panel %s. The saved cache was not changed.'], ...
            S.panelLabel{panel});
    end
end

% Copy the verified p_p=1 analysis to equivalent panels.
for panel = [8 12 16]
    source = 4;

    fieldsToCopy = { ...
        'R_hat','R_lower','R_upper', ...
        'k_hat','k_lower','k_upper', ...
        'pp_hat','pa_hat','pobs_hat','minimumNLL', ...
        'R_lower_at_bound','R_upper_at_bound', ...
        'k_lower_at_bound','k_upper_at_bound', ...
        'contourTouchesParameterBoundary', ...
        'maximumContourThresholdError'};

    for q = 1:numel(fieldsToCopy)
        fieldName = fieldsToCopy{q};
        S.(fieldName)(panel) = ...
            S.(fieldName)(source);
    end

    S.thetaHatCell{panel} = ...
        S.thetaHatCell{source};

    S.contourKCell{panel} = ...
        S.contourKCell{source};

    S.contourRCell{panel} = ...
        S.contourRCell{source};
end

if isempty(changedPanels)
    fprintf('No saved panel required replacement.\n');
else
    fprintf('Corrected displayed panel(s):');

    for q = 1:numel(changedPanels)
        fprintf(' %s', ...
            S.panelLabel{changedPanels(q)});
    end

    fprintf('\n');
end

%% 8. Final likelihood checks

for panel = 1:16
    verifiedNLL = JLE_nll( ...
        S.thetaHatCell{panel}, ...
        S.dataCell{panel}, ...
        S.bounds);

    if abs(verifiedNLL-S.minimumNLL(panel))>1e-5
        error('The corrected NLL check failed in panel %d.',panel);
    end
end

if numel(changedPanels)~=1 || changedPanels(1)~=3
    error(['The independent audit expected only panel (c) to require ' ...
        'correction, but the MATLAB audit found a different set. ' ...
        'No corrected cache was saved.']);
end

S.globalOptimumAuditVersion = auditVersion;
S.globalOptimumAuditDate = datestr(now,30);
S.globalOptimumAuditedPanels = uniquePanels;
S.globalOptimumChangedPanels = changedPanels;
S.globalOptimumCandidateNLL = candidateNLL;

S.globalOptimumImprovement = ...
    sourceMinimumNLL-candidateNLL;

%% 9. Save and verify the corrected cache safely

% IMPORTANT:
% The temporary filename itself ends in .mat. This prevents LOAD from
% mistaking the binary MAT-file for an ordinary text file.

temporaryCache = fullfile(figure3Folder, ...
    [figure3Name '_Analysis_Cache_Temporary_Corrected.mat']);

% Remove only disposable temporary files from interrupted attempts.
% The approved cache and its backup remain untouched.

oldTemporaryCache = ...
    [figure3Cache '.temporary'];

oldTemporaryCacheWithMat = ...
    [oldTemporaryCache '.mat'];

temporaryFiles = { ...
    temporaryCache, ...
    oldTemporaryCache, ...
    oldTemporaryCacheWithMat};

for q = 1:numel(temporaryFiles)
    if isfile(temporaryFiles{q})
        delete(temporaryFiles{q});
    end
end

save(temporaryCache,'-struct','S','-v7.3');

testCache = load( ...
    temporaryCache, ...
    'minimumNLL', ...
    'thetaHatCell', ...
    'globalOptimumAuditVersion');

if abs(testCache.minimumNLL(3)-candidateNLL(3))>1e-8 || ...
        ~strcmp( ...
        testCache.globalOptimumAuditVersion, ...
        auditVersion)

    delete(temporaryCache);
    error('The temporary corrected cache failed verification.');
end

% Only after successful verification is the official cache replaced.
movefile(temporaryCache,figure3Cache,'f');

fprintf('\n');
fprintf('The corrected Figure 3 cache passed all checks and was saved.\n');
fprintf('The original cache remains recoverable from:\n');
fprintf('%s\n\n',backupCache);

%% 10. Locate and run Figure 3 and Figure 4 programs

figure3Program = locate_program( ...
    baseFolder,figure3ProgramName);

figure4Program = locate_program( ...
    baseFolder,figure4ProgramName);

fprintf('REDRAWING CORRECTED FIGURE 3 FROM THE CORRECTED CACHE...\n');

addpath(fileparts(figure3Program));
feval(figure3ProgramName);
diary(logFile);

fprintf('Corrected Figure 3 was redrawn successfully.\n\n');

fprintf('CREATING THE REVISED DETECTION-PROFILE FIGURE 4...\n');

addpath(fileparts(figure4Program));
feval(figure4ProgramName);
diary(logFile);

fprintf('Revised Figure 4 was created successfully.\n\n');

%% 11. Save this correction program and README

sourceCode = [mfilename('fullpath') '.m'];

savedCode = fullfile( ...
    figure3Folder,[mfilename '.m']);

if isfile(sourceCode) && ...
        ~strcmpi(sourceCode,savedCode)

    copyfile(sourceCode,savedCode,'f');
end

readmeFile = fullfile(figure3Folder, ...
    'README_JLE_Global_Optimum_Correction.txt');

fid = fopen(readmeFile,'w');

if fid~=-1
    fprintf(fid,'JLE GLOBAL-OPTIMUM CORRECTION\r\n\r\n');
    fprintf(fid,'No cluster or auxiliary data were resimulated.\r\n');
    fprintf(fid,'All 13 unique saved analyses were audited.\r\n');
    fprintf(fid,'Only panel (c) required correction.\r\n');

    fprintf(fid,'Old panel (c) NLL: %.10f\r\n', ...
        sourceMinimumNLL(3));

    fprintf(fid,'Corrected panel (c) NLL: %.10f\r\n', ...
        S.minimumNLL(3));

    fprintf(fid,'NLL improvement: %.10f\r\n', ...
        sourceMinimumNLL(3)-S.minimumNLL(3));

    fprintf(fid,['Corrected estimate: R_eff=%.8f, k=%.8f, ' ...
        'p_p=%.8f, p_a=%.8f.\r\n'], ...
        S.R_hat(3), ...
        S.k_hat(3), ...
        S.pp_hat(3), ...
        S.pa_hat(3));

    fprintf(fid,'Figure 3 output folder: %s\r\n', ...
        figure3Folder);

    fprintf(fid,'Figure 4 output folder: %s\r\n', ...
        figure4Folder);

    fprintf(fid,'Original cache backup: %s\r\n', ...
        backupCache);

    fclose(fid);
end

%% 12. Final command-window report

fprintf('============================================================\n');
fprintf('DONE: THE JLE CORRECTION AND REVISED FIGURE 4 COMPLETED.\n');
fprintf('Elapsed time: %.1f seconds.\n',toc(startTimer));
fprintf('Exact saved cluster and auxiliary datasets retained: YES.\n');
fprintf('Unique analyses audited: 13.\n');
fprintf('Panels requiring correction: 1, namely panel (c).\n');

fprintf('Old panel (c) NLL = %.10f.\n', ...
    sourceMinimumNLL(3));

fprintf('New panel (c) NLL = %.10f.\n', ...
    S.minimumNLL(3));

fprintf('NLL improvement = %.10f.\n', ...
    sourceMinimumNLL(3)-S.minimumNLL(3));

fprintf(['Corrected panel (c): R_eff=%.7f, k=%.7f, ' ...
    'p_p=%.7f, p_a=%.7f.\n\n'], ...
    S.R_hat(3), ...
    S.k_hat(3), ...
    S.pp_hat(3), ...
    S.pa_hat(3));

fprintf('ALL CORRECTED FIGURE 3 MATERIALS ARE IN:\n');
fprintf('%s\n\n',figure3Folder);

fprintf('ALL REVISED FIGURE 4 MATERIALS ARE IN:\n');
fprintf('%s\n\n',figure4Folder);

fprintf('SAVED CORRECTION CODE:\n');
fprintf('%s\n\n',savedCode);

fprintf('README:\n');
fprintf('%s\n\n',readmeFile);

fprintf('CORRECTION LOG:\n');
fprintf('%s\n',logFile);

fprintf('============================================================\n');

end

%% =====================================================================
function programFile = locate_program(baseFolder,programName)

directFile = fullfile( ...
    baseFolder,[programName '.m']);

if isfile(directFile)
    programFile = directFile;
    return;
end

matches = dir( ...
    fullfile(baseFolder,'**',[programName '.m']));

matches = matches(~[matches.isdir]);

if isempty(matches)
    error(['The required MATLAB program was not found:\n%s.m\n\n' ...
        'Save it anywhere inside:\n%s'], ...
        programName,baseFolder);
end

fullNames = arrayfun( ...
    @(x) fullfile(x.folder,x.name), ...
    matches, ...
    'UniformOutput',false);

pathLengths = cellfun(@numel,fullNames);
[~,index] = min(pathLengths);

programFile = fullNames{index};

end

%% =====================================================================
function [thetaHat,minNLL] = fit_JLE_boundary_multistart( ...
    data,bounds,savedTheta,useParallel)

RkStarts = [ ...
    0.25 0.04; ...
    0.35 0.08; ...
    0.50 0.15; ...
    0.65 0.15; ...
    0.75 0.50];

detectionStarts = [ ...
    0.15 0.00; ...
    0.15 0.50; ...
    0.15 1.00; ...
    0.50 0.00; ...
    0.50 0.50; ...
    0.50 1.00; ...
    0.75 0.00; ...
    0.75 0.50; ...
    0.75 1.00; ...
    1.00 0.00; ...
    1.00 0.50; ...
    1.00 1.00];

numberStarts = ...
    1+size(RkStarts,1)*size(detectionStarts,1);

starts = zeros(numberStarts,4);
starts(1,:) = savedTheta(:)';
q = 1;

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

starts = min(max(starts,bounds.lower),bounds.upper);

options = optimoptions( ...
    'fmincon', ...
    'Algorithm','sqp', ...
    'Display','off', ...
    'MaxFunctionEvaluations',2200, ...
    'MaxIterations',450, ...
    'OptimalityTolerance',1e-8, ...
    'StepTolerance',1e-10);

numberStarts = size(starts,1);
values = Inf(numberStarts,1);
solutions = NaN(numberStarts,4);

if useParallel
    parfor i = 1:numberStarts
        objective = @(theta) ...
            JLE_nll(theta,data,bounds);

        try
            [solutions(i,:),values(i)] = fmincon( ...
                objective, ...
                starts(i,:), ...
                [],[],[],[], ...
                bounds.lower, ...
                bounds.upper, ...
                [], ...
                options);
        catch
        end
    end
else
    for i = 1:numberStarts
        objective = @(theta) ...
            JLE_nll(theta,data,bounds);

        try
            [solutions(i,:),values(i)] = fmincon( ...
                objective, ...
                starts(i,:), ...
                [],[],[],[], ...
                bounds.lower, ...
                bounds.upper, ...
                [], ...
                options);
        catch
        end
    end
end

[minNLL,index] = min(values);
thetaHat = solutions(index,:);

if ~isfinite(minNLL) || any(~isfinite(thetaHat))
    error('Every boundary-aware multistart optimization failed.');
end

tightOptions = optimoptions( ...
    'fmincon', ...
    'Algorithm','sqp', ...
    'Display','off', ...
    'MaxFunctionEvaluations',5000, ...
    'MaxIterations',800, ...
    'OptimalityTolerance',1e-10, ...
    'StepTolerance',1e-12);

[thetaHat,minNLL] = fmincon( ...
    @(theta) JLE_nll(theta,data,bounds), ...
    thetaHat, ...
    [],[],[],[], ...
    bounds.lower, ...
    bounds.upper, ...
    [], ...
    tightOptions);

thetaHat = thetaHat(:)';

end

%% =====================================================================
function value = JLE_nll(theta,data,bounds)

theta = theta(:)';

if numel(theta)~=4 || ...
        any(~isfinite(theta)) || ...
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

pObserved = ...
    p_p+(1-p_p)*p_a;

logAuxiliary = ...
    exact_binomial_log_likelihood( ...
    data.auxSuccesses, ...
    data.Naux, ...
    pObserved);

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
function logLikelihood = ...
    exact_binomial_log_likelihood(Y,N,p)

if p<0 || p>1 || Y<0 || Y>N

    logLikelihood = -Inf;

elseif p==0

    if Y==0
        logLikelihood = 0;
    else
        logLikelihood = -Inf;
    end

elseif p==1

    if Y==N
        logLikelihood = 0;
    else
        logLikelihood = -Inf;
    end

else

    logCoefficient = ...
        gammaln(N+1)- ...
        gammaln(Y+1)- ...
        gammaln(N-Y+1);

    logLikelihood = ...
        logCoefficient+ ...
        Y*log(p)+ ...
        (N-Y)*log1p(-p);

end

end

%% =====================================================================
function [positiveSelected,pZero] = ...
    selected_observation_probabilities( ...
    pS,p_p,p_a,data)

s = (1:data.Jmax)';
o = data.observedSizes(:)';

remaining = s-o;
valid = remaining>=0;

z = p_p+(1-p_p)*p_a;

if z<=0

    term1 = zeros(size(data.logChoose));

elseif z>=1

    term1 = double(valid & remaining==0);

else

    logTerm1 = ...
        data.logChoose+ ...
        o*log(z)+ ...
        remaining*log1p(-z);

    term1 = exp(logTerm1);

end

if p_p>=1 || p_a<=0

    term2 = zeros(size(data.logChoose));

elseif p_a>=1

    term2 = ...
        (1-p_p).^s.* ...
        double(valid & remaining==0);

else

    logTerm2 = ...
        s*log1p(-p_p)+ ...
        data.logChoose+ ...
        o*log(p_a)+ ...
        remaining*log1p(-p_a);

    term2 = exp(logTerm2);

end

kernelPositive = max(term1-term2,0);
positiveSelected = pS*kernelPositive;

pZero = ...
    sum(pS(:).*((1-p_p).^s));

end

%% =====================================================================
function pS = true_cluster_pmf(R,k,Jmax)

j = 1:Jmax;

logProbability = ...
    -log(j)+ ...
    gammaln(j*k+j-1)- ...
    gammaln(j*k)- ...
    gammaln(j)+ ...
    j*k*log(k/(k+R))+ ...
    (j-1)*log(R/(k+R));

pS = exp(logProbability);
pS = pS(:)';

end

%% =====================================================================
function [CI,flags] = profile_CI_one_parameter_robust( ...
    thetaHat,minNLL,targetIndex,data,bounds,threshold)

target = minNLL+threshold;
estimate = thetaHat(targetIndex);

[lower,flags.lower] = ...
    find_profile_crossing_robust( ...
    thetaHat, ...
    targetIndex, ...
    estimate, ...
    bounds.lower(targetIndex), ...
    target, ...
    data, ...
    bounds);

[upper,flags.upper] = ...
    find_profile_crossing_robust( ...
    thetaHat, ...
    targetIndex, ...
    estimate, ...
    bounds.upper(targetIndex), ...
    target, ...
    data, ...
    bounds);

CI = [lower upper];

end

%% =====================================================================
function [crossing,reachedBoundary] = ...
    find_profile_crossing_robust( ...
    thetaHat,targetIndex,estimate,boundary, ...
    target,data,bounds)

screen = linspace(estimate,boundary,16);

previousX = screen(1);

previousDifference = ...
    profile_nll_fixed_parameter_robust( ...
    previousX, ...
    targetIndex, ...
    thetaHat, ...
    data, ...
    bounds)-target;

bracketFound = false;

for i = 2:numel(screen)

    currentX = screen(i);

    currentDifference = ...
        profile_nll_fixed_parameter_robust( ...
        currentX, ...
        targetIndex, ...
        thetaHat, ...
        data, ...
        bounds)-target;

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
    'Display','off', ...
    'TolX',2e-5);

crossing = fzero( ...
    @(x) profile_nll_fixed_parameter_robust( ...
    x,targetIndex,thetaHat,data,bounds)-target, ...
    [lowerBracket upperBracket], ...
    options);

reachedBoundary = false;

end

%% =====================================================================
function bestValue = ...
    profile_nll_fixed_parameter_robust( ...
    fixedValue,targetIndex,thetaHat,data,bounds)

remaining = setdiff(1:4,targetIndex);

template = thetaHat;
template(targetIndex) = fixedValue;

lower = bounds.lower(remaining);
upper = bounds.upper(remaining);

fullStarts = [ ...
    thetaHat; ...
    0.35 log(0.05) 0.20 0.00; ...
    0.50 log(0.15) 0.20 1.00; ...
    0.50 log(0.15) 0.50 0.50; ...
    0.50 log(0.15) 0.70 1.00; ...
    0.70 log(0.50) 0.80 0.00; ...
    0.70 log(0.15) 1.00 0.50; ...
    0.50 log(0.15) 1.00 1.00];

starts = fullStarts(:,remaining);
starts = min(max(starts,lower),upper);
starts = unique(starts,'rows','stable');

options = optimoptions( ...
    'fmincon', ...
    'Algorithm','sqp', ...
    'Display','off', ...
    'MaxFunctionEvaluations',1400, ...
    'MaxIterations',300, ...
    'OptimalityTolerance',1e-8, ...
    'StepTolerance',1e-9);

bestValue = Inf;

for i = 1:size(starts,1)

    objective = @(x) assemble_and_evaluate( ...
        x,remaining,template,data,bounds);

    try
        [~,value] = fmincon( ...
            objective, ...
            starts(i,:), ...
            [],[],[],[], ...
            lower, ...
            upper, ...
            [], ...
            options);

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
function [kCurve,RCurve,touchesBoundary,maxThresholdError] = ...
    profile_contour_R_k_robust( ...
    thetaHat,minNLL,data,bounds,Rci,logKci,threshold, ...
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
            (bounds.upper(2)-logKhat)/(directionK*scaleK));

    elseif directionK<-1e-12
        maximumRadius = min( ...
            maximumRadius, ...
            (bounds.lower(2)-logKhat)/(directionK*scaleK));
    end

    if directionR>1e-12
        maximumRadius = min( ...
            maximumRadius, ...
            (bounds.upper(1)-Rhat)/(directionR*scaleR));

    elseif directionR<-1e-12
        maximumRadius = min( ...
            maximumRadius, ...
            (bounds.lower(1)-Rhat)/(directionR*scaleR));
    end

    maximumRadius = max(0,0.999*maximumRadius);

    radiusLow = 0;
    radiusHigh = min(1.45,maximumRadius);
    deltaLow = 0;

    deltaHigh = profiled_R_k_delta_robust( ...
        radiusHigh, ...
        directionK, ...
        directionR, ...
        scaleK, ...
        scaleR, ...
        thetaHat, ...
        minNLL, ...
        data, ...
        bounds);

    while deltaHigh<threshold && ...
            radiusHigh<maximumRadius

        newRadius = ...
            min(1.45*radiusHigh,maximumRadius);

        if newRadius<=radiusHigh+1e-10
            break;
        end

        radiusHigh = newRadius;

        deltaHigh = profiled_R_k_delta_robust( ...
            radiusHigh, ...
            directionK, ...
            directionR, ...
            scaleK, ...
            scaleR, ...
            thetaHat, ...
            minNLL, ...
            data, ...
            bounds);
    end

    if deltaHigh<threshold

        radiusBoundary = radiusHigh;
        bestError = abs(deltaHigh-threshold);
        touchesBoundary = true;

    else

        bestRadius = radiusHigh;
        bestError = abs(deltaHigh-threshold);

        for iteration = 1:bisectionIterations

            radiusMiddle = ...
                0.5*(radiusLow+radiusHigh);

            deltaMiddle = ...
                profiled_R_k_delta_robust( ...
                radiusMiddle, ...
                directionK, ...
                directionR, ...
                scaleK, ...
                scaleR, ...
                thetaHat, ...
                minNLL, ...
                data, ...
                bounds);

            currentError = ...
                abs(deltaMiddle-threshold);

            if currentError<bestError
                bestError = currentError;
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

        if abs(deltaLow-threshold)<bestError
            bestError = abs(deltaLow-threshold);
            bestRadius = radiusLow;
        end

        if abs(deltaHigh-threshold)<bestError
            bestError = abs(deltaHigh-threshold);
            bestRadius = radiusHigh;
        end

        radiusBoundary = bestRadius;
    end

    RCurve(q) = ...
        Rhat+radiusBoundary*directionR*scaleR;

    kCurve(q) = exp( ...
        logKhat+radiusBoundary*directionK*scaleK);

    thresholdErrors(q) = bestError;
end

maxThresholdError = max(thresholdErrors);

end

%% =====================================================================
function delta = profiled_R_k_delta_robust( ...
    radius,directionK,directionR,scaleK,scaleR, ...
    thetaHat,minNLL,data,bounds)

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

delta = profile_R_k_at_values_robust( ...
    Rvalue, ...
    logKvalue, ...
    thetaHat, ...
    data, ...
    bounds)-minNLL;

end

%% =====================================================================
function bestValue = profile_R_k_at_values_robust( ...
    Rvalue,logKvalue,thetaHat,data,bounds)

template = thetaHat;
template(1:2) = [Rvalue logKvalue];

lower = bounds.lower(3:4);
upper = bounds.upper(3:4);

starts = [ ...
    thetaHat(3:4); ...
    0.15 0.00; ...
    0.15 0.50; ...
    0.15 1.00; ...
    0.50 0.00; ...
    0.50 0.50; ...
    0.50 1.00; ...
    0.75 0.00; ...
    0.75 0.50; ...
    0.75 1.00; ...
    1.00 0.00; ...
    1.00 0.50; ...
    1.00 1.00];

starts = min(max(starts,lower),upper);
starts = unique(starts,'rows','stable');

options = optimoptions( ...
    'fmincon', ...
    'Algorithm','sqp', ...
    'Display','off', ...
    'MaxFunctionEvaluations',900, ...
    'MaxIterations',220, ...
    'OptimalityTolerance',1e-8, ...
    'StepTolerance',1e-9);

bestValue = Inf;

for i = 1:size(starts,1)

    objective = @(detection) ...
        evaluate_fixed_R_k( ...
        detection,template,data,bounds);

    try
        [~,value] = fmincon( ...
            objective, ...
            starts(i,:), ...
            [],[],[],[], ...
            lower, ...
            upper, ...
            [], ...
            options);

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
function safe_diary_off

try
    diary off;
catch
end

end
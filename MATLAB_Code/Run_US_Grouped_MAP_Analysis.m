function Run_US_Grouped_MAP_Analysis
% Run_US_Grouped_MAP_Analysis
%
% publication-ready grouped-data analysis for U.S. measles transmission chains.
%
% IMPORTANT SCIENTIFIC POINTS
% 1. This program never invents the 92 unavailable individual chain sizes.
% 2. It uses only the published yearly grouped summaries:
%       number of chains, total cases, number of size-1 chains,
%       number of size-2 chains, and number of size>=3 chains.
% 3. The main analysis uses completed years 2020-2023.
% 4. The 2024 snapshot is a prespecified sensitivity analysis because one
%    2024 outbreak was still ongoing at the publication cutoff.
% 5. R_eff is restricted to (0,1). Therefore the analysis estimates how
%    close transmission is to one under a subcritical model; it does not
%    prove subcriticality from the data.
% 6. p_p and p_a use transparent priors. No historical study is silently
%    treated as measuring a quantity that it did not measure.
%
% OUTPUTS
% - checkpoint MAT files, so completed calculations are reused
% - Figure 8: primary and full-snapshot posterior summaries
% - Figure 9: prespecified prior-location sensitivity
% - Figure S1: simulation validation
% - Figure S2: posterior-predictive checks
% - CSV, MAT, code copy, README, and command-window log

close all;
clc;
workflowTimer=tic;

%% ------------------------------------------------------------------------
% 1. USER-EDITABLE CONFIGURATION
% -------------------------------------------------------------------------
cfg.version = 'GroupedMAP_20260819_v2_exact_series';

programFolder = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(programFolder);
cfg.baseFolder = fullfile(repositoryRoot,'outputs');

if ~isfolder(cfg.baseFolder)
    mkdir(cfg.baseFolder);
end
cfg.analysisFolderName = 'US_Grouped_MAP_methodologicalFinal';

% The first publication run uses 200 independent datasets per scenario.
% For a very short code test only, change this temporarily to 5. Restore
% it to 200 before using results in the manuscript.
cfg.validationReplicates = 200;
cfg.validationBatchSize = 20;

% Likelihood and probability tolerances. The likelihood itself uses an
% exact truncated power-series recursion through the largest cluster size
% required by the grouped data. SmaxDirectCheck is used only for an
% independent, slower direct-summation verification at final fitted values.
cfg.SmaxDirectCheck = 12000;
cfg.OmaxSamplerInitial = 500;
cfg.maxSamplerO = 4000;
% Across all validation replicates in one scenario, the expected number of
% sampled observed clusters lying beyond the retained PMF must be <=0.05.
% This is a sample-size-based accuracy rule, rather than an unnecessarily
% strict fixed tail probability.
cfg.maximumExpectedOmittedValidationDraws = 0.05;
cfg.probabilityTolerance = 1e-10;
cfg.finalLikelihoodTolerance = 1e-6;
cfg.maximumExpectedOmittedChains = 1e-3;

% Optimization.
cfg.numberOfStarts = 12;
cfg.validationNumberOfStarts = 4;
cfg.optimMaxIter = 1600;
cfg.optimMaxFunEvals = 6000;

% MCMC. Four chains are required for convergence assessment.
cfg.numberOfChains = 4;
cfg.mcmcIterations = 12000;
cfg.mcmcBurnIn = 4000;
cfg.mcmcThin = 2;
cfg.mcmcAdaptEvery = 100;
cfg.minimumESS = 400;
cfg.maximumRhat = 1.01;

% The nine prior-sensitivity posteriors first use exact importance
% reweighting of the primary posterior draws. A separate MCMC fit is run
% automatically only when overlap is inadequate.
cfg.minimumImportanceESS = 1000;
cfg.maximumImportanceWeight = 0.02;

% Posterior-predictive draws. These are checkpointed.
cfg.ppcDraws = 300;

% Reproducibility.
cfg.masterSeed = 260819;

% Use the Parallel Computing Toolbox when available.
cfg.useParallel = true;

% Scientifically prespecified parameter support.
cfg.logKLower = log(1e-3);
cfg.logKUpper = log(1e3);

% Primary weak detection priors.
cfg.primaryPrior.ppAlpha = 2;
cfg.primaryPrior.ppBeta  = 2;
cfg.primaryPrior.paAlpha = 2;
cfg.primaryPrior.paBeta  = 2;
cfg.primaryPrior.name = 'Beta(2,2) for p_p and p_a';

% prespecified sensitivity prior locations. Concentration = 20.
cfg.sensitivityMeans = [0.20 0.50 0.80];
cfg.sensitivityConcentration = 20;

% The program stops before calling results publication-ready if these fail.
cfg.stopOnFailedDiagnostics = true;

%% ------------------------------------------------------------------------
% 2. FOLDERS, DIARY, AND STARTING MESSAGE
% -------------------------------------------------------------------------
analysisRoot = fullfile(cfg.baseFolder,cfg.analysisFolderName);
checkpointFolder = fullfile(analysisRoot,'Checkpoints');
if ~exist(analysisRoot,'dir'); mkdir(analysisRoot); end
if ~exist(checkpointFolder,'dir'); mkdir(checkpointFolder); end

runStamp = datestr(now,'yyyymmdd_HHMMSS');
runFolder = fullfile(analysisRoot,['Run_' runStamp]);
mkdir(runFolder);

logFile = fullfile(runFolder, ...
    'US_Grouped_MAP_methodologicalFinal_Command_Window_Log.txt');
diary(logFile);
diary on;
cleanupDiary = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('\n============================================================\n');
fprintf('STARTING THE final U.S. GROUPED-DATA WORKFLOW\n');
fprintf('No individual U.S. chain sizes will be invented.\n');
fprintf('Primary analysis: completed years 2020-2023.\n');
fprintf('Sensitivity analysis: 2020 through March 28, 2024.\n');
fprintf('R_eff is estimated inside the prespecified subcritical model.\n');
fprintf('Primary detection priors: p_p ~ Beta(2,2), p_a ~ Beta(2,2).\n');
fprintf('Validation datasets per scenario = %d.\n',cfg.validationReplicates);
fprintf('\nCHECKPOINTS ARE STORED IN:\n%s\n',checkpointFolder);
fprintf('\nALL OUTPUTS FROM THIS RUN WILL BE SAVED IN:\n%s\n',runFolder);
fprintf('============================================================\n\n');

% Copy the exact running code immediately when possible.
sourceCode = [mfilename('fullpath') '.m'];
if exist(sourceCode,'file')
    copyfile(sourceCode,fullfile(runFolder,'Run_US_Grouped_MAP_Analysis.m'));
end

%% ------------------------------------------------------------------------
% 3. PUBLISHED GROUPED DATA
% Columns: year, chains, cases, size1, size2, size3plus
% -------------------------------------------------------------------------
allData = [ ...
    2020,  9,  13,  8, 0, 1; ...
    2021,  7,  49,  2, 2, 3; ...
    2022, 15, 121,  9, 1, 5; ...
    2023, 31,  58, 24, 3, 4; ...
    2024, 30,  97, 19, 4, 7];

primaryData = allData(1:4,:);
fullSnapshotData = allData;

largeCases = allData(:,3)-allData(:,4)-2*allData(:,5);
dataTable = array2table([allData largeCases], ...
    'VariableNames',{'Year','ObservedChains','TotalCases','Size1Chains', ...
    'Size2Chains','Size3PlusChains','CasesInSize3PlusChains'});
writetable(dataTable,fullfile(runFolder,'US_Grouped_Data.csv'));

% Save the exact analysis priors in a plain table. These are modeling
% choices, not claims that a historical study directly measured p_p or p_a.
Parameter={'R_eff';'k';'p_p';'p_a'};
PrimaryPrior={'Uniform(0,1)'; ...
    'log(k) Uniform(log(0.001),log(1000))'; ...
    'Beta(2,2)';'Beta(2,2)'};
Interpretation={ ...
    'Prespecified subcritical model support'; ...
    'Scale-neutral dispersion prior over six orders of magnitude'; ...
    'Transparent weak regularization; not an empirical ascertainment prior'; ...
    'Transparent weak regularization; not empirical contact-tracing evidence'};
priorTable=table(Parameter,PrimaryPrior,Interpretation);
writetable(priorTable,fullfile(runFolder,'US_Prior_Specification.csv'));

EvidenceType={ ...
    'Overall case ascertainment or reporting completeness'; ...
    'Passive routine detection probability'; ...
    'Active detection among passively missed infections after activation'; ...
    'Combined probability that an infection is ultimately observed'};
ModelQuantity={'Not automatically p_p';'p_p';'p_a'; ...
    'p_obs = p_p + (1-p_p)p_a'};
RequiredDenominator={ ...
    'All infections in a defined source population'; ...
    'All infections before active follow-up'; ...
    'Passively missed infections eligible for active follow-up'; ...
    'All infections under the stated two-stage surveillance mechanism'};
UseInPrimaryAnalysis={ ...
    'Context only unless the measurement process matches the model'; ...
    'No direct empirical prior available; use Beta(2,2)'; ...
    'No genuine matching contact-tracing prior available; use Beta(2,2)'; ...
    'Derived from p_p and p_a; not supplied as a separate prior'};
evidenceTable=table(EvidenceType,ModelQuantity,RequiredDenominator, ...
    UseInPrimaryAnalysis);
writetable(evidenceTable,fullfile(runFolder, ...
    'Surveillance_Evidence_Classification.csv'));

validate_grouped_data(allData);

fprintf('Published grouped-data arithmetic checks passed.\n');
fprintf('Primary dataset: %d chains and %d cases.\n', ...
    sum(primaryData(:,2)),sum(primaryData(:,3)));
fprintf('Full snapshot: %d chains and %d cases.\n\n', ...
    sum(fullSnapshotData(:,2)),sum(fullSnapshotData(:,3)));

%% ------------------------------------------------------------------------
% 4. START PARALLEL POOL IF AVAILABLE
% -------------------------------------------------------------------------
useParallel = false;
if cfg.useParallel && license('test','Distrib_Computing_Toolbox')
    try
        pool = gcp('nocreate');
        if isempty(pool)
            parpool('local');
        end
        useParallel = true;
        fprintf('Parallel processing is available.\n\n');
    catch ME
        fprintf('Parallel pool was not started: %s\n',ME.message);
        fprintf('The workflow will continue serially.\n\n');
    end
else
    fprintf('The workflow will run serially.\n\n');
end

%% ------------------------------------------------------------------------
% 5. MATHEMATICAL UNIT TESTS
% -------------------------------------------------------------------------
unitCheckpoint = fullfile(checkpointFolder,'Checkpoint_Unit_Tests.mat');
if checkpoint_matches(unitCheckpoint,cfg.version)
    loaded = load(unitCheckpoint,'unitResults');
    unitResults = loaded.unitResults;
    fprintf('Mathematical unit tests loaded from their completed checkpoint.\n\n');
else
    fprintf('RUNNING MATHEMATICAL UNIT TESTS...\n');
    unitResults = run_unit_tests(cfg);
    version = cfg.version; %#ok<NASGU>
    save(unitCheckpoint,'unitResults','version','-v7.3');
    fprintf('All mathematical unit tests passed.\n\n');
end

%% ------------------------------------------------------------------------
% 6. REALISTIC SIMULATION VALIDATION
% -------------------------------------------------------------------------
validationCheckpoint = fullfile(checkpointFolder, ...
    'Checkpoint_Simulation_Validation.mat');

validationScenarios = make_validation_scenarios(cfg);
validation = load_or_initialize_validation(validationCheckpoint, ...
    cfg,validationScenarios);

fprintf('REALISTIC GROUPED-DATA SIMULATION VALIDATION\n');
fprintf('The annual sample sizes are fixed at 9, 7, 15, and 31 chains.\n');
fprintf('No favorable dataset is selected.\n\n');

for sIndex = 1:numel(validationScenarios)
    scenario = validationScenarios(sIndex);
    completed = validation.completed(sIndex);
    if completed >= cfg.validationReplicates
        fprintf('Validation scenario %d/%d already complete: %s\n', ...
            sIndex,numel(validationScenarios),scenario.name);
        continue;
    end

    fprintf('\nScenario %d/%d: %s\n', ...
        sIndex,numel(validationScenarios),scenario.name);
    scenarioTimer=tic;
    fprintf('Truth: R_eff=%.2f, k=%.2f, p_p=%.2f, p_a=%.2f.\n', ...
        scenario.truth);

    validationDrawsInScenario=cfg.validationReplicates*sum(primaryData(:,2));
    samplerQ = build_sampling_distribution(scenario.truth, ...
        cfg.OmaxSamplerInitial,validationDrawsInScenario,cfg);

    while validation.completed(sIndex) < cfg.validationReplicates
        firstRep = validation.completed(sIndex)+1;
        lastRep = min(firstRep+cfg.validationBatchSize-1, ...
            cfg.validationReplicates);
        batchIndices = firstRep:lastRep;
        batchN = numel(batchIndices);
        batchResults = cell(batchN,1);

        if useParallel
            parfor b = 1:batchN
                rep = batchIndices(b);
                repSeed = cfg.masterSeed + 100000*sIndex + rep;
                batchResults{b} = run_one_validation_replicate( ...
                    scenario,primaryData,samplerQ,cfg,repSeed);
            end
        else
            for b = 1:batchN
                rep = batchIndices(b);
                repSeed = cfg.masterSeed + 100000*sIndex + rep;
                batchResults{b} = run_one_validation_replicate( ...
                    scenario,primaryData,samplerQ,cfg,repSeed);
            end
        end

        for b = 1:batchN
            rep = batchIndices(b);
            validation.results{sIndex}(rep) = batchResults{b};
        end
        validation.completed(sIndex) = lastRep;
        version = cfg.version; %#ok<NASGU>
        save(validationCheckpoint,'validation','validationScenarios', ...
            'version','-v7.3');
        fprintf('  checkpoint saved through replicate %d of %d.\n', ...
            lastRep,cfg.validationReplicates);
    end
    fprintf('  scenario elapsed time = %.1f minutes.\n',toc(scenarioTimer)/60);
end

validationSummary = summarize_validation(validation,validationScenarios,cfg);
writetable(validationSummary,fullfile(runFolder, ...
    'Simulation_Validation_Summary.csv'));
[validationDiagnosticsPass,validationRequirements] = ...
    evaluate_validation_requirements(validationSummary,validationScenarios,cfg);
writetable(validationRequirements,fullfile(runFolder, ...
    'Simulation_Validation_Requirements.csv'));
fprintf('Prespecified simulation-validation requirements passed: %s.\n', ...
    yes_no(validationDiagnosticsPass));

%% ------------------------------------------------------------------------
% 7. PRIMARY U.S. FIT: COMPLETED YEARS 2020-2023
% -------------------------------------------------------------------------
primaryCheckpoint = fullfile(checkpointFolder,'Checkpoint_US_Primary_Fit.mat');
fprintf('\nPRIMARY U.S. GROUPED-DATA FIT: 2020-2023\n');
if checkpoint_matches(primaryCheckpoint,cfg.version)
    loaded = load(primaryCheckpoint,'primaryFit');
    primaryFit = loaded.primaryFit;
    fprintf('Completed primary fit loaded from checkpoint.\n');
else
    primaryPrior = cfg.primaryPrior;
    primaryFit = fit_grouped_posterior(primaryData,primaryPrior,cfg, ...
        cfg.masterSeed+700001,useParallel);
    final_probability_check(primaryFit.map,primaryData,cfg);
    version = cfg.version; %#ok<NASGU>
    save(primaryCheckpoint,'primaryFit','primaryPrior','version','-v7.3');
    fprintf('Primary fit completed and checkpointed.\n');
end
display_fit_summary('2020-2023 primary analysis',primaryFit);

%% ------------------------------------------------------------------------
% 8. FULL SNAPSHOT SENSITIVITY: THROUGH MARCH 28, 2024
% -------------------------------------------------------------------------
fullCheckpoint = fullfile(checkpointFolder,'Checkpoint_US_Full_Snapshot_Fit.mat');
fprintf('\nFULL-SNAPSHOT SENSITIVITY FIT: 2020-MARCH 28, 2024\n');
if checkpoint_matches(fullCheckpoint,cfg.version)
    loaded = load(fullCheckpoint,'fullSnapshotFit');
    fullSnapshotFit = loaded.fullSnapshotFit;
    fprintf('Completed full-snapshot fit loaded from checkpoint.\n');
else
    primaryPrior = cfg.primaryPrior;
    fullSnapshotFit = fit_grouped_posterior(fullSnapshotData,primaryPrior, ...
        cfg,cfg.masterSeed+700002,useParallel);
    final_probability_check(fullSnapshotFit.map,fullSnapshotData,cfg);
    version = cfg.version; %#ok<NASGU>
    save(fullCheckpoint,'fullSnapshotFit','primaryPrior','version','-v7.3');
    fprintf('Full-snapshot fit completed and checkpointed.\n');
end
display_fit_summary('2020-March 28, 2024 snapshot',fullSnapshotFit);

%% ------------------------------------------------------------------------
% 9. prespecified PRIOR-LOCATION SENSITIVITY
% -------------------------------------------------------------------------
fprintf('\nprespecified DETECTION-PRIOR SENSITIVITY\n');
fprintf('The 2020-2023 grouped data remain unchanged in all nine fits.\n');
fprintf('Only the p_p and p_a prior means change.\n\n');

means = cfg.sensitivityMeans;
sensitivityFits = cell(numel(means),numel(means));
for iPP = 1:numel(means)
    for iPA = 1:numel(means)
        prior = make_detection_prior(means(iPP),means(iPA), ...
            cfg.sensitivityConcentration);
        sensitivityFile = fullfile(checkpointFolder,sprintf( ...
            'Checkpoint_US_PriorSensitivity_pp%02d_pa%02d.mat', ...
            round(100*means(iPP)),round(100*means(iPA))));

        fprintf('Prior means p_p=%.2f, p_a=%.2f: ',means(iPP),means(iPA));
        if checkpoint_matches(sensitivityFile,cfg.version)
            loaded = load(sensitivityFile,'oneFit');
            oneFit = loaded.oneFit;
            fprintf('loaded from checkpoint.\n');
        else
            oneFit = fit_prior_sensitivity(primaryFit,primaryData,prior,cfg, ...
                cfg.masterSeed+800000+100*iPP+iPA,useParallel);
            final_probability_check(oneFit.map,primaryData,cfg);
            version = cfg.version; %#ok<NASGU>
            save(sensitivityFile,'oneFit','prior','version','-v7.3');
            fprintf('completed and checkpointed.\n');
        end
        sensitivityFits{iPA,iPP} = oneFit;
    end
end

sensitivityTable = make_sensitivity_table(sensitivityFits,means);
writetable(sensitivityTable,fullfile(runFolder, ...
    'US_Detection_Prior_Sensitivity.csv'));

%% ------------------------------------------------------------------------
% 10. POSTERIOR-PREDICTIVE CHECKS
% -------------------------------------------------------------------------
ppcCheckpoint = fullfile(checkpointFolder, ...
    'Checkpoint_Posterior_Predictive_Checks.mat');
fprintf('\nPOSTERIOR-PREDICTIVE CHECKS\n');
if checkpoint_matches(ppcCheckpoint,cfg.version)
    loaded = load(ppcCheckpoint,'primaryPPC','fullPPC');
    primaryPPC = loaded.primaryPPC;
    fullPPC = loaded.fullPPC;
    fprintf('Posterior-predictive checks loaded from checkpoint.\n');
else
    primaryPPC = run_posterior_predictive(primaryFit,primaryData,cfg, ...
        cfg.masterSeed+900001,useParallel);
    fullPPC = run_posterior_predictive(fullSnapshotFit,fullSnapshotData,cfg, ...
        cfg.masterSeed+900002,useParallel);
    version = cfg.version; %#ok<NASGU>
    save(ppcCheckpoint,'primaryPPC','fullPPC','version','-v7.3');
    fprintf('Posterior-predictive checks completed and checkpointed.\n');
end

writetable(ppc_to_table(primaryPPC,'Primary_2020_2023'), ...
    fullfile(runFolder,'US_Primary_Posterior_Predictive_Check.csv'));
writetable(ppc_to_table(fullPPC,'Full_Through_March_2024'), ...
    fullfile(runFolder,'US_FullSnapshot_Posterior_Predictive_Check.csv'));

%% ------------------------------------------------------------------------
% 11. SAVE NUMERICAL SUMMARIES
% -------------------------------------------------------------------------
primarySummaryTable = fit_to_table(primaryFit,'Primary_2020_2023');
fullSummaryTable = fit_to_table(fullSnapshotFit,'Full_Through_March_2024');
writetable(primarySummaryTable,fullfile(runFolder, ...
    'US_Primary_2020_2023_Posterior_Summary.csv'));
writetable(fullSummaryTable,fullfile(runFolder, ...
    'US_FullSnapshot_Posterior_Summary.csv'));

writematrix(primaryFit.samples,fullfile(runFolder, ...
    'US_Primary_2020_2023_Posterior_Samples.csv'));
writematrix(fullSnapshotFit.samples,fullfile(runFolder, ...
    'US_FullSnapshot_Posterior_Samples.csv'));

%% ------------------------------------------------------------------------
% 12. DRAW CLEAN FIGURES
% -------------------------------------------------------------------------
fprintf('\nDRAWING CLEAN FIGURES...\n');

fig8 = draw_figure8(primaryFit,fullSnapshotFit);
save_publication_figure(fig8,runFolder,'Figure 8',4800,3000);
close(fig8);

fig9 = draw_figure9(sensitivityFits,means);
save_publication_figure(fig9,runFolder,'Figure 9',4800,3000);
close(fig9);

figS1 = draw_figureS1(validationSummary);
save_publication_figure(figS1,runFolder,'Figure S1',4200,2600);
close(figS1);

figS2 = draw_figureS2(primaryPPC,fullPPC);
save_publication_figure(figS2,runFolder,'Figure S2',4200,2600);
close(figS2);

%% ------------------------------------------------------------------------
% 13. MASTER MAT, README, AND FINAL DIAGNOSTIC DECISION
% -------------------------------------------------------------------------
masterMat = fullfile(runFolder,'US_Grouped_MAP_methodologicalFinal_All_Results.mat');
save(masterMat,'cfg','allData','primaryData','fullSnapshotData', ...
    'unitResults','validation','validationSummary','primaryFit', ...
    'fullSnapshotFit','sensitivityFits','sensitivityTable', ...
    'validationDiagnosticsPass','validationRequirements', ...
    'primaryPPC','fullPPC','priorTable','evidenceTable','-v7.3');

write_readme(runFolder,cfg,primaryFit,fullSnapshotFit);

allDiagnosticsPass = primaryFit.diagnosticsPass && ...
    fullSnapshotFit.diagnosticsPass && ...
    primaryPPC.pass && fullPPC.pass && validationDiagnosticsPass;
for ii = 1:numel(sensitivityFits)
    allDiagnosticsPass = allDiagnosticsPass && ...
        sensitivityFits{ii}.diagnosticsPass;
end

fprintf('\n============================================================\n');
fprintf('DONE: THE GROUPED-DATA WORKFLOW COMPLETED.\n');
fprintf('No unavailable individual U.S. chain size was invented.\n');
fprintf('Primary analysis used 62 completed-period chains from 2020-2023.\n');
fprintf('The 92-chain dataset was used only as a prespecified snapshot sensitivity.\n');
fprintf('All required convergence and predictive checks passed: %s\n', ...
    yes_no(allDiagnosticsPass));
fprintf('\nALL OUTPUTS WERE SAVED IN:\n%s\n',runFolder);
fprintf('\nFIGURE 8: %s\n',fullfile(runFolder,'Figure 8.pdf'));
fprintf('FIGURE 9: %s\n',fullfile(runFolder,'Figure 9.pdf'));
fprintf('SUPPLEMENTARY FIGURE S1: %s\n',fullfile(runFolder,'Figure S1.pdf'));
fprintf('SUPPLEMENTARY FIGURE S2: %s\n',fullfile(runFolder,'Figure S2.pdf'));
fprintf('MASTER MAT FILE: %s\n',masterMat);
fprintf('COMMAND-WINDOW LOG: %s\n',logFile);
fprintf('TOTAL ELAPSED TIME: %.1f minutes.\n',toc(workflowTimer)/60);
fprintf('============================================================\n');

if cfg.stopOnFailedDiagnostics && ~allDiagnosticsPass
    error(['One or more final diagnostics failed. Numerical outputs were ' ...
        'saved for investigation, but they must not yet be used in the manuscript.']);
elseif ~allDiagnosticsPass
    warning(['One or more diagnostics did not meet the prespecified threshold. ' ...
        'Send the command-window log and CSV files for review before manuscript use.']);
end

end

%% ========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function validate_grouped_data(data)
if any(data(:,2) ~= data(:,4)+data(:,5)+data(:,6))
    error('Grouped data error: chain categories do not sum to total chains.');
end
largeCases = data(:,3)-data(:,4)-2*data(:,5);
if any(data(:,6)==0 & largeCases~=0)
    error('Grouped data error: positive large-chain cases with no large chains.');
end
if any(data(:,6)>0 & largeCases<3*data(:,6))
    error('Grouped data error: size>=3 chains contain too few cases.');
end
end

function tf = checkpoint_matches(fileName,versionExpected)
tf = false;
if ~exist(fileName,'file'); return; end
try
    x = load(fileName,'version');
    tf = isfield(x,'version') && strcmp(x.version,versionExpected);
catch
    tf = false;
end
end

function unit = run_unit_tests(cfg)
unit = struct;

% Test 1: the two-case hand calculation in the manuscript specification.
pp = 0.5; pa = 0.5; s = 2;
p = observation_probabilities_for_one_true_size(s,pp,pa);
expected = [0.25 0.25 0.50];
unit.twoCaseMaximumError = max(abs(p-expected));
if unit.twoCaseMaximumError > cfg.probabilityTolerance
    error('Unit test failed: the two-case surveillance example is wrong.');
end

% Test 2: probabilities sum to one over several settings and true sizes.
settings = [0.25 0.25; 0.25 0.75; 0.50 0.50; 0.75 1.00; 1.00 0.25];
maxSumError = 0;
minProbability = Inf;
for ii = 1:size(settings,1)
    for ss = 1:12
        oneP = observation_probabilities_for_one_true_size( ...
            ss,settings(ii,1),settings(ii,2));
        maxSumError = max(maxSumError,abs(sum(oneP)-1));
        minProbability = min(minProbability,min(oneP));
    end
end
unit.maximumSurveillanceSumError = maxSumError;
unit.minimumSurveillanceProbability = minProbability;
if maxSumError > cfg.probabilityTolerance || minProbability < -cfg.probabilityTolerance
    error('Unit test failed: surveillance probabilities are invalid.');
end

% Test 3: p_p=1 gives perfect observation and makes p_a irrelevant.
for ss = 1:10
    pLow = observation_probabilities_for_one_true_size(ss,1,0.1);
    pHigh = observation_probabilities_for_one_true_size(ss,1,0.9);
    target = zeros(1,ss+1); target(ss+1)=1;
    if max(abs(pLow-target))>cfg.probabilityTolerance || ...
            max(abs(pHigh-target))>cfg.probabilityTolerance
        error('Unit test failed: p_p=1 is not perfect observation.');
    end
end

% Test 4: p_a=0 agrees with ordinary binomial passive thinning.
for ss = 1:10
    pp = 0.37;
    modelP = observation_probabilities_for_one_true_size(ss,pp,0);
    binomP = zeros(1,ss+1);
    for oo = 0:ss
        binomP(oo+1) = exp(gammaln(ss+1)-gammaln(oo+1)- ...
            gammaln(ss-oo+1)+oo*log(pp)+(ss-oo)*log(1-pp));
    end
    if max(abs(modelP-binomP))>1e-9
        error('Unit test failed: p_a=0 does not reduce to binomial thinning.');
    end
end

% Test 5: true total-progeny probabilities are properly normalized.
g = true_cluster_pmf(0.5,0.15,20000);
unit.trueClusterMassAt20000 = sum(g);
if abs(1-sum(g)) > 1e-8
    error('Unit test failed: true cluster probabilities are not normalized.');
end

% Test 6: observed q probabilities sum to one when enough sizes are retained.
theta = [0.5 0.15 0.5 0.5];
[q,diagQ] = observed_cluster_pmf(theta,1000,5000);
unit.observedMassAt1000 = sum(q);
unit.observationNumeratorTail = diagQ.omittedObservedMass;
if abs(1-sum(q)) > 1e-7
    error('Unit test failed: conditional observed-size probabilities do not sum to one.');
end

% Test 7: the fast exact-series result agrees with an independent direct
% sum over hidden true-cluster sizes at the sizes used by the U.S. data.
[qSeries,~]=observed_cluster_pmf(theta,110,0);
[qDirect,dDirect]=observed_cluster_pmf_direct( ...
    theta,110,cfg.SmaxDirectCheck);
unit.maximumSeriesVersusDirectError=max(abs(qSeries-qDirect));
unit.directTrueSizeTail=1-dDirect.trueMass;
if unit.maximumSeriesVersusDirectError>cfg.finalLikelihoodTolerance
    error('Unit test failed: exact series and independent direct sum disagree.');
end

% Test 8: exact tiny grouped likelihood equals 0.0225.
qTiny = zeros(1,8);
qTiny(1)=0.50; qTiny(2)=0.30; qTiny(3)=0.15; qTiny(4)=0.05;
tinyData = [2000 3 8 1 0 2];
logLTiny = grouped_loglik_from_q(tinyData,qTiny);
unit.tinyGroupedLikelihood = exp(logLTiny);
if abs(exp(logLTiny)-0.0225)>1e-12
    error('Unit test failed: tiny grouped convolution is wrong.');
end

unit.pass = true;
end

function p = observation_probabilities_for_one_true_size(s,pp,pa)
p = zeros(1,s+1);
p(1) = (1-pp)^s;
z = pp+(1-pp)*pa;
d = (1-pp)*pa;
a = 1-z;
for o = 1:s
    logChoose = gammaln(s+1)-gammaln(o+1)-gammaln(s-o+1);
    if a==0 && s>o
        value = 0;
    else
        value = exp(logChoose) * a^(s-o) * (z^o-d^o);
    end
    p(o+1)=max(0,value);
end
end

function g = true_cluster_pmf(R,k,Smax)
s = (1:Smax)';
logg = -log(s) + gammaln((k+1).*s-1) - gammaln(k.*s) - ...
    gammaln(s) + (s-1).*log(R/k) - ((k+1).*s-1).*log1p(R/k);
g = exp(logg);
end

function value = total_progeny_pgf(x,R,k)
if x<=0; value=0; return; end
if x>=1; value=1; return; end
lo=0; hi=1;
for iter=1:80
    mid=(lo+hi)/2;
    fmid = mid-x*(1+(R/k)*(1-mid))^(-k);
    if fmid>0
        hi=mid;
    else
        lo=mid;
    end
end
value=(lo+hi)/2;
end

function [q,diagnostics] = observed_cluster_pmf(theta,Omax,unusedSmax) %#ok<INUSD>
R=theta(1); k=theta(2); pp=theta(3); pa=theta(4);
if ~(R>0 && R<1 && k>0 && pp>0 && pp<=1 && pa>=0 && pa<=1)
    q = nan(1,Omax);
    diagnostics = struct('omittedObservedMass',Inf,'trueMass',NaN, ...
        'probabilityObserved',NaN);
    return;
end

% Let T(x) be the PGF of the true final cluster size. The numerator PGF
% for positive observed cluster sizes is
%
%       T(a + z*t) - T(a + d*t),
%
% where z=p_p+(1-p_p)p_a, d=(1-p_p)p_a, and a=1-z. The coefficients are
% obtained through order Omax by formal power-series Newton iteration.
% This is mathematically equivalent to summing over every possible hidden
% true cluster size, but it is faster and uses no hidden-size cutoff.
z = pp+(1-pp)*pa;
d = (1-pp)*pa;
a = 1-z;

coefficientsZ = composed_total_progeny_series(a,z,R,k,Omax);
ratio = d/z;
orders = 0:Omax;
numeratorAll = coefficientsZ .* (1-ratio.^orders);
numerator = numeratorAll(2:end);

m0 = total_progeny_pgf(1-pp,R,k);
pObserved = 1-m0;
q = numerator./pObserved;
q(q<0 & q>-1e-12)=0;

diagnostics.trueMass = 1;
diagnostics.probabilityObserved = pObserved;
diagnostics.retainedObservedMass = sum(q);
diagnostics.omittedObservedMass = max(0,1-sum(q));
diagnostics.algorithm = 'exact truncated PGF power series';
diagnostics.hiddenSizeCutoffUsed = false;
end

function y = composed_total_progeny_series(x0,x1,R,k,N)
% Coefficients of T(x0+x1*t) through t^N, where
% T=x{1+(R/k)(1-T)}^{-k}. Formal Newton iteration doubles the number of
% correct coefficients at each pass. Coefficient n is stored at n+1.
x=zeros(1,N+1); x(1)=x0;
if N>=1; x(2)=x1; end

y=zeros(1,N+1);
y(1)=total_progeny_pgf(x0,R,k);
c=R/k;
numberOfNewtonSteps=ceil(log2(N+1))+3;
for iter=1:numberOfNewtonSteps
    D=zeros(1,N+1); D(1)=1+c; D=D-c*y;
    phi=series_power(D,-k,N);
    F=y-series_multiply(x,phi,N);

    % d/dT of {T-x*phi(T)} is 1-x*R*D^{-k-1}.
    derivative=zeros(1,N+1); derivative(1)=1;
    derivative=derivative-series_multiply( ...
        x,R*series_power(D,-k-1,N),N);
    correction=series_multiply(F,series_power(derivative,-1,N),N);
    y=y-correction;
end
y(abs(y)<1e-15)=0;
end

function c = series_multiply(a,b,N)
c=conv(a,b);
if numel(c)>N+1; c=c(1:N+1); end
if numel(c)<N+1; c=[c zeros(1,N+1-numel(c))]; end %#ok<AGROW>
end

function b = series_power(a,p,N)
% If B(t)=A(t)^p, this recurrence computes B_0,...,B_N up to ordinary
% floating-point rounding. A(0) must be positive.
a=a(:)';
if numel(a)<N+1; a=[a zeros(1,N+1-numel(a))]; end %#ok<AGROW>
if ~(isfinite(a(1)) && a(1)>0)
    b=nan(1,N+1); return;
end
b=zeros(1,N+1); b(1)=a(1)^p;
for n=1:N
    j=1:n;
    weights=((p+1).*j-n).*a(j+1);
    previous=b(n-j+1);
    b(n+1)=sum(weights.*previous)/(n*a(1));
end
end

function [q,diagnostics] = observed_cluster_pmf_direct(theta,Omax,Smax)
% Independent slow calculation used only to verify the fast PGF recursion.
R=theta(1); k=theta(2); pp=theta(3); pa=theta(4);
g=true_cluster_pmf(R,k,Smax);
s=(1:Smax)';
z=pp+(1-pp)*pa; d=(1-pp)*pa; a=1-z;
numerator=zeros(1,Omax);
if a<1e-14
    maximum=min(Omax,Smax);
    for o=1:maximum
        numerator(o)=g(o)*(z^o-d^o);
    end
else
    wz=g.*(a.^s); wd=wz;
    for o=1:Omax
        multiplier=max(s-(o-1),0)./o;
        wz=wz.*multiplier*(z/a);
        wd=wd.*multiplier*(d/a);
        numerator(o)=sum(wz-wd);
    end
end
m0=total_progeny_pgf(1-pp,R,k);
q=numerator/(1-m0);
q(q<0 & q>-1e-12)=0;
diagnostics.trueMass=sum(g);
diagnostics.retainedObservedMass=sum(q);
diagnostics.omittedObservedMass=max(0,1-sum(q));
diagnostics.probabilityObserved=1-m0;
end

function logL = grouped_loglik(theta,data,cfg)
if any(~isfinite(theta)) || theta(1)<=0 || theta(1)>=1 || ...
        theta(2)<=0 || theta(3)<=0 || theta(3)>1 || ...
        theta(4)<0 || theta(4)>1
    logL=-Inf; return;
end
maximumCasesInLargeGroup = max(data(:,3)-data(:,4)-2*data(:,5));
[q,~] = observed_cluster_pmf(theta,maximumCasesInLargeGroup, ...
    0);
if any(~isfinite(q)) || any(q<0)
    logL=-Inf; return;
end
logL = grouped_loglik_from_q(data,q);
end

function logL = grouped_loglik_from_q(data,q)
logL=0;
for yy=1:size(data,1)
    n=data(yy,2); totalCases=data(yy,3);
    n1=data(yy,4); n2=data(yy,5); n3=data(yy,6);
    C=totalCases-n1-2*n2;
    if q(1)<=0 || (n2>0 && q(2)<=0)
        logL=-Inf; return;
    end
    logChoose=gammaln(n+1)-gammaln(n1+1)-gammaln(n2+1)-gammaln(n3+1);
    if n3==0
        if C~=0; logL=-Inf; return; end
        logG=0;
    else
        base=zeros(1,C+1);
        maxO=min(C,numel(q));
        if maxO>=3; base(4:maxO+1)=q(3:maxO); end
        polynomial=1;
        for rr=1:n3
            polynomial=conv(polynomial,base);
            if numel(polynomial)>C+1
                polynomial=polynomial(1:C+1);
            end
        end
        if numel(polynomial)<C+1 || polynomial(C+1)<=0
            logL=-Inf; return;
        end
        logG=log(polynomial(C+1));
    end
    logL=logL+logChoose+n1*log(q(1))+n2*log(max(q(2),realmin))+logG;
end
end

function scenarios = make_validation_scenarios(cfg)
weak = cfg.primaryPrior;
low = make_detection_prior(0.20,0.20,cfg.sensitivityConcentration);
scenarios(1)=struct('name','Moderate transmission; detection known', ...
    'truth',[0.50 0.15 0.50 0.50],'knownDetection',true,'prior',weak);
scenarios(2)=struct('name','Near-threshold transmission; detection known', ...
    'truth',[0.85 0.15 0.50 0.50],'knownDetection',true,'prior',weak);
scenarios(3)=struct('name','Moderate transmission; weak detection priors', ...
    'truth',[0.50 0.15 0.50 0.50],'knownDetection',false,'prior',weak);
scenarios(4)=struct('name','Moderate transmission; misspecified low detection priors', ...
    'truth',[0.50 0.15 0.50 0.50],'knownDetection',false,'prior',low);
end

function validation = load_or_initialize_validation(fileName,cfg,scenarios)
if checkpoint_matches(fileName,cfg.version)
    x=load(fileName,'validation');
    validation=x.validation;
    if numel(validation.completed)~=numel(scenarios) || ...
            validation.targetReplicates~=cfg.validationReplicates
        validation=initialize_validation(cfg,scenarios);
    end
else
    validation=initialize_validation(cfg,scenarios);
end
end

function validation=initialize_validation(cfg,scenarios)
emptyResult=struct('estimate',nan(1,4),'ciLow',nan(1,4), ...
    'ciHigh',nan(1,4),'success',false,'logPosterior',-Inf);
validation.completed=zeros(numel(scenarios),1);
validation.targetReplicates=cfg.validationReplicates;
validation.results=cell(numel(scenarios),1);
for ss=1:numel(scenarios)
    validation.results{ss}=repmat(emptyResult,cfg.validationReplicates,1);
end
end

function q = build_sampling_distribution(theta,Ostart,numberOfDraws,cfg)
Omax=Ostart;
while true
    [q,diagQ]=observed_cluster_pmf(theta,Omax,0);
    expectedOmitted=diagQ.omittedObservedMass*numberOfDraws;
    if expectedOmitted<=cfg.maximumExpectedOmittedValidationDraws
        fprintf(['  validation sampling PMF retained through size %d; ' ...
            'expected omitted draws = %.4g.\n'],Omax,expectedOmitted);
        q=q/sum(q);
        return;
    end
    if Omax>=cfg.maxSamplerO
        error(['Validation sampling PMF still expects %.4g omitted draws ' ...
            'after size %d. Increase cfg.maxSamplerO.'], ...
            expectedOmitted,Omax);
    end
    Omax=min(2*Omax,cfg.maxSamplerO);
end
end

function result = run_one_validation_replicate(scenario,templateData,q,cfg,seed)
rng(seed,'twister');
simData=simulate_grouped_data(templateData,q);
result=fit_grouped_laplace(simData,scenario.prior,cfg, ...
    scenario.knownDetection,scenario.truth(3:4),seed+17);
end

function simulatedData = simulate_grouped_data(templateData,q)
simulatedData=templateData;
cdf=cumsum(q(:)); cdf(end)=1;
for yy=1:size(templateData,1)
    n=templateData(yy,2);
    u=rand(n,1);
    sizes=zeros(n,1);
    for ii=1:n
        sizes(ii)=find(cdf>=u(ii),1,'first');
    end
    simulatedData(yy,3)=sum(sizes);
    simulatedData(yy,4)=sum(sizes==1);
    simulatedData(yy,5)=sum(sizes==2);
    simulatedData(yy,6)=sum(sizes>=3);
end
end

function result=fit_grouped_laplace(data,prior,cfg,knownDetection,fixedDetection,seed)
rng(seed,'twister');
if knownDetection
    dimension=2;
else
    dimension=4;
end
target=@(u) logposterior_u_for_sampling(u,data,prior,cfg, ...
    knownDetection,fixedDetection);
validationCfg=cfg;
validationCfg.numberOfStarts=cfg.validationNumberOfStarts;
[uMode,bestValue,success]=multistart_optimize( ...
    target,dimension,validationCfg,seed);

result.estimate=nan(1,4); result.ciLow=nan(1,4); result.ciHigh=nan(1,4);
result.success=success; result.logPosterior=bestValue;
if ~success; return; end

H=numerical_hessian(@(u) -target(u),uMode);
H=(H+H')/2;
[V,D]=eig(H);
d=diag(D); d(d<1e-6)=1e-6;
covU=V*diag(1./d)*V'; covU=(covU+covU')/2;
sd=sqrt(max(diag(covU),0));

thetaMode=decode_parameters(uMode,knownDetection,fixedDetection);
lowU=uMode(:)-1.96*sd;
highU=uMode(:)+1.96*sd;
thetaLow=decode_parameters(lowU,knownDetection,fixedDetection);
thetaHigh=decode_parameters(highU,knownDetection,fixedDetection);

result.estimate=thetaMode;
result.ciLow=thetaLow;
result.ciHigh=thetaHigh;
if knownDetection
    result.estimate(3:4)=fixedDetection;
    result.ciLow(3:4)=NaN;
    result.ciHigh(3:4)=NaN;
end
end

function [bestU,bestLogPost,success]=multistart_optimize(target,dimension,cfg,seed)
rng(seed,'twister');
starts=zeros(cfg.numberOfStarts,dimension);
baseR=[0.30 0.55 0.80]; baseK=[0.05 0.20 1.0];
for ii=1:cfg.numberOfStarts
    R0=baseR(mod(ii-1,numel(baseR))+1);
    k0=baseK(mod(floor((ii-1)/numel(baseR)),numel(baseK))+1);
    if dimension==2
        starts(ii,:)=[logit(R0),log(k0)];
    else
        pp0=0.15+0.70*rand;
        pa0=0.15+0.70*rand;
        starts(ii,:)=[logit(R0),log(k0),logit(pp0),logit(pa0)];
    end
end
options=optimset('Display','off','MaxIter',cfg.optimMaxIter, ...
    'MaxFunEvals',cfg.optimMaxFunEvals,'TolX',1e-6,'TolFun',1e-6);
bestLogPost=-Inf; bestU=starts(1,:); success=false;
for ii=1:size(starts,1)
    objective=@(u) safe_negative(target,u);
    [u,fval,exitflag]=fminsearch(objective,starts(ii,:),options);
    value=-fval;
    if isfinite(value) && value>bestLogPost
        bestLogPost=value; bestU=u; success=exitflag~=0;
    end
end
success=success && isfinite(bestLogPost);
end

function value=safe_negative(target,u)
lp=target(u);
if ~isfinite(lp); value=1e100; else; value=-lp; end
end

function H=numerical_hessian(fun,u)
u=u(:); d=numel(u); H=zeros(d);
h=1e-3*(1+abs(u)); f0=fun(u);
for i=1:d
    ei=zeros(d,1); ei(i)=h(i);
    H(i,i)=(fun(u+ei)-2*f0+fun(u-ei))/h(i)^2;
    for j=i+1:d
        ej=zeros(d,1); ej(j)=h(j);
        H(i,j)=(fun(u+ei+ej)-fun(u+ei-ej)-fun(u-ei+ej)+fun(u-ei-ej)) ...
            /(4*h(i)*h(j));
        H(j,i)=H(i,j);
    end
end
end

function summary=summarize_validation(validation,scenarios,cfg)
names={}; parameters={}; bias=[]; mae=[]; rmse=[]; coverage=[]; failures=[];
paramNames={'R_eff','k','p_p','p_a'};
for ss=1:numel(scenarios)
    results=validation.results{ss};
    estimates=vertcat(results.estimate);
    lows=vertcat(results.ciLow); highs=vertcat(results.ciHigh);
    success=[results.success]';
    truth=scenarios(ss).truth;
    for pp=1:4
        valid=success & isfinite(estimates(:,pp));
        names{end+1,1}=scenarios(ss).name; %#ok<AGROW>
        parameters{end+1,1}=paramNames{pp}; %#ok<AGROW>
        if any(valid)
            errors=estimates(valid,pp)-truth(pp);
            bias(end+1,1)=mean(errors); %#ok<AGROW>
            mae(end+1,1)=mean(abs(errors)); %#ok<AGROW>
            rmse(end+1,1)=sqrt(mean(errors.^2)); %#ok<AGROW>
            intervalValid=valid & isfinite(lows(:,pp)) & isfinite(highs(:,pp));
            if any(intervalValid)
                coverage(end+1,1)=mean(lows(intervalValid,pp)<=truth(pp) & ...
                    highs(intervalValid,pp)>=truth(pp)); %#ok<AGROW>
            else
                coverage(end+1,1)=NaN; %#ok<AGROW>
            end
        else
            bias(end+1,1)=NaN; mae(end+1,1)=NaN; rmse(end+1,1)=NaN; %#ok<AGROW>
            coverage(end+1,1)=NaN; %#ok<AGROW>
        end
        failures(end+1,1)=sum(~success); %#ok<AGROW>
    end
end
summary=table(names,parameters,bias,mae,rmse,coverage,failures, ...
    repmat(cfg.validationReplicates,numel(names),1), ...
    'VariableNames',{'Scenario','Parameter','Bias','MAE','RMSE', ...
    'Approximate95Coverage','FailedFits','Replicates'});
end

function [pass,T]=evaluate_validation_requirements(summary,scenarios,cfg)
% Hard checks are deliberately limited to things the validation experiment
% is designed to establish. Detection recovery under weak or deliberately
% misspecified priors is reported, not forced to pass, because partial
% identifiability is one of the scientific questions.
Check={}; Requirement={}; Observed=[]; Passed=[];

maximumFailureRate=max(summary.FailedFits./summary.Replicates);
Check{end+1,1}='Optimization completion'; %#ok<AGROW>
Requirement{end+1,1}='At least 95% of fits complete in every scenario'; %#ok<AGROW>
Observed(end+1,1)=1-maximumFailureRate; %#ok<AGROW>
Passed(end+1,1)=maximumFailureRate<=0.05; %#ok<AGROW>

for ss=1:2
    scenarioName=scenarios(ss).name;
    for parameterName={'R_eff','k'}
        mask=strcmp(summary.Scenario,scenarioName) & ...
            strcmp(summary.Parameter,parameterName{1});
        coverageValue=summary.Approximate95Coverage(mask);
        Check{end+1,1}=sprintf('Known-detection coverage: scenario %d, %s', ...
            ss,parameterName{1}); %#ok<AGROW>
        Requirement{end+1,1}='Approximate 95% interval coverage at least 0.85'; %#ok<AGROW>
        Observed(end+1,1)=coverageValue; %#ok<AGROW>
        Passed(end+1,1)=isfinite(coverageValue) && coverageValue>=0.85; %#ok<AGROW>
    end
end

T=table(Check,Requirement,Observed,logical(Passed), ...
    'VariableNames',{'Check','Requirement','ObservedValue','Passed'});
pass=all(T.Passed) && height(T)==5 && cfg.validationReplicates>=200;
end

function prior=make_detection_prior(ppMean,paMean,concentration)
prior.ppAlpha=ppMean*concentration;
prior.ppBeta=(1-ppMean)*concentration;
prior.paAlpha=paMean*concentration;
prior.paBeta=(1-paMean)*concentration;
prior.name=sprintf('p_p mean %.2f; p_a mean %.2f; concentration %.1f', ...
    ppMean,paMean,concentration);
end

function fit=fit_prior_sensitivity(primaryFit,data,newPrior,cfg,seed,useParallel)
% The likelihood and data are unchanged, so posterior samples obtained
% under the primary prior can be reweighted exactly by
% prior_new/prior_primary. A dedicated MCMC analysis is used automatically
% if the importance weights show inadequate posterior overlap.
samples=primaryFit.samples;
pp=samples(:,3); pa=samples(:,4);
oldPrior=cfg.primaryPrior;

logNew=(newPrior.ppAlpha-1).*log(pp) + ...
    (newPrior.ppBeta-1).*log1p(-pp) - ...
    betaln(newPrior.ppAlpha,newPrior.ppBeta) + ...
    (newPrior.paAlpha-1).*log(pa) + ...
    (newPrior.paBeta-1).*log1p(-pa) - ...
    betaln(newPrior.paAlpha,newPrior.paBeta);
logOld=(oldPrior.ppAlpha-1).*log(pp) + ...
    (oldPrior.ppBeta-1).*log1p(-pp) - ...
    betaln(oldPrior.ppAlpha,oldPrior.ppBeta) + ...
    (oldPrior.paAlpha-1).*log(pa) + ...
    (oldPrior.paBeta-1).*log1p(-pa) - ...
    betaln(oldPrior.paAlpha,oldPrior.paBeta);

logWeights=logNew-logOld;
logWeights=logWeights-max(logWeights);
weights=exp(logWeights); weights=weights/sum(weights);
importanceESS=1/sum(weights.^2);
maximumWeight=max(weights);
adequate=importanceESS>=cfg.minimumImportanceESS && ...
    maximumWeight<=cfg.maximumImportanceWeight;

fprintf('importance ESS %.0f; largest weight %.4f; ', ...
    importanceESS,maximumWeight);
if ~adequate
    fprintf('overlap insufficient, running dedicated MCMC...\n');
    fit=fit_grouped_posterior(data,newPrior,cfg,seed,useParallel);
    fit.inferenceMethod='Dedicated MCMC because importance overlap was insufficient';
    fit.importanceESS=importanceESS;
    fit.importanceMaximumWeight=maximumWeight;
    return;
end

fprintf('importance reweighting accepted... ');
knownDetection=false; fixedDetection=[NaN NaN]; dimension=4;
target=@(u) logposterior_mixed_measure(u,data,newPrior,cfg, ...
    knownDetection,fixedDetection);
[uMap,logPostMap,success]=multistart_optimize(target,dimension,cfg,seed);
if ~success
    fprintf('MAP optimization failed, running dedicated MCMC...\n');
    fit=fit_grouped_posterior(data,newPrior,cfg,seed,useParallel);
    fit.inferenceMethod='Dedicated MCMC after sensitivity MAP optimization failed';
    fit.importanceESS=importanceESS;
    fit.importanceMaximumWeight=maximumWeight;
    return;
end

fit.parameterNames={'R_eff','k','p_p','p_a'};
fit.prior=newPrior;
fit.map=decode_parameters(uMap,false,fixedDetection);
fit.logPosteriorAtMap=logPostMap;
fit.samples=[];
fit.posteriorMedian=weighted_column_quantiles(samples,weights,0.50);
fit.ciLow=weighted_column_quantiles(samples,weights,0.025);
fit.ciHigh=weighted_column_quantiles(samples,weights,0.975);
fit.posteriorMean=sum(samples.*weights,1);
fit.rhat=nan(1,4);
fit.ess=repmat(importanceESS,1,4);
fit.acceptance=NaN;
fit.diagnosticsPass=true;
fit.chainSamples=[];
fit.inferenceMethod='Exact importance reweighting of primary posterior draws';
fit.importanceESS=importanceESS;
fit.importanceMaximumWeight=maximumWeight;
fprintf('accepted.\n');
end

function q=weighted_column_quantiles(x,weights,probability)
q=zeros(1,size(x,2));
for jj=1:size(x,2)
    [sorted,index]=sort(x(:,jj));
    cumulative=cumsum(weights(index));
    first=find(cumulative>=probability,1,'first');
    if isempty(first); first=numel(sorted); end
    q(jj)=sorted(first);
end
end

function fit=fit_grouped_posterior(data,prior,cfg,seed,useParallel)
knownDetection=false; fixedDetection=[NaN NaN]; dimension=4;

fitTimer=tic;
fprintf('  multistart posterior optimization is running...\n');
mixedTarget=@(u) logposterior_mixed_measure(u,data,prior,cfg, ...
    knownDetection,fixedDetection);
[uMap,logPostMap,success]=multistart_optimize(mixedTarget,dimension,cfg,seed);
if ~success
    error('All posterior optimizations failed.');
end
thetaMap=decode_parameters(uMap,false,fixedDetection);
fprintf(['  MAP optimization completed. Running %d independent MCMC ' ...
    'chains...\n'],cfg.numberOfChains);

samplingTarget=@(u) logposterior_u_for_sampling(u,data,prior,cfg, ...
    false,fixedDetection);

chainCells=cell(cfg.numberOfChains,1);
acceptance=zeros(cfg.numberOfChains,1);
if useParallel
    parfor cc=1:cfg.numberOfChains
        [chainCells{cc},acceptance(cc)]=run_adaptive_chain( ...
            samplingTarget,uMap,cfg,seed+1000*cc);
    end
else
    for cc=1:cfg.numberOfChains
        [chainCells{cc},acceptance(cc)]=run_adaptive_chain( ...
            samplingTarget,uMap,cfg,seed+1000*cc);
    end
end

nPerChain=size(chainCells{1},1);
chainsTheta=zeros(nPerChain,4,cfg.numberOfChains);
for cc=1:cfg.numberOfChains
    for rr=1:nPerChain
        chainsTheta(rr,:,cc)=decode_parameters(chainCells{cc}(rr,:),false,fixedDetection);
    end
end

rhat=gelman_rhat(chainsTheta);
ess=effective_sample_size(chainsTheta);
samples=reshape(permute(chainsTheta,[1 3 2]),[],4);

fit.parameterNames={'R_eff','k','p_p','p_a'};
fit.prior=prior;
fit.map=thetaMap;
fit.logPosteriorAtMap=logPostMap;
fit.samples=samples;
fit.posteriorMedian=column_quantiles(samples,0.50);
fit.ciLow=column_quantiles(samples,0.025);
fit.ciHigh=column_quantiles(samples,0.975);
fit.posteriorMean=mean(samples,1);
fit.rhat=rhat;
fit.ess=ess;
fit.acceptance=acceptance;
fit.diagnosticsPass=all(rhat<=cfg.maximumRhat) && ...
    all(ess>=cfg.minimumESS) && all(acceptance>0.10 & acceptance<0.60);
fit.chainSamples=chainsTheta;
fit.inferenceMethod='Dedicated four-chain adaptive MCMC';
fit.importanceESS=NaN;
fit.importanceMaximumWeight=NaN;
fprintf('  posterior fit elapsed time = %.1f minutes.\n',toc(fitTimer)/60);
end

function [kept,acceptance]=run_adaptive_chain(target,uCenter,cfg,seed)
rng(seed,'twister');
d=numel(uCenter);
% Independent overdispersed starts help reveal separated posterior regions;
% the convergence diagnostic is not allowed to rely on four nearly
% identical starts at the same optimum.
u=uCenter(:)'+0.65*randn(1,d);
lp=target(u);
tries=0;
while ~isfinite(lp) && tries<100
    u=uCenter(:)'+0.25*randn(1,d);
    lp=target(u); tries=tries+1;
end
if ~isfinite(lp); error('Could not initialize an MCMC chain.'); end

initialSD=[0.18 0.20 0.25 0.25];
proposalCov=diag(initialSD(1:d).^2);
history=zeros(cfg.mcmcIterations,d);
accepted=0; acceptedWindow=0;

for iter=1:cfg.mcmcIterations
    [L,pFlag]=chol(proposalCov,'lower');
    if pFlag~=0
        proposalCov=proposalCov+1e-8*eye(d);
        L=chol(proposalCov,'lower');
    end
    proposal=u+(L*randn(d,1))';
    lpProposal=target(proposal);
    if isfinite(lpProposal) && log(rand)<lpProposal-lp
        u=proposal; lp=lpProposal; accepted=accepted+1; acceptedWindow=acceptedWindow+1;
    end
    history(iter,:)=u;

    if iter<=cfg.mcmcBurnIn && iter>=500 && ...
            mod(iter,cfg.mcmcAdaptEvery)==0
        empirical=cov(history(1:iter,:));
        windowRate=acceptedWindow/cfg.mcmcAdaptEvery;
        scalar=exp(0.8*(windowRate-0.234));
        proposalCov=scalar*(2.38^2/d)*(empirical+1e-6*eye(d));
        acceptedWindow=0;
    elseif iter<=cfg.mcmcBurnIn && mod(iter,cfg.mcmcAdaptEvery)==0
        acceptedWindow=0;
    end
end

indices=(cfg.mcmcBurnIn+1):cfg.mcmcThin:cfg.mcmcIterations;
kept=history(indices,:);
acceptance=accepted/cfg.mcmcIterations;
end

function lp=logposterior_mixed_measure(u,data,prior,cfg,knownDetection,fixedDetection)
theta=decode_parameters(u,knownDetection,fixedDetection);
if theta(2)<exp(cfg.logKLower) || theta(2)>exp(cfg.logKUpper)
    lp=-Inf; return;
end
logL=grouped_loglik(theta,data,cfg);
if ~isfinite(logL); lp=-Inf; return; end
if knownDetection
    logPrior=0;
else
    pp=theta(3); pa=theta(4);
    logPrior=(prior.ppAlpha-1)*log(pp)+(prior.ppBeta-1)*log1p(-pp) ...
        -betaln(prior.ppAlpha,prior.ppBeta) ...
        +(prior.paAlpha-1)*log(pa)+(prior.paBeta-1)*log1p(-pa) ...
        -betaln(prior.paAlpha,prior.paBeta);
end
% This density is defined with respect to dR d(log k) dp_p dp_a.
lp=logL+logPrior;
end

function lp=logposterior_u_for_sampling(u,data,prior,cfg,knownDetection,fixedDetection)
lp=logposterior_mixed_measure(u,data,prior,cfg,knownDetection,fixedDetection);
if ~isfinite(lp); return; end
theta=decode_parameters(u,knownDetection,fixedDetection);
% Jacobians for logit R, logit p_p and logit p_a. No k Jacobian is
% included because the stated prior measure is uniform in log(k).
lp=lp+log(theta(1))+log1p(-theta(1));
if ~knownDetection
    lp=lp+log(theta(3))+log1p(-theta(3)) ...
        +log(theta(4))+log1p(-theta(4));
end
end

function theta=decode_parameters(u,knownDetection,fixedDetection)
u=u(:)';
R=logistic(u(1)); k=exp(u(2));
if knownDetection
    pp=fixedDetection(1); pa=fixedDetection(2);
else
    pp=logistic(u(3)); pa=logistic(u(4));
end
theta=[R k pp pa];
end

function x=logistic(u)
if u>=0; x=1/(1+exp(-u)); else; e=exp(u); x=e/(1+e); end
x=min(max(x,1e-10),1-1e-10);
end

function u=logit(x)
x=min(max(x,1e-10),1-1e-10);
u=log(x/(1-x));
end

function rhat=gelman_rhat(chains)
% Rank-normalized split R-hat, including the folded version. This is more
% sensitive to nonstationarity and scale differences than the older
% unsplit variance-ratio diagnostic.
[n,p,m]=size(chains);
nHalf=floor(n/2);
split=cat(3,chains(1:nHalf,:,:),chains(n-nHalf+1:n,:,:));
rhat=zeros(1,p);
for j=1:p
    one=squeeze(split(:,j,:));
    z=rank_normalize(one);
    rLocation=basic_rhat(z);

    center=median(one(:));
    folded=abs(one-center);
    zFolded=rank_normalize(folded);
    rScale=basic_rhat(zFolded);
    rhat(j)=max(rLocation,rScale);
end
end

function z=rank_normalize(x)
shape=size(x); values=x(:); N=numel(values);
[~,order]=sort(values);
ranks=zeros(N,1); ranks(order)=(1:N)';
probability=(ranks-3/8)/(N+1/4);
z=sqrt(2)*erfinv(2*probability-1);
z=reshape(z,shape);
end

function value=basic_rhat(one)
[n,m]=size(one);
chainMeans=mean(one,1);
B=n*var(chainMeans,0);
W=mean(var(one,0,1));
varHat=((n-1)/n)*W+B/n;
if W<=0; value=Inf; else; value=sqrt(varHat/W); end
end

function ess=effective_sample_size(chains)
[n,p,m]=size(chains); ess=zeros(1,p);
maxLag=min(1000,n-1);
for j=1:p
    rho=zeros(maxLag,1);
    for lag=1:maxLag
        values=zeros(m,1);
        for cc=1:m
            x=chains(:,j,cc); x=x-mean(x);
            denominator=sum(x.^2);
            if denominator>0
                values(cc)=sum(x(1:n-lag).*x(1+lag:n))/denominator;
            end
        end
        rho(lag)=mean(values);
        if lag>1 && rho(lag)+rho(lag-1)<0
            rho=rho(1:lag-2); break;
        end
    end
    ess(j)=min(m*n,m*n/(1+2*sum(rho)));
end
end

function q=column_quantiles(x,probability)
q=zeros(1,size(x,2));
for jj=1:size(x,2)
    q(jj)=simple_quantile(x(:,jj),probability);
end
end

function q=simple_quantile(x,p)
x=sort(x(isfinite(x)));
if isempty(x); q=NaN; return; end
position=1+(numel(x)-1)*p;
lo=floor(position); hi=ceil(position);
if lo==hi; q=x(lo); else; q=x(lo)+(position-lo)*(x(hi)-x(lo)); end
end

function final_probability_check(theta,data,cfg)
Omax=max(data(:,3)-data(:,4)-2*data(:,5));
[qFast,dFast]=observed_cluster_pmf(theta,Omax,0);
[qDirect,dDirect]=observed_cluster_pmf_direct( ...
    theta,Omax,cfg.SmaxDirectCheck);
llFast=grouped_loglik_from_q(data,qFast);
llDirect=grouped_loglik_from_q(data,qDirect);
if max(abs(qFast-qDirect))>cfg.finalLikelihoodTolerance || ...
        abs(llFast-llDirect)>cfg.finalLikelihoodTolerance
    error(['The exact-series likelihood and the independent direct sum do ' ...
        'not agree. Increase cfg.SmaxDirectCheck and rerun.']);
end
if dFast.omittedObservedMass*sum(data(:,2))>cfg.maximumExpectedOmittedChains
    % This diagnostic uses the displayed Omax. Larger observed sizes are
    % possible but do not enter the grouped likelihood directly.
    fprintf(['Note: probability above the largest required grouped size is %.3g; ' ...
        'this is retained in the conditioning denominator.\n'],dFast.omittedObservedMass);
end
if 1-dDirect.trueMass>1e-8
    fprintf('Direct-check true-size tail mass = %.3g.\n',1-dDirect.trueMass);
end
end

function display_fit_summary(label,fit)
fprintf('\n%s\n',upper(label));
for jj=1:4
    fprintf('%s: MAP=%.6g; median=%.6g; 95%% CrI [%.6g, %.6g]; R-hat=%.4f; ESS=%.0f\n', ...
        fit.parameterNames{jj},fit.map(jj),fit.posteriorMedian(jj), ...
        fit.ciLow(jj),fit.ciHigh(jj),fit.rhat(jj),fit.ess(jj));
end
fprintf('Mean chain acceptance = %.3f. Diagnostics pass = %s.\n', ...
    mean(fit.acceptance),yes_no(fit.diagnosticsPass));
end

function tableOut=make_sensitivity_table(fits,means)
rows=[];
methods={}; importanceESS=[]; maximumImportanceWeight=[];
for iPA=1:numel(means)
    for iPP=1:numel(means)
        f=fits{iPA,iPP};
        rows=[rows; means(iPP),means(iPA),f.map,f.posteriorMedian, ...
            f.ciLow,f.ciHigh,max(f.rhat),min(f.ess),f.diagnosticsPass]; %#ok<AGROW>
        methods{end+1,1}=f.inferenceMethod; %#ok<AGROW>
        importanceESS(end+1,1)=f.importanceESS; %#ok<AGROW>
        maximumImportanceWeight(end+1,1)=f.importanceMaximumWeight; %#ok<AGROW>
    end
end
names={'PPPriorMean','PAPriorMean', ...
    'MAP_R','MAP_k','MAP_pp','MAP_pa', ...
    'Median_R','Median_k','Median_pp','Median_pa', ...
    'Low_R','Low_k','Low_pp','Low_pa', ...
    'High_R','High_k','High_pp','High_pa', ...
    'MaximumRhat','MinimumESS','DiagnosticsPass'};
tableOut=array2table(rows,'VariableNames',names);
tableOut.InferenceMethod=methods;
tableOut.ImportanceESS=importanceESS;
tableOut.MaximumImportanceWeight=maximumImportanceWeight;
end

function ppc=run_posterior_predictive(fit,data,cfg,seed,useParallel)
rng(seed,'twister');
nAvailable=size(fit.samples,1);
indices=round(linspace(1,nAvailable,cfg.ppcDraws));
replicates=zeros(cfg.ppcDraws,4);
tailMass=zeros(cfg.ppcDraws,1);

if useParallel
    parfor dd=1:cfg.ppcDraws
        localSeed=seed+dd;
        rng(localSeed,'twister');
        theta=fit.samples(indices(dd),:);
        [q,tailMass(dd)]=ppc_sampling_q(theta,cfg);
        simulated=simulate_grouped_data(data,q);
        replicates(dd,:)=summarize_grouped_dataset(simulated);
    end
else
    for dd=1:cfg.ppcDraws
        rng(seed+dd,'twister');
        theta=fit.samples(indices(dd),:);
        [q,tailMass(dd)]=ppc_sampling_q(theta,cfg);
        simulated=simulate_grouped_data(data,q);
        replicates(dd,:)=summarize_grouped_dataset(simulated);
    end
end

observed=summarize_grouped_dataset(data);
low=column_quantiles(replicates,0.025);
medianValue=column_quantiles(replicates,0.50);
high=column_quantiles(replicates,0.975);
twoSidedP=zeros(1,4);
for jj=1:4
    lower=mean(replicates(:,jj)<=observed(jj));
    upper=mean(replicates(:,jj)>=observed(jj));
    twoSidedP(jj)=min(1,2*min(lower,upper));
end
ppc.names={'Size1Chains','Size2Chains','Size3PlusChains','TotalCases'};
ppc.observed=observed;
ppc.low=low; ppc.median=medianValue; ppc.high=high;
ppc.twoSidedP=twoSidedP;
ppc.replicates=replicates;
ppc.maximumTruncatedTail=max(tailMass);
ppc.pass=all(observed>=low & observed<=high) && max(tailMass)<1e-4;
end

function [q,tail]=ppc_sampling_q(theta,cfg)
Omax=cfg.OmaxSamplerInitial;
while true
    [q,d]=observed_cluster_pmf(theta,Omax,0);
    tail=d.omittedObservedMass;
    if tail<1e-6
        q=q/sum(q); return;
    end
    if Omax>=cfg.maxSamplerO
        q=q/sum(q); return;
    end
    Omax=min(2*Omax,cfg.maxSamplerO);
end
end

function values=summarize_grouped_dataset(data)
values=[sum(data(:,4)),sum(data(:,5)),sum(data(:,6)),sum(data(:,3))];
end

function T=ppc_to_table(ppc,label)
Analysis=repmat({label},4,1);
Summary=ppc.names(:);
Observed=ppc.observed(:); PredictiveMedian=ppc.median(:);
PredictiveLow=ppc.low(:); PredictiveHigh=ppc.high(:);
TwoSidedPosteriorPredictiveP=ppc.twoSidedP(:);
T=table(Analysis,Summary,Observed,PredictiveMedian,PredictiveLow, ...
    PredictiveHigh,TwoSidedPosteriorPredictiveP);
end

function T=fit_to_table(fit,label)
Analysis=repmat({label},4,1);
Parameter=fit.parameterNames(:);
MAP=fit.map(:); PosteriorMean=fit.posteriorMean(:);
PosteriorMedian=fit.posteriorMedian(:); CredibleLow=fit.ciLow(:);
CredibleHigh=fit.ciHigh(:); Rhat=fit.rhat(:); ESS=fit.ess(:);
T=table(Analysis,Parameter,MAP,PosteriorMean,PosteriorMedian, ...
    CredibleLow,CredibleHigh,Rhat,ESS);
end

function fig=draw_figure8(primaryFit,fullFit)
fig=figure('Color','w','Units','pixels','Position',[50 50 1500 900]);
names={'R_{eff}','k','p_p','p_a'};
for jj=1:4
    ax=subplot(2,2,jj,'Parent',fig); hold(ax,'on');
    fits={primaryFit,fullFit};
    colors=[0.05 0.35 0.90; 0.90 0.25 0.10];
    for aa=1:2
        f=fits{aa};
        med=f.posteriorMedian(jj); lo=f.ciLow(jj); hi=f.ciHigh(jj);
        errorbar(ax,aa,med,med-lo,hi-med,'o','Color',colors(aa,:), ...
            'MarkerFaceColor',colors(aa,:),'LineWidth',2,'MarkerSize',8, ...
            'CapSize',12);
        plot(ax,aa,f.map(jj),'d','Color','k','MarkerFaceColor','w', ...
            'MarkerSize',7,'LineWidth',1.2);
    end
    xlim(ax,[0.5 2.5]); xticks(ax,[1 2]);
    xticklabels(ax,{'2020-2023','Through Mar 2024'});
    title(ax,names{jj},'Interpreter','tex','FontWeight','bold');
    grid(ax,'on'); box(ax,'on'); set(ax,'FontSize',12,'LineWidth',1);
    if jj==1; ylim(ax,[0 1]); end
    if jj>=3; ylim(ax,[0 1]); end
end
annotation(fig,'textbox',[0.10 0.94 0.80 0.045], ...
    'String','U.S. grouped-data posterior summaries', ...
    'HorizontalAlignment','center','EdgeColor','none', ...
    'FontSize',20,'FontWeight','bold');
annotation(fig,'textbox',[0.10 0.90 0.80 0.035], ...
    'String','Circles: posterior medians and 95% credible intervals; diamonds: MAP estimates', ...
    'HorizontalAlignment','center','EdgeColor','none','FontSize',12);
end

function fig=draw_figure9(fits,means)
fig=figure('Color','w','Units','pixels','Position',[50 50 1500 900]);
titles={'R_{eff}','k','p_p','p_a'};
for jj=1:4
    values=zeros(numel(means));
    for iPA=1:numel(means)
        for iPP=1:numel(means)
            values(iPA,iPP)=fits{iPA,iPP}.posteriorMedian(jj);
        end
    end
    ax=subplot(2,2,jj,'Parent',fig);
    imagesc(ax,means,means,values); set(ax,'YDir','normal');
    colorbar(ax); title(ax,titles{jj},'Interpreter','tex','FontWeight','bold');
    xlabel(ax,'p_p prior mean','Interpreter','tex');
    ylabel(ax,'p_a prior mean','Interpreter','tex');
    xticks(ax,means); yticks(ax,means); box(ax,'on');
    set(ax,'FontSize',12,'LineWidth',1);
    for rr=1:numel(means)
        for cc=1:numel(means)
            text(ax,means(cc),means(rr),sprintf('%.3f',values(rr,cc)), ...
                'HorizontalAlignment','center','FontWeight','bold','FontSize',11);
        end
    end
end
annotation(fig,'textbox',[0.10 0.94 0.80 0.045], ...
    'String','U.S. detection-prior sensitivity (2020-2023 grouped data)', ...
    'HorizontalAlignment','center','EdgeColor','none', ...
    'FontSize',19,'FontWeight','bold');
annotation(fig,'textbox',[0.10 0.90 0.80 0.032], ...
    'String','Only prior locations change; the observed grouped data are identical in all nine analyses', ...
    'HorizontalAlignment','center','EdgeColor','none','FontSize',12);
end

function fig=draw_figureS1(summary)
fig=figure('Color','w','Units','pixels','Position',[50 50 1400 800]);
parameters={'R_eff','k','p_p','p_a'};
metrics={'MAE','Approximate95Coverage'};
for mm=1:2
    ax=subplot(1,2,mm,'Parent',fig); hold(ax,'on');
    values=[]; labels={};
    for rr=1:height(summary)
        if any(strcmp(summary.Parameter{rr},parameters))
            values(end+1)=summary.(metrics{mm})(rr); %#ok<AGROW>
            labels{end+1}=sprintf('S%d %s',ceil(rr/4),summary.Parameter{rr}); %#ok<AGROW>
        end
    end
    bar(ax,values,'FaceColor',[0.15 0.45 0.80]);
    xticks(ax,1:numel(values)); xticklabels(ax,labels); xtickangle(ax,60);
    grid(ax,'on'); box(ax,'on'); set(ax,'FontSize',9);
    if mm==1; title(ax,'Mean absolute error');
    else; title(ax,'Approximate 95% interval coverage'); ylim(ax,[0 1]);
        yline(ax,0.95,'r--','LineWidth',1.5);
    end
end
annotation(fig,'textbox',[0.10 0.94 0.80 0.045], ...
    'String','Grouped-data simulation validation', ...
    'HorizontalAlignment','center','EdgeColor','none', ...
    'FontSize',18,'FontWeight','bold');
end

function fig=draw_figureS2(primaryPPC,fullPPC)
fig=figure('Color','w','Units','pixels','Position',[50 50 1400 800]);
ppcs={primaryPPC,fullPPC}; labels={'2020-2023','Through Mar 2024'};
for aa=1:2
    ax=subplot(1,2,aa,'Parent',fig); hold(ax,'on');
    p=ppcs{aa}; x=1:4;
    errorbar(ax,x,p.median,p.median-p.low,p.high-p.median,'o', ...
        'Color',[0.05 0.35 0.90],'MarkerFaceColor',[0.05 0.35 0.90], ...
        'LineWidth',2,'MarkerSize',8,'CapSize',10);
    plot(ax,x,p.observed,'rd','MarkerFaceColor','r','MarkerSize',8);
    xticks(ax,x); xticklabels(ax,{'Size 1','Size 2','Size >=3','Total cases'});
    xtickangle(ax,20); title(ax,labels{aa},'FontWeight','bold');
    ylabel(ax,'Count'); grid(ax,'on'); box(ax,'on'); set(ax,'FontSize',11);
end
annotation(fig,'textbox',[0.10 0.94 0.80 0.045], ...
    'String','Posterior-predictive checks', ...
    'HorizontalAlignment','center','EdgeColor','none', ...
    'FontSize',18,'FontWeight','bold');
annotation(fig,'textbox',[0.10 0.90 0.80 0.032], ...
    'String','Blue: predictive median and 95% interval; red: observed grouped summary', ...
    'HorizontalAlignment','center','EdgeColor','none','FontSize',12);
end

function save_publication_figure(fig,folder,baseName,pixelWidth,pixelHeight)
pngFile=fullfile(folder,[baseName '.png']);
pdfFile=fullfile(folder,[baseName '.pdf']);
figFile=fullfile(folder,[baseName '.fig']);
set(fig,'PaperPositionMode','auto');
savefig(fig,figFile);

% Use print rather than exportgraphics for compatibility with older MATLAB.
dpi=300;
widthInches=pixelWidth/dpi; heightInches=pixelHeight/dpi;
set(fig,'PaperUnits','inches','PaperSize',[widthInches heightInches], ...
    'PaperPosition',[0 0 widthInches heightInches]);
print(fig,pngFile,'-dpng',sprintf('-r%d',dpi));
print(fig,pdfFile,'-dpdf','-painters','-bestfit');

if ~exist(pngFile,'file') || ~exist(pdfFile,'file')
    error('A figure export failed for %s.',baseName);
end
end

function write_readme(folder,cfg,primaryFit,fullFit)
fileName=fullfile(folder,'README_US_Grouped_MAP.txt');
fid=fopen(fileName,'w');
if fid<0; return; end
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'U.S. GROUPED-DATA MEASLES ANALYSIS\n\n');
fprintf(fid,'This analysis never reconstructs unavailable individual chain sizes.\n');
fprintf(fid,'Primary dataset: 2020-2023, 62 observed chains, 241 cases.\n');
fprintf(fid,'Sensitivity dataset: 2020-March 28 2024, 92 chains, 338 cases.\n');
fprintf(fid,'The 2024 snapshot contains one outbreak that was ongoing at cutoff.\n\n');
fprintf(fid,'Likelihood:\n');
fprintf(fid,['The likelihood conditions on the observed number of chains in each year. ' ...
    'It uses the exact counts of size-1, size-2, and size>=3 chains, plus the ' ...
    'total cases contributed by the size>=3 chains.\n']);
fprintf(fid,['No individual sizes inside the size>=3 category are reconstructed. ' ...
    'Their possible allocations are summed by dynamic programming.\n']);
fprintf(fid,['Observed-size probabilities are computed by an exact truncated PGF ' ...
    'power series and independently checked by direct summation.\n\n']);
fprintf(fid,'Primary priors:\n');
fprintf(fid,'R_eff: Uniform(0,1), defining the subcritical model.\n');
fprintf(fid,'log(k): Uniform(log(1e-3),log(1e3)).\n');
fprintf(fid,'p_p: Beta(2,2).\n');
fprintf(fid,'p_a: Beta(2,2).\n\n');
fprintf(fid,'Primary posterior medians [95%% credible intervals]:\n');
for jj=1:4
    fprintf(fid,'%s: %.6g [%.6g, %.6g]\n',primaryFit.parameterNames{jj}, ...
        primaryFit.posteriorMedian(jj),primaryFit.ciLow(jj),primaryFit.ciHigh(jj));
end
fprintf(fid,'\nFull-snapshot posterior medians [95%% credible intervals]:\n');
for jj=1:4
    fprintf(fid,'%s: %.6g [%.6g, %.6g]\n',fullFit.parameterNames{jj}, ...
        fullFit.posteriorMedian(jj),fullFit.ciLow(jj),fullFit.ciHigh(jj));
end
fprintf(fid,'\nWorkflow version: %s\n',cfg.version);
end

function answer=yes_no(tf)
if tf; answer='YES'; else; answer='NO'; end
end

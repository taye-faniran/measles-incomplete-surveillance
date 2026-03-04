function Figure1_generate_known_detection
% Generate results for Figure 1: known passive and active detection

clear; clc;

baseFolder = fileparts(mfilename('fullpath'));
load(fullfile(baseFolder, 'Figure1_known_detection_results.mat'), 'XX');

r0 = 0.5;
k = 0.15;

true_CV = sqrt(r0 + r0^2/k) / r0;
fprintf('True CV = %.6f\n', true_CV);

target_observed_chains = 5000;
max_generations = 30000;

num_mm = 4;
num_ii = 4;

RR0 = NaN(num_mm, num_ii);
kk = NaN(num_mm, num_ii);
RR0_lower = NaN(num_mm, num_ii);
RR0_upper = NaN(num_mm, num_ii);
kk_lower = NaN(num_mm, num_ii);
kk_upper = NaN(num_mm, num_ii);
CV_est = NaN(num_mm, num_ii);

for mm = 1:num_mm
    for ii = 1:num_ii
        fprintf('Running cell (row=%d, col=%d)\n', mm, ii);

        active = 0.25 * mm;
        passive = 0.25 * ii;

        [Observed_Cluster_Dist, ~] = simulate_many_chains( ...
            r0, k, passive, active, target_observed_chains, max_generations);

        Observed_Cluster_Distribution = Observed_Cluster_Dist(:, 2);

        est = fit_observed( ...
            Observed_Cluster_Distribution, ...
            struct('R0', [], 'k', [], 'p_p', passive, 'p_a', active), ...
            default_priors());

        RR0(mm,ii) = est.R0;
        kk(mm,ii) = est.k;
        RR0_lower(mm,ii) = est.R0_CI(1);
        RR0_upper(mm,ii) = est.R0_CI(2);
        kk_lower(mm,ii) = est.k_CI(1);
        kk_upper(mm,ii) = est.k_CI(2);
        CV_est(mm,ii) = sqrt(est.R0 + est.R0^2/est.k) / est.R0;
    end
end

XX(:,:,1) = RR0;
XX(:,:,2) = kk;
XX(:,:,3) = RR0_lower;
XX(:,:,4) = RR0_upper;
XX(:,:,5) = kk_lower;
XX(:,:,6) = kk_upper;
XX(:,:,7) = CV_est;

disp('Saved pages in XX:');
disp('1 = R0 estimate');
disp('2 = k estimate');
disp('3 = R0 lower 95% CI');
disp('4 = R0 upper 95% CI');
disp('5 = k lower 95% CI');
disp('6 = k upper 95% CI');
disp('7 = estimated CV');

save(folderPath, 'XX');
fprintf('Saved results to:\n%s\n', folderPath);

end


function est = fit_observed(counts, fix, pri)

counts = counts(:)';
Xmax = numel(counts);

to_est.R0  = ~isfield(fix,'R0')  || isempty(fix.R0);
to_est.k   = ~isfield(fix,'k')   || isempty(fix.k);
to_est.p_p = ~isfield(fix,'p_p') || isempty(fix.p_p);
to_est.p_a = ~isfield(fix,'p_a') || isempty(fix.p_a);

if nargin < 3 || isempty(pri)
    pri = default_priors();
end

R0_0 = pick_default(fix,'R0',0.6);
k_0  = pick_default(fix,'k',2.0);
p1_0 = pick_default(fix,'p_p',0.5);
p2_0 = pick_default(fix,'p_a',0.5);

theta0 = [];
lb = [];
ub = [];

if to_est.R0
    theta0(end+1) = log(R0_0);
    lb(end+1) = -inf;
    ub(end+1) = 0;
end
if to_est.k
    theta0(end+1) = log(k_0);
    lb(end+1) = log(1e-6);
    ub(end+1) = inf;
end
if to_est.p_p
    theta0(end+1) = logit(p1_0);
    lb(end+1) = -inf;
    ub(end+1) = inf;
end
if to_est.p_a
    theta0(end+1) = logit(p2_0);
    lb(end+1) = -inf;
    ub(end+1) = inf;
end

obj = @(th) obj_single(th, counts, fix, to_est, Xmax, pri);
opts = optimoptions('fmincon','Display','off','Algorithm','interior-point');
[th_hat, fval, ~, ~, ~, ~, hessian] = fmincon(obj, theta0, [], [], [], [], lb, ub, [], opts);

[R0,k,p_p,p_a] = unpack_params(th_hat, fix, to_est);

ci.R0 = [NaN, NaN];
ci.k = [NaN, NaN];

if ~isempty(hessian)
    try
        cov_th = inv(hessian);
        if any(diag(cov_th) < 0)
            error('Covariance matrix is not positive definite.');
        end

        se_th = sqrt(diag(cov_th));
        z = 1.96;
        th_lower = th_hat' - z * se_th;
        th_upper = th_hat' + z * se_th;

        idx = 1;
        if to_est.R0
            ci.R0 = exp([th_lower(idx), th_upper(idx)]);
            idx = idx + 1;
        end
        if to_est.k
            ci.k = exp([th_lower(idx), th_upper(idx)]);
        end
    catch ME
        warning('Could not calculate confidence intervals: %s', ME.message);
    end
end

est = struct('R0',R0,'k',k,'p_p',p_p,'p_a',p_a,'obj',fval,'Xmax',Xmax, ...
             'R0_CI',ci.R0,'k_CI',ci.k);
end


function val = obj_single(th, counts, fix, to_est, Xmax, pri)
[R0,k,p_p,p_a] = unpack_params(th, fix, to_est);

if ~(R0>=0 && k>1e-6 && p_p>=0 && p_p<=1 && p_a>=0 && p_a<=1)
    val = 1e12;
    return;
end

px = observed_size_pmf_condpos(R0,k,p_p,p_a,Xmax);
px(px<=0) = realmin;

nll = -sum(counts .* log(px));
val = nll + prior_penalty(p_p,p_a,pri);
end


function pen = prior_penalty(p_p,p_a,pri)
nlogprior_p1 = -((pri.p1.a-1)*log(max(p_p,eps)) + (pri.p1.b-1)*log(max(1-p_p,eps))) + betaln(pri.p1.a,pri.p1.b);
nlogprior_p2 = -((pri.p2.a-1)*log(max(p_a,eps)) + (pri.p2.b-1)*log(max(1-p_a,eps))) + betaln(pri.p2.a,pri.p2.b);
pen = pri.strength * (nlogprior_p1 + nlogprior_p2);
end


function pri = default_priors()
pri.p1.a = 2; 
pri.p1.b = 2;
pri.p2.a = 2; 
pri.p2.b = 2;
pri.strength = 1.0;
end


function val = pick_default(fix, field, def)
if isfield(fix,field) && ~isempty(fix.(field))
    val = fix.(field);
else
    val = def;
end
end

function y = logit(p)
y = log(p) - log(1-p);
end

function p = logistic(y)
p = 1./(1+exp(-y));
end

function [R0,k,p_p,p_a] = unpack_params(th, fix, to_est)
idx = 1;
if to_est.R0,  R0  = exp(th(idx)); idx=idx+1; else, R0  = fix.R0;  end
if to_est.k,   k   = exp(th(idx)); idx=idx+1; else, k   = fix.k;   end
if to_est.p_p, p_p = logistic(th(idx)); idx=idx+1; else, p_p = fix.p_p; end
if to_est.p_a, p_a = logistic(th(idx));           else, p_a = fix.p_a; end
end


function px = observed_size_pmf_condpos(R0,k,p_p,p_a,Xmax)
n = 1;
Jmax = max(2000,10*Xmax);
pJ = true_cluster_pmf(R0,k,n,Jmax);

massX = zeros(1,Xmax+1);
logp1 = log(max(p_p,eps));
logq1 = log(max(1-p_p,eps));
logp2 = log(max(p_a,eps));
logq2 = log(max(1-p_a,eps));

for j = n:Jmax
    pj = pJ(j);
    if pj < 1e-15
        continue;
    end

    massX(1) = massX(1) + pj * exp(j*logq1);

    xmax = min(j,Xmax);
    for x = 1:xmax
        pvec = 1:min(j,x);
        logC1 = gammaln(j+1)-gammaln(pvec+1)-gammaln(j-pvec+1);
        logC2 = gammaln(j-pvec+1)-gammaln(x-pvec+1)-gammaln(j-x+1);
        term = exp(logC1 + pvec*logp1 + (j-pvec)*logq1 + ...
                   logC2 + (x-pvec)*logp2 + (j-x)*logq2);
        massX(x+1) = massX(x+1) + pj * sum(term);
    end
end

p_pos = 1 - massX(1);
if p_pos <= 0
    p_pos = realmin;
end

px = massX(2:end) / p_pos;
px = px / sum(px);
end


function pJ = true_cluster_pmf(R0,k,n,Jmax)
pJ = zeros(1,Jmax);

if R0 == 0
    pJ(n) = 1;
    pJ = pJ / sum(pJ);
    return;
end

jj = n:Jmax;
logRk = log(R0/k);
log1Rk = log(1+R0/k);

val = log(n) - log(jj) + gammaln(k*jj + jj - n) - gammaln(k*jj) - gammaln(jj - n + 1) ...
    + (jj - n).*logRk - (k*jj + jj - n).*log1Rk;

pJ(jj) = exp(val);
s = sum(pJ);
if s > 0
    pJ = pJ / s;
end
end


function [observed_cluster_summary, zero_observed_cases_counter] = simulate_many_chains(reproduction_number, dispersion_parameter, p_passive, p_active, target_observed_chains, max_generations)

observed_outbreak_sizes = {};
zero_observed_cases_counter = 0;

while length(observed_outbreak_sizes) < target_observed_chains
    current_result = simulate_single_chain_with_observation( ...
        reproduction_number, dispersion_parameter, max_generations, p_passive, p_active);

    if current_result > 0
        observed_outbreak_sizes{end+1} = current_result;
    else
        zero_observed_cases_counter = zero_observed_cases_counter + 1;
    end
end

observed_sizes_vector = cell2mat(observed_outbreak_sizes);

if ~isempty(observed_sizes_vector)
    observed_cluster_summary = tabulate(observed_sizes_vector);
    observed_cluster_summary = observed_cluster_summary(:,1:2);
else
    observed_cluster_summary = [];
end
end


function observed_cases_count = simulate_single_chain_with_observation(R, k, max_gens, p_passive, p_active)

case_id = 1;
generation = 1;
cases_table = table(case_id, generation);

current_generation = 1;
last_generation_with_cases = 1;

while current_generation <= last_generation_with_cases && current_generation < max_gens
    num_parents = sum(cases_table.generation == current_generation);

    if num_parents > 0
        offspring_counts = poissrnd(gamrnd(k, R/k, num_parents, 1));
        total_new_cases = sum(offspring_counts);

        if total_new_cases > 0
            new_case_ids = (height(cases_table) + 1):(height(cases_table) + total_new_cases);
            new_generations = repmat(current_generation + 1, total_new_cases, 1);
            new_cases_table = table(new_case_ids', new_generations, ...
                'VariableNames', {'case_id', 'generation'});
            cases_table = [cases_table; new_cases_table];
            last_generation_with_cases = current_generation + 1;
        end
    end

    current_generation = current_generation + 1;
end

num_total_cases = height(cases_table);
is_passively_observed = rand(num_total_cases,1) < p_passive;
is_actively_observed = false(num_total_cases,1);

if any(is_passively_observed)
    missed_indices = find(~is_passively_observed);
    num_missed = length(missed_indices);
    active_detections = rand(num_missed,1) < p_active;
    is_actively_observed(missed_indices) = active_detections;
end

observed_cases_count = sum(is_passively_observed | is_actively_observed);
end
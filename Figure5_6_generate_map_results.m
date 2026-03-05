function Figure5_6_generate_map_results
% Generate MAP simulation results used for:
% Figure 5 = MAP transmission recovery (R_eff and k)
% Figure 6 = MAP detection recovery (p_p and p_a)
%
% This script simulates outbreaks under true R_eff = 0.5 and k = 0.15,
% applies two-stage detection, then estimates R_eff, k, p_p, and p_a
% using MAP with informative priors on p_p and p_a.
%
% Output:
%   Figure5_6_map_results.mat
% saved automatically in the same folder as this script.

clear; clc;

baseFolder = fileparts(mfilename('fullpath'));
outFile    = fullfile(baseFolder, 'Figure5_6_map_results.mat');

true_params.R0 = 0.5;
true_params.k  = 0.15;

target_observed_chains = 5000;
max_generations = 20000;

true_CV = sqrt(true_params.R0 + true_params.R0^2 / true_params.k) / true_params.R0;
fprintf('True Coefficient of Variation (CV): %.6f\n', true_CV);

% Preallocate result arrays
RR0      = NaN(4,4);
RR0_CI_1 = NaN(4,4);
RR0_CI_2 = NaN(4,4);

kk       = NaN(4,4);
kk_CI_1  = NaN(4,4);
kk_CI_2  = NaN(4,4);

pp       = NaN(4,4);
pp_CI_1  = NaN(4,4);
pp_CI_2  = NaN(4,4);

pa       = NaN(4,4);
pa_CI_1  = NaN(4,4);
pa_CI_2  = NaN(4,4);

CV_est   = NaN(4,4);

for mm = 1:4
    for ii = 1:4
        fprintf('Running grid cell: row = %d, col = %d\n', mm, ii);

        true_params.p1 = 0.25 * ii;   % passive detection
        true_params.p2 = 0.25 * mm;   % active detection

        fprintf('   Passive detection p_p = %.2f\n', true_params.p1);
        fprintf('   Active detection  p_a = %.2f\n', true_params.p2);

        [~, observed_counts] = Simulate_True_vs_Observed( ...
            true_params.R0, true_params.k, ...
            true_params.p1, true_params.p2, ...
            target_observed_chains, max_generations);

        % Build intentionally biased priors for MAP
        mean_p1 = 0.75 * true_params.p1;
        ci_p1   = [mean_p1 - 0.1, mean_p1 + 0.25];
        objective_for_alpha = @(alpha_guess) ...
            (betainv(0.975, alpha_guess, alpha_guess * (1/mean_p1 - 1)) - ci_p1(2))^2;
        options_prior = optimset('TolX', 1e-6, 'Display', 'off');
        p1_alpha = fminsearch(objective_for_alpha, 50, options_prior);
        p1_beta  = p1_alpha * (1/mean_p1 - 1);

        pri.p1.a = p1_alpha;
        pri.p1.b = p1_beta;

        act = true_params.p2 / 1.15;
        pri.p2.a = (1 - 1.09 * act) / 0.09;
        pri.p2.b = (1 - act) * (1 - 1.09 * act) / (0.09 * act);
        pri.strength = 1.0;

        fix = struct();

        tic;

        R0_starts = [0.2, 0.5, 0.8];
        k_starts  = [0.2, 0.5, 1.0];
        p1_starts = [0.2, 0.4, 0.7];
        p2_starts = [0.2, 0.4, 0.7];

        [R0_grid, k_grid, p1_grid, p2_grid] = ndgrid(R0_starts, k_starts, p1_starts, p2_starts);
        start_points = [R0_grid(:), k_grid(:), p1_grid(:), p2_grid(:)];
        num_starts = size(start_points, 1);

        all_results = cell(num_starts, 1);

        parfor i = 1:num_starts
            initial_guess = struct( ...
                'R0',  start_points(i,1), ...
                'k',   start_points(i,2), ...
                'p_p', start_points(i,3), ...
                'p_a', start_points(i,4));

            all_results{i} = fit_observed_no_k_prior( ...
                observed_counts, fix, pri, initial_guess, 'off');
        end

        obj_vals = cellfun(@(s) s.obj, all_results);
        [~, best_idx] = min(obj_vals);
        est = all_results{best_idx};

        toc;

        RR0(mm,ii)      = est.R0;
        RR0_CI_1(mm,ii) = est.CI.R0(1);
        RR0_CI_2(mm,ii) = est.CI.R0(2);

        kk(mm,ii)       = est.k;
        kk_CI_1(mm,ii)  = est.CI.k(1);
        kk_CI_2(mm,ii)  = est.CI.k(2);

        pp(mm,ii)       = est.p_p;
        pp_CI_1(mm,ii)  = est.CI.p_p(1);
        pp_CI_2(mm,ii)  = est.CI.p_p(2);

        pa(mm,ii)       = est.p_a;
        pa_CI_1(mm,ii)  = est.CI.p_a(1);
        pa_CI_2(mm,ii)  = est.CI.p_a(2);

        CV_est(mm,ii) = sqrt(est.R0 + est.R0^2 / est.k) / est.R0;

        fprintf('   Estimated R_eff = %.4f\n', est.R0);
        fprintf('   Estimated k     = %.4f\n', est.k);
        fprintf('   Estimated p_p   = %.4f\n', est.p_p);
        fprintf('   Estimated p_a   = %.4f\n', est.p_a);
        fprintf('   Estimated CV    = %.4f\n\n', CV_est(mm,ii));
    end
end

% Pages 1-6: Figure 5
XX(:,:,1)  = RR0;
XX(:,:,2)  = RR0_CI_1;
XX(:,:,3)  = RR0_CI_2;
XX(:,:,4)  = kk;
XX(:,:,5)  = kk_CI_1;
XX(:,:,6)  = kk_CI_2;

% Pages 7-12: Figure 6
XX(:,:,7)  = pp;
XX(:,:,8)  = pp_CI_1;
XX(:,:,9)  = pp_CI_2;
XX(:,:,10) = pa;
XX(:,:,11) = pa_CI_1;
XX(:,:,12) = pa_CI_2;

% Page 13: CV diagnostic
XX(:,:,13) = CV_est;

save(outFile, 'XX');

fprintf('Saved results to:\n%s\n', outFile);

end

function est = fit_observed_no_k_prior(counts, fix, pri, initial_guess, display_option)

if nargin < 5
    display_option = 'iter';
end

counts = counts(:)';
Xmax = numel(counts);

to_est.R0  = ~isfield(fix,'R0')  || isempty(fix.R0);
to_est.k   = ~isfield(fix,'k')   || isempty(fix.k);
to_est.p_p = ~isfield(fix,'p_p') || isempty(fix.p_p);
to_est.p_a = ~isfield(fix,'p_a') || isempty(fix.p_a);

R0_0 = pick_default(fix,'R0',  initial_guess.R0);
k_0  = pick_default(fix,'k',   initial_guess.k);
p1_0 = pick_default(fix,'p_p', initial_guess.p_p);
p2_0 = pick_default(fix,'p_a', initial_guess.p_a);

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
    lb(end+1) = -inf;
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

opts = optimoptions('fmincon', ...
    'Display', display_option, ...
    'Algorithm', 'interior-point', ...
    'MaxFunctionEvaluations', 3000);

obj = @(th) obj_single(th, counts, fix, to_est, Xmax, pri);
[th_hat, fval, ~, ~, ~, ~, hessian] = fmincon(obj, theta0, [], [], [], [], lb, ub, [], opts);

[R0,k,p_p,p_a] = unpack_params(th_hat, fix, to_est);
est = struct('R0',R0,'k',k,'p_p',p_p,'p_a',p_a,'obj',fval,'Xmax',Xmax);

try
    var_cov_matrix_transformed = inv(hessian);
    diag_variances = diag(var_cov_matrix_transformed);

    if any(diag_variances < 0)
        warning('Negative variance detected from Hessian. Using absolute values. CIs may be unreliable.');
        diag_variances = abs(diag_variances);
    end

    std_errors_transformed = sqrt(diag_variances);

    idx = 1;
    if to_est.R0
        ci_log = th_hat(idx) + [-1.96, 1.96] * std_errors_transformed(idx);
        est.CI.R0 = exp(ci_log);
        idx = idx + 1;
    end
    if to_est.k
        ci_log = th_hat(idx) + [-1.96, 1.96] * std_errors_transformed(idx);
        est.CI.k = exp(ci_log);
        idx = idx + 1;
    end
    if to_est.p_p
        ci_logit = th_hat(idx) + [-1.96, 1.96] * std_errors_transformed(idx);
        est.CI.p_p = logistic(ci_logit);
        idx = idx + 1;
    end
    if to_est.p_a
        ci_logit = th_hat(idx) + [-1.96, 1.96] * std_errors_transformed(idx);
        est.CI.p_a = logistic(ci_logit);
    end

catch ME
    warning('Could not calculate CIs. Error: %s', ME.message);
    if to_est.R0,  est.CI.R0  = [NaN, NaN]; end
    if to_est.k,   est.CI.k   = [NaN, NaN]; end
    if to_est.p_p, est.CI.p_p = [NaN, NaN]; end
    if to_est.p_a, est.CI.p_a = [NaN, NaN]; end
end

end

function val = obj_single(th, counts, fix, to_est, Xmax, pri)

[R0,k,p_p,p_a] = unpack_params(th, fix, to_est);

if ~(R0 >= 0 && k > 1e-6 && p_p > 0 && p_p < 1 && p_a > 0 && p_a < 1)
    val = 1e12;
    return
end

px = observed_size_pmf_condpos_FAST(R0,k,p_p,p_a,Xmax);
px(px <= 0) = realmin;

nll = -sum(counts .* log(px));
val = nll + prior_penalty(p_p,p_a,pri);

end

function pen = prior_penalty(p_p,p_a,pri)

nlogprior_p1 = -((pri.p1.a-1)*log(max(p_p,eps)) + (pri.p1.b-1)*log(max(1-p_p,eps))) ...
               + betaln(pri.p1.a,pri.p1.b);
nlogprior_p2 = -((pri.p2.a-1)*log(max(p_a,eps)) + (pri.p2.b-1)*log(max(1-p_a,eps))) ...
               + betaln(pri.p2.a,pri.p2.b);

pen = pri.strength * (nlogprior_p1 + nlogprior_p2);

end

function px = observed_size_pmf_condpos_FAST(R0,k,p_p,p_a,Xmax)

n = 1;
Jmax = max(2000, 10 * Xmax);
pJ = true_cluster_pmf(R0,k,n,Jmax);

p_obs = p_p + (1-p_p)*p_a;

x_vec = (1:Xmax)';
j_vec = 1:Jmax;

j_minus_x_plus_1 = j_vec - x_vec + 1;
invalid_mask = j_minus_x_plus_1 <= 0;
j_minus_x_plus_1(invalid_mask) = 1;

gammaln_j_minus_x = gammaln(j_minus_x_plus_1);
gammaln_j_minus_x(invalid_mask) = Inf;

log_binom_coeff = gammaln(j_vec + 1) - gammaln(x_vec + 1) - gammaln_j_minus_x;

term1_matrix = exp(log_binom_coeff + x_vec * log(p_obs) + (j_vec - x_vec) * log(1 - p_obs));
term2_matrix = exp(j_vec * log(1 - p_p) + log_binom_coeff + x_vec * log(p_a) + (j_vec - x_vec) * log(1 - p_a));

P_O_given_S = term1_matrix - term2_matrix;
P_O_given_S(P_O_given_S < 0) = 0;

massX_positive = P_O_given_S * pJ';
p_O_zero = sum((1-p_p).^j_vec .* pJ);

p_pos = 1 - p_O_zero;
if p_pos <= 0
    p_pos = realmin;
end

px = massX_positive' / p_pos;
if sum(px) > 0
    px = px / sum(px);
end

end

function pJ = true_cluster_pmf(R0,k,n,Jmax)

pJ = zeros(1,Jmax);

if R0 == 0
    pJ(n) = 1;
    pJ = pJ / sum(pJ);
    return
end

jj = n:Jmax;
valid_mask = (k*jj + jj - n) > 0 & (jj - n + 1) > 0;
jj_valid = jj(valid_mask);

if isempty(jj_valid)
    return
end

logRk = log(R0/k);
log1Rk = log(1 + R0/k);

val = log(n) - log(jj_valid) ...
    + gammaln(k*jj_valid + jj_valid - n) ...
    - gammaln(k*jj_valid) ...
    - gammaln(jj_valid - n + 1) ...
    + (jj_valid - n).*logRk ...
    - (k*jj_valid + jj_valid - n).*log1Rk;

pJ(jj_valid) = exp(val);

s = sum(pJ);
if s > 0
    pJ = pJ / s;
end

end

function [true_counts_all, observed_counts] = Simulate_True_vs_Observed( ...
    reproduction_number, dispersion_parameter, p_passive, p_active, ...
    target_observed_chains, max_generations)

all_true_outbreak_sizes = [];
observed_final_sizes = [];

while length(observed_final_sizes) < target_observed_chains
    true_size = simulate_single_true_chain(reproduction_number, dispersion_parameter, max_generations);
    all_true_outbreak_sizes(end+1) = true_size;

    detected_p1 = binornd(true_size, p_passive);

    if detected_p1 > 0
        missed_cases = true_size - detected_p1;
        detected_p2 = binornd(missed_cases, p_active);
        observed_size = detected_p1 + detected_p2;
        observed_final_sizes(end+1) = observed_size;
    end
end

if ~isempty(all_true_outbreak_sizes)
    max_true_size = max(all_true_outbreak_sizes);
    true_counts_all = histcounts(all_true_outbreak_sizes, 1:max_true_size+1)';
else
    true_counts_all = [];
end

if ~isempty(observed_final_sizes)
    max_obs_size = max(observed_final_sizes);
    observed_counts = histcounts(observed_final_sizes, 1:max_obs_size+1)';
else
    observed_counts = [];
end

end

function true_size = simulate_single_true_chain(R, k, max_gens)

popsize = 1;
true_size = 1;
gen = 1;

while popsize > 0 && gen < max_gens
    offspring_counts = poissrnd(gamrnd(k, R/k, popsize, 1));
    popsize = sum(offspring_counts);
    true_size = true_size + popsize;
    gen = gen + 1;
end

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
p = 1 ./ (1 + exp(-y));
end

function [R0,k,p_p,p_a] = unpack_params(th, fix, to_est)
idx = 1;
if to_est.R0
    R0 = exp(th(idx)); idx = idx + 1;
else
    R0 = fix.R0;
end
if to_est.k
    k = exp(th(idx)); idx = idx + 1;
else
    k = fix.k;
end
if to_est.p_p
    p_p = logistic(th(idx)); idx = idx + 1;
else
    p_p = fix.p_p;
end
if to_est.p_a
    p_a = logistic(th(idx));
else
    p_a = fix.p_a;
end
end
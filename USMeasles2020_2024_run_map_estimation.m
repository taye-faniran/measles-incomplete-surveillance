function USMeasles2020_2024_run_map_estimation
% Run MAP estimation for U.S. measles cluster data, 2020-2024.
%
% This script estimates:
%   - R_eff : effective reproduction number
%   - k     : dispersion parameter
%   - p_p   : passive detection probability
%   - p_a   : active detection probability
%
% The observed data are grouped as:
%   [observed cluster size, number of index cases, frequency]
%
% The script uses:
%   - a cluster-size likelihood conditioned on O > 0
%   - informative Beta priors for p_p and p_a
%   - multi-start optimization with fmincon
%
% Output:
%   Prints parameter estimates, 95% confidence intervals, and estimated CV.

clear; clc;

%% Observed U.S. measles cluster data (2020-2024)
% Columns:
%   1 = observed cluster size
%   2 = number of index cases / introductions
%   3 = frequency
observedClusterData = [ ...
     13,  9, 1
     49,  7, 1
    121, 15, 1
     58, 31, 1
     97, 30, 1];

%% Informative priors for detection parameters
% Passive detection prior
pri.p1.a = 1.843;
pri.p1.b = 10.794;

% Active detection prior
pri.p2.a = 3.212;
pri.p2.b = 2.822;

% Prior strength multiplier
pri.strength = 1.0;

%% Define estimation problem
% Empty struct means estimate all parameters
fix = struct();

%% Multi-start optimization
tic;

R0_starts = [0.25, 0.50, 0.75];
k_starts  = [0.50, 1.00, 1.50];
p1_starts = [0.20, 0.40, 0.70];
p2_starts = [0.20, 0.40, 0.70];

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
        observedClusterData, fix, pri, initial_guess, 'off');
end

obj_vals = cellfun(@(s) s.obj, all_results);
[~, best_idx] = min(obj_vals);
est = all_results{best_idx};

elapsedTime = toc;

%% Report results
fprintf('\nMAP estimation completed in %.2f seconds.\n', elapsedTime);

fprintf('\n--- Point Estimates ---\n');
fprintf('R_eff = %.4f\n', est.R0);
fprintf('k     = %.4f\n', est.k);
fprintf('p_p   = %.4f\n', est.p_p);
fprintf('p_a   = %.4f\n', est.p_a);

fprintf('\n--- 95%% Confidence Intervals ---\n');
fprintf('R_eff: [%.4f, %.4f]\n', est.CI.R0(1), est.CI.R0(2));
fprintf('k    : [%.4f, %.4f]\n', est.CI.k(1), est.CI.k(2));
fprintf('p_p  : [%.4f, %.4f]\n', est.CI.p_p(1), est.CI.p_p(2));
fprintf('p_a  : [%.4f, %.4f]\n', est.CI.p_a(1), est.CI.p_a(2));

estimated_CV = sqrt(est.R0 + est.R0^2 / est.k) / est.R0;
fprintf('\n--- Derived Quantity ---\n');
fprintf('Estimated coefficient of variation (CV) = %.4f\n', estimated_CV);

end


%% ========================================================================
% MAP estimation helper
% ========================================================================

function est = fit_observed_no_k_prior(counts, fix, pri, initial_guess, display_option)

if nargin < 5
    display_option = 'iter';
end

to_est.R0  = ~isfield(fix,'R0')  || isempty(fix.R0);
to_est.k   = ~isfield(fix,'k')   || isempty(fix.k);
to_est.p_p = ~isfield(fix,'p_p') || isempty(fix.p_p);
to_est.p_a = ~isfield(fix,'p_a') || isempty(fix.p_a);

R0_0 = pick_default(fix, 'R0',  initial_guess.R0);
k_0  = pick_default(fix, 'k',   initial_guess.k);
p1_0 = pick_default(fix, 'p_p', initial_guess.p_p);
p2_0 = pick_default(fix, 'p_a', initial_guess.p_a);

theta0 = [];
lb = [];
ub = [];

if to_est.R0
    theta0(end+1) = log(R0_0);
    lb(end+1) = -inf;
    ub(end+1) = 0;              % ensures R0 <= 1
end
if to_est.k
    theta0(end+1) = log(k_0);
    lb(end+1) = -inf;
    ub(end+1) = inf;
end
if to_est.p_p
    theta0(end+1) = logit(p1_0);
    lb(end+1) = logit(0.1);     % lower bound on passive detection
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

obj = @(th) obj_single(th, counts, fix, to_est, pri);
[th_hat, fval, ~, ~, ~, ~, hessian] = fmincon(obj, theta0, [], [], [], [], lb, ub, [], opts);

[R0, k, p_p, p_a] = unpack_params(th_hat, fix, to_est);
est = struct('R0', R0, 'k', k, 'p_p', p_p, 'p_a', p_a, 'obj', fval);

try
    var_cov_matrix_transformed = inv(hessian);

    diag_variances = diag(var_cov_matrix_transformed);
    if any(diag_variances < 0)
        warning('Negative variance detected from Hessian. Using absolute values. Confidence intervals may be unreliable.');
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
    warning('Could not calculate confidence intervals: %s', ME.message);
    if to_est.R0,  est.CI.R0  = [NaN, NaN]; end
    if to_est.k,   est.CI.k   = [NaN, NaN]; end
    if to_est.p_p, est.CI.p_p = [NaN, NaN]; end
    if to_est.p_a, est.CI.p_a = [NaN, NaN]; end
end

end


%% ========================================================================
% Objective function
% ========================================================================

function val = obj_single(th, counts, fix, to_est, pri)

[R0, k, p_p, p_a] = unpack_params(th, fix, to_est);

if ~(R0 >= 0 && k > 1e-6 && p_p > 0 && p_p < 1 && p_a > 0 && p_a < 1)
    val = 1e12;
    return;
end

num_clusters = size(counts, 1);
log_px_vector = zeros(num_clusters, 1);

for i = 1:num_clusters
    cluster_size = counts(i, 1);
    num_index    = counts(i, 2);

    px = observed_size_pmf_condpos_FAST(R0, k, p_p, p_a, cluster_size, num_index);

    prob_of_size = px(cluster_size);
    if prob_of_size <= 0
        prob_of_size = realmin;
    end

    log_px_vector(i) = log(prob_of_size);
end

nll = -sum(counts(:, 3) .* log_px_vector);

val = nll + prior_penalty(p_p, p_a, pri);

end


%% ========================================================================
% Prior penalty
% ========================================================================

function pen = prior_penalty(p_p, p_a, pri)

nlogprior_p1 = -((pri.p1.a - 1) * log(max(p_p, eps)) + ...
                 (pri.p1.b - 1) * log(max(1 - p_p, eps))) + ...
                 betaln(pri.p1.a, pri.p1.b);

nlogprior_p2 = -((pri.p2.a - 1) * log(max(p_a, eps)) + ...
                 (pri.p2.b - 1) * log(max(1 - p_a, eps))) + ...
                 betaln(pri.p2.a, pri.p2.b);

pen = pri.strength * (nlogprior_p1 + nlogprior_p2);

end


%% ========================================================================
% Utilities
% ========================================================================

function val = pick_default(fix, field, def)
if isfield(fix, field) && ~isempty(fix.(field))
    val = fix.(field);
else
    val = def;
end
end

function y = logit(p)
y = log(p) - log(1 - p);
end

function p = logistic(y)
p = 1 ./ (1 + exp(-y));
end

function [R0, k, p_p, p_a] = unpack_params(th, fix, to_est)
idx = 1;

if to_est.R0
    R0 = exp(th(idx));
    idx = idx + 1;
else
    R0 = fix.R0;
end

if to_est.k
    k = exp(th(idx));
    idx = idx + 1;
else
    k = fix.k;
end

if to_est.p_p
    p_p = logistic(th(idx));
    idx = idx + 1;
else
    p_p = fix.p_p;
end

if to_est.p_a
    p_a = logistic(th(idx));
else
    p_a = fix.p_a;
end
end


%% ========================================================================
% Conditional observed-size PMF
% ========================================================================

function px = observed_size_pmf_condpos_FAST(R0, k, p_p, p_a, Xmax, n_index)

n = n_index;
Jmax = max(2000, 10 * Xmax);

pJ = true_cluster_pmf(R0, k, n, Jmax);
p_obs = p_p + (1 - p_p) * p_a;

x_vec = (1:Xmax)';
j_vec = n:Jmax;

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

if n == 1
    massX_positive = P_O_given_S * pJ';
    p_O_zero = sum((1 - p_p) .^ j_vec .* pJ);
else
    massX_positive = [zeros(Xmax, n - 1), P_O_given_S] * pJ';
    p_O_zero = sum(([zeros(1, n - 1), (1 - p_p) .^ j_vec]) .* pJ);
end

p_pos = 1 - p_O_zero;
if p_pos <= 0
    p_pos = realmin;
end

px = massX_positive' / p_pos;

end


function pJ = true_cluster_pmf(R0, k, n, Jmax)

pJ = zeros(1, Jmax);

if R0 == 0
    pJ(n) = 1;
    pJ = pJ / sum(pJ);
    return;
end

jj = n:Jmax;
valid_mask = (k * jj + jj - n) > 0 & (jj - n + 1) > 0;
jj_valid = jj(valid_mask);

if isempty(jj_valid)
    return;
end

logRk  = log(R0 / k);
log1Rk = log(1 + R0 / k);

val = log(n) - log(jj_valid) ...
    + gammaln(k * jj_valid + jj_valid - n) ...
    - gammaln(k * jj_valid) ...
    - gammaln(jj_valid - n + 1) ...
    + (jj_valid - n) .* logRk ...
    - (k * jj_valid + jj_valid - n) .* log1Rk;

pJ(jj_valid) = exp(val);

end
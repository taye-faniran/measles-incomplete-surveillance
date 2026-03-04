function Figure8_generate_active_detection_prior
% Generate the hierarchical Beta-Binomial prior for active detection p_a.
%
% This script pools evidence from multiple studies using a hierarchical
% Beta-Binomial model:
%
%   s_i | n_i, p_i ~ Binomial(n_i, p_i)
%   p_i          ~ Beta(alpha_H, beta_H)
%
% The pooled hyperparameters (alpha_H, beta_H) are estimated with a
% random-walk Metropolis-Hastings algorithm on the log scale.
%
% Current input data correspond to the active-detection evidence used in
% the manuscript:
%   1) UK outbreak investigation:      54 detected, 10 missed
%   2) Catalonia comparison study:     41 detected, 130 missed
%   3) Weakly informative Beta(2,2):    2 pseudo-successes, 2 pseudo-failures
%
% Outputs:
%   - posterior mean and 95% credible interval for alpha_H and beta_H
%   - recommended pooled Beta prior for p_a
%   - diagnostic plots
%   - saved MAT file with posterior samples and summaries

clear; clc; close all;
rng(2025);

baseFolder = fileparts(mfilename('fullpath'));
outFile    = fullfile(baseFolder, 'Figure8_active_detection_prior_results.mat');

%% Input data: active-detection evidence
successes = [54, 41, 2];
failures  = [10, 130, 2];
totals    = successes + failures;

%% Hyperpriors for alpha_H and beta_H
prior_shape = 1.0;
prior_rate  = 0.1;

%% MCMC settings
niter = 100000;
burn  = 50000;
thin  = 5;

prop_sd_loga = 0.12;
prop_sd_logb = 0.12;

%% Initial values
loga_curr = log(5);
logb_curr = log(5);

lp_curr = logpost_logab_fun( ...
    loga_curr, logb_curr, successes, totals, prior_shape, prior_rate);

chain = zeros(niter, 2);
lp_store = zeros(niter, 1);
accept = 0;

%% Metropolis-Hastings sampling
for t = 1:niter
    loga_prop = loga_curr + prop_sd_loga * randn;
    logb_prop = logb_curr + prop_sd_logb * randn;

    lp_prop = logpost_logab_fun( ...
        loga_prop, logb_prop, successes, totals, prior_shape, prior_rate);

    if isfinite(lp_prop) && log(rand) < (lp_prop - lp_curr)
        loga_curr = loga_prop;
        logb_curr = logb_prop;
        lp_curr = lp_prop;
        accept = accept + 1;
    end

    chain(t, :) = [loga_curr, logb_curr];
    lp_store(t) = lp_curr;

    if mod(t, 10000) == 0
        fprintf('Iteration %d of %d | acceptance rate = %.3f\n', ...
            t, niter, accept / t);
    end
end

accept_rate = accept / niter;
fprintf('\nFinal acceptance rate = %.3f\n', accept_rate);

%% Posterior samples after burn-in and thinning
chain_post = chain(burn+1:thin:end, :);
a_samples = exp(chain_post(:,1));
b_samples = exp(chain_post(:,2));

%% Posterior summaries
alpha_mean = mean(a_samples);
beta_mean  = mean(b_samples);

alpha_ci = quantile(a_samples, [0.025 0.975]);
beta_ci  = quantile(b_samples, [0.025 0.975]);

fprintf('\nPosterior mean alpha_H = %.3f (95%% CrI [%.3f, %.3f])\n', ...
    alpha_mean, alpha_ci(1), alpha_ci(2));
fprintf('Posterior mean beta_H  = %.3f (95%% CrI [%.3f, %.3f])\n', ...
    beta_mean, beta_ci(1), beta_ci(2));

fprintf('\nRecommended pooled prior for p_a: Beta(%.3f, %.3f)\n', ...
    alpha_mean, beta_mean);

%% Posterior predictive draws for the population-level detection probability
nsamps = 5000;
idx = randi(numel(a_samples), nsamps, 1);
p_draws = betarnd(a_samples(idx), b_samples(idx));

predictive_mean = mean(p_draws);
predictive_ci   = quantile(p_draws, [0.025 0.975]);

fprintf('\nPosterior predictive mean for p_a = %.3f\n', predictive_mean);
fprintf('Posterior predictive 95%% interval for p_a = [%.3f, %.3f]\n', ...
    predictive_ci(1), predictive_ci(2));

%% Diagnostic figure: traces and posterior histograms
fig1 = figure('Color', 'w', 'Position', [100 100 1000 700]);

subplot(2,2,1);
plot(exp(chain(:,1)), 'LineWidth', 1);
title('Trace plot: \alpha_H');
xlabel('Iteration');
ylabel('\alpha_H');

subplot(2,2,2);
plot(exp(chain(:,2)), 'LineWidth', 1);
title('Trace plot: \beta_H');
xlabel('Iteration');
ylabel('\beta_H');

subplot(2,2,3);
histogram(a_samples, 50);
title('Posterior of \alpha_H');
xlabel('\alpha_H');
ylabel('Frequency');

subplot(2,2,4);
histogram(b_samples, 50);
title('Posterior of \beta_H');
xlabel('\beta_H');
ylabel('Frequency');

%% Posterior predictive figure
fig2 = figure('Color', 'w', 'Position', [150 150 700 500]);
histogram(p_draws, 50);
title('Posterior predictive distribution of p_a');
xlabel('p_a');
ylabel('Frequency');

%% Save results
save(outFile, ...
    'successes', 'failures', 'totals', ...
    'prior_shape', 'prior_rate', ...
    'niter', 'burn', 'thin', ...
    'prop_sd_loga', 'prop_sd_logb', ...
    'chain', 'lp_store', 'accept_rate', ...
    'a_samples', 'b_samples', ...
    'alpha_mean', 'beta_mean', ...
    'alpha_ci', 'beta_ci', ...
    'p_draws', 'predictive_mean', 'predictive_ci');

fprintf('\nSaved results to:\n%s\n', outFile);

end

function lp = logpost_logab_fun(loga, logb, successes, totals, prior_shape, prior_rate)
% Log-posterior for log(alpha_H), log(beta_H).

alpha_H = exp(loga);
beta_H  = exp(logb);

% Beta-Binomial marginal likelihood
loglik = sum( ...
    betaln(successes + alpha_H, totals - successes + beta_H) ...
    - betaln(alpha_H, beta_H) ...
    + gammaln(totals + 1) ...
    - gammaln(successes + 1) ...
    - gammaln(totals - successes + 1));

% Gamma hyperpriors on alpha_H and beta_H
logprior = (prior_shape - 1) * log(alpha_H) - prior_rate * alpha_H ...
         + (prior_shape - 1) * log(beta_H)  - prior_rate * beta_H;

% Jacobian for the log transformation
logjac = loga + logb;

lp = loglik + logprior + logjac;
end
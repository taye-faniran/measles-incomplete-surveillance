% Figure7_generate_passive_detection_prior.m
% ------------------------------------------------------------
% Figure 7: Hierarchical Beta-Binomial prior for PASSIVE detection p_p
%
% Model:
%   s_i | n_i, p_i ~ Binomial(n_i, p_i)
%   p_i          ~ Beta(alpha_H, beta_H)
% Hyperpriors:
%   alpha_H ~ Gamma(prior_shape, prior_rate)
%   beta_H  ~ Gamma(prior_shape, prior_rate)
%
% We sample (log(alpha_H), log(beta_H)) via random-walk Metropolis-Hastings.
%
% Data below are the passive-detection evidence used in Appendix B:
%   successes = reported (detected) cases
%   failures  = undetected cases
% ------------------------------------------------------------

clear; clc; close all;
rng(2025);

%% Folder containing this script
baseFolder = fileparts(mfilename('fullpath'));

outMat       = fullfile(baseFolder, 'Figure7_passive_detection_prior_results.mat');
outPdfTraces = fullfile(baseFolder, 'Figure7_passive_detection_prior_traces.pdf');
outPdfPred   = fullfile(baseFolder, 'Figure7_passive_detection_prior_predictive.pdf');

%% Input data: passive-detection evidence (successes, failures)
successes = [148, 745, 10, 2, 0];
failures  = [420, 9255, 25, 66, 8];
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

%% Initial values (on log scale)
loga_curr = log(5);
logb_curr = log(5);

lp_curr = logpost_logab_fun(loga_curr, logb_curr, successes, totals, prior_shape, prior_rate);

chain    = zeros(niter, 2);
lp_store = zeros(niter, 1);
accept   = 0;

%% Metropolis-Hastings loop
for t = 1:niter
    loga_prop = loga_curr + prop_sd_loga * randn;
    logb_prop = logb_curr + prop_sd_logb * randn;

    lp_prop = logpost_logab_fun(loga_prop, logb_prop, successes, totals, prior_shape, prior_rate);

    if isfinite(lp_prop) && log(rand) < (lp_prop - lp_curr)
        loga_curr = loga_prop;
        logb_curr = logb_prop;
        lp_curr   = lp_prop;
        accept    = accept + 1;
    end

    chain(t,:)    = [loga_curr, logb_curr];
    lp_store(t,1) = lp_curr;

    if mod(t, 10000) == 0
        fprintf('Iteration %d of %d | acceptance rate = %.3f\n', t, niter, accept/t);
    end
end

accept_rate = accept / niter;
fprintf('\nFinal acceptance rate = %.3f\n', accept_rate);

%% Posterior samples after burn-in and thinning
chain_post = chain(burn+1:thin:end, :);
a_samples  = exp(chain_post(:,1));
b_samples  = exp(chain_post(:,2));

%% Posterior summaries for alpha_H, beta_H
alpha_mean = mean(a_samples);
beta_mean  = mean(b_samples);
alpha_ci   = quantile(a_samples, [0.025 0.975]);
beta_ci    = quantile(b_samples, [0.025 0.975]);

fprintf('\nPosterior mean alpha_H = %.3f (95%% CrI [%.3f, %.3f])\n', ...
    alpha_mean, alpha_ci(1), alpha_ci(2));
fprintf('Posterior mean beta_H  = %.3f (95%% CrI [%.3f, %.3f])\n', ...
    beta_mean, beta_ci(1), beta_ci(2));

fprintf('\nRecommended pooled prior for p_p: Beta(%.3f, %.3f)\n', alpha_mean, beta_mean);

%% Posterior predictive draws for p_p
nsamps = 5000;
idx = randi(numel(a_samples), nsamps, 1);
p_draws = betarnd(a_samples(idx), b_samples(idx));

pred_mean = mean(p_draws);
pred_ci   = quantile(p_draws, [0.025 0.975]);

fprintf('\nPosterior predictive mean for p_p = %.3f\n', pred_mean);
fprintf('Posterior predictive 95%% interval for p_p = [%.3f, %.3f]\n', pred_ci(1), pred_ci(2));

%% Figure 7A: trace + posterior histograms (alpha_H, beta_H)
fig1 = figure('Color', 'w', 'Position', [100 100 1000 700]);

subplot(2,2,1);
plot(exp(chain(:,1)), 'LineWidth', 1);
title('Trace plot: \alpha_H'); xlabel('Iteration'); ylabel('\alpha_H');

subplot(2,2,2);
plot(exp(chain(:,2)), 'LineWidth', 1);
title('Trace plot: \beta_H'); xlabel('Iteration'); ylabel('\beta_H');

subplot(2,2,3);
histogram(a_samples, 50);
title('Posterior of \alpha_H'); xlabel('\alpha_H'); ylabel('Frequency');

subplot(2,2,4);
histogram(b_samples, 50);
title('Posterior of \beta_H'); xlabel('\beta_H'); ylabel('Frequency');

set(fig1, 'PaperPositionMode', 'auto');
print(fig1, outPdfTraces, '-dpdf', '-bestfit');

%% Figure 7B: posterior predictive distribution for p_p
fig2 = figure('Color', 'w', 'Position', [150 150 700 500]);
histogram(p_draws, 50);
title('Posterior predictive distribution of p_p');
xlabel('p_p'); ylabel('Frequency');

set(fig2, 'PaperPositionMode', 'auto');
print(fig2, outPdfPred, '-dpdf', '-bestfit');

%% Save results to .mat (optional but useful for reproducibility)
save(outMat, ...
    'successes','failures','totals', ...
    'prior_shape','prior_rate', ...
    'niter','burn','thin', ...
    'prop_sd_loga','prop_sd_logb', ...
    'chain','lp_store','accept_rate', ...
    'a_samples','b_samples', ...
    'alpha_mean','beta_mean','alpha_ci','beta_ci', ...
    'p_draws','pred_mean','pred_ci');

fprintf('\nSaved PDFs:\n%s\n%s\n', outPdfTraces, outPdfPred);
fprintf('Saved MAT results:\n%s\n', outMat);

%% ---- local function ----
function lp = logpost_logab_fun(loga, logb, successes, totals, prior_shape, prior_rate)
    alpha_H = exp(loga);
    beta_H  = exp(logb);

    % Beta-Binomial marginal likelihood
    loglik = sum( ...
        betaln(successes + alpha_H, totals - successes + beta_H) ...
        - betaln(alpha_H, beta_H) ...
        + gammaln(totals + 1) ...
        - gammaln(successes + 1) ...
        - gammaln(totals - successes + 1));

    % Gamma hyperpriors
    logprior = (prior_shape - 1) * log(alpha_H) - prior_rate * alpha_H ...
             + (prior_shape - 1) * log(beta_H)  - prior_rate * beta_H;

    % Jacobian for log-transform
    logjac = loga + logb;

    lp = loglik + logprior + logjac;
end
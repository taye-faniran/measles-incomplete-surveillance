function Figure2_generate_unknown_detection
% Generate results for Figure 2:
% unknown detection leads to structural non-identifiability.

clear; clc;

baseFolder = fileparts(mfilename('fullpath'));
load(fullfile(baseFolder, 'Figure2_unknown_detection_results.mat'), 'XX');

index = 1;
r0 = 0.5;
k = 0.15;

RR0 = NaN(4,4);
kk = NaN(4,4);
RR0_lower = NaN(4,4);
RR0_upper = NaN(4,4);
kk_lower = NaN(4,4);
kk_upper = NaN(4,4);
CV_est = NaN(4,4);

for mm = 1:4
    for ii = 1:4
        fprintf('Running grid cell: row = %d, col = %d\n', mm, ii);

        active = 0.25 * mm;
        passive = 0.25 * ii;

        fprintf('   Passive detection p_p = %.2f\n', passive);
        fprintf('   Active detection  p_a = %.2f\n', active);

        % First simulate a reference outbreak-size distribution close to the
        % target truth (R0 = 0.5, k = 0.15).
        opt_arr = [0, 0];
        maxIter = 1e5;
        iter = 0;

        while (opt_arr(1) < r0 - 0.002 || opt_arr(1) > r0 + 0.002 || ...
               opt_arr(2) < (1 - 0.002)*k || opt_arr(2) > 1.002*k)

            iter = iter + 1;

            CD = cluster_distribution_2tier_pobs(r0, k, index, -1, [1 1], 1e6, 10000);
            opt_arr = optimize_cluster_dist043010(CD(:,1)', 0, size(CD,1), ...
                      [0.5 1 1 index 1 1], [1 2]);

            if iter >= maxIter
                warning('Stopped after %d iterations without hitting both target ranges.', maxIter);
                break;
            end
        end

        % Apply the chosen passive/active detection setting.
        Original_Cluster_Distribution = CD(:,1);
        Observed_Cluster_Distribution = cluster_distribution_pobs( ...
            Original_Cluster_Distribution, [passive, active]);

        % Fit all four parameters with detection unknown.
        est = fit_observed(Observed_Cluster_Distribution, ...
              struct('R0',[],'k',[],'p_p',[],'p_a',[]));

        RR0(mm,ii) = est.R0;
        kk(mm,ii) = est.k;
        RR0_lower(mm,ii) = est.R0_CI(1);
        RR0_upper(mm,ii) = est.R0_CI(2);
        kk_lower(mm,ii) = est.k_CI(1);
        kk_upper(mm,ii) = est.k_CI(2);
        CV_est(mm,ii) = sqrt(est.R0 + est.R0^2/est.k) / est.R0;

        fprintf('   Estimated R_eff = %.4f\n', est.R0);
        fprintf('   Estimated k     = %.4f\n', est.k);
        fprintf('   95%% CI for R_eff = [%.4f, %.4f]\n', est.R0_CI(1), est.R0_CI(2));
        fprintf('   95%% CI for k     = [%.4f, %.4f]\n', est.k_CI(1), est.k_CI(2));
        fprintf('   Estimated CV     = %.4f\n\n', CV_est(mm,ii));
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
disp('1 = R_eff estimates');
disp('2 = k estimates');
disp('3 = R_eff lower 95% CI');
disp('4 = R_eff upper 95% CI');
disp('5 = k lower 95% CI');
disp('6 = k upper 95% CI');
disp('7 = estimated CV');

save(folderPath, 'XX');
fprintf('Saved results to:\n%s\n', folderPath);

end


function est = fit_observed(counts, fix)

counts = counts(:)';
Xmax = numel(counts);

to_est.R0  = ~isfield(fix,'R0')  || isempty(fix.R0);
to_est.k   = ~isfield(fix,'k')   || isempty(fix.k);
to_est.p_p = ~isfield(fix,'p_p') || isempty(fix.p_p);
to_est.p_a = ~isfield(fix,'p_a') || isempty(fix.p_a);

R0_0 = pick_default(fix,'R0',0.6);
k_0  = pick_default(fix,'k',1.0);
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

obj = @(th) obj_single(th, counts, fix, to_est, Xmax);
opts = optimoptions('fmincon','Display','off','Algorithm','interior-point');
[th_hat, fval, ~, ~, ~, ~, hessian] = fmincon(obj, theta0, [], [], [], [], lb, ub, [], opts);

[R0,k,p_p,p_a] = unpack_params(th_hat, fix, to_est);

ci.R0 = [NaN, NaN];
ci.k = [NaN, NaN];

if ~isempty(hessian)
    try
        cov_th = inv(hessian);

        if any(diag(cov_th) < 0)
            error('Covariance matrix has negative diagonal elements.');
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
else
    warning('Hessian not returned by fmincon. Cannot compute confidence intervals.');
end

est = struct('R0',R0,'k',k,'p_p',p_p,'p_a',p_a,'obj',fval,'Xmax',Xmax, ...
             'R0_CI',ci.R0,'k_CI',ci.k);
end


function val = obj_single(th, counts, fix, to_est, Xmax)
[R0,k,p_p,p_a] = unpack_params(th, fix, to_est);

if ~(R0>=0 && k>1e-6 && p_p>=0 && p_p<=1 && p_a>=0 && p_a<=1)
    val = 1e12;
    return
end

px = observed_size_pmf_condpos(R0,k,p_p,p_a,Xmax);
px(px<=0) = realmin;
val = -sum(counts .* log(px));
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
        continue
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
    return
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
            new_case_ids = (height(cases_table)+1):(height(cases_table)+total_new_cases);
            new_generations = repmat(current_generation+1, total_new_cases, 1);
            new_cases_table = table(new_case_ids', new_generations, ...
                'VariableNames', {'case_id','generation'});
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


function cluster_dist = cluster_distribution_2tier_pobs(r0,k,n_ps,k_ps,pobs,max_outbreak_size,num_sim)

cluster_dist = zeros(1,2);
n = 0;

if ~(k_ps == -1 || n_ps == 1)
    mu = fzero(@(x) x - n_ps*(1-nbinpdf(0,k_ps,1/(1+x/k_ps))), n_ps - 1);
end

while n < num_sim
    if (k_ps == -1 || n_ps == 1)
        popsize = n_ps;
    else
        popsize = 0;
        while popsize == 0
            popsize = gen_nb_dist(mu,k_ps,1,1);
        end
    end

    outbreak_size = popsize;
    while (popsize > 0 && outbreak_size < max_outbreak_size)
        nb_dist = gen_nb_dist(r0,k,popsize,1);
        popsize = sum(nb_dist);
        outbreak_size = outbreak_size + popsize;
    end

    outbreak_size1 = binornd(outbreak_size,pobs(1));
    if outbreak_size1 > 0
        outbreak_size = outbreak_size1 + binornd(outbreak_size-outbreak_size1,pobs(2));
    else
        outbreak_size = 0;
    end

    if outbreak_size > 0
        if outbreak_size > size(cluster_dist,1)
            cluster_dist(outbreak_size,1) = 0;
        end
        cluster_dist(outbreak_size,1) = cluster_dist(outbreak_size,1)+1;
        if popsize > 0
            cluster_dist(outbreak_size,2) = cluster_dist(outbreak_size,2)+1;
        end
        n = n + 1;
    end
end

XX = cluster_dist(:,1);
XY = cluster_dist(:,2);

XX_cluster = XX(n_ps:length(XX));
XY_cluster = XY(n_ps:length(XX));

cluster_dist = [XX_cluster, XY_cluster];
end


function nb_dist = gen_nb_dist(r0,k,m,n)
nb_dist = poissrnd(gamrnd(k,r0/k,m,n));
end


function opt_par_arr = optimize_cluster_dist043010(outbreak_data, num_nonextinct, max_outbreak_size, init_par_arr, opt_indices)

x = init_par_arr(opt_indices);
OPTIONS = optimset('MaxFunEvals',5000,'MaxIter',5000);
[opt_x,nll_min] = fminsearch(@(x) -mle_call(outbreak_data,num_nonextinct,max_outbreak_size,init_par_arr,x,opt_indices), x, OPTIONS);

opt_par_arr = init_par_arr;
opt_par_arr(opt_indices) = opt_x;
opt_par_arr(end+1) = nll_min;
end


function L = mle_call(outbreak_data,num_nonextinct,max_outbreak_size,init_par_arr,x,opt_indices)
par_arr = init_par_arr;
par_arr(opt_indices) = x;
L = mle_cluster_dist031721(outbreak_data,num_nonextinct,max_outbreak_size,par_arr);
end


function L = mle_cluster_dist031721(outbreak_data,num_nonextinct,max_outbreak_size,par_arr)

if (par_arr(1) < 0 || par_arr(2) < 1e-5 || par_arr(2) > 1e3 || par_arr(3) < 1 ...
        || par_arr(4) > 1e4 || (par_arr(3) < 0 && par_arr(3) ~= -1) || par_arr(6) < 0 || par_arr(6) > 1)
    L = -Inf;
    return
end

XX = length(outbreak_data) + par_arr(4) - 1;

if par_arr(5) == 1
    clust_pdf = calc_cluster_pdf(par_arr, XX);
else
    clust_pdf = calc_cluster_pdf1(par_arr, XX);
end

clust_pdf = clust_pdf(par_arr(4):XX);

if par_arr(5) == 3
else
    L = sum(outbreak_data .* log(clust_pdf(1:length(outbreak_data)))) + ...
        num_nonextinct * log(1 - sum(clust_pdf));
end

if isnan(L) || ~isreal(L) || ~isfinite(L)
    L = -1e300;
end
end


function clust_pdf = calc_cluster_pdf(par_arr,outbreak_size_limit)

r0 = par_arr(1);
k = max(0,par_arr(2));
pobs_mode = par_arr(5);
pobs = par_arr(6);
m = par_arr(4);

if (pobs_mode == 1 && pobs < 1)
    num_calc = min(100*outbreak_size_limit,1000);
else
    num_calc = outbreak_size_limit;
end

if r0 == 0
    true_clust_pdf = zeros(1,num_calc);
    true_clust_pdf(1) = 1;
else
    j = 1:num_calc;
    log_real_clust_pdf = log(m) - log(j) + gammaln(max(k*j+j-m,0)) ...
        - gammaln(max(j-m+1,0)) - gammaln(k*j) + (j-m)*log(r0/k) ...
        - (max(k*j+j-m,0))*log(1+r0/k);
    true_clust_pdf = exp(log_real_clust_pdf);
end

if (pobs_mode == 0 || pobs == 1)
    clust_pdf = true_clust_pdf;
    return
end

j = 1:num_calc;
prob0 = sum(exp(j*log(1-pobs) + log(true_clust_pdf)));
denominator = 1 - prob0;

switch pobs_mode
    case 1
        numerator = zeros(1,outbreak_size_limit);
        for jj = 1:outbreak_size_limit
            l = jj:length(true_clust_pdf);
            numerator(jj) = exp(jj*log(pobs/(1-pobs)) - gammaln(jj+1)) * ...
                sum(exp(log(true_clust_pdf(l)) + l*log(1-pobs) + gammaln(l+1) - gammaln(l+1-jj)));
        end
        clust_pdf = numerator / denominator;

    case 2
        numerator = exp(log(true_clust_pdf(j))) .* (1-(1-pobs).^j);
        clust_pdf = numerator(1:outbreak_size_limit) / denominator;
end
end


function clust_pdf = calc_cluster_pdf1(~,~)
error('calc_cluster_pdf1 is not available in this workflow.');
end


function cluster_dist = cluster_distribution_pobs(Cluster_dist,pobs)

cluster_dist = zeros(1);

FF = (1:length(Cluster_dist))';
GG = [FF, Cluster_dist];
expanded_column = repelem(GG(:,1), GG(:,2));

for jj = 1:length(expanded_column)
    outbreak_size = expanded_column(jj);

    outbreak_size1 = binornd(outbreak_size,pobs(1));
    if outbreak_size1 > 0
        outbreak_size = outbreak_size1 + binornd(outbreak_size-outbreak_size1,pobs(2));
    else
        outbreak_size = 0;
    end

    if outbreak_size > 0
        if outbreak_size > size(cluster_dist,1)
            cluster_dist(outbreak_size,1) = 0;
        end
        cluster_dist(outbreak_size,1) = cluster_dist(outbreak_size,1)+1;
    end
end

XX = cluster_dist(:,1);
cluster_dist = XX(1:length(XX));
end
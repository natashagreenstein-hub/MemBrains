%% NEMO_oxglc.m — Simplified NEMO with metabolic (CMRO2/CMRglc) noise scaling
%                  ** DBS80 version **  (N = 80 regions)
%
%   This is the cleaned-up entry point that *replaces* the original NEMO
%   step-4 (PSO-fit β over neurotransmitter receptor maps) with a single
%   per-region scalar:
%
%       betax_i = beta_baseline * ox_gluc_ratio_i      (i = 1..N regions)
%
%   We pivoted from Schaefer1000 → DBS80 because:
%     - GEC fitting at N=1000 is slow and memory-heavy.
%     - All the upstream group-level pieces (HCP rest timeseries, SC,
%       CMRO2/CMRglc maps) already exist at DBS80 in MemBrain.
%     - The ox_gluc_ratio map will be parcellated to DBS80 before being
%       dropped in here.
%
%   Pipeline:
%     0.  setup paths, atlas size (DBS80, N=80), filter
%     1.  load HCP rest DBS80 timeseries (498 subj × 1200 TR × 80 ROI),
%         compute group-mean FCemp, COVtauemp.
%         f_diff is precomputed in empirical_HCP_rest_DBS80.mat — use it.
%     2.  load SC_dbs80HARDI (80×80) and fit Ceff via RS-GEC
%         [optional toggle, default ON because N=80 is fast]
%     3.  load ox_gluc_ratio (per-region noise weight, length N=80)
%     4.  scale β per region: betax = beta .* ox_gluc_ratio   (← simplified
%         step 4; the PSO + receptor-map block is removed entirely)
%     5.  simulate Hopf with the metabolically-scaled noise; compare to FCemp
%     6.  save results
%
%   Required files in /Users/natashag/Documents/Claude/Projects/MemBrain/:
%       HCP_rest_DBS80_Functional_Timeseries.mat   — `bigmat` (80×1200×498), `TR`
%       empirical_HCP_rest_DBS80.mat               — `f_diff` (1×80) precomputed
%       SC_dbs80HARDIFULL.mat                      — `SC_dbs80HARDI` (80×80)
%       CMRO2_over_CMRglc_DBS80.mat                — fallback for ox_gluc_ratio
%       ox_gluc_ratio.mat                          — variable `ox_gluc_ratio`
%                                                    (1×80, you provide this)

% Clear everything except `metabolic_map` and `ceff_source` so a wrapper
% script (e.g. run_all_NEMO_maps.m) can set them before invoking us;
% standalone runs behave exactly as before because neither var will exist
% yet (defaults are resolved below).
clearvars -except metabolic_map ceff_source;

projDir = '/Users/natashag/Documents/Claude/Projects/MemBrain';
nemoDir = fullfile(projDir, 'NEMO');
gecDir  = fullfile(projDir, 'GEC Scripts');
addpath(genpath(projDir));

% --- atlas / model parameters ---
N      = 80;            % DBS80 parcels (use all 80; for the 62-region subset
                        % uncomment indexN below)
% indexN = [1:31, 50:80];   % optional 62-region subset (matches Ceff_Computation_MemBrain.m)
indexN = 1:N;

beta  = 0.01;          % baseline β (homogeneous noise amplitude)
Tau    = 3;             % lag (TRs) for COV(τ)
TR     = 0.72;          % HCP TR (s) — also stored in the .mat file

% --- bandpass filter (same band as the original NEMO/GEC scripts) ---
fnq = 1/(2*TR);
flp = 0.008;  fhi = 0.08;
Wn  = [flp/fnq, fhi/fnq];
[bfilt, afilt] = butter(2, Wn);

% --- runtime toggles ---
do_gec_fit   = true;    % fit Ceff via legacy RS-GEC (fast at N=80)
                        % Only consulted when ceff_source == 'legacy'.
max_gec_iter = 2000;

% Which algorithm produces Ceff for downstream Step 4 PSO.
%   'legacy'           : original NEMO cooperative-only RS-GEC (Lyapunov-based,
%                        Cnew >= 0, fits group-mean FCemp + COVtauemp). DEFAULT
%                        — byte-identical to the prior pipeline.
%   'gec_competitive'  : Luppi et al 2026 GEC allowing negative weights
%                        (competitive interactions), fit by
%                        fit_competitive_gec_NEMO.m. The fit is cached to a
%                        sibling .mat so re-runs don't redo it.
%   'gec_cooperative'  : same Luppi 2026 algorithm but cooperative-only
%                        (negative GEC values clamped to 0). Useful as a
%                        like-for-like A/B against 'gec_competitive'.
% Respect a value already set by a wrapper; else default to legacy.
if ~exist('ceff_source', 'var') || isempty(ceff_source)
    ceff_source = 'legacy';
end
switch ceff_source
    case {'legacy', 'gec_competitive', 'gec_cooperative'}
        % valid
    otherwise
        error(['ceff_source must be one of ''legacy'', ''gec_competitive'', ' ...
               'or ''gec_cooperative'' (got ''%s'').'], ceff_source);
end

% --- which metabolic map drives betax ---
%   'ox_gluc' : oxphos / glucose ratio   (DEFAULT — unchanged behaviour;
%               loads ox_gluc_ratio.mat or falls back to ratio_values)
%   'ox'      : CMRO2 only   (cmro2_values  in metabolic_maps_DBS80.mat)
%   'glc'     : CMRglc only  (cmrglc_values in metabolic_maps_DBS80.mat)
% In all three cases the loaded vector is normalised to mean 1, so the
% same PSO bounds and seed apply. The output filename gets a matching
% suffix so different runs don't overwrite each other.
% Respect a value already set by a wrapper (run_all_NEMO_maps.m); else default.
if ~exist('metabolic_map', 'var') || isempty(metabolic_map)
    metabolic_map = 'ox_gluc';   % headline map by default (was 'glc', the null —
                                 % a bare run then looked like "no effect / 0.6472").
                                 % run_all_NEMO_maps.m overrides this per map.
end
switch metabolic_map
    case 'ox_gluc', file_tag = 'oxglc';
    case 'ox',      file_tag = 'ox';
    case 'glc',     file_tag = 'glc';
    otherwise, error('metabolic_map must be ''ox_gluc'', ''ox'', or ''glc''.');
end
% Append a ceff-source suffix when non-legacy so different Ceff sources
% don't overwrite each other's results. 'legacy' keeps the original
% filename, preserving byte-identical paths for prior runs.
switch ceff_source
    case 'legacy',          ceff_tag = '';
    case 'gec_competitive', ceff_tag = '_gecCompetitive';
    case 'gec_cooperative', ceff_tag = '_gecCooperative';
end
results_filename = sprintf('results_NEMO_%s_DBS80%s.mat', file_tag, ceff_tag);

%% ===== Step 1: load HCP rest DBS80, compute FCemp / COVtauemp =====
fprintf('[1/8] loading HCP_rest_DBS80_Functional_Timeseries.mat (`bigmat`: %d×T×NSUB)...\n', N);
S = load(fullfile(projDir, 'HCP_rest_DBS80_Functional_Timeseries.mat'));   % `bigmat`, `TR`, etc.
bigmat = S.bigmat;                                  % 80 × T × NSUB
if isfield(S, 'TR'),  TR = S.TR;  end               % use stored TR if present
[nROI, T, NSUB] = size(bigmat);
assert(nROI >= max(indexN), 'bigmat ROI dim (%d) < requested N=%d', nROI, N);
fprintf('     bigmat = %d × %d × %d   (using %d ROIs)\n', nROI, T, NSUB, N);

% load precomputed empirical f_diff (1×80) so we don't redo the spectral peak
Eemp   = load(fullfile(projDir, 'empirical_HCP_rest_DBS80.mat'));
f_diff = Eemp.f_diff(:).';
assert(numel(f_diff) >= max(indexN), 'f_diff length (%d) < N=%d', numel(f_diff), N);
f_diff = f_diff(indexN);
fprintf('     loaded precomputed f_diff (mean %.4f Hz, range [%.4f %.4f])\n', ...
        mean(f_diff), min(f_diff), max(f_diff));

% running accumulators + per-subject FC stash (kept so we can bootstrap later)
FCemp_sub = zeros(N, N, NSUB, 'single');     % N×N×NSUB single (~12 MB at N=80)
keep      = false(1, NSUB);
COVtauemp = zeros(N, N);
n_used    = 0;
for nsub = 1:NSUB
    ts = double(squeeze(bigmat(indexN, :, nsub)));    % N × T
    if any(~isfinite(ts(:))),  continue;  end
    Tsub = size(ts, 2);

    % detrend + bandpass per region
    sig = zeros(N, Tsub);
    for r = 1:N
        x = detrend(ts(r,:) - nanmean(ts(r,:)));
        sig(r,:) = filtfilt(bfilt, afilt, x);
    end
    sig  = sig(:, 100:end-100);                       % drop filter edge transients
    Tred = size(sig, 2);

    % per-subject FC(0); kept for bootstrap
    FCsub                 = corrcoef(sig.');
    FCemp_sub(:, :, nsub) = single(FCsub);
    keep(nsub)            = true;

    % COV(τ): vectorised cross-covariance at lag +Tau, normalised by std·std'
    Sm      = sig - mean(sig, 2);
    cov_tau = (Sm(:, 1+Tau:end) * Sm(:, 1:end-Tau).') / Tred;
    sd      = std(sig, 0, 2);
    COVtauemp = COVtauemp + cov_tau ./ (sd * sd.');

    n_used = n_used + 1;
    if mod(nsub,50)==0 || nsub==NSUB, fprintf('     subj %d/%d\n', nsub, NSUB); end
end
FCemp_sub = FCemp_sub(:, :, keep);                    % N × N × n_used
FCemp     = double(mean(FCemp_sub, 3));               % group-mean target
COVtauemp = COVtauemp / n_used;
fprintf('     FCemp/COVtauemp ready (averaged over %d subjects, per-subject FC kept for bootstrap).\n', n_used);

%% ===== Step 2: load SC, optionally fit Ceff (RS-GEC) =====
fprintf('[2/8] structural connectivity / Ceff...\n');
SCl = load(fullfile(projDir, 'SC_dbs80HARDIFULL.mat'));
if isfield(SCl, 'SC_dbs80HARDI')
    SC = double(SCl.SC_dbs80HARDI);
elseif isfield(SCl, 'SC_dbs80FULL')
    SC = double(SCl.SC_dbs80FULL);
else
    fn = fieldnames(SCl);  SC = double(SCl.(fn{1}));
end
SC = SC(indexN, indexN);
fprintf('     loaded SC (%dx%d)\n', size(SC,1), size(SC,2));

C = SC;  C(C<0) = 0;
if max(C(:)) > 0,  C = C / max(C(:));  end
maxC  = 0.2;
C     = C * maxC;
maskC = C > 0;

switch ceff_source
case 'legacy'
    % ---------- ORIGINAL NEMO cooperative-only RS-GEC (unchanged) ----------
    if do_gec_fit
        epsFC = 0.0004;  epsFCtau = 0.0001;
        Cnew  = C;  olderror = 1e5;
        for iter = 1:max_gec_iter
            [FCsim, COVsim, COVsimtotal, A] = hopf_int(Cnew, f_diff, beta);
            COVtausim = expm((Tau*TR)*A) * COVsimtotal;
            COVtausim = COVtausim(1:N, 1:N);
            dv          = sqrt(diag(COVsim));
            sigratiosim = 1 ./ (dv * dv.');
            COVtausim   = COVtausim .* sigratiosim;
            if mod(iter,100) < 0.1
                errornow = mean(mean((FCemp-FCsim).^2)) + mean(mean((COVtauemp-COVtausim).^2));
                fprintf('     gec iter %4d  err=%.5f\n', iter, errornow);
                if (olderror-errornow)/max(errornow,eps) < 1e-4 || olderror < errornow, break; end
                olderror = errornow;
            end
            dFC  = (FCemp     - FCsim);
            dCOV = (COVtauemp - COVtausim);
            Cnew = Cnew + (epsFC*dFC + epsFCtau*dCOV) .* maskC;
            Cnew(Cnew < 0) = 0;
            Cnew = Cnew / max(Cnew(:)) * maxC;
        end
        Ceff = Cnew;
    else
        Ceff = C;
        fprintf('     do_gec_fit=false; using scaled raw SC as Ceff.\n');
    end

case {'gec_competitive', 'gec_cooperative'}
    % ---------- Luppi et al 2026 GEC (allows competitive / negative weights) ----------
    % Uses NEMO/competitive_hopf/matlab_utils/fcn_Hopf_NonLinear_cooperative_competitive.m
    % via the NEMO-local wrapper fit_competitive_gec_NEMO.m.
    %
    % The fit is per-subject (the upstream function expects a single N x T TS),
    % so by default we fit on subject 1 only (~minutes at N=80). Switch to
    % opts.subject_idx = 1:NSUB (or any subset) to fit per-subject and average
    % at the cost of linearly more compute. Result is cached on disk.
    fprintf('     ceff_source = ''%s'' -> calling fit_competitive_gec_NEMO\n', ceff_source);
    is_coop = strcmp(ceff_source, 'gec_cooperative');
    gec_opts = struct( ...
        'flag_cooperativeOnly_YN', is_coop, ...
        'subject_idx',             1, ...
        'TR',                      TR, ...
        'filter_low',              flp, ...
        'filter_high',             fhi);
    [Ceff, gec_info] = fit_competitive_gec_NEMO(bigmat, SC, indexN, gec_opts);
    fprintf('     Luppi GEC fit done in %.1fs   FC_corr(mean)=%.4f   GEC range [%.4f, %.4f]\n', ...
        gec_info.runtime_sec, gec_info.FC_corr_mean, min(Ceff(:)), max(Ceff(:)));
    if any(Ceff(:) < 0)
        fprintf('     Ceff has %d negative entries (%.1f%% of off-diagonal)\n', ...
            nnz(Ceff < 0), 100*nnz(Ceff<0)/(N*(N-1)));
    end
end

% baseline (homogeneous β) reference simulation
FCsim_base = hopf_int_pert_beta(Ceff, f_diff, beta*ones(1,N));
maxfcfitt  = ssim(single(FCemp), single(FCsim_base));
fprintf('     baseline (β=%.4f) SSIM(FCemp,FCsim) = %.4f\n', beta, maxfcfitt);

%% ===== Step 3: load metabolic map (per-region β scalar) =====
% The variable is named `ox_gluc_ratio` for historical reasons but holds
% whichever map was selected by `metabolic_map` above (ox_gluc / ox / glc).
fprintf('[3/8] loading metabolic map "%s" (DBS80)...\n', metabolic_map);
switch metabolic_map
case 'ox_gluc'
    % --- DEFAULT (unchanged) ---
    ratio_path = fullfile(projDir, 'ox_gluc_ratio.mat');
    if isfile(ratio_path)
        Rload = load(ratio_path);
        if isfield(Rload, 'ox_gluc_ratio')
            ox_gluc_ratio = Rload.ox_gluc_ratio(:).';
        else
            fn = fieldnames(Rload);  ox_gluc_ratio = Rload.(fn{1})(:).';
        end
        fprintf('     loaded ox_gluc_ratio.mat\n');
    else
        warning(['ox_gluc_ratio.mat not found. Falling back to ratio_values ' ...
                 'in CMRO2_over_CMRglc_DBS80.mat normalised to mean 1. ' ...
                 'Drop your DBS80 ox_gluc_ratio.mat into MemBrain to override.']);
        Rload = load(fullfile(projDir, 'CMRO2_over_CMRglc_DBS80.mat'));
        rv    = Rload.ratio_values(:);
        ox_gluc_ratio = (rv ./ mean(rv,'omitnan')).';
    end

case {'ox', 'glc'}
    % CMRO2 or CMRglc alone, normalised to mean 1 (same scale as the
    % ox_gluc_ratio default → PSO bounds & seed don't need to change).
    map_path = fullfile(projDir, 'metabolic_maps_DBS80.mat');
    if ~isfile(map_path)
        error('metabolic_maps_DBS80.mat not found in %s', projDir);
    end
    Mload = load(map_path);
    if strcmp(metabolic_map, 'ox')
        if ~isfield(Mload, 'cmro2_values')
            error('metabolic_maps_DBS80.mat is missing field `cmro2_values`.');
        end
        rv = Mload.cmro2_values(:);
        fprintf('     loaded cmro2_values from metabolic_maps_DBS80.mat\n');
    else  % 'glc'
        if ~isfield(Mload, 'cmrglc_values')
            error('metabolic_maps_DBS80.mat is missing field `cmrglc_values`.');
        end
        rv = Mload.cmrglc_values(:);
        fprintf('     loaded cmrglc_values from metabolic_maps_DBS80.mat\n');
    end
    ox_gluc_ratio = (rv ./ mean(rv, 'omitnan')).';
end
% ensure it lines up with our (possibly subset) indexN
if numel(ox_gluc_ratio) ~= N
    if numel(ox_gluc_ratio) >= max(indexN)
        ox_gluc_ratio = ox_gluc_ratio(indexN);
    else
        error('ox_gluc_ratio length (%d) does not match N=%d', numel(ox_gluc_ratio), N);
    end
end
% guard rail: any NaN region (e.g. subcortex with no metabolic coverage)
% gets a homogeneous-noise weight of 1 so it doesn't poison betax.
nan_mask = ~isfinite(ox_gluc_ratio);
if any(nan_mask)
    warning('NEMO:oxglcNaN', ...
        '%d/%d ox_gluc_ratio entries are NaN/Inf (regions [%s]); filling with 1.0', ...
        sum(nan_mask), numel(ox_gluc_ratio), num2str(find(nan_mask)));
    ox_gluc_ratio(nan_mask) = 1;
end
fprintf('     ox_gluc_ratio: range [%.3f, %.3f], mean %.3f\n', ...
        min(ox_gluc_ratio), max(ox_gluc_ratio), mean(ox_gluc_ratio));

%% ===== Step 4 (simplified): scale β per region by ox_gluc_ratio =====
% --- ORIGINAL NEMO step 4 (REMOVED) ---
%   nvars1 = 2 + Nnr;
%   betax = beta * (BB+AA*neuroreceptormaps);
% --- SIMPLIFIED replacement: ---

% Linear noise model fit by PSO:
%     betax_i = beta * ( BB + AA * ox_gluc_ratio_i )
% The optimiser searches (AA, BB); `beta` enters ONCE, here and inside
% the cost function (FC_prediction_Hopf_Task_ox), so the values of AA/BB
% returned are on the same scale used downstream.

nvars1   = 2;                                    % [BB, AA]
xinit1   = [0 1];                                % BB0 = 0, AA0 = 1
lb1      = [ 0  0];                              % betax_i = beta*(BB+AA*ratio_i) must be >= 0;
ub1      = [ 5  5];                              % since ratio_i >= 0, BB,AA >= 0 is sufficient

% custom initial swarm — clamp to bounds so particleswarm doesn't drop seeds
swarmSize = 40;
initpop1  = 0.01*randn(swarmSize, nvars1) + repmat(xinit1, swarmSize, 1);
initpop1  = max(min(initpop1, ub1), lb1);

options = optimoptions('particleswarm', ...
    'InitialSwarmMatrix', initpop1, ...           % seeded around (BB,AA) = (0,1)
    'SwarmSize',          swarmSize, ...          % match initpop1 height (silences warning)
    'Display',            'iter', ...
    'MaxIterations',      300);

% closure: pass everything the cost function needs explicitly (no globals).
costFn = @(p) FC_prediction_Hopf_Task_ox(p, Ceff, FCemp, f_diff, ox_gluc_ratio, beta);

fprintf('[4/8] PSO fit of (BB, AA) on full-sample FCemp (SwarmSize=%d, MaxIter=300)...\n', swarmSize);
[par, fval] = particleswarm(costFn, nvars1, lb1, ub1, options);
BB = par(1);
AA = par(2);
betax = beta * (BB + AA*ox_gluc_ratio);          % 1 × N per-region noise amplitude
fprintf('     PSO fit: BB = %+.4f, AA = %+.4f   (fval = %.4f, SSIM = %.4f)\n', ...
        BB, AA, fval, 1 - fval);


%% ===== Step 5: simulate + compare =====
fprintf('[5/8] simulating Hopf with metabolically-scaled β...\n');
FCsim_oxglc = hopf_int_pert_beta(Ceff, f_diff, betax);
ssim_oxglc  = ssim(single(FCemp), single(FCsim_oxglc));
fprintf('     metabolic SSIM(FCemp,FCsim) = %.4f   (fixed-β=%.3f baseline = %.4f)\n', ...
        ssim_oxglc, beta, maxfcfitt);

%% ===== Step 6: spatial-permutation null over ox_gluc_ratio =====
%  WHY NOT A HOMOGENEOUS β SWEEP?
%  Under the analytical Hopf solver, FCth = corrcov(Cvth). With a single
%  scalar β, the noise covariance is Q = β²·I, so Cvth ∝ β² everywhere and
%  corrcov divides the scalar out — the correlation matrix (and therefore
%  the SSIM) is *invariant to a uniform β*. A 1D β-sweep would return a
%  flat curve and is meaningless as a control here.
%
%  The right control is a SPATIAL PERMUTATION NULL: keep the same set of
%  ox_gluc_ratio values, shuffle which region gets which one, recompute
%  SSIM. This isolates the contribution of the spatial *pattern* of
%  metabolism (heterogeneous AND aligned to brain anatomy) from the
%  contribution of merely having heterogeneous β values at all.
%  FAIRNESS: each shuffle gets its own (AA, BB) PSO fit, in the same
%  linear model class as the metabolic condition:
%        betax_perm = beta * ( BB_perm + AA_perm * shuffled_ratio )
%  Without this, the metabolic condition gets 2 free fit parameters the
%  null is denied — and the comparison is no longer apples-to-apples.
%  Reduced PSO budget per shuffle (SwarmSize=20, MaxIter=100) seeded
%  near the full-sample (BB, AA) — fast convergence on a 2-D surface.
n_perm = 200;
rng(0, 'twister');                              % reproducible permutations
ssim_perm      = zeros(1, n_perm);
AA_perm        = zeros(1, n_perm);
BB_perm        = zeros(1, n_perm);
FCsim_perm_all = zeros(N, N, n_perm, 'single');  % stash for step 7

swarmSize_perm = 20;
opts_perm = optimoptions('particleswarm', ...
    'SwarmSize',     swarmSize_perm, ...
    'MaxIterations', 100, ...
    'Display',       'off');

fprintf('[6/8] permutation null over ox_gluc_ratio (%d perms, PSO refit per shuffle)...\n', n_perm);
t_perm0 = tic;
for k = 1:n_perm
    perm_ratio = ox_gluc_ratio(randperm(N));

    % seed swarm near full-sample (BB, AA) and clamp to bounds
    initpop_k = 0.05*randn(swarmSize_perm, nvars1) + repmat([BB AA], swarmSize_perm, 1);
    initpop_k = max(min(initpop_k, ub1), lb1);
    opts_k    = optimoptions(opts_perm, 'InitialSwarmMatrix', initpop_k);

    cost_k = @(p) FC_prediction_Hopf_Task_ox(p, Ceff, FCemp, f_diff, perm_ratio, beta);
    [par_k, fval_k] = particleswarm(cost_k, nvars1, lb1, ub1, opts_k);

    BB_perm(k) = par_k(1);
    AA_perm(k) = par_k(2);
    betax_k    = beta * (BB_perm(k) + AA_perm(k) * perm_ratio);
    FCsim_perm_all(:,:,k) = single(hopf_int_pert_beta(Ceff, f_diff, betax_k));
    ssim_perm(k) = 1 - fval_k;       % == ssim(FCemp, FCsim_perm_all(:,:,k))

    if mod(k, 25) == 0
        fprintf('     perm %3d/%d   BB=%+.3f  AA=%+.3f  SSIM=%.4f  (elapsed %.1f min)\n', ...
                k, n_perm, BB_perm(k), AA_perm(k), ssim_perm(k), toc(t_perm0)/60);
    end
end
fprintf('     step 6 done in %.1f min\n', toc(t_perm0)/60);
% one-sided p: P(SSIM_perm >= SSIM_metabolic)
p_perm        = (1 + sum(ssim_perm >= ssim_oxglc)) / (1 + n_perm);
ssim_perm_mu  = mean(ssim_perm);
ssim_perm_sd  = std(ssim_perm);
ssim_perm_max = max(ssim_perm);
ssim_perm_min = min(ssim_perm);

fprintf('\n     ===== fit summary =====\n');
fprintf('     fixed-β baseline (β=%.4f homogeneous)            SSIM = %.4f   [invariant to β in this pipeline]\n', beta, maxfcfitt);
fprintf('     metabolic per-region (β * ox_gluc_ratio)         SSIM = %.4f\n', ssim_oxglc);
fprintf('     permutation null (n=%d, shuffled metabolism)     SSIM = %.4f ± %.4f   [%.4f .. %.4f]\n', ...
        n_perm, ssim_perm_mu, ssim_perm_sd, ssim_perm_min, ssim_perm_max);
fprintf('     Δ(metabolic - null mean)  = %+.4f   (p = %.4f, %d-perm one-sided)\n\n', ...
        ssim_oxglc - ssim_perm_mu, p_perm, n_perm);

%% ===== Step 6b: spatial-autocorrelation-preserving null (OXYGEN ONLY) =====
%  Step 6 shuffles the map UNIFORMLY (randperm), which destroys its spatial
%  smoothness — an easy null. The stronger control keeps the smoothness and
%  only breaks the anatomical alignment. We generate spatially-autocorrelated
%  surrogate maps by Moran spectral randomisation (random sign-flip of the
%  spatial-spectral coefficients) on a spatial weight matrix W, then rank-remap
%  each surrogate onto the empirical values — so every surrogate is a
%  spatially-structured RESHUFFLE of the same numbers. Each surrogate is refit
%  with its own (BB, AA), exactly like Step 6. Run for oxygen only, since that
%  is the result that matters.
%
%  W needs a region geometry. If DBS80_geometry.mat (with `coords` N×3 centroids
%  or a precomputed N×N `D`) is in MemBrain we use Gaussian distance weights;
%  otherwise we fall back to the structural connectivity SC as a proxy adjacency
%  (weaker — flagged — drop a geometry file before publication).
do_spin = strcmp(metabolic_map, 'ox');
if do_spin
    n_spin = 200;
    rng(2, 'twister');

    % --- spatial weight matrix W (N×N, symmetric, zero diagonal) ---
    geom_file = fullfile(projDir, 'DBS80_geometry.mat');
    if isfile(geom_file)
        G = load(geom_file);
        if     isfield(G, 'D'),      Dm = double(G.D);
        elseif isfield(G, 'coords'), Dm = squareform(pdist(double(G.coords)));
        else,  fn = fieldnames(G);   Dm = double(G.(fn{1}));
        end
        if size(Dm,1) > N, Dm = Dm(indexN, indexN); end
        d0 = median(Dm(~eye(N) & isfinite(Dm)));
        W  = exp(-(Dm.^2) / (2*d0^2));  W(1:N+1:end) = 0;
        fprintf('[6b/8] spatial-autocorr null: Gaussian W from geometry file\n');
    else
        warning('NEMO:noGeometry', ...
            ['No DBS80_geometry.mat — using SC as a PROXY spatial weight. ' ...
             'Drop a DBS80_geometry.mat (coords N×3 or D N×N) into %s for publication.'], projDir);
        W = double(SC);  W = 0.5*(W + W.');  W(W < 0) = 0;  W(1:N+1:end) = 0;
        fprintf('[6b/8] spatial-autocorr null: PROXY W from SC\n');
    end

    % --- Moran eigenvectors of the doubly-centred weight matrix (once) ---
    Hc = eye(N) - ones(N)/N;  Bc = Hc * W * Hc;  Bc = 0.5*(Bc + Bc.');
    [Vm, Lm] = eig(Bc);  [~, od] = sort(diag(Lm), 'descend');  Vm = Vm(:, od);
    sorted_vals = sort(ox_gluc_ratio(:));
    coef0       = Vm.' * (ox_gluc_ratio(:) - mean(ox_gluc_ratio));

    ssim_spin = zeros(1, n_spin);
    AA_spin   = zeros(1, n_spin);
    BB_spin   = zeros(1, n_spin);
    opts_spin = optimoptions('particleswarm', ...
        'SwarmSize', 20, 'MaxIterations', 100, 'Display', 'off');

    fprintf('[6b/8] spatial-autocorr null (%d surrogates, PSO refit each)...\n', n_spin);
    t_spin0 = tic;
    for k = 1:n_spin
        sgn  = sign(randn(N, 1));  sgn(sgn == 0) = 1;
        surr = Vm * (sgn .* coef0);                  % spatially-autocorrelated field
        [~, rk]         = sort(surr);
        spin_ratio      = zeros(1, N);
        spin_ratio(rk)  = sorted_vals;               % rank-remap onto empirical values

        initpop_k = max(min(0.05*randn(20, nvars1) + repmat([BB AA], 20, 1), ub1), lb1);
        opts_k    = optimoptions(opts_spin, 'InitialSwarmMatrix', initpop_k);
        cost_k    = @(p) FC_prediction_Hopf_Task_ox(p, Ceff, FCemp, f_diff, spin_ratio, beta);
        [par_k, fval_k] = particleswarm(cost_k, nvars1, lb1, ub1, opts_k);

        BB_spin(k)   = par_k(1);
        AA_spin(k)   = par_k(2);
        ssim_spin(k) = 1 - fval_k;
        if mod(k, 25) == 0
            fprintf('     spin %3d/%d   SSIM=%.4f  (elapsed %.1f min)\n', ...
                    k, n_spin, ssim_spin(k), toc(t_spin0)/60);
        end
    end
    p_spin = (1 + sum(ssim_spin >= ssim_oxglc)) / (1 + n_spin);
    fprintf('     spatial-autocorr null  SSIM = %.4f ± %.4f   metabolic = %.4f   (p = %.4f, n=%d)\n', ...
            mean(ssim_spin), std(ssim_spin), ssim_oxglc, p_spin, n_spin);
    fprintf('     [contrast]  uniform-shuffle p = %.4f   spatial-autocorr p = %.4f\n\n', p_perm, p_spin);
else
    % not oxygen — leave placeholders so the save block always succeeds
    ssim_spin = [];  AA_spin = [];  BB_spin = [];  p_spin = NaN;  n_spin = 0;
end

%% ===== Step 7: bootstrap over subjects (error bars on every condition) =====
%  All three conditions (homogeneous β, metabolic β·ox_gluc_ratio, perm null)
%  use the same analytical Hopf solver, so the only source of variability we
%  can give them is in the SSIM TARGET — i.e., the empirical FC. We resample
%  subjects with replacement (n_used at a time) and average their per-subject
%  FCs to get FCemp_b, then evaluate SSIM under each β specification.
%  Ceff is held fixed (fit once on the full sample) so error bars reflect
%  subject sampling variability, not GEC convergence variability.
n_boot = 200;
rng(1, 'twister');                                    % reproducible
ssim_boot_homo = zeros(1, n_boot);
ssim_boot_meta = zeros(1, n_boot);
ssim_boot_perm = zeros(1, n_boot);
% precompute the homogeneous and metabolic FCsim once — they don't change
FCsim_homo = hopf_int_pert_beta(Ceff, f_diff, beta * ones(1, N));
FCsim_meta = FCsim_oxglc;

% pair each bootstrap iter to one of the fitted-perm FCsims from step 6
% (paired if n_boot==n_perm; otherwise sample with replacement). This
% mirrors how the metabolic condition uses one FIXED fitted FCsim_meta.
if n_boot == n_perm
    perm_idx_boot = 1:n_boot;
else
    perm_idx_boot = randi(n_perm, 1, n_boot);
end

fprintf('[7/8] bootstrap over %d subjects (n_boot=%d)...\n', n_used, n_boot);
for b = 1:n_boot
    idx     = randi(n_used, 1, n_used);                   % sample with replacement
    FCemp_b = double(mean(FCemp_sub(:, :, idx), 3));      % bootstrap empirical target
    ssim_boot_homo(b) = ssim(single(FCemp_b), single(FCsim_homo));
    ssim_boot_meta(b) = ssim(single(FCemp_b), single(FCsim_meta));
    % perm null: use a step-6 PSO-fit shuffle (same model class as metabolic)
    ssim_boot_perm(b) = ssim(single(FCemp_b), FCsim_perm_all(:,:,perm_idx_boot(b)));
    if mod(b, 25) == 0, fprintf('     boot %d/%d\n', b, n_boot); end
end

% bootstrap-aware p-value: fraction of (boot perm SSIM ≥ boot metabolic SSIM)
p_boot = (1 + sum(ssim_boot_perm >= ssim_boot_meta)) / (1 + n_boot);

fprintf('\n     ===== bootstrap summary (n_boot=%d) =====\n', n_boot);
fprintf('     homogeneous β   :  %.4f ± %.4f\n', mean(ssim_boot_homo), std(ssim_boot_homo));
fprintf('     metabolic β·ox  :  %.4f ± %.4f\n', mean(ssim_boot_meta), std(ssim_boot_meta));
fprintf('     perm-null β     :  %.4f ± %.4f   (paired p = %.4f)\n\n', ...
        mean(ssim_boot_perm), std(ssim_boot_perm), p_boot);

%% ===== save =====
fprintf('[8/8] saving...\n');
save(fullfile(projDir, results_filename), ...
     'betax', 'AA', 'BB', 'fval', 'ox_gluc_ratio', 'metabolic_map', ...
     'ceff_source', 'Ceff', ...
     'FCemp', 'FCsim_base', 'FCsim_oxglc', 'FCsim_homo', ...
     'ssim_oxglc', 'maxfcfitt', ...
     'ssim_perm', 'AA_perm', 'BB_perm', 'FCsim_perm_all', ...
     'p_perm', 'n_perm', 'ssim_perm_mu', 'ssim_perm_sd', ...
     'ssim_spin', 'AA_spin', 'BB_spin', 'p_spin', 'n_spin', ...
     'ssim_boot_homo', 'ssim_boot_meta', 'ssim_boot_perm', 'perm_idx_boot', ...
     'n_boot', 'p_boot', 'n_used', ...
     'f_diff', 'beta', 'N', 'indexN', '-v7.3');
fprintf('     wrote %s/%s\n', projDir, results_filename);

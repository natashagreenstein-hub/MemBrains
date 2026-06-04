%% plot_NEMO_oxglc_results.m
%
%  Figure 1 : CMRO2/CMRglc ratio        (results_NEMO_oxglc_DBS80.mat)
%  Figure 2 : CMRO2 only + CMRglc only  (results_NEMO_ox_DBS80.mat, glc)
%  Figure 3 : OXYGEN spatial-autocorrelation control — adds the
%             spatial-autocorr-preserving null next to the metabolic fit,
%             the homogeneous baseline and the uniform-shuffle null.
%             Drawn only if results_NEMO_ox_DBS80.mat carries `ssim_spin`.
%
%  SSIM panels show bootstrap mean ± sd. The spatial-autocorr null is shown
%  as the across-surrogate mean ± sd of ssim_spin. All numeric results
%  (AA, BB, dot values, p-values) are printed to the command window.
%
%  Missing files / fields are handled gracefully (placeholder / skip).

clearvars -except ceff_source; close all;
projDir = '/Users/natashag/Documents/Claude/Projects/MemBrain';

% Which ceff_source run to plot — must match what NEMO_oxglc wrote.
if ~exist('ceff_source', 'var') || isempty(ceff_source)
    ceff_source = 'legacy';
end
switch ceff_source
    case 'legacy',          ceff_tag = '';
    case 'gec_competitive', ceff_tag = '_gecCompetitive';
    case 'gec_cooperative', ceff_tag = '_gecCooperative';
    otherwise, error('unknown ceff_source ''%s''', ceff_source);
end

runs = struct( ...
    'tag',       {'oxglc',          'ox',            'glc'}, ...
    'file',      {sprintf('results_NEMO_oxglc_DBS80%s.mat', ceff_tag), ...
                  sprintf('results_NEMO_ox_DBS80%s.mat',    ceff_tag), ...
                  sprintf('results_NEMO_glc_DBS80%s.mat',   ceff_tag)}, ...
    'map_label', {'CMRO_2\CMRglc\_ratio', 'CMRO_2',        'CMRglc'}, ...
    'panel_ttl', {'CMRO_2 / CMRglc ratio',  'CMRO_2 only',   'CMRglc only'});

%% ===== Figure 1: ratio =====
fig1 = figure('Color', 'w', 'Position', [80 80 700 520]);
plot_one_panel(projDir, runs(1));
f_path1 = fullfile(projDir, sprintf('fig_NEMO_oxglc_SSIM%s.png', ceff_tag));
exportgraphics(fig1, f_path1, 'Resolution', 200);
fprintf('wrote %s\n', f_path1);

%% ===== Figure 2: CMRO2 only + CMRglc only =====
fig2 = figure('Color', 'w', 'Position', [80 80 1300 520]);
subplot(1, 2, 1); plot_one_panel(projDir, runs(2));
subplot(1, 2, 2); plot_one_panel(projDir, runs(3));
f_path2 = fullfile(projDir, sprintf('fig_NEMO_ox_glc_SSIM%s.png', ceff_tag));
exportgraphics(fig2, f_path2, 'Resolution', 200);
fprintf('wrote %s\n', f_path2);

%% ===== Figure 3: oxygen spatial-autocorrelation control =====
if has_spin(projDir, runs(2))
    fig3 = figure('Color', 'w', 'Position', [80 80 760 520]);
    plot_spin_panel(projDir, runs(2));
    f_path3 = fullfile(projDir, sprintf('fig_NEMO_ox_spatialnull%s.png', ceff_tag));
    exportgraphics(fig3, f_path3, 'Resolution', 200);
    fprintf('wrote %s\n', f_path3);
else
    fprintf('\n[Figure 3] no ssim_spin in %s — run NEMO with metabolic_map=''ox'' to generate it.\n', runs(2).file);
end

%% ---------- SSIM per-panel renderer (homogeneous / metabolic / uniform null) ----------
function plot_one_panel(projDir, run)
    f = fullfile(projDir, run.file);
    if ~isfile(f)
        text(0.5, 0.5, sprintf('%s\nnot found —\nrun NEMO with metabolic\\_map = ''%s''', ...
             run.file, run.tag), ...
             'HorizontalAlignment','center','Units','normalized','FontSize',11);
        title(run.panel_ttl); axis off;
        fprintf('\n[%s] %s not found — skipped.\n', run.panel_ttl, run.file);
        return;
    end
    R = load(f);
    vals = [mean(R.ssim_boot_homo), mean(R.ssim_boot_meta), mean(R.ssim_boot_perm)];
    errs = [std(R.ssim_boot_homo),  std(R.ssim_boot_meta),  std(R.ssim_boot_perm)];
    labs = {'homogeneous β', 'metabolic β', 'permutation null'};
    cols = [0.20 0.45 0.70; 0.85 0.10 0.10; 0.50 0.50 0.50];
    draw_points(vals, errs, labs, cols);
    title(run.panel_ttl);

    fprintf('\n================ %s ================\n', run.panel_ttl);
    fprintf('file: %s\n', run.file);
    fprintf('  homogeneous β    : %.4f ± %.4f\n', vals(1), errs(1));
    fprintf('  metabolic β      : %.4f ± %.4f\n', vals(2), errs(2));
    fprintf('  permutation null : %.4f ± %.4f\n', vals(3), errs(3));
    if isfield(R,'AA') && isfield(R,'BB') && isfield(R,'beta')
        fprintf('fit: beta = %.3f , BB = %.4f , AA = %.4f\n', R.beta, R.BB, R.AA);
    end
    if isfield(R, 'p_boot') && isfield(R, 'n_boot')
        fprintf('paired p (permutation null) = %.4f  (n_boot = %d)\n', R.p_boot, R.n_boot);
    end
end

%% ---------- Figure 3 renderer: oxygen + spatial-autocorrelation null ----------
function plot_spin_panel(projDir, run)
    f = fullfile(projDir, run.file);
    R = load(f);
    vals = [mean(R.ssim_boot_homo), mean(R.ssim_boot_meta), mean(R.ssim_boot_perm), mean(R.ssim_spin)];
    errs = [std(R.ssim_boot_homo),  std(R.ssim_boot_meta),  std(R.ssim_boot_perm),  std(R.ssim_spin)];
    labs = {'homogeneous β', 'metabolic β', 'uniform-shuffle null', 'spatial-autocorr null'};
    cols = [0.20 0.45 0.70; 0.85 0.10 0.10; 0.50 0.50 0.50; 0.95 0.55 0.10];
    draw_points(vals, errs, labs, cols);
    title(sprintf(run.panel_ttl));

    fprintf('\n================ %s : spatial-autocorr null ================\n', run.panel_ttl);
    fprintf('  metabolic β           : %.4f\n', mean(R.ssim_boot_meta));
    fprintf('  uniform-shuffle null  : %.4f ± %.4f   (one-sided p = %.4f)\n', ...
            mean(R.ssim_boot_perm), std(R.ssim_boot_perm), R.p_perm);
    fprintf('  spatial-autocorr null : %.4f ± %.4f   (one-sided p = %.4f, n=%d surrogates)\n', ...
            mean(R.ssim_spin), std(R.ssim_spin), R.p_spin, R.n_spin);
end

%% ---------- shared point/errorbar drawer ----------
function draw_points(vals, errs, labs, cols)
    nC = numel(vals);
    hold on;
    for k = 1:nC
        errorbar(k, vals(k), errs(k), 'k', 'LineStyle','none', 'LineWidth',1.4, 'CapSize',14);
    end
    for k = 1:nC
        plot(k, vals(k), 'o', 'MarkerFaceColor', cols(k,:), ...
             'MarkerEdgeColor','k', 'MarkerSize',14, 'LineWidth',1.0);
    end
    xlim([0.5 nC+0.5]);
    xticks(1:nC); xticklabels(labs); xtickangle(20);
    ylabel('SSIM(FC_{emp}, FC_{sim})');
    set(gca, 'FontSize', 11, 'TickDir', 'out');
    ylim([min(vals - errs) - 0.02, max(vals + errs) + 0.025]);
    grid on; box off;
end

%% ---------- helper: does this file carry a spatial-autocorr null? ----------
function tf = has_spin(projDir, run)
    f = fullfile(projDir, run.file);
    tf = false;
    if ~isfile(f), return; end
    w = whos('-file', f);
    if ~any(strcmp({w.name}, 'ssim_spin')), return; end
    S = load(f, 'ssim_spin');
    tf = ~isempty(S.ssim_spin);
end

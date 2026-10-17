%% ========================================================================
% FIGURE 3B — FINAL
% Planned alpha-band decomposition contrasts
%
% Paper 3: PD EEG Rest-to-Walking Network Reconfiguration
%
% Frozen STEP10S planned contrasts:
%
%   C1 = Interhemispheric - mean(Intra-left, Intra-right)
%        Estimate = -0.0038432
%        95% CI   = [-0.018430, 0.010744]
%        Raw p    = 0.60432
%        Holm p   = 0.60432
%
%   C2 = Intra-left - Intra-right
%        Estimate = -0.025702
%        95% CI   = [-0.042546, -0.008859]
%        Raw p    = 0.0029194
%        Holm p   = 0.0058387
%
% Interpretation:
%   These are difference-of-differences contrasts in the PD-HC
%   NetworkType effects.
%
% IMPORTANT:
%   Do NOT interpret C2 as showing that either Intra-left or Intra-right
%   individually has a significant HC-PD effect.
%
% Figure style:
%   Filled circle = CI excludes zero
%   Open circle   = CI includes zero
%
% MATLAB R2023b compatible
%% ========================================================================

clear;
clc;
close all;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% ========================================================================
% 1. OUTPUT DIRECTORY
%% ========================================================================

outDir = ...
    fullfile(cfg.datasetRoot,'PAPER3_PRIMARY_NBS','FIGURES');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ========================================================================
% 2. FROZEN STEP10S RESULTS
%% ========================================================================

% Contrast 1:
% Interhemispheric - mean(Intra-left, Intra-right)

C1_estimate = -0.0038432;
C1_CI_low   = -0.018430;
C1_CI_high  =  0.010744;
C1_rawP     =  0.60432;
C1_holmP    =  0.60432;

% Contrast 2:
% Intra-left - Intra-right

C2_estimate = -0.025702;
C2_CI_low   = -0.042546;
C2_CI_high  = -0.008859;
C2_rawP     =  0.0029194;
C2_holmP    =  0.0058387;

%% ========================================================================
% 3. ARRANGE FOR DISPLAY
%
% Put significant left-right contrast on top.
%% ========================================================================

plotEstimate = [ ...
    C2_estimate, ...
    C1_estimate];

plotLow = [ ...
    C2_CI_low, ...
    C1_CI_low];

plotHigh = [ ...
    C2_CI_high, ...
    C1_CI_high];

plotHolmP = [ ...
    C2_holmP, ...
    C1_holmP];

plotRawP = [ ...
    C2_rawP, ...
    C1_rawP];

plotNames = { ...
    'Intra-left - Intra-right', ...
    'Inter - mean(Intra-left, Intra-right)'};

% Row positions
y = [1 2];

%% ========================================================================
% 4. COLORS
%% ========================================================================

darkGray = [0.26 0.26 0.26];
axisGray = [0.12 0.12 0.12];
zeroGray = [0.72 0.72 0.72];

%% ========================================================================
% 5. CREATE FIGURE
%% ========================================================================

fig = figure( ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [3 3 13.0 7.5]);

ax = axes(fig);

hold(ax, 'on');

%% ========================================================================
% 6. ZERO REFERENCE LINE
%% ========================================================================

xline(ax, 0, '--', ...
    'Color', zeroGray, ...
    'LineWidth', 1.2);

%% ========================================================================
% 7. CONFIDENCE INTERVALS
%% ========================================================================

capHeight = 0.065;

for i = 1:2

    % Main horizontal 95% CI
    plot(ax, ...
        [plotLow(i), plotHigh(i)], ...
        [y(i), y(i)], ...
        '-', ...
        'Color', darkGray, ...
        'LineWidth', 2.0);

    % Left cap
    plot(ax, ...
        [plotLow(i), plotLow(i)], ...
        [y(i)-capHeight, y(i)+capHeight], ...
        '-', ...
        'Color', darkGray, ...
        'LineWidth', 1.6);

    % Right cap
    plot(ax, ...
        [plotHigh(i), plotHigh(i)], ...
        [y(i)-capHeight, y(i)+capHeight], ...
        '-', ...
        'Color', darkGray, ...
        'LineWidth', 1.6);

end

%% ========================================================================
% 8. POINT ESTIMATES
%% ========================================================================

% ------------------------------------------------------------------------
% Significant contrast:
% Intra-left - Intra-right
% Filled marker
% ------------------------------------------------------------------------

scatter(ax, ...
    plotEstimate(1), y(1), ...
    105, ...
    'o', ...
    'MarkerFaceColor', darkGray, ...
    'MarkerEdgeColor', darkGray, ...
    'LineWidth', 1.5);

% ------------------------------------------------------------------------
% Null contrast:
% Inter - mean(Intra-left, Intra-right)
% Open marker
% ------------------------------------------------------------------------

scatter(ax, ...
    plotEstimate(2), y(2), ...
    105, ...
    'o', ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', darkGray, ...
    'LineWidth', 1.8);

%% ========================================================================
% 9. HOLM-ADJUSTED P-VALUE LABELS
%% ========================================================================

% Put annotations inside the plotting area so they cannot be clipped.

pTextX = 0.014;

text(ax, ...
    pTextX, y(1), ...
    sprintf('Holm p = %.4f', plotHolmP(1)), ...
    'FontName', 'Times New Roman', ...
    'FontSize', 9.5, ...
    'FontWeight', 'bold', ...
    'Color', darkGray, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', ...
    'Interpreter', 'none');

text(ax, ...
    pTextX, y(2), ...
    sprintf('Holm p = %.3f', plotHolmP(2)), ...
    'FontName', 'Times New Roman', ...
    'FontSize', 9.5, ...
    'FontWeight', 'bold', ...
    'Color', darkGray, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', ...
    'Interpreter', 'none');

%% ========================================================================
% 10. AXES
%% ========================================================================

xlim(ax, [-0.050 0.045]);
ylim(ax, [0.50 2.50]);

% Significant contrast at top
set(ax, 'YDir', 'reverse');

yticks(ax, y);
yticklabels(ax, plotNames);

% Use concise x-axis label.
xlabel(ax, ...
    'Contrast estimate (\Delta db-wPLI)', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 10.5, ...
    'FontWeight', 'bold');

%% ========================================================================
% 11. X TICKS
%% ========================================================================

xticks(ax, [-0.04 -0.02 0 0.02 0.04]);

xticklabels(ax, { ...
    '-0.04', ...
    '-0.02', ...
    '0', ...
    '0.02', ...
    '0.04'});

%% ========================================================================
% 12. GLOBAL AXES STYLE
%% ========================================================================

set(ax, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 9.5, ...
    'FontWeight', 'bold', ...
    'LineWidth', 1.1, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'TickLength', [0.018 0.018], ...
    'Layer', 'top');

ax.XColor = axisGray;
ax.YColor = axisGray;

%% ========================================================================
% 13. MANUAL AXES POSITION
%
% Large left margin is required for the two contrast labels.
% This prevents labels from compressing the forest plot.
%% ========================================================================

ax.Position = [0.42 0.20 0.55 0.72];

%% ========================================================================
% 14. OPTIONAL ZERO-LINE LABEL
%
% Keep disabled for the final figure because the dashed line is obvious.
%% ========================================================================

% text(ax, 0.001, 2.40, 'No difference', ...
%     'FontName','Times New Roman', ...
%     'FontSize',8.5, ...
%     'Color',zeroGray);

%% ========================================================================
% 15. NO TITLE / PANEL LETTER
%
% Add panel label "b" later during PowerPoint figure assembly.
%% ========================================================================

%% ========================================================================
% 16. EXPORT
%% ========================================================================

pngFile = fullfile( ...
    outDir, ...
    'Figure3B_Alpha_Planned_Contrasts_FINAL_600dpi.png');

pdfFile = fullfile( ...
    outDir, ...
    'Figure3B_Alpha_Planned_Contrasts_FINAL_vector.pdf');

exportgraphics( ...
    fig, ...
    pngFile, ...
    'Resolution', 600, ...
    'BackgroundColor', 'white');

exportgraphics( ...
    fig, ...
    pdfFile, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

%% ========================================================================
% 17. VERIFY NUMBERS IN COMMAND WINDOW
%% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('FIGURE 3B — FROZEN ALPHA PLANNED CONTRASTS\n');
fprintf('============================================================\n\n');

fprintf('Contrast 1:\n');
fprintf('Inter - mean(Intra-left, Intra-right)\n');
fprintf('Estimate  = %.7f\n', C1_estimate);
fprintf('95%% CI    = [%.6f, %.6f]\n', C1_CI_low, C1_CI_high);
fprintf('Raw p     = %.5f\n', C1_rawP);
fprintf('Holm p    = %.5f\n\n', C1_holmP);

fprintf('Contrast 2:\n');
fprintf('Intra-left - Intra-right\n');
fprintf('Estimate  = %.6f\n', C2_estimate);
fprintf('95%% CI    = [%.6f, %.6f]\n', C2_CI_low, C2_CI_high);
fprintf('Raw p     = %.7f\n', C2_rawP);
fprintf('Holm p    = %.7f\n\n', C2_holmP);

fprintf('Saved PNG:\n%s\n\n', pngFile);
fprintf('Saved vector PDF:\n%s\n', pdfFile);

fprintf('\n============================================================\n');
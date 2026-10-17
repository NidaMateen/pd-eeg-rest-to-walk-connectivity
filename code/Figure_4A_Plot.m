%% ========================================================================
% FIGURE 4A — FINAL
% Sensitivity of the primary theta-band group effect to QC exclusions
%
% Outcome:
%   Mean Δ db-wPLI across the 32-edge theta NBS component
%
% Effect estimate:
%   Age- and sex-adjusted PD-HC regression coefficient (beta)
%
% Negative beta:
%   Smaller Rest-to-Walk theta increase in PD than HC
%
% Frozen sensitivity analyses:
%
%   Baseline complete network
%       N = 117 (22 HC / 95 PD)
%       beta = -0.045949
%       95% CI [-0.069387, -0.022511]
%
%   Exclude >=5 bad channels
%       N = 114 (22 HC / 92 PD)
%       beta = -0.044214
%       95% CI [-0.067589, -0.020840]
%
%   Exclude >20% rejected ICs
%       N = 106 (21 HC / 85 PD)
%       beta = -0.044683
%       95% CI [-0.069112, -0.020254]
%
%   Exclude both QC burdens
%       N = 103 (21 HC / 82 PD)
%       beta = -0.042800
%       95% CI [-0.067169, -0.018432]
%
%   Exclude all manual-ICA subjects
%       N = 105 (21 HC / 84 PD)
%       beta = -0.044803
%       95% CI [-0.069273, -0.020332]
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
% 2. FROZEN RESULTS
%% ========================================================================

analysisNames = { ...
    'Baseline complete network', ...
    'Exclude >=5 bad channels', ...
    'Exclude >20% rejected ICs', ...
    'Exclude both QC burdens', ...
    'Exclude all manual-ICA subjects'};

Ntotal = [117, 114, 106, 103, 105];

NHC = [22, 22, 21, 21, 21];

NPD = [95, 92, 85, 82, 84];

beta = [ ...
    -0.045949, ...
    -0.044214, ...
    -0.044683, ...
    -0.042800, ...
    -0.044803];

CI_low = [ ...
    -0.069387, ...
    -0.067589, ...
    -0.069112, ...
    -0.067169, ...
    -0.069273];

CI_high = [ ...
    -0.022511, ...
    -0.020840, ...
    -0.020254, ...
    -0.018432, ...
    -0.020332];

pValue = [ ...
    0.000173, ...
    0.000285, ...
    0.000448, ...
    0.000735, ...
    0.000444];

%% ========================================================================
% 3. EFFECT-MAGNITUDE RETENTION
%
% Absolute beta relative to baseline
%% ========================================================================

retainedPct = abs(beta) ./ abs(beta(1)) .* 100;

%% ========================================================================
% 4. COLORS
%% ========================================================================

baselineColor = [0.18 0.18 0.18];
sensitivityColor = [0.43 0.43 0.43];

zeroColor = [0.72 0.72 0.72];
axisColor = [0.12 0.12 0.12];

%% ========================================================================
% 5. CREATE FIGURE
%% ========================================================================

fig = figure( ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [3 3 13.5 8.5]);

ax = axes(fig);
hold(ax, 'on');

%% ========================================================================
% 6. ROW POSITIONS
%% ========================================================================

y = 1:5;

%% ========================================================================
% 7. ZERO REFERENCE LINE
%% ========================================================================

xline(ax, 0, '--', ...
    'Color', zeroColor, ...
    'LineWidth', 1.2);

%% ========================================================================
% 8. CONFIDENCE INTERVALS
%% ========================================================================

capHeight = 0.065;

for i = 1:5

    if i == 1
        thisColor = baselineColor;
        thisWidth = 2.2;
    else
        thisColor = sensitivityColor;
        thisWidth = 1.8;
    end

    % Horizontal CI
    plot(ax, ...
        [CI_low(i), CI_high(i)], ...
        [y(i), y(i)], ...
        '-', ...
        'Color', thisColor, ...
        'LineWidth', thisWidth);

    % Left CI cap
    plot(ax, ...
        [CI_low(i), CI_low(i)], ...
        [y(i)-capHeight, y(i)+capHeight], ...
        '-', ...
        'Color', thisColor, ...
        'LineWidth', 1.4);

    % Right CI cap
    plot(ax, ...
        [CI_high(i), CI_high(i)], ...
        [y(i)-capHeight, y(i)+capHeight], ...
        '-', ...
        'Color', thisColor, ...
        'LineWidth', 1.4);
end

%% ========================================================================
% 9. POINT ESTIMATES
%% ========================================================================

% ------------------------------------------------------------------------
% Baseline:
% filled diamond
% ------------------------------------------------------------------------

scatter(ax, ...
    beta(1), y(1), ...
    115, ...
    'd', ...
    'MarkerFaceColor', baselineColor, ...
    'MarkerEdgeColor', baselineColor, ...
    'LineWidth', 1.4);

% ------------------------------------------------------------------------
% Sensitivity analyses:
% filled circles
% ------------------------------------------------------------------------

for i = 2:5

    scatter(ax, ...
        beta(i), y(i), ...
        90, ...
        'o', ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', sensitivityColor, ...
        'LineWidth', 1.8);

end

%% ========================================================================
% 10. OPTIONAL SAMPLE-SIZE TEXT
%
% Put N on the right side of the plot.
%% ========================================================================

sampleTextX = -0.0065;

for i = 1:5

    text(ax, ...
        sampleTextX, y(i), ...
        sprintf('N = %d', Ntotal(i)), ...
        'FontName', 'Times New Roman', ...
        'FontSize', 9, ...
        'FontWeight', 'bold', ...
        'Color', sensitivityColor, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');

end

%% ========================================================================
% 11. AXES LIMITS
%% ========================================================================

xlim(ax, [-0.078 0.012]);
ylim(ax, [0.45 5.55]);

% Baseline on top
set(ax, 'YDir', 'reverse');

%% ========================================================================
% 12. Y-AXIS LABELS
%% ========================================================================

yticks(ax, y);

yticklabels(ax, { ...
    'Baseline complete network', ...
    'Exclude ≥ 5 bad channels', ...
    'Exclude > 20% rejected ICs', ...
    'Exclude both QC burdens', ...
    'Exclude all manual-ICA subjects'});

%% ========================================================================
% 13. X-AXIS
%% ========================================================================

xticks(ax, [-0.07 -0.06 -0.05 -0.04 -0.03 -0.02 -0.01 0]);

xticklabels(ax, { ...
    '-0.07', ...
    '-0.06', ...
    '-0.05', ...
    '-0.04', ...
    '-0.03', ...
    '-0.02', ...
    '-0.01', ...
    '0'});

xlabel(ax, ...
    'Adjusted PD-HC \beta (\Delta db-wPLI)', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 10.5, ...
    'FontWeight', 'bold');

%% ========================================================================
% 14. GLOBAL STYLE
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

ax.XColor = axisColor;
ax.YColor = axisColor;

%% ========================================================================
% 15. AXES POSITION
%
% Extra left margin for long sensitivity-analysis labels.
%% ========================================================================

ax.Position = [0.40 0.18 0.56 0.75];

%% ========================================================================
% 16. NO TITLE / PANEL LETTER
%
% Add "a" later when assembling Figure 4 in PowerPoint.
%% ========================================================================

%% ========================================================================
% 17. EXPORT
%% ========================================================================

pngFile = fullfile( ...
    outDir, ...
    'Figure4A_Theta_QC_Sensitivity_FINAL_600dpi.png');

pdfFile = fullfile( ...
    outDir, ...
    'Figure4A_Theta_QC_Sensitivity_FINAL_vector.pdf');

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
% 18. VERIFY VALUES IN COMMAND WINDOW
%% ========================================================================

fprintf('\n');
fprintf('==============================================================\n');
fprintf('FIGURE 4A — THETA COMPONENT QC SENSITIVITY\n');
fprintf('==============================================================\n\n');

for i = 1:length(beta)

    fprintf('%s\n', analysisNames{i});
    fprintf('  N              = %d (%d HC / %d PD)\n', ...
        Ntotal(i), NHC(i), NPD(i));

    fprintf('  beta           = %.6f\n', beta(i));

    fprintf('  95%% CI         = [%.6f, %.6f]\n', ...
        CI_low(i), CI_high(i));

    fprintf('  p              = %.6f\n', pValue(i));

    fprintf('  retained effect = %.2f%% of baseline magnitude\n\n', ...
        retainedPct(i));

end

fprintf('Minimum retained magnitude = %.2f%%\n', min(retainedPct));

fprintf('\nSaved PNG:\n%s\n', pngFile);
fprintf('\nSaved vector PDF:\n%s\n', pdfFile);

fprintf('\n==============================================================\n');
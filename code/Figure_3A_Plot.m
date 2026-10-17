%% ========================================================================
% FIGURE 3A
% Alpha-band NetworkType reconfiguration in HC vs PD
%
% Final frozen STEP10R descriptive summaries:
% Subject-level summary = median Δ db-wPLI within each NetworkType
% Δ db-wPLI = Walk - Rest
%
% Groups:
%   HC = 23
%   PD = 107
%
% Network types:
%   IntraLeft      = 6 edges
%   IntraRight     = 6 edges
%   Interhemispheric = 16 edges
%
% Plot:
%   Group mean of subject medians ± 95% CI
%
% NOTE:
%   This figure visualizes the significant Group × NetworkType interaction.
%   Do NOT interpret individual HC-vs-PD differences at each class as
%   separately significant tests.
%% ========================================================================

clear;
clc;
close all;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% -------------------- OUTPUT DIRECTORY ---------------------------------

outDir = fullfile(cfg.datasetRoot,'PAPER3_PRIMARY_NBS','FIGURES');

if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% -------------------- FROZEN DATA ---------------------------------------

networkNames = {'Intra-left','Intra-right','Interhemispheric'};

nHC = 23;
nPD = 107;

% Mean of subject-level median Δ db-wPLI
HC_mean = [-0.023316, -0.032718, -0.018294];
PD_mean = [-0.035985, -0.019684, -0.021954];

% SD of subject-level median Δ db-wPLI
HC_sd = [0.055626, 0.051432, 0.040811];
PD_sd = [0.053279, 0.043912, 0.036041];

%% -------------------- 95% CONFIDENCE INTERVALS --------------------------

% t-based 95% confidence intervals around descriptive group means

tHC = tinv(0.975, nHC - 1);
tPD = tinv(0.975, nPD - 1);

HC_ci = tHC .* HC_sd ./ sqrt(nHC);
PD_ci = tPD .* PD_sd ./ sqrt(nPD);

%% -------------------- COLORS --------------------------------------------

% Match Figure 2 visual language:
% HC = neutral gray
% PD = blue

HCcolor = [0.42 0.42 0.42];
PDcolor = [0.20 0.43 0.66];

%% -------------------- FIGURE --------------------------------------------

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3 3 10.5 8.5]);

ax = axes(fig);
hold(ax,'on');

x = 1:3;

%% Zero reference line

yline(0,'--', ...
    'Color',[0.60 0.60 0.60], ...
    'LineWidth',1.1);

%% Connecting lines

plot(x, HC_mean, ...
    '-', ...
    'Color',HCcolor, ...
    'LineWidth',2.0);

plot(x, PD_mean, ...
    '-', ...
    'Color',PDcolor, ...
    'LineWidth',2.0);

%% Error bars

errorbar(x, HC_mean, HC_ci, ...
    'LineStyle','none', ...
    'Color',HCcolor, ...
    'LineWidth',1.4, ...
    'CapSize',8);

errorbar(x, PD_mean, PD_ci, ...
    'LineStyle','none', ...
    'Color',PDcolor, ...
    'LineWidth',1.4, ...
    'CapSize',8);

%% Markers
% Use scatter separately for MATLAB R2023b compatibility

scatter(x, HC_mean, ...
    72, ...
    'o', ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor',HCcolor, ...
    'LineWidth',1.8);

scatter(x, PD_mean, ...
    72, ...
    'o', ...
    'MarkerFaceColor',PDcolor, ...
    'MarkerEdgeColor',PDcolor, ...
    'LineWidth',1.8);

%% -------------------- AXES ----------------------------------------------

xlim([0.65 3.35]);
ylim([-0.065 0.015]);

xticks(x);
xticklabels(networkNames);

ylabel('\Delta db-wPLI (Walk - Rest)', ...
    'FontName','Times New Roman', ...
    'FontSize',11, ...
    'FontWeight','bold');

set(ax, ...
    'FontName','Times New Roman', ...
    'FontSize',10.5, ...
    'FontWeight','bold', ...
    'LineWidth',1.1, ...
    'Box','off', ...
    'TickDir','out', ...
    'Layer','top');

%% -------------------- GROUP LABELS --------------------------------------

% Direct labels are cleaner than a legend

text(3.08, HC_mean(3)+0.0038, 'HC', ...
    'Color',HCcolor, ...
    'FontName','Times New Roman', ...
    'FontSize',10.5, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left');

text(3.08, PD_mean(3)-0.0040, 'PD', ...
    'Color',PDcolor, ...
    'FontName','Times New Roman', ...
    'FontSize',10.5, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left');

%% -------------------- OPTIONAL SAMPLE SIZE ------------------------------
% % Keep subtle; delete these two lines if you prefer a cleaner panel.
% 
% text(0.74, 0.009, 'HC, n = 23', ...
%     'Color',HCcolor, ...
%     'FontName','Times New Roman', ...
%     'FontSize',9, ...
%     'FontWeight','bold');
% 
% text(0.74, 0.003, 'PD, n = 107', ...
%     'Color',PDcolor, ...
%     'FontName','Times New Roman', ...
%     'FontSize',9, ...
%     'FontWeight','bold');

%% -------------------- CLEAN APPEARANCE ----------------------------------

ax.XColor = [0.15 0.15 0.15];
ax.YColor = [0.15 0.15 0.15];

% No title and no panel letter.
% Add "a" later during final PowerPoint assembly.

%% -------------------- EXPORT --------------------------------------------

pngFile = fullfile(outDir, ...
    'Figure3A_Alpha_NetworkType_FINAL_600dpi.png');

pdfFile = fullfile(outDir, ...
    'Figure3A_Alpha_NetworkType_FINAL_vector.pdf');

exportgraphics(fig, pngFile, ...
    'Resolution',600, ...
    'BackgroundColor','white');

exportgraphics(fig, pdfFile, ...
    'ContentType','vector', ...
    'BackgroundColor','white');

fprintf('\nSaved Figure 3A:\n');
fprintf('%s\n', pngFile);
fprintf('%s\n\n', pdfFile);

%% -------------------- VERIFY VALUES -------------------------------------

fprintf('Alpha NetworkType descriptive summaries\n');
fprintf('-----------------------------------------\n');

for k = 1:3
    fprintf('%-18s HC = %.6f ± %.6f CI | PD = %.6f ± %.6f CI\n', ...
        networkNames{k}, ...
        HC_mean(k), HC_ci(k), ...
        PD_mean(k), PD_ci(k));
end

fprintf('\nFrozen omnibus Group x NetworkType result:\n');
fprintf('F(2,256) = 4.6496, raw p = 0.010387, Holm p = 0.031162\n');
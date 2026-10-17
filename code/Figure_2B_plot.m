%% ========================================================================
%  PAPER 3 — FIGURE 2B
%  Subject-level theta-band Rest-to-Walk reconfiguration
%  across the 32-edge primary NBS component
%
%  Outcome:
%       Delta db-wPLI = Walk - Rest
%
%  Primary complete-network cohort:
%       HC = 22
%       PD = 95
%
%  Figure shows:
%       1. Individual participants
%       2. Boxplots: median and IQR
%       3. White diamonds: group means
%       4. Black error bars: 95% CI of group mean
%       5. Dashed horizontal line: Delta = 0
%
%  IMPORTANT:
%  db-wPLI is dimensionless.
%
%  No panel letter/title is added here.
%  Add "b" and any panel heading later in PowerPoint.
% ========================================================================

clear;
close all;
clc;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% ========================================================================
% 1. LOAD THE FROZEN PRIMARY NBS RESULTS
% ========================================================================

resultFile = ...
    fullfile(cfg.datasetRoot,'PAPER3_PRIMARY_NBS','Paper3_PRIMARY_NBS_RESULTS.mat');

% If the expected path no longer exists, select the file manually
if ~exist(resultFile,'file')

    warning('Expected primary NBS result file was not found.');

    [fileName,filePath] = uigetfile( ...
        '*.mat', ...
        'Select Paper3_PRIMARY_NBS_RESULTS.mat');

    if isequal(fileName,0)
        error('No MAT file selected.');
    end

    resultFile = fullfile(filePath,fileName);

end

fprintf('\nLoading primary NBS results:\n%s\n\n',resultFile);

R = load(resultFile);

%% ========================================================================
% 2. CHECK REQUIRED VARIABLES
% ========================================================================

requiredFields = { ...
    'subjects', ...
    'groups', ...
    'edgeNames', ...
    'bandNames', ...
    'deltaNBS'};

for k = 1:numel(requiredFields)

    if ~isfield(R,requiredFields{k})

        error( ...
            'Required field "%s" is missing from the MAT file.', ...
            requiredFields{k});

    end

end

subjects  = string(R.subjects(:));
groups    = upper(string(R.groups(:)));
edgeNames = string(R.edgeNames(:));
bandNames = lower(string(R.bandNames(:)));

%% ========================================================================
% 3. VERIFY PRIMARY SAMPLE
% ========================================================================

N   = numel(subjects);
nHC = sum(groups == "HC");
nPD = sum(groups == "PD");

fprintf('Primary cohort:\n');
fprintf('Total N = %d\n',N);
fprintf('HC      = %d\n',nHC);
fprintf('PD      = %d\n\n',nPD);

assert(N == 117, ...
    'Expected primary complete-network N = 117.');

assert(nHC == 22, ...
    'Expected HC n = 22.');

assert(nPD == 95, ...
    'Expected PD n = 95.');

%% ========================================================================
% 4. IDENTIFY THETA BAND
% ========================================================================

thetaIdx = find(bandNames == "theta");

assert(numel(thetaIdx) == 1, ...
    'Theta band could not be uniquely identified.');

%% ========================================================================
% 5. EXACT 32 EDGES OF THE PRIMARY THETA NBS COMPONENT
% ========================================================================

componentEdges = [ ...
    "FrontalRight__CentralRight"
    "ParietalMidline__ParietalRight"
    "FrontalLeft__ParietalMidline"
    "ParietalLeft__ParietalRight"
    "FrontalLeft__FrontalRight"
    "FrontalLeft__TemporalRight"
    "ParietalMidline__TemporalRight"
    "FrontalLeft__CentralRight"
    "CentralRight__ParietalMidline"
    "FrontalMidline__ParietalLeft"
    "FrontalRight__ParietalRight"
    "FrontalLeft__ParietalRight"
    "TemporalRight__Occipital"
    "FrontalLeft__TemporalLeft"
    "FrontalMidline__Occipital"
    "ParietalLeft__TemporalRight"
    "FrontalLeft__ParietalLeft"
    "FrontalRight__ParietalLeft"
    "ParietalRight__Occipital"
    "FrontalRight__TemporalRight"
    "FrontalMidline__FrontalRight"
    "FrontalLeft__Occipital"
    "CentralRight__Occipital"
    "CentralRight__ParietalLeft"
    "FrontalMidline__CentralRight"
    "FrontalRight__CentralLeft"
    "FrontalLeft__FrontalMidline"
    "FrontalRight__TemporalLeft"
    "ParietalLeft__ParietalMidline"
    "CentralLeft__CentralRight"
    "FrontalMidline__ParietalRight"
    "CentralLeft__ParietalRight"];

assert(numel(componentEdges) == 32, ...
    'Expected exactly 32 component edges.');

%% ========================================================================
% 6. MATCH COMPONENT EDGES TO SAVED CONNECTIVITY ORDER
% ========================================================================

componentIdx = zeros(32,1);

for e = 1:32

    idx = find(edgeNames == componentEdges(e));

    % Allow reversed order if the saved edge naming convention differs
    if isempty(idx)

        parts = split(componentEdges(e),"__");

        reversedEdge = parts(2) + "__" + parts(1);

        idx = find(edgeNames == reversedEdge);

    end

    if numel(idx) ~= 1

        error( ...
            'Could not uniquely identify edge: %s', ...
            componentEdges(e));

    end

    componentIdx(e) = idx;

end

assert(numel(unique(componentIdx)) == 32, ...
    'Duplicate edge matching was detected.');

fprintf('All 32 theta component edges matched successfully.\n\n');

%% ========================================================================
% 7. CALCULATE SUBJECT-LEVEL COMPONENT-AVERAGE DELTA db-wPLI
% ========================================================================

% deltaNBS dimensions:
% Subject x Edge x FrequencyBand

thetaDelta = R.deltaNBS(:,:,thetaIdx);

% Mean across the 32 NBS-component edges for each participant
componentDelta = mean( ...
    thetaDelta(:,componentIdx), ...
    2, ...
    'omitnan');

componentDelta = componentDelta(:);

assert(all(isfinite(componentDelta)), ...
    'Nonfinite component-average values detected.');

HC = componentDelta(groups == "HC");
PD = componentDelta(groups == "PD");

%% ========================================================================
% 8. VERIFY AGAINST THE FROZEN RESULTS
% ========================================================================

meanHC = mean(HC);
meanPD = mean(PD);

medianHC = median(HC);
medianPD = median(PD);

fprintf('============================================================\n');
fprintf('FIGURE 2B DATA CHECK\n');
fprintf('============================================================\n');

fprintf('HC: n = %d | mean = %.6f | median = %.6f\n', ...
    numel(HC),meanHC,medianHC);

fprintf('PD: n = %d | mean = %.6f | median = %.6f\n', ...
    numel(PD),meanPD,medianPD);

fprintf('\nFrozen expected means:\n');
fprintf('HC = 0.067641\n');
fprintf('PD = 0.021181\n');

assert(abs(meanHC - 0.067641) < 1e-5, ...
    'HC mean does not reproduce the frozen result.');

assert(abs(meanPD - 0.021181) < 1e-5, ...
    'PD mean does not reproduce the frozen result.');

fprintf('\nPASS: frozen component-average means reproduced.\n');
fprintf('============================================================\n\n');

%% ========================================================================
% 9. GROUP-MEAN 95%% CONFIDENCE INTERVALS
% ========================================================================

semHC = std(HC,0) / sqrt(numel(HC));
semPD = std(PD,0) / sqrt(numel(PD));

ciHC = tinv(0.975,numel(HC)-1) * semHC;
ciPD = tinv(0.975,numel(PD)-1) * semPD;

%% ========================================================================
% 10. VISUAL SETTINGS
% ========================================================================

fontName = 'Times New Roman';

% Neutral control group
HCcolor = [0.30 0.30 0.30];

% Blue consistent with Figure 2A
PDcolor = [0.36 0.55 0.78];

%% ========================================================================
% 11. CREATE FIGURE
% ========================================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 9.5 10.5], ...
    'Renderer','painters');

ax = axes(fig);

hold(ax,'on');

%% ========================================================================
% 12. DETERMINE Y-AXIS LIMITS
%
% Keep all observations visible while creating a small margin.
% ========================================================================

allValues = [HC; PD];

dataMin = min(allValues);
dataMax = max(allValues);

dataRange = dataMax - dataMin;

yMin = dataMin - 0.08*dataRange;
yMax = dataMax + 0.08*dataRange;

% Round outward slightly for cleaner axis limits
yMin = floor(yMin*100)/100;
yMax = ceil(yMax*100)/100;

ylim(ax,[yMin yMax]);

%% ========================================================================
% 13. ZERO REFERENCE LINE
% ========================================================================

yline(ax,0,'--', ...
    'Color',[0.55 0.55 0.55], ...
    'LineWidth',1.2);

%% ========================================================================
% 14. BOXPLOTS
%
% These display median, IQR and whiskers.
% ========================================================================

boxchart(ax, ...
    ones(numel(HC),1), ...
    HC, ...
    'BoxFaceColor',HCcolor, ...
    'BoxFaceAlpha',0.12, ...
    'WhiskerLineColor',HCcolor, ...
    'MarkerStyle','none', ...
    'BoxWidth',0.36, ...
    'LineWidth',1.5);

boxchart(ax, ...
    2*ones(numel(PD),1), ...
    PD, ...
    'BoxFaceColor',PDcolor, ...
    'BoxFaceAlpha',0.12, ...
    'WhiskerLineColor',PDcolor, ...
    'MarkerStyle','none', ...
    'BoxWidth',0.36, ...
    'LineWidth',1.5);

%% ========================================================================
% 15. INDIVIDUAL PARTICIPANT VALUES
%
% swarmchart spreads overlapping observations without arbitrary
% random jitter.
% ========================================================================

swarmchart(ax, ...
    ones(numel(HC),1), ...
    HC, ...
    22, ...
    HCcolor, ...
    'filled', ...
    'MarkerFaceAlpha',0.60, ...
    'MarkerEdgeAlpha',0.60, ...
    'XJitterWidth',0.22);

swarmchart(ax, ...
    2*ones(numel(PD),1), ...
    PD, ...
    22, ...
    PDcolor, ...
    'filled', ...
    'MarkerFaceAlpha',0.48, ...
    'MarkerEdgeAlpha',0.48, ...
    'XJitterWidth',0.22);

%% ========================================================================
% 16. GROUP MEAN AND 95%% CI
%
% White diamond = arithmetic mean
% Black vertical line = 95% CI
% ========================================================================

errorbar(ax, ...
    1,meanHC,ciHC, ...
    'Color','k', ...
    'LineStyle','none', ...
    'LineWidth',2.0, ...
    'CapSize',9);

errorbar(ax, ...
    2,meanPD,ciPD, ...
    'Color','k', ...
    'LineStyle','none', ...
    'LineWidth',2.0, ...
    'CapSize',9);

scatter(ax, ...
    1,meanHC, ...
    85,'d', ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor','k', ...
    'LineWidth',1.8);

scatter(ax, ...
    2,meanPD, ...
    85,'d', ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor','k', ...
    'LineWidth',1.8);

%% ========================================================================
% 17. X-AXIS
%
% Hide MATLAB's default labels entirely.
% We will add controlled custom labels manually.
% ========================================================================

xlim(ax,[0.55 2.45]);

ax.XTick = [1 2];

% Suppress default labels
ax.XTickLabel = {'',''};

%% ========================================================================
% 18. CUSTOM GROUP LABELS
%
% This avoids MATLAB's broken multiline xtick rendering.
% ========================================================================

yl = ylim(ax);
yrange = yl(2) - yl(1);

% Place labels just below plotting region
labelY1 = yl(1) - 0.055*yrange;
labelY2 = yl(1) - 0.105*yrange;

% -------------------------
% HC
% -------------------------

text(ax, ...
    1,labelY1, ...
    'HC', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'FontName',fontName, ...
    'FontSize',15, ...
    'FontWeight','bold', ...
    'Color','k', ...
    'Clipping','off');

text(ax, ...
    1,labelY2, ...
    '(n = 22)', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'FontName',fontName, ...
    'FontSize',12, ...
    'FontWeight','bold', ...
    'Color','k', ...
    'Clipping','off');

% -------------------------
% PD
% -------------------------

text(ax, ...
    2,labelY1, ...
    'PD', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'FontName',fontName, ...
    'FontSize',15, ...
    'FontWeight','bold', ...
    'Color','k', ...
    'Clipping','off');

text(ax, ...
    2,labelY2, ...
    '(n = 95)', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'FontName',fontName, ...
    'FontSize',12, ...
    'FontWeight','bold', ...
    'Color','k', ...
    'Clipping','off');

%% ========================================================================
% 19. Y-AXIS LABEL
%
% db-wPLI is dimensionless.
% ========================================================================

ylabel(ax, ...
    '\Delta db-wPLI (Walk - Rest)', ...
    'FontName',fontName, ...
    'FontSize',15, ...
    'FontWeight','bold', ...
    'Color','k');

%% ========================================================================
% 20. AXIS APPEARANCE
% ========================================================================

set(ax, ...
    'FontName',fontName, ...
    'FontSize',13, ...
    'FontWeight','bold', ...
    'LineWidth',1.4, ...
    'TickDir','out', ...
    'TickLength',[0.018 0.018], ...
    'Box','off', ...
    'XColor','k', ...
    'YColor','k');

% No general gridlines
ax.XGrid = 'off';
ax.YGrid = 'off';

%% ========================================================================
% 21. CLEAN Y-TICK FORMATTING
% ========================================================================

% Use clean decimal formatting
ytickformat(ax,'%.2f');

%% ========================================================================
% 22. COMPACT LAYOUT
%
% Extra bottom room is reserved for HC/PD sample-size labels.
% ========================================================================

ax.Position = [0.20 0.19 0.76 0.77];

%% ========================================================================
% 23. EXPORT
% ========================================================================

outputFolder = fileparts(resultFile);

pdfFile = fullfile( ...
    outputFolder, ...
    'Figure2B_Theta_Component_Delta_FINAL.pdf');

pngFile = fullfile( ...
    outputFolder, ...
    'Figure2B_Theta_Component_Delta_FINAL_600dpi.png');

% Export entire FIGURE rather than only axes,
% because custom x labels extend outside the axes.
exportgraphics(fig, ...
    pdfFile, ...
    'ContentType','vector', ...
    'BackgroundColor','white');

exportgraphics(fig, ...
    pngFile, ...
    'Resolution',600, ...
    'BackgroundColor','white');

fprintf('Figure saved:\n');
fprintf('%s\n',pdfFile);
fprintf('%s\n',pngFile);
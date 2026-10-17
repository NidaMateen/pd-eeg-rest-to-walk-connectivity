%% ========================================================================
% FIGURE 4C — FINAL
% Complete-network vs available-case theta edge-effect concordance
%
% Paper 3:
% PD EEG Rest-to-Walking Network Reconfiguration
%
% Purpose
% -------
% Compare age/sex-adjusted PD-HC beta estimates for all 55 theta edges:
%
%   X-axis:
%       Complete-network cohort
%       N = 117 = 22 HC + 95 PD
%
%   Y-axis:
%       Available-case edgewise cohort
%       Edge-specific N = 117-130
%
% Frozen expected concordance:
%   Pearson r                 = 0.9843
%   Spearman rho              = 0.9866
%   Sign agreement            = 100%
%   Median absolute beta diff = 0.001668 db-wPLI
%
% IMPORTANT
% ---------
% HC3 robust covariance affects SE / CI / p-values, NOT the OLS beta
% coefficients themselves. Therefore ordinary least-squares beta estimates
% are sufficient for this concordance plot.
%
% MATLAB R2023b compatible
%% ========================================================================

clear;
clc;
close all;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% ========================================================================
% 1. PATHS
%% ========================================================================

datasetRoot = ...
    cfg.datasetRoot;

connectivityRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_FINAL_CONNECTIVITY');

statsInputFile = fullfile( ...
    connectivityRoot, ...
    'Paper3_CONNECTIVITY_STATS_INPUT.mat');

participantsFile = fullfile( ...
    datasetRoot, ...
    'participants.tsv');

outDir = fullfile( ...
    datasetRoot, ...
    'PAPER3_PRIMARY_NBS', ...
    'FIGURES');

if ~exist(outDir,'dir')
    mkdir(outDir);
end

if ~exist(statsInputFile,'file')
    error('Statistics-input MAT file not found:\n%s', statsInputFile);
end

if ~exist(participantsFile,'file')
    error('participants.tsv not found:\n%s', participantsFile);
end

%% ========================================================================
% 2. LOAD FROZEN CONNECTIVITY DATA
%% ========================================================================

S = load(statsInputFile);

requiredFields = { ...
    'subjects', ...
    'groups', ...
    'edgeNames', ...
    'bandNames', ...
    'deltaBio', ...
    'completeNetworkMask'};

for k = 1:numel(requiredFields)

    if ~isfield(S,requiredFields{k})
        error('Missing required field: %s', requiredFields{k});
    end

end

subjectsAll = string(S.subjects(:));
groupsAll   = string(S.groups(:));

edgeNames = string(S.edgeNames(:));
bandNames = string(S.bandNames(:));

deltaBio = double(S.deltaBio);

completeMask = logical(S.completeNetworkMask(:));

%% ========================================================================
% 3. STRUCTURAL CHECKS
%% ========================================================================

N_EDGES = 55;

if numel(edgeNames) ~= N_EDGES
    error('Expected 55 edges; found %d.', numel(edgeNames));
end

if size(deltaBio,1) ~= numel(subjectsAll)
    error('Subject dimension mismatch in deltaBio.');
end

if size(deltaBio,2) ~= N_EDGES
    error('Expected deltaBio second dimension = 55 edges.');
end

% Find theta band explicitly
thetaIdx = find(lower(bandNames) == "theta");

if numel(thetaIdx) ~= 1
    error('Could not uniquely identify theta band.');
end

%% ========================================================================
% 4. VERIFY PRIMARY COMPLETE-NETWORK COHORT
%% ========================================================================

nComplete = sum(completeMask);

nHCcomplete = sum( ...
    completeMask & groupsAll == "HC");

nPDcomplete = sum( ...
    completeMask & groupsAll == "PD");

fprintf('\nComplete-network cohort:\n');
fprintf('N  = %d\n', nComplete);
fprintf('HC = %d\n', nHCcomplete);
fprintf('PD = %d\n', nPDcomplete);

if nComplete ~= 117 || ...
   nHCcomplete ~= 22 || ...
   nPDcomplete ~= 95

    error(['Complete-network cohort does not reproduce ' ...
           '117 = 22 HC + 95 PD.']);

end

%% ========================================================================
% 5. LOAD AND ALIGN PARTICIPANT AGE / SEX
%% ========================================================================

P = readtable( ...
    participantsFile, ...
    'FileType','text', ...
    'Delimiter','\t', ...
    'TextType','string');

neededColumns = [ ...
    "participant_id", ...
    "age", ...
    "sex"];

for k = 1:numel(neededColumns)

    if ~ismember(neededColumns(k), ...
            string(P.Properties.VariableNames))

        error('participants.tsv lacks column: %s', ...
            neededColumns(k));

    end

end

participantID = string(P.participant_id);

[tf,loc] = ismember(subjectsAll,participantID);

if ~all(tf)

    fprintf('\nSubjects not found in participants.tsv:\n');
    disp(subjectsAll(~tf));

    error('Participant metadata alignment failed.');

end

%% ========================================================================
% 6. CONVERT AGE ROBUSTLY
%% ========================================================================

ageRaw = P.age(loc);

if isnumeric(ageRaw)

    age = double(ageRaw);

else

    age = str2double(string(ageRaw));

end

if any(~isfinite(age))
    error('Non-finite age values after participant alignment.');
end

%% ========================================================================
% 7. SEX CODING
%
% Frozen coding:
% F = 0
% M = 1
%% ========================================================================

sexRaw = upper(strtrim(string(P.sex(loc))));

if any(~ismember(sexRaw,["F","M"]))
    error('Sex values must be F or M.');
end

Sex = double(sexRaw == "M");

%% ========================================================================
% 8. GROUP CODING
%
% Frozen coding:
% HC = 0
% PD = 1
%% ========================================================================

if any(~ismember(groupsAll,["HC","PD"]))
    error('Unexpected group labels.');
end

Group = double(groupsAll == "PD");

%% ========================================================================
% 9. EXTRACT THETA Δ db-wPLI
%% ========================================================================

thetaDelta = squeeze(deltaBio(:,:,thetaIdx));

% Expected dimensions:
% subjects x 55 edges

if size(thetaDelta,1) ~= numel(subjectsAll) || ...
   size(thetaDelta,2) ~= N_EDGES

    error('Unexpected thetaDelta dimensions.');

end

%% ========================================================================
% 10. ESTIMATE EDGE BETAS
%
% Model for every edge:
%
% Δ db-wPLI =
%     beta0
%   + beta1 * Group
%   + beta2 * Age_c
%   + beta3 * Sex
%   + error
%
% beta1 = adjusted PD-HC effect
%% ========================================================================

betaComplete  = nan(N_EDGES,1);
betaAvailable = nan(N_EDGES,1);

nAvailable = nan(N_EDGES,1);
nHCavail   = nan(N_EDGES,1);
nPDavail   = nan(N_EDGES,1);

for e = 1:N_EDGES

    y = thetaDelta(:,e);

    %% --------------------------------------------------------------------
    % COMPLETE-NETWORK MODEL
    %% --------------------------------------------------------------------

    maskC = completeMask & ...
            isfinite(y) & ...
            isfinite(age) & ...
            isfinite(Sex);

    if sum(maskC) ~= 117

        error(['Edge %d complete-network analysis does not contain ' ...
               'exactly 117 participants.'],e);

    end

    yC = y(maskC);

    groupC = Group(maskC);
    ageC   = age(maskC);
    sexC   = Sex(maskC);

    % Center age within analysis sample.
    % Group coefficient is invariant to this centering.
    ageC = ageC - mean(ageC);

    XC = [ ...
        ones(sum(maskC),1), ...
        groupC, ...
        ageC, ...
        sexC];

    if rank(XC) ~= size(XC,2)
        error('Complete-network design rank deficient at edge %d.',e);
    end

    bC = XC \ yC;

    betaComplete(e) = bC(2);

    %% --------------------------------------------------------------------
    % AVAILABLE-CASE MODEL
    %% --------------------------------------------------------------------

    maskA = isfinite(y) & ...
            isfinite(age) & ...
            isfinite(Sex);

    yA = y(maskA);

    groupA = Group(maskA);
    ageA   = age(maskA);
    sexA   = Sex(maskA);

    ageA = ageA - mean(ageA);

    XA = [ ...
        ones(sum(maskA),1), ...
        groupA, ...
        ageA, ...
        sexA];

    if rank(XA) ~= size(XA,2)
        error('Available-case design rank deficient at edge %d.',e);
    end

    bA = XA \ yA;

    betaAvailable(e) = bA(2);

    nAvailable(e) = sum(maskA);
    nHCavail(e)   = sum(groupA == 0);
    nPDavail(e)   = sum(groupA == 1);

end

%% ========================================================================
% 11. CHECK AVAILABLE-CASE SAMPLE RANGE
%% ========================================================================

fprintf('\nAvailable-case sample sizes:\n');
fprintf('Minimum N = %d\n',min(nAvailable));
fprintf('Maximum N = %d\n',max(nAvailable));

if min(nAvailable) < 117 || max(nAvailable) > 130

    error(['Available-case sample range is outside expected ' ...
           '117-130 subjects.']);

end

%% ========================================================================
% 12. CONCORDANCE STATISTICS
%% ========================================================================

pearsonR = corr( ...
    betaComplete, ...
    betaAvailable, ...
    'Type','Pearson');

spearmanR = corr( ...
    betaComplete, ...
    betaAvailable, ...
    'Type','Spearman');

signAgreement = mean( ...
    sign(betaComplete) == sign(betaAvailable)) * 100;

medianAbsDifference = median( ...
    abs(betaComplete - betaAvailable));

fprintf('\n============================================================\n');
fprintf('COMPLETE vs AVAILABLE-CASE THETA EDGE BETAS\n');
fprintf('============================================================\n');

fprintf('Pearson r             = %.6f\n',pearsonR);
fprintf('Spearman rho          = %.6f\n',spearmanR);
fprintf('Sign agreement        = %.2f%%\n',signAgreement);
fprintf('Median |beta diff|    = %.6f\n',medianAbsDifference);

%% ========================================================================
% 13. VALIDATE AGAINST FROZEN RESULTS
%
% Expected:
% Pearson  = 0.9843
% Spearman = 0.9866
% Sign     = 100%
% Median absolute beta difference = 0.001668
%
% Small tolerance allowed for rounding / implementation precision.
%% ========================================================================

EXPECTED_PEARSON  = 0.9843;
EXPECTED_SPEARMAN = 0.9866;
EXPECTED_SIGN     = 100;
EXPECTED_MAD      = 0.001668;

if abs(pearsonR - EXPECTED_PEARSON) > 0.002
    error(['Pearson correlation does not reproduce frozen result. ' ...
           'Observed %.6f; expected approximately %.4f.'], ...
           pearsonR, EXPECTED_PEARSON);
end

if abs(spearmanR - EXPECTED_SPEARMAN) > 0.002
    error(['Spearman correlation does not reproduce frozen result. ' ...
           'Observed %.6f; expected approximately %.4f.'], ...
           spearmanR, EXPECTED_SPEARMAN);
end

if abs(signAgreement - EXPECTED_SIGN) > 0.01
    error(['Sign agreement does not reproduce frozen result. ' ...
           'Observed %.2f%%.'],signAgreement);
end

if abs(medianAbsDifference - EXPECTED_MAD) > 0.0005
    error(['Median absolute beta difference does not reproduce ' ...
           'frozen result. Observed %.6f.'], ...
           medianAbsDifference);
end

fprintf('\nPASS: frozen concordance metrics reproduced.\n');

%% ========================================================================
% 14. SAVE EDGE-LEVEL RECONSTRUCTION TABLE
%% ========================================================================

T = table( ...
    edgeNames, ...
    betaComplete, ...
    betaAvailable, ...
    nAvailable, ...
    nHCavail, ...
    nPDavail, ...
    'VariableNames',{ ...
        'Edge', ...
        'Beta_Complete117', ...
        'Beta_AvailableCase', ...
        'N_Available', ...
        'N_HC_Available', ...
        'N_PD_Available'});

csvFile = fullfile( ...
    outDir, ...
    'Figure4C_Complete_vs_Available_EdgeBetas.csv');

writetable(T,csvFile);

%% ========================================================================
% 15. FIGURE COLORS
%
% Neutral color: points represent EDGES, not HC or PD.
%% ========================================================================

pointColor = [0.35 0.35 0.35];
lineColor  = [0.62 0.62 0.62];
axisColor  = [0.12 0.12 0.12];

%% ========================================================================
% 16. CREATE FIGURE
%% ========================================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3 3 9.5 9.0]);

ax = axes(fig);
hold(ax,'on');

%% ========================================================================
% 17. AXIS RANGE
%
% Use common symmetric-ish limits based on both vectors.
%% ========================================================================

allBeta = [betaComplete; betaAvailable];

dataMin = min(allBeta);
dataMax = max(allBeta);

padding = 0.006;

lowerLim = floor((dataMin-padding)*100)/100;
upperLim = ceil((dataMax+padding)*100)/100;

% Ensure 0 can be seen if reasonably close.
lowerLim = min(lowerLim,-0.01);
upperLim = max(upperLim,0.01);

%% ========================================================================
% 18. IDENTITY LINE y = x
%% ========================================================================

plot(ax, ...
    [lowerLim upperLim], ...
    [lowerLim upperLim], ...
    '--', ...
    'Color',lineColor, ...
    'LineWidth',1.3);

%% ========================================================================
% 19. EDGE SCATTER
%% ========================================================================

scatter(ax, ...
    betaComplete, ...
    betaAvailable, ...
    46, ...
    'o', ...
    'MarkerFaceColor',pointColor, ...
    'MarkerEdgeColor','w', ...
    'LineWidth',0.6);

%% ========================================================================
% 20. AXES
%% ========================================================================

xlim(ax,[lowerLim upperLim]);
ylim(ax,[lowerLim upperLim]);

axis(ax,'square');

xlabel(ax, ...
    {'Adjusted PD-HC \beta', ...
     'complete-network cohort'}, ...
    'FontName','Times New Roman', ...
    'FontSize',10.5, ...
    'FontWeight','bold');

ylabel(ax, ...
    {'Adjusted PD-HC \beta', ...
     'available-case analysis'}, ...
    'FontName','Times New Roman', ...
    'FontSize',10.5, ...
    'FontWeight','bold');

%% ========================================================================
% 21. METRIC ANNOTATION
%% ========================================================================

xText = lowerLim + 0.055*(upperLim-lowerLim);
yText = upperLim - 0.075*(upperLim-lowerLim);

metricText = sprintf([ ...
    'Pearson r = %.3f\n' ...
    'Spearman \\rho = %.3f\n' ...
    'Sign agreement = %.0f%%\n' ...
    'Median |\\Delta\\beta| = %.4f'], ...
    pearsonR, ...
    spearmanR, ...
    signAgreement, ...
    medianAbsDifference);

text(ax, ...
    xText, ...
    yText, ...
    metricText, ...
    'FontName','Times New Roman', ...
    'FontSize',8, ...
    'FontWeight','bold', ...
    'Color',axisColor, ...
    'VerticalAlignment','top', ...
    'HorizontalAlignment','left', ...
    'Interpreter','tex');

%% ========================================================================
% 22. STYLE
%% ========================================================================

set(ax, ...
    'FontName','Times New Roman', ...
    'FontSize',9.5, ...
    'FontWeight','bold', ...
    'LineWidth',1.1, ...
    'Box','off', ...
    'TickDir','out', ...
    'TickLength',[0.018 0.018], ...
    'Layer','top');

ax.XColor = axisColor;
ax.YColor = axisColor;

%% ========================================================================
% 23. NO TITLE / PANEL LETTER
%
% Add "c" when assembling Figure 4 in PowerPoint.
%% ========================================================================

%% ========================================================================
% 24. EXPORT
%% ========================================================================

pngFile = fullfile( ...
    outDir, ...
    'Figure4C_Complete_vs_Available_FINAL_600dpi.png');

pdfFile = fullfile( ...
    outDir, ...
    'Figure4C_Complete_vs_Available_FINAL_vector.pdf');

exportgraphics( ...
    fig, ...
    pngFile, ...
    'Resolution',600, ...
    'BackgroundColor','white');

exportgraphics( ...
    fig, ...
    pdfFile, ...
    'ContentType','vector', ...
    'BackgroundColor','white');

%% ========================================================================
% 25. FINAL OUTPUT
%% ========================================================================

fprintf('\nSaved Figure 4C PNG:\n%s\n',pngFile);
fprintf('\nSaved Figure 4C vector PDF:\n%s\n',pdfFile);
fprintf('\nSaved edge table:\n%s\n',csvFile);

fprintf('\n============================================================\n');
fprintf('FIGURE 4C COMPLETE\n');
fprintf('============================================================\n');
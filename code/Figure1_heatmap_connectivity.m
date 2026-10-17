%% ========================================================================
% FIGURE — GROUP-MEAN REST-TO-WALK db-wPLI HEATMAPS
%
%                  HC                          PD
%
% Theta       11 x 11 heatmap            11 x 11 heatmap
% Alpha       11 x 11 heatmap            11 x 11 heatmap
% Beta        11 x 11 heatmap            11 x 11 heatmap
%
% Data:
%       regionalDelta(:,:,band)
%
% Delta db-wPLI = Walk - Rest
%
% Diagonal:
%       Gray = no regional edge
%
% MATLAB R2023b
%% ========================================================================

clear;
clc;
close all;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% ========================================================================
% 1. ROOT DIRECTORY
%% ========================================================================

rootDir = ...
    fullfile(cfg.datasetRoot,'PAPER3_FINAL_CONNECTIVITY');

if ~isfolder(rootDir)
    error('Root directory does not exist:\n%s', rootDir);
end

%% ========================================================================
% 2. REGION ABBREVIATIONS
%
% Use a fixed order matching your regional matrices.
%% ========================================================================

regionLabels = { ...
    'FL','FM','FR', ...
    'CL','CR', ...
    'PL','PM','PR', ...
    'TL','TR','O'};

nRegions = 11;

%% ========================================================================
% 3. FIND SUBJECT DIRECTORIES
%% ========================================================================

subjectFolders = dir(fullfile(rootDir,'sub-*'));
subjectFolders = subjectFolders([subjectFolders.isdir]);

fprintf('\nSubject folders found: %d\n',numel(subjectFolders));

%% ========================================================================
% 4. STORAGE
%% ========================================================================

thetaHC = {};
thetaPD = {};

alphaHC = {};
alphaPD = {};

betaHC = {};
betaPD = {};

%% ========================================================================
% 5. LOAD EACH SUBJECT
%% ========================================================================

for s = 1:numel(subjectFolders)

    subjectDir = fullfile( ...
        subjectFolders(s).folder, ...
        subjectFolders(s).name);

    %% Find FINAL_CONNECTIVITY MAT

    matFiles = dir(fullfile( ...
        subjectDir, ...
        '*FINAL_CONNECTIVITY*.mat'));

    if isempty(matFiles)

        fprintf('SKIP %s: connectivity MAT not found.\n', ...
            subjectFolders(s).name);

        continue;

    end

    matFile = fullfile( ...
        matFiles(1).folder, ...
        matFiles(1).name);

    %% Load required variables

    D = load( ...
        matFile, ...
        'Subject', ...
        'Group', ...
        'regionalDelta', ...
        'bandNames');

    %% Check dimensions

    if size(D.regionalDelta,1) ~= 11 || ...
       size(D.regionalDelta,2) ~= 11

        error('%s does not contain an 11 x 11 regional matrix.', ...
            string(D.Subject));

    end

    %% ---------------------------------------------------------------
    % Identify bands from bandNames
    %% ---------------------------------------------------------------

    currentBandNames = ...
        lower(strtrim(string(D.bandNames(:))));

    thetaIdx = find(currentBandNames == "theta",1);
    alphaIdx = find(currentBandNames == "alpha",1);
    betaIdx  = find(currentBandNames == "beta",1);

    if isempty(thetaIdx) || isempty(alphaIdx) || isempty(betaIdx)

        fprintf('\nBand names in %s:\n',string(D.Subject));
        disp(D.bandNames);

        error('Theta/alpha/beta could not be identified.');

    end

    %% ---------------------------------------------------------------
    % Extract Walk-Rest matrices
    %% ---------------------------------------------------------------

    thetaMatrix = double(D.regionalDelta(:,:,thetaIdx));
    alphaMatrix = double(D.regionalDelta(:,:,alphaIdx));
    betaMatrix  = double(D.regionalDelta(:,:,betaIdx));

    %% ---------------------------------------------------------------
    % Force symmetry
    %% ---------------------------------------------------------------

    thetaMatrix = make_symmetric(thetaMatrix);
    alphaMatrix = make_symmetric(alphaMatrix);
    betaMatrix  = make_symmetric(betaMatrix);

    %% ---------------------------------------------------------------
    % Require all 55 unique edges
    %% ---------------------------------------------------------------

    upperMask = triu(true(nRegions),1);

    completeNetwork = ...
        all(isfinite(thetaMatrix(upperMask))) && ...
        all(isfinite(alphaMatrix(upperMask))) && ...
        all(isfinite(betaMatrix(upperMask)));

    if ~completeNetwork

        fprintf('EXCLUDE %s: incomplete 55-edge network.\n', ...
            string(D.Subject));

        continue;

    end

    %% ---------------------------------------------------------------
    % IMPORTANT:
    % Set diagonal to NaN.
    %
    % It will subsequently be shown as gray.
    %% ---------------------------------------------------------------

    thetaMatrix(1:nRegions+1:end) = NaN;
    alphaMatrix(1:nRegions+1:end) = NaN;
    betaMatrix(1:nRegions+1:end)  = NaN;

    %% ---------------------------------------------------------------
    % Store according to group
    %% ---------------------------------------------------------------

    currentGroup = upper(strtrim(string(D.Group)));

    switch currentGroup

        case "HC"

            thetaHC{end+1} = thetaMatrix;
            alphaHC{end+1} = alphaMatrix;
            betaHC{end+1}  = betaMatrix;

        case "PD"

            thetaPD{end+1} = thetaMatrix;
            alphaPD{end+1} = alphaMatrix;
            betaPD{end+1}  = betaMatrix;

        otherwise

            warning('Unknown group for %s.',string(D.Subject));

    end

end

%% ========================================================================
% 6. SAMPLE SIZE CHECK
%% ========================================================================

nHC = numel(thetaHC);
nPD = numel(thetaPD);

fprintf('\n========================================\n');
fprintf('COMPLETE-NETWORK COHORT\n');
fprintf('========================================\n');

fprintf('HC = %d\n',nHC);
fprintf('PD = %d\n',nPD);
fprintf('Total = %d\n',nHC+nPD);

% Expected from the primary cohort:
%
% HC = 22
% PD = 95
% Total = 117

%% ========================================================================
% 7. STACK SUBJECT MATRICES
%% ========================================================================

thetaHCstack = cat(3,thetaHC{:});
thetaPDstack = cat(3,thetaPD{:});

alphaHCstack = cat(3,alphaHC{:});
alphaPDstack = cat(3,alphaPD{:});

betaHCstack = cat(3,betaHC{:});
betaPDstack = cat(3,betaPD{:});

%% ========================================================================
% 8. GROUP-MEAN MATRICES
%% ========================================================================

thetaHCmean = mean(thetaHCstack,3,'omitnan');
thetaPDmean = mean(thetaPDstack,3,'omitnan');

alphaHCmean = mean(alphaHCstack,3,'omitnan');
alphaPDmean = mean(alphaPDstack,3,'omitnan');

betaHCmean = mean(betaHCstack,3,'omitnan');
betaPDmean = mean(betaPDstack,3,'omitnan');

%% ========================================================================
% 9. BAND-SPECIFIC SHARED COLOR LIMITS
%
% IMPORTANT:
%
% HC and PD use SAME limit within a band.
%
% Each limit is symmetric around zero:
%
%       [-L,+L]
%
% Thus:
%       blue  = negative Walk-Rest change
%       white = approximately zero
%       red   = positive Walk-Rest change
%% ========================================================================

thetaLimit = get_symmetric_limit( ...
    thetaHCmean, ...
    thetaPDmean);

alphaLimit = get_symmetric_limit( ...
    alphaHCmean, ...
    alphaPDmean);

betaLimit = get_symmetric_limit( ...
    betaHCmean, ...
    betaPDmean);

fprintf('\nCOLOR BAR RANGES\n');
fprintf('Theta: %.5f to %.5f\n',-thetaLimit,thetaLimit);
fprintf('Alpha: %.5f to %.5f\n',-alphaLimit,alphaLimit);
fprintf('Beta : %.5f to %.5f\n',-betaLimit,betaLimit);

%% ========================================================================
% 10. CREATE FIGURE
%% ========================================================================

fig = figure( ...
    'Color','w', ...
    'Units','pixels', ...
    'Position',[100 40 900 960]);

tl = tiledlayout( ...
    3, ...
    2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% ========================================================================
% ROW 1 — THETA
%% ========================================================================

% ------------------------------------------------------------------------
% Theta HC
% ------------------------------------------------------------------------

ax1 = nexttile;

draw_heatmap( ...
    ax1, ...
    thetaHCmean, ...
    regionLabels, ...
    thetaLimit, ...
    true, ...      % show Y labels
    false);        % show X labels

title(ax1, ...
    sprintf('HC (N = %d)',nHC), ...
    'FontName','Arial', ...
    'FontSize',12, ...
    'FontWeight','bold');

ylabel(ax1, ...
    {'Theta','4–8 Hz'}, ...
    'FontName','Arial', ...
    'FontSize',11, ...
    'FontWeight','bold');

% ------------------------------------------------------------------------
% Theta PD
% ------------------------------------------------------------------------

ax2 = nexttile;

draw_heatmap( ...
    ax2, ...
    thetaPDmean, ...
    regionLabels, ...
    thetaLimit, ...
    false, ...
    false);

title(ax2, ...
    sprintf('PD (N = %d)',nPD), ...
    'FontName','Arial', ...
    'FontSize',12, ...
    'FontWeight','bold');

% Color bar for THETA row

cb1 = colorbar(ax2,'eastoutside');

cb1.Label.String = '\Delta db-wPLI';
cb1.Label.FontName = 'Arial';
cb1.Label.FontSize = 10;

cb1.FontName = 'Arial';
cb1.FontSize = 9;

cb1.Ticks = linspace( ...
    -thetaLimit, ...
     thetaLimit, ...
     5);

%% ========================================================================
% ROW 2 — ALPHA
%% ========================================================================

% ------------------------------------------------------------------------
% Alpha HC
% ------------------------------------------------------------------------

ax3 = nexttile;

draw_heatmap( ...
    ax3, ...
    alphaHCmean, ...
    regionLabels, ...
    alphaLimit, ...
    true, ...
    false);

ylabel(ax3, ...
    {'Alpha','8–13 Hz'}, ...
    'FontName','Arial', ...
    'FontSize',11, ...
    'FontWeight','bold');

% ------------------------------------------------------------------------
% Alpha PD
% ------------------------------------------------------------------------

ax4 = nexttile;

draw_heatmap( ...
    ax4, ...
    alphaPDmean, ...
    regionLabels, ...
    alphaLimit, ...
    false, ...
    false);

% Color bar for ALPHA row

cb2 = colorbar(ax4,'eastoutside');

cb2.Label.String = '\Delta db-wPLI';
cb2.Label.FontName = 'Arial';
cb2.Label.FontSize = 10;

cb2.FontName = 'Arial';
cb2.FontSize = 9;

cb2.Ticks = linspace( ...
    -alphaLimit, ...
     alphaLimit, ...
     5);

%% ========================================================================
% ROW 3 — BETA
%
% Only this row shows X-axis region abbreviations.
%% ========================================================================

% ------------------------------------------------------------------------
% Beta HC
% ------------------------------------------------------------------------

ax5 = nexttile;

draw_heatmap( ...
    ax5, ...
    betaHCmean, ...
    regionLabels, ...
    betaLimit, ...
    true, ...
    true);

ylabel(ax5, ...
    {'Beta','13–30 Hz'}, ...
    'FontName','Arial', ...
    'FontSize',11, ...
    'FontWeight','bold');

xlabel(ax5, ...
    'Regional node', ...
    'FontName','Arial', ...
    'FontSize',10);

% ------------------------------------------------------------------------
% Beta PD
% ------------------------------------------------------------------------

ax6 = nexttile;

draw_heatmap( ...
    ax6, ...
    betaPDmean, ...
    regionLabels, ...
    betaLimit, ...
    false, ...
    true);

xlabel(ax6, ...
    'Regional node', ...
    'FontName','Arial', ...
    'FontSize',10);

% Color bar for BETA row

cb3 = colorbar(ax6,'eastoutside');

cb3.Label.String = '\Delta db-wPLI';
cb3.Label.FontName = 'Arial';
cb3.Label.FontSize = 10;

cb3.FontName = 'Arial';
cb3.FontSize = 9;

cb3.Ticks = linspace( ...
    -betaLimit, ...
     betaLimit, ...
     5);

%% ========================================================================
% 11. EXPORT
%% ========================================================================

pngFile = fullfile( ...
    rootDir, ...
    'Figure1_dbwPLI_GroupMean_Heatmaps.png');

pdfFile = fullfile( ...
    rootDir, ...
    'Figure1_dbwPLI_GroupMean_Heatmaps.pdf');

figFile = fullfile( ...
    rootDir, ...
    'Figure1_dbwPLI_GroupMean_Heatmaps.fig');

exportgraphics( ...
    fig, ...
    pngFile, ...
    'Resolution',600);

exportgraphics( ...
    fig, ...
    pdfFile, ...
    'ContentType','vector');

savefig(fig,figFile);

fprintf('\nFigure saved successfully.\n');

%% ========================================================================
% LOCAL FUNCTION 1
%
% DRAW HEATMAP
%% ========================================================================

function draw_heatmap( ...
    ax, ...
    M, ...
    regionLabels, ...
    maxLimit, ...
    showYLabels, ...
    showXLabels)

    %% ---------------------------------------------------------------
    % Special diagonal color
    %
    % This gray color is NOT part of the blue-white-red color scale.
    %% ---------------------------------------------------------------

    diagonalColor = [0.72 0.72 0.72];

    %% ---------------------------------------------------------------
    % Draw matrix
    %% ---------------------------------------------------------------

    h = imagesc(ax,M);

    axis(ax,'square');

    %% ---------------------------------------------------------------
    % Make NaN diagonal transparent.
    %
    % Transparent cells reveal the gray axes background.
    %% ---------------------------------------------------------------

    h.AlphaData = isfinite(M);

    ax.Color = diagonalColor;

    %% ---------------------------------------------------------------
    % Symmetric color range
    %% ---------------------------------------------------------------

    clim(ax,[-maxLimit maxLimit]);

    %% ---------------------------------------------------------------
    % Diverging color map
    %% ---------------------------------------------------------------

    colormap(ax,blue_white_red(256));

    %% ---------------------------------------------------------------
    % Axis orientation
    %% ---------------------------------------------------------------

    set(ax,'YDir','normal');

    %% ---------------------------------------------------------------
    % Tick positions
    %% ---------------------------------------------------------------

    ax.XTick = 1:11;
    ax.YTick = 1:11;

    %% ---------------------------------------------------------------
    % X labels
    %
    % ONLY bottom row
    %% ---------------------------------------------------------------

    if showXLabels

        ax.XTickLabel = regionLabels;
        ax.XTickLabelRotation = 45;

    else

        ax.XTickLabel = [];

    end

    %% ---------------------------------------------------------------
    % Y labels
    %
    % ONLY left column
    %% ---------------------------------------------------------------

    if showYLabels

        ax.YTickLabel = regionLabels;

    else

        ax.YTickLabel = [];

    end

    %% ---------------------------------------------------------------
    % Appearance
    %% ---------------------------------------------------------------

    ax.FontName = 'Arial';
    ax.FontSize = 9;

    ax.TickLength = [0 0];

    ax.Box = 'on';
    ax.LineWidth = 1;

    %% ---------------------------------------------------------------
    % Draw cell boundaries
    %
    % Gives appearance similar to your example.
    %% ---------------------------------------------------------------

    hold(ax,'on');

    for k = 0.5:1:11.5

        line( ...
            ax, ...
            [0.5 11.5], ...
            [k k], ...
            'Color',[0.30 0.30 0.30], ...
            'LineWidth',0.35);

        line( ...
            ax, ...
            [k k], ...
            [0.5 11.5], ...
            'Color',[0.30 0.30 0.30], ...
            'LineWidth',0.35);

    end

    hold(ax,'off');

    %% Exact matrix bounds

    xlim(ax,[0.5 11.5]);
    ylim(ax,[0.5 11.5]);

end

%% ========================================================================
% LOCAL FUNCTION 2
%
% FIND SHARED HC/PD COLOR LIMIT FOR ONE BAND
%% ========================================================================

function L = get_symmetric_limit(A,B)

    values = [A(:); B(:)];

    values = values(isfinite(values));

    L = max(abs(values));

    if isempty(L) || ~isfinite(L) || L == 0
        L = 1;
    end

end

%% ========================================================================
% LOCAL FUNCTION 3
%
% MAKE MATRIX EXACTLY SYMMETRIC
%% ========================================================================

function M = make_symmetric(M)

    n = size(M,1);

    for i = 1:n-1

        for j = i+1:n

            a = M(i,j);
            b = M(j,i);

            if isfinite(a) && isfinite(b)

                value = mean([a b]);

            elseif isfinite(a)

                value = a;

            elseif isfinite(b)

                value = b;

            else

                value = NaN;

            end

            M(i,j) = value;
            M(j,i) = value;

        end

    end

end

%% ========================================================================
% LOCAL FUNCTION 4
%
% BLUE -> WHITE -> RED
%
% Negative Delta = blue
% Zero           = white
% Positive Delta = red
%
% Gray is deliberately NOT included, so the diagonal is visually distinct.
%% ========================================================================

function cmap = blue_white_red(n)

    if nargin < 1
        n = 256;
    end

    n1 = floor(n/2);
    n2 = n - n1;

    %% Blue -> white

    blue = [0.10 0.30 0.80];

    white = [1 1 1];

    red = [0.80 0.10 0.10];

    firstHalf = [ ...
        linspace(blue(1),white(1),n1)', ...
        linspace(blue(2),white(2),n1)', ...
        linspace(blue(3),white(3),n1)'];

    secondHalf = [ ...
        linspace(white(1),red(1),n2)', ...
        linspace(white(2),red(2),n2)', ...
        linspace(white(3),red(3),n2)'];

    cmap = [firstHalf; secondHalf];

end
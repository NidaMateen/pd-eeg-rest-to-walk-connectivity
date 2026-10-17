%% ========================================================================
% FIGURE 1 — COMPLETE 11-NODE REGIONAL EEG NETWORK
%
% 11 regional nodes
% 55 unique undirected edges
%
% Designed for manuscript Figure 1:
% - darker edges
% - larger nodes
% - more prominent labels
% - publication-friendly export
%
% MATLAB R2023b compatible
%% ========================================================================

clear;
clc;
close all;

%% ========================================================================
% 1. DEFINE REGIONAL NODES
%% ========================================================================

nodeNames = { ...
    'FL', ... % Frontal Left
    'FM', ... % Frontal Midline
    'FR', ... % Frontal Right
    'CL', ... % Central Left
    'CR', ... % Central Right
    'PL', ... % Parietal Left
    'PM', ... % Parietal Midline
    'PR', ... % Parietal Right
    'TL', ... % Temporal Left
    'TR', ... % Temporal Right
    'O'  ...  % Occipital
    };

nNodes = numel(nodeNames);

assert(nNodes == 11, ...
    'Expected exactly 11 regional nodes.');

%% ========================================================================
% 2. DEFINE NODE POSITIONS
%
% Schematic scalp-like layout
%% ========================================================================

x = [ ...
    -1.05, ... % FL
     0.00, ... % FM
     1.05, ... % FR
    -0.78, ... % CL
     0.78, ... % CR
    -0.62, ... % PL
     0.00, ... % PM
     0.62, ... % PR
    -1.35, ... % TL
     1.35, ... % TR
     0.00  ... % O
    ];

y = [ ...
     1.35, ... % FL
     1.65, ... % FM
     1.35, ... % FR
     0.55, ... % CL
     0.55, ... % CR
    -0.15, ... % PL
     0.05, ... % PM
    -0.15, ... % PR
     0.35, ... % TL
     0.35, ... % TR
    -1.10  ... % O
    ];

%% ========================================================================
% 3. GENERATE ALL 55 UNIQUE UNDIRECTED EDGES
%
% For 11 nodes:
%
%       11 choose 2 = 55
%% ========================================================================

edges = nchoosek(1:nNodes, 2);

nEdges = size(edges,1);

fprintf('Number of nodes = %d\n', nNodes);
fprintf('Number of edges = %d\n', nEdges);

assert(nEdges == 55, ...
    'Network must contain exactly 55 undirected edges.');

%% ========================================================================
% 4. CREATE FIGURE
%% ========================================================================

fig = figure( ...
    'Color', 'w', ...
    'Units', 'pixels', ...
    'Position', [100 100 850 850]);

ax = axes(fig);

hold(ax,'on');

%% ========================================================================
% 5. DRAW ALL 55 EDGES
%
% Darker and thicker than previous version
%% ========================================================================

edgeColor = [0.58 0.58 0.58];
edgeWidth = 1.05;

for e = 1:nEdges

    node1 = edges(e,1);
    node2 = edges(e,2);

    line( ...
        ax, ...
        [x(node1) x(node2)], ...
        [y(node1) y(node2)], ...
        'Color', edgeColor, ...
        'LineWidth', edgeWidth);

end

%% ========================================================================
% 6. DRAW NODE CIRCLES
%
% Larger circles with thicker outline
%% ========================================================================

nodeSize = 900;

scatter( ...
    ax, ...
    x, ...
    y, ...
    nodeSize, ...
    'o', ...
    'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', [0.12 0.12 0.12], ...
    'LineWidth', 1.8);

%% ========================================================================
% 7. DRAW NODE LABELS
%
% Larger and bold
%% ========================================================================

for i = 1:nNodes

    text( ...
        ax, ...
        x(i), ...
        y(i), ...
        nodeNames{i}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', 'Arial', ...
        'FontSize', 14, ...
        'FontWeight', 'bold', ...
        'Color', [0.05 0.05 0.05]);

end

%% ========================================================================
% 8. OPTIONAL PANEL LABEL
%
% If this network is part of panel 3 in Figure 1, keep this.
% Otherwise comment it out.
%% ========================================================================

% text( ...
%     ax, ...
%     -1.62, ...
%      1.88, ...
%     '(b)', ...
%     'FontName', 'Arial', ...
%     'FontSize', 16, ...
%     'FontWeight', 'bold');

%% ========================================================================
% 9. OPTIONAL TITLE
%
% Keep only if needed as a standalone network figure.
% For insertion into your larger Figure 1, I would usually remove it.
%% ========================================================================

% title( ...
%     ax, ...
%     '11-node regional network (schematic)', ...
%     'FontName', 'Arial', ...
%     'FontSize', 16, ...
%     'FontWeight', 'bold');

%% ========================================================================
% 10. FORMAT AXES
%% ========================================================================

axis(ax,'equal');
axis(ax,'off');

xlim(ax,[-1.60 1.60]);
ylim(ax,[-1.40 1.95]);

set( ...
    ax, ...
    'Position', [0.04 0.04 0.92 0.92]);

%% ========================================================================
% 11. EXPORT
%% ========================================================================

outputPNG = 'Figure1_11Node_55Edge_Network_Dark.png';
outputPDF = 'Figure1_11Node_55Edge_Network_Dark.pdf';
outputFIG = 'Figure1_11Node_55Edge_Network_Dark.fig';

% High-resolution raster
exportgraphics( ...
    fig, ...
    outputPNG, ...
    'Resolution', 600);

% Vector PDF
exportgraphics( ...
    fig, ...
    outputPDF, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

% Editable MATLAB figure
savefig( ...
    fig, ...
    outputFIG);

fprintf('\nFigure exported successfully.\n');
fprintf('Nodes: %d\n', nNodes);
fprintf('Edges: %d\n', nEdges);

fprintf('\nSaved files:\n');
fprintf('%s\n', outputPNG);
fprintf('%s\n', outputPDF);
fprintf('%s\n', outputFIG);
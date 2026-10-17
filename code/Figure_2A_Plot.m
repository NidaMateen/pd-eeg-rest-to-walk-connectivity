%% ================================================================
% FIGURE 2A — PRIMARY THETA NBS NETWORK
%
% 32-edge significant theta NBS component
% Primary complete-network cohort: N = 117
%
% Blue edges = membership in the significant NBS component ONLY.
% Edge width/color do NOT encode effect magnitude or edge-level p-value.
%
% Regional locations are schematic and are NOT source-localized anatomy.
% ================================================================

clear;
close all;
clc;

%% ---------------------------------------------------------------
% OUTPUT FILES
% ---------------------------------------------------------------

outDir = pwd;

pdfFile = fullfile(outDir, ...
    'Figure2A_Theta_NBS_Network.pdf');

pngFile = fullfile(outDir, ...
    'Figure2A_Theta_NBS_Network_600dpi.png');

%% ---------------------------------------------------------------
% 1. REGIONAL NODES
% ---------------------------------------------------------------

nodes = { ...
    'FrontalLeft'
    'FrontalMidline'
    'FrontalRight'
    'CentralLeft'
    'CentralRight'
    'ParietalLeft'
    'ParietalMidline'
    'ParietalRight'
    'TemporalLeft'
    'TemporalRight'
    'Occipital'};

% Labels displayed inside nodes
shortLabels = { ...
    'FL'
    'FM'
    'FR'
    'CL'
    'CR'
    'PL'
    'PM'
    'PR'
    'TL'
    'TR'
    'O'};

%% ---------------------------------------------------------------
% 2. COMPACT SCHEMATIC POSITIONS
%
% Left hemisphere  = negative x
% Right hemisphere = positive x
% Midline          = x = 0
% ---------------------------------------------------------------

pos = [ ...
   -0.82,  1.25;    % FL
    0.00,  1.43;    % FM
    0.82,  1.25;    % FR

   -0.72,  0.48;    % CL
    0.72,  0.48;    % CR

   -0.70, -0.25;    % PL
    0.00, -0.10;    % PM
    0.70, -0.25;    % PR

   -1.12,  0.30;    % TL
    1.12,  0.30;    % TR

    0.00, -0.95];   % O

%% ---------------------------------------------------------------
% 3. EXACT 32-EDGE PRIMARY THETA NBS COMPONENT
% ---------------------------------------------------------------

edgeNames = { ...
    'FrontalRight',    'CentralRight';
    'ParietalMidline', 'ParietalRight';
    'FrontalLeft',     'ParietalMidline';
    'ParietalLeft',    'ParietalRight';
    'FrontalLeft',     'FrontalRight';
    'FrontalLeft',     'TemporalRight';
    'ParietalMidline', 'TemporalRight';
    'FrontalLeft',     'CentralRight';
    'CentralRight',    'ParietalMidline';
    'FrontalMidline',  'ParietalLeft';
    'FrontalRight',    'ParietalRight';
    'FrontalLeft',     'ParietalRight';
    'TemporalRight',   'Occipital';
    'FrontalLeft',     'TemporalLeft';
    'FrontalMidline',  'Occipital';
    'ParietalLeft',    'TemporalRight';
    'FrontalLeft',     'ParietalLeft';
    'FrontalRight',    'ParietalLeft';
    'ParietalRight',   'Occipital';
    'FrontalRight',    'TemporalRight';
    'FrontalMidline',  'FrontalRight';
    'FrontalLeft',     'Occipital';
    'CentralRight',    'Occipital';
    'CentralRight',    'ParietalLeft';
    'FrontalMidline',  'CentralRight';
    'FrontalRight',    'CentralLeft';
    'FrontalLeft',     'FrontalMidline';
    'FrontalRight',    'TemporalLeft';
    'ParietalLeft',    'ParietalMidline';
    'CentralLeft',     'CentralRight';
    'FrontalMidline',  'ParietalRight';
    'CentralLeft',     'ParietalRight'};

assert(size(edgeNames,1) == 32, ...
    'Expected exactly 32 NBS component edges.');

%% ---------------------------------------------------------------
% 4. CONVERT NODE NAMES TO INDICES
% ---------------------------------------------------------------

nEdges = size(edgeNames,1);
edgeIdx = zeros(nEdges,2);

for e = 1:nEdges

    edgeIdx(e,1) = find(strcmp(nodes,edgeNames{e,1}));
    edgeIdx(e,2) = find(strcmp(nodes,edgeNames{e,2}));

end

%% ---------------------------------------------------------------
% 5. VALIDATE NODE DEGREE
% ---------------------------------------------------------------

degreeCount = zeros(numel(nodes),1);

for e = 1:nEdges

    degreeCount(edgeIdx(e,1)) = ...
        degreeCount(edgeIdx(e,1)) + 1;

    degreeCount(edgeIdx(e,2)) = ...
        degreeCount(edgeIdx(e,2)) + 1;

end

fprintf('\nNode participation:\n');

for n = 1:numel(nodes)

    fprintf('%-18s %d\n', ...
        nodes{n}, degreeCount(n));

end

% Expected node degrees:
expectedDegree = [9 6 8 3 7 7 5 7 2 5 5]';

assert(isequal(degreeCount,expectedDegree), ...
    'Node-degree validation failed.');

%% ---------------------------------------------------------------
% 6. CREATE COMPACT FIGURE
% ---------------------------------------------------------------

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 11.5 10.5]);

ax = axes(fig);

hold(ax,'on');
axis(ax,'equal');
axis(ax,'off');

% Tight limits = less unused white space
xlim(ax,[-1.30 1.30]);
ylim(ax,[-1.12 1.58]);

%% ---------------------------------------------------------------
% 7. DRAW COMPONENT EDGES
% ---------------------------------------------------------------

% Muted scientific blue.
% Meaning: membership in significant theta NBS component.
edgeColor = [0.8 0.8 0.8];

for e = 1:nEdges

    i = edgeIdx(e,1);
    j = edgeIdx(e,2);

    plot(ax, ...
        [pos(i,1) pos(j,1)], ...
        [pos(i,2) pos(j,2)], ...
        '-', ...
        'Color',edgeColor, ...
        'LineWidth',1.35);

end

%% ---------------------------------------------------------------
% 8. DRAW NODES
% ---------------------------------------------------------------

nodeSize = 680;

scatter(ax, ...
    pos(:,1), ...
    pos(:,2), ...
    nodeSize, ...
    'o', ...
    'MarkerFaceColor',[1 1 1], ...
    'MarkerEdgeColor',[0.08 0.08 0.08], ...
    'LineWidth',1.6);

%% ---------------------------------------------------------------
% 9. PROMINENT NODE LABELS
% ---------------------------------------------------------------

for n = 1:numel(nodes)

    text(ax, ...
        pos(n,1), ...
        pos(n,2), ...
        shortLabels{n}, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontName','Arial', ...
        'FontSize',12, ...
        'FontWeight','bold', ...
        'Color',[0.05 0.05 0.05]);

end

%% ---------------------------------------------------------------
% 10. MINIMIZE AXIS MARGINS
% ---------------------------------------------------------------

ax.Position = [0.02 0.02 0.96 0.96];

%% ---------------------------------------------------------------
% 11. EXPORT
% ---------------------------------------------------------------

exportgraphics(ax, ...
    pdfFile, ...
    'ContentType','vector', ...
    'BackgroundColor','white');

exportgraphics(ax, ...
    pngFile, ...
    'Resolution',600, ...
    'BackgroundColor','white');

fprintf('\nSaved:\n');
fprintf('%s\n',pdfFile);
fprintf('%s\n',pngFile);
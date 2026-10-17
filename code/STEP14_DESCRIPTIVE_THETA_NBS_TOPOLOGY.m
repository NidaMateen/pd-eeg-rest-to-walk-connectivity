%% =========================================================================
% PAPER 3 — STEP 14
% DESCRIPTIVE TOPOLOGY OF THE SIGNIFICANT PRIMARY THETA NBS COMPONENT
%
% File:
%   STEP14_DESCRIPTIVE_THETA_NBS_TOPOLOGY.m
%
% PURPOSE
% -------
% Describe the already-established 32-edge theta NBS component without
% introducing any new inferential testing.
%
% This step reports:
%   1) exact component edge membership
%   2) node degree / participation
%   3) lobe-pair distribution
%   4) broad hemispheric/anatomical distribution
%   5) homologous L-R edge membership
%   6) component density
%   7) adjusted beta and v-statistic ranges
%   8) optional linkage to Step-11 edge-level q-values
%
% IMPORTANT
% ---------
% This is DESCRIPTIVE topology only.
% No new p-values are generated.
% No nodes/edges are selected using a second threshold.
% The 32-edge component is taken exactly from the frozen primary NBS result.
%% =========================================================================

clear;
clc;
close all;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% =========================================================================
% 1. PATHS
%% =========================================================================

datasetRoot = ...
    cfg.datasetRoot;

primaryNBSRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_PRIMARY_NBS');

edgeLocRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_EDGE_LOCALIZATION');

primaryNBSFile = fullfile( ...
    primaryNBSRoot, ...
    'Paper3_PRIMARY_NBS_RESULTS.mat');

edgeLocFile = fullfile( ...
    edgeLocRoot, ...
    'Paper3_EDGE_LOCALIZATION_RESULTS.mat');

outRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_THETA_NBS_TOPOLOGY');

if ~exist(outRoot,'dir')
    mkdir(outRoot);
end

if ~exist(primaryNBSFile,'file')
    error('Missing primary NBS result:\n%s',primaryNBSFile);
end

%% =========================================================================
% 2. LOAD PRIMARY NBS RESULT
%% =========================================================================

R = load(primaryNBSFile);

requiredFields = { ...
    'regionNames', ...
    'edgeNames', ...
    'edgeRegion1', ...
    'edgeRegion2', ...
    'edgeIndex1', ...
    'edgeIndex2', ...
    'betaGroup', ...
    'vObserved', ...
    'observedComponents', ...
    'ComponentTable'};

for k = 1:numel(requiredFields)
    if ~isfield(R,requiredFields{k})
        error('Primary NBS MAT lacks field: %s',requiredFields{k});
    end
end

regionNames = string(R.regionNames(:));

edgeNames   = string(R.edgeNames(:));
edgeRegion1 = string(R.edgeRegion1(:));
edgeRegion2 = string(R.edgeRegion2(:));

edgeIndex1 = double(R.edgeIndex1(:));
edgeIndex2 = double(R.edgeIndex2(:));

betaGroup = double(R.betaGroup);
vObserved = double(R.vObserved);

CT = R.ComponentTable;

if numel(regionNames) ~= 11
    error('Expected 11 regional nodes.');
end

if numel(edgeNames) ~= 55
    error('Expected 55 regional edges.');
end

%% =========================================================================
% 3. LOCK THE SIGNIFICANT THETA-NEGATIVE COMPONENT
%% =========================================================================

sigMask = ...
    CT.Significant_FWER_0p05 & ...
    lower(string(CT.Band))=="theta" & ...
    contains(lower(string(CT.Direction)),"negative");

if sum(sigMask) ~= 1
    error('Expected exactly one significant theta-negative primary component.');
end

sigRow = CT(sigMask,:);

componentID = double(sigRow.ComponentID);
componentExtent = double(sigRow.ExtentEdges);
componentP = double(sigRow.JointFWER_P);

thetaNegativeComponents = R.observedComponents{1,2};

if componentID < 1 || componentID > numel(thetaNegativeComponents)
    error('Invalid component ID.');
end

componentEdges = double( ...
    thetaNegativeComponents(componentID).edgeIdx(:));

componentNodes = double( ...
    thetaNegativeComponents(componentID).nodeIdx(:));

if numel(componentEdges) ~= componentExtent
    error('Component edge count does not match ComponentTable.');
end

if componentExtent ~= 32
    error('Expected frozen primary theta component extent = 32.');
end

fprintf('\n============================================================\n');
fprintf('PRIMARY THETA COMPONENT LOCK\n');
fprintf('============================================================\n');

fprintf('Component edges = %d / 55\n',componentExtent);
fprintf('Component nodes = %d / 11\n',numel(componentNodes));
fprintf('Joint FWER p    = %.8f\n',componentP);

%% =========================================================================
% 4. OPTIONAL LOAD STEP-11 EDGE LOCALIZATION
%% =========================================================================

haveEdgeLoc = exist(edgeLocFile,'file')==2;

thetaEdgeLoc = table();

if haveEdgeLoc

    E = load(edgeLocFile);

    if isfield(E,'EdgeResults')

        thetaEdgeLoc = E.EdgeResults( ...
            lower(string(E.EdgeResults.Band))=="theta",:);

        if height(thetaEdgeLoc) ~= 55
            warning('Step-11 theta table does not contain 55 rows; q-values omitted.');
            haveEdgeLoc = false;
        end

    else

        warning('Step-11 result MAT lacks EdgeResults; q-values omitted.');
        haveEdgeLoc = false;

    end

end

%% =========================================================================
% 5. DEFINE FIXED DESCRIPTIVE ANATOMICAL LABELS
%% =========================================================================

% Hemisphere / laterality label for each node.
nodeHemisphere = strings(11,1);

% Lobe/system label for each node.
nodeLobe = strings(11,1);

for n = 1:11

    rn = regionNames(n);

    if contains(rn,"Left")
        nodeHemisphere(n) = "Left";
    elseif contains(rn,"Right")
        nodeHemisphere(n) = "Right";
    elseif contains(rn,"Midline")
        nodeHemisphere(n) = "Midline";
    elseif rn=="Occipital"
        nodeHemisphere(n) = "BilateralCombined";
    else
        error('Unknown hemisphere label for node: %s',rn);
    end

    if startsWith(rn,"Frontal")
        nodeLobe(n) = "Frontal";
    elseif startsWith(rn,"Central")
        nodeLobe(n) = "Central";
    elseif startsWith(rn,"Parietal")
        nodeLobe(n) = "Parietal";
    elseif startsWith(rn,"Temporal")
        nodeLobe(n) = "Temporal";
    elseif rn=="Occipital"
        nodeLobe(n) = "Occipital";
    else
        error('Unknown lobe label for node: %s',rn);
    end

end

%% =========================================================================
% 6. BUILD EXACT 32-EDGE COMPONENT TABLE
%% =========================================================================

nComp = numel(componentEdges);

EdgeIndex = componentEdges;
Edge = edgeNames(componentEdges);

Region1 = edgeRegion1(componentEdges);
Region2 = edgeRegion2(componentEdges);

Node1 = edgeIndex1(componentEdges);
Node2 = edgeIndex2(componentEdges);

Beta_PDminusHC_Adjusted = betaGroup(componentEdges,1);
VStatistic = vObserved(componentEdges,1);

AbsBeta = abs(Beta_PDminusHC_Adjusted);
AbsV = abs(VStatistic);

Lobe1 = nodeLobe(Node1);
Lobe2 = nodeLobe(Node2);

Hemisphere1 = nodeHemisphere(Node1);
Hemisphere2 = nodeHemisphere(Node2);

LobePair = strings(nComp,1);
HemisphereClass = strings(nComp,1);

InvolvesMidline = false(nComp,1);
InvolvesOccipital = false(nComp,1);
IsHomologousLR = false(nComp,1);

for e = 1:nComp

    LobePair(e) = canonical_pair_label(Lobe1(e),Lobe2(e));

    InvolvesMidline(e) = ...
        Hemisphere1(e)=="Midline" || Hemisphere2(e)=="Midline";

    InvolvesOccipital(e) = ...
        Region1(e)=="Occipital" || Region2(e)=="Occipital";

    HemisphereClass(e) = classify_edge_anatomy( ...
        Hemisphere1(e), ...
        Hemisphere2(e), ...
        InvolvesMidline(e), ...
        InvolvesOccipital(e));

    IsHomologousLR(e) = is_homologous_lr( ...
        Region1(e),Region2(e));

end

% Rank by absolute v and absolute beta for descriptive display.
[~,ordV] = sort(AbsV,'descend');
RankByAbsV = nan(nComp,1);
RankByAbsV(ordV) = (1:nComp)';

[~,ordB] = sort(AbsBeta,'descend');
RankByAbsBeta = nan(nComp,1);
RankByAbsBeta(ordB) = (1:nComp)';

ComponentEdgeTable = table( ...
    EdgeIndex, ...
    Edge, ...
    Region1, ...
    Region2, ...
    Lobe1, ...
    Lobe2, ...
    LobePair, ...
    Hemisphere1, ...
    Hemisphere2, ...
    HemisphereClass, ...
    InvolvesMidline, ...
    InvolvesOccipital, ...
    IsHomologousLR, ...
    Beta_PDminusHC_Adjusted, ...
    VStatistic, ...
    AbsBeta, ...
    AbsV, ...
    RankByAbsBeta, ...
    RankByAbsV);

% Add Step-11 available-case edge localization statistics if available.
if haveEdgeLoc

    [tf,loc] = ismember(Edge, string(thetaEdgeLoc.Edge));

    if all(tf)

        ComponentEdgeTable.N_HC_AvailableCase = ...
            double(thetaEdgeLoc.N_HC(loc));

        ComponentEdgeTable.N_PD_AvailableCase = ...
            double(thetaEdgeLoc.N_PD(loc));

        ComponentEdgeTable.Edge_P_Raw = ...
            double(thetaEdgeLoc.P_Raw(loc));

        ComponentEdgeTable.Edge_Q_BH_FDR = ...
            double(thetaEdgeLoc.Q_BH_FDR(loc));

        ComponentEdgeTable.Edge_FDR_Significant = ...
            logical(thetaEdgeLoc.Significant_BH_FDR_0p05(loc));

        ComponentEdgeTable.Edge_Hedges_g_Adjusted = ...
            double(thetaEdgeLoc.HedgesType_gAdjusted(loc));

    else

        warning('Could not match all 32 component edges to Step-11 table.');

    end

end

%% =========================================================================
% 7. NODE DEGREE / PARTICIPATION
%% =========================================================================

Node = regionNames;
Lobe = nodeLobe;
Hemisphere = nodeHemisphere;

ComponentDegree = zeros(11,1);
PossibleDegree = repmat(10,11,1);

IncidentBetaMean = nan(11,1);
IncidentBetaMedian = nan(11,1);
IncidentAbsBetaMean = nan(11,1);

IncidentVMean = nan(11,1);
IncidentAbsVMean = nan(11,1);

for n = 1:11

    inc = ...
        Node1==n | Node2==n;

    ComponentDegree(n) = sum(inc);

    if any(inc)

        IncidentBetaMean(n) = mean( ...
            Beta_PDminusHC_Adjusted(inc));

        IncidentBetaMedian(n) = median( ...
            Beta_PDminusHC_Adjusted(inc));

        IncidentAbsBetaMean(n) = mean( ...
            AbsBeta(inc));

        IncidentVMean(n) = mean( ...
            VStatistic(inc));

        IncidentAbsVMean(n) = mean( ...
            AbsV(inc));

    end

end

DegreeFractionOfPossible = ...
    ComponentDegree ./ PossibleDegree;

NodeTable = table( ...
    Node, ...
    Lobe, ...
    Hemisphere, ...
    ComponentDegree, ...
    PossibleDegree, ...
    DegreeFractionOfPossible, ...
    IncidentBetaMean, ...
    IncidentBetaMedian, ...
    IncidentAbsBetaMean, ...
    IncidentVMean, ...
    IncidentAbsVMean);

NodeTable = sortrows( ...
    NodeTable, ...
    {'ComponentDegree','IncidentAbsVMean'}, ...
    {'descend','descend'});

%% =========================================================================
% 8. LOBE-PAIR DISTRIBUTION
%% =========================================================================

lobeCats = unique(LobePair,'stable');

LobePairCategory = lobeCats;
EdgeCount = nan(numel(lobeCats),1);
PercentOfComponent = nan(numel(lobeCats),1);

MeanAdjustedBeta = nan(numel(lobeCats),1);
MedianAdjustedBeta = nan(numel(lobeCats),1);

for i = 1:numel(lobeCats)

    idx = LobePair==lobeCats(i);

    EdgeCount(i) = sum(idx);

    PercentOfComponent(i) = ...
        100*sum(idx)/nComp;

    MeanAdjustedBeta(i) = ...
        mean(Beta_PDminusHC_Adjusted(idx));

    MedianAdjustedBeta(i) = ...
        median(Beta_PDminusHC_Adjusted(idx));

end

LobeDistributionTable = table( ...
    LobePairCategory, ...
    EdgeCount, ...
    PercentOfComponent, ...
    MeanAdjustedBeta, ...
    MedianAdjustedBeta);

LobeDistributionTable = sortrows( ...
    LobeDistributionTable, ...
    'EdgeCount', ...
    'descend');

%% =========================================================================
% 9. HEMISPHERIC / ANATOMICAL DISTRIBUTION
%% =========================================================================

hemCats = [ ...
    "Left_intrahemispheric"
    "Right_intrahemispheric"
    "Interhemispheric_lateralized"
    "Midline_involving"
    "Occipital_involving"];

AnatomicalClass = hemCats;
EdgeCountH = zeros(numel(hemCats),1);
PercentOfComponentH = zeros(numel(hemCats),1);
MeanAdjustedBetaH = nan(numel(hemCats),1);

for i = 1:numel(hemCats)

    idx = HemisphereClass==hemCats(i);

    EdgeCountH(i) = sum(idx);

    PercentOfComponentH(i) = ...
        100*sum(idx)/nComp;

    if any(idx)
        MeanAdjustedBetaH(i) = ...
            mean(Beta_PDminusHC_Adjusted(idx));
    end

end

HemisphereDistributionTable = table( ...
    AnatomicalClass, ...
    EdgeCountH, ...
    PercentOfComponentH, ...
    MeanAdjustedBetaH, ...
    'VariableNames', ...
    {'AnatomicalClass', ...
     'EdgeCount', ...
     'PercentOfComponent', ...
     'MeanAdjustedBeta_PDminusHC'});

%% =========================================================================
% 10. HOMOLOGOUS LEFT-RIGHT MEMBERSHIP
%% =========================================================================

HomologousPair = [ ...
    "Frontal_LR"
    "Central_LR"
    "Parietal_LR"
    "Temporal_LR"];

HomRegion1 = [ ...
    "FrontalLeft"
    "CentralLeft"
    "ParietalLeft"
    "TemporalLeft"];

HomRegion2 = [ ...
    "FrontalRight"
    "CentralRight"
    "ParietalRight"
    "TemporalRight"];

IncludedInPrimaryComponent = false(4,1);
AdjustedBeta = nan(4,1);
V = nan(4,1);

for h = 1:4

    m = ...
        ((Region1==HomRegion1(h)) & (Region2==HomRegion2(h))) | ...
        ((Region1==HomRegion2(h)) & (Region2==HomRegion1(h)));

    IncludedInPrimaryComponent(h) = any(m);

    if any(m)

        AdjustedBeta(h) = ...
            Beta_PDminusHC_Adjusted(find(m,1));

        V(h) = ...
            VStatistic(find(m,1));

    end

end

HomologousMembershipTable = table( ...
    HomologousPair, ...
    HomRegion1, ...
    HomRegion2, ...
    IncludedInPrimaryComponent, ...
    AdjustedBeta, ...
    V);

%% =========================================================================
% 11. COMPONENT MATRICES
%% =========================================================================

BinaryAdjacency = zeros(11,11);
BetaMatrix = nan(11,11);
VMatrix = nan(11,11);

for k = 1:nComp

    i = Node1(k);
    j = Node2(k);

    BinaryAdjacency(i,j) = 1;
    BinaryAdjacency(j,i) = 1;

    BetaMatrix(i,j) = Beta_PDminusHC_Adjusted(k);
    BetaMatrix(j,i) = Beta_PDminusHC_Adjusted(k);

    VMatrix(i,j) = VStatistic(k);
    VMatrix(j,i) = VStatistic(k);

end

BinaryAdjacencyTable = array2table( ...
    BinaryAdjacency, ...
    'VariableNames',cellstr(regionNames), ...
    'RowNames',cellstr(regionNames));

BetaMatrixTable = array2table( ...
    BetaMatrix, ...
    'VariableNames',cellstr(regionNames), ...
    'RowNames',cellstr(regionNames));

VMatrixTable = array2table( ...
    VMatrix, ...
    'VariableNames',cellstr(regionNames), ...
    'RowNames',cellstr(regionNames));

%% =========================================================================
% 12. COMPONENT-LEVEL DESCRIPTIVE SUMMARY
%% =========================================================================

ComponentDensityPercent = ...
    100*nComp/55;

NNodesInComponent = numel(unique([Node1;Node2]));

NAllNegativeBeta = sum(Beta_PDminusHC_Adjusted<0);

BetaMin = min(Beta_PDminusHC_Adjusted);
BetaMedian = median(Beta_PDminusHC_Adjusted);
BetaMax = max(Beta_PDminusHC_Adjusted);

VMin = min(VStatistic);
VMedian = median(VStatistic);
VMax = max(VStatistic);

NMidlineEdges = sum(InvolvesMidline);
NOccipitalEdges = sum(InvolvesOccipital);
NHomologousIncluded = sum(IsHomologousLR);

NComponentEdgesFDR = NaN;

if ismember('Edge_FDR_Significant', ...
        ComponentEdgeTable.Properties.VariableNames)

    NComponentEdgesFDR = ...
        sum(ComponentEdgeTable.Edge_FDR_Significant);

end

SummaryMetric = [ ...
    "Component_edges"
    "Possible_11node_edges"
    "Component_density_percent"
    "Nodes_in_component"
    "Possible_nodes"
    "Edges_with_negative_adjusted_beta"
    "Midline_involving_edges"
    "Occipital_involving_edges"
    "Homologous_LR_edges_in_component"
    "Component_edges_FDR_significant_in_Step11"
    "Adjusted_beta_min"
    "Adjusted_beta_median"
    "Adjusted_beta_max"
    "V_min"
    "V_median"
    "V_max"
    "Primary_joint_FWER_p"];

SummaryValue = [ ...
    nComp
    55
    ComponentDensityPercent
    NNodesInComponent
    11
    NAllNegativeBeta
    NMidlineEdges
    NOccipitalEdges
    NHomologousIncluded
    NComponentEdgesFDR
    BetaMin
    BetaMedian
    BetaMax
    VMin
    VMedian
    VMax
    componentP];

SummaryTable = table( ...
    SummaryMetric, ...
    SummaryValue);

%% =========================================================================
% 13. PRINT CONSOLE SUMMARY
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('DESCRIPTIVE THETA NBS TOPOLOGY\n');
fprintf('============================================================\n');

disp(SummaryTable);

fprintf('\nNode participation, sorted by component degree:\n');
disp(NodeTable);

fprintf('\nLobe-pair distribution:\n');
disp(LobeDistributionTable);

fprintf('\nBroad anatomical distribution:\n');
disp(HemisphereDistributionTable);

fprintf('\nHomologous L-R edges in the primary component:\n');
disp(HomologousMembershipTable);

fprintf('\nExact 32-edge component membership, sorted by |v|:\n');

ComponentEdgeTableForDisplay = sortrows( ...
    ComponentEdgeTable, ...
    'AbsV', ...
    'descend');

displayVars = { ...
    'Edge', ...
    'LobePair', ...
    'HemisphereClass', ...
    'Beta_PDminusHC_Adjusted', ...
    'VStatistic', ...
    'IsHomologousLR'};

if ismember('Edge_Q_BH_FDR', ...
        ComponentEdgeTable.Properties.VariableNames)
    displayVars = [displayVars,{'Edge_Q_BH_FDR'}];
end

disp(ComponentEdgeTableForDisplay(:,displayVars));

%% =========================================================================
% 14. SAVE TABLES
%% =========================================================================

edgeFile = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_component_32_edges.csv');

nodeFile = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_node_degree.csv');

lobeFile = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_lobe_pair_distribution.csv');

hemiFile = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_anatomical_distribution.csv');

homFile = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_homologous_membership.csv');

summaryFile = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_topology_summary.csv');

adjFile = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_binary_adjacency.csv');

betaMatrixFile = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_beta_matrix.csv');

vMatrixFile = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_v_matrix.csv');

writetable(ComponentEdgeTable,edgeFile);
writetable(NodeTable,nodeFile);
writetable(LobeDistributionTable,lobeFile);
writetable(HemisphereDistributionTable,hemiFile);
writetable(HomologousMembershipTable,homFile);
writetable(SummaryTable,summaryFile);

writetable( ...
    BinaryAdjacencyTable, ...
    adjFile, ...
    'WriteRowNames',true);

writetable( ...
    BetaMatrixTable, ...
    betaMatrixFile, ...
    'WriteRowNames',true);

writetable( ...
    VMatrixTable, ...
    vMatrixFile, ...
    'WriteRowNames',true);

%% =========================================================================
% 15. SAVE MAT
%% =========================================================================

resultFile = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_TOPOLOGY_RESULTS.mat');

T = struct();

T.componentEdgeIdx = componentEdges;
T.componentNodeIdx = componentNodes;

T.ComponentEdgeTable = ComponentEdgeTable;
T.NodeTable = NodeTable;
T.LobeDistributionTable = LobeDistributionTable;
T.HemisphereDistributionTable = HemisphereDistributionTable;
T.HomologousMembershipTable = HomologousMembershipTable;
T.SummaryTable = SummaryTable;

T.BinaryAdjacency = BinaryAdjacency;
T.BetaMatrix = BetaMatrix;
T.VMatrix = VMatrix;

save(resultFile,'-struct','T','-v7.3');

%% =========================================================================
% 16. FIGURE 1 — COMPONENT BETA MATRIX
%% =========================================================================

fig1 = figure( ...
    'Color','w', ...
    'Position',[100 100 1000 850]);

imagesc(BetaMatrix);

axis square;

cb = colorbar;
ylabel(cb,'Adjusted PD-HC beta for Walk-Rest db-wPLI');

set(gca, ...
    'XTick',1:11, ...
    'XTickLabel',regionNames, ...
    'XTickLabelRotation',45, ...
    'YTick',1:11, ...
    'YTickLabel',regionNames, ...
    'FontSize',10);

title({ ...
    'Primary theta NBS component', ...
    'Adjusted PD-HC beta on the 32 component edges'});

xlabel('Region');
ylabel('Region');

betaFigPng = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_component_beta_matrix.png');

betaFigPdf = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_component_beta_matrix.pdf');

exportgraphics(fig1,betaFigPng,'Resolution',300);
exportgraphics(fig1,betaFigPdf,'ContentType','vector');

close(fig1);

%% =========================================================================
% 17. FIGURE 2 — NODE DEGREE
%% =========================================================================

% Restore original anatomical node order for plotting.
[~,nodeOrder] = ismember(regionNames,NodeTable.Node);
degreeOriginalOrder = NodeTable.ComponentDegree(nodeOrder);

fig2 = figure( ...
    'Color','w', ...
    'Position',[100 100 1100 650]);

bar(degreeOriginalOrder);

set(gca, ...
    'XTick',1:11, ...
    'XTickLabel',regionNames, ...
    'XTickLabelRotation',45, ...
    'FontSize',10);

ylabel('Number of incident component edges');
xlabel('Regional node');

title('Node participation in the 32-edge theta NBS component');

ylim([0 10.5]);

degreeFigPng = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_node_degree.png');

degreeFigPdf = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_node_degree.pdf');

exportgraphics(fig2,degreeFigPng,'Resolution',300);
exportgraphics(fig2,degreeFigPdf,'ContentType','vector');

close(fig2);

%% =========================================================================
% 18. TEXT SUMMARY FOR MANUSCRIPT PREPARATION
%% =========================================================================

textFile = fullfile( ...
    outRoot, ...
    'Paper3_THETA_NBS_topology_manuscript_notes.txt');

fid = fopen(textFile,'w');

if fid==-1
    error('Could not open manuscript-notes text file.');
end

fprintf(fid,'PAPER 3 — DESCRIPTIVE TOPOLOGY OF PRIMARY THETA NBS COMPONENT\n\n');

fprintf(fid,'Primary result:\n');
fprintf(fid,'- Theta negative PD<HC Rest-to-Walk component\n');
fprintf(fid,'- Extent = %d edges\n',nComp);
fprintf(fid,'- Joint FWER p = %.8f\n\n',componentP);

fprintf(fid,'Topology:\n');
fprintf(fid,'- Nodes represented = %d/11\n',NNodesInComponent);
fprintf(fid,'- Component density = %.1f%% of all 55 possible regional edges\n', ...
    ComponentDensityPercent);

fprintf(fid,'- Negative adjusted beta edges = %d/%d\n', ...
    NAllNegativeBeta,nComp);

fprintf(fid,'- Midline-involving edges = %d\n',NMidlineEdges);
fprintf(fid,'- Occipital-involving edges = %d\n',NOccipitalEdges);
fprintf(fid,'- Homologous L-R edges represented = %d/4\n\n',NHomologousIncluded);

fprintf(fid,'Adjusted beta across component edges:\n');
fprintf(fid,'- min = %.6f\n',BetaMin);
fprintf(fid,'- median = %.6f\n',BetaMedian);
fprintf(fid,'- max = %.6f\n\n',BetaMax);

fprintf(fid,'Highest-degree nodes:\n');

topN = min(5,height(NodeTable));

for i = 1:topN
    fprintf(fid,'- %s: degree %d/10\n', ...
        NodeTable.Node(i), ...
        NodeTable.ComponentDegree(i));
end

fprintf(fid,'\nBroad anatomical edge-class counts:\n');

for i = 1:height(HemisphereDistributionTable)

    fprintf(fid,'- %s: %d edges (%.1f%%)\n', ...
        HemisphereDistributionTable.AnatomicalClass(i), ...
        HemisphereDistributionTable.EdgeCount(i), ...
        HemisphereDistributionTable.PercentOfComponent(i));

end

fprintf(fid,'\nLargest lobe-pair classes:\n');

topL = min(8,height(LobeDistributionTable));

for i = 1:topL

    fprintf(fid,'- %s: %d edges (%.1f%%)\n', ...
        LobeDistributionTable.LobePairCategory(i), ...
        LobeDistributionTable.EdgeCount(i), ...
        LobeDistributionTable.PercentOfComponent(i));

end

fprintf(fid,'\nInterpretive guardrail:\n');
fprintf(fid,[ ...
    'These topology summaries are descriptive characteristics of the already-\n' ...
    'significant NBS component. They do not constitute independent inferential\n' ...
    'tests, and node degree / lobe counts should not be assigned p-values.\n']);

fclose(fid);

%% =========================================================================
% 19. FINAL CONSOLE MESSAGE
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('STEP 14 DESCRIPTIVE TOPOLOGY COMPLETE\n');
fprintf('============================================================\n');

fprintf('Saved:\n');
fprintf('32 edges            : %s\n',edgeFile);
fprintf('Node degree         : %s\n',nodeFile);
fprintf('Lobe distribution   : %s\n',lobeFile);
fprintf('Anatomical classes  : %s\n',hemiFile);
fprintf('Homologous membership: %s\n',homFile);
fprintf('Summary             : %s\n',summaryFile);
fprintf('Adjacency matrix    : %s\n',adjFile);
fprintf('Beta matrix         : %s\n',betaMatrixFile);
fprintf('V matrix            : %s\n',vMatrixFile);
fprintf('Full MAT            : %s\n',resultFile);
fprintf('Manuscript notes    : %s\n',textFile);
fprintf('Beta matrix figure  : %s\n',betaFigPng);
fprintf('Node-degree figure  : %s\n',degreeFigPng);

fprintf('\nNo new inferential tests were performed.\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function label = canonical_pair_label(a,b)

    order = [ ...
        "Frontal"
        "Central"
        "Parietal"
        "Temporal"
        "Occipital"];

    ia = find(order==a,1);
    ib = find(order==b,1);

    if isempty(ia) || isempty(ib)
        error('Unknown lobe label.');
    end

    if ia<=ib
        label = a + "__" + b;
    else
        label = b + "__" + a;
    end

end

%% =========================================================================

function label = classify_edge_anatomy( ...
    h1,h2,involvesMidline,involvesOccipital)

    if involvesOccipital
        label = "Occipital_involving";
        return;
    end

    if involvesMidline
        label = "Midline_involving";
        return;
    end

    if h1=="Left" && h2=="Left"
        label = "Left_intrahemispheric";
        return;
    end

    if h1=="Right" && h2=="Right"
        label = "Right_intrahemispheric";
        return;
    end

    if (h1=="Left" && h2=="Right") || ...
       (h1=="Right" && h2=="Left")

        label = "Interhemispheric_lateralized";
        return;
    end

    error('Unclassified anatomical edge type: %s / %s',h1,h2);

end

%% =========================================================================

function tf = is_homologous_lr(r1,r2)

    pairs = [ ...
        "FrontalLeft","FrontalRight"
        "CentralLeft","CentralRight"
        "ParietalLeft","ParietalRight"
        "TemporalLeft","TemporalRight"];

    tf = false;

    for i = 1:size(pairs,1)

        if (r1==pairs(i,1) && r2==pairs(i,2)) || ...
           (r1==pairs(i,2) && r2==pairs(i,1))

            tf = true;
            return;

        end

    end

end

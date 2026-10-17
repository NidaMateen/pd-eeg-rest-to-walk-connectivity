%% =========================================================================
% PAPER 3 — STEP 13
% PRESPECIFIED ARTIFACT / SENSITIVITY / ROBUSTNESS ANALYSIS
%
% File:
%   STEP13_ARTIFACT_SENSITIVITY_QC.m
%
% PURPOSE
% -------
% Stress-test the primary theta Rest->Walk reconfiguration finding against
% preprocessing burden, residual high-frequency artifact, regional
% missingness, and cohort-definition sensitivity.
%
% This script is a SENSITIVITY/QC step. It does not redefine the primary
% NBS result or create a new primary endpoint.
%
% FROZEN PRIMARY RESULT BEING STRESS-TESTED
% ------------------------------------------
% Complete-network cohort:
%   N = 117 = 22 HC + 95 PD
%
% Primary outcome:
%   Delta db-wPLI = Walk - Rest
%
% Primary significant component:
%   theta, negative PD<HC change, 32 edges
%
% QC / SENSITIVITY QUESTIONS
% --------------------------
% A. Is the subject-level primary theta component signal associated with:
%       - NBadChannels
%       - PctRejectedIC
%       - HighPowerRatio = Walk/Rest 30-45 Hz global power
%
% B. Does the adjusted PD-HC effect remain after:
%       - excluding subjects with >5% bad EEG channels
%         (NBadChannels >= 5 bad scalp channels; frozen protocol rule)
%       - excluding subjects with >20% rejected ICs
%       - applying both restrictions
%
% C. Does the 30-45 Hz artifact-control connectivity band reproduce the
%    biological theta pattern?
%       - same 32 theta-component edges, but high-band Delta
%       - whole 55-edge high-band mean
%       - 55 high-band edge models + BH-FDR
%       - edge-effect topology similarity theta vs high
%       - high-band NBS using the SAME |v|>2.5 threshold and 10,000
%         Freedman-Lane/Rademacher resamples, corrected across high +/- only
%
% D. Does complete-network restriction materially change theta edge effects?
%       - compare 55 theta edge beta estimates using:
%           1) complete N=117
%           2) available-case per edge
%
% INTERPRETATION
% --------------
% Reassuring:
%   - primary theta Group beta remains materially similar with QC adjustment
%   - exclusion sensitivities preserve direction/magnitude
%   - high-band artifact control does NOT reproduce widespread theta topology
%   - complete-vs-available theta edge effects are highly concordant
%
% Concerning:
%   - primary theta signal strongly tracks HighPowerRatio / ICA burden
%   - Group effect collapses after reasonable QC adjustment/restriction
%   - 30-45 Hz shows a comparable widespread negative topology or NBS network
%
% IMPORTANT
% ---------
% Do not alter thresholds after seeing results.
%% =========================================================================

clear;
clc;
close all;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% =========================================================================
% 1. FROZEN SETTINGS
%% =========================================================================

EXPECTED_COMPLETE_N  = 117;
EXPECTED_COMPLETE_HC = 22;
EXPECTED_COMPLETE_PD = 95;

EXPECTED_FULL_N  = 130;
EXPECTED_FULL_HC = 23;
EXPECTED_FULL_PD = 107;

NBS_THRESHOLD = 2.5;
N_RESAMPLES = 10000;
RNG_SEED = 20260916;

FDR_ALPHA = 0.05;

% High bad-channel burden:
% Frozen written protocol: exclude subjects with >=5 bad channels.
BAD_CHANNEL_EXCLUDE_THRESHOLD = 5;

% This threshold was already part of the frozen ICA manual-QC safety rule.
HIGH_IC_REJECTION_THRESHOLD_PCT = 20;

%% =========================================================================
% 2. PATHS
%% =========================================================================

datasetRoot = ...
    cfg.datasetRoot;

connectivityRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_FINAL_CONNECTIVITY');

primaryNBSRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_PRIMARY_NBS');

statsInputFile = fullfile( ...
    connectivityRoot, ...
    'Paper3_CONNECTIVITY_STATS_INPUT.mat');

subjectSummaryFile = fullfile( ...
    connectivityRoot, ...
    'Paper3_CONNECTIVITY_subject_summary.csv');

participantsFile = fullfile( ...
    datasetRoot, ...
    'participants.tsv');

primaryNBSFile = fullfile( ...
    primaryNBSRoot, ...
    'Paper3_PRIMARY_NBS_RESULTS.mat');

outRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_ARTIFACT_SENSITIVITY');

if ~exist(outRoot,'dir')
    mkdir(outRoot);
end

filesRequired = { ...
    statsInputFile, ...
    subjectSummaryFile, ...
    participantsFile, ...
    primaryNBSFile};

for k = 1:numel(filesRequired)
    if ~exist(filesRequired{k},'file')
        error('Missing required file:\n%s',filesRequired{k});
    end
end

%% =========================================================================
% 3. LOAD CONNECTIVITY + PRIMARY NBS
%% =========================================================================

S = load(statsInputFile);
R = load(primaryNBSFile);

requiredS = { ...
    'subjects','groups','edgeNames','edgeRegion1','edgeRegion2', ...
    'edgeIndex1','edgeIndex2','bandNames','deltaBio','highDelta', ...
    'completeNetworkMask'};

for k = 1:numel(requiredS)
    if ~isfield(S,requiredS{k})
        error('Statistics input lacks field: %s',requiredS{k});
    end
end

subjectsAll = string(S.subjects(:));
groupsAll   = string(S.groups(:));

edgeNames   = string(S.edgeNames(:));
edgeRegion1 = string(S.edgeRegion1(:));
edgeRegion2 = string(S.edgeRegion2(:));

edgeIndex1 = double(S.edgeIndex1(:));
edgeIndex2 = double(S.edgeIndex2(:));

bandNames = lower(string(S.bandNames(:)));

deltaBio = double(S.deltaBio);
highDelta = double(S.highDelta);

completeMask = logical(S.completeNetworkMask(:));

if numel(subjectsAll) ~= EXPECTED_FULL_N || ...
   sum(groupsAll=="HC") ~= EXPECTED_FULL_HC || ...
   sum(groupsAll=="PD") ~= EXPECTED_FULL_PD

    error('Full connectivity cohort mismatch.');
end

if sum(completeMask) ~= EXPECTED_COMPLETE_N || ...
   sum(groupsAll(completeMask)=="HC") ~= EXPECTED_COMPLETE_HC || ...
   sum(groupsAll(completeMask)=="PD") ~= EXPECTED_COMPLETE_PD

    error('Complete-network cohort mismatch.');
end

if size(deltaBio,2) ~= 55 || size(highDelta,2) ~= 55
    error('Expected 55 regional edges.');
end

thetaIdx = find(bandNames=="theta");

if numel(thetaIdx) ~= 1
    error('Could not identify theta band.');
end

%% =========================================================================
% 4. IDENTIFY THE PRIMARY SIGNIFICANT THETA NBS COMPONENT
%% =========================================================================

if ~isfield(R,'ComponentTable') || ~isfield(R,'observedComponents')
    error('Primary NBS MAT does not contain expected result fields.');
end

CT = R.ComponentTable;

sig = CT.Significant_FWER_0p05;

thetaNeg = ...
    lower(string(CT.Band))=="theta" & ...
    contains(lower(string(CT.Direction)),"negative") & ...
    sig;

if sum(thetaNeg) ~= 1
    error('Expected exactly one significant theta-negative NBS component.');
end

sigRow = CT(thetaNeg,:);

compID = double(sigRow.ComponentID);

thetaNegComps = R.observedComponents{1,2};

if compID < 1 || compID > numel(thetaNegComps)
    error('Invalid theta component ID.');
end

primaryEdgeIdx = double(thetaNegComps(compID).edgeIdx(:));

if numel(primaryEdgeIdx) ~= 32
    error('Expected primary theta component extent of 32 edges.');
end

fprintf('\n============================================================\n');
fprintf('PRIMARY COMPONENT LOCK\n');
fprintf('============================================================\n');
fprintf('Theta-negative component edges = %d\n',numel(primaryEdgeIdx));
fprintf('Primary NBS joint FWER p       = %.8f\n',double(sigRow.JointFWER_P));

%% =========================================================================
% 5. LOAD SUBJECT QC SUMMARY AND ALIGN
%% =========================================================================

Q = readtable( ...
    subjectSummaryFile, ...
    'TextType','string');

requiredQ = [ ...
    "Subject"
    "Group"
    "NBadChannels"
    "PctRejectedIC"
    "HighPowerRatio"
    "ManualICAFlag"];

for k = 1:numel(requiredQ)
    if ~ismember(requiredQ(k),string(Q.Properties.VariableNames))
        error('Subject summary lacks column: %s',requiredQ(k));
    end
end

[tfQ,locQ] = ismember(subjectsAll,string(Q.Subject));

if ~all(tfQ)
    error('Could not align all connectivity subjects to QC summary.');
end

if ~isequal(string(Q.Group(locQ)),groupsAll)
    error('Group mismatch between connectivity and QC summary.');
end

NBadChannelsAll = column_to_double(Q.NBadChannels(locQ));
PctRejectedICAll = column_to_double(Q.PctRejectedIC(locQ));
HighPowerRatioAll = column_to_double(Q.HighPowerRatio(locQ));
ManualICAFlagAll = column_to_logical(Q.ManualICAFlag(locQ));

if any(~isfinite(NBadChannelsAll)) || ...
   any(~isfinite(PctRejectedICAll)) || ...
   any(~isfinite(HighPowerRatioAll)) || ...
   any(HighPowerRatioAll<=0)

    error('QC variables contain invalid values.');
end

%% =========================================================================
% 6. LOAD / ALIGN DEMOGRAPHICS
%% =========================================================================

P = readtable( ...
    participantsFile, ...
    'FileType','text', ...
    'Delimiter','\t', ...
    'TextType','string');

neededCols = ["participant_id","group","age","sex"];

for k = 1:numel(neededCols)
    if ~ismember(neededCols(k),string(P.Properties.VariableNames))
        error('participants.tsv lacks required column: %s',neededCols(k));
    end
end

[tfP,locP] = ismember(subjectsAll,string(P.participant_id));

if ~all(tfP)
    error('Participant alignment failed.');
end

if ~isequal(string(P.group(locP)),groupsAll)
    error('Group mismatch between connectivity and participants.tsv.');
end

ageAll = column_to_double(P.age(locP));

sexRawAll = upper(strtrim(string(P.sex(locP))));

if any(~isfinite(ageAll)) || any(~ismember(sexRawAll,["F","M"]))
    error('Invalid age/sex data.');
end

GroupAll = double(groupsAll=="PD");
SexAll   = double(sexRawAll=="M");

%% =========================================================================
% 7. PRIMARY COMPLETE-NETWORK SUBJECT-LEVEL SIGNAL
%% =========================================================================

subjects = subjectsAll(completeMask);
groups   = groupsAll(completeMask);

Group = GroupAll(completeMask);
Sex   = SexAll(completeMask);
age   = ageAll(completeMask);

NBadChannels = NBadChannelsAll(completeMask);
PctRejectedIC = PctRejectedICAll(completeMask);
HighPowerRatio = HighPowerRatioAll(completeMask);
ManualICAFlag = ManualICAFlagAll(completeMask);

thetaComplete = deltaBio(completeMask,:,thetaIdx);
highComplete  = highDelta(completeMask,:);

if any(~isfinite(thetaComplete),'all') || ...
   any(~isfinite(highComplete),'all')
    error('Complete cohort contains nonfinite theta/high connectivity.');
end

primaryThetaSignal = mean(thetaComplete(:,primaryEdgeIdx),2);
primaryHighSignal  = mean(highComplete(:,primaryEdgeIdx),2);

thetaWholeSignal = mean(thetaComplete,2);
highWholeSignal  = mean(highComplete,2);

Age_c = age-mean(age);

Xbase = [ones(numel(age),1),Group,Age_c,Sex];

if rank(Xbase) ~= 4
    error('Primary sensitivity design is rank deficient.');
end

%% =========================================================================
% 8. QC BURDEN DESCRIPTIVES
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('A. QC BURDEN DESCRIPTIVES — PRIMARY N=117\n');
fprintf('============================================================\n');

fprintf('NBadChannels: HC median %.1f | PD median %.1f | max %d\n', ...
    median(NBadChannels(Group==0)), ...
    median(NBadChannels(Group==1)), ...
    max(NBadChannels));

fprintf('PctRejectedIC: HC median %.2f%% | PD median %.2f%% | max %.2f%%\n', ...
    median(PctRejectedIC(Group==0)), ...
    median(PctRejectedIC(Group==1)), ...
    max(PctRejectedIC));

fprintf('HighPowerRatio: HC median %.4f | PD median %.4f | range %.4f to %.4f\n', ...
    median(HighPowerRatio(Group==0)), ...
    median(HighPowerRatio(Group==1)), ...
    min(HighPowerRatio),max(HighPowerRatio));

fprintf('Subjects with NBadChannels >= %d: %d\n', ...
    BAD_CHANNEL_EXCLUDE_THRESHOLD, ...
    sum(NBadChannels>=BAD_CHANNEL_EXCLUDE_THRESHOLD));

fprintf('Subjects with ManualICAFlag = true: %d\n', ...
    sum(ManualICAFlag));

fprintf('Subjects with PctRejectedIC > %.1f%%: %d\n', ...
    HIGH_IC_REJECTION_THRESHOLD_PCT, ...
    sum(PctRejectedIC>HIGH_IC_REJECTION_THRESHOLD_PCT));

%% =========================================================================
% 9. BASELINE PRIMARY COMPONENT EFFECT
%% =========================================================================

baseStats = univariate_group_effect_hc3(primaryThetaSignal,Xbase);

fprintf('\nBaseline 32-edge theta component effect:\n');
fprintf('Adjusted PD-HC beta = %.6f\n',baseStats.betaGroup);
fprintf('95%% CI              = [%.6f, %.6f]\n',baseStats.ciLow,baseStats.ciHigh);
fprintf('HC3 p               = %.8g\n',baseStats.p);

%% =========================================================================
% 10. QC COVARIATE ADJUSTMENT / ASSOCIATION
%
% Each QC metric is added separately to the frozen Group+Age+Sex model.
% Report:
%   - QC coefficient
%   - QC p-value
%   - adjusted Group coefficient
%   - percent change in Group beta magnitude vs baseline
%% =========================================================================

QCName = [ ...
    "NBadChannels"
    "PctRejectedIC"
    "HighPowerRatio"];

QCData = [ ...
    NBadChannels, ...
    PctRejectedIC, ...
    HighPowerRatio];

nQC = size(QCData,2);

QC_Beta = nan(nQC,1);
QC_SE = nan(nQC,1);
QC_CI_Low = nan(nQC,1);
QC_CI_High = nan(nQC,1);
QC_P = nan(nQC,1);

GroupBeta_Adjusted = nan(nQC,1);
GroupCI_Low = nan(nQC,1);
GroupCI_High = nan(nQC,1);
GroupP_Adjusted = nan(nQC,1);
PctChange_GroupBetaMagnitude = nan(nQC,1);

for q = 1:nQC

    z = QCData(:,q);

    % Scale for numerical interpretability, but preserve direction.
    z = (z-mean(z))/std(z);

    X = [Xbase,z];

    st = univariate_named_effect_hc3(primaryThetaSignal,X,5);

    QC_Beta(q) = st.beta;
    QC_SE(q) = st.se;
    QC_CI_Low(q) = st.ciLow;
    QC_CI_High(q) = st.ciHigh;
    QC_P(q) = st.p;

    gst = univariate_named_effect_hc3(primaryThetaSignal,X,2);

    GroupBeta_Adjusted(q) = gst.beta;
    GroupCI_Low(q) = gst.ciLow;
    GroupCI_High(q) = gst.ciHigh;
    GroupP_Adjusted(q) = gst.p;

    PctChange_GroupBetaMagnitude(q) = ...
        100*(abs(gst.beta)-abs(baseStats.betaGroup))/abs(baseStats.betaGroup);

end

QCAssociationTable = table( ...
    QCName, ...
    QC_Beta, ...
    QC_SE, ...
    QC_CI_Low, ...
    QC_CI_High, ...
    QC_P, ...
    GroupBeta_Adjusted, ...
    GroupCI_Low, ...
    GroupCI_High, ...
    GroupP_Adjusted, ...
    PctChange_GroupBetaMagnitude);

fprintf('\n============================================================\n');
fprintf('B. PRIMARY THETA SIGNAL vs QC METRICS\n');
fprintf('============================================================\n');
fprintf('QC coefficients are per +1 SD of the QC metric.\n');
disp(QCAssociationTable);

% Joint adjustment for all three QC metrics.
Zall = zscore_columns(QCData);

XallQC = [Xbase,Zall];

jointGroup = univariate_named_effect_hc3(primaryThetaSignal,XallQC,2);

fprintf('\nJoint adjustment for all 3 QC metrics:\n');
fprintf('Group beta = %.6f\n',jointGroup.beta);
fprintf('95%% CI     = [%.6f, %.6f]\n',jointGroup.ciLow,jointGroup.ciHigh);
fprintf('p          = %.8g\n',jointGroup.p);
fprintf('Change in |Group beta| vs baseline = %.2f%%\n', ...
    100*(abs(jointGroup.beta)-abs(baseStats.betaGroup))/abs(baseStats.betaGroup));

%% =========================================================================
% 11. EXCLUSION SENSITIVITY
%% =========================================================================

Sensitivity = [ ...
    "Baseline_complete117"
    "Exclude_badChannels_ge5"
    "Exclude_rejectedIC_gt20pct"
    "Exclude_both_QC_burdens"
    "Exclude_all_manualICA_subjects"];

keepMasks = cell(5,1);

keepMasks{1} = true(numel(Group),1);

keepMasks{2} = ...
    NBadChannels < BAD_CHANNEL_EXCLUDE_THRESHOLD;

keepMasks{3} = ...
    PctRejectedIC <= HIGH_IC_REJECTION_THRESHOLD_PCT;

keepMasks{4} = ...
    keepMasks{2} & keepMasks{3};

% Supplementary robustness:
% exclude every subject that triggered manual ICA review, regardless of
% whether borderline ICs were subsequently restored.
keepMasks{5} = ...
    ~ManualICAFlag;

Sens_N = nan(5,1);
Sens_HC = nan(5,1);
Sens_PD = nan(5,1);

Sens_Beta = nan(5,1);
Sens_SE = nan(5,1);
Sens_CI_Low = nan(5,1);
Sens_CI_High = nan(5,1);
Sens_P = nan(5,1);
Sens_g = nan(5,1);

for s = 1:5

    keep = keepMasks{s};

    y = primaryThetaSignal(keep);
    g = Group(keep);
    a = age(keep);
    sx = Sex(keep);

    X = [ones(sum(keep),1),g,a-mean(a),sx];

    st = univariate_group_effect_hc3(y,X);

    Sens_N(s) = sum(keep);
    Sens_HC(s) = sum(g==0);
    Sens_PD(s) = sum(g==1);

    Sens_Beta(s) = st.betaGroup;
    Sens_SE(s) = st.se;
    Sens_CI_Low(s) = st.ciLow;
    Sens_CI_High(s) = st.ciHigh;
    Sens_P(s) = st.p;
    Sens_g(s) = st.gAdjusted;

end

SensitivityTable = table( ...
    Sensitivity, ...
    Sens_N, ...
    Sens_HC, ...
    Sens_PD, ...
    Sens_Beta, ...
    Sens_SE, ...
    Sens_CI_Low, ...
    Sens_CI_High, ...
    Sens_P, ...
    Sens_g, ...
    'VariableNames', ...
    {'Sensitivity', ...
     'N', ...
     'N_HC', ...
     'N_PD', ...
     'AdjustedBeta_PDminusHC', ...
     'HC3_SE', ...
     'CI95_Low', ...
     'CI95_High', ...
     'P_Raw', ...
     'HedgesType_gAdjusted'});

fprintf('\n============================================================\n');
fprintf('C. PRIMARY THETA QC-EXCLUSION SENSITIVITY\n');
fprintf('============================================================\n');
disp(SensitivityTable);

%% =========================================================================
% 12. HIGH-BAND CONTROL — COMPONENT AND WHOLE NETWORK
%% =========================================================================

highCompStats = univariate_group_effect_hc3(primaryHighSignal,Xbase);
highWholeStats = univariate_group_effect_hc3(highWholeSignal,Xbase);

HighSummary = table( ...
    ["Theta_primary_component";"High_same_32_edges";"Theta_whole55";"High_whole55"], ...
    [mean(primaryThetaSignal(Group==0)); ...
     mean(primaryHighSignal(Group==0)); ...
     mean(thetaWholeSignal(Group==0)); ...
     mean(highWholeSignal(Group==0))], ...
    [mean(primaryThetaSignal(Group==1)); ...
     mean(primaryHighSignal(Group==1)); ...
     mean(thetaWholeSignal(Group==1)); ...
     mean(highWholeSignal(Group==1))], ...
    [baseStats.betaGroup; ...
     highCompStats.betaGroup; ...
     univariate_group_effect_hc3(thetaWholeSignal,Xbase).betaGroup; ...
     highWholeStats.betaGroup], ...
    [baseStats.ciLow; ...
     highCompStats.ciLow; ...
     univariate_group_effect_hc3(thetaWholeSignal,Xbase).ciLow; ...
     highWholeStats.ciLow], ...
    [baseStats.ciHigh; ...
     highCompStats.ciHigh; ...
     univariate_group_effect_hc3(thetaWholeSignal,Xbase).ciHigh; ...
     highWholeStats.ciHigh], ...
    [baseStats.p; ...
     highCompStats.p; ...
     univariate_group_effect_hc3(thetaWholeSignal,Xbase).p; ...
     highWholeStats.p], ...
    'VariableNames', ...
    {'Signal', ...
     'MeanDelta_HC', ...
     'MeanDelta_PD', ...
     'AdjustedBeta_PDminusHC', ...
     'CI95_Low', ...
     'CI95_High', ...
     'P_Raw'});

fprintf('\n============================================================\n');
fprintf('D. 30-45 Hz ARTIFACT-CONTROL SUMMARY\n');
fprintf('============================================================\n');
disp(HighSummary);

%% =========================================================================
% 13. HIGH-BAND EDGE-LEVEL AVAILABLE-CASE CONTROL
%% =========================================================================

HighEdge = edge_level_models( ...
    highDelta, ...
    ageAll, ...
    SexAll, ...
    GroupAll, ...
    edgeNames);

HighEdge.Q_BH_FDR = bh_adjust(HighEdge.P_Raw);
HighEdge.Significant_BH_FDR_0p05 = HighEdge.Q_BH_FDR < FDR_ALPHA;

ThetaEdgeAvailable = edge_level_models( ...
    deltaBio(:,:,thetaIdx), ...
    ageAll, ...
    SexAll, ...
    GroupAll, ...
    edgeNames);

ThetaEdgeAvailable.Q_BH_FDR = bh_adjust(ThetaEdgeAvailable.P_Raw);
ThetaEdgeAvailable.Significant_BH_FDR_0p05 = ...
    ThetaEdgeAvailable.Q_BH_FDR < FDR_ALPHA;

nHighSig = sum(HighEdge.Significant_BH_FDR_0p05);
nThetaSig = sum(ThetaEdgeAvailable.Significant_BH_FDR_0p05);

overlapSig = sum( ...
    HighEdge.Significant_BH_FDR_0p05 & ...
    ThetaEdgeAvailable.Significant_BH_FDR_0p05);

pearsonThetaHigh = corr_manual( ...
    ThetaEdgeAvailable.AdjustedBeta_PDminusHC, ...
    HighEdge.AdjustedBeta_PDminusHC);

spearmanThetaHigh = corr_manual( ...
    tied_rank(ThetaEdgeAvailable.AdjustedBeta_PDminusHC), ...
    tied_rank(HighEdge.AdjustedBeta_PDminusHC));

sameSignThetaHigh = mean( ...
    sign(ThetaEdgeAvailable.AdjustedBeta_PDminusHC) == ...
    sign(HighEdge.AdjustedBeta_PDminusHC));

fprintf('\nHigh-band edge-level control:\n');
fprintf('Theta FDR-significant edges = %d / 55\n',nThetaSig);
fprintf('High  FDR-significant edges = %d / 55\n',nHighSig);
fprintf('FDR-significant overlap     = %d edges\n',overlapSig);
fprintf('Theta-vs-high edge beta Pearson r  = %.4f\n',pearsonThetaHigh);
fprintf('Theta-vs-high edge beta Spearman r = %.4f\n',spearmanThetaHigh);
fprintf('Theta-vs-high beta sign agreement  = %.1f%%\n',100*sameSignThetaHigh);

%% =========================================================================
% 14. HIGH-BAND NBS ARTIFACT CONTROL
%
% Same primary complete cohort, same model, same |v| threshold, same number
% of resamples and same seed. Family here is high+ and high- ONLY because
% this is a separate artifact-control sensitivity test.
%% =========================================================================

Yhigh = highComplete;

vHighObs = compute_signed_v_all( ...
    Yhigh, ...
    Xbase, ...
    2, ...
    Group==0, ...
    Group==1);

[obsHighPosExtent,obsHighPosEdges] = max_component_extent( ...
    vHighObs(:)>NBS_THRESHOLD, ...
    edgeIndex1,edgeIndex2,11);

[obsHighNegExtent,obsHighNegEdges] = max_component_extent( ...
    vHighObs(:)<-NBS_THRESHOLD, ...
    edgeIndex1,edgeIndex2,11);

% Freedman-Lane reduced model: Intercept + Age + Sex
Z = [ones(numel(Group),1),Age_c,Sex];

Hred = Z*pinv(Z);
FittedReduced = Hred*Yhigh;
ResidualReduced = Yhigh-FittedReduced;

rng(RNG_SEED,'twister');

signMatrix = int8(ones(N_RESAMPLES,numel(Group)));

flipMask = rand(N_RESAMPLES,numel(Group))<0.5;

signMatrix(flipMask) = -1;

nullHighMax = zeros(N_RESAMPLES,1);

fprintf('\nRunning high-band artifact-control NBS (%d resamples)...\n',N_RESAMPLES);

for b = 1:N_RESAMPLES

    signs = double(signMatrix(b,:))';

    Ystar = FittedReduced + ResidualReduced.*signs;

    vStar = compute_signed_v_all( ...
        Ystar, ...
        Xbase, ...
        2, ...
        Group==0, ...
        Group==1);

    posExt = max_component_extent_only( ...
        vStar(:)>NBS_THRESHOLD, ...
        edgeIndex1,edgeIndex2,11);

    negExt = max_component_extent_only( ...
        vStar(:)<-NBS_THRESHOLD, ...
        edgeIndex1,edgeIndex2,11);

    nullHighMax(b) = max(posExt,negExt);

end

pHighPos = NaN;
pHighNeg = NaN;

if obsHighPosExtent > 0
    pHighPos = ...
        (1+sum(nullHighMax>=obsHighPosExtent))/(N_RESAMPLES+1);
end

if obsHighNegExtent > 0
    pHighNeg = ...
        (1+sum(nullHighMax>=obsHighNegExtent))/(N_RESAMPLES+1);
end

HighNBS = table( ...
    ["Positive_PD_gt_HC_change";"Negative_PD_lt_HC_change"], ...
    [obsHighPosExtent;obsHighNegExtent], ...
    [pHighPos;pHighNeg], ...
    'VariableNames', ...
    {'Direction','ExtentEdges','FWER_P_HighBandPlusMinus'});

fprintf('\nHigh-band NBS artifact-control result:\n');
disp(HighNBS);

%% =========================================================================
% 15. COMPLETE-NETWORK vs AVAILABLE-CASE THETA EDGE EFFECTS
%% =========================================================================

ThetaEdgeComplete = edge_level_models_complete( ...
    thetaComplete, ...
    age, ...
    Sex, ...
    Group, ...
    edgeNames);

rCompleteAvailablePearson = corr_manual( ...
    ThetaEdgeComplete.AdjustedBeta_PDminusHC, ...
    ThetaEdgeAvailable.AdjustedBeta_PDminusHC);

rCompleteAvailableSpearman = corr_manual( ...
    tied_rank(ThetaEdgeComplete.AdjustedBeta_PDminusHC), ...
    tied_rank(ThetaEdgeAvailable.AdjustedBeta_PDminusHC));

signAgreementCompleteAvailable = mean( ...
    sign(ThetaEdgeComplete.AdjustedBeta_PDminusHC) == ...
    sign(ThetaEdgeAvailable.AdjustedBeta_PDminusHC));

medianAbsBetaDifference = median(abs( ...
    ThetaEdgeComplete.AdjustedBeta_PDminusHC - ...
    ThetaEdgeAvailable.AdjustedBeta_PDminusHC));

fprintf('\n============================================================\n');
fprintf('E. COMPLETE vs AVAILABLE-CASE THETA EDGE ROBUSTNESS\n');
fprintf('============================================================\n');

fprintf('Pearson beta concordance  = %.4f\n',rCompleteAvailablePearson);
fprintf('Spearman beta concordance = %.4f\n',rCompleteAvailableSpearman);
fprintf('Sign agreement            = %.1f%%\n',100*signAgreementCompleteAvailable);
fprintf('Median |beta difference|  = %.6f db-wPLI\n',medianAbsBetaDifference);

%% =========================================================================
% 16. SUMMARY FLAGS
%% =========================================================================

baselineAbs = abs(baseStats.betaGroup);

jointQCAttenuationPct = ...
    100*(1-abs(jointGroup.beta)/baselineAbs);

allSensitivitySameSign = all(sign(Sens_Beta)==sign(baseStats.betaGroup));

minSensitivityRatio = min(abs(Sens_Beta))/baselineAbs;

highNBSAnySig = any(HighNBS.FWER_P_HighBandPlusMinus < 0.05);

SummaryMetric = [ ...
    "Primary_theta_beta"
    "Joint_QC_adjusted_theta_beta"
    "Joint_QC_beta_attenuation_percent"
    "Minimum_sensitivity_beta_ratio"
    "All_sensitivity_effects_same_sign"
    "High_band_FDR_edges"
    "High_band_NBS_any_FWER_lt_0p05"
    "ThetaHigh_edge_beta_Pearson_r"
    "CompleteAvailable_theta_beta_Pearson_r"
    "CompleteAvailable_theta_sign_agreement_percent"];

SummaryValue = [ ...
    baseStats.betaGroup
    jointGroup.beta
    jointQCAttenuationPct
    minSensitivityRatio
    double(allSensitivitySameSign)
    nHighSig
    double(highNBSAnySig)
    pearsonThetaHigh
    rCompleteAvailablePearson
    100*signAgreementCompleteAvailable];

SummaryTable = table(SummaryMetric,SummaryValue);

%% =========================================================================
% 17. SAVE OUTPUTS
%% =========================================================================

qcAssocFile = fullfile( ...
    outRoot, ...
    'Paper3_QC_primary_theta_associations.csv');

sensitivityFile = fullfile( ...
    outRoot, ...
    'Paper3_QC_exclusion_sensitivity.csv');

highSummaryFile = fullfile( ...
    outRoot, ...
    'Paper3_HIGHBAND_component_wholenetwork_summary.csv');

highEdgesFile = fullfile( ...
    outRoot, ...
    'Paper3_HIGHBAND_edge_control_BHFDR.csv');

highNBSFile = fullfile( ...
    outRoot, ...
    'Paper3_HIGHBAND_NBS_artifact_control.csv');

completeVsAvailableFile = fullfile( ...
    outRoot, ...
    'Paper3_THETA_complete_vs_available_edge_effects.csv');

summaryFile = fullfile( ...
    outRoot, ...
    'Paper3_ARTIFACT_SENSITIVITY_summary.csv');

writetable(QCAssociationTable,qcAssocFile);
writetable(SensitivityTable,sensitivityFile);
writetable(HighSummary,highSummaryFile);
writetable(HighEdge,highEdgesFile);
writetable(HighNBS,highNBSFile);
writetable(SummaryTable,summaryFile);

CompleteVsAvailable = table( ...
    edgeNames, ...
    ThetaEdgeComplete.AdjustedBeta_PDminusHC, ...
    ThetaEdgeAvailable.AdjustedBeta_PDminusHC, ...
    ThetaEdgeAvailable.N_Total, ...
    'VariableNames', ...
    {'Edge', ...
     'Beta_CompleteN117', ...
     'Beta_AvailableCase', ...
     'N_AvailableCase'});

writetable(CompleteVsAvailable,completeVsAvailableFile);

resultFile = fullfile( ...
    outRoot, ...
    'Paper3_ARTIFACT_SENSITIVITY_RESULTS.mat');

save( ...
    resultFile, ...
    'QCAssociationTable', ...
    'SensitivityTable', ...
    'HighSummary', ...
    'HighEdge', ...
    'ThetaEdgeAvailable', ...
    'ThetaEdgeComplete', ...
    'HighNBS', ...
    'nullHighMax', ...
    'obsHighPosEdges', ...
    'obsHighNegEdges', ...
    'SummaryTable', ...
    'primaryEdgeIdx', ...
    'primaryThetaSignal', ...
    'primaryHighSignal', ...
    'subjects', ...
    'groups', ...
    '-v7.3');

%% =========================================================================
% 18. LOG
%% =========================================================================

logFile = fullfile( ...
    outRoot, ...
    'Paper3_ARTIFACT_SENSITIVITY_analysis_log.txt');

fid = fopen(logFile,'w');

if fid==-1
    error('Could not open analysis log.');
end

fprintf(fid,'PAPER 3 ARTIFACT / SENSITIVITY ANALYSIS\n');
fprintf(fid,'Run time: %s\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));

fprintf(fid,'Primary complete cohort: 117 = 22 HC + 95 PD\n');
fprintf(fid,'Primary theta component: 32 edges\n');

fprintf(fid,'Bad-channel sensitivity: exclude NBadChannels >= %d\n', ...
    BAD_CHANNEL_EXCLUDE_THRESHOLD);
fprintf(fid,'Manual-ICA sensitivity: exclude every ManualICAFlag subject\n');

fprintf(fid,'ICA sensitivity: exclude PctRejectedIC > %.1f%%\n', ...
    HIGH_IC_REJECTION_THRESHOLD_PCT);

fprintf(fid,'High band: 30-45 Hz artifact control\n');
fprintf(fid,'High-band NBS: |v|>%.2f, B=%d, seed=%d, high+/- max family\n', ...
    NBS_THRESHOLD,N_RESAMPLES,RNG_SEED);

fprintf(fid,'High-band FDR edges: %d/55\n',nHighSig);
fprintf(fid,'High-band NBS any p<0.05: %d\n',highNBSAnySig);

fprintf(fid,'Theta/high edge beta Pearson r: %.6f\n',pearsonThetaHigh);
fprintf(fid,'Complete/available theta beta Pearson r: %.6f\n', ...
    rCompleteAvailablePearson);

fclose(fid);

%% =========================================================================
% 19. FINAL CONSOLE SUMMARY
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('ARTIFACT / SENSITIVITY ANALYSIS COMPLETE\n');
fprintf('============================================================\n');

disp(SummaryTable);

fprintf('\nInterpretation guardrails:\n');

if highNBSAnySig || nHighSig>0
    fprintf([ ...
        'ATTENTION: the 30-45 Hz artifact-control band shows corrected\n' ...
        'structure. Inspect topology and QC associations carefully before\n' ...
        'making a biological theta interpretation.\n']);
else
    fprintf([ ...
        'The 30-45 Hz artifact-control band does not show a corrected\n' ...
        'edge/NBS effect under these prespecified checks.\n']);
end

if allSensitivitySameSign && minSensitivityRatio>=0.70
    fprintf([ ...
        'The primary theta Group effect retains direction and at least 70%%\n' ...
        'of its original magnitude across the predefined QC exclusions.\n']);
else
    fprintf([ ...
        'The primary theta effect materially changes under at least one QC\n' ...
        'restriction. Inspect sensitivity rows before interpretation.\n']);
end

fprintf('\nSaved files:\n');
fprintf('QC associations       : %s\n',qcAssocFile);
fprintf('QC exclusions         : %s\n',sensitivityFile);
fprintf('High summary          : %s\n',highSummaryFile);
fprintf('High edge control     : %s\n',highEdgesFile);
fprintf('High NBS control      : %s\n',highNBSFile);
fprintf('Complete vs available : %s\n',completeVsAvailableFile);
fprintf('Summary               : %s\n',summaryFile);
fprintf('Full MAT              : %s\n',resultFile);
fprintf('Analysis log          : %s\n',logFile);

fprintf('\n============================================================\n');
fprintf('DO NOT CHANGE QC THRESHOLDS AFTER VIEWING THESE RESULTS.\n');
fprintf('============================================================\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function st = univariate_group_effect_hc3(y,X)

    st0 = univariate_named_effect_hc3(y,X,2);

    residual = double(y(:))-X*(X\double(y(:)));

    df = numel(y)-size(X,2);

    rmse = sqrt(sum(residual.^2)/df);

    dAdj = st0.beta/rmse;

    J = 1-3/(4*df-1);

    st = st0;

    st.betaGroup = st0.beta;
    st.gAdjusted = J*dAdj;

end

%% =========================================================================

function st = univariate_named_effect_hc3(y,X,coefIdx)

    y = double(y(:));

    N = numel(y);
    p = size(X,2);

    beta = X\y;

    residual = y-X*beta;

    XtXinv = pinv(X'*X);

    Hdiag = sum((X*XtXinv).*X,2);

    if any(Hdiag>=1)
        error('Invalid leverage in HC3.');
    end

    uAdj2 = residual.^2./(1-Hdiag).^2;

    meat = X'*(X.*uAdj2);

    V = XtXinv*meat*XtXinv;

    se = sqrt(V(coefIdx,coefIdx));

    t = beta(coefIdx)/se;

    df = N-p;

    pValue = t_two_sided_p(t,df);

    tc = student_t_critical_975(df);

    st = struct();

    st.beta = beta(coefIdx);
    st.se = se;
    st.ciLow = beta(coefIdx)-tc*se;
    st.ciHigh = beta(coefIdx)+tc*se;
    st.t = t;
    st.df = df;
    st.p = pValue;

end

%% =========================================================================

function T = edge_level_models(Y,ageAll,SexAll,GroupAll,edgeNames)

    nEdges = size(Y,2);

    N_Total = nan(nEdges,1);
    N_HC = nan(nEdges,1);
    N_PD = nan(nEdges,1);

    AdjustedBeta_PDminusHC = nan(nEdges,1);
    CI95_Low = nan(nEdges,1);
    CI95_High = nan(nEdges,1);
    P_Raw = nan(nEdges,1);

    for e = 1:nEdges

        y0 = Y(:,e);

        use = ...
            isfinite(y0) & ...
            isfinite(ageAll) & ...
            isfinite(SexAll) & ...
            isfinite(GroupAll);

        y = y0(use);
        g = GroupAll(use);
        a = ageAll(use);
        sx = SexAll(use);

        X = [ones(sum(use),1),g,a-mean(a),sx];

        st = univariate_group_effect_hc3(y,X);

        N_Total(e) = sum(use);
        N_HC(e) = sum(g==0);
        N_PD(e) = sum(g==1);

        AdjustedBeta_PDminusHC(e) = st.betaGroup;
        CI95_Low(e) = st.ciLow;
        CI95_High(e) = st.ciHigh;
        P_Raw(e) = st.p;

    end

    T = table( ...
        edgeNames(:), ...
        N_Total, ...
        N_HC, ...
        N_PD, ...
        AdjustedBeta_PDminusHC, ...
        CI95_Low, ...
        CI95_High, ...
        P_Raw, ...
        'VariableNames', ...
        {'Edge','N_Total','N_HC','N_PD', ...
         'AdjustedBeta_PDminusHC','CI95_Low','CI95_High','P_Raw'});

end

%% =========================================================================

function T = edge_level_models_complete(Y,age,Sex,Group,edgeNames)

    nEdges = size(Y,2);

    beta = nan(nEdges,1);
    lo = nan(nEdges,1);
    hi = nan(nEdges,1);
    p = nan(nEdges,1);

    X = [ones(numel(age),1),Group,age-mean(age),Sex];

    for e = 1:nEdges

        st = univariate_group_effect_hc3(Y(:,e),X);

        beta(e) = st.betaGroup;
        lo(e) = st.ciLow;
        hi(e) = st.ciHigh;
        p(e) = st.p;

    end

    T = table( ...
        edgeNames(:), ...
        beta, ...
        lo, ...
        hi, ...
        p, ...
        'VariableNames', ...
        {'Edge','AdjustedBeta_PDminusHC','CI95_Low','CI95_High','P_Raw'});

end

%% =========================================================================
% Signed Aspin-Welch-type v used in frozen primary NBS implementation.
%% =========================================================================

function v = compute_signed_v_all(Y,M,groupColumn,vgHC,vgPD)

    Mpinv = pinv(M);

    Bhat = Mpinv*Y;

    E = Y-M*Bhat;

    H = M*Mpinv;
    Rdiag = diag(eye(size(M,1))-H);

    sumR_HC = sum(Rdiag(vgHC));
    sumR_PD = sum(Rdiag(vgPD));

    sseHC = sum(E(vgHC,:).^2,1);
    ssePD = sum(E(vgPD,:).^2,1);

    wHC = sumR_HC./sseHC;
    wPD = sumR_PD./ssePD;

    A_HC = M(vgHC,:)'*M(vgHC,:);
    A_PD = M(vgPD,:)'*M(vgPD,:);

    C = zeros(size(M,2),1);
    C(groupColumn) = 1;

    betaGroup = Bhat(groupColumn,:);

    K = size(Y,2);

    v = nan(1,K);

    for k = 1:K

        MtWM = wHC(k)*A_HC + wPD(k)*A_PD;

        q = C'*(MtWM\C);

        if ~isfinite(q) || q<=0
            error('Invalid contrast variance in v statistic.');
        end

        v(k) = betaGroup(k)/sqrt(q);

    end

end

%% =========================================================================

function [maxExt,edgeIdxBest] = max_component_extent( ...
    edgeMask,edgeIndex1,edgeIndex2,nNodes)

    selectedEdges = find(edgeMask);

    if isempty(selectedEdges)
        maxExt = 0;
        edgeIdxBest = zeros(0,1);
        return;
    end

    A = false(nNodes,nNodes);

    for q = 1:numel(selectedEdges)

        e = selectedEdges(q);

        i = edgeIndex1(e);
        j = edgeIndex2(e);

        A(i,j) = true;
        A(j,i) = true;

    end

    activeNodes = find(any(A,2));
    visited = false(nNodes,1);

    maxExt = 0;
    edgeIdxBest = zeros(0,1);

    for seed = activeNodes(:)'

        if visited(seed)
            continue;
        end

        queue = seed;
        visited(seed) = true;

        nodesThis = zeros(0,1);

        while ~isempty(queue)

            node = queue(1);
            queue(1) = [];

            nodesThis(end+1,1) = node; %#ok<AGROW>

            nbr = find(A(node,:));

            for nb = nbr
                if ~visited(nb)
                    visited(nb) = true;
                    queue(end+1) = nb; %#ok<AGROW>
                end
            end

        end

        nodeMask = false(nNodes,1);
        nodeMask(nodesThis) = true;

        edgeThis = ...
            edgeMask & ...
            nodeMask(edgeIndex1) & ...
            nodeMask(edgeIndex2);

        ext = sum(edgeThis);

        if ext > maxExt
            maxExt = ext;
            edgeIdxBest = find(edgeThis);
        end

    end

end

%% =========================================================================

function maxExt = max_component_extent_only( ...
    edgeMask,edgeIndex1,edgeIndex2,nNodes)

    [maxExt,~] = max_component_extent( ...
        edgeMask,edgeIndex1,edgeIndex2,nNodes);

end

%% =========================================================================

function q = bh_adjust(p)

    p = double(p(:));

    m = numel(p);

    [pSort,ord] = sort(p,'ascend');

    raw = pSort.*m./(1:m)';

    qSort = nan(m,1);

    runningMin = 1;

    for i = m:-1:1
        runningMin = min(runningMin,raw(i));
        qSort(i) = min(1,runningMin);
    end

    q = nan(m,1);
    q(ord) = qSort;

end

%% =========================================================================

function Z = zscore_columns(X)

    Z = nan(size(X));

    for j = 1:size(X,2)

        mu = mean(X(:,j));
        sd = std(X(:,j));

        if sd<=0
            error('QC covariate has zero variance.');
        end

        Z(:,j) = (X(:,j)-mu)/sd;

    end

end

%% =========================================================================

function r = corr_manual(x,y)

    x = double(x(:));
    y = double(y(:));

    use = isfinite(x)&isfinite(y);

    x = x(use);
    y = y(use);

    x = x-mean(x);
    y = y-mean(y);

    r = sum(x.*y)/sqrt(sum(x.^2)*sum(y.^2));

end

%% =========================================================================

function r = tied_rank(x)

    x = double(x(:));

    [xs,ord] = sort(x);

    rSort = nan(size(x));

    i = 1;

    while i<=numel(x)

        j = i;

        while j<numel(x) && xs(j+1)==xs(i)
            j = j+1;
        end

        rSort(i:j) = mean(i:j);

        i = j+1;

    end

    r = nan(size(x));
    r(ord) = rSort;

end

%% =========================================================================

function p = t_two_sided_p(t,df)

    x = df/(df+t^2);

    p = betainc(x,df/2,0.5);

    p = max(0,min(1,p));

end

%% =========================================================================

function tcrit = student_t_critical_975(df)

    x = betaincinv(0.05,df/2,0.5);

    tcrit = sqrt(df*(1-x)/x);

end

%% =========================================================================

function x = column_to_double(col)

    if isnumeric(col)
        x = double(col(:));
    else
        x = str2double(string(col(:)));
    end

end


%% =========================================================================

function x = column_to_logical(col)

    if islogical(col)
        x = col(:);
        return;
    end

    if isnumeric(col)
        v = double(col(:));
        if any(~ismember(v,[0 1]))
            error('Numeric logical column contains values other than 0/1.');
        end
        x = logical(v);
        return;
    end

    s = lower(strtrim(string(col(:))));

    trueMask = ismember(s,["true","1","yes","y"]);
    falseMask = ismember(s,["false","0","no","n"]);

    if any(~(trueMask|falseMask))
        bad = find(~(trueMask|falseMask),1);
        error('Cannot parse logical value: %s',s(bad));
    end

    x = false(size(s));
    x(trueMask) = true;
    x(falseMask) = false;

end

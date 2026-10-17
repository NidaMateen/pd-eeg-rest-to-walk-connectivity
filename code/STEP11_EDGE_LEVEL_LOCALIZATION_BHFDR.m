%% =========================================================================
% PAPER 3 — STEP 11
% PRESPECIFIED EDGE-LEVEL LOCALIZATION
%
% File:
%   STEP11_EDGE_LEVEL_LOCALIZATION_BHFDR.m
%
% PURPOSE
% -------
% Localize Rest->Walk connectivity reconfiguration differences between
% Parkinson's disease (PD) and healthy controls (HC) at the 55 regional
% edges, separately for:
%
%   theta 4-8 Hz
%   alpha 8-13 Hz
%   beta  13-30 Hz
%
% This is the prespecified Level-3 localization analysis. It does NOT
% replace or redefine the primary NBS result.
%
% OUTCOME
% -------
%   Delta db-wPLI = Walk - Rest
%
% EDGE MODEL
% ----------
% For each edge and each band:
%
%   Delta db-wPLI =
%       beta0 + beta1*Group + beta2*Age_c + beta3*Sex + error
%
% Coding:
%   Group: HC = 0, PD = 1
%   Sex:   F  = 0, M  = 1
%
% Therefore:
%   beta1 > 0 : larger Walk-Rest change in PD than HC
%   beta1 < 0 : smaller Walk-Rest change in PD than HC
%
% COHORT
% ------
% AVAILABLE CASE PER EDGE.
%
% A participant contributes to an edge only if that edge is finite for that
% band. This intentionally preserves the frozen regional-coverage plan:
%
%   - most edges: N=130 = 23 HC + 107 PD
%   - Occipital edges: typically N=126 = 23 HC + 103 PD
%   - ParietalMidline edges: typically N=118 = 22 HC + 96 PD
%   - ParietalMidline<->Occipital: N=117 = 22 HC + 95 PD
%
% INFERENCE
% ---------
% OLS coefficients with HC3 heteroscedasticity-robust covariance.
%
% For each edge report:
%   - HC / PD sample sizes
%   - HC / PD mean Delta db-wPLI
%   - adjusted Group beta (PD-HC)
%   - HC3 robust SE
%   - 95% CI
%   - two-sided p-value
%   - BH-FDR q-value
%
% MULTIPLE COMPARISON CONTROL
% ---------------------------
% Benjamini-Hochberg FDR is applied SEPARATELY within each band across the
% 55 regional edges:
%
%   theta: 55 p-values -> 55 q-values
%   alpha: 55 p-values -> 55 q-values
%   beta : 55 p-values -> 55 q-values
%
% FDR threshold:
%   q < 0.05
%
% IMPORTANT
% ---------
% - Do not pool the 165 tests into one BH family.
% - Do not change the q threshold after seeing results.
% - Do not use nominal p<0.05 edges as "significant" if q>=0.05.
% - NBS component membership does not imply individual-edge significance.
%% =========================================================================

clear;
clc;
close all;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% =========================================================================
% 1. FROZEN SETTINGS
%% =========================================================================

BIO_BANDS = [ ...
    "theta"
    "alpha"
    "beta"];

FDR_ALPHA = 0.05;

EXPECTED_TOTAL_SUBJECTS = 130;
EXPECTED_TOTAL_HC = 23;
EXPECTED_TOTAL_PD = 107;

%% =========================================================================
% 2. PATHS
%% =========================================================================

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

outRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_EDGE_LOCALIZATION');

if ~exist(outRoot,'dir')
    mkdir(outRoot);
end

if ~exist(statsInputFile,'file')
    error('Missing STEP 8 statistics input:\n%s',statsInputFile);
end

if ~exist(participantsFile,'file')
    error('Missing participants.tsv:\n%s',participantsFile);
end

%% =========================================================================
% 3. LOAD CONNECTIVITY DATA
%% =========================================================================

S = load(statsInputFile);

requiredFields = { ...
    'subjects', ...
    'groups', ...
    'edgeNames', ...
    'edgeRegion1', ...
    'edgeRegion2', ...
    'bandNames', ...
    'deltaBio'};

for k = 1:numel(requiredFields)

    if ~isfield(S,requiredFields{k})
        error('STEP 8 MAT lacks required field: %s',requiredFields{k});
    end

end

subjectsAll = string(S.subjects(:));
groupsAll   = string(S.groups(:));

edgeNames   = string(S.edgeNames(:));
edgeRegion1 = string(S.edgeRegion1(:));
edgeRegion2 = string(S.edgeRegion2(:));

bandNamesAll = lower(string(S.bandNames(:)));

deltaBio = double(S.deltaBio);

if numel(subjectsAll) ~= EXPECTED_TOTAL_SUBJECTS
    error('Expected 130 connectivity subjects; found %d.',numel(subjectsAll));
end

if sum(groupsAll=="HC") ~= EXPECTED_TOTAL_HC || ...
   sum(groupsAll=="PD") ~= EXPECTED_TOTAL_PD

    error('Expected full connectivity cohort 23 HC + 107 PD.');

end

if numel(edgeNames) ~= 55 || size(deltaBio,2) ~= 55
    error('Expected exactly 55 regional edges.');
end

%% =========================================================================
% 4. IDENTIFY BIOLOGICAL BAND INDICES
%% =========================================================================

bandIdx = nan(3,1);

for b = 1:3

    idx = find(bandNamesAll == BIO_BANDS(b));

    if numel(idx) ~= 1
        error('Could not uniquely identify band: %s',BIO_BANDS(b));
    end

    bandIdx(b) = idx;

end

%% =========================================================================
% 5. LOAD AND ALIGN AGE / SEX FOR ALL 130 SUBJECTS
%% =========================================================================

P = readtable( ...
    participantsFile, ...
    'FileType','text', ...
    'Delimiter','\t', ...
    'TextType','string');

neededCols = [ ...
    "participant_id"
    "group"
    "age"
    "sex"];

for k = 1:numel(neededCols)

    if ~ismember(neededCols(k),string(P.Properties.VariableNames))
        error('participants.tsv lacks required column: %s',neededCols(k));
    end

end

[tf,loc] = ismember(subjectsAll,string(P.participant_id));

if ~all(tf)

    fprintf('Subjects missing from participants.tsv:\n');
    disp(subjectsAll(~tf));

    error('Participant alignment failed.');

end

participantGroup = string(P.group(loc));

if ~isequal(participantGroup,groupsAll)

    bad = participantGroup ~= groupsAll;

    disp(table( ...
        subjectsAll(bad), ...
        groupsAll(bad), ...
        participantGroup(bad), ...
        'VariableNames', ...
        {'Subject','ConnectivityGroup','ParticipantsGroup'}));

    error('Group mismatch between STEP 8 and participants.tsv.');

end

ageAll = column_to_double(P.age(loc));

if any(~isfinite(ageAll))
    error('Missing/nonfinite age in connectivity cohort.');
end

sexRawAll = upper(strtrim(string(P.sex(loc))));

if any(~ismember(sexRawAll,["F","M"]))
    error('Sex must be F or M for all subjects.');
end

GroupAll = double(groupsAll=="PD");
SexAll   = double(sexRawAll=="M");

%% =========================================================================
% 6. PREALLOCATE 165-ROW RESULT TABLE
%% =========================================================================

nEdges = 55;
nBands = 3;
nRows = nEdges*nBands;

Band = strings(nRows,1);
EdgeIndex = nan(nRows,1);

Region1 = strings(nRows,1);
Region2 = strings(nRows,1);
Edge = strings(nRows,1);

N_Total = nan(nRows,1);
N_HC = nan(nRows,1);
N_PD = nan(nRows,1);

MeanDelta_HC = nan(nRows,1);
MeanDelta_PD = nan(nRows,1);

AdjustedBeta_PDminusHC = nan(nRows,1);
HC3_SE = nan(nRows,1);

CI95_Low = nan(nRows,1);
CI95_High = nan(nRows,1);

T_HC3 = nan(nRows,1);
DF = nan(nRows,1);
P_Raw = nan(nRows,1);

Q_BH_FDR = nan(nRows,1);
Significant_BH_FDR_0p05 = false(nRows,1);

HedgesType_gAdjusted = nan(nRows,1);

%% =========================================================================
% 7. FIT AVAILABLE-CASE EDGE MODELS
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('PRESPECIFIED EDGE-LEVEL LOCALIZATION\n');
fprintf('============================================================\n');

row = 0;

for b = 1:nBands

    bb = bandIdx(b);

    fprintf('\n%s band\n',upper(BIO_BANDS(b)));

    for e = 1:nEdges

        row = row+1;

        yAll = deltaBio(:,e,bb);

        use = ...
            isfinite(yAll) & ...
            isfinite(ageAll) & ...
            isfinite(SexAll) & ...
            isfinite(GroupAll);

        y = yAll(use);

        group = GroupAll(use);
        age = ageAll(use);
        sex = SexAll(use);

        groupLabel = groupsAll(use);

        n = numel(y);
        nhc = sum(group==0);
        npd = sum(group==1);

        if nhc < 5 || npd < 5
            error( ...
                'Insufficient group sample size for %s / %s: HC=%d PD=%d.', ...
                BIO_BANDS(b),edgeNames(e),nhc,npd);
        end

        Age_c = age-mean(age);

        X = [ ...
            ones(n,1), ...
            group, ...
            Age_c, ...
            sex];

        if rank(X) ~= size(X,2)
            error( ...
                'Rank-deficient design at %s / %s.', ...
                BIO_BANDS(b),edgeNames(e));
        end

        st = univariate_group_effect_hc3(y,X);

        Band(row) = BIO_BANDS(b);
        EdgeIndex(row) = e;

        Region1(row) = edgeRegion1(e);
        Region2(row) = edgeRegion2(e);
        Edge(row) = edgeNames(e);

        N_Total(row) = n;
        N_HC(row) = nhc;
        N_PD(row) = npd;

        MeanDelta_HC(row) = mean(y(group==0));
        MeanDelta_PD(row) = mean(y(group==1));

        AdjustedBeta_PDminusHC(row) = st.betaGroup;
        HC3_SE(row) = st.se;

        CI95_Low(row) = st.ciLow;
        CI95_High(row) = st.ciHigh;

        T_HC3(row) = st.t;
        DF(row) = st.df;
        P_Raw(row) = st.p;

        HedgesType_gAdjusted(row) = st.gAdjusted;

    end

end

%% =========================================================================
% 8. BH-FDR SEPARATELY WITHIN EACH BAND
%% =========================================================================

for b = 1:nBands

    idx = find(Band==BIO_BANDS(b));

    if numel(idx) ~= 55
        error('Expected 55 %s-band edge tests.',BIO_BANDS(b));
    end

    q = bh_adjust(P_Raw(idx));

    Q_BH_FDR(idx) = q;

    Significant_BH_FDR_0p05(idx) = q < FDR_ALPHA;

end

%% =========================================================================
% 9. BUILD MASTER RESULT TABLE
%% =========================================================================

EdgeResults = table( ...
    Band, ...
    EdgeIndex, ...
    Region1, ...
    Region2, ...
    Edge, ...
    N_Total, ...
    N_HC, ...
    N_PD, ...
    MeanDelta_HC, ...
    MeanDelta_PD, ...
    AdjustedBeta_PDminusHC, ...
    HC3_SE, ...
    CI95_Low, ...
    CI95_High, ...
    T_HC3, ...
    DF, ...
    P_Raw, ...
    Q_BH_FDR, ...
    Significant_BH_FDR_0p05, ...
    HedgesType_gAdjusted);

%% =========================================================================
% 10. QC SAMPLE-SIZE PATTERN
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('AVAILABLE-CASE SAMPLE-SIZE QC\n');
fprintf('============================================================\n');

for b = 1:nBands

    idx = Band==BIO_BANDS(b);

    fprintf('%-5s unique N_Total = %s\n', ...
        BIO_BANDS(b), ...
        mat2str(unique(N_Total(idx))'));

end

% Frozen sample-size counts by edge class should be among:
expectedN = [117 118 126 130];

unexpectedN = setdiff(unique(N_Total),expectedN);

if ~isempty(unexpectedN)

    fprintf('\nUnexpected edge sample sizes:\n');
    disp(unexpectedN);

    error('Available-case edge N differs from frozen regional-coverage audit.');

end

%% =========================================================================
% 11. SAVE ALL 165 EDGE TESTS
%% =========================================================================

allFile = fullfile( ...
    outRoot, ...
    'Paper3_EDGE_LOCALIZATION_all_165_tests.csv');

writetable(EdgeResults,allFile);

%% =========================================================================
% 12. SAVE FDR-SIGNIFICANT EDGES ONLY
%% =========================================================================

SignificantEdges = EdgeResults( ...
    EdgeResults.Significant_BH_FDR_0p05, ...
    :);

if ~isempty(SignificantEdges)

    SignificantEdges = sortrows( ...
        SignificantEdges, ...
        {'Band','Q_BH_FDR','P_Raw'}, ...
        {'ascend','ascend','ascend'});

end

sigFile = fullfile( ...
    outRoot, ...
    'Paper3_EDGE_LOCALIZATION_FDR_significant_edges.csv');

writetable(SignificantEdges,sigFile);

%% =========================================================================
% 13. SAVE NOMINAL P<0.05 FOR TRANSPARENT SUPPLEMENTARY REPORTING
%
% These are explicitly labelled NOMINAL and must not be called significant.
%% =========================================================================

NominalEdges = EdgeResults(P_Raw<0.05,:);

if ~isempty(NominalEdges)

    NominalEdges = sortrows( ...
        NominalEdges, ...
        {'Band','P_Raw'}, ...
        {'ascend','ascend'});

end

nominalFile = fullfile( ...
    outRoot, ...
    'Paper3_EDGE_LOCALIZATION_nominal_p_lt_0p05_NOT_FDR.csv');

writetable(NominalEdges,nominalFile);

%% =========================================================================
% 14. BAND SUMMARY
%% =========================================================================

BandSummary_Band = BIO_BANDS;
N_FDR_Significant = nan(3,1);
N_Nominal_Plt0p05 = nan(3,1);

Min_P = nan(3,1);
Min_Q = nan(3,1);

MedianBeta = nan(3,1);

for b = 1:3

    idx = Band==BIO_BANDS(b);

    N_FDR_Significant(b) = sum(Significant_BH_FDR_0p05(idx));
    N_Nominal_Plt0p05(b) = sum(P_Raw(idx)<0.05);

    Min_P(b) = min(P_Raw(idx));
    Min_Q(b) = min(Q_BH_FDR(idx));

    MedianBeta(b) = median(AdjustedBeta_PDminusHC(idx));

end

BandSummary = table( ...
    BandSummary_Band, ...
    N_FDR_Significant, ...
    N_Nominal_Plt0p05, ...
    Min_P, ...
    Min_Q, ...
    MedianBeta, ...
    'VariableNames', ...
    {'Band', ...
     'N_FDR_Significant', ...
     'N_Nominal_Plt0p05', ...
     'Min_P', ...
     'Min_Q', ...
     'Median_AdjustedBeta_PDminusHC'});

summaryFile = fullfile( ...
    outRoot, ...
    'Paper3_EDGE_LOCALIZATION_band_summary.csv');

writetable(BandSummary,summaryFile);

%% =========================================================================
% 15. SAVE MAT
%% =========================================================================

resultFile = fullfile( ...
    outRoot, ...
    'Paper3_EDGE_LOCALIZATION_RESULTS.mat');

R = struct();

R.subjectsAll = subjectsAll;
R.groupsAll = groupsAll;

R.ageAll = ageAll;
R.sexRawAll = sexRawAll;

R.bandNames = BIO_BANDS;

R.edgeNames = edgeNames;
R.edgeRegion1 = edgeRegion1;
R.edgeRegion2 = edgeRegion2;

R.EdgeResults = EdgeResults;
R.SignificantEdges = SignificantEdges;
R.NominalEdges = NominalEdges;
R.BandSummary = BandSummary;

save(resultFile,'-struct','R','-v7.3');

%% =========================================================================
% 16. ANALYSIS LOG
%% =========================================================================

logFile = fullfile( ...
    outRoot, ...
    'Paper3_EDGE_LOCALIZATION_analysis_log.txt');

fid = fopen(logFile,'w');

if fid == -1
    error('Could not open edge-localization analysis log.');
end

fprintf(fid,'PAPER 3 PRESPECIFIED EDGE-LEVEL LOCALIZATION\n');
fprintf(fid,'Run time: %s\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));

fprintf(fid,'Outcome: Delta db-wPLI = Walk - Rest\n');
fprintf(fid,'Model per edge: Delta ~ Group + Age_c + Sex\n');
fprintf(fid,'Coding: HC=0, PD=1; F=0, M=1\n');
fprintf(fid,'Inference: OLS beta with HC3 robust covariance\n');

fprintf(fid,'Cohort rule: available case per edge\n');
fprintf(fid,'Multiplicity: BH-FDR separately within each band across 55 edges\n');
fprintf(fid,'FDR alpha: %.3f\n',FDR_ALPHA);

fprintf(fid,'Bands: theta, alpha, beta\n');
fprintf(fid,'Total edge tests: 165\n\n');

for b = 1:3

    fprintf(fid,'%s: FDR significant = %d; nominal p<0.05 = %d; min q = %.8g\n', ...
        BIO_BANDS(b), ...
        N_FDR_Significant(b), ...
        N_Nominal_Plt0p05(b), ...
        Min_Q(b));

end

fclose(fid);

%% =========================================================================
% 17. CONSOLE REPORT
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('EDGE-LEVEL LOCALIZATION COMPLETE\n');
fprintf('============================================================\n');

disp(BandSummary);

fprintf('\nFDR-significant regional edges (q < 0.05):\n');

if isempty(SignificantEdges)

    fprintf('None.\n');

else

    disp(SignificantEdges(:,{ ...
        'Band', ...
        'Edge', ...
        'N_HC', ...
        'N_PD', ...
        'MeanDelta_HC', ...
        'MeanDelta_PD', ...
        'AdjustedBeta_PDminusHC', ...
        'CI95_Low', ...
        'CI95_High', ...
        'P_Raw', ...
        'Q_BH_FDR', ...
        'HedgesType_gAdjusted'}));

end

fprintf('\nNominal p<0.05 edges are saved separately and are NOT FDR-significant\n');
fprintf('unless their q-value is also <0.05.\n');

fprintf('\nSaved files:\n');
fprintf('All 165 tests       : %s\n',allFile);
fprintf('FDR-significant     : %s\n',sigFile);
fprintf('Nominal p<0.05      : %s\n',nominalFile);
fprintf('Band summary        : %s\n',summaryFile);
fprintf('Full MAT            : %s\n',resultFile);
fprintf('Analysis log        : %s\n',logFile);

fprintf('\n============================================================\n');
fprintf('DO NOT REDEFINE THE BH-FDR FAMILY AFTER VIEWING RESULTS.\n');
fprintf('============================================================\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function st = univariate_group_effect_hc3(y,X)

    y = double(y(:));

    N = numel(y);
    p = size(X,2);

    beta = X\y;

    residual = y-X*beta;

    XtXinv = pinv(X'*X);

    Hdiag = sum((X*XtXinv).*X,2);

    if any(Hdiag >= 1)
        error('Invalid leverage encountered in HC3 calculation.');
    end

    uAdj2 = residual.^2 ./ (1-Hdiag).^2;

    meat = X'*(X.*uAdj2);

    Vhc3 = XtXinv*meat*XtXinv;

    varGroup = Vhc3(2,2);

    if ~isfinite(varGroup) || varGroup<=0
        error('Invalid HC3 group variance.');
    end

    se = sqrt(varGroup);

    tStat = beta(2)/se;

    df = N-p;

    pValue = t_two_sided_p(tStat,df);

    tcrit = student_t_critical_975(df);

    ciLow  = beta(2)-tcrit*se;
    ciHigh = beta(2)+tcrit*se;

    % Descriptive standardized adjusted contrast:
    % adjusted Group beta / full-model residual SD with small-sample J.
    rmse = sqrt(sum(residual.^2)/df);

    dAdj = beta(2)/rmse;

    J = 1-3/(4*df-1);

    gAdj = J*dAdj;

    st = struct();

    st.betaGroup = beta(2);
    st.se = se;

    st.t = tStat;
    st.df = df;
    st.p = pValue;

    st.ciLow = ciLow;
    st.ciHigh = ciHigh;

    st.gAdjusted = gAdj;

end

%% =========================================================================
% BENJAMINI-HOCHBERG ADJUSTED P-VALUES
%% =========================================================================

function q = bh_adjust(p)

    p = double(p(:));

    m = numel(p);

    if any(~isfinite(p)) || any(p<0) || any(p>1)
        error('Invalid p-values supplied to BH-FDR.');
    end

    [pSort,ord] = sort(p,'ascend');

    qSort = nan(m,1);

    % Raw BH adjusted values.
    raw = pSort .* m ./ (1:m)';

    % Enforce monotonicity from largest to smallest rank.
    runningMin = 1;

    for i = m:-1:1

        runningMin = min(runningMin,raw(i));

        qSort(i) = min(1,runningMin);

    end

    q = nan(m,1);

    q(ord) = qSort;

end

%% =========================================================================
% TWO-SIDED STUDENT-t P WITHOUT STATISTICS TOOLBOX
%% =========================================================================

function p = t_two_sided_p(t,df)

    if ~isfinite(t) || df<=0
        error('Invalid t statistic or degrees of freedom.');
    end

    x = df/(df+t^2);

    % For t>=0:
    % two-sided P = I_x(df/2, 1/2)
    p = betainc(x,df/2,0.5);

    p = max(0,min(1,p));

end

%% =========================================================================
% 95% t CRITICAL VALUE WITHOUT STATISTICS TOOLBOX
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

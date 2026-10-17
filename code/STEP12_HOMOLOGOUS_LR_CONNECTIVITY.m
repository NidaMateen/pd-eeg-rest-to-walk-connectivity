%% =========================================================================
% PAPER 3 — STEP 12
% PRESPECIFIED HOMOLOGOUS LEFT-RIGHT CONNECTIVITY ANALYSIS
%
% File:
%   STEP12_HOMOLOGOUS_LR_CONNECTIVITY.m
%
% PURPOSE
% -------
% Test the prespecified homologous left-right regional connections:
%
%   1) FrontalLeft  <-> FrontalRight
%   2) CentralLeft  <-> CentralRight
%   3) ParietalLeft <-> ParietalRight
%   4) TemporalLeft <-> TemporalRight
%
% across theta, alpha, and beta Rest->Walk reconfiguration.
%
% OUTCOME
% -------
%   Delta db-wPLI = Walk - Rest
%
% INFERENTIAL STRUCTURE
% ---------------------
% For each band separately, the four homologous edges are modeled jointly:
%
%   [FL-FR, CL-CR, PL-PR, TL-TR] ~ Group + Age_c + Sex
%
% Coding:
%   Group: HC=0, PD=1
%   Sex:   F=0, M=1
%
% A 4-df subject-level multivariate HC3 robust Wald test evaluates:
%
%   H0: all four adjusted PD-HC homologous-edge effects = 0
%
% producing ONE omnibus p-value per band.
%
% MULTIPLE-COMPARISON CONTROL
% ---------------------------
% The three band-level omnibus p-values (theta, alpha, beta) are corrected
% together using Holm correction.
%
% The four individual homologous-edge adjusted effects within each band are
% reported with HC3 95% CIs for anatomical interpretation. They are not used
% to redefine the omnibus test.
%
% COHORT
% ------
% All four homologous edges involve only lateralized regions that were usable
% in all 130 connectivity participants, so expected:
%
%   N = 130 = 23 HC + 107 PD
%
% IMPORTANT
% ---------
% This is the frozen Level-4 targeted homologous L-R analysis.
% Do not add/remove homologous edges after results are viewed.
%% =========================================================================

clear;
clc;
close all;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% =========================================================================
% 1. FROZEN SETTINGS
%% =========================================================================

EXPECTED_N  = 130;
EXPECTED_HC = 23;
EXPECTED_PD = 107;

BIO_BANDS = ["theta";"alpha";"beta"];
ALPHA = 0.05;

PAIR_REGION1 = [ ...
    "FrontalLeft"
    "CentralLeft"
    "ParietalLeft"
    "TemporalLeft"];

PAIR_REGION2 = [ ...
    "FrontalRight"
    "CentralRight"
    "ParietalRight"
    "TemporalRight"];

PAIR_LABEL = [ ...
    "Frontal_LR"
    "Central_LR"
    "Parietal_LR"
    "Temporal_LR"];

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
    'PAPER3_HOMOLOGOUS_LR');

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

if numel(subjectsAll) ~= EXPECTED_N || ...
   sum(groupsAll=="HC") ~= EXPECTED_HC || ...
   sum(groupsAll=="PD") ~= EXPECTED_PD

    error('Expected full cohort 130 = 23 HC + 107 PD.');
end

if size(deltaBio,2) ~= 55
    error('Expected 55 regional edges.');
end

%% =========================================================================
% 4. IDENTIFY BAND INDICES
%% =========================================================================

bandIdx = nan(3,1);

for b = 1:3
    idx = find(bandNamesAll==BIO_BANDS(b));
    if numel(idx) ~= 1
        error('Could not uniquely identify band %s.',BIO_BANDS(b));
    end
    bandIdx(b) = idx;
end

%% =========================================================================
% 5. IDENTIFY EXACT HOMOLOGOUS EDGES
%% =========================================================================

homEdgeIdx = nan(4,1);

for h = 1:4

    mask = ...
        (edgeRegion1==PAIR_REGION1(h) & edgeRegion2==PAIR_REGION2(h)) | ...
        (edgeRegion1==PAIR_REGION2(h) & edgeRegion2==PAIR_REGION1(h));

    idx = find(mask);

    if numel(idx) ~= 1
        error( ...
            'Expected one edge for %s <-> %s; found %d.', ...
            PAIR_REGION1(h),PAIR_REGION2(h),numel(idx));
    end

    homEdgeIdx(h) = idx;

end

fprintf('\n============================================================\n');
fprintf('PRESPECIFIED HOMOLOGOUS LEFT-RIGHT EDGES\n');
fprintf('============================================================\n');

for h = 1:4
    fprintf('%-12s : %s\n',PAIR_LABEL(h),edgeNames(homEdgeIdx(h)));
end

%% =========================================================================
% 6. COHORT CHECK
%% =========================================================================

finiteMask = true(numel(subjectsAll),1);

for b = 1:3
    tmp = deltaBio(:,homEdgeIdx,bandIdx(b));
    finiteMask = finiteMask & all(isfinite(tmp),2);
end

subjects = subjectsAll(finiteMask);
groups   = groupsAll(finiteMask);

N   = numel(subjects);
nHC = sum(groups=="HC");
nPD = sum(groups=="PD");

fprintf('\n============================================================\n');
fprintf('HOMOLOGOUS L-R COHORT\n');
fprintf('============================================================\n');
fprintf('N total = %d\n',N);
fprintf('N HC    = %d\n',nHC);
fprintf('N PD    = %d\n',nPD);

if N ~= EXPECTED_N || nHC ~= EXPECTED_HC || nPD ~= EXPECTED_PD
    fprintf('Excluded subjects:\n');
    disp(subjectsAll(~finiteMask));
    error('Homologous L-R cohort does not match expected 130/23/107.');
end

delta = deltaBio(finiteMask,:,:);

%% =========================================================================
% 7. LOAD / ALIGN AGE AND SEX
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

[tf,loc] = ismember(subjects,string(P.participant_id));

if ~all(tf)
    error('Participant alignment failed.');
end

participantGroup = string(P.group(loc));

if ~isequal(participantGroup,groups)
    error('Group mismatch between connectivity and participants.tsv.');
end

age = column_to_double(P.age(loc));

if any(~isfinite(age))
    error('Age contains missing/nonfinite values.');
end

sexRaw = upper(strtrim(string(P.sex(loc))));

if any(~ismember(sexRaw,["F","M"]))
    error('Sex must be F or M.');
end

Group = double(groups=="PD");
Sex   = double(sexRaw=="M");

Age_c = age-mean(age);

X = [ones(N,1),Group,Age_c,Sex];

if rank(X) ~= size(X,2)
    error('Design matrix is rank deficient.');
end

fprintf('Age centering mean = %.6f years\n',mean(age));
fprintf('Coding: HC=0, PD=1; F=0, M=1\n');

%% =========================================================================
% 8. BUILD N x 4 x 3 HOMOLOGOUS-EDGE DATA
%% =========================================================================

Yhom = nan(N,4,3);

for b = 1:3
    Yhom(:,:,b) = squeeze(delta(:,homEdgeIdx,bandIdx(b)));
end

if any(~isfinite(Yhom),'all')
    error('Unexpected nonfinite homologous-edge values.');
end

%% =========================================================================
% 9. BAND-LEVEL 4-DF OMNIBUS TESTS
%% =========================================================================

Band = BIO_BANDS;

WaldChi2 = nan(3,1);
Frobust  = nan(3,1);
DF1      = repmat(4,3,1);
DF2      = nan(3,1);
PRaw     = nan(3,1);

for b = 1:3

    Y = Yhom(:,:,b);

    st = multivariate_group_hc3(Y,X);

    WaldChi2(b) = st.Wald;
    Frobust(b)  = st.F;
    DF2(b)      = st.df2;
    PRaw(b)     = st.p;

end

PHolm = holm_adjust(PRaw);
SignificantHolm = PHolm < ALPHA;

OmnibusTable = table( ...
    Band, ...
    WaldChi2, ...
    Frobust, ...
    DF1, ...
    DF2, ...
    PRaw, ...
    PHolm, ...
    SignificantHolm);

%% =========================================================================
% 10. EDGE-SPECIFIC ADJUSTED EFFECTS FOR INTERPRETATION
%% =========================================================================

nRows = 3*4;

BandE = strings(nRows,1);
HomologousPair = strings(nRows,1);
Edge = strings(nRows,1);

MeanDelta_HC = nan(nRows,1);
MeanDelta_PD = nan(nRows,1);

AdjustedBeta_PDminusHC = nan(nRows,1);
HC3_SE = nan(nRows,1);
CI95_Low = nan(nRows,1);
CI95_High = nan(nRows,1);

P_Raw_Descriptive = nan(nRows,1);
HedgesType_gAdjusted = nan(nRows,1);

row = 0;

for b = 1:3

    for h = 1:4

        row = row+1;

        y = Yhom(:,h,b);

        st = univariate_group_effect_hc3(y,X);

        BandE(row) = BIO_BANDS(b);
        HomologousPair(row) = PAIR_LABEL(h);
        Edge(row) = edgeNames(homEdgeIdx(h));

        MeanDelta_HC(row) = mean(y(Group==0));
        MeanDelta_PD(row) = mean(y(Group==1));

        AdjustedBeta_PDminusHC(row) = st.betaGroup;
        HC3_SE(row) = st.se;
        CI95_Low(row) = st.ciLow;
        CI95_High(row) = st.ciHigh;
        P_Raw_Descriptive(row) = st.p;
        HedgesType_gAdjusted(row) = st.gAdjusted;

    end

end

EdgeEffectTable = table( ...
    BandE, ...
    HomologousPair, ...
    Edge, ...
    MeanDelta_HC, ...
    MeanDelta_PD, ...
    AdjustedBeta_PDminusHC, ...
    HC3_SE, ...
    CI95_Low, ...
    CI95_High, ...
    P_Raw_Descriptive, ...
    HedgesType_gAdjusted, ...
    'VariableNames', ...
    {'Band', ...
     'HomologousPair', ...
     'Edge', ...
     'MeanDelta_HC', ...
     'MeanDelta_PD', ...
     'AdjustedBeta_PDminusHC', ...
     'HC3_SE', ...
     'CI95_Low', ...
     'CI95_High', ...
     'P_Raw_Descriptive', ...
     'HedgesType_gAdjusted'});

%% =========================================================================
% 11. SAVE SUBJECT-LEVEL LONG TABLE
%% =========================================================================

nLong = N*4*3;

SubjectCol = strings(nLong,1);
GroupCol = strings(nLong,1);
BandCol = strings(nLong,1);
PairCol = strings(nLong,1);
EdgeCol = strings(nLong,1);

AgeCol = nan(nLong,1);
SexCol = strings(nLong,1);
DeltaCol = nan(nLong,1);

rr = 0;

for b = 1:3
    for h = 1:4
        for i = 1:N

            rr = rr+1;

            SubjectCol(rr) = subjects(i);
            GroupCol(rr) = groups(i);
            BandCol(rr) = BIO_BANDS(b);
            PairCol(rr) = PAIR_LABEL(h);
            EdgeCol(rr) = edgeNames(homEdgeIdx(h));

            AgeCol(rr) = age(i);
            SexCol(rr) = sexRaw(i);
            DeltaCol(rr) = Yhom(i,h,b);

        end
    end
end

SubjectTable = table( ...
    SubjectCol, ...
    GroupCol, ...
    AgeCol, ...
    SexCol, ...
    BandCol, ...
    PairCol, ...
    EdgeCol, ...
    DeltaCol, ...
    'VariableNames', ...
    {'Subject','Group','Age','Sex','Band','HomologousPair','Edge','Delta_dbWPLI'});

%% =========================================================================
% 12. SAVE FILES
%% =========================================================================

omnibusFile = fullfile( ...
    outRoot, ...
    'Paper3_HOMOLOGOUS_LR_band_omnibus_Holm.csv');

effectFile = fullfile( ...
    outRoot, ...
    'Paper3_HOMOLOGOUS_LR_edge_effects.csv');

subjectFile = fullfile( ...
    outRoot, ...
    'Paper3_HOMOLOGOUS_LR_subject_data.csv');

resultFile = fullfile( ...
    outRoot, ...
    'Paper3_HOMOLOGOUS_LR_RESULTS.mat');

writetable(OmnibusTable,omnibusFile);
writetable(EdgeEffectTable,effectFile);
writetable(SubjectTable,subjectFile);

R = struct();

R.subjects = subjects;
R.groups = groups;

R.age = age;
R.Age_c = Age_c;
R.sexRaw = sexRaw;
R.Group = Group;
R.Sex = Sex;

R.bandNames = BIO_BANDS;

R.homologousPairLabels = PAIR_LABEL;
R.homologousEdgeIdx = homEdgeIdx;
R.homologousEdgeNames = edgeNames(homEdgeIdx);

R.Yhom = Yhom;

R.OmnibusTable = OmnibusTable;
R.EdgeEffectTable = EdgeEffectTable;

save(resultFile,'-struct','R','-v7.3');

%% =========================================================================
% 13. LOG
%% =========================================================================

logFile = fullfile( ...
    outRoot, ...
    'Paper3_HOMOLOGOUS_LR_analysis_log.txt');

fid = fopen(logFile,'w');

if fid == -1
    error('Could not open log file.');
end

fprintf(fid,'PAPER 3 PRESPECIFIED HOMOLOGOUS LEFT-RIGHT ANALYSIS\n');
fprintf(fid,'Run time: %s\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));

fprintf(fid,'Cohort: N=%d; HC=%d; PD=%d\n',N,nHC,nPD);
fprintf(fid,'Age center: %.10f\n',mean(age));

fprintf(fid,'Edges:\n');
for h = 1:4
    fprintf(fid,'  %s = %s\n',PAIR_LABEL(h),edgeNames(homEdgeIdx(h)));
end

fprintf(fid,'\nOutcome: Delta db-wPLI = Walk - Rest\n');
fprintf(fid,'Band omnibus model: 4 homologous outcomes ~ Group + Age_c + Sex\n');
fprintf(fid,'Inference: 4-df subject-level multivariate HC3 robust Wald F approximation\n');
fprintf(fid,'Multiplicity: Holm across theta/alpha/beta omnibus tests\n');
fprintf(fid,'Individual edge estimates: descriptive/anatomical interpretation\n');

fclose(fid);

%% =========================================================================
% 14. CONSOLE OUTPUT
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('PRESPECIFIED HOMOLOGOUS L-R RESULTS\n');
fprintf('============================================================\n');

disp(OmnibusTable);

fprintf('\nHomologous-edge adjusted PD-HC effects:\n');
disp(EdgeEffectTable);

fprintf('\nInterpretation rule:\n');

if any(SignificantHolm)

    fprintf([ ...
        'At least one frequency band shows a Holm-significant joint PD-HC\n' ...
        'effect across the four prespecified homologous L-R connections.\n' ...
        'Use the individual edge estimates to describe anatomical pattern.\n']);

else

    fprintf([ ...
        'No frequency band shows a Holm-significant joint homologous L-R\n' ...
        'effect. This does not negate the primary NBS result.\n']);

end

fprintf('\nSaved files:\n');
fprintf('Band omnibus       : %s\n',omnibusFile);
fprintf('Edge effects       : %s\n',effectFile);
fprintf('Subject data       : %s\n',subjectFile);
fprintf('Full MAT           : %s\n',resultFile);
fprintf('Analysis log       : %s\n',logFile);

fprintf('\n============================================================\n');
fprintf('DO NOT ALTER HOMOLOGOUS EDGE SET AFTER VIEWING RESULTS.\n');
fprintf('============================================================\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function st = multivariate_group_hc3(Y,X)

    Y = double(Y);

    [N,m] = size(Y);
    p = size(X,2);

    XtXinv = pinv(X'*X);

    B = XtXinv*X'*Y;

    E = Y-X*B;

    Hdiag = sum((X*XtXinv).*X,2);

    if any(Hdiag >= 1)
        error('Invalid leverage in HC3 multivariate calculation.');
    end

    % Influence scalar for Group coefficient (row 2 of B).
    a = (XtXinv(2,:)*X')';

    Vg = zeros(m,m);

    for i = 1:N

        eAdj = E(i,:)/(1-Hdiag(i));

        Vg = Vg + ...
            (a(i)^2)*(eAdj'*eAdj);

    end

    betaGroup = B(2,:)';

    Wald = betaGroup'*pinv(Vg)*betaGroup;

    q = m;

    F = Wald/q;

    df1 = q;
    df2 = N-p;

    pValue = f_upper_tail(F,df1,df2);

    st = struct();

    st.betaGroup = betaGroup;
    st.Wald = Wald;
    st.F = F;
    st.df1 = df1;
    st.df2 = df2;
    st.p = pValue;

end

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
        error('Invalid leverage in HC3 calculation.');
    end

    uAdj2 = residual.^2 ./ (1-Hdiag).^2;

    meat = X'*(X.*uAdj2);

    Vhc3 = XtXinv*meat*XtXinv;

    se = sqrt(Vhc3(2,2));

    df = N-p;

    tStat = beta(2)/se;
    pValue = t_two_sided_p(tStat,df);

    tcrit = student_t_critical_975(df);

    ciLow = beta(2)-tcrit*se;
    ciHigh = beta(2)+tcrit*se;

    rmse = sqrt(sum(residual.^2)/df);

    dAdj = beta(2)/rmse;

    J = 1-3/(4*df-1);

    gAdj = J*dAdj;

    st = struct();

    st.betaGroup = beta(2);
    st.se = se;
    st.ciLow = ciLow;
    st.ciHigh = ciHigh;
    st.p = pValue;
    st.gAdjusted = gAdj;

end

%% =========================================================================

function pAdj = holm_adjust(p)

    p = double(p(:));
    m = numel(p);

    [pSort,ord] = sort(p,'ascend');

    adjSort = nan(m,1);
    runningMax = 0;

    for i = 1:m

        rawAdj = (m-i+1)*pSort(i);

        runningMax = max(runningMax,rawAdj);

        adjSort(i) = min(1,runningMax);

    end

    pAdj = nan(m,1);
    pAdj(ord) = adjSort;

end

%% =========================================================================

function p = f_upper_tail(F,df1,df2)

    if ~isfinite(F) || F<0
        error('Invalid F statistic.');
    end

    x = (df1*F)/(df1*F+df2);

    cdf = betainc(x,df1/2,df2/2);

    p = max(0,min(1,1-cdf));

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

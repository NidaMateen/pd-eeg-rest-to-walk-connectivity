%% =========================================================================
% PAPER 3 — STEP 15B
% CROSS-COHORT CLINICAL-STATUS ASSOCIATIONS (HC + PD)
%
% File:
%   STEP15B_CROSSCOHORT_CLINICAL_STATUS_ASSOCIATIONS.m
%
% PURPOSE
% -------
% Test whether subject-level theta Rest->Walk network reconfiguration is
% associated with clinical status ACROSS the combined HC + PD cohort.
%
% This analysis is deliberately separated from STEP 15 (PD-only severity).
%
% PRIMARY CROSS-COHORT CLINICAL VARIABLES
% ---------------------------------------
% Only variables that are meaningfully measured in BOTH HC and PD:
%
%   1) MDS-UPDRS Part III
%   2) MoCA
%
% Disease duration and LEDD are NOT analyzed across HC + PD because they are
% PD-specific exposures/characteristics and are undefined in healthy controls.
%
% PRIMARY NEURAL OUTCOME
% ----------------------
% Mean theta Delta db-wPLI across ALL 55 regional edges:
%
%       Delta theta = Walk - Rest
%
% using the complete-network cohort:
%
%       Expected N = 117 = 22 HC + 95 PD
%
% before clinical-score missingness.
%
% ANALYSES FOR EACH CLINICAL VARIABLE
% -----------------------------------
% A) Pooled Spearman correlation across HC + PD
%    - descriptive "clinical spectrum" association
%    - 10,000-permutation two-sided P value
%    - BH-FDR across the 2 clinical variables
%
% B) Group-adjusted HC3 regression
%
%       ThetaDelta ~ ClinicalScore + Group + Age_c + Sex
%
%    Group: HC=0, PD=1
%    Sex:   F=0, M=1
%
%    This asks whether the clinical score is associated with theta
%    reconfiguration BEYOND the mean PD-vs-HC group difference.
%
% C) Group x ClinicalScore interaction
%
%       ThetaDelta ~ ClinicalScore_z + Group
%                    + ClinicalScore_z*Group
%                    + Age_c + Sex
%
%    This asks whether the score-reconfiguration slope differs between HC
%    and PD.
%
% IMPORTANT:
% For UPDRS-III, HC and PD score ranges may have little or no overlap.
% Therefore the script explicitly quantifies score-range overlap and VIF.
% If common support is absent, the pooled Spearman correlation may largely
% reflect group separation, and the interaction should be interpreted with
% caution.
%
% D) Within-group Spearman correlations
%    - HC separately
%    - PD separately
%    These are descriptive context only.
%
% MULTIPLICITY
% ------------
% Separate BH-FDR correction across the 2 clinical variables for:
%   1) pooled Spearman tests
%   2) group-adjusted clinical-score coefficients
%   3) Group x ClinicalScore interaction coefficients
%
% SENSITIVITY OUTCOME
% -------------------
% The same analysis is repeated using mean theta Delta across the already-
% significant 32-edge primary NBS component.
%
% This NBS32 analysis is explicitly POST-SELECTION / sensitivity only.
%
% STATUS
% ------
% Exploratory cross-cohort clinical-status validation.
% It does not alter the frozen primary PD-vs-HC NBS hypothesis test.
%% =========================================================================

clear;
clc;
close all;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% =========================================================================
% 1. FROZEN SETTINGS
%% =========================================================================

EXPECTED_FULL_N  = 130;
EXPECTED_FULL_HC = 23;
EXPECTED_FULL_PD = 107;

EXPECTED_COMPLETE_N  = 117;
EXPECTED_COMPLETE_HC = 22;
EXPECTED_COMPLETE_PD = 95;

EXPECTED_PRIMARY_COMPONENT_EDGES = 32;

N_PERM_SPEARMAN = 10000;
RNG_SEED = 20260917;
FDR_ALPHA = 0.05;

clinicalVars = [ ...
    "updrs_part_iii"
    "moca"];

clinicalLabels = [ ...
    "MDS-UPDRS Part III"
    "MoCA"];

clinicalUnits = [ ...
    "points"
    "points"];

rng(RNG_SEED,'twister');

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

outRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_CROSSCOHORT_CLINICAL_ASSOCIATIONS');

statsInputFile = fullfile( ...
    connectivityRoot, ...
    'Paper3_CONNECTIVITY_STATS_INPUT.mat');

participantsFile = fullfile( ...
    datasetRoot, ...
    'participants.tsv');

primaryNBSFile = fullfile( ...
    primaryNBSRoot, ...
    'Paper3_PRIMARY_NBS_RESULTS.mat');

if ~exist(outRoot,'dir')
    mkdir(outRoot);
end

requiredFiles = { ...
    statsInputFile, ...
    participantsFile, ...
    primaryNBSFile};

for k = 1:numel(requiredFiles)

    if ~exist(requiredFiles{k},'file')
        error('Missing required file:\n%s',requiredFiles{k});
    end

end

%% =========================================================================
% 3. LOAD CONNECTIVITY INPUT
%% =========================================================================

S = load(statsInputFile);

requiredS = { ...
    'subjects', ...
    'groups', ...
    'edgeNames', ...
    'bandNames', ...
    'deltaBio', ...
    'completeNetworkMask'};

for k = 1:numel(requiredS)

    if ~isfield(S,requiredS{k})
        error('Statistics input lacks field: %s',requiredS{k});
    end

end

subjectsAll = string(S.subjects(:));
groupsAll   = string(S.groups(:));

edgeNames = string(S.edgeNames(:));
bandNames = lower(string(S.bandNames(:)));

deltaBio = double(S.deltaBio);
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

if size(deltaBio,2) ~= 55
    error('Expected 55 regional connectivity edges.');
end

thetaIdx = find(bandNames=="theta");

if numel(thetaIdx) ~= 1
    error('Could not uniquely identify theta band.');
end

%% =========================================================================
% 4. LOAD PRIMARY 32-EDGE THETA NBS COMPONENT
%% =========================================================================

R = load(primaryNBSFile);

if ~isfield(R,'ComponentTable') || ...
   ~isfield(R,'observedComponents')

    error('Primary NBS MAT lacks required component fields.');

end

CT = R.ComponentTable;

sigThetaNeg = ...
    logical(CT.Significant_FWER_0p05) & ...
    lower(string(CT.Band))=="theta" & ...
    contains(lower(string(CT.Direction)),"negative");

if sum(sigThetaNeg) ~= 1
    error('Expected exactly one significant theta-negative NBS component.');
end

sigRow = CT(sigThetaNeg,:);

componentID = double(sigRow.ComponentID);
componentP  = double(sigRow.JointFWER_P);

thetaNegativeComponents = R.observedComponents{1,2};

if componentID < 1 || componentID > numel(thetaNegativeComponents)
    error('Invalid primary theta component ID.');
end

primaryEdgeIdx = double( ...
    thetaNegativeComponents(componentID).edgeIdx(:));

if numel(primaryEdgeIdx) ~= EXPECTED_PRIMARY_COMPONENT_EDGES
    error('Expected 32 edges in the primary theta component.');
end

fprintf('\n============================================================\n');
fprintf('PRIMARY NBS COMPONENT LOCK\n');
fprintf('============================================================\n');

fprintf('Primary theta-negative component = %d edges\n', ...
    numel(primaryEdgeIdx));

fprintf('Primary NBS joint FWER p         = %.8f\n', ...
    componentP);

%% =========================================================================
% 5. LOAD PARTICIPANTS.TSV
%% =========================================================================

opts = detectImportOptions( ...
    participantsFile, ...
    'FileType','text', ...
    'Delimiter','\t');

allVars = string(opts.VariableNames);

for k = 1:numel(allVars)

    try

        opts = setvaropts( ...
            opts, ...
            allVars(k), ...
            'TreatAsMissing', ...
            {'n/a','N/A','NA','na','NaN','nan',''});

    catch
        % Ignore variable-option classes that do not expose TreatAsMissing.
    end

end

P = readtable(participantsFile,opts);

requiredP = [ ...
    "participant_id"
    "group"
    "age"
    "sex"
    clinicalVars];

for k = 1:numel(requiredP)

    if ~ismember(requiredP(k),string(P.Properties.VariableNames))

        error('participants.tsv lacks required column: %s', ...
            requiredP(k));

    end

end

%% =========================================================================
% 6. ALIGN PARTICIPANT METADATA TO CONNECTIVITY ORDER
%% =========================================================================

participantIDs = string(P.participant_id);

[tf,loc] = ismember(subjectsAll,participantIDs);

if ~all(tf)

    missingIDs = subjectsAll(~tf);

    error( ...
        'Could not match %d connectivity subjects. First missing: %s', ...
        sum(~tf), ...
        missingIDs(1));

end

P = P(loc,:);

if any(string(P.group) ~= groupsAll)

    bad = find(string(P.group) ~= groupsAll,1);

    error('Group mismatch for subject %s.',subjectsAll(bad));

end

ageAll = column_to_double(P.age);

sexLabelAll = upper(strtrim(string(P.sex)));

if any(~ismember(sexLabelAll,["M","F"]))

    bad = find(~ismember(sexLabelAll,["M","F"]),1);

    error( ...
        'Unexpected sex coding for %s: %s', ...
        subjectsAll(bad), ...
        sexLabelAll(bad));

end

sexIntAll = double(sexLabelAll=="M");

clinicalAll = nan(EXPECTED_FULL_N,numel(clinicalVars));

for c = 1:numel(clinicalVars)

    clinicalAll(:,c) = column_to_double( ...
        P.(clinicalVars(c)));

end

%% =========================================================================
% 7. COMPLETE-NETWORK HC + PD COHORT
%% =========================================================================

subjects = subjectsAll(completeMask);
groups = groupsAll(completeMask);

groupInt = double(groups=="PD");

age = ageAll(completeMask);
sexInt = sexIntAll(completeMask);
sexLabel = sexLabelAll(completeMask);

clinical = clinicalAll(completeMask,:);

if numel(subjects) ~= EXPECTED_COMPLETE_N || ...
   sum(groups=="HC") ~= EXPECTED_COMPLETE_HC || ...
   sum(groups=="PD") ~= EXPECTED_COMPLETE_PD

    error('Complete-network cohort mismatch after metadata alignment.');

end

%% =========================================================================
% 8. SUBJECT-LEVEL THETA FEATURES
%% =========================================================================

thetaComplete = squeeze( ...
    deltaBio(completeMask,:,thetaIdx));

if size(thetaComplete,1) ~= EXPECTED_COMPLETE_N || ...
   size(thetaComplete,2) ~= 55

    error('Unexpected theta complete-network matrix dimensions.');

end

if any(~isfinite(thetaComplete(:)))

    error('Complete-network theta matrix contains nonfinite values.');

end

whole55ThetaDelta = ...
    mean(thetaComplete,2);

nbs32ThetaDelta = ...
    mean(thetaComplete(:,primaryEdgeIdx),2);

fprintf('\n============================================================\n');
fprintf('CROSS-COHORT ANALYSIS COHORT\n');
fprintf('============================================================\n');

fprintf('Complete-network N  = %d\n',numel(subjects));
fprintf('HC = %d | PD = %d\n',sum(groups=="HC"),sum(groups=="PD"));

fprintf('\nWhole55 theta Delta means:\n');
fprintf('HC: %.6f +/- %.6f\n', ...
    mean(whole55ThetaDelta(groups=="HC")), ...
    std(whole55ThetaDelta(groups=="HC"),0));

fprintf('PD: %.6f +/- %.6f\n', ...
    mean(whole55ThetaDelta(groups=="PD")), ...
    std(whole55ThetaDelta(groups=="PD"),0));

fprintf('\nNBS32 theta Delta means:\n');
fprintf('HC: %.6f +/- %.6f\n', ...
    mean(nbs32ThetaDelta(groups=="HC")), ...
    std(nbs32ThetaDelta(groups=="HC"),0));

fprintf('PD: %.6f +/- %.6f\n', ...
    mean(nbs32ThetaDelta(groups=="PD")), ...
    std(nbs32ThetaDelta(groups=="PD"),0));

%% =========================================================================
% 9. SCORE DISTRIBUTION / COMMON-SUPPORT DIAGNOSTICS
%% =========================================================================

ScoreDistribution = build_score_distribution_table( ...
    clinical, ...
    groups, ...
    clinicalLabels, ...
    clinicalUnits);

scoreDistFile = fullfile( ...
    outRoot, ...
    'Paper3_CROSSCOHORT_clinical_score_distributions.csv');

writetable(ScoreDistribution,scoreDistFile);

%% =========================================================================
% 10. PRIMARY OUTCOME:
%     WHOLE-NETWORK 55-EDGE THETA DELTA
%% =========================================================================

PrimaryResults = run_crosscohort_family( ...
    whole55ThetaDelta, ...
    clinical, ...
    age, ...
    sexInt, ...
    groups, ...
    groupInt, ...
    clinicalVars, ...
    clinicalLabels, ...
    clinicalUnits, ...
    "Whole-network mean theta Delta (55 edges)", ...
    N_PERM_SPEARMAN, ...
    FDR_ALPHA);

primaryFile = fullfile( ...
    outRoot, ...
    'Paper3_CROSSCOHORT_PRIMARY_whole55_theta.csv');

writetable(PrimaryResults,primaryFile);

%% =========================================================================
% 11. SENSITIVITY OUTCOME:
%     POST-SELECTION 32-EDGE NBS COMPONENT
%% =========================================================================

SensitivityResults = run_crosscohort_family( ...
    nbs32ThetaDelta, ...
    clinical, ...
    age, ...
    sexInt, ...
    groups, ...
    groupInt, ...
    clinicalVars, ...
    clinicalLabels, ...
    clinicalUnits, ...
    "Primary NBS component mean theta Delta (32 edges)", ...
    N_PERM_SPEARMAN, ...
    FDR_ALPHA);

sensitivityFile = fullfile( ...
    outRoot, ...
    'Paper3_CROSSCOHORT_SENSITIVITY_NBS32_theta.csv');

writetable(SensitivityResults,sensitivityFile);

%% =========================================================================
% 12. SUBJECT-LEVEL ANALYSIS FILE
%% =========================================================================

Subject = subjects;
Group = groups;
Age = age;
Sex = sexLabel;

UPDRS_III = clinical(:,1);
MoCA = clinical(:,2);

Whole55Theta_Delta = whole55ThetaDelta;
NBS32Theta_Delta = nbs32ThetaDelta;

SubjectLevelTable = table( ...
    Subject, ...
    Group, ...
    Age, ...
    Sex, ...
    UPDRS_III, ...
    MoCA, ...
    Whole55Theta_Delta, ...
    NBS32Theta_Delta);

subjectFile = fullfile( ...
    outRoot, ...
    'Paper3_CROSSCOHORT_subject_level_data.csv');

writetable(SubjectLevelTable,subjectFile);

%% =========================================================================
% 13. OUTCOME CONCORDANCE
%% =========================================================================

Pearson_r = pearson_corr( ...
    whole55ThetaDelta, ...
    nbs32ThetaDelta);

Spearman_rho = spearman_rho( ...
    whole55ThetaDelta, ...
    nbs32ThetaDelta);

Outcome1 = "Whole-network 55-edge theta Delta";
Outcome2 = "NBS32 theta Delta";
N = numel(subjects);

OutcomeConcordance = table( ...
    Outcome1, ...
    Outcome2, ...
    N, ...
    Pearson_r, ...
    Spearman_rho);

concordanceFile = fullfile( ...
    outRoot, ...
    'Paper3_CROSSCOHORT_outcome_concordance.csv');

writetable(OutcomeConcordance,concordanceFile);

%% =========================================================================
% 14. MANUSCRIPT-READY SUMMARY
%% =========================================================================

ClinicalVariable = PrimaryResults.ClinicalVariable;
N = PrimaryResults.N;

PooledSpearman_rho = PrimaryResults.PooledSpearman_rho;
PooledSpearman_Pperm = PrimaryResults.PooledSpearman_Pperm;
PooledSpearman_Q = PrimaryResults.PooledSpearman_Q_BH_FDR;

AdjustedClinical_StdBeta = ...
    PrimaryResults.AdjustedClinical_StdBeta;

AdjustedClinical_CI95_Low = ...
    PrimaryResults.AdjustedClinical_StdCI95_Low;

AdjustedClinical_CI95_High = ...
    PrimaryResults.AdjustedClinical_StdCI95_High;

AdjustedClinical_P = ...
    PrimaryResults.AdjustedClinical_P;

AdjustedClinical_Q = ...
    PrimaryResults.AdjustedClinical_Q_BH_FDR;

Interaction_StdBeta = ...
    PrimaryResults.Interaction_StdBeta;

Interaction_CI95_Low = ...
    PrimaryResults.Interaction_StdCI95_Low;

Interaction_CI95_High = ...
    PrimaryResults.Interaction_StdCI95_High;

Interaction_P = ...
    PrimaryResults.Interaction_P;

Interaction_Q = ...
    PrimaryResults.Interaction_Q_BH_FDR;

HC_Spearman_rho = ...
    PrimaryResults.HC_Spearman_rho;

PD_Spearman_rho = ...
    PrimaryResults.PD_Spearman_rho;

CommonSupport = ...
    PrimaryResults.CommonSupport;

Clinical_VIF = ...
    PrimaryResults.AdjustedClinical_VIF;

Group_VIF = ...
    PrimaryResults.AdjustedGroup_VIF;

ManuscriptSummary = table( ...
    ClinicalVariable, ...
    N, ...
    PooledSpearman_rho, ...
    PooledSpearman_Pperm, ...
    PooledSpearman_Q, ...
    AdjustedClinical_StdBeta, ...
    AdjustedClinical_CI95_Low, ...
    AdjustedClinical_CI95_High, ...
    AdjustedClinical_P, ...
    AdjustedClinical_Q, ...
    Interaction_StdBeta, ...
    Interaction_CI95_Low, ...
    Interaction_CI95_High, ...
    Interaction_P, ...
    Interaction_Q, ...
    HC_Spearman_rho, ...
    PD_Spearman_rho, ...
    CommonSupport, ...
    Clinical_VIF, ...
    Group_VIF);

manuscriptFile = fullfile( ...
    outRoot, ...
    'Paper3_CROSSCOHORT_MANUSCRIPT_SUMMARY.csv');

writetable(ManuscriptSummary,manuscriptFile);

%% =========================================================================
% 15. SAVE FULL MAT
%% =========================================================================

resultFile = fullfile( ...
    outRoot, ...
    'Paper3_CROSSCOHORT_CLINICAL_ASSOCIATIONS_RESULTS.mat');

C = struct();

C.settings = struct();
C.settings.clinicalVars = clinicalVars;
C.settings.clinicalLabels = clinicalLabels;
C.settings.N_PERM_SPEARMAN = N_PERM_SPEARMAN;
C.settings.RNG_SEED = RNG_SEED;
C.settings.FDR_ALPHA = FDR_ALPHA;

C.settings.primaryOutcome = ...
    "Whole-network mean theta Delta across all 55 edges";

C.settings.sensitivityOutcome = ...
    "Mean theta Delta across 32 primary NBS component edges";

C.settings.adjustedModel = ...
    "ThetaDelta ~ ClinicalScore + Group + Age_c + Sex, HC3";

C.settings.interactionModel = ...
    "ThetaDelta ~ ClinicalScore_z + Group + ClinicalScore_z*Group + Age_c + Sex, HC3";

C.settings.status = ...
    "Exploratory cross-cohort clinical-status analysis";

C.subjects = subjects;
C.groups = groups;
C.groupInt = groupInt;
C.age = age;
C.sexInt = sexInt;

C.clinical = clinical;

C.whole55ThetaDelta = whole55ThetaDelta;
C.nbs32ThetaDelta = nbs32ThetaDelta;
C.primaryEdgeIdx = primaryEdgeIdx;

C.ScoreDistribution = ScoreDistribution;
C.PrimaryResults = PrimaryResults;
C.SensitivityResults = SensitivityResults;
C.SubjectLevelTable = SubjectLevelTable;
C.OutcomeConcordance = OutcomeConcordance;
C.ManuscriptSummary = ManuscriptSummary;

save(resultFile,'-struct','C','-v7.3');

%% =========================================================================
% 16. WRITE ANALYSIS LOG
%% =========================================================================

logFile = fullfile( ...
    outRoot, ...
    'Paper3_CROSSCOHORT_clinical_analysis_log.txt');

fid = fopen(logFile,'w');

if fid==-1
    error('Could not create analysis log.');
end

fprintf(fid,'PAPER 3 — STEP 15B CROSS-COHORT CLINICAL-STATUS ASSOCIATIONS\n');
fprintf(fid,'Run time: %s\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));

fprintf(fid,'STATUS\n');
fprintf(fid,'Exploratory cross-cohort clinical-status validation.\n');
fprintf(fid,'Does not modify the frozen primary NBS inference.\n\n');

fprintf(fid,'COHORT\n');
fprintf(fid,'Complete network: N=%d = %d HC + %d PD before score missingness.\n\n', ...
    numel(subjects),sum(groups=="HC"),sum(groups=="PD"));

fprintf(fid,'PRIMARY NEURAL OUTCOME\n');
fprintf(fid,'Mean theta Delta db-wPLI across all 55 regional edges.\n');
fprintf(fid,'Delta = Walk - Rest.\n\n');

fprintf(fid,'CLINICAL VARIABLES\n');
fprintf(fid,'1. MDS-UPDRS Part III\n');
fprintf(fid,'2. MoCA\n');
fprintf(fid,'Disease duration and LEDD excluded because they are PD-specific.\n\n');

fprintf(fid,'ANALYTIC HIERARCHY FOR EACH SCORE\n');
fprintf(fid,'A. Pooled HC+PD Spearman rho with 10,000-permutation P.\n');
fprintf(fid,'B. HC3 model: ThetaDelta ~ Score + Group + Age_c + Sex.\n');
fprintf(fid,'C. HC3 interaction: ThetaDelta ~ Score_z + Group + Score_z*Group + Age_c + Sex.\n');
fprintf(fid,'D. HC and PD within-group Spearman rho as descriptive context.\n\n');

fprintf(fid,'MULTIPLICITY\n');
fprintf(fid,'Separate BH-FDR across UPDRS-III and MoCA for pooled Spearman,\n');
fprintf(fid,'group-adjusted score effects, and score-by-group interactions.\n\n');

fprintf(fid,'COMMON-SUPPORT WARNING\n');
fprintf(fid,['A pooled correlation can be driven by disease-group separation. ' ...
    'Range overlap and VIF are therefore reported explicitly.\n\n']);

fprintf(fid,'PRIMARY WHOLE55 RESULTS\n');

for c = 1:height(PrimaryResults)

    fprintf(fid,'\n%s\n',PrimaryResults.ClinicalVariable(c));

    fprintf(fid,'N=%d (HC=%d, PD=%d)\n', ...
        PrimaryResults.N(c), ...
        PrimaryResults.N_HC(c), ...
        PrimaryResults.N_PD(c));

    fprintf(fid,'HC range: %.4f to %.4f\n', ...
        PrimaryResults.HC_Min(c), ...
        PrimaryResults.HC_Max(c));

    fprintf(fid,'PD range: %.4f to %.4f\n', ...
        PrimaryResults.PD_Min(c), ...
        PrimaryResults.PD_Max(c));

    fprintf(fid,'Common support: %d | overlap width=%.4f\n', ...
        PrimaryResults.CommonSupport(c), ...
        PrimaryResults.OverlapWidth(c));

    fprintf(fid,'Pooled Spearman rho=%.4f; perm p=%.6g; q=%.6g\n', ...
        PrimaryResults.PooledSpearman_rho(c), ...
        PrimaryResults.PooledSpearman_Pperm(c), ...
        PrimaryResults.PooledSpearman_Q_BH_FDR(c));

    fprintf(fid,['Adjusted clinical std beta=%.4f; 95%% CI [%.4f, %.4f]; ' ...
        'p=%.6g; q=%.6g\n'], ...
        PrimaryResults.AdjustedClinical_StdBeta(c), ...
        PrimaryResults.AdjustedClinical_StdCI95_Low(c), ...
        PrimaryResults.AdjustedClinical_StdCI95_High(c), ...
        PrimaryResults.AdjustedClinical_P(c), ...
        PrimaryResults.AdjustedClinical_Q_BH_FDR(c));

    fprintf(fid,'Adjusted model VIF: clinical=%.3f; group=%.3f\n', ...
        PrimaryResults.AdjustedClinical_VIF(c), ...
        PrimaryResults.AdjustedGroup_VIF(c));

    fprintf(fid,['Interaction std beta=%.4f; 95%% CI [%.4f, %.4f]; ' ...
        'p=%.6g; q=%.6g\n'], ...
        PrimaryResults.Interaction_StdBeta(c), ...
        PrimaryResults.Interaction_StdCI95_Low(c), ...
        PrimaryResults.Interaction_StdCI95_High(c), ...
        PrimaryResults.Interaction_P(c), ...
        PrimaryResults.Interaction_Q_BH_FDR(c));

    fprintf(fid,'Within HC Spearman rho=%.4f | Within PD rho=%.4f\n', ...
        PrimaryResults.HC_Spearman_rho(c), ...
        PrimaryResults.PD_Spearman_rho(c));

end

fprintf(fid,'\n\nINTERPRETIVE GUARDRAILS\n');
fprintf(fid,'- Pooled Spearman is not evidence of within-PD severity tracking.\n');
fprintf(fid,'- If the adjusted score effect disappears after adding Group, the pooled\n');
fprintf(fid,'  association is primarily explained by PD-vs-HC separation.\n');
fprintf(fid,'- UPDRS interaction is weakly interpretable when HC/PD score ranges do not overlap.\n');
fprintf(fid,'- NBS32 results are post-selection sensitivity analyses only.\n');
fprintf(fid,'- Do not add further clinical variables after viewing these results.\n');

fclose(fid);

%% =========================================================================
% 17. CONSOLE REPORT
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('STEP 15B — CROSS-COHORT CLINICAL-STATUS ANALYSIS COMPLETE\n');
fprintf('============================================================\n');

fprintf('\nClinical score distributions / common support:\n');
disp(ScoreDistribution);

fprintf('\n------------------------------------------------------------\n');
fprintf('PRIMARY: WHOLE-NETWORK 55-EDGE THETA DELTA\n');
fprintf('------------------------------------------------------------\n');

disp(PrimaryResults(:,{ ...
    'ClinicalVariable', ...
    'N', ...
    'N_HC', ...
    'N_PD', ...
    'PooledSpearman_rho', ...
    'PooledSpearman_Pperm', ...
    'PooledSpearman_Q_BH_FDR', ...
    'AdjustedClinical_StdBeta', ...
    'AdjustedClinical_StdCI95_Low', ...
    'AdjustedClinical_StdCI95_High', ...
    'AdjustedClinical_P', ...
    'AdjustedClinical_Q_BH_FDR', ...
    'Interaction_StdBeta', ...
    'Interaction_StdCI95_Low', ...
    'Interaction_StdCI95_High', ...
    'Interaction_P', ...
    'Interaction_Q_BH_FDR', ...
    'HC_Spearman_rho', ...
    'PD_Spearman_rho', ...
    'CommonSupport', ...
    'AdjustedClinical_VIF', ...
    'AdjustedGroup_VIF'}));

fprintf('\n------------------------------------------------------------\n');
fprintf('SENSITIVITY: POST-SELECTION NBS32 THETA DELTA\n');
fprintf('------------------------------------------------------------\n');

disp(SensitivityResults(:,{ ...
    'ClinicalVariable', ...
    'N', ...
    'PooledSpearman_rho', ...
    'PooledSpearman_Pperm', ...
    'PooledSpearman_Q_BH_FDR', ...
    'AdjustedClinical_StdBeta', ...
    'AdjustedClinical_StdCI95_Low', ...
    'AdjustedClinical_StdCI95_High', ...
    'AdjustedClinical_P', ...
    'AdjustedClinical_Q_BH_FDR', ...
    'Interaction_StdBeta', ...
    'Interaction_P', ...
    'Interaction_Q_BH_FDR', ...
    'HC_Spearman_rho', ...
    'PD_Spearman_rho'}));

fprintf('\nOutcome concordance:\n');
disp(OutcomeConcordance);

fprintf('\nSaved:\n');
fprintf('Score distributions : %s\n',scoreDistFile);
fprintf('Primary results     : %s\n',primaryFile);
fprintf('Sensitivity results : %s\n',sensitivityFile);
fprintf('Subject-level data  : %s\n',subjectFile);
fprintf('Manuscript summary  : %s\n',manuscriptFile);
fprintf('Outcome concordance : %s\n',concordanceFile);
fprintf('Full MAT            : %s\n',resultFile);
fprintf('Analysis log        : %s\n',logFile);

fprintf('\n============================================================\n');
fprintf('STEP 15B FAMILY CLOSED.\n');
fprintf('Do not add further score-by-feature analyses after viewing results.\n');
fprintf('============================================================\n');

%% =========================================================================
%% LOCAL FUNCTION: CROSS-COHORT ANALYSIS FAMILY
%% =========================================================================

function T = run_crosscohort_family( ...
    yAll, ...
    clinicalMatrix, ...
    ageAll, ...
    sexAll, ...
    groupsAll, ...
    groupIntAll, ...
    clinicalVars, ...
    clinicalLabels, ...
    clinicalUnits, ...
    outcomeLabel, ...
    nPerm, ...
    fdrAlpha)

    nClinical = numel(clinicalVars);

    Outcome = repmat(string(outcomeLabel),nClinical,1);
    ClinicalVariable = clinicalLabels;
    VariableName = clinicalVars;
    Unit = clinicalUnits;

    N = nan(nClinical,1);
    N_HC = nan(nClinical,1);
    N_PD = nan(nClinical,1);

    HC_Min = nan(nClinical,1);
    HC_Max = nan(nClinical,1);
    PD_Min = nan(nClinical,1);
    PD_Max = nan(nClinical,1);

    OverlapLow = nan(nClinical,1);
    OverlapHigh = nan(nClinical,1);
    OverlapWidth = nan(nClinical,1);
    CommonSupport = false(nClinical,1);

    PooledSpearman_rho = nan(nClinical,1);
    PooledSpearman_Pperm = nan(nClinical,1);

    HC_Spearman_rho = nan(nClinical,1);
    PD_Spearman_rho = nan(nClinical,1);

    AdjustedClinical_RawSlope = nan(nClinical,1);
    AdjustedClinical_HC3_SE = nan(nClinical,1);

    AdjustedClinical_RawCI95_Low = nan(nClinical,1);
    AdjustedClinical_RawCI95_High = nan(nClinical,1);

    AdjustedClinical_StdBeta = nan(nClinical,1);
    AdjustedClinical_StdCI95_Low = nan(nClinical,1);
    AdjustedClinical_StdCI95_High = nan(nClinical,1);

    AdjustedClinical_T = nan(nClinical,1);
    AdjustedClinical_DF = nan(nClinical,1);
    AdjustedClinical_P = nan(nClinical,1);

    AdjustedClinical_VIF = nan(nClinical,1);
    AdjustedGroup_VIF = nan(nClinical,1);

    Interaction_StdBeta = nan(nClinical,1);
    Interaction_StdCI95_Low = nan(nClinical,1);
    Interaction_StdCI95_High = nan(nClinical,1);

    Interaction_T = nan(nClinical,1);
    Interaction_DF = nan(nClinical,1);
    Interaction_P = nan(nClinical,1);

    InteractionClinical_VIF = nan(nClinical,1);
    InteractionGroup_VIF = nan(nClinical,1);
    InteractionTerm_VIF = nan(nClinical,1);

    HC_InteractionSlope_Std = nan(nClinical,1);
    PD_InteractionSlope_Std = nan(nClinical,1);

    for c = 1:nClinical

        x = double(clinicalMatrix(:,c));

        good = ...
            isfinite(yAll) & ...
            isfinite(x) & ...
            isfinite(ageAll) & ...
            isfinite(sexAll) & ...
            isfinite(groupIntAll);

        y = double(yAll(good));
        x = double(x(good));
        age = double(ageAll(good));
        sex = double(sexAll(good));
        groupInt = double(groupIntAll(good));
        groups = string(groupsAll(good));

        if numel(y) < 30
            error('Too few observations for %s: N=%d', ...
                clinicalLabels(c),numel(y));
        end

        if std(x,0)==0 || std(y,0)==0
            error('Zero variance for %s.',clinicalLabels(c));
        end

        N(c) = numel(y);
        N_HC(c) = sum(groups=="HC");
        N_PD(c) = sum(groups=="PD");

        xHC = x(groups=="HC");
        xPD = x(groups=="PD");

        yHC = y(groups=="HC");
        yPD = y(groups=="PD");

        HC_Min(c) = min(xHC);
        HC_Max(c) = max(xHC);

        PD_Min(c) = min(xPD);
        PD_Max(c) = max(xPD);

        OverlapLow(c) = max(HC_Min(c),PD_Min(c));
        OverlapHigh(c) = min(HC_Max(c),PD_Max(c));

        OverlapWidth(c) = ...
            OverlapHigh(c)-OverlapLow(c);

        CommonSupport(c) = ...
            OverlapHigh(c) >= OverlapLow(c);

        % -------------------------------------------------------------
        % A. Pooled Spearman across HC + PD
        % -------------------------------------------------------------

        PooledSpearman_rho(c) = ...
            spearman_rho(x,y);

        PooledSpearman_Pperm(c) = ...
            spearman_permutation_p( ...
                x,y,nPerm,PooledSpearman_rho(c));

        % -------------------------------------------------------------
        % D. Within-group descriptive Spearman
        % -------------------------------------------------------------

        HC_Spearman_rho(c) = ...
            spearman_rho(xHC,yHC);

        PD_Spearman_rho(c) = ...
            spearman_rho(xPD,yPD);

        % -------------------------------------------------------------
        % B. Group-adjusted HC3 regression:
        %    y ~ x + Group + Age_c + Sex
        % -------------------------------------------------------------

        ageC = age-mean(age);

        Xadj = [ ...
            ones(numel(y),1), ...
            x, ...
            groupInt, ...
            ageC, ...
            sex];

        stAdj = hc3_predictor_effect( ...
            y,Xadj,2);

        AdjustedClinical_RawSlope(c) = ...
            stAdj.beta;

        AdjustedClinical_HC3_SE(c) = ...
            stAdj.se;

        AdjustedClinical_RawCI95_Low(c) = ...
            stAdj.ciLow;

        AdjustedClinical_RawCI95_High(c) = ...
            stAdj.ciHigh;

        AdjustedClinical_T(c) = ...
            stAdj.t;

        AdjustedClinical_DF(c) = ...
            stAdj.df;

        AdjustedClinical_P(c) = ...
            stAdj.p;

        scaleFactor = std(x,0)/std(y,0);

        AdjustedClinical_StdBeta(c) = ...
            stAdj.beta*scaleFactor;

        AdjustedClinical_StdCI95_Low(c) = ...
            stAdj.ciLow*scaleFactor;

        AdjustedClinical_StdCI95_High(c) = ...
            stAdj.ciHigh*scaleFactor;

        vifAdj = predictor_vifs(Xadj(:,2:end));

        AdjustedClinical_VIF(c) = vifAdj(1);
        AdjustedGroup_VIF(c) = vifAdj(2);

        % -------------------------------------------------------------
        % C. Group x Score interaction
        %
        % Use globally standardized x and y for directly standardized
        % interaction coefficients.
        % -------------------------------------------------------------

        xZ = (x-mean(x))/std(x,0);
        yZ = (y-mean(y))/std(y,0);

        xGroup = xZ.*groupInt;

        Xint = [ ...
            ones(numel(y),1), ...
            xZ, ...
            groupInt, ...
            xGroup, ...
            ageC, ...
            sex];

        stInt = hc3_predictor_effect( ...
            yZ,Xint,4);

        Interaction_StdBeta(c) = ...
            stInt.beta;

        Interaction_StdCI95_Low(c) = ...
            stInt.ciLow;

        Interaction_StdCI95_High(c) = ...
            stInt.ciHigh;

        Interaction_T(c) = ...
            stInt.t;

        Interaction_DF(c) = ...
            stInt.df;

        Interaction_P(c) = ...
            stInt.p;

        betaIntAll = Xint\yZ;

        HC_InteractionSlope_Std(c) = ...
            betaIntAll(2);

        PD_InteractionSlope_Std(c) = ...
            betaIntAll(2)+betaIntAll(4);

        vifInt = predictor_vifs(Xint(:,2:end));

        InteractionClinical_VIF(c) = vifInt(1);
        InteractionGroup_VIF(c) = vifInt(2);
        InteractionTerm_VIF(c) = vifInt(3);

    end

    % Separate correction families.
    PooledSpearman_Q_BH_FDR = ...
        bh_adjust(PooledSpearman_Pperm);

    AdjustedClinical_Q_BH_FDR = ...
        bh_adjust(AdjustedClinical_P);

    Interaction_Q_BH_FDR = ...
        bh_adjust(Interaction_P);

    PooledSpearman_FDR_Significant = ...
        PooledSpearman_Q_BH_FDR < fdrAlpha;

    AdjustedClinical_FDR_Significant = ...
        AdjustedClinical_Q_BH_FDR < fdrAlpha;

    Interaction_FDR_Significant = ...
        Interaction_Q_BH_FDR < fdrAlpha;

    T = table( ...
        Outcome, ...
        ClinicalVariable, ...
        VariableName, ...
        Unit, ...
        N, ...
        N_HC, ...
        N_PD, ...
        HC_Min, ...
        HC_Max, ...
        PD_Min, ...
        PD_Max, ...
        OverlapLow, ...
        OverlapHigh, ...
        OverlapWidth, ...
        CommonSupport, ...
        PooledSpearman_rho, ...
        PooledSpearman_Pperm, ...
        PooledSpearman_Q_BH_FDR, ...
        PooledSpearman_FDR_Significant, ...
        HC_Spearman_rho, ...
        PD_Spearman_rho, ...
        AdjustedClinical_RawSlope, ...
        AdjustedClinical_HC3_SE, ...
        AdjustedClinical_RawCI95_Low, ...
        AdjustedClinical_RawCI95_High, ...
        AdjustedClinical_StdBeta, ...
        AdjustedClinical_StdCI95_Low, ...
        AdjustedClinical_StdCI95_High, ...
        AdjustedClinical_T, ...
        AdjustedClinical_DF, ...
        AdjustedClinical_P, ...
        AdjustedClinical_Q_BH_FDR, ...
        AdjustedClinical_FDR_Significant, ...
        AdjustedClinical_VIF, ...
        AdjustedGroup_VIF, ...
        Interaction_StdBeta, ...
        Interaction_StdCI95_Low, ...
        Interaction_StdCI95_High, ...
        Interaction_T, ...
        Interaction_DF, ...
        Interaction_P, ...
        Interaction_Q_BH_FDR, ...
        Interaction_FDR_Significant, ...
        InteractionClinical_VIF, ...
        InteractionGroup_VIF, ...
        InteractionTerm_VIF, ...
        HC_InteractionSlope_Std, ...
        PD_InteractionSlope_Std);

end

%% =========================================================================
%% LOCAL FUNCTION: SCORE DISTRIBUTION TABLE
%% =========================================================================

function T = build_score_distribution_table( ...
    clinicalMatrix,groups,labels,units)

    nC = numel(labels);

    ClinicalVariable = labels;
    Unit = units;

    N_HC = nan(nC,1);
    HC_Mean = nan(nC,1);
    HC_SD = nan(nC,1);
    HC_Median = nan(nC,1);
    HC_Min = nan(nC,1);
    HC_Max = nan(nC,1);

    N_PD = nan(nC,1);
    PD_Mean = nan(nC,1);
    PD_SD = nan(nC,1);
    PD_Median = nan(nC,1);
    PD_Min = nan(nC,1);
    PD_Max = nan(nC,1);

    OverlapLow = nan(nC,1);
    OverlapHigh = nan(nC,1);
    OverlapWidth = nan(nC,1);
    CommonSupport = false(nC,1);

    for c = 1:nC

        x = double(clinicalMatrix(:,c));

        xHC = x(groups=="HC" & isfinite(x));
        xPD = x(groups=="PD" & isfinite(x));

        N_HC(c) = numel(xHC);
        HC_Mean(c) = mean(xHC);
        HC_SD(c) = std(xHC,0);
        HC_Median(c) = median(xHC);
        HC_Min(c) = min(xHC);
        HC_Max(c) = max(xHC);

        N_PD(c) = numel(xPD);
        PD_Mean(c) = mean(xPD);
        PD_SD(c) = std(xPD,0);
        PD_Median(c) = median(xPD);
        PD_Min(c) = min(xPD);
        PD_Max(c) = max(xPD);

        OverlapLow(c) = max(HC_Min(c),PD_Min(c));
        OverlapHigh(c) = min(HC_Max(c),PD_Max(c));
        OverlapWidth(c) = OverlapHigh(c)-OverlapLow(c);

        CommonSupport(c) = ...
            OverlapHigh(c) >= OverlapLow(c);

    end

    T = table( ...
        ClinicalVariable, ...
        Unit, ...
        N_HC, ...
        HC_Mean, ...
        HC_SD, ...
        HC_Median, ...
        HC_Min, ...
        HC_Max, ...
        N_PD, ...
        PD_Mean, ...
        PD_SD, ...
        PD_Median, ...
        PD_Min, ...
        PD_Max, ...
        OverlapLow, ...
        OverlapHigh, ...
        OverlapWidth, ...
        CommonSupport);

end

%% =========================================================================
%% LOCAL FUNCTION: SPEARMAN PERMUTATION P VALUE
%% =========================================================================

function p = spearman_permutation_p(x,y,nPerm,rhoObs)

    x = double(x(:));
    y = double(y(:));

    good = isfinite(x) & isfinite(y);

    x = x(good);
    y = y(good);

    n = numel(x);

    if n<3
        p = NaN;
        return;
    end

    rx = average_ranks(x);
    ry = average_ranks(y);

    rhoObsAbs = abs(rhoObs);

    exceed = 0;

    for b = 1:nPerm

        rp = ry(randperm(n));

        rhoPerm = pearson_corr(rx,rp);

        if abs(rhoPerm) >= rhoObsAbs
            exceed = exceed+1;
        end

    end

    p = (1+exceed)/(nPerm+1);

end

%% =========================================================================
%% LOCAL FUNCTION: HC3 ROBUST COEFFICIENT INFERENCE
%% =========================================================================

function st = hc3_predictor_effect(y,X,predictorColumn)

    y = double(y(:));
    X = double(X);

    N = numel(y);
    p = size(X,2);

    if size(X,1) ~= N
        error('X/y dimension mismatch.');
    end

    if predictorColumn < 1 || predictorColumn > p
        error('Invalid predictor column.');
    end

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

    varB = Vhc3(predictorColumn,predictorColumn);

    if ~isfinite(varB) || varB<=0
        error('Invalid HC3 predictor variance.');
    end

    se = sqrt(varB);

    tStat = beta(predictorColumn)/se;

    df = N-p;

    if df<=0
        error('Nonpositive residual degrees of freedom.');
    end

    pValue = t_two_sided_p(tStat,df);

    tcrit = student_t_critical_975(df);

    ciLow = beta(predictorColumn)-tcrit*se;
    ciHigh = beta(predictorColumn)+tcrit*se;

    st = struct();

    st.beta = beta(predictorColumn);
    st.se = se;
    st.t = tStat;
    st.df = df;
    st.p = pValue;
    st.ciLow = ciLow;
    st.ciHigh = ciHigh;

end

%% =========================================================================
%% LOCAL FUNCTION: PREDICTOR VIFS
%% =========================================================================

function vifs = predictor_vifs(Z)

    Z = double(Z);

    [n,p] = size(Z);

    vifs = nan(p,1);

    for j = 1:p

        y = Z(:,j);

        others = setdiff(1:p,j);

        X = [ones(n,1),Z(:,others)];

        b = X\y;

        resid = y-X*b;

        sst = sum((y-mean(y)).^2);
        sse = sum(resid.^2);

        if sst<=0

            vifs(j) = Inf;

        else

            R2 = 1-sse/sst;

            R2 = max(0,min(1,R2));

            if R2>=1
                vifs(j) = Inf;
            else
                vifs(j) = 1/(1-R2);
            end

        end

    end

end

%% =========================================================================
%% LOCAL FUNCTION: BH-FDR
%% =========================================================================

function q = bh_adjust(p)

    p = double(p(:));

    m = numel(p);

    if any(~isfinite(p)) || ...
       any(p<0) || ...
       any(p>1)

        error('Invalid P values supplied to BH-FDR.');

    end

    [pSort,ord] = sort(p,'ascend');

    raw = pSort .* m ./ (1:m)';

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
%% LOCAL FUNCTION: SPEARMAN RHO
%% =========================================================================

function rho = spearman_rho(x,y)

    x = double(x(:));
    y = double(y(:));

    good = isfinite(x) & isfinite(y);

    x = x(good);
    y = y(good);

    if numel(x)<3
        rho = NaN;
        return;
    end

    rx = average_ranks(x);
    ry = average_ranks(y);

    rho = pearson_corr(rx,ry);

end

%% =========================================================================
%% LOCAL FUNCTION: AVERAGE RANKS WITH TIES
%% =========================================================================

function r = average_ranks(x)

    x = double(x(:));

    n = numel(x);

    [xs,ord] = sort(x,'ascend');

    rSorted = nan(n,1);

    i = 1;

    while i<=n

        j = i;

        while j<n && xs(j+1)==xs(i)
            j = j+1;
        end

        meanRank = mean(i:j);

        rSorted(i:j) = meanRank;

        i = j+1;

    end

    r = nan(n,1);
    r(ord) = rSorted;

end

%% =========================================================================
%% LOCAL FUNCTION: PEARSON CORRELATION
%% =========================================================================

function r = pearson_corr(x,y)

    x = double(x(:));
    y = double(y(:));

    good = isfinite(x) & isfinite(y);

    x = x(good);
    y = y(good);

    if numel(x)<3 || std(x,0)==0 || std(y,0)==0
        r = NaN;
        return;
    end

    x = x-mean(x);
    y = y-mean(y);

    r = sum(x.*y) / sqrt(sum(x.^2)*sum(y.^2));

    r = max(-1,min(1,r));

end

%% =========================================================================
%% LOCAL FUNCTION: TWO-SIDED STUDENT-t P
%% =========================================================================

function p = t_two_sided_p(t,df)

    if ~isfinite(t) || df<=0
        error('Invalid t statistic or degrees of freedom.');
    end

    x = df/(df+t^2);

    p = betainc(x,df/2,0.5);

    p = max(0,min(1,p));

end

%% =========================================================================
%% LOCAL FUNCTION: 95% t CRITICAL VALUE
%% =========================================================================

function tcrit = student_t_critical_975(df)

    x = betaincinv(0.05,df/2,0.5);

    tcrit = sqrt(df*(1-x)/x);

end

%% =========================================================================
%% LOCAL FUNCTION: TABLE COLUMN -> DOUBLE
%% =========================================================================

function x = column_to_double(col)

    if isnumeric(col)

        x = double(col(:));
        return;

    end

    s = strtrim(string(col(:)));

    missingTokens = ...
        ismissing(s) | ...
        lower(s)=="n/a" | ...
        lower(s)=="na" | ...
        lower(s)=="nan" | ...
        s=="";

    x = str2double(s);

    x(missingTokens) = NaN;

end

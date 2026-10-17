%% ========================================================================
% PAPER 3 — STEP 10S
% FROZEN-PROTOCOL DECOMPOSITION OF SIGNIFICANT ALPHA
% GROUP × NETWORKTYPE INTERACTION
%
% PLANNED CONTRASTS ONLY:
%
% C1:
%   Inter - mean(IntraLeft, IntraRight)
%
% C2:
%   IntraLeft - IntraRight
%
% These contrasts test the NETWORKTYPE dependence of the PD-HC
% rest-to-walking reconfiguration difference.
%
% Holm correction is applied across THESE TWO TESTS ONLY.
%
% No additional pairwise comparisons.
% No alpha edge-level follow-up.
% ========================================================================

clear;
clc;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% ========================================================================
% 1. PATHS
% ========================================================================

baseDir = ...
    fullfile(cfg.datasetRoot,'PAPER3_SECONDARY_HEMISPHERIC_PROTOCOL_ALIGNED');

inputFile = fullfile(baseDir, ...
    'Paper3_SECONDARY_PROTOCOL_ALIGNED_subject_networktype_medians_long.csv');

outTable = fullfile(baseDir, ...
    'Paper3_STEP10S_ALPHA_FrozenProtocol_Decomposition_Holm2.csv');

outMAT = fullfile(baseDir, ...
    'Paper3_STEP10S_ALPHA_FrozenProtocol_Decomposition_RESULTS.mat');

outLog = fullfile(baseDir, ...
    'Paper3_STEP10S_ALPHA_FrozenProtocol_Decomposition_log.txt');


fprintf('\n============================================================\n');
fprintf('PAPER 3 — STEP 10S\n');
fprintf('FROZEN-PROTOCOL ALPHA INTERACTION DECOMPOSITION\n');
fprintf('============================================================\n\n');


%% ========================================================================
% 2. LOAD DATA
% ========================================================================

T = readtable(inputFile, 'VariableNamingRule','preserve');

fprintf('Loaded:\n%s\n\n', inputFile);


%% ========================================================================
% 3. STANDARDIZE VARIABLES
% ========================================================================

D = table;

D.Subject = categorical(string(T.Subject));

D.Group = categorical(upper(strtrim(string(T.Group))));

D.Band = lower(strtrim(string(T.Band)));

D.NetworkType = categorical(string(T.NetworkType));

D.Age_c = double(T.Age_c);

D.Sex = categorical(string(T.Sex));

D.Y = double(T.DeltaFC_Median);


%% ========================================================================
% 4. KEEP ALPHA ONLY
% ========================================================================

A = D(D.Band == "alpha", :);

A = rmmissing(A);

A.Subject     = removecats(A.Subject);
A.Group       = removecats(A.Group);
A.NetworkType = removecats(A.NetworkType);
A.Sex         = removecats(A.Sex);


% Force the same reference categories used in Step 10R
A.Group = reordercats(A.Group, {'HC','PD'});

A.NetworkType = reordercats(A.NetworkType, ...
    {'IntraLeft','IntraRight','Inter'});


%% ========================================================================
% 5. COHORT CHECK
% ========================================================================

nTotal = numel(unique(A.Subject));
nHC    = numel(unique(A.Subject(A.Group == 'HC')));
nPD    = numel(unique(A.Subject(A.Group == 'PD')));

fprintf('Cohort:\n');
fprintf('N total = %d\n', nTotal);
fprintf('N HC    = %d\n', nHC);
fprintf('N PD    = %d\n\n', nPD);

if nTotal ~= 130 || nHC ~= 23 || nPD ~= 107
    warning('Cohort differs from protocol-aligned Step 10R.');
end


%% ========================================================================
% 6. FIT IDENTICAL ALPHA LME
% ========================================================================

formula = ...
    'Y ~ Group*NetworkType + Age_c + Sex + (1|Subject)';

lme = fitlme(A, formula, ...
    'FitMethod','REML');


%% ========================================================================
% 7. VERIFY ORIGINAL INTERACTION
% ========================================================================

anovaTbl = anova(lme, ...
    'DFMethod','satterthwaite');

disp(anovaTbl);

terms = string(anovaTbl.Term);

idxOmnibus = ...
    contains(terms,'Group:NetworkType') | ...
    contains(terms,'NetworkType:Group');

Fomni  = anovaTbl.FStat(idxOmnibus);
DF1omni = anovaTbl.DF1(idxOmnibus);
DF2omni = anovaTbl.DF2(idxOmnibus);
Pomni   = anovaTbl.pValue(idxOmnibus);

fprintf('\nOriginal alpha Group x NetworkType interaction:\n');
fprintf('F(%g,%g) = %.6f, p = %.8f\n\n', ...
    DF1omni, DF2omni, Fomni, Pomni);


% Verify against the current Step 10R output, not hard-coded historical values.
% This avoids false failures when reproducibility runs differ slightly
% from the original manuscript reference statistics.
refFile = fullfile(baseDir, ...
    'Paper3_SECONDARY_PROTOCOL_ALIGNED_GroupxNetworkType_LME_Holm.csv');
if ~isfile(refFile)
    error('Missing Step 10R reference table: %s', refFile);
end
ref = readtable(refFile, 'VariableNamingRule','preserve');
idxAlpha = find(strcmpi(strtrim(string(ref.Band)), 'alpha'));
if numel(idxAlpha) ~= 1
    error('Expected exactly one alpha row in Step 10R results.');
end
refF = ref.F_GroupxNetworkType(idxAlpha);
refP = ref.P_Raw(idxAlpha);
refDF1 = ref.DF1(idxAlpha);
refDF2 = ref.DF2(idxAlpha);
fprintf('Step 10R reference: F(%g,%g) = %.6f, p = %.8f\n', ...
    refDF1, refDF2, refF, refP);
if abs(Fomni - refF) > 1e-4 || ...
        abs(Pomni - refP) > 1e-5 || ...
        abs(DF1omni - refDF1) > 1e-4 || ...
        abs(DF2omni - refDF2) > 1e-3
    error(['Alpha model does not reproduce current Step 10R output. ' ...
           'Stop before interpreting decomposition.']);
end
fprintf('STEP 10R REPRODUCTION CHECK PASSED.\n\n');


%% ========================================================================
% 8. FIXED EFFECT COEFFICIENTS
% ========================================================================

[beta,betaNames] = fixedEffects(lme);

coefNames = string(betaNames.Name);

fprintf('Fixed-effect coefficient order:\n');

for k = 1:numel(coefNames)
    fprintf('%2d : %s\n', k, coefNames(k));
end

fprintf('\n');


%% ========================================================================
% 9. IDENTIFY INTERACTION COEFFICIENTS
% ========================================================================

idxRight = find( ...
    contains(coefNames,'Group_PD') & ...
    contains(coefNames,'NetworkType_IntraRight') & ...
    contains(coefNames,':'));

idxInter = find( ...
    contains(coefNames,'Group_PD') & ...
    contains(coefNames,'NetworkType_Inter') & ...
    contains(coefNames,':'));


if numel(idxRight) ~= 1
    error('Could not uniquely identify PD x IntraRight coefficient.');
end

if numel(idxInter) ~= 1
    error('Could not uniquely identify PD x Inter coefficient.');
end


%% ========================================================================
% 10. DEFINE THE TWO FROZEN-PROTOCOL CONTRASTS
% ========================================================================
%
% Let:
%
% G_L = PD-HC effect for IntraLeft
% G_R = PD-HC effect for IntraRight
% G_I = PD-HC effect for Inter
%
%
% Model coding gives:
%
% bRight = G_R - G_L
% bInter = G_I - G_L
%
%
% --------------------------------------------------------------
% C1:
%
% G_I - (G_L + G_R)/2
%
% = bInter - 0.5*bRight
%
%
% --------------------------------------------------------------
% C2:
%
% G_L - G_R
%
% = -bRight
%
% ========================================================================

p = numel(beta);

H = zeros(2,p);


% C1:
% Inter - average(IntraLeft, IntraRight)
H(1,idxInter) =  1;
H(1,idxRight) = -0.5;


% C2:
% IntraLeft - IntraRight
H(2,idxRight) = -1;


contrastNames = [
    "Inter - mean(IntraLeft, IntraRight)"
    "IntraLeft - IntraRight"
    ];


%% ========================================================================
% 11. TEST CONTRASTS
% ========================================================================

Vbeta = lme.CoefficientCovariance;

nTests = 2;

Estimate = nan(nTests,1);
SE       = nan(nTests,1);
tStat    = nan(nTests,1);
DF       = nan(nTests,1);
FStat    = nan(nTests,1);
P_Raw    = nan(nTests,1);
CI95_Low = nan(nTests,1);
CI95_High = nan(nTests,1);


for j = 1:nTests

    h = H(j,:);

    % Contrast estimate
    Estimate(j) = h * beta;

    % Contrast SE
    varContrast = h * Vbeta * h';
    SE(j) = sqrt(varContrast);

    tStat(j) = Estimate(j) / SE(j);


    % Satterthwaite inference
    [P_Raw(j), FStat(j), ~, DF(j)] = ...
        coefTest(lme, h, 0, ...
        'DFMethod','satterthwaite');


    % 95% confidence interval
    tcrit = tinv(0.975, DF(j));

    CI95_Low(j) = ...
        Estimate(j) - tcrit*SE(j);

    CI95_High(j) = ...
        Estimate(j) + tcrit*SE(j);

end


%% ========================================================================
% 12. HOLM CORRECTION ACROSS EXACTLY TWO TESTS
% ========================================================================

[Psorted, order] = sort(P_Raw);

m = 2;

PadjSorted = nan(m,1);

for k = 1:m
    PadjSorted(k) = ...
        (m-k+1) * Psorted(k);
end

% Enforce Holm monotonicity
PadjSorted = cummax(PadjSorted);

PadjSorted(PadjSorted > 1) = 1;

P_Holm_2 = nan(m,1);

P_Holm_2(order) = PadjSorted;


Significant_Holm_0p05 = ...
    P_Holm_2 < 0.05;


%% ========================================================================
% 13. RESULT TABLE
% ========================================================================

Results = table( ...
    contrastNames, ...
    Estimate, ...
    SE, ...
    tStat, ...
    DF, ...
    FStat, ...
    CI95_Low, ...
    CI95_High, ...
    P_Raw, ...
    P_Holm_2, ...
    Significant_Holm_0p05, ...
    'VariableNames', { ...
        'PlannedContrast', ...
        'Estimate', ...
        'SE', ...
        'tStat', ...
        'DF', ...
        'FStat', ...
        'CI95_Low', ...
        'CI95_High', ...
        'P_Raw', ...
        'P_Holm_2', ...
        'Significant_Holm_0p05'});


fprintf('\n============================================================\n');
fprintf('FROZEN-PROTOCOL ALPHA DECOMPOSITION\n');
fprintf('HOLM CORRECTION ACROSS TWO PLANNED TESTS\n');
fprintf('============================================================\n\n');

disp(Results);


%% ========================================================================
% 14. SAVE
% ========================================================================

writetable(Results, outTable);

save(outMAT, ...
    'Results', ...
    'lme', ...
    'anovaTbl', ...
    'H', ...
    'beta', ...
    'betaNames', ...
    'A');


%% ========================================================================
% 15. LOG
% ========================================================================

fid = fopen(outLog,'w');

fprintf(fid, ...
    'PAPER 3 STEP 10S — FROZEN PROTOCOL ALPHA DECOMPOSITION\n\n');

fprintf(fid, ...
    'Primary omnibus interaction reproduced:\n');

fprintf(fid, ...
    'F(%g,%g) = %.8f, p = %.10f\n\n', ...
    DF1omni, DF2omni, Fomni, Pomni);


fprintf(fid, ...
    'ONLY TWO PRE-SPECIFIED DECOMPOSITION CONTRASTS TESTED:\n\n');

fprintf(fid, ...
    '1. Inter - mean(IntraLeft, IntraRight)\n');

fprintf(fid, ...
    '2. IntraLeft - IntraRight\n\n');

fprintf(fid, ...
    'Both contrasts operate on the PD-HC reconfiguration effect.\n\n');

fprintf(fid, ...
    'Multiplicity correction: Holm across exactly 2 contrasts.\n\n');

fprintf(fid, ...
    'No additional pairwise NetworkType comparisons performed.\n');

fprintf(fid, ...
    'No alpha edge-level follow-up performed.\n\n');

fprintf(fid, ...
    'Results:\n');

for j = 1:height(Results)

    fprintf(fid, ...
        '\n%s\n', Results.PlannedContrast(j));

    fprintf(fid, ...
        'Estimate = %.10f\n', Results.Estimate(j));

    fprintf(fid, ...
        'SE = %.10f\n', Results.SE(j));

    fprintf(fid, ...
        '95%% CI = [%.10f, %.10f]\n', ...
        Results.CI95_Low(j), ...
        Results.CI95_High(j));

    fprintf(fid, ...
        'p_raw = %.10f\n', ...
        Results.P_Raw(j));

    fprintf(fid, ...
        'p_Holm_2 = %.10f\n', ...
        Results.P_Holm_2(j));

end

fclose(fid);


fprintf('\n============================================================\n');
fprintf('STEP 10S COMPLETE\n');
fprintf('============================================================\n');

fprintf('\nSaved:\n');
fprintf('Results : %s\n', outTable);
fprintf('MAT     : %s\n', outMAT);
fprintf('Log     : %s\n', outLog);

fprintf('\nNo additional alpha decomposition should be performed.\n');
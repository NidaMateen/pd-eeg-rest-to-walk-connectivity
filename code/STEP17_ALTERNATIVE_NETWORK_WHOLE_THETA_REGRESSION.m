%% =========================================================================
% PAPER 3 — STEP 17
% ALTERNATIVE-NETWORK WHOLE-THETA REGRESSION
%
% PURPOSE
% -------
% Test whether the PD-HC difference in Rest-to-Walk theta reconfiguration
% persists when the network representation is changed from the primary
% regional median-aggregation network to the frozen one-electrode-per-node
% alternative network.
%
% MAIN sensitivity:
%   Matched cohort N = 117 = 22 HC + 95 PD
%   SAME subjects as original primary complete-network analysis.
%
% SUPPORTIVE sensitivity:
%   Expanded alternative cohort N = 124 = 23 HC + 101 PD
%
% Outcome:
%   Mean theta Delta db-wPLI across all 55 alternative-network edges.
%
% Model:
%   WholeThetaDelta ~ Group + Age_c + Sex
%
% Group:
%   HC = 0
%   PD = 1
%
% Sex:
%   F = 0
%   M = 1
%
% Inference:
%   OLS coefficient + HC3 heteroscedasticity-consistent covariance.
%
% PRIMARY NETWORK REFERENCE:
%   Original whole-network theta adjusted PD-HC beta = -0.042697
%   95% CI = [-0.066857, -0.018536]
%
% IMPORTANT:
%   The matched N=117 analysis is the principal network-construction
%   robustness analysis because only network representation changes.
%% =========================================================================

clear;
clc;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% =========================================================================
% 1. PATHS
%% =========================================================================

datasetRoot = ...
    cfg.datasetRoot;

altRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_ALTERNATIVE_NETWORK');

matchedFile = fullfile( ...
    altRoot, ...
    'Paper3_ALTERNATIVE11_MATCHED_N117.mat');

expandedFile = fullfile( ...
    altRoot, ...
    'Paper3_ALTERNATIVE11_EXPANDED_N124.mat');

participantsFile = fullfile( ...
    datasetRoot, ...
    'participants.tsv');

if ~exist(matchedFile,'file')
    error('Matched N=117 file not found:\n%s',matchedFile);
end

if ~exist(expandedFile,'file')
    error('Expanded N=124 file not found:\n%s',expandedFile);
end

if ~exist(participantsFile,'file')
    error('participants.tsv not found:\n%s',participantsFile);
end

%% =========================================================================
% 2. LOAD MATCHED N=117 ALTERNATIVE DATASET
%% =========================================================================

M = load(matchedFile);

requiredMatched = { ...
    'subjectsMatched', ...
    'groupsMatched', ...
    'deltaMatched', ...
    'thetaIdx', ...
    'EdgeTable', ...
    'MappingTable'};

for k = 1:numel(requiredMatched)

    if ~isfield(M,requiredMatched{k})
        error('Matched file missing field: %s',requiredMatched{k});
    end

end

subjectsMatched = string(M.subjectsMatched(:));
groupsMatched   = string(M.groupsMatched(:));

deltaMatched = double(M.deltaMatched);

thetaIdxMatched = double(M.thetaIdx);

if size(deltaMatched,1) ~= 117
    error('Expected matched N=117.');
end

if size(deltaMatched,2) ~= 55
    error('Expected 55 alternative-network edges.');
end

if sum(groupsMatched=="HC") ~= 22 || ...
   sum(groupsMatched=="PD") ~= 95

    error('Expected matched cohort 22 HC + 95 PD.');

end

%% =========================================================================
% 3. LOAD EXPANDED N=124 ALTERNATIVE DATASET
%% =========================================================================

E = load(expandedFile);

requiredExpanded = { ...
    'subjectsExpanded', ...
    'groupsExpanded', ...
    'deltaExpanded', ...
    'thetaIdx'};

for k = 1:numel(requiredExpanded)

    if ~isfield(E,requiredExpanded{k})
        error('Expanded file missing field: %s',requiredExpanded{k});
    end

end

subjectsExpanded = string(E.subjectsExpanded(:));
groupsExpanded   = string(E.groupsExpanded(:));

deltaExpanded = double(E.deltaExpanded);

thetaIdxExpanded = double(E.thetaIdx);

if size(deltaExpanded,1) ~= 124
    error('Expected expanded N=124.');
end

if size(deltaExpanded,2) ~= 55
    error('Expected 55 alternative-network edges.');
end

if sum(groupsExpanded=="HC") ~= 23 || ...
   sum(groupsExpanded=="PD") ~= 101

    error('Expected expanded cohort 23 HC + 101 PD.');

end

if thetaIdxMatched ~= thetaIdxExpanded
    error('Theta-band index differs between datasets.');
end

thetaIdx = thetaIdxMatched;

%% =========================================================================
% 4. VALIDATE THETA EDGE COMPLETENESS
%% =========================================================================

thetaMatchedEdges = ...
    deltaMatched(:,:,thetaIdx);

thetaExpandedEdges = ...
    deltaExpanded(:,:,thetaIdx);

if any(~isfinite(thetaMatchedEdges(:)))
    error('Matched theta data contains missing/nonfinite values.');
end

if any(~isfinite(thetaExpandedEdges(:)))
    error('Expanded theta data contains missing/nonfinite values.');
end

fprintf('\n============================================================\n');
fprintf('STEP 17 — ALTERNATIVE WHOLE-NETWORK THETA REGRESSION\n');
fprintf('============================================================\n');

fprintf('Matched cohort  : %d = %d HC + %d PD\n', ...
    numel(subjectsMatched), ...
    sum(groupsMatched=="HC"), ...
    sum(groupsMatched=="PD"));

fprintf('Expanded cohort : %d = %d HC + %d PD\n', ...
    numel(subjectsExpanded), ...
    sum(groupsExpanded=="HC"), ...
    sum(groupsExpanded=="PD"));

fprintf('Edges per subject = %d\n', ...
    size(thetaMatchedEdges,2));

%% =========================================================================
% 5. CALCULATE WHOLE-NETWORK THETA DELTA
%% =========================================================================

wholeThetaMatched = ...
    mean(thetaMatchedEdges,2);

wholeThetaExpanded = ...
    mean(thetaExpandedEdges,2);

wholeThetaMatched = ...
    wholeThetaMatched(:);

wholeThetaExpanded = ...
    wholeThetaExpanded(:);

if any(~isfinite(wholeThetaMatched))
    error('Matched whole-theta outcome contains nonfinite values.');
end

if any(~isfinite(wholeThetaExpanded))
    error('Expanded whole-theta outcome contains nonfinite values.');
end

%% =========================================================================
% 6. LOAD PARTICIPANT DEMOGRAPHICS
%% =========================================================================

P = readtable( ...
    participantsFile, ...
    'FileType','text', ...
    'Delimiter','\t', ...
    'TextType','string');

if ~ismember('participant_id',P.Properties.VariableNames)
    error('participants.tsv lacks participant_id.');
end

if ~ismember('age',P.Properties.VariableNames)
    error('participants.tsv lacks age.');
end

if ~ismember('sex',P.Properties.VariableNames)
    error('participants.tsv lacks sex.');
end

%% =========================================================================
% 7. ALIGN MATCHED COHORT DEMOGRAPHICS
%% =========================================================================

[tfM,locM] = ...
    ismember( ...
        subjectsMatched, ...
        string(P.participant_id));

if ~all(tfM)

    fprintf('\nMissing matched subjects in participants.tsv:\n');
    disp(subjectsMatched(~tfM));

    error('Matched participant alignment failed.');

end

ageMatched = ...
    column_to_double(P.age(locM));

sexMatchedRaw = ...
    upper(strtrim(string(P.sex(locM))));

if any(~ismember(sexMatchedRaw,["F","M"]))
    error('Unexpected sex coding in matched cohort.');
end

GroupMatched = ...
    double(groupsMatched=="PD");

SexMatched = ...
    double(sexMatchedRaw=="M");

AgeMatched_c = ...
    ageMatched - mean(ageMatched);

%% =========================================================================
% 8. ALIGN EXPANDED COHORT DEMOGRAPHICS
%% =========================================================================

[tfE,locE] = ...
    ismember( ...
        subjectsExpanded, ...
        string(P.participant_id));

if ~all(tfE)

    fprintf('\nMissing expanded subjects in participants.tsv:\n');
    disp(subjectsExpanded(~tfE));

    error('Expanded participant alignment failed.');

end

ageExpanded = ...
    column_to_double(P.age(locE));

sexExpandedRaw = ...
    upper(strtrim(string(P.sex(locE))));

if any(~ismember(sexExpandedRaw,["F","M"]))
    error('Unexpected sex coding in expanded cohort.');
end

GroupExpanded = ...
    double(groupsExpanded=="PD");

SexExpanded = ...
    double(sexExpandedRaw=="M");

AgeExpanded_c = ...
    ageExpanded - mean(ageExpanded);

%% =========================================================================
% 9. MATCHED N=117 DESCRIPTIVE VALUES
%% =========================================================================

matchedHC = ...
    wholeThetaMatched(GroupMatched==0);

matchedPD = ...
    wholeThetaMatched(GroupMatched==1);

fprintf('\n============================================================\n');
fprintf('MATCHED N=117 — DESCRIPTIVE VALUES\n');
fprintf('============================================================\n');

fprintf('HC mean Delta = %.6f\n', ...
    mean(matchedHC));

fprintf('HC SD         = %.6f\n', ...
    std(matchedHC));

fprintf('PD mean Delta = %.6f\n', ...
    mean(matchedPD));

fprintf('PD SD         = %.6f\n', ...
    std(matchedPD));

fprintf('Raw PD-HC difference = %.6f\n', ...
    mean(matchedPD)-mean(matchedHC));

%% =========================================================================
% 10. MATCHED N=117 HC3 REGRESSION
%% =========================================================================

Xmatched = [ ...
    ones(numel(wholeThetaMatched),1), ...
    GroupMatched, ...
    AgeMatched_c, ...
    SexMatched];

Rmatched = ...
    hc3_regression( ...
        wholeThetaMatched, ...
        Xmatched);

matchedBeta = ...
    Rmatched.beta(2);

matchedSE = ...
    Rmatched.SE(2);

matchedCI_low = ...
    Rmatched.CI_low(2);

matchedCI_high = ...
    Rmatched.CI_high(2);

matchedP = ...
    Rmatched.p(2);

fprintf('\n============================================================\n');
fprintf('MATCHED N=117 — ALTERNATIVE WHOLE-THETA MODEL\n');
fprintf('Theta Delta ~ Group + Age_c + Sex\n');
fprintf('HC3 heteroscedasticity-consistent inference\n');
fprintf('============================================================\n');

fprintf('Adjusted PD-HC beta = %.6f\n', ...
    matchedBeta);

fprintf('HC3 SE              = %.6f\n', ...
    matchedSE);

fprintf('95%% CI               = [%.6f, %.6f]\n', ...
    matchedCI_low, ...
    matchedCI_high);

fprintf('p                    = %.8f\n', ...
    matchedP);

%% =========================================================================
% 11. COMPARE MATCHED ALTERNATIVE BETA WITH PRIMARY NETWORK
%% =========================================================================

primaryBeta = ...
    -0.042697;

primaryCI_low = ...
    -0.066857;

primaryCI_high = ...
    -0.018536;

betaDifference = ...
    matchedBeta - primaryBeta;

magnitudeRatio = ...
    abs(matchedBeta) / abs(primaryBeta);

percentMagnitudeDifference = ...
    100 * ...
    (abs(matchedBeta)-abs(primaryBeta)) / ...
    abs(primaryBeta);

sameDirection = ...
    sign(matchedBeta)==sign(primaryBeta);

fprintf('\n============================================================\n');
fprintf('MATCHED N=117 — COMPARISON WITH PRIMARY NETWORK\n');
fprintf('============================================================\n');

fprintf('Primary-network beta     = %.6f\n', ...
    primaryBeta);

fprintf('Primary-network 95%% CI   = [%.6f, %.6f]\n', ...
    primaryCI_low, ...
    primaryCI_high);

fprintf('Alternative-network beta = %.6f\n', ...
    matchedBeta);

fprintf('Alternative-network 95%% CI = [%.6f, %.6f]\n', ...
    matchedCI_low, ...
    matchedCI_high);

fprintf('Beta difference (Alt - Primary) = %.6f\n', ...
    betaDifference);

fprintf('|Alt beta| / |Primary beta| = %.3f\n', ...
    magnitudeRatio);

fprintf('Percent difference in magnitude = %.2f%%\n', ...
    percentMagnitudeDifference);

fprintf('Same effect direction = %d\n', ...
    sameDirection);

%% =========================================================================
% 12. EXPANDED N=124 DESCRIPTIVE VALUES
%% =========================================================================

expandedHC = ...
    wholeThetaExpanded(GroupExpanded==0);

expandedPD = ...
    wholeThetaExpanded(GroupExpanded==1);

fprintf('\n============================================================\n');
fprintf('EXPANDED N=124 — DESCRIPTIVE VALUES\n');
fprintf('============================================================\n');

fprintf('HC mean Delta = %.6f\n', ...
    mean(expandedHC));

fprintf('HC SD         = %.6f\n', ...
    std(expandedHC));

fprintf('PD mean Delta = %.6f\n', ...
    mean(expandedPD));

fprintf('PD SD         = %.6f\n', ...
    std(expandedPD));

fprintf('Raw PD-HC difference = %.6f\n', ...
    mean(expandedPD)-mean(expandedHC));

%% =========================================================================
% 13. EXPANDED N=124 HC3 REGRESSION
%% =========================================================================

Xexpanded = [ ...
    ones(numel(wholeThetaExpanded),1), ...
    GroupExpanded, ...
    AgeExpanded_c, ...
    SexExpanded];

Rexpanded = ...
    hc3_regression( ...
        wholeThetaExpanded, ...
        Xexpanded);

expandedBeta = ...
    Rexpanded.beta(2);

expandedSE = ...
    Rexpanded.SE(2);

expandedCI_low = ...
    Rexpanded.CI_low(2);

expandedCI_high = ...
    Rexpanded.CI_high(2);

expandedP = ...
    Rexpanded.p(2);

fprintf('\n============================================================\n');
fprintf('EXPANDED N=124 — ALTERNATIVE WHOLE-THETA MODEL\n');
fprintf('Theta Delta ~ Group + Age_c + Sex\n');
fprintf('HC3 heteroscedasticity-consistent inference\n');
fprintf('============================================================\n');

fprintf('Adjusted PD-HC beta = %.6f\n', ...
    expandedBeta);

fprintf('HC3 SE              = %.6f\n', ...
    expandedSE);

fprintf('95%% CI               = [%.6f, %.6f]\n', ...
    expandedCI_low, ...
    expandedCI_high);

fprintf('p                    = %.8f\n', ...
    expandedP);

%% =========================================================================
% 14. COMPARE MATCHED VERSUS EXPANDED ALTERNATIVE NETWORK
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('MATCHED VS EXPANDED ALTERNATIVE NETWORK\n');
fprintf('============================================================\n');

fprintf('Matched N=117 beta = %.6f\n', ...
    matchedBeta);

fprintf('Expanded N=124 beta = %.6f\n', ...
    expandedBeta);

fprintf('Difference expanded - matched = %.6f\n', ...
    expandedBeta-matchedBeta);

fprintf('Same direction = %d\n', ...
    sign(expandedBeta)==sign(matchedBeta));

%% =========================================================================
% 15. OPTIONAL STANDARDIZED DESCRIPTIVE EFFECT SIZE
%
% Hedges-type g using raw group difference.
% This is descriptive only.
%% =========================================================================

gMatched = ...
    hedges_g( ...
        matchedHC, ...
        matchedPD);

gExpanded = ...
    hedges_g( ...
        expandedHC, ...
        expandedPD);

fprintf('\n============================================================\n');
fprintf('DESCRIPTIVE STANDARDIZED EFFECTS\n');
fprintf('============================================================\n');

fprintf('Matched N=117 Hedges g  = %.4f\n', ...
    gMatched);

fprintf('Expanded N=124 Hedges g = %.4f\n', ...
    gExpanded);

fprintf('(Descriptive only; primary inferential quantity is adjusted beta.)\n');

%% =========================================================================
% 16. CREATE SUMMARY TABLE
%% =========================================================================

Cohort = [ ...
    "Primary regional network reference"
    "Alternative matched N117"
    "Alternative expanded N124"];

Ntotal = [ ...
    117
    117
    124];

NHC = [ ...
    22
    22
    23];

NPD = [ ...
    95
    95
    101];

GroupBeta = [ ...
    primaryBeta
    matchedBeta
    expandedBeta];

CI95_Low = [ ...
    primaryCI_low
    matchedCI_low
    expandedCI_low];

CI95_High = [ ...
    primaryCI_high
    matchedCI_high
    expandedCI_high];

Pvalue = [ ...
    0.00066426
    matchedP
    expandedP];

SummaryTable = table( ...
    Cohort, ...
    Ntotal, ...
    NHC, ...
    NPD, ...
    GroupBeta, ...
    CI95_Low, ...
    CI95_High, ...
    Pvalue);

fprintf('\n============================================================\n');
fprintf('FINAL SUMMARY TABLE\n');
fprintf('============================================================\n');

disp(SummaryTable);

%% =========================================================================
% 17. SAVE MATCHED SUBJECT-LEVEL DATA
%% =========================================================================

MatchedSubjectTable = table( ...
    subjectsMatched, ...
    groupsMatched, ...
    ageMatched, ...
    SexMatched, ...
    wholeThetaMatched, ...
    'VariableNames',{ ...
        'Subject', ...
        'Group', ...
        'Age', ...
        'Sex_M1', ...
        'AlternativeWholeThetaDelta'});

writetable( ...
    MatchedSubjectTable, ...
    fullfile( ...
        altRoot, ...
        'Paper3_ALT11_whole_theta_MATCHED_N117_subjects.csv'));

%% =========================================================================
% 18. SAVE EXPANDED SUBJECT-LEVEL DATA
%% =========================================================================

ExpandedSubjectTable = table( ...
    subjectsExpanded, ...
    groupsExpanded, ...
    ageExpanded, ...
    SexExpanded, ...
    wholeThetaExpanded, ...
    'VariableNames',{ ...
        'Subject', ...
        'Group', ...
        'Age', ...
        'Sex_M1', ...
        'AlternativeWholeThetaDelta'});

writetable( ...
    ExpandedSubjectTable, ...
    fullfile( ...
        altRoot, ...
        'Paper3_ALT11_whole_theta_EXPANDED_N124_subjects.csv'));

%% =========================================================================
% 19. SAVE SUMMARY CSV / MAT
%% =========================================================================

summaryCSV = fullfile( ...
    altRoot, ...
    'Paper3_ALT11_whole_theta_regression_summary.csv');

writetable( ...
    SummaryTable, ...
    summaryCSV);

save( ...
    fullfile( ...
        altRoot, ...
        'Paper3_ALT11_whole_theta_regression_results.mat'), ...
    'SummaryTable', ...
    'MatchedSubjectTable', ...
    'ExpandedSubjectTable', ...
    'Rmatched', ...
    'Rexpanded', ...
    'matchedBeta', ...
    'expandedBeta', ...
    'primaryBeta', ...
    'gMatched', ...
    'gExpanded', ...
    'percentMagnitudeDifference', ...
    'magnitudeRatio');

%% =========================================================================
% 20. FINAL INTERPRETATION FLAGS
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('ROBUSTNESS FLAGS\n');
fprintf('============================================================\n');

fprintf('Matched alternative beta negative?       %d\n', ...
    matchedBeta < 0);

fprintf('Matched alternative CI excludes zero?    %d\n', ...
    matchedCI_high < 0 || matchedCI_low > 0);

fprintf('Matched alternative p < 0.05?            %d\n', ...
    matchedP < 0.05);

fprintf('Same direction as primary network?        %d\n', ...
    sameDirection);

fprintf('Expanded alternative beta negative?      %d\n', ...
    expandedBeta < 0);

fprintf('Expanded alternative CI excludes zero?   %d\n', ...
    expandedCI_high < 0 || expandedCI_low > 0);

fprintf('Expanded alternative p < 0.05?           %d\n', ...
    expandedP < 0.05);

fprintf('\nSaved:\n%s\n', ...
    summaryCSV);

fprintf('\n============================================================\n');
fprintf('STEP 17 COMPLETE\n');
fprintf('============================================================\n');

%% =========================================================================
%% LOCAL FUNCTION — HC3 REGRESSION
%% =========================================================================

function R = hc3_regression(y,X)

    y = double(y(:));
    X = double(X);

    [N,p] = size(X);

    if numel(y) ~= N
        error('X/y dimension mismatch.');
    end

    if rank(X) < p
        error('Design matrix is rank deficient.');
    end

    % OLS coefficient estimates.
    beta = X\y;

    fitted = X*beta;

    residual = y-fitted;

    % (X'X)^-1
    XtX_inv = pinv(X'*X);

    % Hat diagonal.
    h = ...
        sum( ...
            (X*XtX_inv).*X, ...
            2);

    if any(h >= 1)
        error('Invalid leverage >= 1.');
    end

    % HC3 residual inflation.
    omega = ...
        residual.^2 ./ ...
        (1-h).^2;

    % Sandwich meat.
    meat = ...
        X' * ...
        (X .* omega);

    % HC3 covariance.
    V_HC3 = ...
        XtX_inv * ...
        meat * ...
        XtX_inv;

    SE = ...
        sqrt(diag(V_HC3));

    % Finite-sample t inference.
    df = ...
        N-p;

    tStatistic = ...
        beta ./ SE;

    pValue = ...
        2 * ...
        tcdf( ...
            -abs(tStatistic), ...
            df);

    tCritical = ...
        tinv( ...
            0.975, ...
            df);

    CI_low = ...
        beta - ...
        tCritical*SE;

    CI_high = ...
        beta + ...
        tCritical*SE;

    % Model description.
    SST = ...
        sum( ...
            (y-mean(y)).^2);

    SSE = ...
        sum(residual.^2);

    R2 = ...
        1-SSE/SST;

    R = struct();

    R.beta = beta;
    R.SE = SE;
    R.CI_low = CI_low;
    R.CI_high = CI_high;
    R.t = tStatistic;
    R.p = pValue;
    R.df = df;
    R.V_HC3 = V_HC3;
    R.fitted = fitted;
    R.residual = residual;
    R.leverage = h;
    R.R2 = R2;

end

%% =========================================================================
%% LOCAL FUNCTION — COLUMN TO DOUBLE
%% =========================================================================

function x = column_to_double(x)

    if isnumeric(x)

        x = double(x(:));

    else

        x = ...
            str2double( ...
                string(x(:)));

    end

    if any(~isfinite(x))
        error('Numeric conversion produced NaN/Inf.');
    end

end

%% =========================================================================
%% LOCAL FUNCTION — HEDGES G
%% =========================================================================

function g = hedges_g(HC,PD)

    HC = double(HC(:));
    PD = double(PD(:));

    HC = HC(isfinite(HC));
    PD = PD(isfinite(PD));

    nHC = numel(HC);
    nPD = numel(PD);

    varHC = var(HC,0);
    varPD = var(PD,0);

    pooledSD = ...
        sqrt( ...
            ((nHC-1)*varHC + ...
             (nPD-1)*varPD) / ...
            (nHC+nPD-2));

    if pooledSD==0
        g = NaN;
        return;
    end

    d = ...
        (mean(PD)-mean(HC)) / ...
        pooledSD;

    df = ...
        nHC+nPD-2;

    J = ...
        1 - ...
        3/(4*df-1);

    g = ...
        J*d;

end
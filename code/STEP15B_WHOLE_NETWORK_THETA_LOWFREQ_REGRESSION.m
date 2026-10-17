%% =========================================================================
% PAPER 3 — STEP 15B
% WHOLE-NETWORK THETA REGRESSION
% LOW-FREQUENCY ARTIFACT SENSITIVITY
%
% Purpose
% -------
% Test whether the non-edge-selected whole-network theta Rest-to-Walk
% group difference remains after adjustment for the subject-level
% 1–3 Hz Walk/Rest power ratio derived in STEP 15.
%
% Outcome:
%   Mean theta Delta db-wPLI across all 55 regional edges
%
% Baseline model:
%   WholeThetaDelta ~ Group + Age_c + Sex
%
% Low-frequency-adjusted model:
%   WholeThetaDelta ~ Group + Age_c + Sex + LowFreqRatio_z
%
% Group coding:
%   HC = 0
%   PD = 1
%
% Sex coding:
%   F = 0
%   M = 1
%
% Inference:
%   OLS coefficients with HC3 heteroscedasticity-consistent covariance.
%
% IMPORTANT:
% This is a supportive post-hoc robustness analysis.
% It does not replace the primary NBS inference.
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

connectivityFile = fullfile( ...
    datasetRoot, ...
    'PAPER3_FINAL_CONNECTIVITY', ...
    'Paper3_CONNECTIVITY_STATS_INPUT.mat');

lowFreqFile = fullfile( ...
    datasetRoot, ...
    'PAPER3_LOWFREQ_ARTIFACT_SENSITIVITY', ...
    'Paper3_low_frequency_artifact_metric.mat');

participantsFile = fullfile( ...
    datasetRoot, ...
    'participants.tsv');

outRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_LOWFREQ_ARTIFACT_SENSITIVITY');

if ~exist(outRoot,'dir')
    mkdir(outRoot);
end

if ~exist(connectivityFile,'file')
    error('Connectivity input file not found:\n%s',connectivityFile);
end

if ~exist(lowFreqFile,'file')
    error('STEP 15 low-frequency file not found:\n%s',lowFreqFile);
end

if ~exist(participantsFile,'file')
    error('participants.tsv not found:\n%s',participantsFile);
end

%% =========================================================================
% 2. LOAD CONNECTIVITY DATA
%% =========================================================================

S = load(connectivityFile);

requiredFields = { ...
    'subjects', ...
    'groups', ...
    'deltaBio', ...
    'completeNetworkMask', ...
    'bandNames'};

for k = 1:numel(requiredFields)

    if ~isfield(S,requiredFields{k})
        error('Missing connectivity field: %s',requiredFields{k});
    end

end

subjectsAll = string(S.subjects(:));
groupsAll   = string(S.groups(:));

deltaBio = double(S.deltaBio);

completeMask = ...
    logical(S.completeNetworkMask(:));

bandNames = ...
    string(S.bandNames(:));

if ~isequal(size(deltaBio),[130 55 3])
    error('Expected deltaBio dimensions 130 x 55 x 3.');
end

if sum(completeMask) ~= 117
    error('Expected complete-network cohort N=117.');
end

thetaIndex = ...
    find(lower(bandNames)=="theta",1);

if isempty(thetaIndex)
    error('Theta band not found.');
end

%% =========================================================================
% 3. WHOLE-NETWORK THETA DELTA
%% =========================================================================

% Mean across all 55 regional theta edges.
wholeThetaAll = ...
    mean( ...
        deltaBio(:,:,thetaIndex), ...
        2, ...
        'omitnan');

wholeThetaAll = ...
    squeeze(wholeThetaAll);

% Complete-network cohort only.
subjects = subjectsAll(completeMask);
groups   = groupsAll(completeMask);

wholeTheta = ...
    wholeThetaAll(completeMask);

if any(~isfinite(wholeTheta))
    error('Nonfinite whole-network theta values in N=117 cohort.');
end

N = numel(wholeTheta);

if N ~= 117
    error('Expected N=117.');
end

%% =========================================================================
% 4. LOAD LOW-FREQUENCY ARTIFACT METRIC
%% =========================================================================

L = load(lowFreqFile);

if isfield(L,'LowFreqTable')

    LF = L.LowFreqTable;

    lfSubjects = ...
        string(LF.Subject);

    lfRatio = ...
        double(LF.WalkRestRatio_1to3Hz);

elseif isfield(L,'lowFreqRatio')

    lfSubjects = subjectsAll;

    lfRatio = ...
        double(L.lowFreqRatio(:));

else

    error([ ...
        'Could not find LowFreqTable or lowFreqRatio ' ...
        'in STEP 15 output.']);

end

[tf,locLF] = ...
    ismember(subjects,lfSubjects);

if ~all(tf)

    fprintf('\nMissing low-frequency subjects:\n');
    disp(subjects(~tf));

    error('Low-frequency data could not be aligned.');
end

lowFreqRatio = ...
    lfRatio(locLF);

if any(~isfinite(lowFreqRatio))
    error('Low-frequency ratio contains NaN/Inf.');
end

% Same standardization as STEP 15:
% within the N=117 primary complete-network cohort.
LowFreqRatio_z = ...
    (lowFreqRatio - mean(lowFreqRatio)) ./ ...
    std(lowFreqRatio);

%% =========================================================================
% 5. LOAD AGE AND SEX
%% =========================================================================

P = readtable( ...
    participantsFile, ...
    'FileType','text', ...
    'Delimiter','\t', ...
    'TextType','string');

[tf,loc] = ...
    ismember( ...
        subjects, ...
        string(P.participant_id));

if ~all(tf)
    error('Could not align all subjects to participants.tsv.');
end

participantGroup = ...
    string(P.group(loc));

if ~isequal(participantGroup,groups)
    error('Group mismatch with participants.tsv.');
end

age = ...
    column_to_double(P.age(loc));

sexRaw = ...
    upper( ...
        strtrim( ...
            string(P.sex(loc))));

if any(~ismember(sexRaw,["F","M"]))
    error('Unexpected sex coding.');
end

Group = ...
    double(groups=="PD");

Sex = ...
    double(sexRaw=="M");

Age_c = ...
    age - mean(age);

nHC = sum(Group==0);
nPD = sum(Group==1);

if nHC ~= 22 || nPD ~= 95
    error('Expected 22 HC + 95 PD.');
end

%% =========================================================================
% 6. DESCRIPTIVE WHOLE-NETWORK THETA VALUES
%% =========================================================================

meanHC = ...
    mean(wholeTheta(Group==0));

meanPD = ...
    mean(wholeTheta(Group==1));

sdHC = ...
    std(wholeTheta(Group==0));

sdPD = ...
    std(wholeTheta(Group==1));

fprintf('\n============================================================\n');
fprintf('WHOLE-NETWORK THETA DESCRIPTIVE VALUES\n');
fprintf('============================================================\n');

fprintf('N = %d = %d HC + %d PD\n',N,nHC,nPD);

fprintf('HC mean Delta = %.6f (SD %.6f)\n', ...
    meanHC,sdHC);

fprintf('PD mean Delta = %.6f (SD %.6f)\n', ...
    meanPD,sdPD);

fprintf('Raw PD-HC difference = %.6f\n', ...
    meanPD-meanHC);

%% =========================================================================
% 7. BASELINE MODEL
%
% WholeTheta ~ Group + Age_c + Sex
%% =========================================================================

X0 = [ ...
    ones(N,1), ...
    Group, ...
    Age_c, ...
    Sex];

baseline = ...
    hc3_regression( ...
        wholeTheta, ...
        X0);

betaGroup0 = ...
    baseline.beta(2);

seGroup0 = ...
    baseline.SE(2);

ciLow0 = ...
    baseline.CI_low(2);

ciHigh0 = ...
    baseline.CI_high(2);

pGroup0 = ...
    baseline.p(2);

%% =========================================================================
% 8. LOW-FREQUENCY-ADJUSTED MODEL
%
% WholeTheta ~ Group + Age_c + Sex + LowFreqRatio_z
%% =========================================================================

X1 = [ ...
    ones(N,1), ...
    Group, ...
    Age_c, ...
    Sex, ...
    LowFreqRatio_z];

adjusted = ...
    hc3_regression( ...
        wholeTheta, ...
        X1);

betaGroup1 = ...
    adjusted.beta(2);

seGroup1 = ...
    adjusted.SE(2);

ciLow1 = ...
    adjusted.CI_low(2);

ciHigh1 = ...
    adjusted.CI_high(2);

pGroup1 = ...
    adjusted.p(2);

% Low-frequency nuisance coefficient.
betaLowFreq = ...
    adjusted.beta(5);

ciLowFreqLow = ...
    adjusted.CI_low(5);

ciLowFreqHigh = ...
    adjusted.CI_high(5);

pLowFreq = ...
    adjusted.p(5);

%% =========================================================================
% 9. EFFECT-STABILITY METRICS
%% =========================================================================

absoluteBetaChange = ...
    betaGroup1-betaGroup0;

percentMagnitudeChange = ...
    100 * ...
    (abs(betaGroup1)-abs(betaGroup0)) / ...
    abs(betaGroup0);

percentAttenuation = ...
    100 * ...
    (abs(betaGroup0)-abs(betaGroup1)) / ...
    abs(betaGroup0);

%% =========================================================================
% 10. RESULTS
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('BASELINE WHOLE-NETWORK THETA MODEL\n');
fprintf('Delta ~ Group + Age_c + Sex\n');
fprintf('HC3 heteroscedasticity-consistent inference\n');
fprintf('============================================================\n');

fprintf('Adjusted PD-HC beta = %.6f\n',betaGroup0);

fprintf('HC3 SE              = %.6f\n',seGroup0);

fprintf('95%% CI               = [%.6f, %.6f]\n', ...
    ciLow0,ciHigh0);

fprintf('p                    = %.8f\n',pGroup0);

fprintf('\n============================================================\n');
fprintf('LOW-FREQUENCY-ADJUSTED WHOLE-NETWORK THETA MODEL\n');
fprintf('Delta ~ Group + Age_c + Sex + LowFreqRatio_z\n');
fprintf('HC3 heteroscedasticity-consistent inference\n');
fprintf('============================================================\n');

fprintf('Adjusted PD-HC beta = %.6f\n',betaGroup1);

fprintf('HC3 SE              = %.6f\n',seGroup1);

fprintf('95%% CI               = [%.6f, %.6f]\n', ...
    ciLow1,ciHigh1);

fprintf('p                    = %.8f\n',pGroup1);

fprintf('\nLow-frequency nuisance coefficient:\n');

fprintf('beta per 1 SD increase = %.6f\n',betaLowFreq);

fprintf('95%% CI                 = [%.6f, %.6f]\n', ...
    ciLowFreqLow,ciLowFreqHigh);

fprintf('p                       = %.8f\n',pLowFreq);

fprintf('\n============================================================\n');
fprintf('EFFECT STABILITY\n');
fprintf('============================================================\n');

fprintf('Baseline Group beta       = %.6f\n',betaGroup0);
fprintf('Low-frequency-adjusted beta = %.6f\n',betaGroup1);

fprintf('Absolute beta change      = %.6f\n',absoluteBetaChange);

fprintf('Percent attenuation in |beta| = %.2f%%\n', ...
    percentAttenuation);

fprintf('Percent change in |beta|      = %.2f%%\n', ...
    percentMagnitudeChange);

%% =========================================================================
% 11. OPTIONAL CORRELATION / GROUP DIFFERENCE FOR LOW-FREQUENCY METRIC
%% =========================================================================

[rPearson,pPearson] = ...
    corr( ...
        lowFreqRatio, ...
        wholeTheta, ...
        'Type','Pearson', ...
        'Rows','complete');

[rSpearman,pSpearman] = ...
    corr( ...
        lowFreqRatio, ...
        wholeTheta, ...
        'Type','Spearman', ...
        'Rows','complete');

fprintf('\n============================================================\n');
fprintf('DESCRIPTIVE LOW-FREQUENCY ASSOCIATION\n');
fprintf('============================================================\n');

fprintf('Pearson r(low-frequency ratio, theta Delta) = %.4f, p = %.6g\n', ...
    rPearson,pPearson);

fprintf('Spearman rho = %.4f, p = %.6g\n', ...
    rSpearman,pSpearman);

fprintf('HC median low-frequency ratio = %.4f\n', ...
    median(lowFreqRatio(Group==0)));

fprintf('PD median low-frequency ratio = %.4f\n', ...
    median(lowFreqRatio(Group==1)));

%% =========================================================================
% 12. SAVE TABLE
%% =========================================================================

Model = [ ...
    "Baseline_AgeSex"
    "Adjusted_AgeSex_LowFreq"];

GroupBeta = [ ...
    betaGroup0
    betaGroup1];

HC3_SE = [ ...
    seGroup0
    seGroup1];

CI95_Low = [ ...
    ciLow0
    ciLow1];

CI95_High = [ ...
    ciHigh0
    ciHigh1];

P = [ ...
    pGroup0
    pGroup1];

ResultsTable = table( ...
    Model, ...
    GroupBeta, ...
    HC3_SE, ...
    CI95_Low, ...
    CI95_High, ...
    P);

disp(ResultsTable);

outCSV = fullfile( ...
    outRoot, ...
    'Paper3_whole_network_theta_lowfreq_regression.csv');

writetable( ...
    ResultsTable, ...
    outCSV);

%% =========================================================================
% 13. SAVE SUBJECT-LEVEL INPUTS
%% =========================================================================

SubjectTable = table( ...
    subjects, ...
    groups, ...
    age, ...
    Sex, ...
    wholeTheta, ...
    lowFreqRatio, ...
    LowFreqRatio_z, ...
    'VariableNames',{ ...
        'Subject', ...
        'Group', ...
        'Age', ...
        'Sex_M1', ...
        'WholeNetworkThetaDelta', ...
        'WalkRestPowerRatio_1to3Hz', ...
        'WalkRestPowerRatio_1to3Hz_Z'});

writetable( ...
    SubjectTable, ...
    fullfile( ...
        outRoot, ...
        'Paper3_whole_network_theta_lowfreq_subject_data.csv'));

%% =========================================================================
% 14. SAVE MAT
%% =========================================================================

save( ...
    fullfile( ...
        outRoot, ...
        'Paper3_whole_network_theta_lowfreq_regression.mat'), ...
    'SubjectTable', ...
    'ResultsTable', ...
    'baseline', ...
    'adjusted', ...
    'betaGroup0', ...
    'betaGroup1', ...
    'percentAttenuation', ...
    'betaLowFreq', ...
    'pLowFreq');

fprintf('\nSaved:\n%s\n',outCSV);

fprintf('\n============================================================\n');
fprintf('STEP 15B COMPLETE\n');
fprintf('============================================================\n');

%% =========================================================================
%% LOCAL FUNCTIONS
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

    % --------------------------------------------------------------
    % OLS coefficient estimates.
    % --------------------------------------------------------------

    beta = X\y;

    fitted = X*beta;

    residual = y-fitted;

    % --------------------------------------------------------------
    % Hat-matrix diagonal.
    % --------------------------------------------------------------

    XtX_inv = ...
        pinv(X'*X);

    h = ...
        sum( ...
            (X*XtX_inv).*X, ...
            2);

    if any(h >= 1)
        error('Invalid leverage >= 1.');
    end

    % --------------------------------------------------------------
    % HC3 covariance:
    %
    % residual_i^2 / (1-h_i)^2
    % --------------------------------------------------------------

    omega = ...
        residual.^2 ./ ...
        (1-h).^2;

    meat = ...
        X' * ...
        (X .* omega);

    V_HC3 = ...
        XtX_inv * ...
        meat * ...
        XtX_inv;

    SE = ...
        sqrt( ...
            diag(V_HC3));

    % --------------------------------------------------------------
    % Robust t statistics.
    %
    % Use residual df = N-p for finite-sample inference.
    % --------------------------------------------------------------

    df = N-p;

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

    % --------------------------------------------------------------
    % Conventional R-squared retained only as model description.
    % --------------------------------------------------------------

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

function x = column_to_double(x)

    if isnumeric(x)

        x = double(x(:));

    else

        x = ...
            str2double( ...
                string(x(:)));

    end

end
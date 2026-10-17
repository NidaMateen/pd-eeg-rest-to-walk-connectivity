%% =========================================================================
% PAPER 3 — STEP 15C
% WHOLE-NETWORK THETA REGRESSION
% LOG-TRANSFORMED LOW-FREQUENCY POWER-RATIO SENSITIVITY
%
% Purpose
% -------
% Check whether the whole-network theta PD-HC effect remains after
% adjustment for:
%
%   log(Walk 1-3 Hz power / Rest 1-3 Hz power)
%
% The log transformation is useful because power ratios are positive,
% multiplicative, and often right-skewed.
%
% Outcome:
%   Mean theta Delta db-wPLI across all 55 regional edges
%
% Models:
%
%   Baseline:
%     WholeThetaDelta ~ Group + Age_c + Sex
%
%   Log-low-frequency adjusted:
%     WholeThetaDelta ~ Group + Age_c + Sex + LogLowFreqRatio_z
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
%   OLS coefficients with HC3 heteroscedasticity-consistent covariance.
%
% This is a post-hoc robustness analysis only.
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
    error('Connectivity file not found:\n%s',connectivityFile);
end

if ~exist(lowFreqFile,'file')
    error('Low-frequency file not found:\n%s',lowFreqFile);
end

if ~exist(participantsFile,'file')
    error('participants.tsv not found:\n%s',participantsFile);
end

%% =========================================================================
% 2. LOAD CONNECTIVITY DATA
%% =========================================================================

S = load(connectivityFile);

subjectsAll = string(S.subjects(:));
groupsAll   = string(S.groups(:));

deltaBio = double(S.deltaBio);

completeMask = ...
    logical(S.completeNetworkMask(:));

bandNames = ...
    string(S.bandNames(:));

if ~isequal(size(deltaBio),[130 55 3])
    error('Expected deltaBio = 130 x 55 x 3.');
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

wholeThetaAll = ...
    squeeze( ...
        mean( ...
            deltaBio(:,:,thetaIndex), ...
            2, ...
            'omitnan'));

subjects = ...
    subjectsAll(completeMask);

groups = ...
    groupsAll(completeMask);

wholeTheta = ...
    wholeThetaAll(completeMask);

if any(~isfinite(wholeTheta))
    error('Nonfinite whole-network theta values.');
end

N = numel(wholeTheta);

if N ~= 117
    error('Expected N=117.');
end

%% =========================================================================
% 4. LOAD LOW-FREQUENCY POWER RATIO
%% =========================================================================

L = load(lowFreqFile);

if isfield(L,'LowFreqTable')

    LF = L.LowFreqTable;

    lfSubjects = ...
        string(LF.Subject);

    lfRatio = ...
        double(LF.WalkRestRatio_1to3Hz);

elseif isfield(L,'lowFreqRatio')

    lfSubjects = ...
        subjectsAll;

    lfRatio = ...
        double(L.lowFreqRatio(:));

else

    error('Low-frequency ratio not found.');

end

[tf,locLF] = ...
    ismember(subjects,lfSubjects);

if ~all(tf)

    disp(subjects(~tf));

    error('Could not align low-frequency data.');

end

lowFreqRatio = ...
    lfRatio(locLF);

if any(~isfinite(lowFreqRatio))
    error('Low-frequency ratio contains NaN/Inf.');
end

% Ratios must be positive before log transform.
if any(lowFreqRatio <= 0)

    error([ ...
        'Low-frequency Walk/Rest ratios must all be >0 ' ...
        'before logarithmic transformation.']);

end

%% =========================================================================
% 5. LOG TRANSFORM + STANDARDIZE
%% =========================================================================

logLowFreqRatio = ...
    log(lowFreqRatio);

LogLowFreqRatio_z = ...
    ( ...
    logLowFreqRatio - ...
    mean(logLowFreqRatio) ...
    ) ./ ...
    std(logLowFreqRatio);

if any(~isfinite(LogLowFreqRatio_z))
    error('Log-transformed artifact covariate contains NaN/Inf.');
end

fprintf('\n============================================================\n');
fprintf('LOW-FREQUENCY RATIO DISTRIBUTION\n');
fprintf('============================================================\n');

fprintf('Raw ratio median = %.4f\n', ...
    median(lowFreqRatio));

fprintf('Raw ratio IQR    = [%.4f, %.4f]\n', ...
    percentile_local(lowFreqRatio,0.25), ...
    percentile_local(lowFreqRatio,0.75));

fprintf('Log ratio mean   = %.4f\n', ...
    mean(logLowFreqRatio));

fprintf('Log ratio SD     = %.4f\n', ...
    std(logLowFreqRatio));

%% =========================================================================
% 6. AGE / SEX
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
    error('Could not align participants.tsv.');
end

participantGroup = ...
    string(P.group(loc));

if ~isequal(participantGroup,groups)
    error('Group mismatch.');
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
% 7. DESCRIPTIVE THETA VALUES
%% =========================================================================

meanHC = ...
    mean(wholeTheta(Group==0));

meanPD = ...
    mean(wholeTheta(Group==1));

fprintf('\n============================================================\n');
fprintf('WHOLE-NETWORK THETA DESCRIPTIVE VALUES\n');
fprintf('============================================================\n');

fprintf('N = %d = %d HC + %d PD\n', ...
    N,nHC,nPD);

fprintf('HC mean Delta = %.6f\n',meanHC);
fprintf('PD mean Delta = %.6f\n',meanPD);

fprintf('Raw PD-HC difference = %.6f\n', ...
    meanPD-meanHC);

%% =========================================================================
% 8. BASELINE MODEL
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

ciLow0 = ...
    baseline.CI_low(2);

ciHigh0 = ...
    baseline.CI_high(2);

pGroup0 = ...
    baseline.p(2);

%% =========================================================================
% 9. LOG-RATIO-ADJUSTED MODEL
%
% WholeTheta ~ Group + Age_c + Sex + LogLowFreqRatio_z
%% =========================================================================

Xlog = [ ...
    ones(N,1), ...
    Group, ...
    Age_c, ...
    Sex, ...
    LogLowFreqRatio_z];

logAdjusted = ...
    hc3_regression( ...
        wholeTheta, ...
        Xlog);

betaGroupLog = ...
    logAdjusted.beta(2);

seGroupLog = ...
    logAdjusted.SE(2);

ciLowLog = ...
    logAdjusted.CI_low(2);

ciHighLog = ...
    logAdjusted.CI_high(2);

pGroupLog = ...
    logAdjusted.p(2);

% Nuisance coefficient.
betaLogFreq = ...
    logAdjusted.beta(5);

ciLogFreqLow = ...
    logAdjusted.CI_low(5);

ciLogFreqHigh = ...
    logAdjusted.CI_high(5);

pLogFreq = ...
    logAdjusted.p(5);

%% =========================================================================
% 10. EFFECT STABILITY
%% =========================================================================

absoluteBetaChange = ...
    betaGroupLog-betaGroup0;

percentAttenuation = ...
    100 * ...
    ( ...
    abs(betaGroup0) - ...
    abs(betaGroupLog) ...
    ) / ...
    abs(betaGroup0);

%% =========================================================================
% 11. RESULTS
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('BASELINE MODEL\n');
fprintf('Delta ~ Group + Age_c + Sex\n');
fprintf('HC3 inference\n');
fprintf('============================================================\n');

fprintf('Adjusted PD-HC beta = %.6f\n', ...
    betaGroup0);

fprintf('95%% CI               = [%.6f, %.6f]\n', ...
    ciLow0,ciHigh0);

fprintf('p                    = %.8f\n', ...
    pGroup0);

fprintf('\n============================================================\n');
fprintf('LOG-LOW-FREQUENCY-ADJUSTED MODEL\n');
fprintf('Delta ~ Group + Age_c + Sex + z(log(Walk/Rest 1-3 Hz power))\n');
fprintf('HC3 inference\n');
fprintf('============================================================\n');

fprintf('Adjusted PD-HC beta = %.6f\n', ...
    betaGroupLog);

fprintf('HC3 SE              = %.6f\n', ...
    seGroupLog);

fprintf('95%% CI               = [%.6f, %.6f]\n', ...
    ciLowLog,ciHighLog);

fprintf('p                    = %.8f\n', ...
    pGroupLog);

fprintf('\nLog low-frequency nuisance coefficient:\n');

fprintf('beta per 1 SD increase = %.6f\n', ...
    betaLogFreq);

fprintf('95%% CI                 = [%.6f, %.6f]\n', ...
    ciLogFreqLow,ciLogFreqHigh);

fprintf('p                       = %.8f\n', ...
    pLogFreq);

fprintf('\n============================================================\n');
fprintf('EFFECT STABILITY\n');
fprintf('============================================================\n');

fprintf('Baseline Group beta       = %.6f\n', ...
    betaGroup0);

fprintf('Log-adjusted Group beta   = %.6f\n', ...
    betaGroupLog);

fprintf('Absolute beta change      = %.6f\n', ...
    absoluteBetaChange);

fprintf('Percent attenuation in |beta| = %.2f%%\n', ...
    percentAttenuation);

%% =========================================================================
% 12. CORRELATIONS
%% =========================================================================

[rRawPearson,pRawPearson] = ...
    corr( ...
        lowFreqRatio, ...
        wholeTheta, ...
        'Type','Pearson', ...
        'Rows','complete');

[rLogPearson,pLogPearson] = ...
    corr( ...
        logLowFreqRatio, ...
        wholeTheta, ...
        'Type','Pearson', ...
        'Rows','complete');

[rLogSpearman,pLogSpearman] = ...
    corr( ...
        logLowFreqRatio, ...
        wholeTheta, ...
        'Type','Spearman', ...
        'Rows','complete');

fprintf('\n============================================================\n');
fprintf('CORRELATION CHECK\n');
fprintf('============================================================\n');

fprintf('Raw ratio Pearson r = %.4f, p = %.6g\n', ...
    rRawPearson,pRawPearson);

fprintf('Log ratio Pearson r = %.4f, p = %.6g\n', ...
    rLogPearson,pLogPearson);

fprintf('Log ratio Spearman rho = %.4f, p = %.6g\n', ...
    rLogSpearman,pLogSpearman);

%% =========================================================================
% 13. SAVE RESULTS
%% =========================================================================

Model = [ ...
    "Baseline_AgeSex"
    "Adjusted_AgeSex_LogLowFreq"];

GroupBeta = [ ...
    betaGroup0
    betaGroupLog];

CI95_Low = [ ...
    ciLow0
    ciLowLog];

CI95_High = [ ...
    ciHigh0
    ciHighLog];

Pvalue = [ ...
    pGroup0
    pGroupLog];

ResultsTable = table( ...
    Model, ...
    GroupBeta, ...
    CI95_Low, ...
    CI95_High, ...
    Pvalue);

disp(ResultsTable);

outCSV = fullfile( ...
    outRoot, ...
    'Paper3_whole_network_theta_LOG_lowfreq_regression.csv');

writetable( ...
    ResultsTable, ...
    outCSV);

%% =========================================================================
% 14. SAVE SUBJECT-LEVEL INPUTS
%% =========================================================================

SubjectTable = table( ...
    subjects, ...
    groups, ...
    age, ...
    Sex, ...
    wholeTheta, ...
    lowFreqRatio, ...
    logLowFreqRatio, ...
    LogLowFreqRatio_z, ...
    'VariableNames',{ ...
        'Subject', ...
        'Group', ...
        'Age', ...
        'Sex_M1', ...
        'WholeNetworkThetaDelta', ...
        'WalkRestPowerRatio_1to3Hz', ...
        'LogWalkRestPowerRatio_1to3Hz', ...
        'LogWalkRestPowerRatio_1to3Hz_Z'});

writetable( ...
    SubjectTable, ...
    fullfile( ...
        outRoot, ...
        'Paper3_whole_network_theta_LOG_lowfreq_subject_data.csv'));

save( ...
    fullfile( ...
        outRoot, ...
        'Paper3_whole_network_theta_LOG_lowfreq_regression.mat'), ...
    'SubjectTable', ...
    'ResultsTable', ...
    'baseline', ...
    'logAdjusted', ...
    'betaGroup0', ...
    'betaGroupLog', ...
    'percentAttenuation', ...
    'betaLogFreq', ...
    'pLogFreq');

fprintf('\nSaved:\n%s\n',outCSV);

fprintf('\n============================================================\n');
fprintf('STEP 15C COMPLETE\n');
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
        error('Design matrix rank deficient.');
    end

    beta = X\y;

    fitted = X*beta;

    residual = y-fitted;

    XtX_inv = ...
        pinv(X'*X);

    h = ...
        sum( ...
            (X*XtX_inv).*X, ...
            2);

    if any(h >= 1)
        error('Invalid leverage >=1.');
    end

    % HC3 covariance.
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
        sqrt(diag(V_HC3));

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

%% =========================================================================

function q = percentile_local(x,p)

    x = ...
        sort( ...
            double(x(:)));

    x = ...
        x(isfinite(x));

    if isempty(x)

        q = NaN;
        return;

    end

    idx = ...
        1 + ...
        (numel(x)-1)*p;

    lo = floor(idx);
    hi = ceil(idx);

    if lo==hi

        q = x(lo);

    else

        w = idx-lo;

        q = ...
            (1-w)*x(lo) + ...
            w*x(hi);

    end

end
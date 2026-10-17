%% =========================================================================
% PAPER 3 — STEP 9B
% VALIDATE PRIMARY NBS EXECUTION + ESTIMATE EFFECT MAGNITUDE
%
% PURPOSE
% -------
% 1) Verify that the very fast 10,000-resample NBS run actually generated
%    a non-degenerate null distribution and that selected resamples can be
%    reproduced exactly from the frozen seed/sign-flip scheme.
%
% 2) Quantify effect magnitude without changing the primary inference:
%    A. Primary NBS component-average effect (DESCRIPTIVE / POST-SELECTION)
%    B. Whole-network band-average effects for theta/alpha/beta
%
% IMPORTANT
% ---------
% The NBS component p-value is the confirmatory primary inference.
% NBS itself does not have a single conventional Cohen-d/Hedges-g effect
% size or confidence interval. Any scalar effect computed by averaging
% edges inside the data-selected significant component is post-selection
% and should be labelled descriptive, not an independent confirmatory CI.
%
% This script DOES NOT:
%   - change the NBS threshold
%   - rerun feature selection
%   - alter the cohort
%   - change bands
%   - alter multiplicity correction
%
% It is a read-only validation/effect-estimation step.
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

connectivityRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_FINAL_CONNECTIVITY');

nbsRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_PRIMARY_NBS');

statsInputFile = fullfile( ...
    connectivityRoot, ...
    'Paper3_CONNECTIVITY_STATS_INPUT.mat');

resultFile = fullfile( ...
    nbsRoot, ...
    'Paper3_PRIMARY_NBS_RESULTS.mat');

checkpointFile = fullfile( ...
    nbsRoot, ...
    'Paper3_PRIMARY_NBS_resampling_checkpoint.mat');

participantsFile = fullfile( ...
    datasetRoot, ...
    'participants.tsv');

outFile = fullfile( ...
    nbsRoot, ...
    'Paper3_PRIMARY_NBS_effect_size_summary.csv');

%% =========================================================================
% 2. LOAD
%% =========================================================================

if ~exist(statsInputFile,'file'), error('Missing %s',statsInputFile); end
if ~exist(resultFile,'file'), error('Missing %s',resultFile); end
if ~exist(checkpointFile,'file'), error('Missing %s',checkpointFile); end
if ~exist(participantsFile,'file'), error('Missing %s',participantsFile); end

S = load(statsInputFile);
R = load(resultFile);
CP = load(checkpointFile);

%% =========================================================================
% 3. BASIC PRIMARY RESULT CHECK
%% =========================================================================

subjects = string(R.subjects(:));
groups   = string(R.groups(:));

N = numel(subjects);

if N ~= 117 || sum(groups=="HC") ~= 22 || sum(groups=="PD") ~= 95
    error('Frozen primary cohort mismatch.');
end

if R.N_RESAMPLES ~= 10000 || R.RNG_SEED ~= 20260916 || R.NBS_THRESHOLD ~= 2.5
    error('Frozen NBS parameters do not match protocol.');
end

nullMax = double(R.nullMaxExtent(:));

if numel(nullMax) ~= 10000 || any(~isfinite(nullMax))
    error('Null maximum distribution is incomplete/nonfinite.');
end

fprintf('\n============================================================\n');
fprintf('A. NBS EXECUTION VALIDATION\n');
fprintf('============================================================\n');

fprintf('Null resamples = %d\n',numel(nullMax));
fprintf('Unique null-max extents = %d\n',numel(unique(nullMax)));
fprintf('Null max: min=%g | median=%g | 95th=%g | max=%g\n', ...
    min(nullMax),median(nullMax),empirical_percentile(nullMax,0.95),max(nullMax));

%% =========================================================================
% 4. VERIFY OBSERVED COMPONENT / MONTE CARLO P
%% =========================================================================

CT = R.ComponentTable;

sig = CT.Significant_FWER_0p05;

if sum(sig) ~= 1
    error('Expected exactly one significant primary NBS component.');
end

sigRow = CT(sig,:);

fprintf('\nSignificant primary component:\n');
disp(sigRow(:,{'Band','Direction','ExtentEdges','JointFWER_P'}));

obsExtent = double(sigRow.ExtentEdges);

nExceed = sum(nullMax >= obsExtent);

pRecomputed = (1+nExceed)/(numel(nullMax)+1);

fprintf('Observed component extent = %d edges\n',obsExtent);
fprintf('Null maxima >= observed   = %d / %d\n',nExceed,numel(nullMax));
fprintf('Recomputed p_FWER         = %.8f\n',pRecomputed);

if abs(pRecomputed-double(sigRow.JointFWER_P)) > 1e-12
    error('Stored and recomputed joint FWER p-values do not match.');
end

% Monte Carlo uncertainty of the exceedance probability.
% Exact Clopper-Pearson interval for x exceedances among B resamples.
B = numel(nullMax);
x = nExceed;

if x == 0
    mcLo = 0;
else
    mcLo = betaincinv(0.025,x,B-x+1);
end

if x == B
    mcHi = 1;
else
    mcHi = betaincinv(0.975,x+1,B-x);
end

fprintf('Monte Carlo exceedance-rate 95%% CI = [%.6f, %.6f]\n',mcLo,mcHi);
fprintf(['NOTE: this CI describes Monte Carlo uncertainty of the resampling p-value,\n' ...
         'not uncertainty in the biological effect magnitude.\n']);

%% =========================================================================
% 5. RECONSTRUCT DESIGN AND FREEDMAN-LANE NULL INPUT
%% =========================================================================

deltaNBS = double(R.deltaNBS);

Y = reshape(deltaNBS,N,55*3);

Group = double(groups=="PD");
age    = double(R.age(:));
Age_c  = double(R.Age_c(:));
Sex    = double(R.Sex(:));

M = [ones(N,1),Group,Age_c,Sex];
Z = [ones(N,1),Age_c,Sex];

C = [0;1;0;0];

Mpinv = pinv(M);
Zpinv = pinv(Z);

Hfull = M*Mpinv;
Rfull = eye(N)-Hfull;
Rdiag = diag(Rfull);

vgHC = Group==0;
vgPD = Group==1;

sumR_HC = sum(Rdiag(vgHC));
sumR_PD = sum(Rdiag(vgPD));

A_HC = M(vgHC,:)'*M(vgHC,:);
A_PD = M(vgPD,:)'*M(vgPD,:);

FittedReduced   = Z*(Zpinv*Y);
ResidualReduced = Y-FittedReduced;

%% =========================================================================
% 6. VERIFY RADEMACHER MATRIX IS NON-DEGENERATE
%% =========================================================================

rng(R.RNG_SEED,'twister');

signMatrix = int8(ones(R.N_RESAMPLES,N));
flipMask = rand(R.N_RESAMPLES,N) < 0.5;
signMatrix(flipMask) = -1;

nNegative = sum(signMatrix==-1,2);

fprintf('\nRademacher sign-pattern diagnostics:\n');
fprintf('Negative signs per resample: min=%d | median=%.1f | max=%d\n', ...
    min(nNegative),median(nNegative),max(nNegative));
fprintf('Unique sign rows among 10,000 resamples = %d\n', ...
    size(unique(signMatrix,'rows'),1));

if size(unique(signMatrix,'rows'),1) < 9990
    warning('Unexpectedly many duplicated Rademacher sign patterns.');
end

%% =========================================================================
% 7. EXACTLY RECOMPUTE A FIXED SET OF SAVED NULL ITERATIONS
%
% This verifies that the saved parallel output corresponds to the intended
% resampling algorithm rather than merely trusting runtime.
%% =========================================================================

checkIdx = [ ...
    1 2 3 10 25 50 100 250 500 750 ...
    1000 1500 2000 3000 4000 5000 ...
    6000 7000 8000 9000 9500 9999 10000];

recomputed = nan(numel(checkIdx),1);

fprintf('\nRecomputing %d selected null iterations sequentially...\n',numel(checkIdx));

tCheck = tic;

for q = 1:numel(checkIdx)

    b = checkIdx(q);

    signs = double(signMatrix(b,:))';

    Ystar = FittedReduced + ResidualReduced.*signs;

    vStar = compute_signed_v_all_only( ...
        Ystar,M,Mpinv,C,2,vgHC,vgPD, ...
        sumR_HC,sumR_PD,A_HC,A_PD);

    recomputed(q) = joint_max_component_extent( ...
        vStar,55,3,R.NBS_THRESHOLD, ...
        double(R.edgeIndex1(:)), ...
        double(R.edgeIndex2(:)), ...
        11);

end

checkElapsed = toc(tCheck);

savedSelected = nullMax(checkIdx);

ValidationTable = table( ...
    checkIdx(:), ...
    savedSelected(:), ...
    recomputed(:), ...
    savedSelected(:)==recomputed(:), ...
    'VariableNames', ...
    {'Resample','SavedNullMax','RecomputedNullMax','ExactMatch'});

disp(ValidationTable);

if ~all(ValidationTable.ExactMatch)
    error('Selected null iterations do NOT reproduce exactly.');
end

fprintf('PASS: every selected null iteration reproduced exactly.\n');
fprintf('Sequential validation time = %.4f s for %d iterations\n', ...
    checkElapsed,numel(checkIdx));

%% =========================================================================
% 8. LOCATE SIGNIFICANT THETA NEGATIVE COMPONENT
%% =========================================================================

if ~strcmpi(string(sigRow.Band),"theta") || ...
   ~contains(string(sigRow.Direction),"Negative")

    error('Significant component is not the expected theta-negative component.');

end

compID = double(sigRow.ComponentID);

thetaNegComps = R.observedComponents{1,2};

if compID < 1 || compID > numel(thetaNegComps)
    error('Invalid significant component ID.');
end

edgeIdx = double(thetaNegComps(compID).edgeIdx(:));

if numel(edgeIdx) ~= obsExtent
    error('Component edge count does not match ComponentTable.');
end

fprintf('\n============================================================\n');
fprintf('B. EFFECT MAGNITUDE ESTIMATION\n');
fprintf('============================================================\n');

fprintf('Significant theta component contains %d edges.\n',numel(edgeIdx));

%% =========================================================================
% 9. ALIGN REST/WALK DATA TO PRIMARY COHORT
%% =========================================================================

subjectsAll = string(S.subjects(:));

[tf,loc] = ismember(subjects,subjectsAll);

if ~all(tf)
    error('Could not align primary NBS subjects to connectivity file.');
end

restBio = double(S.restBio(loc,:,:));
walkBio = double(S.walkBio(loc,:,:));

if any(~isfinite(restBio(:,edgeIdx,1)),'all') || ...
   any(~isfinite(walkBio(:,edgeIdx,1)),'all')
    error('Nonfinite theta component data in primary cohort.');
end

%% =========================================================================
% 10. COMPONENT-AVERAGE SUBJECT SUMMARY
%
% IMPORTANT:
% The component was selected from these same data by NBS.
% Therefore its scalar effect estimate is descriptive/post-selection.
%% =========================================================================

compRest  = mean(restBio(:,edgeIdx,1),2);
compWalk  = mean(walkBio(:,edgeIdx,1),2);
compDelta = compWalk-compRest;

compStats = adjusted_effect_summary(compDelta,M,Group);

fprintf('\n------------------------------------------------------------\n');
fprintf('32-EDGE THETA COMPONENT AVERAGE — DESCRIPTIVE/POST-SELECTION\n');
fprintf('------------------------------------------------------------\n');

fprintf('HC Rest mean   = %.6f\n',mean(compRest(vgHC)));
fprintf('HC Walk mean   = %.6f\n',mean(compWalk(vgHC)));
fprintf('HC Delta mean  = %.6f\n',mean(compDelta(vgHC)));

fprintf('PD Rest mean   = %.6f\n',mean(compRest(vgPD)));
fprintf('PD Walk mean   = %.6f\n',mean(compWalk(vgPD)));
fprintf('PD Delta mean  = %.6f\n',mean(compDelta(vgPD)));

fprintf('\nAdjusted Group beta (PD-HC) = %.6f db-wPLI\n',compStats.betaGroup);
fprintf('HC3 robust SE                = %.6f\n',compStats.seHC3);
fprintf('HC3 robust 95%% CI            = [%.6f, %.6f]\n', ...
    compStats.ciLow,compStats.ciHigh);
fprintf('Hedges-type adjusted g       = %.3f\n',compStats.gAdjusted);

fprintf(['\nCAUTION: because the 32 edges were selected by NBS in this same dataset,\n' ...
         'this component-average beta/CI/g is DESCRIPTIVE and may be\n' ...
         'optimistically biased by post-selection. It must not replace the\n' ...
         'joint-FWER NBS p-value as the primary inferential result.\n']);

%% =========================================================================
% 11. WHOLE-NETWORK BAND-AVERAGE EFFECTS
%
% Report all 3 biological bands together to avoid selectively reporting
% only the band that yielded the significant NBS component.
%% =========================================================================

bandNames = ["theta";"alpha";"beta"];

EffectType = strings(4,1);
Band = strings(4,1);

MeanDelta_HC = nan(4,1);
MeanDelta_PD = nan(4,1);

AdjustedBeta_PDminusHC = nan(4,1);
HC3_SE = nan(4,1);
CI95_Low = nan(4,1);
CI95_High = nan(4,1);
HedgesType_gAdjusted = nan(4,1);

% First row = selected component descriptive.
EffectType(1) = "NBS_component_average_POST_SELECTION";
Band(1) = "theta";

MeanDelta_HC(1) = mean(compDelta(vgHC));
MeanDelta_PD(1) = mean(compDelta(vgPD));

AdjustedBeta_PDminusHC(1) = compStats.betaGroup;
HC3_SE(1) = compStats.seHC3;
CI95_Low(1) = compStats.ciLow;
CI95_High(1) = compStats.ciHigh;
HedgesType_gAdjusted(1) = compStats.gAdjusted;

for b = 1:3

    yBand = mean(deltaNBS(:,:,b),2);

    st = adjusted_effect_summary(yBand,M,Group);

    rr = b+1;

    EffectType(rr) = "Whole_55_edge_network_average";
    Band(rr) = bandNames(b);

    MeanDelta_HC(rr) = mean(yBand(vgHC));
    MeanDelta_PD(rr) = mean(yBand(vgPD));

    AdjustedBeta_PDminusHC(rr) = st.betaGroup;
    HC3_SE(rr) = st.seHC3;
    CI95_Low(rr) = st.ciLow;
    CI95_High(rr) = st.ciHigh;
    HedgesType_gAdjusted(rr) = st.gAdjusted;

end

EffectSummary = table( ...
    EffectType, ...
    Band, ...
    MeanDelta_HC, ...
    MeanDelta_PD, ...
    AdjustedBeta_PDminusHC, ...
    HC3_SE, ...
    CI95_Low, ...
    CI95_High, ...
    HedgesType_gAdjusted);

fprintf('\n------------------------------------------------------------\n');
fprintf('EFFECT-SIZE SUMMARY\n');
fprintf('------------------------------------------------------------\n');

disp(EffectSummary);

writetable(EffectSummary,outFile);

fprintf('\nSaved:\n%s\n',outFile);

%% =========================================================================
% 12. FINAL VALIDATION STATEMENT
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('STEP 9B COMPLETE\n');
fprintf('============================================================\n');

fprintf(['The fast primary run is computationally plausible because the\n' ...
         'resampling stage operates only on a 117 x 165 statistics matrix,\n' ...
         'not on raw EEG or Fourier data. The exact-reproduction check above\n' ...
         'provides the more important validation of execution correctness.\n']);

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function st = adjusted_effect_summary(y,M,Group)

    y = double(y(:));

    N = numel(y);
    p = size(M,2);

    beta = M\y;

    residual = y-M*beta;

    XtXinv = pinv(M'*M);

    Hdiag = sum((M*XtXinv).*M,2);

    % HC3 sandwich covariance:
    uAdj2 = residual.^2 ./ (1-Hdiag).^2;

    meat = M'*(M.*uAdj2);

    Vhc3 = XtXinv*meat*XtXinv;

    se = sqrt(Vhc3(2,2));

    df = N-p;
    tcrit = student_t_critical_975(df);

    ciLow = beta(2)-tcrit*se;
    ciHigh = beta(2)+tcrit*se;

    % Hedges-type adjusted standardized contrast:
    % beta_group divided by full-model residual SD, then small-sample J.
    rmse = sqrt(sum(residual.^2)/df);

    dAdj = beta(2)/rmse;

    J = 1 - 3/(4*df-1);

    gAdj = J*dAdj;

    st = struct();

    st.betaGroup = beta(2);
    st.seHC3 = se;
    st.ciLow = ciLow;
    st.ciHigh = ciHigh;
    st.gAdjusted = gAdj;

end

%% =========================================================================

function tcrit = student_t_critical_975(df)

    % Exact positive t critical value for CDF=0.975 using betaincinv.
    % For t>0:
    % P(T>t) = 0.5*I_x(df/2,1/2),
    % x = df/(df+t^2).

    x = betaincinv(0.05,df/2,0.5);

    tcrit = sqrt(df*(1-x)/x);

end

%% =========================================================================

function v = compute_signed_v_all_only( ...
    Y,M,Mpinv,C,groupColumn, ...
    vgHC,vgPD,sumR_HC,sumR_PD,A_HC,A_PD)

    Bhat = Mpinv*Y;

    E = Y-M*Bhat;

    sseHC = sum(E(vgHC,:).^2,1);
    ssePD = sum(E(vgPD,:).^2,1);

    if any(sseHC<=0) || any(ssePD<=0)
        error('Nonpositive variance-group residual SSE.');
    end

    wHC = sumR_HC./sseHC;
    wPD = sumR_PD./ssePD;

    betaGroup = Bhat(groupColumn,:);

    K = size(Y,2);

    v = nan(1,K);

    for k = 1:K

        MtWM = wHC(k)*A_HC + wPD(k)*A_PD;

        q = C'*(MtWM\C);

        if ~isfinite(q) || q<=0
            error('Invalid contrast variance.');
        end

        v(k) = betaGroup(k)/sqrt(q);

    end

end

%% =========================================================================

function maxExtent = joint_max_component_extent( ...
    vVec,nEdges,nBands,threshold,edgeIndex1,edgeIndex2,nNodes)

    V = reshape(vVec,nEdges,nBands);

    maxExtent = 0;

    for b = 1:nBands

        extPos = max_component_extent( ...
            V(:,b)>threshold, ...
            edgeIndex1,edgeIndex2,nNodes);

        extNeg = max_component_extent( ...
            V(:,b)<-threshold, ...
            edgeIndex1,edgeIndex2,nNodes);

        maxExtent = max([maxExtent,extPos,extNeg]);

    end

end

%% =========================================================================

function maxExt = max_component_extent( ...
    edgeMask,edgeIndex1,edgeIndex2,nNodes)

    selectedEdges = find(edgeMask);

    if isempty(selectedEdges)
        maxExt = 0;
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

        maxExt = max(maxExt,sum(edgeThis));

    end

end

%% =========================================================================

function q = empirical_percentile(x,p)

    x = sort(double(x(:)));

    if isempty(x)
        q = NaN;
        return;
    end

    pos = 1+(numel(x)-1)*p;

    lo = floor(pos);
    hi = ceil(pos);

    if lo==hi
        q = x(lo);
    else
        w = pos-lo;
        q = (1-w)*x(lo)+w*x(hi);
    end

end

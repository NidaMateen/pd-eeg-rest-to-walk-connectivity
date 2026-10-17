%% =========================================================================
% PAPER 3 — STEP 7C
% VERIFY ALL MANUAL ICA DECISIONS BEFORE FINAL db-wPLI CONNECTIVITY
%
% File:
%   STEP7C_VERIFY_FINAL_ICA_DECISIONS.m
%
% PURPOSE
% -------
% This is a READ-ONLY preflight audit. It does NOT modify ICA decisions.
%
% It verifies that:
%   1) The full ICA master summary contains 130 subjects (23 HC, 107 PD).
%   2) The >20% manual-review trigger identifies exactly the flagged subjects.
%   3) Every flagged subject has an authoritative FINAL_ICA_reject.mat file.
%   4) Saved autoRejectIC matches the CURRENT full-batch automatic decision.
%   5) restoredIC is a subset of autoRejectIC.
%   6) finalRejectIC = autoRejectIC minus restoredIC.
%   7) All IC indices are valid for that subject's current ICA decomposition.
%   8) The saved manual decisions match the decisions finalized during review.
%   9) The cumulative manual-review audit contains all flagged subjects.
%
% If everything passes, the script prints:
%
%   >>> MANUAL ICA GATE: PASS <<<
%
% and it is safe to move to the final db-wPLI connectivity pipeline.
%% =========================================================================

clear;
clc;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();
addpath(cfg.fieldtripRoot);

%% =========================================================================
% 1. PATHS
%% =========================================================================

datasetRoot = ...
    cfg.datasetRoot;

icaRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_FULL130_ICA');

finalRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_FINAL_CONNECTIVITY');

reviewRoot = fullfile( ...
    datasetRoot, ...
    'PAPER3_MANUAL_ICA_REVIEW');

masterSummaryFile = fullfile( ...
    icaRoot, ...
    'Paper3_FULL130_ICA_summary.csv');

manualAuditFile = fullfile( ...
    reviewRoot, ...
    'Paper3_MANUAL_ICA_FINAL_DECISIONS.csv');

%% =========================================================================
% 2. EXPECTED MANUAL DECISIONS FROM THE COMPLETED REVIEW
%
% These are used only as an integrity check against accidental wrong entry.
% They do NOT change any file.
%% =========================================================================

expectedSubject = [ ...
    "sub-002"
    "sub-006"
    "sub-034"
    "sub-064"
    "sub-072"
    "sub-086"
    "sub-087"
    "sub-089"
    "sub-097"
    "sub-128"
    "sub-130"
    "sub-134"
    "sub-138"
    "sub-140"];

expectedRestore = { ...
    [1 3 5 14]    % sub-002
    []            % sub-006
    [1 2]         % sub-034
    []            % sub-064
    []            % sub-072
    [12]          % sub-086
    []            % sub-087
    [8 15]        % sub-089
    [3]           % sub-097
    []            % sub-128
    [1 2]         % sub-130
    []            % sub-134
    [3]           % sub-138
    []};          % sub-140

expectedMap = containers.Map('KeyType','char','ValueType','any');

for k = 1:numel(expectedSubject)
    expectedMap(char(expectedSubject(k))) = expectedRestore{k};
end

%% =========================================================================
% 3. BASIC FILE CHECKS
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('PAPER 3 — FINAL MANUAL ICA DECISION PREFLIGHT\n');
fprintf('============================================================\n');

if ~exist(masterSummaryFile,'file')
    error('Missing master ICA summary:\n%s',masterSummaryFile);
end

if ~exist(finalRoot,'dir')
    error('Missing final-decision root:\n%s',finalRoot);
end

%% =========================================================================
% 4. LOAD MASTER SUMMARY
%% =========================================================================

M = readtable( ...
    masterSummaryFile, ...
    'TextType','string');

requiredCols = ["Subject","Group","NIC","NRejectedIC","PctRejectedIC"];

for k = 1:numel(requiredCols)

    if ~ismember(requiredCols(k),string(M.Properties.VariableNames))
        error('Master summary is missing required column: %s',requiredCols(k));
    end

end

nAll = height(M);
nHC = sum(M.Group=="HC");
nPD = sum(M.Group=="PD");

fprintf('Master subjects = %d\n',nAll);
fprintf('HC              = %d\n',nHC);
fprintf('PD              = %d\n',nPD);

globalPass = true;

if nAll ~= 130 || nHC ~= 23 || nPD ~= 107

    warning('Expected 130 subjects = 23 HC + 107 PD.');
    globalPass = false;

end

%% =========================================================================
% 5. IDENTIFY >20%% FLAGGED SUBJECTS FROM MASTER SUMMARY
%% =========================================================================

F = M(M.PctRejectedIC > 20,:);

F = sortrows(F,'Subject');

flagged = string(F.Subject);

fprintf('\nFlagged >20%% subjects = %d\n',height(F));

expectedSorted = sort(expectedSubject);

if height(F) ~= numel(expectedSubject) || ...
   ~isequal(sort(flagged),expectedSorted)

    fprintf('\nERROR: Current >20%% flagged set does not match expected review set.\n');

    fprintf('\nCurrent flagged subjects:\n');
    disp(sort(flagged));

    fprintf('Expected flagged subjects:\n');
    disp(expectedSorted);

    globalPass = false;

else

    fprintf('Flagged-subject set matches the completed manual-review set.\n');

end

%% =========================================================================
% 6. VERIFY EACH FLAGGED SUBJECT
%% =========================================================================

N = height(F);

Subject       = strings(N,1);
Group         = strings(N,1);
AutoPct       = nan(N,1);
NIC           = nan(N,1);
NAutoReject   = nan(N,1);
RestoredIC    = strings(N,1);
NFinalReject  = nan(N,1);
FinalPct      = nan(N,1);
Status        = strings(N,1);

for i = 1:N

    subj = string(F.Subject(i));
    grp  = string(F.Group(i));

    Subject(i) = subj;
    Group(i) = grp;
    AutoPct(i) = F.PctRejectedIC(i);
    NIC(i) = F.NIC(i);

    subjectPass = true;
    reasons = strings(0,1);

    icaDir = fullfile( ...
        icaRoot, ...
        char(subj));

    autoFile = fullfile( ...
        icaDir, ...
        char(subj + "_automaticICA_decision.mat"));

    icaFile = fullfile( ...
        icaDir, ...
        char(subj + "_commonICA.mat"));

    finalFile = fullfile( ...
        finalRoot, ...
        char(subj), ...
        char(subj + "_FINAL_ICA_reject.mat"));

    %% --------------------------------------------------------------------
    % Required files
    %% --------------------------------------------------------------------

    if ~exist(autoFile,'file')
        subjectPass = false;
        reasons(end+1) = "missing automatic decision"; %#ok<SAGROW>
    end

    if ~exist(icaFile,'file')
        subjectPass = false;
        reasons(end+1) = "missing common ICA"; %#ok<SAGROW>
    end

    if ~exist(finalFile,'file')
        subjectPass = false;
        reasons(end+1) = "missing FINAL decision"; %#ok<SAGROW>
    end

    if ~subjectPass

        Status(i) = "FAIL: " + strjoin(reasons,"; ");
        globalPass = false;
        continue;

    end

    %% --------------------------------------------------------------------
    % Load CURRENT automatic decision and decomposition
    %% --------------------------------------------------------------------

    A = load(autoFile);
    I = load(icaFile);
    D = load(finalFile);

    if ~isfield(A,'rejectIC')

        subjectPass = false;
        reasons(end+1) = "auto file lacks rejectIC"; %#ok<SAGROW>
        autoReject = [];

    else

        autoReject = unique(double(A.rejectIC(:)'),'stable');

    end

    if ~isfield(I,'comp')

        subjectPass = false;
        reasons(end+1) = "ICA file lacks comp"; %#ok<SAGROW>
        currentNIC = NaN;

    else

        currentNIC = size(I.comp.topo,2);

    end

    %% --------------------------------------------------------------------
    % Required fields in FINAL file
    %% --------------------------------------------------------------------

    requiredFinalFields = { ...
        'finalRejectIC', ...
        'autoRejectIC', ...
        'restoredIC', ...
        'NIC'};

    for q = 1:numel(requiredFinalFields)

        if ~isfield(D,requiredFinalFields{q})

            subjectPass = false;
            reasons(end+1) = ...
                "FINAL lacks " + string(requiredFinalFields{q}); %#ok<SAGROW>

        end

    end

    if ~subjectPass

        Status(i) = "FAIL: " + strjoin(reasons,"; ");
        globalPass = false;
        continue;

    end

    savedAuto = unique(double(D.autoRejectIC(:)'),'stable');
    restored  = unique(double(D.restoredIC(:)'),'stable');
    finalRej  = unique(double(D.finalRejectIC(:)'),'stable');

    NAutoReject(i) = numel(autoReject);
    RestoredIC(i) = vecstr(restored);
    NFinalReject(i) = numel(finalRej);

    if isfinite(currentNIC)
        FinalPct(i) = 100*numel(finalRej)/currentNIC;
    end

    %% --------------------------------------------------------------------
    % A. Current automatic decision must equal saved autoRejectIC
    %% --------------------------------------------------------------------

    if ~isequal(autoReject,savedAuto)

        subjectPass = false;
        reasons(end+1) = "saved autoRejectIC != current auto reject"; %#ok<SAGROW>

    end

    %% --------------------------------------------------------------------
    % B. NIC consistency
    %% --------------------------------------------------------------------

    if ~isequal(double(D.NIC),double(currentNIC)) || ...
       ~isequal(double(F.NIC(i)),double(currentNIC))

        subjectPass = false;
        reasons(end+1) = "NIC mismatch"; %#ok<SAGROW>

    end

    %% --------------------------------------------------------------------
    % C. Valid IC indices
    %% --------------------------------------------------------------------

    allIndices = [savedAuto restored finalRej];

    if any(~isfinite(allIndices)) || ...
       any(mod(allIndices,1)~=0) || ...
       any(allIndices<1) || ...
       any(allIndices>currentNIC)

        subjectPass = false;
        reasons(end+1) = "invalid IC index"; %#ok<SAGROW>

    end

    %% --------------------------------------------------------------------
    % D. restoredIC must be subset of autoRejectIC
    %% --------------------------------------------------------------------

    if any(~ismember(restored,savedAuto))

        subjectPass = false;
        reasons(end+1) = "restored IC not in auto reject set"; %#ok<SAGROW>

    end

    %% --------------------------------------------------------------------
    % E. finalRejectIC must equal autoRejectIC minus restoredIC
    %% --------------------------------------------------------------------

    mathematicallyExpectedFinal = setdiff( ...
        savedAuto, ...
        restored, ...
        'stable');

    if ~isequal(finalRej,mathematicallyExpectedFinal)

        subjectPass = false;
        reasons(end+1) = "final reject != auto minus restored"; %#ok<SAGROW>

    end

    %% --------------------------------------------------------------------
    % F. Master NRejectedIC must match current automatic set
    %% --------------------------------------------------------------------

    if double(F.NRejectedIC(i)) ~= numel(autoReject)

        subjectPass = false;
        reasons(end+1) = "master NRejectedIC mismatch"; %#ok<SAGROW>

    end

    %% --------------------------------------------------------------------
    % G. Check against the decisions finalized during manual review
    %% --------------------------------------------------------------------

    key = char(subj);

    if ~isKey(expectedMap,key)

        subjectPass = false;
        reasons(end+1) = "subject absent from expected decision map"; %#ok<SAGROW>

    else

        tmpExpected = expectedMap(key);

        expectedRestored = unique( ...
            double(tmpExpected(:))', ...
            'stable');

        if ~isequal(restored,expectedRestored)

            subjectPass = false;
            reasons(end+1) = ...
                "saved restore " + vecstr(restored) + ...
                " != reviewed restore " + vecstr(expectedRestored); %#ok<SAGROW>

        end

    end

    %% --------------------------------------------------------------------
    % Final subject status
    %% --------------------------------------------------------------------

    if subjectPass

        Status(i) = "PASS";

    else

        Status(i) = "FAIL: " + strjoin(reasons,"; ");
        globalPass = false;

    end

end

%% =========================================================================
% 7. SHOW SUBJECT-LEVEL PREFLIGHT TABLE
%% =========================================================================

T = table( ...
    Subject, ...
    Group, ...
    NIC, ...
    AutoPct, ...
    NAutoReject, ...
    RestoredIC, ...
    NFinalReject, ...
    FinalPct, ...
    Status);

fprintf('\n============================================================\n');
fprintf('SUBJECT-LEVEL MANUAL ICA PREFLIGHT\n');
fprintf('============================================================\n');

disp(T);

%% =========================================================================
% 8. VERIFY CUMULATIVE MANUAL AUDIT FILE
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('CUMULATIVE MANUAL-REVIEW AUDIT\n');
fprintf('============================================================\n');

if ~exist(manualAuditFile,'file')

    warning('Manual audit CSV is missing:\n%s',manualAuditFile);
    globalPass = false;

else

    A = readtable( ...
        manualAuditFile, ...
        'TextType','string');

    if ~ismember("Subject",string(A.Properties.VariableNames))

        warning('Manual audit CSV lacks Subject column.');
        globalPass = false;

    else

        auditSubjects = unique(string(A.Subject));

        fprintf('Unique subjects in manual audit = %d\n',numel(auditSubjects));

        missingFromAudit = setdiff(expectedSubject,auditSubjects);
        unexpectedAudit = setdiff(auditSubjects,expectedSubject);

        if ~isempty(missingFromAudit)

            fprintf('Missing from manual audit:\n');
            disp(missingFromAudit);

            globalPass = false;

        end

        if ~isempty(unexpectedAudit)

            fprintf('Unexpected subjects in manual audit:\n');
            disp(unexpectedAudit);

            globalPass = false;

        end

        if isempty(missingFromAudit) && isempty(unexpectedAudit)

            fprintf('Manual audit contains exactly all 14 flagged subjects.\n');

        end

    end

end

%% =========================================================================
% 9. OPTIONAL FIELDTRIP PATH WARNING
%
% This does NOT alter the MATLAB path.
% It only reports what MATLAB would currently call.
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('CURRENT FIELDTRIP RESOLUTION — INFORMATION ONLY\n');
fprintf('============================================================\n');

ftComp = which('ft_componentanalysis');
ftConn = which('ft_connectivityanalysis');
ftDefs = which('ft_defaults');

fprintf('ft_defaults             : %s\n',ftDefs);
fprintf('ft_componentanalysis    : %s\n',ftComp);
fprintf('ft_connectivityanalysis : %s\n',ftConn);

if contains(string(ftComp),'fieldtrip-20240214','IgnoreCase',true) || ...
   contains(string(ftConn),'fieldtrip-20240214','IgnoreCase',true)

    fprintf(['\nWARNING: MATLAB currently resolves at least one FieldTrip function\n' ...
             'to the older fieldtrip-20240214 installation.\n' ...
             'Do NOT start final connectivity until the intended FieldTrip\n' ...
             'installation is the one resolved on the MATLAB path.\n']);

else

    fprintf('\nNo fieldtrip-20240214 path collision detected in these functions.\n');

end

%% =========================================================================
% 10. FINAL GATE
%% =========================================================================

fprintf('\n============================================================\n');

if globalPass

    fprintf('>>> MANUAL ICA GATE: PASS <<<\n');
    fprintf('All 14 flagged subjects have valid, internally consistent final decisions.\n');
    fprintf('The final manual ICA decision files are ready for connectivity.\n');

else

    fprintf('>>> MANUAL ICA GATE: FAIL <<<\n');
    fprintf('Do NOT start final db-wPLI connectivity yet.\n');
    fprintf('Resolve every FAIL line above first.\n');

end

fprintf('============================================================\n');

%% =========================================================================
%% LOCAL FUNCTION
%% =========================================================================

function s = vecstr(x)

    x = double(x(:)');

    if isempty(x)
        s = "[]";
    else
        s = string(mat2str(x));
    end

end

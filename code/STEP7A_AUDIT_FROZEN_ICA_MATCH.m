%% STEP7A_AUDIT_FROZEN_ICA_MATCH
% READ-ONLY audit of a fresh STEP7 ICA run against the original validated
% 14-subject ICA decision record (Paper3 Data QC dossier, 17 Sep 2026,
% Appendix A). Does not modify any data or ICA decisions.
%
% Run from MATLAB with the repository code/ folder on the MATLAB path and
% PAPER3_DATASET_ROOT + FIELDTRIP_ROOT set as documented in README.
% IMPORTANT: Equal component indices/rejection sets do not, by themselves,
% prove identical ICA spatial maps. Compare original and new decompositions
% before transferring manual decisions when original MAT files are available.

cfg = paper3_config();
root = fullfile(cfg.datasetRoot,'PAPER3_FULL130_ICA');
summaryPath = fullfile(root,'Paper3_FULL130_ICA_summary.csv');
assert(isfile(summaryPath),'Missing ICA summary: %s',summaryPath);
T = readtable(summaryPath,'TextType','string');
req = ["Subject","Group","NIC","NRejectedIC","PctRejectedIC"];
assert(all(ismember(req,string(T.Properties.VariableNames))), ...
    'Unexpected summary columns. Please inspect the master summary.');
assert(height(T)==130 && sum(T.Group=="HC")==23 && sum(T.Group=="PD")==107, ...
    'Cohort differs from expected 130 (23 HC, 107 PD).');

IDs = ["sub-002";"sub-006";"sub-034";"sub-064";"sub-072"; ...
    "sub-086";"sub-087";"sub-089";"sub-097";"sub-128"; ...
    "sub-130";"sub-134";"sub-138";"sub-140"];
NExpectedIC = [57;59;59;59;59;59;59;57;59;59;58;57;59;59];
% Frozen decisions from historical final rejected-component list
ExpectedFinal = { ...
    [7 8 9 10 11 12 18 20 28 37]; ...           % sub-002
    [1 3 4 8 10 11 14 17 19 23 26 30 36 46]; ... % sub-006
    [5 11 13 15 18 24 26 27 28 31 32 37 41 50]; ... % sub-034
    [5 6 7 10 13 15 18 20 23 26 28 29 35 50]; ... % sub-064
    [1 2 4 6 11 12 17 18 21 24 25 33 36 56]; ... % sub-072
    [1 2 3 4 5 7 15 16 26 33 36 46]; ...       % sub-086
    [1 2 3 6 10 12 17 20 24 27 32 33 35 38 55 56]; ... % sub-087
    [4 5 13 18 19 26 27 28 36 37 44 50]; ...   % sub-089
    [6 11 20 21 26 29 42 43 44 45 47 54]; ...  % sub-097
    [1 2 9 10 12 22 25 29 30 40 46 48 59]; ... % sub-128
    [5 8 9 15 17 18 19 26 31 34]; ...          % sub-130
    [1 2 4 9 10 11 12 21 24 29 33 37]; ...    % sub-134
    [1 4 5 10 12 16 19 20 21 23 24 26 28 31 36 39 42 49]; ... % sub-138
    [2 4 7 9 12 15 18 20 22 23 24 26 28 31 35 36 37 41 42 46 54 59]}; % sub-140
ExpectedRestore = { ...
    [1 3 5 14]; []; [1 2]; []; []; [12]; []; [8 15]; [3]; []; [1 2]; []; [3]; []};

flagged = sort(string(T.Subject(T.PctRejectedIC>20)));
assert(numel(flagged)==14 && isequal(flagged,sort(IDs)), ...
    'Flagged-subject IDs differ from the original manually reviewed set.');

status = strings(numel(IDs),1);
actualNIC = nan(numel(IDs),1);
actualAutoCount = nan(numel(IDs),1);
expectedAutoCount = nan(numel(IDs),1);

fprintf('\n*** READ-ONLY ICA FROZEN-DECISION AUDIT ***\n');
for k=1:numel(IDs)
    subj = IDs(k);
    idx = find(string(T.Subject)==subj);
    reasons = strings(0,1);
    if numel(idx)~=1
        status(k) = "FAIL: duplicate/missing subject in master";
        continue;
    end
    autoPath = fullfile(root,char(subj),char(subj+"_automaticICA_decision.mat"));
    compPath = fullfile(root,char(subj),char(subj+"_commonICA.mat"));
    if ~isfile(autoPath) || ~isfile(compPath)
        status(k) = "FAIL: missing auto decision or common ICA";
        continue;
    end
    A = load(autoPath,'rejectIC');
    I = load(compPath,'comp','nIC');
    if ~isfield(A,'rejectIC') || ~isfield(I,'comp') || ~isfield(I,'nIC')
        status(k) = "FAIL: missing required variables";
        continue;
    end
    observed = sort(unique(double(A.rejectIC(:)')));
    historical = sort(unique([ExpectedFinal{k},ExpectedRestore{k}]));
    nIC = size(I.comp.topo,2);
    actualNIC(k)=nIC;
    actualAutoCount(k)=numel(observed);
    expectedAutoCount(k)=numel(historical);

    if nIC ~= NExpectedIC(k) || nIC~=double(I.nIC) || nIC~=double(T.NIC(idx))
        reasons(end+1) = "NIC mismatch";
    end
    if numel(observed)~=double(T.NRejectedIC(idx))
        reasons(end+1) = "master rejected-count mismatch";
    end
    if ~isequal(observed,historical)
        reasons(end+1) = "automatic IC index-set differs";
        fprintf('%s actual automatic: %s\n',subj,mat2str(observed));
        fprintf('%s original expected: %s\n',subj,mat2str(historical));
    end
    if isempty(reasons)
        status(k)="PASS_INDEX_SET";
    else
        status(k)="FAIL: "+strjoin(reasons,'; ');
    end
end
Audit = table(IDs,actualNIC,NExpectedIC,actualAutoCount,expectedAutoCount,status, ...
    'VariableNames',{'Subject','NewNIC','OriginalNIC','NewAutoRejectCount', ...
    'OriginalAutoRejectCount','Status'});
disp(Audit);
if all(status=="PASS_INDEX_SET")
    fprintf('INDEX-SET GATE: PASS (14/14). No files were written.\n');
    fprintf('NEXT: verify original/new ICA component identity, then freeze final decisions.\n');
else
    fprintf('INDEX-SET GATE: FAIL. Do not apply frozen indices yet.\n');
end

%% ========================================================================
% PAPER 3 — STEP 6B
% EXACT SAMPLE SIZE FOR ALL 55 REGIONAL EDGES
%
% Uses already-generated regional coverage audit.
% NO EEG processing.
%
% FROZEN REGION USABILITY RULE:
%   NGood >= 2
%   AND
%   CoverageFraction >= 0.50
%
% A subject contributes to an edge only when BOTH regions forming
% that edge satisfy the regional usability rule.
%% ========================================================================

clear;
clc;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% ========================================================================
% 1. PATH
%% ========================================================================

auditDir = ...
fullfile(cfg.datasetRoot,'PAPER3_FULL130_BADCHANNEL_AUDIT');

regionFile = fullfile( ...
    auditDir, ...
    'Paper3_130_region_coverage_long.csv');

if ~exist(regionFile,'file')
    error('Region audit file not found:\n%s',regionFile);
end

R = readtable( ...
    regionFile, ...
    'TextType','string');

%% ========================================================================
% 2. FROZEN 11-REGION ORDER
%% ========================================================================

regions = [ ...
    "FrontalLeft"
    "FrontalMidline"
    "FrontalRight"
    "CentralLeft"
    "CentralRight"
    "ParietalLeft"
    "ParietalMidline"
    "ParietalRight"
    "TemporalLeft"
    "TemporalRight"
    "Occipital"
    ];

%% ========================================================================
% 3. APPLY FROZEN REGIONAL COVERAGE RULE
%% ========================================================================

R.RegionUsable = ...
    R.NGood >= 2 & ...
    R.CoverageFraction >= 0.50;

fprintf('\n============================================================\n');
fprintf('REGION USABILITY COUNTS\n');
fprintf('============================================================\n');

for r = 1:numel(regions)

    reg = regions(r);

    idx = R.Region == reg;

    fprintf('%-20s : %3d / %3d usable\n', ...
        reg, ...
        sum(R.RegionUsable(idx)), ...
        sum(idx));

end

%% ========================================================================
% 4. SUBJECT LIST
%% ========================================================================

subjects = unique(R.Subject,'stable');

Nsub = numel(subjects);

fprintf('\nSubjects found = %d\n',Nsub);

if Nsub ~= 130
    warning('Expected 130 subjects, found %d.',Nsub);
end

%% ========================================================================
% 5. CREATE SUBJECT x REGION USABILITY MATRIX
%% ========================================================================

usableMatrix = false(Nsub,numel(regions));

groupVector = strings(Nsub,1);

for s = 1:Nsub

    subj = subjects(s);

    rowsSub = R.Subject == subj;

    grp = unique(R.Group(rowsSub));

    if numel(grp) ~= 1
        error('Unexpected group assignment for %s.',subj);
    end

    groupVector(s) = grp;

    for r = 1:numel(regions)

        reg = regions(r);

        row = ...
            rowsSub & ...
            R.Region == reg;

        if sum(row) ~= 1
            error( ...
                'Expected one row for %s / %s, found %d.', ...
                subj,reg,sum(row));
        end

        usableMatrix(s,r) = ...
            R.RegionUsable(row);

    end

end

%% ========================================================================
% 6. CALCULATE ALL 55 EDGE SAMPLE SIZES
%% ========================================================================

Region1 = strings(0,1);
Region2 = strings(0,1);

N_Total = [];
N_HC = [];
N_PD = [];

Pct_Total = [];
Pct_HC = [];
Pct_PD = [];

HCmask = groupVector == "HC";
PDmask = groupVector == "PD";

NHCtotal = sum(HCmask);
NPDtotal = sum(PDmask);

for i = 1:numel(regions)-1

    for j = i+1:numel(regions)

        usableEdge = ...
            usableMatrix(:,i) & ...
            usableMatrix(:,j);

        Region1(end+1,1) = regions(i);
        Region2(end+1,1) = regions(j);

        nAll = sum(usableEdge);
        nHC  = sum(usableEdge & HCmask);
        nPD  = sum(usableEdge & PDmask);

        N_Total(end+1,1) = nAll;
        N_HC(end+1,1) = nHC;
        N_PD(end+1,1) = nPD;

        Pct_Total(end+1,1) = ...
            100*nAll/Nsub;

        Pct_HC(end+1,1) = ...
            100*nHC/NHCtotal;

        Pct_PD(end+1,1) = ...
            100*nPD/NPDtotal;

    end

end

EdgeSamples = table( ...
    Region1, ...
    Region2, ...
    N_Total, ...
    N_HC, ...
    N_PD, ...
    Pct_Total, ...
    Pct_HC, ...
    Pct_PD);

%% ========================================================================
% 7. SORT BY LOWEST AVAILABLE SAMPLE SIZE
%% ========================================================================

EdgeSamplesSorted = ...
    sortrows(EdgeSamples,'N_Total','ascend');

fprintf('\n============================================================\n');
fprintf('ALL 55 REGIONAL EDGES — SAMPLE SIZE\n');
fprintf('SORTED FROM LOWEST TO HIGHEST N\n');
fprintf('============================================================\n');

disp(EdgeSamplesSorted);

%% ========================================================================
% 8. CRITICAL PARIETAL MIDLINE ↔ OCCIPITAL EDGE
%% ========================================================================

critical = ...
    ((EdgeSamples.Region1 == "ParietalMidline") & ...
     (EdgeSamples.Region2 == "Occipital")) | ...
    ((EdgeSamples.Region2 == "ParietalMidline") & ...
     (EdgeSamples.Region1 == "Occipital"));

fprintf('\n============================================================\n');
fprintf('CRITICAL EDGE: PARIETAL MIDLINE ↔ OCCIPITAL\n');
fprintf('============================================================\n');

disp(EdgeSamples(critical,:));

%% ========================================================================
% 9. EDGES THAT DO NOT RETAIN ALL 130 SUBJECTS
%% ========================================================================

reduced = ...
    EdgeSamples.N_Total < Nsub;

ReducedEdges = ...
    sortrows( ...
        EdgeSamples(reduced,:), ...
        'N_Total', ...
        'ascend');

fprintf('\n============================================================\n');
fprintf('EDGES WITH N < 130\n');
fprintf('============================================================\n');

disp(ReducedEdges);

%% ========================================================================
% 10. OVERALL NETWORK SAMPLE-SIZE RANGE
%% ========================================================================

fprintf('\n============================================================\n');
fprintf('55-EDGE SAMPLE SIZE SUMMARY\n');
fprintf('============================================================\n');

fprintf('Minimum edge N = %d\n', ...
    min(EdgeSamples.N_Total));

fprintf('Maximum edge N = %d\n', ...
    max(EdgeSamples.N_Total));

fprintf('Median edge N  = %.1f\n', ...
    median(EdgeSamples.N_Total));

fprintf('\nEdges with N = 130 : %d / 55\n', ...
    sum(EdgeSamples.N_Total == 130));

fprintf('Edges with N < 130 : %d / 55\n', ...
    sum(EdgeSamples.N_Total < 130));

fprintf('Edges with N < 125 : %d / 55\n', ...
    sum(EdgeSamples.N_Total < 125));

fprintf('Edges with N < 120 : %d / 55\n', ...
    sum(EdgeSamples.N_Total < 120));

%% ========================================================================
% 11. WHICH SUBJECTS FAIL PARIETAL MIDLINE / OCCIPITAL
%% ========================================================================

idxPM = find(regions=="ParietalMidline");
idxOcc = find(regions=="Occipital");

PMusable = usableMatrix(:,idxPM);
Occusable = usableMatrix(:,idxOcc);

FailureType = strings(Nsub,1);

for s = 1:Nsub

    if PMusable(s) && Occusable(s)
        FailureType(s) = "PASS_BOTH";

    elseif ~PMusable(s) && Occusable(s)
        FailureType(s) = "FAIL_PAR_MID_ONLY";

    elseif PMusable(s) && ~Occusable(s)
        FailureType(s) = "FAIL_OCC_ONLY";

    else
        FailureType(s) = "FAIL_BOTH";
    end

end

CriticalSubjects = table( ...
    subjects, ...
    groupVector, ...
    PMusable, ...
    Occusable, ...
    FailureType, ...
    'VariableNames', ...
    {'Subject', ...
     'Group', ...
     'ParietalMidlineUsable', ...
     'OccipitalUsable', ...
     'Status'});

fprintf('\n============================================================\n');
fprintf('PARIETAL MIDLINE / OCCIPITAL FAILURE BREAKDOWN\n');
fprintf('============================================================\n');

disp(groupsummary(CriticalSubjects, ...
    {'Group','Status'}));

fprintf('\nSubjects failing one or both critical regions:\n');

disp( ...
    CriticalSubjects( ...
        CriticalSubjects.Status ~= "PASS_BOTH", ...
        :));

%% ========================================================================
% 12. SAVE RESULTS
%% ========================================================================

edgeFile = fullfile( ...
    auditDir, ...
    'Paper3_55edge_available_sample_sizes.csv');

writetable( ...
    EdgeSamplesSorted, ...
    edgeFile);

criticalFile = fullfile( ...
    auditDir, ...
    'Paper3_ParietalMidline_Occipital_subject_coverage.csv');

writetable( ...
    CriticalSubjects, ...
    criticalFile);

fprintf('\n============================================================\n');
fprintf('FILES SAVED\n');
fprintf('============================================================\n');

fprintf('\n55-edge sample sizes:\n%s\n',edgeFile);

fprintf('\nCritical-region subject table:\n%s\n',criticalFile);
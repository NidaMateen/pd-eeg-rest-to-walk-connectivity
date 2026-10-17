%% =========================================================================
% PAPER 3 — STEP 10R
% PROTOCOL-ALIGNED SECONDARY HEMISPHERIC ANALYSIS
%
% Exact implementation of the originally documented secondary protocol:
%
%   Subject-level NetworkType summary = MEDIAN Delta connectivity
%   Model, separately for theta/alpha/beta:
%
%   DeltaFC_summary ~ Group*NetworkType + Age_c + Sex + (1|Subject)
%
%   Group: HC / PD
%   NetworkType: IntraLeft / IntraRight / Inter
%   Subject-specific random intercept
%
% Main test:
%   omnibus Group x NetworkType interaction within each band
%
% Multiplicity:
%   Holm correction across exactly 3 band-level interaction P values
%
% IMPORTANT
% ---------
% This script exists to cleanly reproduce the FROZEN WRITTEN PROTOCOL.
% The earlier STEP10 mean-summary + multivariate-HC3 implementation should
% be retained only as an alternative/sensitivity implementation, not used
% as the protocol-concordant secondary analysis.
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

ALPHA = 0.05;

LEFT_REGIONS = [ ...
    "FrontalLeft"
    "CentralLeft"
    "ParietalLeft"
    "TemporalLeft"];

RIGHT_REGIONS = [ ...
    "FrontalRight"
    "CentralRight"
    "ParietalRight"
    "TemporalRight"];

BIO_BANDS = [ ...
    "theta"
    "alpha"
    "beta"];

NETWORK_TYPES = [ ...
    "IntraLeft"
    "IntraRight"
    "Inter"];

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
    'PAPER3_SECONDARY_HEMISPHERIC_PROTOCOL_ALIGNED');

if ~exist(outRoot,'dir')
    mkdir(outRoot);
end

if ~exist(statsInputFile,'file')
    error('Missing STEP 8 statistics input:\n%s',statsInputFile);
end

if ~exist(participantsFile,'file')
    error('Missing participants.tsv:\n%s',participantsFile);
end

if exist('fitlme','file') ~= 2
    error(['fitlme is unavailable. The Statistics and Machine Learning ' ...
           'Toolbox is required for the frozen LME implementation.']);
end

%% =========================================================================
% 3. LOAD CONNECTIVITY DATA
%% =========================================================================

S = load(statsInputFile);

requiredFields = { ...
    'subjects', ...
    'groups', ...
    'regionNames', ...
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

regionNames = string(S.regionNames(:));

edgeNames   = string(S.edgeNames(:));
edgeRegion1 = string(S.edgeRegion1(:));
edgeRegion2 = string(S.edgeRegion2(:));

bandNamesAll = lower(string(S.bandNames(:)));

deltaBio = double(S.deltaBio);

if size(deltaBio,2) ~= 55
    error('Expected 55 regional edges; found %d.',size(deltaBio,2));
end

%% =========================================================================
% 4. VALIDATE REQUIRED REGIONS / BANDS
%% =========================================================================

requiredRegions = [LEFT_REGIONS;RIGHT_REGIONS];

missingRegions = setdiff(requiredRegions,regionNames);

if ~isempty(missingRegions)
    fprintf('Missing required lateralized regions:\n');
    disp(missingRegions);
    error('Cannot construct frozen hemispheric analysis.');
end

bandIdx = nan(3,1);

for b = 1:3
    idx = find(bandNamesAll == BIO_BANDS(b));
    if numel(idx) ~= 1
        error('Could not uniquely identify band: %s',BIO_BANDS(b));
    end
    bandIdx(b) = idx;
end

%% =========================================================================
% 5. BUILD EXACT PRESPECIFIED EDGE SETS
%% =========================================================================

intraLeftEdges = all_pair_edges( ...
    LEFT_REGIONS,edgeRegion1,edgeRegion2);

intraRightEdges = all_pair_edges( ...
    RIGHT_REGIONS,edgeRegion1,edgeRegion2);

interEdges = cross_pair_edges( ...
    LEFT_REGIONS,RIGHT_REGIONS,edgeRegion1,edgeRegion2);

if numel(intraLeftEdges) ~= 6
    error('Expected 6 IntraLeft edges.');
end

if numel(intraRightEdges) ~= 6
    error('Expected 6 IntraRight edges.');
end

if numel(interEdges) ~= 16
    error('Expected 16 Inter edges.');
end

allSecondaryEdges = [ ...
    intraLeftEdges(:)
    intraRightEdges(:)
    interEdges(:)];

if numel(unique(allSecondaryEdges)) ~= 28
    error('Expected 28 unique lateralized edges.');
end

fprintf('\n============================================================\n');
fprintf('PROTOCOL-ALIGNED HEMISPHERIC EDGE SETS\n');
fprintf('============================================================\n');
fprintf('Intra-left  = %d edges\n',numel(intraLeftEdges));
fprintf('Intra-right = %d edges\n',numel(intraRightEdges));
fprintf('Inter       = %d edges\n',numel(interEdges));

%% =========================================================================
% 6. DEFINE SECONDARY COHORT
%
% All 28 lateralized edges must be finite for theta/alpha/beta.
%% =========================================================================

secondaryMask = true(size(deltaBio,1),1);

for b = 1:3
    tmp = deltaBio(:,allSecondaryEdges,bandIdx(b));
    secondaryMask = secondaryMask & all(isfinite(tmp),2);
end

subjects = subjectsAll(secondaryMask);
groups   = groupsAll(secondaryMask);
delta    = deltaBio(secondaryMask,:,:);

N   = numel(subjects);
nHC = sum(groups=="HC");
nPD = sum(groups=="PD");

fprintf('\n============================================================\n');
fprintf('SECONDARY HEMISPHERIC COHORT\n');
fprintf('============================================================\n');
fprintf('N total = %d\n',N);
fprintf('N HC    = %d\n',nHC);
fprintf('N PD    = %d\n',nPD);

if N ~= EXPECTED_N || nHC ~= EXPECTED_HC || nPD ~= EXPECTED_PD
    error('Expected secondary cohort 130 = 23 HC + 107 PD.');
end

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

if ~isequal(string(P.group(loc)),groups)
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

Age_c = age - mean(age);

fprintf('Age centering mean = %.6f years\n',mean(age));

%% =========================================================================
% 8. COMPUTE EXACT PROTOCOL MEDIAN NETWORK-TYPE SUMMARIES
%
% Dimensions:
%   networkSummary(subject, networkType, band)
%
% NetworkType:
%   1 = IntraLeft
%   2 = IntraRight
%   3 = Inter
%% =========================================================================

networkSummary = nan(N,3,3);

for b = 1:3

    bb = bandIdx(b);

    networkSummary(:,1,b) = ...
        median(delta(:,intraLeftEdges,bb),2,'omitnan');

    networkSummary(:,2,b) = ...
        median(delta(:,intraRightEdges,bb),2,'omitnan');

    networkSummary(:,3,b) = ...
        median(delta(:,interEdges,bb),2,'omitnan');

end

if any(~isfinite(networkSummary),'all')
    error('Unexpected nonfinite median NetworkType summaries.');
end

%% =========================================================================
% 9. SAVE LONG-FORM SUBJECT DATA USED BY THE LME
%% =========================================================================

nRows = N*3*3;

SubjectLong = strings(nRows,1);
GroupLong = strings(nRows,1);
SexLong = strings(nRows,1);
BandLong = strings(nRows,1);
NetworkTypeLong = strings(nRows,1);

Age_c_Long = nan(nRows,1);
DeltaFC_Median = nan(nRows,1);

row = 0;

for b = 1:3
    for nt = 1:3
        idx = row + (1:N);

        SubjectLong(idx) = subjects;
        GroupLong(idx) = groups;
        SexLong(idx) = sexRaw;

        BandLong(idx) = BIO_BANDS(b);
        NetworkTypeLong(idx) = NETWORK_TYPES(nt);

        Age_c_Long(idx) = Age_c;
        DeltaFC_Median(idx) = networkSummary(:,nt,b);

        row = row + N;
    end
end

LongTable = table( ...
    SubjectLong, ...
    GroupLong, ...
    SexLong, ...
    BandLong, ...
    NetworkTypeLong, ...
    Age_c_Long, ...
    DeltaFC_Median, ...
    'VariableNames', ...
    {'Subject','Group','Sex','Band','NetworkType','Age_c','DeltaFC_Median'});

longFile = fullfile( ...
    outRoot, ...
    'Paper3_SECONDARY_PROTOCOL_ALIGNED_subject_networktype_medians_long.csv');

writetable(LongTable,longFile);

%% =========================================================================
% 10. FIT THE FROZEN LME SEPARATELY BY BAND
%% =========================================================================

Band = BIO_BANDS;

N_Subjects = repmat(N,3,1);
N_HC = repmat(nHC,3,1);
N_PD = repmat(nPD,3,1);

F_GroupxNetworkType = nan(3,1);
DF1 = nan(3,1);
DF2 = nan(3,1);
P_Raw = nan(3,1);

lmeObjects = cell(3,1);
anovaTables = cell(3,1);
coefficientTables = cell(3,1);

for b = 1:3

    y = networkSummary(:,:,b);

    % Long format for this band only.
    Subject = repmat(subjects,3,1);
    Group   = repmat(groups,3,1);
    Sex     = repmat(sexRaw,3,1);
    AgeC    = repmat(Age_c,3,1);

    NetworkType = [ ...
        repmat("IntraLeft",N,1)
        repmat("IntraRight",N,1)
        repmat("Inter",N,1)];

    DeltaFC = [ ...
        y(:,1)
        y(:,2)
        y(:,3)];

    T = table( ...
        categorical(Subject), ...
        categorical(Group,["HC","PD"]), ...
        categorical(NetworkType,["IntraLeft","IntraRight","Inter"]), ...
        AgeC, ...
        categorical(Sex,["F","M"]), ...
        DeltaFC, ...
        'VariableNames', ...
        {'Subject','Group','NetworkType','Age_c','Sex','DeltaFC'});

    % Exact documented fixed/random-effects structure.
    %
    % REML is used for variance-component estimation. The fixed-effect
    % omnibus term is assessed from the fitted LME using Satterthwaite
    % denominator degrees of freedom.
    lme = fitlme( ...
        T, ...
        'DeltaFC ~ Group*NetworkType + Age_c + Sex + (1|Subject)', ...
        'FitMethod','REML');

    a = anova(lme,'DFMethod','satterthwaite');

    [fVal,df1Val,df2Val,pVal] = ...
        extract_interaction_test(a,"Group:NetworkType");

    F_GroupxNetworkType(b) = fVal;
    DF1(b) = df1Val;
    DF2(b) = df2Val;
    P_Raw(b) = pVal;

    lmeObjects{b} = lme;
    anovaTables{b} = a;
    coefficientTables{b} = lme.Coefficients;

    % Save exact model outputs per band.
    writetable( ...
        anova_output_to_table(a), ...
        fullfile(outRoot, ...
        sprintf('Paper3_SECONDARY_PROTOCOL_ALIGNED_%s_LME_ANOVA.csv', ...
        char(BIO_BANDS(b)))));

    writetable( ...
        generic_tabular_to_table(lme.Coefficients), ...
        fullfile(outRoot, ...
        sprintf('Paper3_SECONDARY_PROTOCOL_ALIGNED_%s_LME_coefficients.csv', ...
        char(BIO_BANDS(b)))));

end

P_Holm = holm_adjust(P_Raw);
Significant_Holm_0p05 = P_Holm < ALPHA;

InteractionTable = table( ...
    Band, ...
    N_Subjects, ...
    N_HC, ...
    N_PD, ...
    F_GroupxNetworkType, ...
    DF1, ...
    DF2, ...
    P_Raw, ...
    P_Holm, ...
    Significant_Holm_0p05);

interactionFile = fullfile( ...
    outRoot, ...
    'Paper3_SECONDARY_PROTOCOL_ALIGNED_GroupxNetworkType_LME_Holm.csv');

writetable(InteractionTable,interactionFile);

%% =========================================================================
% 11. DESCRIPTIVE GROUP SUMMARIES OF SUBJECT-LEVEL MEDIANS
%
% These are descriptive only. No extra inferential family is added.
%% =========================================================================

nDesc = 3*3;

BandD = strings(nDesc,1);
NetworkTypeD = strings(nDesc,1);

HC_Mean = nan(nDesc,1);
HC_SD = nan(nDesc,1);
HC_Median = nan(nDesc,1);

PD_Mean = nan(nDesc,1);
PD_SD = nan(nDesc,1);
PD_Median = nan(nDesc,1);

RawMeanDifference_PDminusHC = nan(nDesc,1);

row = 0;

for b = 1:3
    for nt = 1:3

        row = row+1;

        vals = networkSummary(:,nt,b);

        vHC = vals(groups=="HC");
        vPD = vals(groups=="PD");

        BandD(row) = BIO_BANDS(b);
        NetworkTypeD(row) = NETWORK_TYPES(nt);

        HC_Mean(row) = mean(vHC);
        HC_SD(row) = std(vHC,0);
        HC_Median(row) = median(vHC);

        PD_Mean(row) = mean(vPD);
        PD_SD(row) = std(vPD,0);
        PD_Median(row) = median(vPD);

        RawMeanDifference_PDminusHC(row) = ...
            PD_Mean(row)-HC_Mean(row);

    end
end

DescriptiveTable = table( ...
    BandD, ...
    NetworkTypeD, ...
    HC_Mean, ...
    HC_SD, ...
    HC_Median, ...
    PD_Mean, ...
    PD_SD, ...
    PD_Median, ...
    RawMeanDifference_PDminusHC, ...
    'VariableNames', ...
    {'Band','NetworkType', ...
     'HC_MeanOfSubjectMedians','HC_SD','HC_MedianOfSubjectMedians', ...
     'PD_MeanOfSubjectMedians','PD_SD','PD_MedianOfSubjectMedians', ...
     'RawMeanDifference_PDminusHC'});

descFile = fullfile( ...
    outRoot, ...
    'Paper3_SECONDARY_PROTOCOL_ALIGNED_descriptive_networktype_medians.csv');

writetable(DescriptiveTable,descFile);

%% =========================================================================
% 12. SAVE MAT + LOG
%% =========================================================================

resultFile = fullfile( ...
    outRoot, ...
    'Paper3_SECONDARY_PROTOCOL_ALIGNED_RESULTS.mat');

save( ...
    resultFile, ...
    'subjects','groups','age','sexRaw','Age_c', ...
    'intraLeftEdges','intraRightEdges','interEdges', ...
    'networkSummary','InteractionTable','DescriptiveTable', ...
    'lmeObjects','anovaTables','coefficientTables', ...
    '-v7.3');

logFile = fullfile( ...
    outRoot, ...
    'Paper3_SECONDARY_PROTOCOL_ALIGNED_analysis_log.txt');

fid = fopen(logFile,'w');

if fid==-1
    error('Could not create analysis log.');
end

fprintf(fid,'PAPER 3 — PROTOCOL-ALIGNED HEMISPHERIC SECONDARY ANALYSIS\n\n');
fprintf(fid,'Subject summary: MEDIAN Delta connectivity within each NetworkType.\n');
fprintf(fid,'Model: DeltaFC ~ Group*NetworkType + Age_c + Sex + (1|Subject)\n');
fprintf(fid,'Fit: REML\n');
fprintf(fid,'Denominator DF: Satterthwaite\n');
fprintf(fid,'Multiplicity: Holm across theta/alpha/beta interaction tests.\n');
fprintf(fid,'Cohort: %d = %d HC + %d PD.\n\n',N,nHC,nPD);

for b = 1:3
    fprintf(fid, ...
        '%s: F=%.8g, df1=%.8g, df2=%.8g, raw p=%.8g, Holm p=%.8g, sig=%d\n', ...
        BIO_BANDS(b), ...
        F_GroupxNetworkType(b), ...
        DF1(b), ...
        DF2(b), ...
        P_Raw(b), ...
        P_Holm(b), ...
        Significant_Holm_0p05(b));
end

fprintf(fid,'\nThis implementation supersedes STEP10 for protocol provenance.\n');
fprintf(fid,['The earlier mean-summary + multivariate-HC3 result may be retained ' ...
             'as an alternative robustness analysis only.\n']);

fclose(fid);

%% =========================================================================
% 13. CONSOLE REPORT
%% =========================================================================

fprintf('\n============================================================\n');
fprintf('PROTOCOL-ALIGNED STEP 10R COMPLETE\n');
fprintf('============================================================\n');

fprintf('\nGroup x NetworkType omnibus LME tests:\n');
disp(InteractionTable);

fprintf('\nDescriptive NetworkType summaries based on SUBJECT MEDIANS:\n');
disp(DescriptiveTable);

fprintf('\nSaved:\n');
fprintf('Interaction table : %s\n',interactionFile);
fprintf('Descriptives      : %s\n',descFile);
fprintf('Long data         : %s\n',longFile);
fprintf('Full MAT          : %s\n',resultFile);
fprintf('Analysis log      : %s\n',logFile);

fprintf('\nUse this result, not the original STEP10 result, for the manuscript\n');
fprintf('secondary hemispheric analysis.\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function edgeIdx = all_pair_edges(regions,edgeRegion1,edgeRegion2)

    edgeIdx = zeros(nchoosek(numel(regions),2),1);

    k = 0;

    for i = 1:numel(regions)-1
        for j = i+1:numel(regions)

            k = k+1;

            edgeIdx(k) = find_unique_edge( ...
                regions(i),regions(j),edgeRegion1,edgeRegion2);

        end
    end

end

function edgeIdx = cross_pair_edges(leftRegions,rightRegions,edgeRegion1,edgeRegion2)

    edgeIdx = zeros(numel(leftRegions)*numel(rightRegions),1);

    k = 0;

    for i = 1:numel(leftRegions)
        for j = 1:numel(rightRegions)

            k = k+1;

            edgeIdx(k) = find_unique_edge( ...
                leftRegions(i),rightRegions(j),edgeRegion1,edgeRegion2);

        end
    end

end

function idx = find_unique_edge(a,b,edgeRegion1,edgeRegion2)

    m = ...
        (edgeRegion1==a & edgeRegion2==b) | ...
        (edgeRegion1==b & edgeRegion2==a);

    idx = find(m);

    if numel(idx) ~= 1
        error('Could not uniquely identify edge %s -- %s.',a,b);
    end

end

function [Fval,df1,df2,pval] = extract_interaction_test(a,targetTerm)

    terms = get_anova_term_labels(a);

    targetCanonical = canonical_term(targetTerm);

    canonicalTerms = strings(numel(terms),1);

    for ii = 1:numel(terms)
        canonicalTerms(ii) = canonical_term(terms(ii));
    end

    idx = find(canonicalTerms==targetCanonical);

    if numel(idx) ~= 1
        fprintf('ANOVA terms returned by MATLAB:\n');
        disp(terms);
        error('Could not uniquely locate Group x NetworkType interaction.');
    end

    Fval = get_numeric_anova_column(a,idx,["FStat","F"]);
    df1  = get_numeric_anova_column(a,idx,["DF1","NumDF"]);
    df2  = get_numeric_anova_column(a,idx,["DF2","DenDF"]);
    pval = get_numeric_anova_column(a,idx,["pValue","pValueF","PValue"]);

end

function terms = get_anova_term_labels(a)

    % Modern MATLAB table
    if istable(a)

        if ~isempty(a.Properties.RowNames)
            terms = string(a.Properties.RowNames);
            return;
        end

        vars = string(a.Properties.VariableNames);

        hit = find(strcmpi(vars,"Term"),1);

        if ~isempty(hit)
            terms = string(a.(char(vars(hit))));
            return;
        end

        error('Cannot identify ANOVA term labels from table output.');

    end

    % Legacy Statistics Toolbox dataset array
    if isa(a,'dataset')

        if isprop(a.Properties,'ObsNames') && ...
           ~isempty(a.Properties.ObsNames)

            terms = string(a.Properties.ObsNames);
            return;
        end

        vars = string(a.Properties.VarNames);

        hit = find(strcmpi(vars,"Term"),1);

        if ~isempty(hit)
            terms = string(a.(char(vars(hit))));
            return;
        end

        error('Cannot identify ANOVA term labels from dataset output.');

    end

    error('Unsupported ANOVA output class: %s',class(a));

end

function c = canonical_term(s)

    c = lower(string(s));

    c = erase(c," ");
    c = erase(c,"(");
    c = erase(c,")");

    % Normalize common interaction separators / orderings.
    c = replace(c,"*","");

    if c=="networktype:group"
        c = "group:networktype";
    end

end

function vars = get_tabular_var_names(a)

    if istable(a)
        vars = string(a.Properties.VariableNames);
        return;
    end

    if isa(a,'dataset')
        vars = string(a.Properties.VarNames);
        return;
    end

    error('Unsupported tabular class: %s',class(a));

end

function x = get_numeric_anova_column(a,row,candidates)

    vars = get_tabular_var_names(a);

    idx = [];

    for k = 1:numel(candidates)

        hit = find(strcmpi(vars,candidates(k)),1);

        if ~isempty(hit)
            idx = hit;
            break;
        end

    end

    if isempty(idx)
        fprintf('ANOVA columns returned by MATLAB:\n');
        disp(vars);
        error('Required ANOVA column not found.');
    end

    col = a.(char(vars(idx)));

    if iscell(col)
        val = col{row};
    else
        val = col(row);
    end

    x = double(val);

end

function Tout = anova_output_to_table(a)

    terms = get_anova_term_labels(a);

    Tout = generic_tabular_to_table(a);

    % Add term labels explicitly for reproducible CSV output.
    if ismember('Term',Tout.Properties.VariableNames)
        Tout.Term = string(Tout.Term);
    else
        Tout = addvars( ...
            Tout, ...
            terms(:), ...
            'Before',1, ...
            'NewVariableNames','Term');
    end

end

function Tout = generic_tabular_to_table(a)

    if istable(a)
        Tout = a;
        return;
    end

    if ~isa(a,'dataset')
        error('Unsupported tabular class: %s',class(a));
    end

    vars = string(a.Properties.VarNames);

    Tout = table();

    for k = 1:numel(vars)

        name = char(vars(k));
        col = a.(name);

        % Preserve ordinary numeric/logical/categorical/string data.
        % Convert char matrices/cellstr only where needed for CSV safety.
        if ischar(col)
            col = string(cellstr(col));
        elseif iscellstr(col)
            col = string(col);
        end

        Tout.(name) = col;

    end

end

function pAdj = holm_adjust(p)

    p = double(p(:));
    m = numel(p);

    if any(~isfinite(p)) || any(p<0) || any(p>1)
        error('Invalid p-values for Holm adjustment.');
    end

    [pSort,ord] = sort(p,'ascend');

    adjSort = nan(m,1);
    runningMax = 0;

    for i = 1:m

        candidate = (m-i+1)*pSort(i);
        runningMax = max(runningMax,candidate);
        adjSort(i) = min(1,runningMax);

    end

    pAdj = nan(m,1);
    pAdj(ord) = adjSort;

end

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

%% PAPER 3 — STEP 7B
% Manual ICA review for subjects with >20% automatic IC rejection.
% >20% is a REVIEW TRIGGER, not an exclusion criterion or quota.
%
% Review ONLY automatically rejected ICs. You may restore false positives.
% Do NOT add new rejected ICs in this step.
%
% Frozen automatic rules shown for reference:
% Eye: max |corr(IC,EOG)| across Rest/Walk >= 0.40
% Muscle/electrode: HF ratio P30-45/P4-30 >= 0.60 in either condition
%                   AND topo focality max(abs(topo))/RMS(topo) >= 4.0
%
% Authoritative output used by STEP 8:
%   PAPER3_FINAL_CONNECTIVITY\sub-XXX\sub-XXX_FINAL_ICA_reject.mat
%
% Diagnostic figures/audit:
%   PAPER3_MANUAL_ICA_REVIEW\...

clear;
clc;
close all;


% Repository path configuration (machine-specific paths are not stored in this file).
cfg = paper3_config();

%% 1) PATHS / SETTINGS

datasetRoot = cfg.datasetRoot;
icaRoot = fullfile(datasetRoot,'PAPER3_FULL130_ICA');
connectivityRoot = fullfile(datasetRoot,'PAPER3_FINAL_CONNECTIVITY');
reviewRoot = fullfile(datasetRoot,'PAPER3_MANUAL_ICA_REVIEW');
masterSummaryFile = fullfile(icaRoot,'Paper3_FULL130_ICA_summary.csv');

REVIEW_EXISTING = false;   % true only if you intentionally want to redo a decision
IC_PER_PAGE = 4;
MAX_PLOT_POINTS = 4000;
ZOOM_SECONDS = 20;

EYE_CORR_THRESHOLD = 0.40;
HF_RATIO_THRESHOLD = 0.60;
FOCALITY_THRESHOLD = 4.00;

if ~exist(masterSummaryFile,'file')
    error('Master ICA summary not found:\n%s',masterSummaryFile);
end
if ~exist(connectivityRoot,'dir'), mkdir(connectivityRoot); end
if ~exist(reviewRoot,'dir'), mkdir(reviewRoot); end

%% 2) LOAD THE EXACT >20% REVIEW SET

M = readtable(masterSummaryFile,'TextType','string');
needVars = ["Subject","Group","PctRejectedIC"];
for k = 1:numel(needVars)
    if ~ismember(needVars(k),string(M.Properties.VariableNames))
        error('Missing master-summary column: %s',needVars(k));
    end
end

R = M(M.PctRejectedIC > 20,:);
R = sortrows(R,'PctRejectedIC','descend');

fprintf('\n============================================================\n');
fprintf('PAPER 3 MANUAL ICA REVIEW\n');
fprintf('============================================================\n');
fprintf('Flagged subjects = %d / %d\n',height(R),height(M));
fprintf('Criterion = automatic rejection >20%% (REVIEW TRIGGER ONLY)\n\n');

showVars = intersect({'Subject','Group','NBadChannels','NIC','NEyeIC', ...
    'NMuscleIC','NRejectedIC','PctRejectedIC','RejectedICLabels'}, ...
    R.Properties.VariableNames,'stable');
disp(R(:,showVars));

auditFile = fullfile(reviewRoot,'Paper3_MANUAL_ICA_FINAL_DECISIONS.csv');

%% 3) REVIEW ONE SUBJECT AT A TIME

for ss = 1:height(R)

    subj = string(R.Subject(ss));
    grp  = string(R.Group(ss));

    icaDir = fullfile(icaRoot,char(subj));
    subjConnDir = fullfile(connectivityRoot,char(subj));
    subjReviewDir = fullfile(reviewRoot,char(subj));
    if ~exist(subjConnDir,'dir'), mkdir(subjConnDir); end
    if ~exist(subjReviewDir,'dir'), mkdir(subjReviewDir); end

    finalDecisionFile = fullfile(subjConnDir,char(subj + "_FINAL_ICA_reject.mat"));
    if exist(finalDecisionFile,'file') && ~REVIEW_EXISTING
        fprintf('\n[SKIP] %s already finalized.\n',subj);
        continue;
    end

    icaFile  = fullfile(icaDir,char(subj + "_commonICA.mat"));
    autoFile = fullfile(icaDir,char(subj + "_automaticICA_decision.mat"));

    if ~exist(icaFile,'file'),  error('Missing %s',icaFile); end
    if ~exist(autoFile,'file'), error('Missing %s',autoFile); end

    I = load(icaFile);
    A = load(autoFile);
    if ~isfield(I,'comp'), error('%s: comp missing.',subj); end
    if ~isfield(A,'rejectIC'), error('%s: rejectIC missing.',subj); end

    comp = I.comp;
    autoRejectIC = unique(A.rejectIC(:)','stable');
    nIC = size(comp.topo,2);

    if numel(comp.trial) < 2
        error('%s: expected common ICA with Rest and Walk trials.',subj);
    end

    fs = get_sampling_rate(comp);
    restAct = double(comp.trial{1});
    walkAct = double(comp.trial{2});

    %% 3A) Recompute display metrics from CURRENT full-batch ICA
    focality = nan(nIC,1);
    hfRest = nan(nIC,1);
    hfWalk = nan(nIC,1);

    for ic = 1:nIC
        topo = double(comp.topo(:,ic));
        topoRMS = sqrt(mean(topo.^2,'omitnan'));
        if topoRMS > 0
            focality(ic) = max(abs(topo))/topoRMS;
        end
        hfRest(ic) = component_hf_ratio(restAct(ic,:),fs);
        hfWalk(ic) = component_hf_ratio(walkAct(ic,:),fs);
    end

    % Recover saved eye/muscle labels if the MAT contains recognizable fields.
    eyeIC = find_saved_ic_list(A,nIC,["eye"]);
    muscleIC = find_saved_ic_list(A,nIC,["muscle"]);
    eogMetric = find_saved_metric(A,nIC,["eog","corr"]);

    [xy,coordLabels] = load_subject_xy(datasetRoot,subj,string(comp.topolabel));

    fprintf('\n\n============================================================\n');
    fprintf('MANUAL REVIEW: %s | %s\n',subj,grp);
    fprintf('============================================================\n');
    fprintf('NIC = %d\n',nIC);
    fprintf('Automatic rejection = %d / %d (%.2f%%)\n', ...
        numel(autoRejectIC),nIC,100*numel(autoRejectIC)/nIC);
    fprintf('Auto rejected ICs: %s\n',mat2str(autoRejectIC));
    fprintf('\nRULE: enter ONLY ICs that should be RESTORED.\n');
    fprintf('Do not add new rejected ICs. Do not target a desired percentage.\n\n');

    reasons = strings(numel(autoRejectIC),1);
    for q = 1:numel(autoRejectIC)
        reasons(q) = reason_text(autoRejectIC(q),eyeIC,muscleIC);
    end

    D = table(autoRejectIC(:),reasons,hfRest(autoRejectIC),hfWalk(autoRejectIC), ...
        focality(autoRejectIC),eogMetric(autoRejectIC), ...
        'VariableNames',{'IC','AutoReason','HF_Rest','HF_Walk','Focality','MaxAbsEOGCorr'});
    disp(D);
    writetable(D,fullfile(subjReviewDir,char(subj + "_REJECTED_IC_DIAGNOSTICS.csv")));

    %% 3B) Diagnostic pages
    nPages = ceil(numel(autoRejectIC)/IC_PER_PAGE);

    for pg = 1:nPages
        firstIdx = (pg-1)*IC_PER_PAGE + 1;
        lastIdx = min(pg*IC_PER_PAGE,numel(autoRejectIC));
        pageIC = autoRejectIC(firstIdx:lastIdx);

        fig = figure('Color','w','Units','normalized','Position',[0.01 0.04 0.98 0.90], ...
            'Name',sprintf('%s page %d/%d',subj,pg,nPages));
        tl = tiledlayout(numel(pageIC),5,'TileSpacing','compact','Padding','compact');
        title(tl,sprintf('%s | %s | automatically rejected ICs | page %d/%d', ...
            subj,grp,pg,nPages),'Interpreter','none','FontWeight','bold');

        for rr = 1:numel(pageIC)
            ic = pageIC(rr);

            nexttile;
            plot_component_scalp(xy,coordLabels,string(comp.topolabel),double(comp.topo(:,ic)));
            title(sprintf('IC %d | %s',ic,reason_text(ic,eyeIC,muscleIC)),'Interpreter','none');

            nexttile;
            plot_full_trace(restAct(ic,:),fs,MAX_PLOT_POINTS);
            title('Rest activation'); xlabel('Time (s)'); ylabel('a.u.');

            nexttile;
            plot_full_trace(walkAct(ic,:),fs,MAX_PLOT_POINTS);
            title('Walk activation'); xlabel('Time (s)'); ylabel('a.u.');

            nexttile;
            plot_peak_zoom(restAct(ic,:),walkAct(ic,:),fs,ZOOM_SECONDS);
            title(sprintf('%g-s strongest-event zoom',ZOOM_SECONDS));
            xlabel('Relative time (s)'); ylabel('robust z');

            nexttile;
            plot_component_psd(restAct(ic,:),walkAct(ic,:),fs);
            if isfinite(eogMetric(ic))
                eogTxt = sprintf('EOG=%.2f',eogMetric(ic));
            else
                eogTxt = 'EOG=n/a';
            end
            title(sprintf('HF R/W=%.2f/%.2f | Foc=%.2f | %s', ...
                hfRest(ic),hfWalk(ic),focality(ic),eogTxt),'Interpreter','none');
        end

        pageFile = fullfile(subjReviewDir,sprintf('%s_manualICA_page_%02d.png',subj,pg));
        exportgraphics(fig,pageFile,'Resolution',160);

        fprintf('\n%s page %d/%d displayed: ICs %s\n',subj,pg,nPages,mat2str(pageIC));
        fprintf('Inspect the figure. Press ENTER when ready for the next page.\n');
        input('','s');
        close(fig);
    end

    %% 3C) Manual decision
    confirmed = false;
    while ~confirmed
        fprintf('\nAUTO REJECT = %s\n',mat2str(autoRejectIC));
        txt = input(['Enter ICs to RESTORE (example [11 42]); [] = restore none; ' ...
                     'q = stop and resume later:\n> '],'s');

        if strcmpi(strtrim(txt),'q')
            fprintf('Stopped before finalizing %s. Rerun later to resume.\n',subj);
            rebuild_manual_audit(R,connectivityRoot,auditFile);
            return;
        end
        if isempty(strtrim(txt)), txt = '[]'; end

        restoreIC = str2num(txt); %#ok<ST2NM>
        if isempty(restoreIC) && ~strcmp(strtrim(txt),'[]')
            fprintf('Could not parse entry. Try again.\n');
            continue;
        end
        restoreIC = unique(restoreIC(:)','stable');

        if any(~ismember(restoreIC,autoRejectIC))
            fprintf('ERROR: restore list must be a subset of AUTO REJECT.\n');
            continue;
        end

        finalRejectIC = setdiff(autoRejectIC,restoreIC,'stable');

        fprintf('\n------------------------------------------------------------\n');
        fprintf('SUBJECT      : %s\n',subj);
        fprintf('AUTO REJECT  : %s\n',mat2str(autoRejectIC));
        fprintf('RESTORE      : %s\n',mat2str(restoreIC));
        fprintf('FINAL REJECT : %s\n',mat2str(finalRejectIC));
        fprintf('FINAL %%      : %.2f%%\n',100*numel(finalRejectIC)/nIC);
        fprintf('------------------------------------------------------------\n');

        answ = input('Confirm FINAL decision? [y/n]: ','s');
        confirmed = strcmpi(strtrim(answ),'y');
    end

    reviewerNote = string(input('Optional note (ENTER for none):\n> ','s'));

    %% 3D) Save authoritative decision
    decisionSource = "MANUAL_REVIEW_GT20_TRIGGER";
    manualFlag = true;
    reviewDate = string(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));

    S = struct();
    S.Subject = subj;
    S.Group = grp;
    S.finalRejectIC = finalRejectIC;
    S.autoRejectIC = autoRejectIC;
    S.restoredIC = restoreIC;
    S.decisionSource = decisionSource;
    S.manualFlag = manualFlag;
    S.reviewDate = reviewDate;
    S.reviewerNote = reviewerNote;
    S.NIC = nIC;
    S.NAutoRejected = numel(autoRejectIC);
    S.NRestored = numel(restoreIC);
    S.NFinalRejected = numel(finalRejectIC);
    S.PctAutoRejected = 100*numel(autoRejectIC)/nIC;
    S.PctFinalRejected = 100*numel(finalRejectIC)/nIC;
    S.eyeIC_fromSavedDecision = eyeIC;
    S.muscleIC_fromSavedDecision = muscleIC;
    S.HFRatioRest = hfRest;
    S.HFRatioWalk = hfWalk;
    S.TopoFocality = focality;
    S.MaxAbsEOGCorr_fromSavedDecision = eogMetric;
    S.frozenThresholds = struct('EyeCorr',EYE_CORR_THRESHOLD, ...
        'HFRatio',HF_RATIO_THRESHOLD,'Focality',FOCALITY_THRESHOLD);

    atomic_save_mat(finalDecisionFile,S);

    % Text audit
    txtFile = fullfile(subjReviewDir,char(subj + "_MANUAL_ICA_DECISION.txt"));
    fid = fopen(txtFile,'w');
    fprintf(fid,'Subject: %s\n',subj);
    fprintf(fid,'Group: %s\n',grp);
    fprintf(fid,'Review date: %s\n',reviewDate);
    fprintf(fid,'NIC: %d\n',nIC);
    fprintf(fid,'Automatic rejection: %s\n',mat2str(autoRejectIC));
    fprintf(fid,'Restored: %s\n',mat2str(restoreIC));
    fprintf(fid,'FINAL rejection: %s\n',mat2str(finalRejectIC));
    fprintf(fid,'Automatic rejection %%: %.4f\n',100*numel(autoRejectIC)/nIC);
    fprintf(fid,'Final rejection %%: %.4f\n',100*numel(finalRejectIC)/nIC);
    fprintf(fid,'Reviewer note: %s\n',reviewerNote);
    fclose(fid);

    rebuild_manual_audit(R,connectivityRoot,auditFile);

    fprintf('\n[SAVED] %s final ICA decision:\n%s\n',subj,finalDecisionFile);
end

%% 4) FINAL STATUS
fprintf('\n============================================================\n');
fprintf('MANUAL ICA REVIEW STATUS\n');
fprintf('============================================================\n');

nFinal = 0;
remaining = strings(0,1);
for ss = 1:height(R)
    subj = string(R.Subject(ss));
    f = fullfile(connectivityRoot,char(subj),char(subj + "_FINAL_ICA_reject.mat"));
    if exist(f,'file')
        nFinal = nFinal + 1;
    else
        remaining(end+1,1) = subj; %#ok<SAGROW>
    end
end

fprintf('Finalized = %d / %d flagged subjects\n',nFinal,height(R));
if isempty(remaining)
    fprintf('ALL >20%% SUBJECTS ARE FINALIZED. Manual ICA gate CLOSED.\n');
else
    fprintf('Still requiring review:\n');
    disp(remaining);
end

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function fs = get_sampling_rate(comp)
    if isfield(comp,'fsample') && isfinite(comp.fsample)
        fs = double(comp.fsample); return;
    end
    if isfield(comp,'time') && ~isempty(comp.time) && numel(comp.time{1})>=2
        fs = 1/median(diff(double(comp.time{1}))); return;
    end
    error('Could not determine component sampling rate.');
end

function ratio = component_hf_ratio(x,fs)
    x = double(x(:)); x = x(isfinite(x));
    if numel(x) < round(4*fs), ratio = NaN; return; end
    win = min(round(2*fs),numel(x));
    if mod(win,2)==1, win = win-1; end
    overlap = floor(win/2);
    nfft = max(1024,2^nextpow2(win));
    [P,f] = pwelch(x,win,overlap,nfft,fs);
    lo = f>=4 & f<30; hi = f>=30 & f<=45;
    ratio = trapz(f(hi),P(hi))/(trapz(f(lo),P(lo))+eps);
end

function idx = find_saved_ic_list(S,nIC,keywords)
    idx = [];
    fn = fieldnames(S);
    for k = 1:numel(fn)
        nm = lower(string(fn{k}));
        hit = false;
        for q = 1:numel(keywords), hit = hit || contains(nm,lower(keywords(q))); end
        if ~hit, continue; end
        x = S.(fn{k});
        if ~(isnumeric(x)||islogical(x)), continue; end
        x = double(x(:)');
        if numel(x)==nIC && all(ismember(unique(x),[0 1]))
            idx = unique([idx find(x~=0)],'stable');
        elseif numel(x)<=nIC && all(isfinite(x)) && all(mod(x,1)==0) && ...
                all(x>=1) && all(x<=nIC)
            idx = unique([idx x],'stable');
        end
    end
end

function metric = find_saved_metric(S,nIC,keywords)
    C = [];
    fn = fieldnames(S);
    for k = 1:numel(fn)
        nm = lower(string(fn{k}));
        hit = true;
        for q = 1:numel(keywords), hit = hit && contains(nm,lower(keywords(q))); end
        if ~hit, continue; end
        x = S.(fn{k});
        if ~isnumeric(x), continue; end
        x = double(x);
        if isvector(x) && numel(x)==nIC
            C = [C abs(x(:))]; %#ok<AGROW>
        elseif size(x,1)==nIC
            C = [C max(abs(x),[],2,'omitnan')]; %#ok<AGROW>
        elseif size(x,2)==nIC
            C = [C max(abs(x),[],1,'omitnan')']; %#ok<AGROW>
        end
    end
    if isempty(C), metric = nan(nIC,1); else, metric = max(C,[],2,'omitnan'); end
end

function [xy,labels] = load_subject_xy(datasetRoot,subj,topoLabels)
    eegDir = fullfile(datasetRoot,char(subj),'eeg');
    f = fullfile(eegDir,char(subj + "_space-CapTrak_electrodes.tsv"));
    topoLabels = string(topoLabels(:)); labels = topoLabels;
    xy = nan(numel(topoLabels),2);

    if exist(f,'file')
        E = readtable(f,'FileType','text','Delimiter','\t','TextType','string');
        vars = lower(string(E.Properties.VariableNames));
        cName = find(vars=="name",1); cX = find(vars=="x",1); cY = find(vars=="y",1);
        if ~isempty(cName)&&~isempty(cX)&&~isempty(cY)
            en = string(E{:,cName}); ex = double(E{:,cX}); ey = double(E{:,cY});
            for k=1:numel(topoLabels)
                ii=find(strcmpi(en,topoLabels(k)),1);
                if ~isempty(ii), xy(k,:)=[ex(ii) ey(ii)]; end
            end
        end
    end

    good=all(isfinite(xy),2);
    if sum(good)>=max(10,round(.7*numel(topoLabels)))
        xy(good,:)=xy(good,:)-mean(xy(good,:),1,'omitnan');
        rad=max(sqrt(sum(xy(good,:).^2,2)));
        if rad>0, xy(good,:)=xy(good,:)/rad; end
        miss=find(~good);
        if ~isempty(miss)
            a=linspace(0,2*pi,numel(miss)+1); a(end)=[];
            xy(miss,:)=.95*[cos(a(:)) sin(a(:))];
        end
    else
        a=linspace(0,2*pi,numel(topoLabels)+1); a(end)=[];
        xy=[cos(a(:)) sin(a(:))];
    end
end

function plot_component_scalp(xy,coordLabels,topoLabels,topo)
    [tf,loc]=ismember(lower(string(topoLabels(:))),lower(string(coordLabels(:))));
    if ~all(tf), loc=(1:numel(topoLabels))'; end
    X=xy(loc,1); Y=xy(loc,2); C=double(topo(:));
    lim=max(abs(C),[],'omitnan'); if ~isfinite(lim)||lim==0, lim=1; end
    scatter(X,Y,85,C,'filled','MarkerEdgeColor',[.25 .25 .25]); hold on;
    th=linspace(0,2*pi,300); plot(1.08*cos(th),1.08*sin(th),'k-','LineWidth',1);
    plot([0 -.08 0 .08 0],[1.08 1.16 1.12 1.16 1.08],'k-','LineWidth',1);
    axis equal; axis([-1.2 1.2 -1.2 1.2]); axis off; caxis([-lim lim]); colorbar;
end

function plot_full_trace(x,fs,maxPts)
    x=double(x(:)); N=numel(x); step=max(1,ceil(N/maxPts)); idx=1:step:N;
    t=(idx-1)/fs; plot(t,x(idx),'k-','LineWidth',.6); grid on; xlim([t(1) t(end)]);
end

function plot_peak_zoom(restX,walkX,fs,zoomSec)
    rz=robust_z(double(restX(:))); wz=robust_z(double(walkX(:)));
    [~,ir]=max(abs(rz)); [~,iw]=max(abs(wz)); h=round((zoomSec/2)*fs);
    rr=max(1,ir-h):min(numel(rz),ir+h); ww=max(1,iw-h):min(numel(wz),iw+h);
    tr=((1:numel(rr))-ceil(numel(rr)/2))/fs; tw=((1:numel(ww))-ceil(numel(ww)/2))/fs;
    plot(tr,rz(rr),'LineWidth',.8); hold on; plot(tw,wz(ww),'LineWidth',.8); yline(0,'k:'); grid on;
    legend({'Rest','Walk'},'Location','best','Box','off');
end

function z = robust_z(x)
    med=median(x,'omitnan'); md=median(abs(x-med),'omitnan'); sc=1.4826*md;
    if ~isfinite(sc)||sc==0, sc=std(x,'omitnan'); end
    if ~isfinite(sc)||sc==0, sc=1; end
    z=(x-med)/sc;
end

function plot_component_psd(restX,walkX,fs)
    r=double(restX(:)); w=double(walkX(:)); win=min([round(2*fs),numel(r),numel(w)]);
    if mod(win,2)==1, win=win-1; end
    ov=floor(win/2); nfft=max(1024,2^nextpow2(win));
    [PR,f]=pwelch(r,win,ov,nfft,fs); [PW,~]=pwelch(w,win,ov,nfft,fs);
    ii=f>=1 & f<=45; plot(f(ii),10*log10(PR(ii)+eps),'LineWidth',.9); hold on;
    plot(f(ii),10*log10(PW(ii)+eps),'LineWidth',.9); xline(4,':'); xline(30,':'); xline(45,':');
    grid on; xlim([1 45]); xlabel('Hz'); ylabel('dB');
    legend({'Rest','Walk'},'Location','best','Box','off');
end

function txt = reason_text(ic,eyeIC,muscleIC)
    e=ismember(ic,eyeIC); m=ismember(ic,muscleIC);
    if e&&m, txt="EYE+MUSCLE"; elseif e, txt="EYE"; elseif m, txt="MUSCLE"; else, txt="AUTO-REJECT"; end
end

function atomic_save_mat(finalFile,S)
    folder=fileparts(finalFile); if ~exist(folder,'dir'), mkdir(folder); end
    tmp=[tempname(folder) '.mat']; save(tmp,'-struct','S','-v7.3');
    [ok,msg]=movefile(tmp,finalFile,'f');
    if ~ok, if exist(tmp,'file'), delete(tmp); end; error('Could not finalize MAT: %s',msg); end
end

function rebuild_manual_audit(R,connectivityRoot,auditFile)
    rows=cell(height(R),1); keep=false(height(R),1);
    for i=1:height(R)
        subj=string(R.Subject(i)); grp=string(R.Group(i));
        f=fullfile(connectivityRoot,char(subj),char(subj + "_FINAL_ICA_reject.mat"));
        if ~exist(f,'file'), continue; end
        F=load(f);
        rows{i}=table(subj,grp,F.NIC,F.NAutoRejected,F.PctAutoRejected,F.NRestored, ...
            F.NFinalRejected,F.PctFinalRejected,string(mat2str(F.autoRejectIC)), ...
            string(mat2str(F.restoredIC)),string(mat2str(F.finalRejectIC)), ...
            string(F.reviewDate),string(F.reviewerNote), ...
            'VariableNames',{'Subject','Group','NIC','NAutoRejected','PctAutoRejected', ...
            'NRestored','NFinalRejected','PctFinalRejected','AutoRejectedIC','RestoredIC', ...
            'FinalRejectedIC','ReviewDate','ReviewerNote'});
        keep(i)=true;
    end
    if any(keep), writetable(vertcat(rows{keep}),auditFile); end
end

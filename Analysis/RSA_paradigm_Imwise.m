%% RSA Paradigm Analysis

%Georgia Milne, 2026
%Using RSA Toolbox (Nili et al, 2014)
%Read subject RSAs from Analysis/Data/derivatives
%Read subject behavioural responses from Analysis/Data/Results

clear;clc; close all hidden;

addpath(genpath("rsatoolbox-develop"))
import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*
import rsa.fsl.*

labelSize = 24;
tickSize = 22;
titleSize = 28;
plotcolors = [136 204 238; 68 170 153; 221 204 119; 204 102 119; 170 68 153; 17 119 51; 136 34 85; 51 34 136];
plotcolors = plotcolors./255;
Colors = [249 150 69; 163 239 163]./255;
Colors2 = [255 181 199; 159 219 243; 255 231 157]./255;
Colors3 = [255, 206, 178; 189, 204, 242; 162, 229, 203; 209, 235, 160]./255;
Colors4 = [51, 34, 136; 17, 119, 51; 204, 102, 119; 136, 34, 85]./255;

Colorsp = [250 230 235; 248 211 221; 245 195 210; 244 180 199; 235 165 192]./255;
Colorsb = [215 238 250; 186 229 246; 172 224 245; 159 219 243; 144 203 238]./255;
Colorsy = [252 243 220; 252 241 195; 252 236 177; 251 231 157; 249 227 131]./255;
Colorsg = [240 240 240; 235 235 235; 230 230 230; 200 200 200; 180 180 180]./255;

%% Load Participants

Within = {'13' '14' '15' '16'};
MixedWithin = {'03' '04' '41' '42'};
MixedAcross = {'05' '06' '09' '10'};
Across = {'17' '18' '19' '20'};
groupNames = {'Within' 'MixedWithin' 'MixedAcross' 'Across'};
subNumsPlus = [Within MixedWithin MixedAcross Across];
maskNames = {"V1_0-9_vol" "PostcentralGyrus"};
subjectNamesPlus = "sub-" + subNumsPlus;

ImN = 20;
n = 60;

%condition indexes
naiveTTidx = logical([ones(1,20) zeros(1,64)]);
GSidx = logical([zeros(1,30) ones(1,20) zeros(1,34)]);
cuedTTidx = logical([zeros(1,54) ones(1,20) zeros(1,10)]);
TTGSidx = logical(naiveTTidx + GSidx + cuedTTidx);
ObIdx = logical([ones(1,10) zeros(1,20) ones(1,10) zeros(1,14) ones(1,10) zeros(1,20)]);
AnIdx = logical([zeros(1,10) ones(1,10) zeros(1,20) ones(1,10) zeros(1,14) ones(1,10) zeros(1,10)]);
matAnIdx = [zeros(1,10), ones(1,10), zeros(1,10), ones(1,10), zeros(1,10), ones(1,10)];
matIdIdx = (((1:60)-1)*(n+1)+1);

mat = [];
for ss = 1:length(subNumsPlus)
    userOptions=projectOptions_fsl_demoGG(subjectNamesPlus(ss));
    si = size(mat,1)+1;
    for mm = 1:length(maskNames)
        mat{si,mm} = load(['Data/derivatives/' subjectNamesPlus{ss} '/mat_' char(maskNames{mm}) '.csv']);
        mat_o{si,mm} = load(['Data/derivatives/' subjectNamesPlus{ss} '/mat_o_' char(maskNames{mm}) '.csv']);
        umat{si,mm} = load(['Data/derivatives/' subjectNamesPlus{ss} '/umat_' char(maskNames{mm}) '.csv']);
        umat_o{si,mm} = load(['Data/derivatives/' subjectNamesPlus{ss} '/umat_o_' char(maskNames{mm}) '.csv']);
    end

    %behavioural
    if contains(subjectNamesPlus{ss}, 'sub-03') || contains(subjectNamesPlus{ss}, 'sub-04') ...
            || contains(subjectNamesPlus{ss}, 'sub-09') || contains(subjectNamesPlus{ss}, 'sub-10')
        userOptions.run_names = userOptions.run_names(1:7);
    elseif contains(subjectNamesPlus{ss}, 'sub-05') || contains(subjectNamesPlus{ss}, 'sub-06')
        userOptions.run_names = userOptions.run_names(1:8);
    else
        userOptions.run_names = userOptions.run_names(1:6);
    end

    if contains(subjectNamesPlus{ss}, 'sub-03')
        files = dir('Data/Results/002');
    else
        files = dir(['Data/Results/0' subNumsPlus{ss}]);
    end
    clearvars runs
    for f = 1:length(files)
        for r = 1:length(userOptions.run_names)
            rr = userOptions.run_names{r}(2);
            if strfind(files(f).name,"events_run-" + rr  + "_scored.tsv")
                runs{r} =  readtable(fullfile(files(f).folder,files(f).name),"FileType","text");
            end
        end
    end

    eventsSheet = runs{1};
    for r = 1:length(runs)-1
        eventsSheet = [eventsSheet; runs{r+1}];
    end

    %find experimental order
    order = sortrows(eventsSheet,2).ImNum;
    GSorder{ss} = unique(order(order<90),'stable')';
    newGSorder{ss} = GSorder{ss};

    if contains(subjectNamesPlus{ss}, 'sub-03') || contains(subjectNamesPlus{ss}, 'sub-04')...
            || contains(subjectNamesPlus{ss}, 'sub-05') || contains(subjectNamesPlus{ss}, 'sub-06')
        newGSorder{ss}(ismember(GSorder{ss},[3 10 15 19 23])) = 0;
        newGSorder{ss}(ismember(newGSorder{ss},[4:9 11:14])) = newGSorder{ss}(ismember(newGSorder{ss},[4:9 11:14])) -1;
        newGSorder{ss}(ismember(newGSorder{ss},16:18)) = newGSorder{ss}(ismember(newGSorder{ss},16:18)) -2;
        newGSorder{ss}(ismember(newGSorder{ss},21)) = newGSorder{ss}(ismember(newGSorder{ss},21)) -4;
        newGSorder{ss}(ismember(newGSorder{ss},22)) = newGSorder{ss}(ismember(newGSorder{ss},22)) -3;
        newGSorder{ss}(ismember(newGSorder{ss},24)) = newGSorder{ss}(ismember(newGSorder{ss},24)) -6;
    end

    if contains(subjectNamesPlus{ss}, 'sub-03') || contains(subjectNamesPlus{ss}, 'sub-04')...
            || contains(subjectNamesPlus{ss}, 'sub-05') || contains(subjectNamesPlus{ss}, 'sub-06')
        naiveTTidx = logical([1:24 zeros(1,60)]);
        GSidx = logical([zeros(1,30) 1:24 zeros(1,30)]);
        cuedTTidx = logical([zeros(1,54) 1:24 zeros(1,6)]);
        TTGSidx = logical(naiveTTidx + GSidx + cuedTTidx);
        ObIdx = [];
        AnIdx = [];

        %reduce matrices to 20 ims (leave out last run)
        for mm = 1:length(maskNames)
            mat_o{ss,mm} = mat_o{ss,mm}([1:20 25:44 49:68],[1:20 25:44 49:68]);
            umat_o{ss,mm} = umat_o{ss,mm}([1:20 25:44 49:68],[1:20 25:44 49:68]);

            mat{ss,mm} = mat{ss,mm}([sort(GSorder{ss}(1:20)) sort(GSorder{ss}(1:20))+24 sort(GSorder{ss}(1:20))+48],[sort(GSorder{ss}(1:20)) sort(GSorder{ss}(1:20))+24 sort(GSorder{ss}(1:20))+48]);
            umat{ss,mm} = umat{ss,mm}([sort(GSorder{ss}(1:20)) sort(GSorder{ss}(1:20))+24 sort(GSorder{ss}(1:20))+48],[sort(GSorder{ss}(1:20)) sort(GSorder{ss}(1:20))+24 sort(GSorder{ss}(1:20))+48]);

            GSorder{ss} = GSorder{ss}(1:20);
        end
    else
        %condition indexes
        naiveTTidx = logical([ismember(1:20, newGSorder{ss}) zeros(1,64)]);
        GSidx = logical([zeros(1,30) ismember(1:20, newGSorder{ss}) zeros(1,34)]);
        cuedTTidx = logical([zeros(1,54) ismember(1:20, newGSorder{ss}) zeros(1,10)]);
        TTGSidx = logical(naiveTTidx + GSidx + cuedTTidx);
        ObIdx = logical([ones(1,10) zeros(1,20) ones(1,10) zeros(1,14) ones(1,10) zeros(1,20)]);
        AnIdx = logical([zeros(1,10) ones(1,10) zeros(1,20) ones(1,10) zeros(1,14) ones(1,10) zeros(1,10)]);
    end

    meanScores{si} = array2table(sort(GSorder{ss}(1:20))',"VariableNames", {'ImNum'});
    perTrialScores{si} = array2table([sort(GSorder{ss}(1:20))' nan(ImN,12)],"VariableNames", {'ImNum','NaiveCat1','NaiveCat2','NaiveCat3','NaiveCat4','NaiveCat5','NaiveCat6'...
        'CuedCat1','CuedCat2','CuedCat3','CuedCat4','CuedCat5','CuedCat6'});
    for ii = 1:ImN
        meanScores{si}.naiveVerb(ii) = mean(eventsSheet.VerbCorrect(eventsSheet.ImNum == ii & contains(eventsSheet.ImCon, 'a')), 'omitnan');
        NaiveCat = eventsSheet.CatCorrect(eventsSheet.ImNum == ii & contains(eventsSheet.ImCon, 'a'));
        meanScores{si}.naiveCat(ii) = mean(NaiveCat, 'omitnan');

        meanScores{si}.greyVerb(ii) = mean(eventsSheet.VerbCorrect(eventsSheet.ImNum == ii & contains(eventsSheet.ImCon, 'b')), 'omitnan');
        meanScores{si}.greyCat(ii) = mean(eventsSheet.CatCorrect(eventsSheet.ImNum == ii & contains(eventsSheet.ImCon, 'b')), 'omitnan');

        meanScores{si}.cuedVerb(ii) = mean(eventsSheet.VerbCorrect(eventsSheet.ImNum == ii & contains(eventsSheet.ImCon, 'c')), 'omitnan');
        CuedCat = eventsSheet.CatCorrect(eventsSheet.ImNum == ii & contains(eventsSheet.ImCon, 'c'));
        meanScores{si}.cuedCat(ii) = mean(CuedCat, 'omitnan');

        for r = 1:length(CuedCat)
            perTrialScores{si}(ii,r+1) = {NaiveCat(r)};
            perTrialScores{si}(ii,r+7) = {CuedCat(r)};
        end
    end
end

%% Verbal Behaviour (Figure S1)
    clearvars meandata groupEbars groupdata
    figure('Position', [476 360 700 560]);
    hold on
    perImData = [];
for g = 1:length(groupNames)
    group = groupNames{g};
    idx = contains(subNumsPlus,eval(group));
    clearvars meandata 
    for s = 1:sum(idx)
        ss = max(find(idx,s));
        meandata(s,:) = [mean(meanScores{ss}.naiveVerb,'omitnan') mean(meanScores{ss}.cuedVerb,'omitnan')];
        scatter(1-0.025*length(groupNames) + 0.025*g + 0.025*s, meandata(s,1),150,'filled','MarkerFaceColor', Colors4(g,:),'MarkerFaceAlpha',.4);
        scatter(2-0.025*length(groupNames) + 0.025*g + 0.025*s, meandata(s,2),150,'filled','MarkerFaceColor', Colors4(g,:),'MarkerFaceAlpha',.4);

        %restructure table
        naivetmp = meanScores{ss}(:,2:3); cuedtmp = meanScores{ss}(:,6:7);
        naivetmp = renamevars(naivetmp,["naiveVerb" "naiveCat"],["Verb" "Cat"]);
        naivetmp.Condition = repmat('pre',height(naivetmp),1);
        cuedtmp = renamevars(cuedtmp,["cuedVerb" "cuedCat"],["Verb" "Cat"]);
        cuedtmp.Condition = repmat('pos',height(cuedtmp),1);

       
        tmptable = [meanScores{ss}(:,1) naivetmp; meanScores{ss}(:,1) cuedtmp];
        tmptable.pID = repmat(subjectNamesPlus{ss},height(tmptable),1);
        tmptable.Paradigm = repmat({group},height(tmptable),1);

        perImData = [perImData; tmptable];
    end
    groupdata(g,:) = mean(meandata,1,'omitnan');
    alldata{g} = meandata;
    groupEbars(g,:) = std(meandata)./sqrt(sum(idx));
    
    %stats
    disp(['Diff = ' num2str(mean(meandata(:,2),'omitnan') - mean(meandata(:,1),'omitnan'))])
    [H,p,~,STATS]=ttest2(meandata(:,1),meandata(:,2));
    disp('p = ')
    fprintf('%.10f\n',p)
    disp(STATS)

    %per paradigm model
    disp([group ' model'])
    subset = perImData(contains(perImData.Paradigm,group), :);
    model = fitglme(subset, 'Verb ~ Condition + (1 | pID) + (1 | ImNum)', 'Distribution', 'Binomial', 'Link', 'logit', 'FitMethod', 'Laplace');
    disp(model)
end

for g = 1:length(groupNames)
plot([1 2],groupdata(g,:), 'Color',Colors4(g,:),'LineWidth',6);
end 

set(gca,'xticklabels', {'pre' 'post'}, 'FontSize', labelSize);
    set(gca,'LineWidth',3)

    set(gca,'xlim', [0.5 2.5]);
    set(gca,'xtick', [1 2]);
    set(gca,'ylim', [0 1]);
    set(gca,'ytick', 0:0.5:1, 'yticklabels', {'0' '50' '100'},'FontSize', tickSize)
    ylabel('Image recognition (%)', 'FontSize', labelSize);

    ax = gca;
    ax.LineWidth = 3;
    box off
    axis square

    %stats mixed-across vs all others
    disp(['Mixed-Across Diff = ' num2str(mean(alldata{3}(:,2) - alldata{3}(:,1)))])
    disp(['Others Diff = ' num2str(mean([alldata{1}(:,2) - alldata{1}(:,1);alldata{2}(:,2) - alldata{2}(:,1);alldata{4}(:,2) - alldata{4}(:,1)] ))])
    [~,p,~,STATS]=ttest2(alldata{3}(:,2) - alldata{3}(:,1),[alldata{1}(:,2) - alldata{1}(:,1);alldata{2}(:,2) - alldata{2}(:,1);alldata{4}(:,2) - alldata{4}(:,1)]);
    disp('p = ')
    fprintf('%.10f\n',p)
    disp(STATS)

    %stats across vs within & mixed-within
    disp(['Across Diff = ' num2str(mean(alldata{4}(:,2) - alldata{4}(:,1)))])
    disp(['Others Diff = ' num2str(mean([alldata{1}(:,2) - alldata{1}(:,1);alldata{2}(:,2) - alldata{2}(:,1)] ))])
    [~,p,~,STATS]=ttest2(alldata{4}(:,2) - alldata{4}(:,1),[alldata{1}(:,2) - alldata{1}(:,1);alldata{2}(:,2) - alldata{2}(:,1)]);
    disp('p = ')
    fprintf('%.10f\n',p)
    disp(STATS)

    %all paradigms model
    model_mixed = fitglme(perImData, 'Verb ~ Condition * Paradigm + (1 | pID) + (1 | ImNum)', 'Distribution', 'binomial', 'Link', 'logit')
%% Category Behaviour (Figure S1)
    clearvars meandata groupEbars groupdata alldata
    fig = figure('Position', [476 360 700 560]);
    hold on
for g = 1:length(groupNames)
    group = groupNames{g};
    idx = contains(subNumsPlus,eval(group));
    clearvars meandata 
    for s = 1:sum(idx)
        ss = max(find(idx,s));
        meandata(s,:) = [mean(meanScores{ss}.naiveCat,'omitnan') mean(meanScores{ss}.cuedCat,'omitnan')];
        scatter(1-0.025*length(groupNames) + 0.025*g + 0.025*s, meandata(s,1),150,'filled','MarkerFaceColor', Colors4(g,:),'MarkerFaceAlpha',.4);
        scatter(2-0.025*length(groupNames) + 0.025*g + 0.025*s, meandata(s,2),150,'filled','MarkerFaceColor', Colors4(g,:),'MarkerFaceAlpha',.4);
    end
        groupdata(g,:) = mean(meandata,1,'omitnan');
        alldata{g} = meandata;
    groupEbars(g,:) = std(meandata)./sqrt(sum(idx));
   
    %stats
    disp(['Diff = ' num2str(mean(meandata(:,2),'omitnan') - mean(meandata(:,1),'omitnan'))])
    [~,p]=ttest2(meandata(:,1),meandata(:,2));
    disp('p = ')
    fprintf('%.10f\n',p)

    %per paradigm model
    disp([group ' model'])
    subset = perImData(contains(perImData.Paradigm,group), :);
    model = fitlme(subset, 'Cat ~ Condition + (1 | pID) + (1 | ImNum)');
    disp(model)
end

for g = 1:length(groupNames)
plot([1 2],groupdata(g,:), 'Color',Colors4(g,:),'LineWidth',6);
end 

set(gca,'xticklabels', {'pre' 'post'}, 'FontSize', labelSize);
    set(gca,'LineWidth',3)

    set(gca,'xlim', [0.5 2.5]);
    set(gca,'xtick', [1 2]);
    set(gca,'ylim', [0.5 1]);
    set(gca,'ytick', 0.5:0.25:1, 'yticklabels', {'50' '75' '100'}, 'FontSize', tickSize)
    ylabel('Image categorisation (%)', 'FontSize', labelSize);

    ax = gca;
    ax.LineWidth = 3;
    box off
    axis square

    %stats mixed-across vs all others
    disp(['Mixed-Across Diff = ' num2str(mean(alldata{3}(:,2) - alldata{3}(:,1)))])
    disp(['Others Diff = ' num2str(mean([alldata{1}(:,2) - alldata{1}(:,1);alldata{2}(:,2) - alldata{2}(:,1);alldata{4}(:,2) - alldata{4}(:,1)] ))])
    [~,p,~,STATS]=ttest2(alldata{3}(:,2) - alldata{3}(:,1),[alldata{1}(:,2) - alldata{1}(:,1);alldata{2}(:,2) - alldata{2}(:,1);alldata{4}(:,2) - alldata{4}(:,1)]);
    disp('p = ')
    fprintf('%.10f\n',p)
    disp(STATS)

    %stats across vs within & mixed-within
    disp(['Across Diff = ' num2str(mean(alldata{4}(:,2) - alldata{4}(:,1)))])
    disp(['Others Diff = ' num2str(mean([alldata{1}(:,2) - alldata{1}(:,1);alldata{2}(:,2) - alldata{2}(:,1)] ))])
    [~,p,~,STATS]=ttest2(alldata{4}(:,2) - alldata{4}(:,1),[alldata{1}(:,2) - alldata{1}(:,1);alldata{2}(:,2) - alldata{2}(:,1)]);
    disp('p = ')
    fprintf('%.10f\n',p)
    disp(STATS)

%% Image order RSMs noise-normalised, V1 & S1 (Figs 2A & 3A)

matrix = umat;
for maska = 1:length(maskNames)
    for g = 1:length(groupNames)
        clearvars groupdata
        group = groupNames{g};
        idx = contains(subNumsPlus,eval(group));

        for s = 1:sum(idx)
            ss = find(idx,s);
            groupdata(:,:,s) = 1 - matrix{max(ss),maska};
        end
        groupdata = mean(groupdata,3,'omitnan');

        fig = figure;
        imagesc(groupdata)

        ax1 = gca;
        hold on
        colormap(flip(redblue))
        clim manual
        clim([-1 1]);

        axis(ax1, 'square');

        cb = colorbar(ax1,"eastoutside","FontSize",tickSize);

        set(cb,'YTick',[-1 0 1],'LineWidth',2,'FontSize',tickSize);

        ax2 = ax1;

        set(ax2, 'XAxisLocation', 'top','YAxisLocation','left', 'Color', 'none');
        set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
        set(ax2,'xtick', []);
        set(ax2,'ytick', []);
        xline(ax1,[ 20.5, 40.5 ],'Color', 'k', "LineWidth", 2 );
        yline(ax1,[ 20.5, 40.5 ],'Color', 'k', "LineWidth", 2 );
        ax2.Box = 'on';
        ax2.XColor = [0 0 0 ];
        ax2.LineWidth = 2;
    end
end

%% Image order Similarity barplots, noise-normalised, V1 & S1, imwise
matrix = umat;
for mm = 1:length(maskNames)
    for ss = 1:length(subNumsPlus)

        toplot = 1 - matrix{ss,mm};
        ImN = size(toplot,1)/3;
        demomat = zeros(size(toplot));

        naiveGSOn{ss,mm} = [];
        for i = 1:ImN
            naiveGSOn{ss,mm} = [naiveGSOn{ss,mm} toplot(i,i+ImN)];
            demomat(i,i+ImN) = 1/3;
        end

        cuedGSOn{ss,mm} = [];
        for i = 1:ImN
            cuedGSOn{ss,mm} = [cuedGSOn{ss,mm} toplot(i+ImN,i+(2*ImN))];
            demomat(i+ImN,i+(2*ImN)) = 2/3;
        end

        naivecuedOn{ss,mm} = [];
        for i = 1:ImN
            naivecuedOn{ss,mm} = [naivecuedOn{ss,mm} toplot(i,i+(ImN*2))];
            demomat(i,i+(ImN*2)) = 3/3;
        end

    end

    %plot group averages
    for gg = 1:length(groupNames)
        group = groupNames{gg};
        groupIdx = contains(subNumsPlus,eval(group));
        gIdx =  find(groupIdx);

        imwisemeans(gg,mm,:,:) = [mean(cell2mat(naiveGSOn(groupIdx,mm)),'omitnan'); mean(cell2mat(cuedGSOn(groupIdx,mm)),'omitnan');  mean(cell2mat(naivecuedOn(groupIdx,mm)),'omitnan')];
        imwiseebars(gg,mm,:) = std(squeeze(imwisemeans(gg,mm,:,:)),[],2,'omitnan');

        disp("Per participant stats, " + maskNames{mm})
        for s = 1:length(gIdx)
            ss = gIdx(s);

            disp(' ')
            disp('cue-pre vs cue-post:')
            disp([groupNames{gg} ': ' num2str(ss)])
            [~,p,~,STATS]=ttest(cell2mat(naiveGSOn(ss,mm)),cell2mat(cuedGSOn(ss,mm)));
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(cell2mat(naiveGSOn(ss,mm))));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

            ps(ss,mm,1) = p;
        end

        for s = 1:length(gIdx)
            ss = gIdx(s);
            disp(' ')
            disp('cue-post vs pre-post:')
            disp([groupNames{gg} ': ' num2str(ss)])
            [~,p,~,STATS]=ttest(cell2mat(cuedGSOn(ss,mm)),cell2mat(naivecuedOn(ss,mm)));
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(cell2mat(cuedGSOn(ss,mm))));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

            ps(ss,mm,2) = p;
        end

        for s = 1:length(gIdx)
            ss = gIdx(s);
            disp(' ')
            disp('cue-pre vs pre-post:')
            disp([groupNames{gg} ': ' num2str(ss)])
            [~,p,~,STATS]=ttest(cell2mat(naiveGSOn(ss,mm)),cell2mat(naivecuedOn(ss,mm)));
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2,'significant');
            d = STATS.tstat / sqrt(length(cell2mat(naiveGSOn(ss,mm))));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

            ps(ss,mm,3) = p;
        end

    end
end

fdrPs = mafdr(ps(:),'BHFDR', 'true');
fdrPs = reshape(fdrPs, size(ps));

ps_rounded = round(ps, 2, 'significant');
fdrs_rounded = round(fdrPs, 2, 'significant');

%% Temporally ordered RSMs, noise-normalised, V1 & S1 (Figs 2B & 3B)

matrix = umat_o;
for maska = 1:length(maskNames)
    for g = 1:length(groupNames)
        clearvars groupdata
        group = groupNames{g};
        idx = contains(subNumsPlus,eval(group));

        for s = 1:sum(idx)
            ss = find(idx,s);
            groupdata(:,:,s) = 1 -matrix{max(ss),maska};
        end
        groupdata = mean(groupdata,3,'omitnan');

        figure;
        imagesc(groupdata)

        ax1 = gca;
        hold on
        colormap(flip(redblue))
        clim manual
        clim([-1 1]);

        axis(ax1, 'square');

        cb = colorbar(ax1,"eastoutside");

        set(cb,'YTick',[-1 0 1],'LineWidth',2,'FontSize',tickSize);

        ax2 = ax1;

        set(ax2, 'XAxisLocation', 'top','YAxisLocation','left', 'Color', 'none');
        set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
        set(ax2,'xtick', []);
        set(ax2,'ytick', []);
        xline(ax1,[ 20.5, 40.5 ],'Color', 'k', "LineWidth", 2 );
        yline(ax1,[ 20.5, 40.5 ],'Color', 'k', "LineWidth", 2 );
        ax2.Box = 'on';
        ax2.XColor = [0 0 0 ];
        ax2.LineWidth = 2;

    end
end

%% Temporally ordered Uncorrected & Corrected Similarity barplots, noise-normalised, V1 & S1 imwise

for mm = 1:length(maskNames)
    for ss = 1:length(subNumsPlus)

        toplot = 1 - matrix{ss,mm};
        ImN = size(toplot,1)/3;
        demomat = zeros(size(toplot));

        %remove corresponding images
        for i = 1:ImN
            toplot(i,i+(ImN*2)) = NaN;
            toplot(i,i+(ImN)) = NaN;
            toplot(i+ImN,i+(ImN*2)) = NaN;
            demomat(i,i+(ImN)) = 1/3;
            demomat(i,i+(ImN*2)) = 1/3;
            demomat(i+ImN,i+(ImN*2)) = 1/3;
        end

        %control noncorresponding images
        if ismember(subNumsPlus{ss},Across)
            bb = 10; %run size
        else
            bb = 4; %run size
        end

        naivecuedBlocks{ss,mm} = [];
        naivecuedMeans{ss,mm} = [];
        for r = 1:ImN/bb %per run
            naivecuedBlocks{ss,mm} = [naivecuedBlocks{ss,mm}; toplot(1+(r-1)*bb:r*bb,41+(r-1)*bb:40+r*bb)];
            naivecuedMeans{ss,mm} = [naivecuedMeans{ss,mm}; mean(mean(toplot(1+(r-1)*bb:r*bb,41+(r-1)*bb:40+r*bb),'omitnan'),'omitnan')];
            toplot(1+(r-1)*bb:r*bb,41+(r-1)*bb:40+r*bb) = NaN;
            demomat(1+(r-1)*bb:r*bb,41+(r-1)*bb:40+r*bb) = 2/3*r;
        end

        cuedGSBlocks{ss,mm} = [];
        cuedGSMeans{ss,mm} = [];
        for r = 1:ImN/bb %per run
            cuedGSBlocks{ss,mm} = [cuedGSBlocks{ss,mm}; toplot(1+(r-1)*bb+ImN:r*bb+ImN,41+(r-1)*bb:40+r*bb)];
            cuedGSMeans{ss,mm} = [cuedGSMeans{ss,mm}; mean(mean(toplot(1+(r-1)*bb+ImN:r*bb+ImN,41+(r-1)*bb:40+r*bb),'omitnan'),'omitnan')];
            toplot(1+(r-1)*bb+ImN:r*bb+ImN,41+(r-1)*bb:40+r*bb) = NaN;
            demomat(1+(r-1)*bb+ImN:r*bb+ImN,41+(r-1)*bb:40+r*bb)= 2/3*r;
        end

        naiveGSBlocks{ss,mm} = [];
        naiveGSMeans{ss,mm} = [];
        for r = 1:ImN/bb %per run
            naiveGSBlocks{ss,mm} = [naiveGSBlocks{ss,mm}; toplot(1+(r-1)*bb:r*bb,21+(r-1)*bb:20+r*bb)];
            naiveGSMeans{ss,mm} = [naiveGSMeans{ss,mm}; mean(mean(toplot(1+(r-1)*bb:r*bb,21+(r-1)*bb:20+r*bb),'omitnan'),'omitnan')];
            toplot(1+(r-1)*bb:r*bb,21+(r-1)*bb:20+r*bb) = NaN;
            demomat(1+(r-1)*bb:r*bb,21+(r-1)*bb:20+r*bb)= 2/3*r;
        end

        background{ss,mm,:,:} = toplot;

        %correct
        tmp = cell2mat(naiveGSBlocks(ss,mm));
        confound_cuepre{ss,:,mm} =  tmp(:);

        tmp = cell2mat(cuedGSBlocks(ss,mm));
        confound_cuepost{ss,:,mm} =  tmp(:);

        tmp = cell2mat(naivecuedBlocks(ss,mm));
        confound_prepost{ss,:,mm} =  tmp(:);

        %convert corresponding ims to GS order
        [~, order] = sort(GSorder{ss}(1:20));
        ranks(order) = 1:ImN;
        naiveGSOn_o{ss,mm} = naiveGSOn{ss,mm}(ranks);
        cuedGSOn_o{ss,mm} = cuedGSOn{ss,mm}(ranks);
        naivecuedOn_o{ss,mm} = naivecuedOn{ss,mm}(ranks);

        %correct
        correction_cuepre(ss,:,mm) =  repelem(cell2mat(naiveGSMeans(ss,mm)),bb,1);
        correction_cuepost(ss,:,mm) =  repelem(cell2mat(cuedGSMeans(ss,mm)),bb,1);
        correction_prepost(ss,:,mm) =  repelem(cell2mat(naivecuedMeans(ss,mm)),bb,1);

        corrected_cuepre(ss,:,mm) =  cell2mat(naiveGSOn_o(ss,mm)) - correction_cuepre(ss,:,mm);
        corrected_cuepost(ss,:,mm) = cell2mat(cuedGSOn_o(ss,mm)) - correction_cuepost(ss,:,mm);
        corrected_prepost(ss,:,mm) = cell2mat(naivecuedOn_o(ss,mm)) - correction_prepost(ss,:,mm);
    end

    %plot figures
    for gg = 1:length(groupNames)
        group = groupNames{gg};
        groupIdx = contains(subNumsPlus,eval(group));
        gIdx =  find(groupIdx);


        corrected_imwisemeans(gg,mm,:,:) = [mean(corrected_cuepre(groupIdx,:,mm) ,'omitnan'); mean(corrected_cuepost(groupIdx,:,mm),'omitnan');   mean(corrected_prepost(groupIdx,:,mm),'omitnan')];
        corrected_imwiseebars(gg,mm,:) = std(squeeze(corrected_imwisemeans(gg,mm,:,:)),[],2,'omitnan');

        correction_imwisemeans(gg,mm,:,:) = [mean(correction_cuepre(groupIdx,:,mm) ,'omitnan'); mean(correction_cuepost(groupIdx,:,mm),'omitnan');   mean(correction_prepost(groupIdx,:,mm),'omitnan')];
        correction_imwiseebars(gg,mm,:) = std(squeeze(correction_imwisemeans(gg,mm,:,:)),[],2,'omitnan');

        %uncorrected
        figure('Position', [476 360 400 300]);
        hold on

        meandata = mean(squeeze(imwisemeans(gg,mm,:,:)),2,'omitnan');
        b = bar(meandata,'facecolor', 'flat','CData', Colors2, 'LineWidth',3, 'BarWidth', 0.6);
        e = errorbar([1 2 3], meandata, squeeze(imwiseebars(gg,mm,:)), 'Color', [0 0 0], 'LineWidth', 3, 'LineStyle', 'none','CapSize',10);

        set(gca,'LineWidth',3)
        set(gca,'xlim', [0.4 3.6], 'xticklabels', ' ');
        set(gca,'xticklabels', {'pre v cue' 'post v cue'  'pre v post'}, 'FontSize', labelSize);

        if mean(mean(matrix{1} == umat_o{1})) == 1
            set(gca,'ylim', [0 0.5]);
            set(gca,'ytick', 0:0.1:0.5, 'FontSize', tickSize)
        elseif mean(mean(matrix{1} == mat_o{1})) == 1
            set(gca,'ylim', [0 0.6]);
            set(gca,'ytick', 0:0.1:0.6, 'FontSize', tickSize)
        end

        ax = gca;
        ax.LineWidth = 3;
        box off

        %correction
        figure('Position', [476 360 400 300]);
        hold on

        meandata = mean(squeeze(correction_imwisemeans(gg,mm,:,:)),2,'omitnan');
        b = bar(meandata,'facecolor', 'flat', 'CData', [.9 .9 .9], 'LineWidth',3, 'BarWidth', 0.5,'FaceAlpha',.6);
        errorbar([1 2 3], meandata, squeeze(correction_imwiseebars(gg,mm,:)), 'Color', [0 0 0], 'LineWidth', 3, 'LineStyle', 'none','CapSize',10);

        set(gca,'xticklabels', {'pre v cue' 'post v cue'  'pre v post'}, 'FontSize', labelSize);
        set(gca,'LineWidth',3)
        set(gca,'xlim', [0.4 3.6]);
        h=gca; h.XAxis.TickLength = [0 0];

        if mean(mean(matrix{1} == umat_o{1})) == 1
            set(gca,'ylim', [0 0.5]);
            set(gca,'ytick', 0:0.1:0.5, 'FontSize', tickSize)
        elseif mean(mean(matrix{1} == mat_o{1})) == 1
            set(gca,'ylim', [0 0.5]);
            set(gca,'ytick', 0:0.1:0.5, 'FontSize', tickSize)
        end

        ax = gca;
        ax.LineWidth = 3;
        box off

        %corrected
        fig = figure('Position', [476 360 400 300]);
        hold on

        meandata = mean(squeeze(corrected_imwisemeans(gg,mm,:,:)),2,'omitnan');
        for bb = 1:3
            b = bar(bb, meandata(bb),'facecolor', 'flat', 'CData', [1 1 1 ], 'EdgeColor', Colors2(bb,:), 'LineWidth',6, 'BarWidth', 0.5);
        end
        errorbar([1 2 3], meandata, squeeze(corrected_imwiseebars(gg,mm,:)), 'Color', [0 0 0], 'LineWidth', 3, 'LineStyle', 'none','CapSize',10);

        set(gca,'xticklabels', {'pre v cue' 'post v cue'  'pre v post'}, 'FontSize', labelSize);
        set(gca,'LineWidth',3)
        set(gca,'xlim', [0.4 3.6]);
        h=gca; h.XAxis.TickLength = [0 0];

        if mean(mean(matrix{1} == umat_o{1})) == 1
            set(gca,'ylim', [-0.1 0.2]);
            set(gca,'ytick', -0.1:0.1:0.2, 'FontSize', tickSize)
        elseif mean(mean(matrix{1} == mat_o{1})) == 1
            set(gca,'ylim', [-0.2 0.2]);
            set(gca,'ytick', -0.2:0.1:0.2, 'FontSize', tickSize)
        end

        ax = gca;
        ax.LineWidth = 3;
        box off

        disp("Per participant stats, " + maskNames{mm})
        for s = 1:length(gIdx)
            ss = gIdx(s);

            disp(['%%%%%%%%%%%%%%%%% CONFOUND STATS %%%%%%%%%%%%%%%%%', newline, 'cue-pre vs cue-post:'])
            disp([groupNames{gg} ': ' num2str(s)])
            [~,p,~,STATS]=ttest(confound_cuepre{ss,:,mm},confound_cuepost{ss,:,mm});
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(confound_cuepre{ss,:,mm}));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

            ps1(ss,mm,1) = p;

            disp([newline, 'cue-post vs pre-post:'])
            disp([groupNames{gg} ': ' num2str(s)])
            [~,p,~,STATS]=ttest(confound_cuepost{ss,:,mm},confound_prepost{ss,:,mm});
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(confound_cuepost{ss,:,mm}));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

            ps1(ss,mm,2) = p;

            disp([newline, 'cue-pre vs pre-post:'])
            disp([groupNames{gg} ': ' num2str(s)])
            [~,p,~,STATS]=ttest(confound_cuepre{ss,:,mm},confound_prepost{ss,:,mm});
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(confound_cuepre{ss,:,mm}));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

            ps1(ss,mm,3) = p;
        end

        for s = 1:length(gIdx)
            ss = gIdx(s);

            disp(['%%%%%%%%%%%%%%%%% CORRECTED STATS %%%%%%%%%%%%%%%%%', newline, 'cue-pre vs cue-post:'])
            disp([groupNames{gg} ': ' num2str(s)])
            [~,p,~,STATS]=ttest(corrected_cuepre(ss,:,mm),corrected_cuepost(ss,:,mm));
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(corrected_cuepre(ss,:,mm)));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

            ps2(ss,mm,1) = p;

            disp([newline, 'cue-post vs pre-post:'])
            disp([groupNames{gg} ': ' num2str(s)])
            [~,p,~,STATS]=ttest(corrected_cuepost(ss,:,mm),corrected_prepost(ss,:,mm));
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(corrected_cuepost(ss,:,mm)));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

            ps2(ss,mm,2) = p;

            disp([newline, 'cue-pre vs pre-post:'])
            disp([groupNames{gg} ': ' num2str(s)])
            [~,p,~,STATS]=ttest(corrected_cuepre(ss,:,mm),corrected_prepost(ss,:,mm));
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(corrected_cuepre(ss,:,mm)));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);
            ps2(ss,mm,3) = p;
        end
    end
end

%for supplementary comparison
for mm = 1:2
    for ss = 1:12
        umeanSim(ss,mm) = mean([corrected_cuepre(ss,:,mm) corrected_cuepost(ss,:,mm) corrected_prepost(ss,:,mm)]);
        uAllSim(ss,mm,:) = [corrected_cuepre(ss,:,mm) corrected_cuepost(ss,:,mm) corrected_prepost(ss,:,mm)];
    end
end
%confound fdrs
fdrPs1 = mafdr(ps1(:),'BHFDR', 'true');
fdrPs1 = reshape(fdrPs1, size(ps1));

ps1_rounded = round(ps1, 2, 'significant');
fdrs1_rounded = round(fdrPs1, 2, 'significant');

%corrected fdrs
fdrPs2 = mafdr(ps2(:),'BHFDR', 'true');
fdrPs2 = reshape(fdrPs2, size(ps2));

ps2_rounded = round(ps2, 2, 'significant');
fdrs2_rounded = round(fdrPs2, 2, 'significant');
%% Temporally ordered Corrected Similarity barplots, noise-normalised, V1 & S1 imwise (Figs 2 & 3)

matrix = umat_o;
for mm = 1:length(maskNames)

    %plot figures
    for gg = 1:length(groupNames)
        group = groupNames{gg};
        groupIdx = contains(subNumsPlus,eval(group));
        gIdx =  find(groupIdx);

        %uncorrected
        figure('Position', [476 360 450 350]);
        hold on

        N = length(gIdx);
        a = 17.5;
        x = [linspace(1-(N/a), 1+(N/a), N); linspace(2-(N/a), 2+(N/a), N); linspace(3-(N/a), 3+(N/a), N)];

        tmpmeans = mean([mean(cell2mat(naiveGSOn(gIdx,mm)),2) mean(cell2mat(cuedGSOn(gIdx,mm)),2) mean(cell2mat(naivecuedOn(gIdx,mm)),2)]);

        plot(x(1,[1 end])+[-.1; .1]', [tmpmeans(1) tmpmeans(1)] ,'LineWidth',9,'Color',Colorsp(end,:), 'LineStyle','-')
        plot(x(2,[1 end])+[-.1; .1]', [tmpmeans(2) tmpmeans(2)] ,'LineWidth',9,'Color',Colorsb(end,:), 'LineStyle','-')
        plot(x(3,[1 end])+[-.1; .1]', [tmpmeans(3) tmpmeans(3)] ,'LineWidth',9,'Color',Colorsy(end,:), 'LineStyle','-')

        for s = 1:length(gIdx)
            hold on

            c = s;
            ss = gIdx(s);

            tmpmeans = mean([cell2mat(naiveGSOn(ss,mm))' cell2mat(cuedGSOn(ss,mm))' cell2mat(naivecuedOn(ss,mm))']);
            tmpCIs = bootci(2000, @mean, [cell2mat(naiveGSOn(ss,mm))' cell2mat(cuedGSOn(ss,mm))' cell2mat(naivecuedOn(ss,mm))']);
            NegDi = tmpmeans- tmpCIs(1,:);
            PosDi=tmpCIs(2,:)-tmpmeans;
            errorbar(x(:,s), tmpmeans, NegDi, PosDi, 'Color', [0 0 0], 'LineWidth', 4, 'LineStyle', 'none','CapSize',9)
        end

        set(gca,'ylim', [0 0.5]);
        set(gca,'ytick', 0:.1:.5 , 'FontSize', tickSize)
        set(gca,'xlim', [0.6 3.4]);
        set(gca,'xtick', [], 'FontSize', tickSize)
        set(gca,'LineWidth',3)

        %correction
        figure('Position', [476 360 450 350]);
        hold on


        tmpmeans = [mean(cell2mat(confound_cuepre(gIdx,:,mm)),'omitnan') mean(cell2mat(confound_cuepost(gIdx,:,mm)),'omitnan') mean(cell2mat(confound_prepost(gIdx,:,mm)),'omitnan')];

        plot(x(1,[1 end])+[-.1; .1]', [tmpmeans(1) tmpmeans(1)] ,'LineWidth',9,'Color',Colorsg(end,:), 'LineStyle','-')
        plot(x(2,[1 end])+[-.1; .1]', [tmpmeans(2) tmpmeans(2)] ,'LineWidth',9,'Color',Colorsg(end,:), 'LineStyle','-')
        plot(x(3,[1 end])+[-.1; .1]', [tmpmeans(3) tmpmeans(3)] ,'LineWidth',9,'Color',Colorsg(end,:), 'LineStyle','-')

        for s = 1:length(gIdx)
            hold on
            c = s;
            ss = gIdx(s);

            tmpmeans = mean([confound_cuepre{ss,:,mm} confound_cuepost{ss,:,mm} confound_prepost{ss,:,mm}],'omitnan');
            tmpCIs = bootci(2000, @mean, [confound_cuepre{ss,:,mm}(~isnan(confound_cuepre{ss,:,mm})) confound_cuepost{ss,:,mm}(~isnan(confound_cuepost{ss,:,mm})) confound_prepost{ss,:,mm}(~isnan(confound_prepost{ss,:,mm}))]);
            NegDi = tmpmeans- tmpCIs(1,:);
            PosDi=tmpCIs(2,:)-tmpmeans;
            errorbar(x(:,s), tmpmeans, NegDi, PosDi, 'Color', [0 0 0], 'LineWidth', 4, 'LineStyle', 'none','CapSize',9)
        end

        set(gca,'ylim', [0 0.5]);
        set(gca,'ytick', -0.1:.1:.5 , 'FontSize', tickSize)
        set(gca,'xlim', [0.6 3.4]);
        set(gca,'xtick', [], 'FontSize', tickSize)
        set(gca,'LineWidth',3)

        %corrected
        fig = figure('Position', [476 360 450 350]);
        hold on
        plot([0 4], [0 0],'LineWidth',3,'Color',[.5 .5 .5], 'LineStyle',':')

        tmpmeans = mean([mean(corrected_cuepre(gIdx,:,mm),2) mean(corrected_cuepost(gIdx,:,mm),2) mean(corrected_prepost(gIdx,:,mm),2)]);

        plot(x(1,[1 end])+[-.1; .1]', [tmpmeans(1) tmpmeans(1)] ,'LineWidth',12,'Color',Colorsp(end,:), 'LineStyle','-')
        plot(x(1,[1 end])+[-.07; .07]', [tmpmeans(1) tmpmeans(1)] ,'LineWidth',4,'Color',[1 1 1], 'LineStyle','-')

        plot(x(2,[1 end])+[-.1; .1]', [tmpmeans(2) tmpmeans(2)] ,'LineWidth',12,'Color',Colorsb(end,:), 'LineStyle','-')
        plot(x(2,[1 end])+[-.07; .07]', [tmpmeans(2) tmpmeans(2)] ,'LineWidth',4,'Color',[1 1 1], 'LineStyle','-')

        plot(x(3,[1 end])+[-.1; .1]', [tmpmeans(3) tmpmeans(3)] ,'LineWidth',12,'Color',Colorsy(end,:), 'LineStyle','-')
        plot(x(3,[1 end])+[-.07; .07]', [tmpmeans(3) tmpmeans(3)] ,'LineWidth',4,'Color',[1 1 1], 'LineStyle','-')

        for s = 1:length(gIdx)
            hold on
            c = s;
            ss = gIdx(s);

            tmpmeans = mean([corrected_cuepre(ss,:,mm)' corrected_cuepost(ss,:,mm)' corrected_prepost(ss,:,mm)']);
            tmpCIs = bootci(2000, @mean, [corrected_cuepre(ss,:,mm)' corrected_cuepost(ss,:,mm)' corrected_prepost(ss,:,mm)']);
            NegDi = tmpmeans- tmpCIs(1,:);
            PosDi=tmpCIs(2,:)-tmpmeans;
            errorbar(x(:,s), tmpmeans, NegDi, PosDi, 'Color', [0 0 0], 'LineWidth', 4, 'LineStyle', 'none','CapSize',9)
        end

        set(gca,'ylim', [-0.05 0.2]);
        set(gca,'ytick', -0.1:.1:.2 , 'FontSize', tickSize)
        set(gca,'xlim', [0.6 3.4]);
        set(gca,'xtick', [], 'FontSize', tickSize)
        set(gca,'LineWidth',3)
    end
end

%% SNR comparisons in V1, normalised
figure
hold on
mm = 1;
allData = [];
allPositions = [];
gap = 3;  % controls the spacing between session blocks
pos = 1;  % initial x-axis position
for ss = 1:16
    tmp = [cell2mat(naiveGSOn(ss,mm)) cell2mat(cuedGSOn(ss,mm)) cell2mat(naivecuedOn(ss,mm))]';
    tmp1{ss,:} = tmp(~isnan(tmp));

    tmp = [confound_cuepre{ss,:,mm}; confound_cuepost{ss,:,mm}; confound_prepost{ss,:,mm}];
    tmp2{ss,:} = tmp(~isnan(tmp));

    [~,p,~,STATS]=ttest2(tmp1{ss,:},tmp2{ss,:});
    p_rounded = round(p, 2, 'significant');
    t_rounded = round(STATS.tstat, 2, 'significant');
    sd_rounded = round(STATS.sd, 2, 'significant');
    d = STATS.tstat / sqrt(1/length(tmp1{ss,:})+1/length(tmp2{ss,:}));
    d_rounded = round(d, 2, 'significant');
    fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

    % Store data and x-axis positions
    allData = [allData; tmp1{ss,:}; tmp2{ss,:}];
    allPositions = [allPositions; ...
        repmat(pos,     numel(tmp1{ss,:}), 1);  % tmp1
        repmat(pos + 1, numel(tmp2{ss,:}), 1)]; % tmp2

    pos = pos + gap;  % move to next session block
end

%% Anova of within vs across comparisons

within = [cell2mat(naiveGSOn(ismember(subNumsPlus,Within),1)); cell2mat(cuedGSOn(ismember(subNumsPlus,[Within MixedWithin]),1)); cell2mat(naivecuedOn(ismember(subNumsPlus,Within),1))];
across = [cell2mat(naiveGSOn(~ismember(subNumsPlus,Within),1)); cell2mat(cuedGSOn(~ismember(subNumsPlus,[Within MixedWithin]),1)); cell2mat(naivecuedOn(~ismember(subNumsPlus,Within),1))];

data = [across; within];

groupLabels = [repmat({'A'}, size(across, 1), 1); ...
    repmat({'B'}, size(within, 1), 1)];

groupLabels = repmat(groupLabels,20,1);

subIDs = [MixedWithin MixedAcross Across MixedAcross Across MixedWithin MixedAcross Across Within Within MixedWithin Within]';

ImIds = repmat([1:20],1,length(subIDs))';
subIDs = repmat(subIDs,20,1);

conds = [repmat("PreCue",length([MixedWithin MixedAcross Across]),1); repmat("PostCue",length([MixedAcross Across]),1); repmat("PrePost",length([MixedWithin MixedAcross Across]),1);...
    repmat("PreCue",length(Within),1); repmat("PostCue",length([Within MixedWithin]),1); repmat("PrePost",length(Within),1)];

conds =repmat(conds,20,1);

T = table(groupLabels, subIDs, conds, ImIds, data(:));

T.subIDs = categorical(T.subIDs);
T.groupLabels = categorical(T.groupLabels);
T.conds = categorical(T.conds);
T.ImIds = categorical(T.ImIds);
T.z_corr = atanh(T.Var5);

lme = fitlme(T, 'z_corr ~ groupLabels * conds + (1|subIDs) + (1|ImIds)')

%% RSM correlations

%image order
matrix = umat;
for maska = 1:length(maskNames)
    for g = 1:length(groupNames)
        groupname = groupNames{g};
        group = eval(groupname);
        p = length(group);
        ii = 1;
        for i = 1:p-1
            idx1 = group{i};
            idx1 = contains(subNumsPlus,idx1);
            for j = i+1:p
                idx2 = group{j};
                idx2 = contains(subNumsPlus,idx2);

                [rho{g}(ii,maska),pval{g}(ii,maska)] = corr(matrix{idx1,maska}(:), matrix{idx2,maska}(:));
                ii = ii + 1;

            end
        end

    end
end

corr_means = mean([rho{1}; rho{2}; rho{3}; rho{4}])

%run-order
matrix = umat_o;
for maska = 1:length(maskNames)
    for g = 1:length(groupNames)
        groupname = groupNames{g};
        group = eval(groupname);
        p = length(group);
        ii = 1;
        for i = 1:p-1
            idx1 = group{i};
            idx1 = contains(subNumsPlus,idx1);
            for j = i+1:p
                idx2 = group{j};
                idx2 = contains(subNumsPlus,idx2);

                [rho_o{g}(ii,maska),pval_o{g}(ii,maska)] = corr(matrix{idx1,maska}(:), matrix{idx2,maska}(:));
                ii = ii + 1;

            end
        end

    end
end

corr_means_o = mean([rho_o{1}; rho_o{2}; rho_o{3}; rho_o{4}])


disp([newline, 'im-order vs run-order RSM correlations'])
[~,p,~,STATS]=ttest([rho{1}; rho{2}; rho{3}; rho{4}],[rho_o{1}; rho_o{2}; rho_o{3}; rho_o{4}]);
p_rounded = round(p, 2, 'significant');
t_rounded = round(STATS.tstat, 2, 'significant');
sd_rounded = round(STATS.sd, 2, 'significant');
d = STATS.tstat / sqrt(length([rho{1}; rho{2}; rho{3}; rho{4}]));
d_rounded = round(d, 2, 'significant');
fprintf('p = %.3g, t(%d) = %.3g, Cohen’s d = %.3g, SD: %.3g\n', p_rounded(1), STATS.df(1), t_rounded(1), d_rounded(1), sd_rounded(1));
fprintf('p = %.3g, t(%d) = %.3g, Cohen’s d = %.3g, SD: %.3g\n', p_rounded(2), STATS.df(2), t_rounded(2), d_rounded(2), sd_rounded(2))

%% Image ordered RSM unnormalised, V1 & S1

matrix = mat;
for maska = 1:length(maskNames)
    for g = 1:length(groupNames)
        clearvars groupdata
        group = groupNames{g};
        idx = contains(subNumsPlus,eval(group));

        for s = 1:sum(idx)
            ss = find(idx,s);
            groupdata(:,:,s) = 1 - matrix{max(ss),maska};
        end
        groupdata = mean(groupdata,3,'omitnan');

        figure;
        imagesc(groupdata)

        ax1 = gca;
        hold on
        colormap(flip(redblue))
        clim manual
        clim([-1 1]);

        axis(ax1, 'square');
        cb = colorbar(ax1,"eastoutside","FontSize",tickSize);
        set(cb,'YTick',[-1 0 1],'LineWidth',2,'FontSize',tickSize);
        ax2 = ax1;

        set(ax2, 'XAxisLocation', 'top','YAxisLocation','left', 'Color', 'none');
        set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
        set(ax2,'xtick', []);
        set(ax2,'ytick', []);
        xline(ax1,[ 20.5, 40.5 ],'Color', 'k', "LineWidth", 2 );
        yline(ax1,[ 20.5, 40.5 ],'Color', 'k', "LineWidth", 2 );
        ax2.Box = 'on';
        ax2.XColor = [0 0 0 ];
        ax2.LineWidth = 2;

    end
end

%% Image ordered Similarity unnormalised, V1 & S1

matrix = mat;
for mm = 1:length(maskNames)
    for ss = 1:length(subNumsPlus)

        toplot = 1 - matrix{ss,mm};
        ImN = size(toplot,1)/3;
        demomat = zeros(size(toplot));

        naiveGSOn{ss,mm} = [];
        for i = 1:ImN
            naiveGSOn{ss,mm} = [naiveGSOn{ss,mm} toplot(i,i+ImN)];
            demomat(i,i+ImN) = 1/3;
        end

        cuedGSOn{ss,mm} = [];
        for i = 1:ImN
            cuedGSOn{ss,mm} = [cuedGSOn{ss,mm} toplot(i+ImN,i+(2*ImN))];
            demomat(i+ImN,i+(2*ImN)) = 2/3;
        end

        naivecuedOn{ss,mm} = [];
        for i = 1:ImN
            naivecuedOn{ss,mm} = [naivecuedOn{ss,mm} toplot(i,i+(ImN*2))];
            demomat(i,i+(ImN*2)) = 3/3;
        end
    end

    %plot group averages
    for gg = 1:length(groupNames)
        group = groupNames{gg};
        groupIdx = contains(subNumsPlus,eval(group));
        gIdx =  find(groupIdx);

        imwisemeans(gg,mm,:,:) = [mean(cell2mat(naiveGSOn(groupIdx,mm)),'omitnan'); mean(cell2mat(cuedGSOn(groupIdx,mm)),'omitnan');  mean(cell2mat(naivecuedOn(groupIdx,mm)),'omitnan')];
        imwiseebars(gg,mm,:) = std(squeeze(imwisemeans(gg,mm,:,:)),[],2,'omitnan');

        disp("Per participant stats, " + maskNames{mm} )
        for s = 1:length(gIdx)
            ss = gIdx(s);

            disp(' ')
            disp('cue-pre vs cue-post:')
            disp([groupNames{gg} ': ' num2str(ss)])
            [h,p,CI,STATS]=ttest(cell2mat(naiveGSOn(ss,mm)),cell2mat(cuedGSOn(ss,mm)));
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(cell2mat(naiveGSOn(ss,mm))));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);
            ps(ss,mm,1) = p;
        end

        for s = 1:length(gIdx)
            ss = gIdx(s);

            disp(' ')
            disp('cue-post vs pre-post:')
            disp([groupNames{gg} ': ' num2str(ss)])
            [h,p,CI,STATS]=ttest(cell2mat(cuedGSOn(ss,mm)),cell2mat(naivecuedOn(ss,mm)));
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(cell2mat(cuedGSOn(ss,mm))));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);
            ps(ss,mm,2) = p;
        end

        for s = 1:length(gIdx)
            ss = gIdx(s);

            disp(' ')
            disp('cue-pre vs pre-post:')
            disp([groupNames{gg} ': ' num2str(ss)])
            [h,p,CI,STATS]=ttest(cell2mat(naiveGSOn(ss,mm)),cell2mat(naivecuedOn(ss,mm)));
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2,'significant');
            d = STATS.tstat / sqrt(length(cell2mat(naiveGSOn(ss,mm))));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);
            ps(ss,mm,3) = p;
        end
    end
end

fdrPs = mafdr(ps(:),'BHFDR', 'true');
fdrPs = reshape(fdrPs, size(ps));

ps_rounded = round(ps, 2, 'significant');
fdrs_rounded = round(fdrPs, 2, 'significant');
%% Temporally ordered RSM unnormalised, V1 & S1 (Fig S3)

matrix = mat_o;
for maska = 1:length(maskNames)
    for g = 1:length(groupNames)
        clearvars groupdata
        group = groupNames{g};
        idx = contains(subNumsPlus,eval(group));

        for s = 1:sum(idx)
            ss = find(idx,s);
            groupdata(:,:,s) = 1 -matrix{max(ss),maska};
        end
        groupdata = mean(groupdata,3,'omitnan');

        figure;
        imagesc(groupdata)

        ax1 = gca;
        hold on
        colormap(flip(redblue))
        clim manual
        clim([-1 1]);

        axis(ax1, 'square');
        cb = colorbar(ax1,"eastoutside");
        set(cb,'YTick',[-1 0 1],'LineWidth',2,'FontSize',tickSize);
        ax2 = ax1;

        set(ax2, 'XAxisLocation', 'top','YAxisLocation','left', 'Color', 'none');
        set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
        set(ax2,'xtick', []);
        set(ax2,'ytick', []);
        xline(ax1,[ 20.5, 40.5 ],'Color', 'k', "LineWidth", 2 );
        yline(ax1,[ 20.5, 40.5 ],'Color', 'k', "LineWidth", 2 );
        ax2.Box = 'on';
        ax2.XColor = [0 0 0 ];
        ax2.LineWidth = 2;
    end
end

%% Temporally ordered Corrected Similarity barplots unnormalised, V1 & S1

matrix = mat_o;
for mm = 1:length(maskNames)
    for ss = 1:length(subNumsPlus)

        toplot = 1 - matrix{ss,mm};
        ImN = size(toplot,1)/3;
        demomat = zeros(size(toplot));

        %remove corresponding images
        for i = 1:ImN
            toplot(i,i+(ImN*2)) = NaN;
            toplot(i,i+(ImN)) = NaN;
            toplot(i+ImN,i+(ImN*2)) = NaN;
            demomat(i,i+(ImN)) = 1/3;
            demomat(i,i+(ImN*2)) = 1/3;
            demomat(i+ImN,i+(ImN*2)) = 1/3;
        end

        %control noncorresponding images
        if ismember(subNumsPlus{ss},Across)
            bb = 10; %run size
        else
            bb = 4; %run size
        end

        naivecuedBlocks{ss,mm} = [];
        naivecuedMeans{ss,mm} = [];
        for r = 1:ImN/bb %per run
            naivecuedBlocks{ss,mm} = [naivecuedBlocks{ss,mm}; toplot(1+(r-1)*bb:r*bb,41+(r-1)*bb:40+r*bb)];
            naivecuedMeans{ss,mm} = [naivecuedMeans{ss,mm}; mean(mean(toplot(1+(r-1)*bb:r*bb,41+(r-1)*bb:40+r*bb),'omitnan'),'omitnan')];
            demomat(1+(r-1)*bb:r*bb,41+(r-1)*bb:40+r*bb) = 2/3*r;
        end

        cuedGSBlocks{ss,mm} = [];
        cuedGSMeans{ss,mm} = [];
        for r = 1:ImN/bb %per run
            cuedGSBlocks{ss,mm} = [cuedGSBlocks{ss,mm}; toplot(1+(r-1)*bb+ImN:r*bb+ImN,41+(r-1)*bb:40+r*bb)];
            cuedGSMeans{ss,mm} = [cuedGSMeans{ss,mm}; mean(mean(toplot(1+(r-1)*bb+ImN:r*bb+ImN,41+(r-1)*bb:40+r*bb),'omitnan'),'omitnan')];
            demomat(1+(r-1)*bb+ImN:r*bb+ImN,41+(r-1)*bb:40+r*bb)= 2/3*r;
        end

        naiveGSBlocks{ss,mm} = [];
        naiveGSMeans{ss,mm} = [];
        for r = 1:ImN/bb %per run
            naiveGSBlocks{ss,mm} = [naiveGSBlocks{ss,mm}; toplot(1+(r-1)*bb:r*bb,21+(r-1)*bb:20+r*bb)];
            naiveGSMeans{ss,mm} = [naiveGSMeans{ss,mm}; mean(mean(toplot(1+(r-1)*bb:r*bb,21+(r-1)*bb:20+r*bb),'omitnan'),'omitnan')];
            demomat(1+(r-1)*bb:r*bb,21+(r-1)*bb:20+r*bb)= 2/3*r;
        end

        %correct
        tmp = cell2mat(naiveGSBlocks(ss,mm));
        confound_cuepre{ss,:,mm} =  tmp(:);

        tmp = cell2mat(cuedGSBlocks(ss,mm));
        confound_cuepost{ss,:,mm} =  tmp(:);

        tmp = cell2mat(naivecuedBlocks(ss,mm));
        confound_prepost{ss,:,mm} =  tmp(:);

        %convert corresponding ims to GS order
        [~, order] = sort(GSorder{ss}(1:20));
        ranks(order) = 1:ImN;
        naiveGSOn_o{ss,mm} = naiveGSOn{ss,mm}(ranks);
        cuedGSOn_o{ss,mm} = cuedGSOn{ss,mm}(ranks);
        naivecuedOn_o{ss,mm} = naivecuedOn{ss,mm}(ranks);

        %correct
        correction_cuepre(ss,:,mm) =  repelem(cell2mat(naiveGSMeans(ss,mm)),bb,1);
        correction_cuepost(ss,:,mm) =  repelem(cell2mat(cuedGSMeans(ss,mm)),bb,1);
        correction_prepost(ss,:,mm) =  repelem(cell2mat(naivecuedMeans(ss,mm)),bb,1);

        corrected_cuepre(ss,:,mm) =  cell2mat(naiveGSOn_o(ss,mm)) - correction_cuepre(ss,:,mm);
        corrected_cuepost(ss,:,mm) = cell2mat(cuedGSOn_o(ss,mm)) - correction_cuepost(ss,:,mm);
        corrected_prepost(ss,:,mm) = cell2mat(naivecuedOn_o(ss,mm)) - correction_prepost(ss,:,mm);
    end

    %plot figures
    for gg = 1:length(groupNames)
        group = groupNames{gg};
        groupIdx = contains(subNumsPlus,eval(group));
        gIdx =  find(groupIdx);


        corrected_imwisemeans(gg,mm,:,:) = [mean(corrected_cuepre(groupIdx,:,mm) ,'omitnan'); mean(corrected_cuepost(groupIdx,:,mm),'omitnan');   mean(corrected_prepost(groupIdx,:,mm),'omitnan')];
        corrected_imwiseebars(gg,mm,:) = std(squeeze(corrected_imwisemeans(gg,mm,:,:)),[],2,'omitnan');

        correction_imwisemeans(gg,mm,:,:) = [mean(correction_cuepre(groupIdx,:,mm) ,'omitnan'); mean(correction_cuepost(groupIdx,:,mm),'omitnan');   mean(correction_prepost(groupIdx,:,mm),'omitnan')];
        correction_imwiseebars(gg,mm,:) = std(squeeze(correction_imwisemeans(gg,mm,:,:)),[],2,'omitnan');

        disp("Per participant stats, " + maskNames{mm})
        for s = 1:length(gIdx)
            ss = gIdx(s);

            disp(['%%%%%%%%%%%%%%%%% CONFOUND STATS %%%%%%%%%%%%%%%%%', newline, 'cue-pre vs cue-post:'])
            disp([groupNames{gg} ': ' num2str(s)])
            [~,p,~,STATS]=ttest(confound_cuepre{ss,:,mm},confound_cuepost{ss,:,mm});
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(confound_cuepre{ss,:,mm}));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

            ps1(ss,mm,1) = p;

            disp([newline, 'cue-post vs pre-post:'])
            disp([groupNames{gg} ': ' num2str(s)])
            [~,p,~,STATS]=ttest(confound_cuepost{ss,:,mm},confound_prepost{ss,:,mm});
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(confound_cuepost{ss,:,mm}));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

            ps1(ss,mm,2) = p;

            disp([newline, 'cue-pre vs pre-post:'])
            disp([groupNames{gg} ': ' num2str(s)])
            [~,p,~,STATS]=ttest(confound_cuepre{ss,:,mm},confound_prepost{ss,:,mm});
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(confound_cuepre{ss,:,mm}));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

            ps1(ss,mm,3) = p;
        end

        for s = 1:length(gIdx)
            ss = gIdx(s);

            disp(['%%%%%%%%%%%%%%%%% CORRECTED STATS %%%%%%%%%%%%%%%%%', newline, 'cue-pre vs cue-post:'])
            disp([groupNames{gg} ': ' num2str(s)])
            [~,p,~,STATS]=ttest(corrected_cuepre(ss,:,mm),corrected_cuepost(ss,:,mm));
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(corrected_cuepre(ss,:,mm)));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

            ps2(ss,mm,1) = p;

            disp([newline, 'cue-post vs pre-post:'])
            disp([groupNames{gg} ': ' num2str(s)])
            [~,p,~,STATS]=ttest(corrected_cuepost(ss,:,mm),corrected_prepost(ss,:,mm));
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(corrected_cuepost(ss,:,mm)));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

            ps2(ss,mm,2) = p;

            disp([newline, 'cue-pre vs pre-post:'])
            disp([groupNames{gg} ': ' num2str(s)])
            [~,p,~,STATS]=ttest(corrected_cuepre(ss,:,mm),corrected_prepost(ss,:,mm));
            p_rounded = round(p, 2, 'significant');
            t_rounded = round(STATS.tstat, 2, 'significant');
            sd_rounded = round(STATS.sd, 2, 'significant');
            d = STATS.tstat / sqrt(length(corrected_cuepre(ss,:,mm)));
            d_rounded = round(d, 2, 'significant');
            fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);
            ps2(ss,mm,3) = p;
        end

    end
end

%confound fdrs
fdrPs1 = mafdr(ps1(:),'BHFDR', 'true');
fdrPs1 = reshape(fdrPs1, size(ps1));

ps1_rounded = round(ps1, 2, 'significant');
fdrs1_rounded = round(fdrPs1, 2, 'significant');

%corrected fdrs
fdrPs2 = mafdr(ps2(:),'BHFDR', 'true');
fdrPs2 = reshape(fdrPs2, size(ps2));

ps2_rounded = round(ps2, 2, 'significant');
fdrs2_rounded = round(fdrPs2, 2, 'significant');

%% Temporally ordered Corrected Similarity barplots normalised, V1 & S1 imwise (Fig S3)

for mm = 1:length(maskNames)
    %plot figures
    for gg = 1:length(groupNames)
        group = groupNames{gg};
        groupIdx = contains(subNumsPlus,eval(group));
        gIdx =  find(groupIdx);

        %uncorrected
        figure('Position', [476 360 450 350]);
        hold on
        if mm == 2
            plot([0 4], [0 0],'LineWidth',3,'Color',[.5 .5 .5], 'LineStyle',':')
        end

        N = length(gIdx);
        a = 17.5;
        x = [linspace(1-(N/a), 1+(N/a), N); linspace(2-(N/a), 2+(N/a), N); linspace(3-(N/a), 3+(N/a), N)];

        tmpmeans = mean([mean(cell2mat(naiveGSOn(gIdx,mm)),2) mean(cell2mat(cuedGSOn(gIdx,mm)),2) mean(cell2mat(naivecuedOn(gIdx,mm)),2)]);

        plot(x(1,[1 end])+[-.1; .1]', [tmpmeans(1) tmpmeans(1)] ,'LineWidth',9,'Color',Colorsp(end,:), 'LineStyle','-')
        plot(x(2,[1 end])+[-.1; .1]', [tmpmeans(2) tmpmeans(2)] ,'LineWidth',9,'Color',Colorsb(end,:), 'LineStyle','-')
        plot(x(3,[1 end])+[-.1; .1]', [tmpmeans(3) tmpmeans(3)] ,'LineWidth',9,'Color',Colorsy(end,:), 'LineStyle','-')

        for s = 1:length(gIdx)
            hold on
            c = s;
            ss = gIdx(s);

            tmpmeans = mean([cell2mat(naiveGSOn(ss,mm))' cell2mat(cuedGSOn(ss,mm))' cell2mat(naivecuedOn(ss,mm))']);
            tmpCIs = bootci(2000, @mean, [cell2mat(naiveGSOn(ss,mm))' cell2mat(cuedGSOn(ss,mm))' cell2mat(naivecuedOn(ss,mm))']);
            NegDi = tmpmeans- tmpCIs(1,:);
            PosDi=tmpCIs(2,:)-tmpmeans;
            errorbar(x(:,s), tmpmeans, NegDi, PosDi, 'Color', [0 0 0], 'LineWidth', 4, 'LineStyle', 'none','CapSize',9)
        end

        if mm == 1
            set(gca,'ylim', [0 0.8]);
            set(gca,'ytick', 0:.2:0.8 , 'FontSize', tickSize)
        else
            set(gca,'ylim', [-0.1 0.8]);
            set(gca,'ytick', -0.2:.2:0.8 , 'FontSize', tickSize)
        end
        set(gca,'xlim', [0.6 3.4]);
        set(gca,'xtick', [], 'FontSize', tickSize)
        set(gca,'LineWidth',3)

        %correction
        figure('Position', [476 360 450 350]);
        hold on
        if mm == 2
            plot([0 4], [0 0],'LineWidth',3,'Color',[.5 .5 .5], 'LineStyle',':')
        end
        tmpmeans = [mean(cell2mat(confound_cuepre(gIdx,:,mm)),'omitnan') mean(cell2mat(confound_cuepost(gIdx,:,mm)),'omitnan') mean(cell2mat(confound_prepost(gIdx,:,mm)),'omitnan')];

        plot(x(1,[1 end])+[-.1; .1]', [tmpmeans(1) tmpmeans(1)] ,'LineWidth',9,'Color',Colorsg(end,:), 'LineStyle','-')
        plot(x(2,[1 end])+[-.1; .1]', [tmpmeans(2) tmpmeans(2)] ,'LineWidth',9,'Color',Colorsg(end,:), 'LineStyle','-')
        plot(x(3,[1 end])+[-.1; .1]', [tmpmeans(3) tmpmeans(3)] ,'LineWidth',9,'Color',Colorsg(end,:), 'LineStyle','-')

        for s = 1:length(gIdx)
            hold on
            %c = 4 - length(gIdx) + s;
            c = s;
            ss = gIdx(s);

            tmpmeans = mean([confound_cuepre{ss,:,mm} confound_cuepost{ss,:,mm} confound_prepost{ss,:,mm}],'omitnan');
            tmpCIs = bootci(2000, @mean, [confound_cuepre{ss,:,mm}(~isnan(confound_cuepre{ss,:,mm})) confound_cuepost{ss,:,mm}(~isnan(confound_cuepost{ss,:,mm})) confound_prepost{ss,:,mm}(~isnan(confound_prepost{ss,:,mm}))]);
            NegDi = tmpmeans- tmpCIs(1,:);
            PosDi=tmpCIs(2,:)-tmpmeans;
            errorbar(x(:,s), tmpmeans, NegDi, PosDi, 'Color', [0 0 0], 'LineWidth', 4, 'LineStyle', 'none','CapSize',9)
        end

        if mm == 1
            set(gca,'ylim', [0 0.8]);
            set(gca,'ytick', 0:.2:0.8 , 'FontSize', tickSize)
        else
            set(gca,'ylim', [-0.1 0.8]);
            set(gca,'ytick', -0.2:.2:0.8 , 'FontSize', tickSize)
        end
        set(gca,'xlim', [0.6 3.4]);
        set(gca,'xtick', [], 'FontSize', tickSize)
        set(gca,'LineWidth',3)

        %corrected
        fig = figure('Position', [476 360 450 350]);
        hold on
        plot([0 4], [0 0],'LineWidth',3,'Color',[.5 .5 .5], 'LineStyle',':')


        tmpmeans = mean([mean(corrected_cuepre(gIdx,:,mm),2) mean(corrected_cuepost(gIdx,:,mm),2) mean(corrected_prepost(gIdx,:,mm),2)]);

        plot(x(1,[1 end])+[-.1; .1]', [tmpmeans(1) tmpmeans(1)] ,'LineWidth',12,'Color',Colorsp(end,:), 'LineStyle','-')
        plot(x(1,[1 end])+[-.07; .07]', [tmpmeans(1) tmpmeans(1)] ,'LineWidth',4,'Color',[1 1 1], 'LineStyle','-')

        plot(x(2,[1 end])+[-.1; .1]', [tmpmeans(2) tmpmeans(2)] ,'LineWidth',12,'Color',Colorsb(end,:), 'LineStyle','-')
        plot(x(2,[1 end])+[-.07; .07]', [tmpmeans(2) tmpmeans(2)] ,'LineWidth',4,'Color',[1 1 1], 'LineStyle','-')

        plot(x(3,[1 end])+[-.1; .1]', [tmpmeans(3) tmpmeans(3)] ,'LineWidth',12,'Color',Colorsy(end,:), 'LineStyle','-')
        plot(x(3,[1 end])+[-.07; .07]', [tmpmeans(3) tmpmeans(3)] ,'LineWidth',4,'Color',[1 1 1], 'LineStyle','-')

        for s = 1:length(gIdx)
            hold on

            c = s;
            ss = gIdx(s);

            tmpmeans = mean([corrected_cuepre(ss,:,mm)' corrected_cuepost(ss,:,mm)' corrected_prepost(ss,:,mm)']);
            tmpCIs = bootci(2000, @mean, [corrected_cuepre(ss,:,mm)' corrected_cuepost(ss,:,mm)' corrected_prepost(ss,:,mm)']);
            NegDi = tmpmeans- tmpCIs(1,:);
            PosDi=tmpCIs(2,:)-tmpmeans;
            errorbar(x(:,s), tmpmeans, NegDi, PosDi, 'Color', [0 0 0], 'LineWidth', 4, 'LineStyle', 'none','CapSize',9)
        end

        set(gca,'ylim', [-0.2 0.2]);
        set(gca,'ytick', -0.2:.1:.2 , 'FontSize', tickSize)
        set(gca,'xlim', [0.6 3.4]);
        set(gca,'xtick', [], 'FontSize', tickSize)
        set(gca,'LineWidth',3)

        disp([char(maskNames{mm}), ', ', groupNames{gg}])
        [~,p]=ttest2(squeeze(corrected_imwisemeans(gg,mm,1,:)),squeeze(corrected_imwisemeans(gg,mm,2,:)));
        disp('cue-pre vs cue-post:')
        fprintf('%.30f\n',p)

        [~,p]=ttest2(squeeze(corrected_imwisemeans(gg,mm,2,:)),squeeze(corrected_imwisemeans(gg,mm,3,:)));
        disp('cue-post vs pre-post:')
        fprintf('%.30f\n',p)

        [~,p]=ttest2(squeeze(corrected_imwisemeans(gg,mm,1,:)),squeeze(corrected_imwisemeans(gg,mm,3,:)));
        disp('cue-pre vs pre-post:')
        fprintf('%.30f\n',p)
    end
end

%% SNR comparisons in V1, unnormalised

mm = 1;
allData = [];
allPositions = [];
gap = 3;  % controls the spacing between session blocks
pos = 1;  % initial x-axis position
for ss = 1:16
    tmp = [cell2mat(naiveGSOn(ss,mm)) cell2mat(cuedGSOn(ss,mm)) cell2mat(naivecuedOn(ss,mm))]';
    tmp1{ss,:} = tmp(~isnan(tmp));

    tmp = [confound_cuepre{ss,:,mm}; confound_cuepost{ss,:,mm}; confound_prepost{ss,:,mm}];
    tmp2{ss,:} = tmp(~isnan(tmp));

    [h,p,CI,STATS]=ttest2(tmp1{ss,:},tmp2{ss,:});
    p_rounded = round(p, 2, 'significant');
    t_rounded = round(STATS.tstat, 2, 'significant');
    sd_rounded = round(STATS.sd, 2, 'significant');
    d = STATS.tstat / sqrt(1/length(tmp1{ss,:})+1/length(tmp2{ss,:}));
    d_rounded = round(d, 2, 'significant');
    fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

    % Store data and x-axis positions
    allData = [allData; tmp1{ss,:}; tmp2{ss,:}];
    allPositions = [allPositions; ...
        repmat(pos,     numel(tmp1{ss,:}), 1);  % tmp1
        repmat(pos + 1, numel(tmp2{ss,:}), 1)]; % tmp2

    pos = pos + gap;  % move to next session block
end

%% Anova of within vs across comparisons

within = [cell2mat(naiveGSOn(ismember(subNumsPlus,Within),1)); cell2mat(cuedGSOn(ismember(subNumsPlus,[Within MixedWithin]),1)); cell2mat(naivecuedOn(ismember(subNumsPlus,Within),1))];
across = [cell2mat(naiveGSOn(~ismember(subNumsPlus,Within),1)); cell2mat(cuedGSOn(~ismember(subNumsPlus,[Within MixedWithin]),1)); cell2mat(naivecuedOn(~ismember(subNumsPlus,Within),1))];

data = [across; within];

groupLabels = [repmat({'A'}, size(across, 1), 1); ...
    repmat({'B'}, size(within, 1), 1)];

groupLabels = repmat(groupLabels,20,1);

subIDs = [MixedWithin MixedAcross Across MixedAcross Across MixedWithin MixedAcross Across Within Within MixedWithin Within]';

ImIds = repmat(1:20,1,length(subIDs))';
subIDs = repmat(subIDs,20,1);

conds = [repmat("PreCue",length([MixedWithin MixedAcross Across]),1); repmat("PostCue",length([MixedAcross Across]),1); repmat("PrePost",length([MixedWithin MixedAcross Across]),1);...
    repmat("PreCue",length(Within),1); repmat("PostCue",length([Within MixedWithin]),1); repmat("PrePost",length(Within),1)];

conds =repmat(conds,20,1);

%[stats] = manova(groupLabels,data)

T = table(groupLabels, subIDs, conds, ImIds, data(:));

T.subIDs = categorical(T.subIDs);
T.groupLabels = categorical(T.groupLabels);
T.conds = categorical(T.conds);
T.ImIds = categorical(T.ImIds);
T.z_corr = atanh(T.Var5);

lme = fitlme(T, 'z_corr ~ groupLabels * conds + (1|subIDs) + (1|ImIds)')
%% Mean correlations
for mm = 1
    for ss = 1:16
        meanRSM(ss,mm) = mean(1- mat_o{ss,mm}(:));
        stdRSM(ss,mm) = std(1- mat_o{ss,mm}(:));

        umeanRSM(ss,mm) = mean(1- umat_o{ss,mm}(:));
        ustdRSM(ss,mm) = std(1- umat_o{ss,mm}(:));

    end

    figure
    boxplot([meanRSM(:,mm); umeanRSM(:,mm)], [ones(16,1); ones(16,1)*2])
    title('Mean');
    xticklabels({'unnormalised', 'normalised'})
    ylabel('pearsons r');

    disp("comparison of means: unnormalised vs normalised, " + maskNames{mm})
    [~,p,~,STATS]=ttest(meanRSM(:,mm),umeanRSM(:,mm));
    p_rounded = round(p, 2, 'significant');
    t_rounded = round(STATS.tstat, 2, 'significant');
    sd_rounded = round(STATS.sd, 2, 'significant');
    d = STATS.tstat / sqrt(length(meanRSM(:,mm)));
    d_rounded = round(d, 2, 'significant');
    fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

    figure
    boxplot([stdRSM(:); ustdRSM(:)], [ones(16,1); ones(16,1)*2])
    title('Std');

    xticklabels({'unnormalised', 'normalised'})
    ylabel('pearsons r');

    disp("comparison of stds: unnormalised vs normalised, " + maskNames{mm})
    [~,p,~,STATS]=ttest(stdRSM(:,mm),ustdRSM(:,mm));
    p_rounded = round(p, 2, 'significant');
    t_rounded = round(STATS.tstat, 2, 'significant');
    sd_rounded = round(STATS.sd, 2, 'significant');
    d = STATS.tstat / sqrt(length(stdRSM(:,mm)));
    d_rounded = round(d, 2, 'significant');
    fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

end

%% Mean within image similarities
for mm = 1:2
    for ss = 1:12
        meanSim(ss,mm) = mean([corrected_cuepre(ss,:,mm) corrected_cuepost(ss,:,mm) corrected_prepost(ss,:,mm)]);
        AllSim(ss,mm,:) = [corrected_cuepre(ss,:,mm) corrected_cuepost(ss,:,mm) corrected_prepost(ss,:,mm)];

        disp("comparison of means: unnormalised vs 0, " + maskNames{mm})
        [~,p,~,STATS]=ttest(AllSim(ss,mm,:));
        p_rounded = round(p, 2, 'significant');
        t_rounded = round(STATS.tstat, 2, 'significant');
        sd_rounded = round(STATS.sd, 2, 'significant');
        d = STATS.tstat / sqrt(length(AllSim(ss,mm,:)));
        d_rounded = round(d, 2, 'significant');
        fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

        disp("comparison of means: normalised vs 0, " + maskNames{mm})
        [~,p,~,STATS]=ttest(uAllSim(ss,mm,:));
        p_rounded = round(p, 2, 'significant');
        t_rounded = round(STATS.tstat, 2, 'significant');
        sd_rounded = round(STATS.sd, 2, 'significant');
        d = STATS.tstat / sqrt(length(uAllSim(ss,mm,:)));
        d_rounded = round(d, 2, 'significant');
        fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);

    end

    figure
    boxplot([meanSim(:,mm); umeanSim(:,mm)], [ones(12,1); ones(12,1)*2])
    title('Mean corrected similarity');
    xticklabels({'unnormalised', 'normalised'})
    ylabel('pearsons r');

    disp("comparison of means: unnormalised vs normalised, " + maskNames{mm})
    [~,p,~,STATS]=ttest(meanSim(:,mm),umeanSim(:,mm));
    p_rounded = round(p, 2, 'significant');
    t_rounded = round(STATS.tstat, 2, 'significant');
    sd_rounded = round(STATS.sd, 2, 'significant');
    d = STATS.tstat / sqrt(length(meanSim(:,mm)));
    d_rounded = round(d, 2, 'significant');
    fprintf('p = %.3g, t(%d) = %.3g, Cohen’s dz = %.3g, SD: %.3g\n', p_rounded, STATS.df, t_rounded, d_rounded, sd_rounded);
end
%% Image ordered RSM normalised single subject representative, V1 & S1 (Fig S2)

matrix = umat;
for maska = [1 2]
    for g = 1:length(groupNames)
        clearvars groupdata
        group = groupNames{g};
        idx = contains(subNumsPlus,eval(group));

        for s = 1:sum(idx)
            ss = find(idx,s);
            groupdata(:,:,s) = 1 - matrix{max(ss),maska};
        end

        for s = 1:size(groupdata,3)
            fig = figure;
            imagesc(groupdata(:,:,s))

            ax1 = gca;
            hold on
            colormap(flip(redblue))
            clim manual
            clim([-1 1]);
            axis(ax1, 'square');

            cb = colorbar(ax1,"eastoutside");

            set(cb,'YTick',[-1 0 1],'LineWidth',2,'FontSize',tickSize);

            ax2 = ax1;

            set(ax2, 'XAxisLocation', 'top','YAxisLocation','left', 'Color', 'none');
            set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
            set(ax2,'xtick', []);
            set(ax2,'ytick', []);
            xline(ax1,[ 20.5, 40.5 ],'Color', 'k', "LineWidth", 2 );
            yline(ax1,[ 20.5, 40.5 ],'Color', 'k', "LineWidth", 2 );
            ax2.Box = 'on';
            ax2.XColor = [0 0 0 ];
            ax2.LineWidth = 2;
        end
    end
end

%% Temporally ordered RDM normalised single subject, V1 & S1 (Fig S2)

matrix = umat_o;
for maska = [1 2]
    for g = 1:length(groupNames)
        clearvars groupdata
        group = groupNames{g};
        idx = contains(subNumsPlus,eval(group));

        for s = 1:sum(idx)
            ss = find(idx,s);
            groupdata(:,:,s) = 1 - matrix{max(ss),maska};
        end

        for s = 1:size(groupdata,3)
            fig = figure;
            imagesc(groupdata(:,:,s))

            ax1 = gca;
            hold on
            colormap(flip(redblue))
            clim manual
            clim([-1 1]);

            axis(ax1, 'square');

            cb = colorbar(ax1,"eastoutside");

            set(cb,'YTick',[-1 0 1],'LineWidth',2,'FontSize',tickSize);

            ax2 = ax1;

            set(ax2, 'XAxisLocation', 'top','YAxisLocation','left', 'Color', 'none');
            set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
            set(ax2,'xtick', []);
            set(ax2,'ytick', []);
            xline(ax1,[ 20.5, 40.5 ],'Color', 'k', "LineWidth", 2 );
            yline(ax1,[ 20.5, 40.5 ],'Color', 'k', "LineWidth", 2 );
            ax2.Box = 'on';
            ax2.XColor = [0 0 0 ];
            ax2.LineWidth = 2;

        end
    end
end

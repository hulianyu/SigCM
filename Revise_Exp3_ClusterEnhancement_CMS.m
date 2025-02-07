addpath([cd '/']);
addpath([cd '/Datasets']);
%% Load Data sets
filename = char('lenses','lung-cancer','soybean-small','zoo','dna-promoter',...
    'hayes-roth','lymphography','heart-disease','solar-flare','primary-tumor',...
    'dermatology','house-votes','balance-scale','credit-approval','breast-cancer-wisconsin',...
    'mammographic-mass','tic-tac-toe','car');
rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf',...
    'Pt','De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce'};
Alpha = 0.05;
%% Partitions: CMS
load('Partitions_CMS.mat');
RT = 100;
Enhancement_ACC_CMS = zeros(18,RT);
Enhancement_NMI_CMS = zeros(18,RT);
Enhancement_ARI_CMS = zeros(18,RT);
Enhancement_FSC_CMS = zeros(18,RT);
for I=1:18
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1);
    GT = X_data(:,1); %Ground Truth
    %
    parfor run = 1:RT
        disp(I);
        disp(run);
        pi = Partitions_CMS{I,1}(:,run);
        [~, ~, pi] = unique(pi);
        pi_new = pi;
        maxK = max(pi);
        [Nk, ~] = histcounts(pi, 0.5:maxK+0.5);
        wait = [];
        for k =1:maxK
            clusterk = find(pi==k);
            list = SigCM_intra(X,clusterk);
            setK = setdiff(1:maxK,k);
            pv_point = Assign_Control(list,X,k,pi);
            pv_point = [list(:,2) pv_point];
            for th = 1:size(pv_point,1)
                [minpv,b] = min(pv_point(th,:)); 
                if b~=1
                    N_that_K = Nk(setK(b-1));
                    if minpv <= Alpha/N_that_K
                        wait = [wait;list(th,1) setK(b-1)];
                    end
                end
            end
        end
        if ~isempty(wait)
            pi_new(wait(:,1)) = wait(:,2);
            all_metrics = ClusteringMeasure(GT, pi_new);
            Enhancement_ACC_CMS(I,run) = all_metrics(1);
            Enhancement_NMI_CMS(I,run) = all_metrics(2);
            Enhancement_ARI_CMS(I,run) = all_metrics(4);
            Enhancement_FSC_CMS(I,run) = all_metrics(7);
        else
            all_metrics = ClusteringMeasure(GT, pi);
            Enhancement_ACC_CMS(I,run) = all_metrics(1);
            Enhancement_NMI_CMS(I,run) = all_metrics(2);
            Enhancement_ARI_CMS(I,run) = all_metrics(4);
            Enhancement_FSC_CMS(I,run) = all_metrics(7);
        end
    end
end
load('Exp1_CMIvsACC_NMI_ARI_Fscore.mat')
%% ACC 
Enhancement_ACC_CMS_mean = mean(Enhancement_ACC_CMS,2);
Enhancement_ACC_CMS_mean_up = 100*(Enhancement_ACC_CMS_mean-com_ACC_CMS_mean)./com_ACC_CMS_mean;
%
threeBars_ACC = zeros(18,3);
for I=1:18
    PLOT = [com_ACC_CMS(I,:)' Enhancement_ACC_CMS(I,:)'];
    threeBars_ACC(I,:) = [sum(PLOT(:,2)>PLOT(:,1)) sum(PLOT(:,2)==PLOT(:,1)) sum(PLOT(:,2)<PLOT(:,1))];
end
%% NMI
Enhancement_NMI_CMS_mean = mean(Enhancement_NMI_CMS,2);
Enhancement_NMI_CMS_mean_up = 100*(Enhancement_NMI_CMS_mean-com_NMI_CMS_mean)./com_NMI_CMS_mean;
%
threeBars_NMI = zeros(18,3);
for I=1:18
    PLOT = [com_NMI_CMS(I,:)' Enhancement_NMI_CMS(I,:)'];
    threeBars_NMI(I,:) = [sum(PLOT(:,2)>PLOT(:,1)) sum(PLOT(:,2)==PLOT(:,1)) sum(PLOT(:,2)<PLOT(:,1))];
end
%% ARI
Enhancement_ARI_CMS_mean = mean(Enhancement_ARI_CMS,2);
Enhancement_ARI_CMS_mean_up = 100*(Enhancement_ARI_CMS_mean-com_ARI_CMS_mean)./com_ARI_CMS_mean;
%
threeBars_ARI = zeros(18,3);
for I=1:18
    PLOT = [com_ARI_CMS(I,:)' Enhancement_ARI_CMS(I,:)'];
    threeBars_ARI(I,:) = [sum(PLOT(:,2)>PLOT(:,1)) sum(PLOT(:,2)==PLOT(:,1)) sum(PLOT(:,2)<PLOT(:,1))];
end
%% FSC
Enhancement_FSC_CMS_mean = mean(Enhancement_FSC_CMS,2);
Enhancement_FSC_CMS_mean_up = 100*(Enhancement_FSC_CMS_mean-com_FSC_CMS_mean)./com_FSC_CMS_mean;
%
threeBars_FSC = zeros(18,3);
for I=1:18
    PLOT = [com_FSC_CMS(I,:)' Enhancement_FSC_CMS(I,:)'];
    threeBars_FSC(I,:) = [sum(PLOT(:,2)>PLOT(:,1)) sum(PLOT(:,2)==PLOT(:,1)) sum(PLOT(:,2)<PLOT(:,1))];
end
%% Plot
figure
Plot_Promote_Bars(threeBars_ACC, 'ACC (CMS) @Enhancement');
figure
Plot_Promote_Bars(threeBars_NMI, 'NMI (CMS) @Enhancement');
figure
Plot_Promote_Bars(threeBars_ARI, 'ARI (CMS) @Enhancement');
figure
Plot_Promote_Bars(threeBars_FSC, 'FSC (CMS) @Enhancement');

%% save 'Enhancement_' vars 
vars = who('Enhancement_*');
% save('Revise_Exp3_ClusterEnhancement_CMS.mat', vars{:});
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
%% Partitions: COForest
load('Partitions_COForest.mat');
RT = 100;
Enhancement_ACC_COForest = zeros(18,RT);
Enhancement_NMI_COForest = zeros(18,RT);
Enhancement_ARI_COForest = zeros(18,RT);
Enhancement_FSC_COForest = zeros(18,RT);
for I=1:18
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1);
    GT = X_data(:,1); %Ground Truth
    %
    parfor run = 1:RT
        disp(I);
        disp(run);
        pi = Partitions_COForest{I,1}(:,run);
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
            Enhancement_ACC_COForest(I,run) = all_metrics(1);
            Enhancement_NMI_COForest(I,run) = all_metrics(2);
            Enhancement_ARI_COForest(I,run) = all_metrics(4);
            Enhancement_FSC_COForest(I,run) = all_metrics(7);
        else
            all_metrics = ClusteringMeasure(GT, pi);
            Enhancement_ACC_COForest(I,run) = all_metrics(1);
            Enhancement_NMI_COForest(I,run) = all_metrics(2);
            Enhancement_ARI_COForest(I,run) = all_metrics(4);
            Enhancement_FSC_COForest(I,run) = all_metrics(7);
        end
    end
end
load('Revise2_Exp1_CMIvsACC_NMI_ARI_Fscore.mat')
%% ACC 
Enhancement_ACC_COForest_mean = mean(Enhancement_ACC_COForest,2);
Enhancement_ACC_COForest_mean_up = 100*(Enhancement_ACC_COForest_mean-com_ACC_COForest_mean)./com_ACC_COForest_mean;
%
threeBars_ACC = zeros(18,3);
for I=1:18
    PLOT = [com_ACC_COForest(I,:)' Enhancement_ACC_COForest(I,:)'];
    threeBars_ACC(I,:) = [sum(PLOT(:,2)>PLOT(:,1)) sum(PLOT(:,2)==PLOT(:,1)) sum(PLOT(:,2)<PLOT(:,1))];
end
%% NMI
Enhancement_NMI_COForest_mean = mean(Enhancement_NMI_COForest,2);
Enhancement_NMI_COForest_mean_up = 100*(Enhancement_NMI_COForest_mean-com_NMI_COForest_mean)./com_NMI_COForest_mean;
%
threeBars_NMI = zeros(18,3);
for I=1:18
    PLOT = [com_NMI_COForest(I,:)' Enhancement_NMI_COForest(I,:)'];
    threeBars_NMI(I,:) = [sum(PLOT(:,2)>PLOT(:,1)) sum(PLOT(:,2)==PLOT(:,1)) sum(PLOT(:,2)<PLOT(:,1))];
end
%% ARI
Enhancement_ARI_COForest_mean = mean(Enhancement_ARI_COForest,2);
Enhancement_ARI_COForest_mean_up = 100*(Enhancement_ARI_COForest_mean-com_ARI_COForest_mean)./com_ARI_COForest_mean;
%
threeBars_ARI = zeros(18,3);
for I=1:18
    PLOT = [com_ARI_COForest(I,:)' Enhancement_ARI_COForest(I,:)'];
    threeBars_ARI(I,:) = [sum(PLOT(:,2)>PLOT(:,1)) sum(PLOT(:,2)==PLOT(:,1)) sum(PLOT(:,2)<PLOT(:,1))];
end
%% FSC
Enhancement_FSC_COForest_mean = mean(Enhancement_FSC_COForest,2);
Enhancement_FSC_COForest_mean_up = 100*(Enhancement_FSC_COForest_mean-com_FSC_COForest_mean)./com_FSC_COForest_mean;
%
threeBars_FSC = zeros(18,3);
for I=1:18
    PLOT = [com_FSC_COForest(I,:)' Enhancement_FSC_COForest(I,:)'];
    threeBars_FSC(I,:) = [sum(PLOT(:,2)>PLOT(:,1)) sum(PLOT(:,2)==PLOT(:,1)) sum(PLOT(:,2)<PLOT(:,1))];
end
%% Plot
figure
Plot_Promote_Bars(threeBars_ACC, 'ACC (COForest) @Enhancement');
figure
Plot_Promote_Bars(threeBars_NMI, 'NMI (COForest) @Enhancement');
figure
Plot_Promote_Bars(threeBars_ARI, 'ARI (COForest) @Enhancement');
figure
Plot_Promote_Bars(threeBars_FSC, 'FSC (COForest) @Enhancement');

%% save 'Enhancement_' vars 
vars = who('Enhancement_*');
%save('Revise2_Exp3_ClusterEnhancement_COForest.mat', vars{:});
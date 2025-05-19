addpath([cd '/']);
addpath([cd '/Datasets']);
%% Load Data sets
filename = char('lenses','lung-cancer','soybean-small','zoo','dna-promoter',...
    'hayes-roth','lymphography','heart-disease','solar-flare','primary-tumor',...
    'dermatology','house-votes','balance-scale','credit-approval','breast-cancer-wisconsin',...
    'mammographic-mass','tic-tac-toe','car');
rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf',...
    'Pt','De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce'};
%% Partitions: ADC
load('Partitions_ADC.mat');
RT = 100;
Refinement_ACC_ADC = zeros(18,RT);
Refinement_NMI_ADC = zeros(18,RT);
Refinement_ARI_ADC = zeros(18,RT);
Refinement_FSC_ADC = zeros(18,RT);
for I=1:18
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1);
    GT = X_data(:,1); %Ground Truth
    %
    parfor run = 1:RT
        disp(I);
        disp(run);
        pi = Partitions_ADC{I,1}(:,run);
        [~, ~, pi] = unique(pi);
        maxK = max(pi);
        remain = [];
        for k =1:maxK
            clusterk = find(pi==k);
            pvals = SigCM_intra(X,clusterk);
            [~,~,list] = FWER_Control(pvals);
            if ~isempty(list)
                remain = [remain;list(:,1)];
            end
        end
        if length(remain) > N/2
            GTr = GT(remain);
            [~, ~, GTr] = unique(GTr);
            pir = pi(remain);
            [~, ~, pir] = unique(pir);
            try
                all_metrics = ClusteringMeasure(GTr, pir);
                Refinement_ACC_ADC(I,run) = all_metrics(1);
                Refinement_NMI_ADC(I,run) = all_metrics(2);
                Refinement_ARI_ADC(I,run) = all_metrics(4);
                Refinement_FSC_ADC(I,run) = all_metrics(7);
                % Refinement_ACC_ADC(I,run) = ACC(GTr, pir);
            catch exception
                all_metrics = ClusteringMeasure(pir, GTr);
                Refinement_ACC_ADC(I,run) = all_metrics(1);
                Refinement_NMI_ADC(I,run) = all_metrics(2);
                Refinement_ARI_ADC(I,run) = all_metrics(4);
                Refinement_FSC_ADC(I,run) = all_metrics(7);
                % Refinement_ACC_ADC(I,run) = ACC(pir,GTr);
            end
        else
            all_metrics = ClusteringMeasure(GT, pi);
            Refinement_ACC_ADC(I,run) = all_metrics(1);
            Refinement_NMI_ADC(I,run) = all_metrics(2);
            Refinement_ARI_ADC(I,run) = all_metrics(4);
            Refinement_FSC_ADC(I,run) = all_metrics(7);
            % Refinement_ACC_ADC(I,run) = ACC(GT,pi);
        end
    end
end
load('Revise2_Exp1_CMIvsACC_NMI_ARI_Fscore.mat')
%% ACC 
Refinement_ACC_ADC_mean = mean(Refinement_ACC_ADC,2);
Refinement_ACC_ADC_mean_up = 100*(Refinement_ACC_ADC_mean-com_ACC_ADC_mean)./com_ACC_ADC_mean;
%
threeBars_ACC = zeros(18,3);
for I=1:18
    PLOT = [com_ACC_ADC(I,:)' Refinement_ACC_ADC(I,:)'];
    threeBars_ACC(I,:) = [sum(PLOT(:,2)>PLOT(:,1)) sum(PLOT(:,2)==PLOT(:,1)) sum(PLOT(:,2)<PLOT(:,1))];
end
%% NMI
Refinement_NMI_ADC_mean = mean(Refinement_NMI_ADC,2);
Refinement_NMI_ADC_mean_up = 100*(Refinement_NMI_ADC_mean-com_NMI_ADC_mean)./com_NMI_ADC_mean;
%
threeBars_NMI = zeros(18,3);
for I=1:18
    PLOT = [com_NMI_ADC(I,:)' Refinement_NMI_ADC(I,:)'];
    threeBars_NMI(I,:) = [sum(PLOT(:,2)>PLOT(:,1)) sum(PLOT(:,2)==PLOT(:,1)) sum(PLOT(:,2)<PLOT(:,1))];
end
%% ARI
Refinement_ARI_ADC_mean = mean(Refinement_ARI_ADC,2);
Refinement_ARI_ADC_mean_up = 100*(Refinement_ARI_ADC_mean-com_ARI_ADC_mean)./com_ARI_ADC_mean;
%
threeBars_ARI = zeros(18,3);
for I=1:18
    PLOT = [com_ARI_ADC(I,:)' Refinement_ARI_ADC(I,:)'];
    threeBars_ARI(I,:) = [sum(PLOT(:,2)>PLOT(:,1)) sum(PLOT(:,2)==PLOT(:,1)) sum(PLOT(:,2)<PLOT(:,1))];
end
%% FSC
Refinement_FSC_ADC_mean = mean(Refinement_FSC_ADC,2);
Refinement_FSC_ADC_mean_up = 100*(Refinement_FSC_ADC_mean-com_FSC_ADC_mean)./com_FSC_ADC_mean;
%
threeBars_FSC = zeros(18,3);
for I=1:18
    PLOT = [com_FSC_ADC(I,:)' Refinement_FSC_ADC(I,:)'];
    threeBars_FSC(I,:) = [sum(PLOT(:,2)>PLOT(:,1)) sum(PLOT(:,2)==PLOT(:,1)) sum(PLOT(:,2)<PLOT(:,1))];
end
%% Plot
figure
Plot_Promote_Bars(threeBars_ACC, 'ACC (ADC) @Refinement');
figure
Plot_Promote_Bars(threeBars_NMI, 'NMI (ADC) @Refinement');
figure
Plot_Promote_Bars(threeBars_ARI, 'ARI (ADC) @Refinement');
figure
Plot_Promote_Bars(threeBars_FSC, 'FSC (ADC) @Refinement');

%% save 'Enhancement_' vars 
vars = who('Refinement_*');
%save('Revise2_Exp2_Refinement_ADC.mat', vars{:});
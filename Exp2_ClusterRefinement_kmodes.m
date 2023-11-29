addpath([cd '/']);
addpath([cd '/Datasets']);
%% Load Data sets
filename = char('lenses','lung-cancer','soybean-small','zoo','dna-promoter',...
    'hayes-roth','lymphography','heart-disease','solar-flare','primary-tumor',...
    'dermatology','house-votes','balance-scale','credit-approval','breast-cancer-wisconsin',...
    'mammographic-mass','tic-tac-toe','car');
rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf',...
    'Pt','De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce'};
%% Partitions: kmodes
load('Partitions_kmodes.mat');
RT = 100;
Refinement_ACC_kmodes = zeros(18,RT);
for I=1:18
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1);
    GT = X_data(:,1); %Ground Truth
    %
    parfor run = 1:RT
        disp(I);
        disp(run);
        pi = Partitions_kmodes{I,1}(:,run);
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
                Refinement_ACC_kmodes(I,run) = ACC(GTr, pir);
            catch exception
                Refinement_ACC_kmodes(I,run) = ACC(pir,GTr);
            end
        else
            Refinement_ACC_kmodes(I,run) = ACC(GT,pi);
        end
    end
end
load('Exp1_CMIvsACC.mat')
Refinement_ACC_kmodes_mean = mean(Refinement_ACC_kmodes,2);
Refinement_ACC_kmodes_mean_up = 100*(Refinement_ACC_kmodes_mean-com_ACC_kmodes_mean)./com_ACC_kmodes_mean;
%
threeBars = zeros(18,3);
for I=1:18
    PLOT = [com_ACC_kmodes(I,:)' Refinement_ACC_kmodes(I,:)'];
    threeBars(I,:) = [sum(PLOT(:,2)>PLOT(:,1)) sum(PLOT(:,2)==PLOT(:,1)) sum(PLOT(:,2)<PLOT(:,1))];
end
Plot_Promote_Bars(threeBars);
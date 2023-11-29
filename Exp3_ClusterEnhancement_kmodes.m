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
Enhancement_ACC_kmodes = zeros(18,RT);
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
        pi_new = pi;
        maxK = max(pi);
        wait = [];
        for k =1:maxK
            clusterk = find(pi==k);
            list = SigCM_intra(X,clusterk);
            setK = setdiff(1:maxK,k);
            pv_point = Assign_Control(list,X,k,pi);
            pv_point = [list(:,2) pv_point];
            for th = 1:size(pv_point,1)
                if pv_point(th,1)~=min(pv_point(th,:))
                    [~,b] = min(pv_point(th,:));
                    wait = [wait;list(th,1) setK(b-1)];
                end
            end
        end
        if ~isempty(wait)
            pi_new(wait(:,1)) = wait(:,2);
            Enhancement_ACC_kmodes(I,run) = ACC(GT, pi_new);
        else
            Enhancement_ACC_kmodes(I,run) = ACC(GT,pi);
        end
    end
end
load('Exp1_CMIvsACC.mat')
Enhancement_ACC_kmodes_mean = mean(Enhancement_ACC_kmodes,2);
Enhancement_ACC_kmodes_mean_up = 100*(Enhancement_ACC_kmodes_mean-com_ACC_kmodes_mean)./com_ACC_kmodes_mean;
%
threeBars = zeros(18,3);
for I=1:18
    PLOT = [com_ACC_kmodes(I,:)' Enhancement_ACC_kmodes(I,:)'];
    threeBars(I,:) = [sum(PLOT(:,2)>PLOT(:,1)) sum(PLOT(:,2)==PLOT(:,1)) sum(PLOT(:,2)<PLOT(:,1))];
end
Plot_Promote_Bars(threeBars);
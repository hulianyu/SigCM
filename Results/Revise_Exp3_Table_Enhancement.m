load('Revise_Exp3_ClusterEnhancement_kmodes.mat')
load('Revise_Exp3_ClusterEnhancement_CDCDR.mat')
load('Revise_Exp3_ClusterEnhancement_Het2Hom.mat')
load('Revise_Exp3_ClusterEnhancement_CMS.mat')
load('Exp1_CMIvsACC_NMI_ARI_Fscore.mat')
%% k-modes
Table_kmodes_ACC = [Enhancement_ACC_kmodes_mean';
    (Enhancement_ACC_kmodes_mean_up./100)'];
Table_kmodes_NMI = [Enhancement_NMI_kmodes_mean';
    (Enhancement_NMI_kmodes_mean_up./100)'];
Table_kmodes_ARI = [Enhancement_ARI_kmodes_mean';
    (Enhancement_ARI_kmodes_mean_up./100)'];
Table_kmodes_FSC = [Enhancement_FSC_kmodes_mean';
    (Enhancement_FSC_kmodes_mean_up./100)'];
Table_kmodes = [Table_kmodes_ACC;Table_kmodes_NMI;
Table_kmodes_ARI;Table_kmodes_FSC];
Table_kmodes_Overall = [mean(Enhancement_ACC_kmodes_mean);(mean(Enhancement_ACC_kmodes_mean)-mean(com_ACC_kmodes_mean))/mean(com_ACC_kmodes_mean);
    mean(Enhancement_NMI_kmodes_mean);(mean(Enhancement_NMI_kmodes_mean)-mean(com_NMI_kmodes_mean))/mean(com_NMI_kmodes_mean);
        mean(Enhancement_ARI_kmodes_mean);(mean(Enhancement_ARI_kmodes_mean)-mean(com_ARI_kmodes_mean))/mean(com_ARI_kmodes_mean);
        mean(Enhancement_FSC_kmodes_mean);(mean(Enhancement_FSC_kmodes_mean)-mean(com_FSC_kmodes_mean))/mean(com_FSC_kmodes_mean)];
Table_kmodes = [Table_kmodes Table_kmodes_Overall];
%% CDCDR
Table_CDCDR_ACC = [Enhancement_ACC_CDCDR_mean';
    (Enhancement_ACC_CDCDR_mean_up./100)'];
Table_CDCDR_NMI = [Enhancement_NMI_CDCDR_mean';
    (Enhancement_NMI_CDCDR_mean_up./100)'];
Table_CDCDR_ARI = [Enhancement_ARI_CDCDR_mean';
    (Enhancement_ARI_CDCDR_mean_up./100)'];
Table_CDCDR_FSC = [Enhancement_FSC_CDCDR_mean';
    (Enhancement_FSC_CDCDR_mean_up./100)'];
Table_CDCDR = [Table_CDCDR_ACC;Table_CDCDR_NMI;
Table_CDCDR_ARI;Table_CDCDR_FSC];
Table_CDCDR_Overall = [mean(Enhancement_ACC_CDCDR_mean);(mean(Enhancement_ACC_CDCDR_mean)-mean(com_ACC_CDCDR_mean))/mean(com_ACC_CDCDR_mean);
    mean(Enhancement_NMI_CDCDR_mean);(mean(Enhancement_NMI_CDCDR_mean)-mean(com_NMI_CDCDR_mean))/mean(com_NMI_CDCDR_mean);
        mean(Enhancement_ARI_CDCDR_mean);(mean(Enhancement_ARI_CDCDR_mean)-mean(com_ARI_CDCDR_mean))/mean(com_ARI_CDCDR_mean);
        mean(Enhancement_FSC_CDCDR_mean);(mean(Enhancement_FSC_CDCDR_mean)-mean(com_FSC_CDCDR_mean))/mean(com_FSC_CDCDR_mean)];
Table_CDCDR = [Table_CDCDR Table_CDCDR_Overall];
%% Het2Hom
Table_Het2Hom_ACC = [Enhancement_ACC_Het2Hom_mean';
    (Enhancement_ACC_Het2Hom_mean_up./100)'];
Table_Het2Hom_NMI = [Enhancement_NMI_Het2Hom_mean';
    (Enhancement_NMI_Het2Hom_mean_up./100)'];
Table_Het2Hom_ARI = [Enhancement_ARI_Het2Hom_mean';
    (Enhancement_ARI_Het2Hom_mean_up./100)'];
Table_Het2Hom_FSC = [Enhancement_FSC_Het2Hom_mean';
    (Enhancement_FSC_Het2Hom_mean_up./100)'];
Table_Het2Hom = [Table_Het2Hom_ACC;Table_Het2Hom_NMI;
Table_Het2Hom_ARI;Table_Het2Hom_FSC];
Table_Het2Hom_Overall = [mean(Enhancement_ACC_Het2Hom_mean);(mean(Enhancement_ACC_Het2Hom_mean)-mean(com_ACC_Het2Hom_mean))/mean(com_ACC_Het2Hom_mean);
    mean(Enhancement_NMI_Het2Hom_mean);(mean(Enhancement_NMI_Het2Hom_mean)-mean(com_NMI_Het2Hom_mean))/mean(com_NMI_Het2Hom_mean);
        mean(Enhancement_ARI_Het2Hom_mean);(mean(Enhancement_ARI_Het2Hom_mean)-mean(com_ARI_Het2Hom_mean))/mean(com_ARI_Het2Hom_mean);
        mean(Enhancement_FSC_Het2Hom_mean);(mean(Enhancement_FSC_Het2Hom_mean)-mean(com_FSC_Het2Hom_mean))/mean(com_FSC_Het2Hom_mean)];
Table_Het2Hom = [Table_Het2Hom Table_Het2Hom_Overall];
%% CMS
Table_CMS_ACC = [Enhancement_ACC_CMS_mean';
    (Enhancement_ACC_CMS_mean_up./100)'];
Table_CMS_NMI = [Enhancement_NMI_CMS_mean';
    (Enhancement_NMI_CMS_mean_up./100)'];
Table_CMS_ARI = [Enhancement_ARI_CMS_mean';
    (Enhancement_ARI_CMS_mean_up./100)'];
Table_CMS_FSC = [Enhancement_FSC_CMS_mean';
    (Enhancement_FSC_CMS_mean_up./100)'];
Table_CMS = [Table_CMS_ACC;Table_CMS_NMI;
Table_CMS_ARI;Table_CMS_FSC];
Table_CMS_Overall = [mean(Enhancement_ACC_CMS_mean);(mean(Enhancement_ACC_CMS_mean)-mean(com_ACC_CMS_mean))/mean(com_ACC_CMS_mean);
    mean(Enhancement_NMI_CMS_mean);(mean(Enhancement_NMI_CMS_mean)-mean(com_NMI_CMS_mean))/mean(com_NMI_CMS_mean);
        mean(Enhancement_ARI_CMS_mean);(mean(Enhancement_ARI_CMS_mean)-mean(com_ARI_CMS_mean))/mean(com_ARI_CMS_mean);
        mean(Enhancement_FSC_CMS_mean);(mean(Enhancement_FSC_CMS_mean)-mean(com_FSC_CMS_mean))/mean(com_FSC_CMS_mean)];
Table_CMS = [Table_CMS Table_CMS_Overall];
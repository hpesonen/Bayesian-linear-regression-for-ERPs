function bayesfuncreg

% Analysis for one general amplitude voltage curve

close all

load('cong.mat');
load('incong.mat');

[N_channel,N,N_trial,N_subject] = size(Y0);

h = waitbar(0/N_subject);
t = linspace(0,800,N);
nb = 21;
Hfda = ones(N,nb);
omega = pi/800;
for i = 1:((nb-1)/2)
    Hfda(:,(i*2):(i*2+1)) = [cos(i*omega*t)',sin(i*omega*t)'];
end

% common timewindow length
tw1 = 75;  

% N400 timewindow length
tw2 = 126; 
H_t0 = blkdiag(zeros(tw1),eye(tw2));

% Common element for higher level of hierarchy
Xblock1 = Hfda; 

% Additive incongruent response in time window [300,800]ms
Xblock2 = H_t0*Hfda; 

X0 = [Xblock1,zeros(N,nb)];
X1 = [Xblock1,Xblock2];

M = size(X0,2);

inspectwindow = 76:151; % Timewindow 300-600
K = length(inspectwindow);
tmp = X1-X0;
H = tmp(inspectwindow,:);
erp_alpha = zeros(N_subject,M,N_channel);
erp_Sigma = zeros(M,M,N_subject,N_channel);

%% Channel level analysis
for iSubject = 1:N_subject

    N_trial1 = nnz(Y1(1,1,:,iSubject));
    N_trial0 = nnz(Y0(1,1,:,iSubject));
    N_trialpre = nnz(epsilon_pre(1,1,:,iSubject));
    N_trialpost = nnz(epsilon_post(1,1,:,iSubject));

    
    for iChannel = 1:N_channel
        alpha_bar = zeros(M,1);
        Sigma = 10^2*eye(M);
    
        % Background EEG covariance
        R = zeros(201);

        N_trial = N_trialpre;
        for iTrial = 1:N_trial
           [~,Ri] = corrmtx(epsilon_pre(iChannel,:,iTrial,iSubject)',N-1);
           R = R + Ri; 
        end 
        N_trial = N_trialpost; 
        for iTrial = 1:N_trial 
           [~,Ri] = corrmtx(epsilon_post(iChannel,:,iTrial,iSubject)',N-1);
           R = R + Ri;
        end 
        R = R/(N_trialpre+N_trialpost);     

        N_trial = N_trial0;
        for iTrial = 1:N_trial
            S = X0*Sigma*X0'+R;
            Kgain = Sigma*X0'/S;
            alpha_bar = alpha_bar+Kgain*(Y0(iChannel,:,iTrial,iSubject)' - X0*alpha_bar);
            Sigma = (eye(M) - Kgain*X0)*Sigma;
        end

        N_trial = N_trial1;
        for iTrial = 1:N_trial
            S = X1*Sigma*X1'+R;
            Kgain = Sigma*X1'/S;
            alpha_bar = alpha_bar+Kgain*(Y1(iChannel,:,iTrial,iSubject)' - X1*alpha_bar);
            Sigma = (eye(M) - Kgain*X1)*Sigma;
        end
                
        erp_alpha(iSubject,:,iChannel) = alpha_bar;
        erp_Sigma(:,:,iSubject,iChannel) = Sigma;

    end
    
    waitbar(iSubject/N_subject,h)
end
%%
% Subject and group level analyses

alpha_channel = zeros(M,N_subject);
Sigma_channel = repmat(10^2*eye(M),[1,1,N_subject]);
alpha_group = zeros(M,1);
Sigma_group = 10^2*eye(M);

erp_prob_subject = zeros(N_subject,1);
for iSubject = 1:N_subject
    for iChannel = 1:N_channel
        S = Sigma_channel(:,:,iSubject)+erp_Sigma(:,:,iSubject,iChannel);
        Kgain = Sigma_channel(:,:,iSubject)/S;
        alpha_channel(:,iSubject) = alpha_channel(:,iSubject)+Kgain*(erp_alpha(iSubject,:,iChannel)' - alpha_channel(:,iSubject));
        Sigma_channel(:,:,iSubject) = (eye(M) - Kgain)*Sigma_channel(:,:,iSubject);        
    end
    erp_prob_subject(iSubject) = normcdf(0,mean(H*alpha_channel(:,iSubject)),1/K*sqrt(sum(sum(H*Sigma_channel(:,:,iSubject)*H'))));
    
    S = Sigma_group+Sigma_channel(:,:,iSubject);
    Kgain = Sigma_group/S;
    alpha_group = alpha_group+Kgain*(alpha_channel(:,iSubject) - alpha_group);
    Sigma_group = (eye(2*nb) - Kgain)*Sigma_group;
end
erp_prob_group = normcdf(0,mean(H*alpha_group),1/K*sqrt(sum(sum(H*Sigma_group*H'))));

plot(t(inspectwindow),H*alpha_channel,'b',t(inspectwindow),H*alpha_group,'r','LineWidth',3)
set(gca,'TickDir','out','FontSize',16)
box off
title('\color{blue}Subject level ERP effects \color{black}with \color{red} group level ERP effect')
disp('=======================================')
disp('Subject level ERP effect probabilities:')
disp(erp_prob_subject)
disp('=======================================')
disp('Group level ERP effect probability:')
disp(erp_prob_group)

clear;
clc;

addpath('functions')
addpath('data')

SMALL  = [33, 110, 82];
MEDIUM = [33, 110, 82, 8, 1, 101, 21, 73, 74, 50, 78, 96, 60, 25, 107, 130, 55, 68, 83];
LARGE  = [33, 110, 82, 46 1 2 3 4 5 77 76 83 102 100 99 101 74 68 69 70 71 73 72 80 78 79 81 89 90 84 85 86 88 87 6 21 22 23 25 62 61 63 19 106 104 103 105];
%HUGE   = [33, 110, 82, 1:32, 34:81, 83:109, 111:130];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HUGE   = [33, 110, 82, 1:32, 34:81, 83:109, 111:130];
% Drop Avg hrs: mfg, all house starts and permit variables
HUGE   = [33, 110, 82, 1:32, 34:47 49 60:81, 83:109, 111:130];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A,B] = xlsread('FREDMD.xlsx');
DATA   = A(2:end,1:end);
tcode = A(1,1:end); %% Vector of transformations
series = B(1,2:end); %% The mnemonic for the variables in the panel
dates = B(3:end,1);  %% Dates of observations

% Transform variables
Y = zeros(660,130);
for i = 1:130
    if tcode(i) == 5 || tcode(i) == 6
        Y(:,i) = 100*transxf(DATA(:,i),tcode(i));
    else
        Y(:,i) = transxf(DATA(:,i),tcode(i));
    end
end

Y = Y(3:end,HUGE);
Y = standard(Y);
[T, M] = size(Y);

% Take lags, and correct observations
p = 12;
n_psi = 10;
DIM_TOT = 500;

Ylag = mlag2(Y,p);
Ylag = Ylag(p+1:end,:);
T = T-p;

% Final VAR matrices/vectors
% 1/ VAR matrices for traditional matrix form
X = [Ylag];
Y_all = Y;
Y = Y(p+1:end,:);

% pairs_index = combnk(1:T,2);   % All pairs for the T observations
% T_pairs = size(pairs_index,1); % Total number of pairs is T_pairs, but this is a large number
% k = randperm(T_pairs);         % Permute these pairs randomly...
% j = k(1:2000);                 % and choose 2000 out of all T_pairs (the paper takes only 100)
% dist_original = sqrt(sum((Y(pairs_index(j,1),:) - Y(pairs_index(j,2),:)).^2,2));

for dim = 1:DIM_TOT
    disp(['Dim = ',num2str(dim)]);
    % 0) Do AR(1)
    if dim == 1
        for i=1:M
            forc_err_AR(:,i) = Y(:,i) - ((Ylag(:,i)'*Ylag(:,i))\(Ylag(:,i)'*Y(:,i)))*Ylag(:,i);
        end
    end
    
%     % 1) Do RP of Ylag (v.1, a single draw of PHI, no averaging involved)
%     PHI = genRP(1,dim,size(Ylag,2));
%     Y_RP = Ylag*PHI';
%     
%     for i=1:M
%         forc_err_RP(:,i,dim) = Y(:,i) - Y_RP*((Y_RP'*Y_RP)\(Y_RP'*Y(:,i)));
%     end
    
    % 1) Do RP of Ylag (v2, averaging over n_psi different draws of PHI,
    % using BIC to compute posterior model probs
    forc_err_tmp  = [];
    for s=1:n_psi
       
        PHI = genRP(1,dim,size(Ylag,2));
        Y_RP = Ylag*PHI';
        for i=1:M
            forc_err_tmp(:,i,s) = Y(:,i) - Y_RP*((Y_RP'*Y_RP)\(Y_RP'*Y(:,i)));
        end
        BIC(s) = log(det(forc_err_tmp(:,:,s)'*forc_err_tmp(:,:,s)/size(forc_err_tmp,1))) + (log(size(forc_err_tmp,1))/size(forc_err_tmp,1))*(size(forc_err_tmp,2)*dim);
        tmp = BIC - min(BIC);
        Post_weights = exp(-.5*tmp) / sum(exp(-.5*tmp)); % - min(BIC) in both numerator and denominator is for stability
    end
    for i=1:M
        forc_err_RP(:,i,dim) = (squeeze(forc_err_tmp(:,i,:)))*Post_weights';
    end
    
    % 2) Do PCA of Ylag
    [Y_PCA,L,~,~] = pc_factor(Ylag,dim);
    
    for i=1:M
        forc_err_PCA(:,i,dim) = Y(:,i) - Y_PCA*((Y_PCA'*Y_PCA)\(Y_PCA'*Y(:,i)));
    end
    
    % 3) PCA on 1 factor only
    
    if dim == 1
        % Extract factors and create VAR matric
        [F,~,~,~] = pc_factor(Y_all,dim);
        Flag = mlag2(F,p);
        Flag = Flag(p+1:end,:);
        
        for i=1:M
            forc_err_PCA_1(:,i) = Y(:,i) - Flag*((Flag'*Flag)\(Flag'*Y(:,i)));
        end
    end
    
    % 3) PCA on 3 factor only
    
    if dim == 3
        % Extract factors and create VAR matric
        [F,~,~,~] = pc_factor(Y_all,dim);
        Flag = mlag2(F,p);
        Flag = Flag(p+1:end,:);
        
        for i=1:M
            forc_err_PCA_3(:,i) = Y(:,i) - Flag*((Flag'*Flag)\(Flag'*Y(:,i)));
        end
    end
end

% Plot MSFE
MSFE_AR     = repmat(mean(forc_err_AR.^2)',1,DIM_TOT);
MSFE_PCA    = squeeze(mean(forc_err_PCA.^2));
MSFE_RP     = squeeze(mean(forc_err_RP.^2));
MSFE_PCA_1  = repmat(mean(forc_err_PCA_1.^2)',1,DIM_TOT);
MSFE_PCA_3  = repmat(mean(forc_err_PCA_3.^2)',1,DIM_TOT);

% First plot only means for all models
fullscreen = get(0,'ScreenSize');
h = figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
plot(mean(MSFE_PCA),'r','linewidth',4);
hold on
plot(mean(MSFE_RP),'k','linewidth',4);
plot(mean(MSFE_AR),'b','linewidth',4);
plot(mean(MSFE_PCA_1),'g','linewidth',4);
plot(mean(MSFE_PCA_3),'c','linewidth',4);

legend('PCA','RP','AR(1)','PCA (1 factor)','PCA (3 factors)');
title('Average MSFE across all series in the panel (solid lines)');
xlabel('Projection dimension');
ylabel('MSFE');
grid on;
set(gcf,'PaperPositionMode','auto')
print('-depsc','-tiff','-r600',['Figure_MSFE_1 (full sample).eps'])
close all;

% Next plot top models with bands
% First plot only means for all models
fullscreen = get(0,'ScreenSize');
h = figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
plot(mean(MSFE_PCA),'r','linewidth',4);
hold on
plot(mean(MSFE_RP),'k','linewidth',4);
plot(mean(MSFE_AR),'b','linewidth',4);
plot(prctile(MSFE_PCA,2.5),'r--','linewidth',2);
plot(prctile(MSFE_PCA,97.5),'r--','linewidth',2);
plot(prctile(MSFE_RP,2.5),'k--','linewidth',2);
plot(prctile(MSFE_RP,97.5),'k--','linewidth',2);
plot(prctile(MSFE_AR,2.5),'b--','linewidth',2);
plot(prctile(MSFE_AR,97.5),'b--','linewidth',2);

legend('PCA','RP','AR(1)');
title('Average MSFE across all series in the panel (solid lines), and 95% coverage probabilities (dotted lines)');
xlabel('Projection dimension');
ylabel('MSFE');
grid on;
set(gcf,'PaperPositionMode','auto')
print('-depsc','-tiff','-r600',['Figure_MSFE_2 (full sample).eps'])
close all;

% Plot only for first three series
% First plot only means for all models
fullscreen = get(0,'ScreenSize');
h = figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
for p=1:3
    subplot(3,1,p);
    plot(MSFE_PCA(p,:),'r','linewidth',4);
    hold on
    plot(MSFE_RP(p,:),'k','linewidth',4);
    plot(MSFE_AR(p,:),'b','linewidth',4);
    plot(MSFE_PCA_1(p,:),'g','linewidth',4);
    plot(MSFE_PCA_3(p,:),'c','linewidth',4);
    if p==1
        legend('PCA','RP','AR(1)','PCA (1 factor)','PCA (3 factors)');
    end
    title(['MSFE for series # ',num2str(p)]);
    xlabel('Projection dimension');
    ylabel('MSFE');
    grid on;

end

set(gcf,'PaperPositionMode','auto')
print('-depsc','-tiff','-r600',['Figure_MSFE_3 (full sample).eps'])
close all;
keyboard

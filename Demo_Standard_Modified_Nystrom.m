
% Farhad Pourkamali-Anaraki, farhad.pourkamali@colorado.edu

% This matlab code compares the modified Nystrom method (Nystrom via QR decomposition)
% with standard Nystrom on a synthetic data set. The target rank is fixed r=2 
% and varying numbers of landmark points are used. 


close all;clc;
rng(1) % reproducibility

% generate data
number_of_clusters = 2; % number of clusters
num_sample_cluster = 2000; % number of samples per cluster
n = number_of_clusters * num_sample_cluster; % total number of samples

r = (0:number_of_clusters-1)*3+1; 
R = zeros(number_of_clusters,num_sample_cluster); % radius
for i = 1 : number_of_clusters
    R(i,:) = r(i) + 0.4.*randn(1,num_sample_cluster);
end

Theta = pi.*rand(number_of_clusters,num_sample_cluster); % angle

X = []; % data matrix 
for i = 1 : number_of_clusters
    X = [X,[R(i,:).*cos(Theta(i,:));R(i,:).*sin(Theta(i,:))]];
end

% kernel matrix
d = 2; % polynomial kernel with d=2 (no constant)
KernelMatrix = (X'*X).^d;
rank = 2;
m_all = (1:20)*rank;

% approximation error
error = @(K) norm(KernelMatrix - K,'fro')./norm(KernelMatrix,'fro'); 

% direct eigenvalue decomposition
[UK,SK] = eigs(KernelMatrix, rank);
err_SVD = error( UK * SK * UK' );

num_trial = 10;

err_standard = zeros(num_trial,size(m_all,2));
err_new      = zeros(num_trial,size(m_all,2)); % QR

param.type   = 'uni-sample';
kernel.type  = 'Poly'; 
kernel.par   = [2,0];

for mm = 1 : size(m_all,2)
    m = m_all(mm);
    for trial = 1 : num_trial
        Z = FindRep(X , m , param); % same landmark set for both
           
        L_standard = NysLowRank(X , Z , rank , kernel); % standard Nystrom method
        
        L_new      = NysDecom(X , Z , rank , kernel); % Our approach: QR factorization
        
        err_standard(trial,mm) = error(L_standard * L_standard');
        err_new(trial,mm) = error(L_new * L_new');
    end
end


% plot error
jitter = 0.2;
figure(4);clf;
errorbar( m_all-jitter , err_SVD.*ones(1,size(m_all,2)) , zeros(1,size(m_all,2)) ,'-o','LineWidth',2,'Color',[0,0,1]);
hold all
errorbar( m_all , mean(err_standard) , std(err_standard),':^','LineWidth',2,'Color',[1,0,0]);
errorbar( m_all+jitter , mean(err_new) , std(err_new),'-.p','LineWidth',2,'Color',[0.7490 0 0.7490]);
loc = get(gcf,'position');
set( gcf, 'position', [loc(1), loc(2), 480*1.8,360*1.8] );
hl = legend('EVD',...
    'Standard Nystr\"om','Nystr\"om via QR','location','best');
set(hl,'Interpreter','latex')
set(hl,'fontsize',22);
myList = findobj(gca,'-property','markerfacecolor');
for obj = myList'
    clr = get(obj,'color');
    set(obj,'markerfacecolor',clr);
end
myList = findobj(gca,'-property','markerfacecolor');
for obj = myList'
    set(obj,'markersize',10);
end
set(gca,'fontsize',22);
set(gca,'fontname','times')  % Set it to times
xlabel('number of landmark points m','Interpreter','LaTex');
ylabel('approximation error','Interpreter','LaTex');
ax = gca;
xlim([min(m_all)-.8,max(m_all)]+0.5);
ylim([0.39,0.61])
ax.XTick = m_all;


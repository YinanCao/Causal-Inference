clear; clc; close all;
load('Data_15subjs_22Trls_MEGextract.mat')
Data = Data_sum(:,:,1);

Nsample = 10000; % increase # of samples for smooth curves

pcMean=0.3795; % mean pcom of CI (MA) across participants (empirical values)

pc_vec=[1,pcMean,0];

a1=2.5;
a2=4;
v=3.75;

for whichpc=1:length(pc_vec)

    X = [log(50),  16,   log(a1), log(a2), log(v), log(a1), log(v), 0, 0, pc_vec(whichpc)];
    ep = 10^-32;
    Utility_flag = 'MA'; % decision strategy, model averaging

    data_flag = 'mul';

    do_nonCausalPred=0;

if length(X) ~= 10
    warning('Seriously wrong!! use the right function')
    warning('Seriously wrong!! use the right function')
    warning('Seriously wrong!! use the right function')
    warning('Seriously wrong!! use the right function')
end

if size(Data, 1) ~= 76 % Data must be 76 lines
    warning('Error: incorrect data dimension!')
    warning('Error: incorrect data dimension!')
    warning('Error: incorrect data dimension!')
    warning('Error: incorrect data dimension!')
end

Vlevel = Data(:,3)';
Alevel = Data(:,4)';
Arel = Data(:,2)';
Arel(isnan(Arel)) = 3;

r4 = 11/0.55; % highest physical rate
r1 = 5/0.55;  % lowest physical rate

% parameters:
% log(SDprior)   1
% mu_prior       2
% log(SDa1_high) 3
% log(SDa1_low)  4
% log(SDv1)      5

% log(SDa4_high) 6
% log(SDa4_low)  6
% log(SDv4)      7
% powera         8
% powerv         9
% pc(aud)       10
% pc(vis)       11

% sd of noise was sampled in log-space for multi-starts
% so exponentiated here, and squared for variance.
Varp = exp(X(1)).^2; % variance of prior (perceptual bias)
mup  = X(2);    % mean of prior
var1 = exp( X([3,4,5]) ).^2; % Var of lowest rate % Order: aud_high, aud_low, vis
% var4 = exp( X([6,6,7]) ).^2; % Var of highest rate

var4 = exp( X([3,4,5]) ).^2; % Var of highest rate

% var4 = exp( X([6,6,7]) ).^2; % Var of highest rate

var4(2) = var4(2) + var1(2) - var1(1);
p    = X([8,8,9]); % power exponent

% VarA = var1(Arel) + (var4(Arel)-var1(Arel))./((r4).^p(Arel)-(r1).^p(Arel))...
%     .*(Alevel.^p(Arel) - (r1).^p(Arel));

VarA=var1(Arel);

% VarV = var1(3) + (var4(3)-var1(3))./((r4)^p(3)-(r1).^p(3))...
%     .*(Vlevel.^p(3) - (r1).^p(3));

VarV=var1(3)*ones(1,76);

% [Arel', Vlevel', Alevel', VarV', VarA']

% plot(VarA(65:68), 'ro-')
% hold on;
% plot(VarA(69:72), 'ro--')
% hold on;
% plot(VarV(73:76), 'bo-')
%
pc = X([10, 10]);
task_idx = [ones(32,1); ones(32,1)*2];
pc = pc(task_idx);


N = Nsample; % number of simulated samples
seed1 = randn(N, 1); % Simulated samples
seed2 = randn(N, 1); % Simulated samples

idx = 1:64;
M_n = length(idx);
p_com = repmat(pc, N, 1);    % 10000 * 64
var_p = repmat(Varp, N, M_n);
mu_p  = repmat(mup, N, M_n);

var_a = repmat(VarA(idx), N, 1);
S_a = repmat(Data(idx, 4)', N, 1);
x_a = bsxfun(@times,seed1,sqrt(var_a)) + S_a;

var_v = repmat(VarV(idx), N, 1);
S_v = repmat(Data(idx, 3)', N, 1);
x_v = bsxfun(@times,seed2,sqrt(var_v)) + S_v;

% Causal inference:
p_AVc1 = 1./(2.*pi.*sqrt(var_v.*var_a+var_v.*var_p+var_a.*var_p)).*...
    exp(  -1/2.*((x_v-x_a).^2.*(var_p) + (x_v-mu_p).^2.*(var_a)+(x_a-mu_p).^2.*(var_v))./...
    (var_v.*var_a+var_v.*var_p+var_a.*var_p)  );

p_AVc2 = 1./(2.*pi.*sqrt((var_v+var_p).*(var_a+var_p))).*...
    exp(-1/2.*((x_v-mu_p).^2)./(var_v+var_p)+((x_a-mu_p).^2)./(var_a+var_p));

p_c1AV = (p_AVc1.*p_com) ./...
         (p_AVc1.*(p_com) + p_AVc2.*(1-(p_com)));

% Sest_AVc1 = (x_v./var_v+x_a./var_a+mu_p./var_p)./...
%             (1./var_v+1./var_a+1./var_p);
        
        Sest_AVc1 = (x_v./var_v+x_a./var_a)./...
            (1./var_v+1./var_a);
        
        
        
% Sest_Vc2 = (x_v./var_v + mu_p./var_p)./(1./var_v+1./var_p);
% Sest_Ac2 = (x_a./var_a + mu_p./var_p)./(1./var_a+1./var_p);

Sest_Vc2 = x_v;
Sest_Ac2 = x_a;

switch Utility_flag % decision strategy
    case 'MA' % Model averaging
        Sest_A = p_c1AV.*Sest_AVc1 + (1-p_c1AV).*Sest_Ac2; 
        Sest_V = p_c1AV.*Sest_AVc1 + (1-p_c1AV).*Sest_Vc2;
    case 'PM' % Probability matching
%         thiseta = unifrnd(0,1,N,M_n); % stochastic
        thiseta = rand(N, M_n); % stochastic
        Sest_AVc1(p_c1AV <= thiseta) = 0;
        Sest_Ac2(p_c1AV > thiseta) = 0;
        Sest_Vc2(p_c1AV > thiseta) = 0;
        Sest_A = Sest_AVc1 + Sest_Ac2;
        Sest_V = Sest_AVc1 + Sest_Vc2;
    case 'MS' % Model selection
        Sest_AVc1(p_c1AV <= 0.5) = 0;
        Sest_Ac2(p_c1AV > 0.5) = 0;
        Sest_Vc2(p_c1AV > 0.5) = 0;
        Sest_A = Sest_AVc1 + Sest_Ac2;
        Sest_V = Sest_AVc1 + Sest_Vc2;
    case 'Fusion'
        Sest_A = Sest_AVc1;
        Sest_V = Sest_AVc1;
    case 'Seg'
        Sest_A = Sest_Ac2;
        Sest_V = Sest_Vc2;
end


    switch data_flag
    case 'all'
    % unisensory aud
    idx = 65:72; 
    M_n = length(idx);
    var_a = repmat(VarA(idx), N, 1);
    var_p = repmat(Varp, N, M_n);
    mu_p = repmat(mup, N, M_n);
    S_a = repmat(Data(idx, 4)',N,1);
    x_a = bsxfun(@times,seed,sqrt(var_a)) + S_a;
    Sest_A_uni = (x_a./var_a + mu_p./var_p)./(1./var_a+1./var_p);
    


    
    idx = 73:76; % Visual-only conditions (4)
    M_n = length(idx);
    var_v = repmat(VarV(idx), N, 1);
    var_p = repmat(Varp, N, M_n);
    mu_p = repmat(mup, N, M_n);
    S_v = repmat(Data(idx, 3)',N,1);
    x_v = bsxfun(@times,seed,sqrt(var_v)) + S_v;
    Sest_V_uni = (x_v./var_v + mu_p./var_p)./(1./var_v+1./var_p);
        
    cond_sel = 1:76;
    
    case 'mul'
    Sest_A_uni = [];
    Sest_V_uni = [];
    cond_sel = 1:64;
    end
    
    
Sest = [Sest_A(:,1:32) Sest_V(:,33:64) Sest_A_uni Sest_V_uni]; % 10000 * 76

% Binning predicted samples into 4 bins
% then, predict response probabilities:
% x = Sest';
% level = unique(Data(1:16,3)); % four true rate levels
% xc = mat2cell(x, ones(1, size(x, 1)), size(x, 2));
% edge = [-inf, mean(level(1:2)), mean(level(2:3)), mean(level(3:4)), inf];
% hcell = cellfun(@(x) histcounts(x, edge), xc, 'Uni', 0);
% hmtx = cell2mat(hcell);
% hmtx = hmtx + ep;
% p_pred = hmtx/(N + 4*ep); % model predicted response probabilities [76 * 4]
% 
% error = -sum(sum(Data(:, 5:8).*log(p_pred))); % -Log-likelihood

level = unique(Data(1:16, 3)); % four true rate levels
edge = [-inf, mean(level(1:2)), mean(level(2:3)), mean(level(3:4)), inf];

% Binning predicted samples into 4 bins
% then, predict response probabilities:
% tic
% x = Sest';
% xc = mat2cell(x, ones(1, size(x, 1)), size(x, 2));
% hcell = cellfun(@(x) histcounts(x, edge), xc, 'Uni', 0);
% hmtx = cell2mat(hcell);
% t1 = toc;

[nx1, nx2] = size(Sest);
s1 = sum(Sest < edge(2));
s2 = sum(Sest < edge(3));
s3 = sum(Sest < edge(4));
s4 = nx1*ones(1,nx2);
hmtx = diff(cat(1,zeros(1,nx2),s1,s2,s3,s4))';

hmtx = hmtx + ep;
p_pred = hmtx/(N + 4*ep); % model predicted response probabilities [76 * 4]

% size(Data(attention_idx, 5:8))
% 
% size(p_pred)


error = -sum(sum(Data(cond_sel, 5:8).*log(p_pred))); % -Log-likelihood
Ntot = sum(Data(cond_sel , end));


error_A = -sum(sum(Data(1:32, 5:8).*log(p_pred(1:32, :))));
error_V = -sum(sum(Data(33:64, 5:8).*log(p_pred(33:64, :))));

out.Ntot = Ntot;
out.error_A = error_A;
out.error_V = error_V;
out.p_pred = p_pred;
out.Sest = mean(Sest)';
CI=out.Sest;

% plot(CI,'ko-'); axis tight

pred_rate_all(1,:,whichpc)=CI;


out.p_AVc1 = p_AVc1;
out.p_AVc2 = p_AVc2;
out.p_c1AV = p_c1AV;

out.Sest_raw = Sest;
out.Sest_AVc1 = mean(Sest_AVc1);
out.Sest_Vc2 = mean(Sest_Vc2);
out.Sest_Ac2 = mean(Sest_Ac2);
out.negLL = NaN;


if do_nonCausalPred
    Seg=[Sest_Ac2(:,1:32), Sest_Vc2(:,33:64)];
    Ssum=cat(3,Sest_AVc1,Seg);
    for kk=1:2
        Sest=Ssum(:,:,kk);
        [nx1, nx2] = size(Sest);
        s1 = sum(Sest < edge(2));
        s2 = sum(Sest < edge(3));
        s3 = sum(Sest < edge(4));
        s4 = nx1*ones(1,nx2);
        hmtx = diff(cat(1,zeros(1,nx2),s1,s2,s3,s4))';
        hmtx = hmtx + ep;
        p_pred = hmtx/(N + 4*ep);
        negLL(kk) = -sum(sum(Data(cond_sel, 5:8).*log(p_pred)));
        out.negLL = negLL;
    end

end

yheight=0.2;

for rel=[0,16]
whichc=4+32 + rel; % visual task
CI=Sest(:,whichc);
V=Sest_Vc2(:,whichc);
A=Sest_Ac2(:,whichc);

X=[A,V,CI];
if whichpc==1 % fusion
    if rel==0
        Fusion(:,1)=CI;
    else
        Fusion(:,2)=CI;
    end
end

if whichpc==1
    if rel==0
        f1=figure('position',[381,509,298,184],'Name',num2str(whichpc));
    else
        figure(f1)
    end
    ccc=[0.7,0.7,0.7;0.7,0.7,0.7;0,0,1];
    for i=1:3
        [f,xi]=ksdensity(X(:,i));
        if i==2
            plot(xi,f,'color',ccc(i,:),'linewidth',5);
        else
            if rel==0
                 plot(xi,f,'color',ccc(i,:),'linewidth',2);
            else
                 plot(xi,f,'--','color',ccc(i,:),'linewidth',2);
            end
        end
        hold on;
        plot([-10,30],[0,0],'k-','linewidth',2)
        xlim([-10,30])
        ylim([0,yheight])
%         axis tight
    end
%     f1=get(f1,'children');
    set(gca,'box','off','linewidth',2,'xtick',[],'ytick',[])
elseif whichpc==2
    if rel==0
        f2=figure('position',[381,509,298,184],'Name',num2str(whichpc));
    else
        figure(f2)
    end
%     X=[X,Fusion];
    ccc=[0.7,0.7,0.7;0.7,0.7,0.7;1,0,0;0,0,1];
    for i=1:3
        [f,xi]=ksdensity(X(:,i));
        mm=mean(X);
        if i==2
            plot(xi,f,'color',ccc(i,:),'linewidth',5);
        else
            if rel==0
                 plot(xi,f,'color',ccc(i,:),'linewidth',2);
            else
                 plot(xi,f,'--','color',ccc(i,:),'linewidth',2);
            end
        end
        hold on;
%         plot(mm,0,'go-')
        plot([-10,30],[0,0],'k-','linewidth',2)
        xlim([-10,30])
        ylim([0,yheight])
    end
    set(gca,'box','off','linewidth',2,'xtick',[],'ytick',[])
elseif whichpc==3
    if rel==0
        f3=figure('position',[381,509,298,184],'Name',num2str(whichpc));
    else
        figure(f3)
    end
    ccc=[0.7,0.7,0.7;0.7,0.7,0.7;0,1,0];
    for i=1:3
        [f,xi]=ksdensity(X(:,i));
        if i==3
           plot(xi,f,'color',ccc(i,:),'linewidth',5);
        else
            if rel==0
                 plot(xi,f,'color',ccc(i,:),'linewidth',2);
            else
                 plot(xi,f,'--','color',ccc(i,:),'linewidth',2);
            end
        end
        hold on;
        plot([-10,30],[0,0],'k-','linewidth',2)
        xlim([-10,30])
       ylim([0,yheight])
    end
    set(gca,'box','off','linewidth',2,'xtick',[],'ytick',[])
    
end
end

end % end of pc vec



% A:fusion, a, v grey, fused in blue
% B:
% C:seg, 


% histogram(A); hold on;
% histogram(V); hold on;
% histogram(CI); hold on;

figure(f2)
ccc=[0.7,0.7,0.7;0.7,0.7,0.7;1,0,0;0,0,1];
for i=1:2
    [f,xi]=ksdensity(Fusion(:,i));
    if i==1
         plot(xi,f,'color',[0,0,1],'linewidth',2);
    else
         plot(xi,f,'--','color',[0,0,1],'linewidth',2);
    end
    hold on;
    xlim([-10,30])
    ylim([0,yheight])
end



clc; 

all_model=repmat(pred_rate_all,[15,1,1]);
size(all_model)

all_model;



% all_model=pred_rate_all;
% % tmp=all_model';
% for sub=1:15
%     length(unique(tmp(:,sub)))
% end

% all_model=cat(3,pred_rate_all(:,:,1),SEG);
% D = Info.Data(1:64, 1:4);
% open D

D=[
         0    1.0000    9.0909    9.0909
         0    1.0000    9.0909   12.7273
         0    1.0000    9.0909   16.3636
         0    1.0000    9.0909   20.0000
         0    1.0000   12.7273    9.0909
         0    1.0000   12.7273   12.7273
         0    1.0000   12.7273   16.3636
         0    1.0000   12.7273   20.0000
         0    1.0000   16.3636    9.0909
         0    1.0000   16.3636   12.7273
         0    1.0000   16.3636   16.3636
         0    1.0000   16.3636   20.0000
         0    1.0000   20.0000    9.0909
         0    1.0000   20.0000   12.7273
         0    1.0000   20.0000   16.3636
         0    1.0000   20.0000   20.0000
         0    2.0000    9.0909    9.0909
         0    2.0000    9.0909   12.7273
         0    2.0000    9.0909   16.3636
         0    2.0000    9.0909   20.0000
         0    2.0000   12.7273    9.0909
         0    2.0000   12.7273   12.7273
         0    2.0000   12.7273   16.3636
         0    2.0000   12.7273   20.0000
         0    2.0000   16.3636    9.0909
         0    2.0000   16.3636   12.7273
         0    2.0000   16.3636   16.3636
         0    2.0000   16.3636   20.0000
         0    2.0000   20.0000    9.0909
         0    2.0000   20.0000   12.7273
         0    2.0000   20.0000   16.3636
         0    2.0000   20.0000   20.0000
    1.0000    1.0000    9.0909    9.0909
    1.0000    1.0000    9.0909   12.7273
    1.0000    1.0000    9.0909   16.3636
    1.0000    1.0000    9.0909   20.0000
    1.0000    1.0000   12.7273    9.0909
    1.0000    1.0000   12.7273   12.7273
    1.0000    1.0000   12.7273   16.3636
    1.0000    1.0000   12.7273   20.0000
    1.0000    1.0000   16.3636    9.0909
    1.0000    1.0000   16.3636   12.7273
    1.0000    1.0000   16.3636   16.3636
    1.0000    1.0000   16.3636   20.0000
    1.0000    1.0000   20.0000    9.0909
    1.0000    1.0000   20.0000   12.7273
    1.0000    1.0000   20.0000   16.3636
    1.0000    1.0000   20.0000   20.0000
    1.0000    2.0000    9.0909    9.0909
    1.0000    2.0000    9.0909   12.7273
    1.0000    2.0000    9.0909   16.3636
    1.0000    2.0000    9.0909   20.0000
    1.0000    2.0000   12.7273    9.0909
    1.0000    2.0000   12.7273   12.7273
    1.0000    2.0000   12.7273   16.3636
    1.0000    2.0000   12.7273   20.0000
    1.0000    2.0000   16.3636    9.0909
    1.0000    2.0000   16.3636   12.7273
    1.0000    2.0000   16.3636   16.3636
    1.0000    2.0000   16.3636   20.0000
    1.0000    2.0000   20.0000    9.0909
    1.0000    2.0000   20.0000   12.7273
    1.0000    2.0000   20.0000   16.3636
    1.0000    2.0000   20.0000   20.0000];

cond=[];
for task=[0,1]
    for ar=1:2
        for v=1:4
            for a=1:4
                cond=[cond;task,ar,v,a];
            end
        end
    end
end
visrate=cond(:,3);
audrate=cond(:,4);

docorrect=1
h = figure('Position',[391   397   300   300]);
hold;
for whichmodel = [1:3]
    pred = all_model(:,:,whichmodel);
    BB = [];
    for whichsub = 1:15
        clear thispred
        thispred = [D, pred(whichsub, :)']; % predicted resp
        AVd = thispred(:, 3) - thispred(:, 4);
        cong=thispred(thispred(:,3)==thispred(:,4),:);
        XX=[];
        ss=[]; 
        for whichtask=[0,1]
            for whichar=1:2
            ind=cong(:,1)==whichtask&cong(:,2)==whichar;
            DV=cong(ind,5);
            physical=cong(ind,4); py=unique(physical);
            predR=[ones(size(DV)),cong(ind,4),cong(ind,4).^2]; % physical
            [B,~,R]=regress(DV,predR);
            tmptmp=DV-R; % regressed rates (corrected physical);
            newr=[];
            for i=1:4
                newr=[newr,tmptmp(physical==py(i))];
            end
            percrate=newr(1,:);
            
            DVreg=newr(1,:)';
%             DVreg=DV;
            rep_idx=cond(:,1)==whichtask&cond(:,2)==whichar;
            
            repRate=[DVreg(visrate(rep_idx)),DVreg(audrate(rep_idx))];
            ss=[repmat([whichtask,whichar],size(repRate,1),1),repRate];
            XX=[XX;ss];
            end
        end

        RegRate=XX(:,3:4); % corrected physical rates. Vis, Aud
        
        if docorrect
            thispred(:,3:4)=RegRate;
        end
        
        thispred = [thispred, zeros(64,1)];
        
        thispred(thispred(:,1)==0, 6) = thispred(thispred(:,1)==0, 4);
        thispred(thispred(:,1)==1, 6) = thispred(thispred(:,1)==1, 3); % target signal
        thispred(thispred(:,1)==0, 10) = thispred(thispred(:,1)==0, 3);
        thispred(thispred(:,1)==1, 10) = thispred(thispred(:,1)==1, 4); 
        thispred(:, 7) = thispred(:, 5) - thispred(:, 6); %bias
        
%         thispred(:,7)=thispred(:,7)./(thispred(:,10)-thispred(:,6));
        
%         AVd = thispred(:, 3) - thispred(:, 4);
       
        % group data:
        AVd=round(AVd,3);
        uAVd = unique(AVd)
        thisBias = [];
        for task = [0,1]
            for aRel = 1:2
                bias = [];
                for i = 1:length(uAVd)
                    bias(1,i) = nanmean(thispred(AVd==uAVd(i)&thispred(:,1)==task&thispred(:,2)==aRel, 7));
                end
                thisBias = [thisBias; task,aRel,bias];
                
            end
        end
        BB(:,:,whichsub) = thisBias;
    end
    M = mean(BB, 3);
    SE = squeeze(std(permute(BB, [3,1,2])))./sqrt(15);
    
    SE(:,1:2) = M(:,1:2);
    ccc = [1,0,0; 0,0,0];
    cmodel={'b';'r';'g'};
    scale = 2/0.55; 
    pp=[1,1];
    for task = [1]
        for AR = 1:2
            idx = M(:,1)==task&M(:,2)==AR;
        if AR == 1
            plot((-3:3)*pp(task+1)*scale,M(idx,3:end),'-o','linewidth',2.4,...
               'color',cmodel{whichmodel},'markersize',7,'markerfacecolor',cmodel{whichmodel})
%            errorbar((-3:3)*pp(task+1)*scale,M(idx,3:end),SE(idx,3:end),'o-','linewidth',2.4,...
%                'color',ccc(task+1,:),'markersize',7,'markerfacecolor',ccc(task+1,:))
        else
%            errorbar((-3:3)*pp(task+1)*scale,M(idx,3:end),SE(idx,3:end),'s--',...
%                'color',ccc(task+1,:),'markersize',9,'linewidth',1.5);
plot((-3:3)*pp(task+1)*scale,M(idx,3:end),'s--',...
               'color',cmodel{whichmodel},'markersize',9,'linewidth',1.5);
        end
        hold on;
        end
    end
% plot([-40,40], [0,0],'k--')
hold on;

end
% xlim([-3.5,3.5])
ylim([-2.8,2.8])
% axis tight
xlabel('Disparity')
ylabel({'Crossmodal bias'})
axis square
xlim([-12,12])
leg = {'A report, high A reliability';
       'A report, low  A reliability';
       'V report, high A reliability';
       'V report, low  A reliability';
       };
% legend(leg, 'Location', 'northoutside')
h = get(h, 'children');
set(h, 'linewidth', 2, 'box', 'off', 'fontsize', 15,'xtick',[],'ytick',[])
% ylim(range)

range=[-2.5,2];


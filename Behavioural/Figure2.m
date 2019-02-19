%% perceptual bias
%%
clear; clc; close all;
load('Data_15subjs_22Trls_MEGextract.mat')
addpath('YC_private')
docorrect=1;
g_bias = [];
g_rt = [];

winS=[392   389   300   308];

subj_bias=[];
for whichsub = 1:size(Data_runs,3)
    D = Data_runs(:,:,whichsub);
    multi_idx = D(:,3)<5 & D(:,4)<5;
    uniD=D(~multi_idx,:);
    for task=0
        for arel=1:2
            s=uniD(uniD(:,1)==task&uniD(:,2)==arel,3:5);
            b=regress(s(:,3),[ones(size(s,1),1),s(:,2),s(:,2).^2]);
            pred=[ones(size(s,1),1),s(:,2),s(:,2).^2]*b;
            Auni(:,arel)=unique(pred);
            
%             for i=1:4
%                 Auni(i,arel)=mean(s(s(:,2)==i,3));
%             end
            
        end
    end
    
    for task=1
        for arel=3
            s=uniD(uniD(:,1)==task&uniD(:,2)==arel,3:5);
            b=regress(s(:,3),[ones(size(s,1),1),s(:,1),s(:,1).^2]);
            pred=[ones(size(s,1),1),s(:,1),s(:,1).^2]*b;
            Vuni(:,1)=unique(pred);
%             for i=1:4
%                 Vuni(i,1)=mean(s(s(:,1)==i,3));
%             end
        end
    end
    

%     open D % task,audR,vis,aud,rsp,RT
    D = D(multi_idx,:); % select only multisensory trials
    task = D(:,1); AR = D(:,2);
    
    % congruent trials, irrespective of reliability
    cong=D(D(:,3)==D(:,4),:);
    % regression:
    
    for whichtask=[0,1]
        for whichar=1:2
        ind=cong(:,1)==whichtask&cong(:,2)==whichar;
        DV=cong(ind,5); % response
        physical=cong(ind,4);
        percrate=[];
        for i=1:4
            percrate(i,1)=mean(DV(physical==i));
        end
       
        DVreg(:,whichtask+1,whichar)=percrate;
        end
    end
    
    Disp = (D(:,3)-D(:,4)); % Vis - Aud, physical
    disp_level = unique(Disp)';
%     disp_class{1} = find(Disp<2);
%     disp_class{2} = find(Disp>=2);
    
    all_bias = []; all_rt = [];
     for whichtask = [0,1] % aud vis
      for whichar = 1:2    % high low
          bias_cond = []; rt_cond = [];
          
          aud_r=Auni(D(task==whichtask&AR==whichar,4),whichar);
          vis_r=Vuni(D(task==whichtask&AR==whichar,3));
              
        for whichdisp = disp_level
           S = D(Disp==whichdisp & task==whichtask & AR==whichar,3:6);
           
          tmp_rate=S(:,1:2);
           
          aud_r=Auni(D(task==whichtask&AR==whichar&Disp==whichdisp,4),whichar);
          vis_r=Vuni(D(task==whichtask&AR==whichar&Disp==whichdisp,3));
          
           
           if whichtask == 0 % aud attention
               phy=S(:,2);
               if docorrect
               thisDVreg=DVreg(:,whichtask+1,whichar);
               phy=thisDVreg(phy);
               VA=thisDVreg(S(:,1))-thisDVreg(S(:,2));
%                VA(VA==0)=NaN;
               end
               bias = S(:,3)-phy;
%                bias = bias./VA;
           else
               phy=S(:,1);
               if docorrect
               thisDVreg=DVreg(:,whichtask+1,whichar);
               phy=thisDVreg(phy);
               VA=thisDVreg(S(:,2))-thisDVreg(S(:,1));
%                VA(VA==0)=NaN;
               end
               bias = S(:,3)-phy;
%                bias = bias./VA;
           end
           size([aud_r,vis_r,bias]);
           n=size(bias,1);
           
               
           x=[whichsub*ones(n,1),...
               whichtask*ones(n,1),whichar*ones(n,1),...
               whichdisp*ones(n,1),abs(aud_r-vis_r),bias,abs(VA),tmp_rate,S(:,4),log(S(:,4))];
%            x=mean(x,1);
           subj_bias=[subj_bias;x];
               
%            if whichsub==4
%                bias=NaN;
%            end
           bias_cond = [bias_cond, nanmean(bias)];
%            rt_cond = [rt_cond, mean(S(:,4))];
        end
        all_bias = [all_bias; [whichtask,whichar,bias_cond]];
%         all_rt = [all_rt; [whichtask,whichar,rt_cond]];
      end
     end
     g_bias(:,:,whichsub) = all_bias;
%      g_rt(:,:,whichsub) = all_rt;
end

Beta=[];  RT_stack=[];
BehavModelFree=[];
for whichsub=1:15
    Y=subj_bias(subj_bias(:,1)==whichsub,:);
    
    Ys=sortrows(Y,[2,3,8,9]);
    tmp=YC_transpose(Ys,22);
    tmp2=squeeze(mean(tmp))';
    Y=tmp2;
    RT_stack(:,:,whichsub)=Y;
    
    measure=1;
    if measure==1
        thisbias=abs(Y(:,6));
    elseif measure==2
        thisbias=Y(:,10); % raw RT
    else
        thisbias=Y(:,11); % log RT
    end
    
    task=Y(:,2);
    Arel=Y(:,3);
    uni=0;
    if uni
        AVd=abs(Y(:,5));
    else
        AVd=abs(Y(:,7));
    end
    task=zscore(task);
    Arel=zscore(Arel);
    AVd=zscore(AVd);
    AVd2=AVd.^2;
    AVd2=zscore(AVd2);
    
    AVd3=zscore(AVd.^3);
    AVd4=zscore(AVd.^4);
    
    doFull=0;
    if doFull
    Xreg=[task,Arel,task.*Arel,...
        AVd,AVd2,...
        task.*AVd,task.*AVd2,...
        Arel.*AVd,Arel.*AVd2,...
        task.*Arel.*AVd,...
        task.*Arel.*AVd2
        ];
%         Xreg=[AVd,AVd2,...
%         ];
    
%      Xreg=[task,Arel,task.*Arel,...
%          AVd,AVd2,...
%         ];
    
    else
    Xreg=[task,Arel,task.*Arel,AVd,AVd2];
    end
    
    Nreg=size(Xreg,2);
    
    Yreg=thisbias;
    Yreg=zscore(Yreg);
    Xreg=zscore(Xreg);
    Beta(whichsub,:)=regress(Yreg,[ones(size(Yreg)),(Xreg)]);
    
    tmp=[Yreg,Xreg];
    size(tmp)
    BehavModelFree(:,:,whichsub)=tmp;
    
    
    
    b=glmfit(Xreg,Yreg);
    CMbiasBs(whichsub,:)=b(1:end)';
end


Beta=CMbiasBs

T=BLG_Ttest(Beta)

[h,p]=ttest(Beta)

d=mean(Beta)./std(Beta)

out=[T',h',p',d'];

if doFull
    varName={'int','task','Arel','task.*Arel',...
            'AVd','AVd2',...
            'task.*AVd','task.*AVd2',...
            'Arel.*AVd','Arel.*AVd2',...
            'task.*Arel.*AVd',...
            'task.*Arel.*AVd2'};
else
    varName={'int','task','Arel','task.*Arel',...
        'AVd','AVd2'};
end
    
table(varName',out(:,1),out(:,3),out(:,4),mean(Beta)',std(Beta)'./sqrt(15),...
    'VariableNames',{'effect','T14','p','Cohend','m','se'})

%        effect         T14          p         Cohend 
%     ____________    _______    __________    _______
% 
%     'int'            2.6632      0.018545    0.68763
%     'task'           1.8367       0.08757    0.47424
%     'Arel'           1.9108      0.076725    0.49337
%     'task.*Arel'    -6.3643    1.7558e-05    -1.6433
%     'AVd'            22.456    2.2253e-12      5.798
%     'AVd2'          -9.2767    2.3479e-07    -2.3952


%% Permutation test maximum-statistics
nu=14;
T=BLG_Ttest(CMbiasBs,14);
p=(1-tcdf(abs(T),nu))*2; %two sided parametric pvalue
p=p';
T=T';
BetaP=p;
BetaT=T;

CMbiasDat=repmat(BehavModelFree,[1,1,1,2]);
DatPerm=CMbiasDat; %nobjs nvars nsubjs nrois 
          % for behaviour you will nrois=1. As simple as that

          %%% this is how input data needs to be structured for behaviour
          %%% analysed with this code

[N,nvars,nsubjs,nrois]=size(CMbiasDat); %N = 64 conditions; nvars is the zscored regressors (including the intercept = first regressor)
                                        %nrois = 1 for behavior
DatPerm=permute(DatPerm,[1 2 4 3]);%nobjs nvars nrois nsubjs
X=DatPerm;
X(:,1,:,:)=ones(size(X(:,1,:,:))); % regressors
Y=DatPerm(:,1,:,:);
 
nperms=10000;
% N=64;

%%% a different, equivalent way to do permutation used here (matrix
%%% multiplication approach)
%same set of permutations for each coefficient for all ROIs, but different
%across subjects.
[~,perms]=sort(rand([N,nsubjs,nperms-1]));
perms2=(perms-1)*N+repmat([1:N]',1,nsubjs,nperms-1);
PPs=zeros([N N nsubjs nperms-1]);
z=zeros(N);
for whichsubj=1:nsubjs
    for j=1:nperms-1
        P=z;
        P(perms2(:,whichsubj,j))=1;
        PPs(:,:,whichsubj,j)=P; %permutation matrices
    end
end

r=nvars-1;%size(XRegr,2);
NY=nrois;
Nsubjs=nsubjs;
whichX=[2 3 4 5 6]; %which GLM coefficients are to be tested with permutations


sX=size(X);
X=reshape(X,[sX(1:2) prod(sX(3:4))]);
Mplus=BLG_pinv3D(X,0);
Mplus=reshape(Mplus,[sX(2) sX(1) sX(3:4)]);%Moore-Penrose pseudo inverse of design matrix
PlugBetas=mtimesx(Mplus,'n',Y,'n'); %regression coefficients
PlugBetas=permute(PlugBetas(2:end,:,:,:),[4 1 3 2]); %same as CMbiasBs

%close all %sanity check: same as CMbiasBs
%plot(CMbiasBs(:),CMbiasBs(:)-PlugBetas(:),'.')


PermBetas=[];

for isubj=1:nsubjs
    thisperms=PPs(:,:,isubj,:);
    thisYperm=mtimesx(thisperms,Y(:,:,:,isubj)); %another way to do permutation
               %through matrix multiplication
    thisBetas=mtimesx(Mplus(:,:,:,isubj),'n',thisYperm);
    PermBetas(isubj,:,:,:)=permute(thisBetas(2:end,:,:,:),[2 1 3 4]);
    disp(num2str(isubj))
end
PermBetas=cat(4,PermBetas,PlugBetas);
PermTs=BLG_Ttest(PermBetas,nsubjs-1);
PermTs=permute(PermTs,[2 1 3]); %rois, effects, perms


%%% FWE = 0.05 with max stat only across ROIs AND effects
% thissigperm1=zeros(size(thissig));

%%
thisTperm=abs(PermTs);
maxT=prctile(thisTperm,100,2); % across effects
Tthr=prctile(maxT,95,3)

% symmetricl, 2 sided, take abs, then 95

thissigperm=double(abs(BetaT)-Tthr(1)>=0)

% thisTperm(removerois,:,:)=NaN;
% maxstatStim=prctile(thisTperm(stimresp==1,:,:),100,2); 
% maxstatResp=prctile(thisTperm(stimresp==2,:,:),100,2);
% maxstatStim=prctile(maxstatStim,100,1); 
% maxstatResp=prctile(maxstatResp,100,1);
% TthrStim=prctile(maxstatStim,95,3);
% TthrResp=prctile(maxstatResp,95,3);

table(varName',thissigperm,out(:,1),out(:,3),out(:,4),mean(Beta)',std(Beta)'./sqrt(15),...
    'VariableNames',{'effect','sig','T14','p','Cohend','m','se'})




%% verify distribution of permutation Ts for different effects
close all
betanams={'Task' 'Reliability' 'TaskxRel' 'Linear bias' 'Nonlinear bias'};
xlims=prctile(PermTs(:),[0.05 99.95]);
for i=1:5
    tmp=PermTs(:,i,:);
    tmp=tmp(:);
    subplot(2,3,i)
    prcs=prctile(tmp,[1 50 99]);
    tmp(tmp<xlims(1)| tmp>xlims(2))=[]; %for display purposes
    hist(tmp,100)
    title([betanams(i) {num2str(prcs,'%0.2f ')}],'fontsize',10) %excellent, prctiles same for all Ts; bias almost zero!
    xlim(xlims)
    xlabel('Permutation Ts')
end


%% Figure 2F/G

close all;
scale = 2/0.55;

% % % perceptual bias: plot
h = figure('Position',[654   389   300   308]);
M = nanmean(g_bias,3);
Nsubj = size(g_bias,3);
SE = squeeze(nanstd(permute(g_bias,[3,1,2]))./sqrt(15));
SE(:,1:2) = M(:,1:2);
% ccc = [1,0,0; 0.5,0.5,0.5];
ccc = [.7,0.7,0.7; 0,0,0];
pp=[1,1];
s=[-1,1];
for task = [0,1] % aud, vis
    for AR = 1:2 %only high reliability
        idx = M(:,1)==task&M(:,2)==AR;
        
        x=[-3:3]*pp(task+1)*scale;
        x = x+0.2*s(AR);
        y=M(idx,3:end)*scale;
        
        se=SE(idx,3:end)*scale;
        w=1.8;
        hold on;
%           plot([-11,11],[0 0],'r:','linewidth',2)
        if AR == 1
%            errorbar([-3:3]*pp(task+1)*scale,M(idx,3:end)*scale,SE(idx,3:end)*scale,'o-','linewidth',2.4,...
%                'color',ccc(task+1,:),'markersize',7,'markerfacecolor',ccc(task+1,:))
           
           plot(x,y,'o-','linewidth',2,...
               'color',ccc(task+1,:),'markersize',7,'markerfacecolor',ccc(task+1,:))
           hold on;
           errorbar_YC(x,y,se,w,2,ccc(task+1,:))
           
           hold on;
           
           % fusion prediction:
%            plot(x,yfu,'-','linewidth',2,...
%                'color',ccc(task+1,:),'markersize',7,'markerfacecolor',ccc(task+1,:))
           
           
           
        else
%            errorbar([-3:3]*pp(task+1)*scale,M(idx,3:end)*scale,SE(idx,3:end)*scale,'s-',...
%                'color',ccc(task+1,:),'markersize',9,'linewidth',1.5);

                       
           hold on;
           errorbar_YC(x,y,se,w,2,ccc(task+1,:))
           
%            plot(x,yfu,'--','linewidth',2,...
%                'color',ccc(task+1,:),'markersize',7,'markerfacecolor',ccc(task+1,:))
%            
           hold on;
           plot(x,y,'s--',...
              'color',ccc(task+1,:),'markersize',12,'linewidth',1.5,'markerfacecolor','w');
          
          hold on;

           

           
        end
        hold on;
    end
end
% plot([-40,40], [0,0],'k--')
% xlim([-3.5,3.5])
xlabel('Vis. - Aud. rate (Hz)')
% ylim([-1.5,1.5])
% xlabel('Vis. - [Aud.] rate (Hz)','color',[.7,0.7,0.7])
% ylabel({'Crossmodal bias (Hz)'; '(reported - task relevant)'})
ylabel('Crossmodal bias (Hz)')
% axis square
xlim([-12,12])
leg = {'A report, high A reliability';
       'A report, low  A reliability';
       'V report, high A reliability';
       'V report, low  A reliability';
       };
leg = {'A report, high A reliability';
       'A report, low  A reliability';
   'V report, high A reliability';
       'V report, low  A reliability';
   };
% legend(leg, 'Location', 'northoutside')
h = get(h, 'children');
set(h, 'linewidth', 2, 'box', 'off', 'fontsize', 15,'tickdir','out','ytick',[-2,0,2])
range=[-2.5,2];
ylim(range)
%

%
%
%
% causal inference occurs only when the rate from the other modality is
% higher than the rate from its own modality

% when the rate from the other modality is lower than the rate from its own
% modality


% do statisticsal test:
% effect of reliability:
% absAV, task, 
h = figure('Position',[392   389   320   308]);
G=abs(g_bias);
B=[];
% for subj=1:15
%     thisg=G(:,3:end,subj)';
%     Yreg=thisg(:);
%     task=[zeros(14,1);ones(14,1)];
%     Arel=[ones(7,1);ones(7,1)*2;
%           ones(7,1);ones(7,1)*2];
%     AVd=repmat([3,2,1,0,1,2,3]',4,1);
%     task=zscore(task);
%     Arel=zscore(Arel);
%     AVd=zscore(AVd);
%     AVd2=AVd.^2;
%     AVd2=zscore(AVd2);
%     Xreg=[task,Arel,task.*Arel,AVd,AVd2];
%     Nreg=size(Xreg,2);
%     
%     Yreg=zscore(Yreg);
%     Xreg=zscore(Xreg);
%     B(subj,:)=regress(Yreg,[ones(size(Yreg)),(Xreg)]);
%     
% end

B=Beta
B=B(:,2:end); % remove intercept
M=mean(B);
SE=std(B)/sqrt(15);

% bar([1,2],M([1:2]),'facecolor','w','linewidth',2,'edgecolor','k');
hold on;
bar(1:Nreg,M,'facecolor','w','linewidth',2,'edgecolor','k');
hold on;
% errorbar(M,SE,'k.','linewidth',2)
% errorbar_YC([1:5],M,SE,0.3,2,'k')
% errorbar_YC([1:2],M(1:2),SE(1:2),0.3,2,'k')
% errorbar_YC([3,4,5],M([3,4,5]),SE(3:5),0.3,2,'k')
errorbar_YC(1:Nreg,M,SE,0.3,2,'k')
plot([-1,10],[0,0],'k-')

ylabel({'Standardised beta'})
hold on;
for k=1:Nreg
    % plot([k]-0.3+randn(15,1)*0.1,B(:,k),'ko','markerfacecolor','k','markeredgecolor','w')
    [f,xi,bw]=ksdensity(B(:,k),B(:,k));
    f=f./max(f)/7;
    rng
    for i=1:length(f)
%     plot([k]-0.3+randn*f(i),B(i,k),'ko','markerfacecolor','k','markeredgecolor','w');hold on;
    plot([k]-0.3,B(i,k),'ko','markeredgecolor','k');hold on;
    end
end
axis tight
xlim([0.5,Nreg+.5])
str={'T','AR','TxAR','Lin.D','Qua.D'};


[H,P,CI,STATS]=ttest(B,0,'tail','both');
P
P=P*Nreg

h = get(h, 'children');
set(h, 'linewidth', 2, 'box', 'off', 'fontsize', 15)
set(h,'xtick',[1:Nreg],'ytick',[-0.3,0,0.3,0.6],'tickdir','out');
xlabel('Effect')
set(h,'xticklabel',str,'xticklabelrotation',45,...
    'xtick',[1:5],'ytick',[-0.3,0,0.6],'tickdir','out');

%%
clc;


beta_log=[];error_log=[];
rt_log=[];
var_log=[];
for whichsub=1:15
    clear X D RT Y perf
    D = Data_sum(:,:,whichsub);
    X=D(D(:,3)==D(:,4),:);
    tmp=X(:,5:8);
    var=[];
    for k=1:size(tmp,1)
        out=[];
        for j=1:4
            out=[out,repmat(j,1,tmp(k,j))];
        end
        var=[var;std(out)];
    end
    multi_var=var;
    
    tmp=D(65:end,5:8); % unisensory std
    var=[];
    for k=1:size(tmp,1)
        out=[];
        for j=1:4
            out=[out,repmat(j,1,tmp(k,j))];
        end
        var=[var;std(out)];
    end
    uni_var=[var;var(9:12)];
    var_log=[var_log, (multi_var-uni_var)];
    
    
    X(:,5:8)=X(:,5:8)/22;
    p=X(:,5:8);
    pp=YC_transpose(p,4);
    px=[];
    for i=1:4
        px=[px;diag(pp(:,:,i))];
    end
    meanR=p*[1,2,3,4]';
    targ=repmat([1:4]',4,1);
    
    puni=D(65:end,5:8)/22;
    
    meanRu=puni*[1,2,3,4]';
    meanRu=[meanRu;meanRu(9:12)];
    puni=YC_transpose(puni,4);
    pxu=[];
    for i=1:size(puni,3)
        pxu=[pxu;diag(puni(:,:,i))];
    end
    vis=pxu(9:12);
    pxu=[pxu;vis];
    
    di=px-pxu;
    di=reshape(di,4,4);
    
    beta_log=[beta_log;mean(di)];
    
    error=abs(meanR-targ)-abs(meanRu-targ);
    error2=reshape(error,4,4)*2/0.55;
    error_log=[error_log;mean(error2)];
    
    doLog=1;
    D = Data_runs(:,:,whichsub);
    M=D(D(:,3)==D(:,4),1:6);
    tmp=[];
    for task=[0,1]
        for AR=[1,2]
            rt=M(M(:,1)==task&M(:,2)==AR,end);
            if doLog
               tmp=[tmp;mean(log(rt))];
            else
               tmp=[tmp;mean((rt))];
            end
        end
    end
    mul=tmp;
    UniA=D(D(:,3)==5,1:6);
    tmp=[];
    for AR=1:2
        if doLog
            tmp=[tmp;mean(log(UniA(UniA(:,2)==AR,end)))];
        else
            tmp=[tmp;mean((UniA(UniA(:,2)==AR,end)))];
        end
    end
    aud=tmp;
    UniV=D(D(:,4)==5,1:6);
    if doLog
        vis=mean(log(UniV(:,end)));
    else
        vis=mean((UniV(:,end)));
    end
    uni=[aud;vis;vis];
    d_rt=mul-uni;
    rt_log=[rt_log;d_rt'];
    
end



str={'Atask/high';'Atask/low';';Vtask/high';'Vtask/low';'mean'};

h=figure('position',[218   473   895   320]);hold

subplot(1,3,1)
plotX=error_log;
plotX=[plotX,mean(plotX,2)];
boxplot(plotX); hold on;
plot([0,5],[0,0],'r--');
hold on;
plot([1:size(plotX,2)]+0.2,plotX,'ko');
errorbar_YC([1:size(plotX,2)]-0.2,mean(plotX),std(plotX)/sqrt(15),...
    0.2,2,'r');
plot([1:size(plotX,2)]-0.2,mean(plotX),'ro','markersize',10)
title('error (Mul-Uni)')
axis tight
[~,p1,~,st1]=ttest(plotX(:,end));
[p2,~,st2]=signrank(plotX(:,end));
table(p1,p2,'variableNames',{'t','sign'})

subplot(1,3,2)
plotX=beta_log;
plotX=[plotX,mean(plotX,2)];
boxplot(plotX); hold on;
plot([0,5],[0,0],'r--');
hold on;
plot([1:size(plotX,2)]+0.2,plotX,'ko');
errorbar_YC([1:size(plotX,2)]-0.2,mean(plotX),std(plotX)/sqrt(15),...
    0.2,2,'r');
plot([1:size(plotX,2)]-0.2,mean(plotX),'ro','markersize',10)
title('% corr (Mul-Uni)')
axis tight
[~,p1,~,st1]=ttest(plotX(:,end));
[p2,~,st2]=signrank(plotX(:,end));
table(p1,p2,'variableNames',{'t','sign'})


subplot(1,3,3)
plotX=rt_log;
plotX=[plotX,mean(plotX,2)];
boxplot(plotX); hold on;
plot([0,5],[0,0],'r--');
hold on;
plot([1:size(plotX,2)]+0.2,plotX,'ko');
errorbar_YC([1:size(plotX,2)]-0.2,mean(plotX),std(plotX)/sqrt(15),...
    0.2,2,'r');
plot([1:size(plotX,2)]-0.2,mean(plotX),'ro','markersize',10)
if doLog
title('Log RT (Mul-Uni)')
else
    title('RT (Mul-Uni)')
end
axis tight
[~,p1,~,st1]=ttest(plotX(:,end));
[p2,~,st2]=signrank(plotX(:,end),[],'tail','both','alpha',0.05,'method','approximate');
table(p1,p2,'variableNames',{'t','sign'})

h=get(h,'children');
set(h,'xtick',[1:5],'xticklabel',str,'xticklabelrotation',45);


%%
close all;clc;
figure;
RTgain=mean(plotX,2);
Vargain=mean(var_log)';
y=[RTgain, Vargain*2/0.55];
% plotyy(1,mean(y(:,1)),2,mean(y(:,2)),'ko')
errorbar_YC([1],mean(y(:,1)),std(y(:,1))./sqrt(15),0.2,2,'k')

yyaxis right
errorbar_YC([2],mean(y(:,2)),std(y(:,2))./sqrt(15),0.2,2,'k')

%%
[h,p,~,st]=ttest(y,0)

%%
clc;
AV=[2,2,2,4];
A=[1,3,3,3];
V=[1,3,3,3];
R=[0,0,0,10];
X=[AV;A;V;R]

rep_rate = mean(bsxfun(@times, X, [1:4]),2)'
D2=squareform(pdist(X))
D2(1:3,end)







%% pyschometric function
% also RT , new analysis 
clc; close all;
beta_log=[];
for whichsub=1:15
    clear X D RT Y perf
    D = Data_runs(:,:,whichsub);
    task = D(:,1);
    aRel = D(:,2);
    visR = D(:,3);
    audR = D(:,4);
    
    Aonly=D(visR==5,1:6);
    taskr=Aonly(:,4);
    choice=Aonly(:,5);
    ordL=mnrfit(zscore(taskr),ordinal(-choice),'model','ordinal','link','logit');
    ordP=mnrfit(zscore(taskr),ordinal(-choice),'model','ordinal','link','probit');
    auni=[0,0,ordL(4),ordP(4)];
    
    Vonly=D(audR==5,1:6);
    taskr=Vonly(:,3);
    choice=Vonly(:,5);
    ordL=mnrfit(zscore(taskr),ordinal(-choice),'model','ordinal','link','logit');
    ordP=mnrfit(zscore(taskr),ordinal(-choice),'model','ordinal','link','probit');
    vuni=[1,0,ordL(4),ordP(4)];
    
    M=D(audR<5 & visR<5,1:6);
    d=abs(M(:,3)-M(:,4));
    
    % group by task, disparity
    reg_beta=[];
    for thistask=[0,1]
        for thisd=1:3
            if thisd==1
                X=M(M(:,1)==thistask&d==0,:);
            elseif thisd==2
                X=M(M(:,1)==thistask&d==1,:);
            else
                X=M(M(:,1)==thistask&d>=2,:);
            end
            if thistask==0
                taskr=X(:,4);
            else
                taskr=X(:,3);
            end
            choice=X(:,5);
            ordL=mnrfit(zscore(taskr),ordinal(-choice),'model','ordinal','link','logit');
            ordP=mnrfit(zscore(taskr),ordinal(-choice),'model','ordinal','link','probit');
            reg_beta=[reg_beta;...
                      thistask,thisd,ordL(4),ordP(4)];
        end
    end
    
    beta_log(whichsub,:,:)=[reg_beta; auni; vuni];
    
end

[T,p,M,SE]=BLG_Ttest(beta_log);
T
M




%%
beta_log=[];
for whichsub=1:15
    clear X D RT Y perf
    D = Data_runs(:,:,whichsub);
    task = D(:,1);
    aRel = D(:,2);
    visR = D(:,3);
    audR = D(:,4);
    
    Aonly=D(visR==5,1:6);
    Vonly=D(audR==5,1:6);
    taskr=[Aonly(:,4);Vonly(:,3)];
    choice=[Aonly(:,5);Vonly(:,5)];
    ordL=mnrfit(zscore(taskr),ordinal(-choice),'model','ordinal','link','logit');
    ordP=mnrfit(zscore(taskr),ordinal(-choice),'model','ordinal','link','probit');
% ordP=[1,1,1,corr(zscore(tiedrank(taskr)), zscore(tiedrank(choice)))];

accd=ordP(4);
accd=mean(abs(taskr-choice));
% accd=mean(taskr==choice);

    uni=[NaN,ordL(4),accd];
    
    uniRT=nanmean(log([Aonly(:,6);Vonly(:,6)]));
    
    
    
    M=D(audR<5 & visR<5,1:6);
    d=abs(M(:,3)-M(:,4));
    
    % group by task, disparity
    reg_beta=[];
    
    for thisd=1:3
        taskr_2=[];
        choice_2=[];
        mRT=[];
        for thistask=[0,1]
            if thisd==1
                X=M(M(:,1)==thistask&d==0,:);
            elseif thisd==2
                X=M(M(:,1)==thistask&d==1,:);
            else
                X=M(M(:,1)==thistask&d>=2,:);
            end
            if thistask==0
                taskr=X(:,4);
            else
                taskr=X(:,3);
            end
            choice=X(:,5);
            
            taskr_2=[taskr_2;taskr];
            choice_2=[choice_2;choice];
            mRT=[mRT;log(X(:,6))];
        end
        
        
        ordL=mnrfit(zscore(taskr_2),ordinal(-choice_2),'model','ordinal','link','logit');
        ordP=mnrfit(zscore(taskr_2),ordinal(-choice_2),'model','ordinal','link','probit');
%         ordP=[1,1,1,corr(zscore(tiedrank(taskr_2)), zscore(tiedrank(choice_2)))];
accd=ordP(4);
accd=mean(abs(taskr_2-choice_2));

% accd=mean(taskr_2==choice_2);

        reg_beta=[reg_beta;...
                  thisd,ordL(4),accd];
              rt_all(whichsub,thisd)=nanmean(mRT)-uniRT;
    end
    
    beta_log(whichsub,:,:)=bsxfun(@minus,reg_beta,uni);
    
end



%
% close all; clc;
%
% h=figure('Position',[392   389   600   308]);
close all;
h = figure('Position',[567   517   342   276]);
hold
[T,p,M,SE]=BLG_Ttest(rt_all);
d=M./(SE*sqrt(15))
[~,p]=ttest(rt_all)
subplot(1,2,1)
x=1:3;
y=M;
errorbar_YC(1-0.2,y(1),SE(1),.3,2,'k');
hold on;
boxplot(rt_all(:,1));
plot(1-0.2,y(1),'ko-','markersize',9,'markerfacecolor','w','linewidth',2); hold on
plot(1+0.2,rt_all(:,1),'ko')


hold on;
plot([-1,5],[0,0],'k--')
% ylabel({'Multisensory Log RT (s)' 'relative to Unisensory'})
axis tight
% ylim([-0.04, 0.04])
xlim([0.5,1.5])
% xlabel('Abs(audiovisual disparity)')
hold on;


%let's plot accraucy
% subj,AVd,[avdlevel,corr], type
beta_log=beta_log*(2/0.55);
[T,p,M,SE]=BLG_Ttest(beta_log);
subplot(1,2,2)
x=1:3;
y=M(:,3);
d=M./(SE*sqrt(15))
[~,p,~,st]=ttest(beta_log(:,:,3))

[p,h,sts]=signrank(squeeze(beta_log(:,1,3)))
% errorbar_YC(x,y,SE(:,3),.3,2,'k');
errorbar_YC(1-0.2,y(1),SE(1,3),.3,2,'k');
boxplot(beta_log(:,1,3))
hold on;
plot(1-0.2,y(1),'ko-','markersize',9,'markerfacecolor','w','linewidth',2)
hold on;
plot(1+0.2,beta_log(:,1,3),'ko')
hold on;
plot([-1,5],[0,0],'k--')
% ylabel({'Multisensory accuracy' 'relative to Unisensory'})
axis tight
xlim([0.5,1.5])
% ylim([-0.01,0.03])

absAVd=2/0.55*[0,1,2]

h=get(h,'children');
set(h,'fontsize',15,'linewidth',2,'box','off','tickdir','out',...
    'xtick',[1],'xticklabel',{},'xticklabelrotation',20)
% set(h,'fontsize',15,'linewidth',2,'box','off','tickdir','out',...
%     'xtick',[1,2,3],'xticklabel',{'0';'3.6';['>=7.3']},'xticklabelrotation',20)

%%


%
for whichsub = 1:15
    clear X D RT Y perf
    D = Data_runs(:,:,whichsub);
    task = D(:,1);
    aRel = D(:,2);
    visR = D(:,3);
    audR = D(:,4);
    
    A_uni_high = D(task==0 & aRel==1 & visR==5, [4, 5, 6]);
    A_uni_low  = D(task==0 & aRel==2 & visR==5, [4, 5, 6]);
    V_uni      = D(task==1 & audR==5, [3, 5, 6]); % target, resp, RT
    
    A_mul_high = D(task==0 & aRel==1 & visR==audR,  [4, 5, 6]);
    A_mul_low  = D(task==0 & aRel==2 & visR==audR,  [4, 5, 6]);
    V_mul_high = D(task==1 & aRel==1 & visR==audR,  [3, 5, 6]);
    V_mul_low  = D(task==1 & aRel==2 & visR==audR,  [3, 5, 6]); % congruent
    
    X = {A_uni_high,A_uni_low,V_uni,A_mul_high,A_mul_low,V_mul_high,V_mul_low};
    for i = 1:length(X)
        Y = X{i};
        perf(i, 1) = corr(Y(:,1), Y(:,2), 'Type', 'Pearson');
        RT(i, :) = [median(Y(:,3)), mean(Y(:,3)), median(log(Y(:,3))), mean(log(Y(:,3))), mean(log10(Y(:,3)))];
    end
%     perf;
    % the contrast:
    % 1 vs. 4  uniAhigh-mulAhigh
    % 2 vs. 5  uniAlow-mulAlow
    % 3 vs. 6  uniV - mulV_ahigh
    % 3 vs. 7  uniV - mulV_alow

    conPerf(:, whichsub) = mean([perf(1)-perf(4); perf(2)-perf(5); perf(3)-perf(6); perf(3)-perf(7)]);
    conRT(:,:,whichsub)  = [RT(1,:)-RT(4,:); RT(2,:)-RT(5,:); RT(3,:)-RT(6,:); RT(3,:)-RT(7,:)];
    
    Uni = D(visR==5 | audR==5, [3, 4, 5, 6, 1]);
    Uni(Uni==5) = 0;
    uniTarg = Uni(:,1) + Uni(:,2);
    Uni = [uniTarg, Uni(:,3:5)];
    
    
    
    
    uniRT=Uni;
    aT=Uni(:,4)==0;
    vT=Uni(:,4)==1;
    Uni_subj(whichsub,:,:)=[[mean([uniRT(aT,3),log(uniRT(aT,3))]),median(uniRT(aT,3))];
    [mean([uniRT(vT,3),log(uniRT(vT,3))]),median(uniRT(vT,3))]];
    % 12*22=264 trials in total for this subjec
    
    Multi=D(audR<5 & visR<5,1:6);
    d=abs(Multi(:,3)-Multi(:,4));
    RTall=[];
    acc=[];
    
    for task=[0,1]
      for AR=[1,2]
        for AVd=[0,1,2,3]
            res=Multi(Multi(:,1)==task&Multi(:,2)==AR&d==AVd, [5:6]);
            sss=Multi(Multi(:,1)==task&Multi(:,2)==AR&d==AVd, [1,2,3,4,5]);
            if task==0
                targ=sss(:,[4,5]);
            else
                targ=sss(:,[3,5]);
            end
            
            useCorr=0;
            if useCorr
                acc=[acc; task,AR,AVd,corr(targ(:,1),targ(:,2))];
            else
                tmp=abs(targ(:,1)-targ(:,2));
                acc=[acc; task,AR,AVd,length(find(tmp==0))./length(tmp)];
%                 acc=[acc; task,AR,AVd, abs(mean(targ(:,1)-targ(:,2)))];
            end
            
            RT=[res(:,2), log(res(:,2))];
            RTall=[RTall; task,AR,AVd,mean(RT), median(res(:,2))];
        end
      end
    end
    
    acc2=[];
    for task=[0,1]
        for AVd=[0,1,2,3]
            res=Multi(Multi(:,1)==task& d==AVd, [5:6]);
            sss=Multi(Multi(:,1)==task& d==AVd, [1,2,3,4,5]);
            if task==0
                targ=sss(:,[4,5]);
            else
                targ=sss(:,[3,5]);
            end
            acc2=[acc2; task,AVd,corr(targ(:,1),targ(:,2))];
        end
    end
    
    acc3=[];
        for AVd=[0,1,2,3]
            res=Multi(d==AVd, [5:6]);
            sss=Multi(d==AVd, [1,2,3,4,5]);

            targ=[sss(sss(:,1)==0,[4,5]);sss(sss(:,1)==1,[3,5])];

            acc3=[acc3; AVd,corr(targ(:,1),targ(:,2))];
        end
    
    
    
    
    
    % average for acc:
    acc=sortrows(acc,3);
    acc1=[];
    for i=[0,1,2,3]
        acc1=[acc1; i,mean(acc(acc(:,3)==i,end))];
    end
    
    % accruacy measure, pearson corelaiton, varies with disparity.
    if useCorr
        accuracy(whichsub,:,:,1)=[-1,corr(Uni(:,1),Uni(:,2)); acc3];
        accuracy(whichsub,:,:,2)=[-1,corr(Uni(:,1),Uni(:,2)); acc1];
    else
        tmp=abs(Uni(:,1)-Uni(:,2));
        tmp=length(find(tmp==0))./length(tmp);
        accuracy(whichsub,:,:,1)=[-1,tmp; acc3];
        accuracy(whichsub,:,:,2)=[-1,tmp; acc1];

    end
    
    
    RTall_subj(whichsub,:,:)=RTall;
    
%         whichArel = 2; % in Vtask, choose high arel, or low arel, to match N across uni vs. mul
%     Mul = D(visR==audR & ( (task==0) | (task==1&aRel==whichArel) ), [4, 5, 6]);
%     
%     corrStr = 'Spearman';
% %     corrStr = 'Pearson';
%     tmp(whichsub, 1) = corr(Mul(:,1), Mul(:,2),'Type',corrStr) - corr(Uni(:,1), Uni(:,2),'Type',corrStr);
%     % we use true signal rate as regressor, and response as DV,
%     
%     diff_RT(whichsub,1) = (median(Mul(:,3)) - median(Uni(:,3)))*1000;
    
end

close all;

h=figure('Position',[392   389   600   308]);
hold
for whicht=[0,1]

X=(mean(Uni_subj(:,whicht+1,:),2)); % mean acros tasks
size(X) % all subj
[squeeze(mean(X)), squeeze(std(X)/sqrt(15))] % diff types of measure, meanRT, meanlogRT, mean(medRT)
unisensoryRT=X;


task=RTall(:,1);
AR=RTall(:,2);
AVd=RTall(:,3);

tmp=[];
for i=[0,1,2,3]
     tmp(:,i+1,:)=mean(RTall_subj(:,AVd==i & task==whicht,4:6), 2);
end
size(tmp)
squeeze(mean(tmp))
squeeze(std(tmp)./sqrt(15))
multisensoryRT=tmp;

RTtest=cat(2,X,tmp);
size(RTtest)

MvU=squeeze(RTtest(:,1,:)-RTtest(:,2,:));
size(MvU)
[H,P,CI,STATS]=ttest(MvU,0)
% [p, h, stats] = signrank(MvU,0)


% let's use logRT

subplot(1,2,1)
finalRT=bsxfun(@minus,RTtest(:,2:end,2),RTtest(:,1,2));
x=1:size(finalRT,2);
y=mean(finalRT);

if whicht==0
    plot(x,y,'ko-','markersize',9,'markerfacecolor','w','linewidth',2); hold on
    errorbar_YC(x,y,std(finalRT)/sqrt(15),.3,2,'k');
hold on;
else
    plot(x,y,'ro-','markersize',9,'markerfacecolor','w','linewidth',2); hold on
    errorbar_YC(x,y,std(finalRT)/sqrt(15),.3,2,'r');
hold on;
end
hold on;
plot([-1,5],[0,0],'k--')
ylabel({'Multisensory Log RT (s)' 'relative to Unisensory'})
axis tight
ylim([-0.04, 0.04])
xlim([0.5,4.5])
xlabel('Abs(audiovisual disparity)')
hold on;


%let's plot accraucy
% subj,AVd,[avdlevel,corr], type
x=squeeze(mean(accuracy));
se=squeeze(std(accuracy))/sqrt(15);
for type=2
%     errorbar(x(:,1,type),x(:,2,type),se(:,2,type))
%     hold on;
    Y=squeeze(accuracy(:,:,2,type));
    Ytest=bsxfun(@minus,Y(:,2:end),Y(:,1));
    [H,Pacc,CI,STATS]=ttest(Ytest,0)
end
finalACC=Ytest;

subplot(1,2,2)
x=1:size(finalACC,2);
y=mean(finalACC);

errorbar_YC(x,y,std(finalACC)/sqrt(15),.3,2,'k');
hold on;
plot(x,y,'ko-','markersize',9,'markerfacecolor','w','linewidth',2)
hold on;
plot([-1,5],[0,0],'k--')
ylabel({'Multisensory accuracy' 'relative to Unisensory'})
xlim([0.5,4.5])
% ylim([-0.15,0.05])

absAVd=2/0.55*[0,3];


end
h=get(h,'children');
set(h,'fontsize',15,'linewidth',2,'box','off','tickdir','out',...
    'xtick',[1,4],'xticklabel',round(absAVd,1))


% 
% h = figure('Position',[391   397   491   576]);
% M = mean(g_bias,3);
% Nsubj = size(g_bias,3);
% SE = squeeze(std(permute(g_bias,[3,1,2]))./sqrt(15));
% SE(:,1:2) = M(:,1:2);
% 
% XX=abs(M(:,3:end)) + fliplr(abs(M(:,3:end)));
% M=[M(:,1:2),XX(:,4:end)/2];
% 
% XX=abs(SE(:,3:end)) + fliplr(abs(SE(:,3:end)));
% SE=[SE(:,1:2),XX(:,4:end)/2];
% 
% ccc = [1,0,0; 0,0,0];
% pp=[1,-1];
% for task = [0,1]
%     for AR = 1:2
%         idx = M(:,1)==task&M(:,2)==AR;
%         if AR == 1
%            errorbar([0:3]*pp(task+1)*scale,M(idx,3:end)*scale,SE(idx,3:end)*scale,'o-','linewidth',2.4,...
%                'color',ccc(task+1,:),'markersize',7,'markerfacecolor',ccc(task+1,:))
%         else
%            errorbar([0:3]*pp(task+1)*scale,M(idx,3:end)*scale,SE(idx,3:end)*scale,'s--',...
%                'color',ccc(task+1,:),'markersize',9,'linewidth',1.5);
%         end
%         hold on;
%     end
% end
% plot([-40,40], [0,0],'k--')
% % xlim([-3.5,3.5])
% % ylim([-1.5,1.5])
% xlabel('Distractor rate - Target rate (Hz)')
% ylabel({'Crossmodal bias'; '(reported rate - signal rate) Hz'})
% axis tight
% xlim([-12,12])
% leg = {'A report, high A reliability';
%        'A report, low  A reliability';
%        'V report, high A reliability';
%        'V report, low  A reliability';
%        };
% legend(leg, 'Location', 'northoutside')
% h = get(h, 'children');
% set(h, 'linewidth', 2, 'box', 'off', 'fontsize', 15)
% ylim([-3.5,4.5])



%%



% model prediciton


clc;

rootmain = '/analyse/Project0127/Cao.Y.17/';
behavdir = [rootmain, '/RSA_scripts/'];

load(['/analyse/Project0127/Cao.Y.17/Behav_modelling/',...
    'A_Modelfit2018/Mulconditions_fixpcfixnoise/modelpred/',...
    'ModelPred_Fixpc_FixNoise_final2018.mat']);



outdir=['/analyse/Project0127/Cao.Y.17/Behav_modelling/A_Modelfit2018',...
'/Mulconditions_fixpcfixnoise_corrected_tmp/modelpred/'];
new=load([outdir, 'ModelPred_Fixpc_FixNoise_final2018_MAFu.mat']);



% load(['/analyse/Project0127/Cao.Y.17/Behav_modelling/',...
% 'A_Modelfit2018/Mulconditions_fixpcfixnoise_noRelWeightFusion/',...
% 'modelpred/ModelPred_Fixpc_FixNoise_FusionNoWeighting.mat'])

% load(['/analyse/Project0127/Cao.Y.17/Behav_modelling/',...
% 'A_Modelfit2018/Mulconditions_fixpcfixnoise_fixAV/',...
% 'modelpred/ModelPred_Fixpc_FixNoise_FusionfixAV.mat'])

% load(['/analyse/Project0127/Cao.Y.17/Behav_modelling/',...
%     'A_Modelfit2018/Mulconditions_fixpcfixnoise_perturbPosLogOdds/modelpred/',...
%     'ModelPred_Fixpc_FixNoise_final2018_noise.mat']);

% load(['/analyse/Project0127/Cao.Y.17/Behav_modelling/',...
%     'A_Modelfit2018/Mulconditions_fixpcfixnoise_fullpower/modelpred/',...
%     'ModelPred_Fixpc_FixNoise_final2018_fullpower.mat']);
% % 



clc; 
% load('/analyse/Project0127/Cao.Y.17/Behav_modelling/AttentionModulation/constrained_powerFunc/model3_plugin/model_parm_model1_fixNoise_freepc/modelpred/ModelPred_Freepc_FixNoise.mat');
% % pred_rate_all

SA=[];
SV=[];
SEG=[];
% 
SA=pred_rate_all(:,:,4);
SV=pred_rate_all(:,:,5);
% % 
% % SA=FusionSeg_fromCI(:,:,2);
% % SV=FusionSeg_fromCI(:,:,3);
% 
SEG=[SA(:,1:32), SV(:,33:64)];
% 
all_model=cat(3,pred_rate_all(:,:,1:3),SEG,new.pred_rate_all);
size(all_model)

%%
CIcomp=all_model(:,:,[1,5]);
Fucomp=all_model(:,:,[3,6]);

RSQCI=[];
RSQFu=[];
for subj=1:15
    
tmp=squeeze(CIcomp(subj,:,:));
Preds = {tmp(:,1), tmp(:,2)};
ncond = size(Preds{1},1);
nreg = length(Preds);
f1 = @(x) x*ones(1,ncond);
f2 = @(x) (vectorform(abs(f1(x)-f1(x)')))';
DPred = cellfun(f2,Preds,'un',0);
DPred = tiedrank(cell2mat(DPred));
RSQCI(subj)=(corr(DPred(:,1),DPred(:,2)))^2;

tmp=squeeze(Fucomp(subj,:,:));
Preds = {tmp(:,1), tmp(:,2)};
ncond = size(Preds{1},1);
nreg = length(Preds);
f1 = @(x) x*ones(1,ncond);
f2 = @(x) (vectorform(abs(f1(x)-f1(x)')))';
DPred = cellfun(f2,Preds,'un',0);
DPred = tiedrank(cell2mat(DPred));
RSQFu(subj)=(corr(DPred(:,1),DPred(:,2)))^2;

end
bar(RSQFu)
close all;
figure;
subplot(2,2,1)
A=CIcomp(:,:,1);
B=CIcomp(:,:,2);
plot(A(:),B(:),'.'); hold on;
plot([10,20],[10,20],'r-','linewidth',2)
title('CI old vs. new')
axis tight;
subplot(2,2,2)
bar(RSQCI); title('CI RSQ old new')
axis tight; ylim([0,1])

subplot(2,2,3)
A=Fucomp(:,:,1);
B=Fucomp(:,:,2);
plot(A(:),B(:),'.'); title('Fusion old vs. new');
hold on;
plot([10,26],[10,26],'r-','linewidth',2)
axis tight;
subplot(2,2,4)
bar(RSQFu); title('Fusion RSQ old new')
axis tight; ylim([0,1])


figure;
plot(new.Info.fitted_pc(:,1), Info.fitted_pc(:,1),'ko')



%%



% all_model=pred_rate_all;
% % tmp=all_model';
% for sub=1:15
%     length(unique(tmp(:,sub)))
% end

% all_model=cat(3,pred_rate_all(:,:,1),SEG);
D = Info.Data(1:64, 1:4);
% open D

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
close all;
docorrect=1
h = figure('Position',[391   397   300   300]);
hold;
for whichmodel = [3,6]
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
        uAVd = unique(AVd);
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
    scale = 2/0.55; 
    pp=[1,1];
    for task = [0,1]
        for AR = 1:2
            idx = M(:,1)==task&M(:,2)==AR;
        if AR == 1
            plot((-3:3)*pp(task+1)*scale,M(idx,3:end),'o-','linewidth',2.4,...
               'color',ccc(task+1,:),'markersize',7,'markerfacecolor',ccc(task+1,:))
%            errorbar((-3:3)*pp(task+1)*scale,M(idx,3:end),SE(idx,3:end),'o-','linewidth',2.4,...
%                'color',ccc(task+1,:),'markersize',7,'markerfacecolor',ccc(task+1,:))
        else
%            errorbar((-3:3)*pp(task+1)*scale,M(idx,3:end),SE(idx,3:end),'s--',...
%                'color',ccc(task+1,:),'markersize',9,'linewidth',1.5);
plot((-3:3)*pp(task+1)*scale,M(idx,3:end),'s--',...
               'color',ccc(task+1,:),'markersize',9,'linewidth',1.5);
        end
        hold on;
        end
    end
plot([-40,40], [0,0],'k--')
hold on;

end
% xlim([-3.5,3.5])
% ylim([-1.5,1.5])
% axis tight
xlabel('Vis. - Aud. rate (Hz)')
ylabel({'Crossmodal bias'; '(reported rate - signal rate) Hz'})
axis square
xlim([-12,12])
leg = {'A report, high A reliability';
       'A report, low  A reliability';
       'V report, high A reliability';
       'V report, low  A reliability';
       };
% legend(leg, 'Location', 'northoutside')
h = get(h, 'children');
set(h, 'linewidth', 2, 'box', 'off', 'fontsize', 15)
% ylim(range)

range=[-2.5,2];
% ylim(range)




%%
% check probabilities

clear; clc; close all;
cd('/analyse/Project0127/Cao.Y.17/Behav_modelling/')
load('Data_15subjs_22Trls_MEGextract.mat')

rootmain = '/analyse/Project0127/Cao.Y.17/';
behavdir = [rootmain, '/RSA_scripts/'];
load(['/analyse/Project0127/Cao.Y.17/Behav_modelling/',...
    'A_Modelfit2018/Mulconditions_fixpcfixnoise/modelpred/',...
    'ModelPred_Fixpc_FixNoise_final2018_X.mat']);

XX=[];
data=[];
for whichsub=1:15
    D=Data_sum(:,:,whichsub);
    pOb=D(1:64,5:8)/22;
    
    error=[];
    for whichmodel=[1]
        pMod=squeeze(Info.probabilities(whichsub,:,:,whichmodel));
%         error=[error;sum(sum((pOb-pMod).^2))];
    end
%     XX(:,whichsub)=error;

    data=[data; pMod(:), pOb(:)];
end
plot(data(:,2),data(:,1),'ko')
axis square
%%
% plot(XX([1,3],:),'ko-')







%% Create figures 6-10
% old file was figures2024rev1
load ExportShifts.mat
Ntable=ExportShifts;
n=sum(Ntable{:,:},"all");
I=size(Ntable,1);
%% Correspondence analysis based on all the obs
CorAna(Ntable)
title('')
if prin==1
    % print to postscript
    print -depsc figs2024\SIca.eps;
end

%%  under H0
findEmp=struct;
% Generate nsimul contingency tables
findEmp.nsimul=1000; % 0;
% Under the null hypothesis of independence
findEmp.underH0=true;
% Store the nsimul robust distance sorted (for each row)
findEmp.StoreSim=true;

% Detect outlying rows using a simultaneous confidence level
conflev=1-0.01/I;
RAWH0=mcdCorAna(Ntable,'plots',0,'findEmpiricalEnvelope',findEmp, ...
    'conflev',conflev,'bdp',0.5);


%%  under H1
% Example of findEmpiricalEnvelope a struct
findEmp=struct;
% Generate nsimul contingency tables
findEmp.nsimul=1000; % 0;
% Under the null hypothesis of independence
findEmp.underH0=false;
% Store the nsimul robust distance sorted (for each row)
findEmp.StoreSim=true;

% Detect outlying rows using a simultaneous confidence level
conflev=1-0.01/I;
RAWH1=mcdCorAna(Ntable,'plots',0,'findEmpiricalEnvelope',findEmp, ...
    'conflev',conflev,'bdp',0.5);


%%  Create left panel of Figure 7
plots=struct;
plots.minMD = true;
plots.inertia = false;
plots.addRowNames = true;

outFSH0=FSCorAna(RAWH0,'plots',plots,'init',round(n*0.25),'conflev',[0.05 0.95],'resc',false);
title('')
ylabel('Minimum MD')
h = findobj('Type','figure','Tag','pl_in');   % find figure handles with that tag
if ~isempty(h)
    close(h)                                  
end

if prin==1
    % print to postscript
    print -depsc figs2024\SiFSH0.eps;
end

%% Create right panel of Figure 7 (Under H0 rescaled)
plots=struct;
plots.minMD = true;
plots.inertia = false;
plots.addRowNames = true;

outFSH0=FSCorAna(RAWH0,'plots',plots,'init',round(n*0.25),'conflev',[0.05 0.95],'resc',true);
title('')
ylabel('Minimum MD')
if prin==1
    % print to postscript
    print -depsc figs2024\SiFSH0resc.eps;
end

%% Create left panel of Figure 8: Under H1 not rescaled 
FSCorAna(RAWH1,'plots',plots,'init',round(n*0.25),'conflev',[0.05 0.95],'resc',false);
title('')
ylabel('Minimum MD')
% ylabel('')
if prin==1
    % print to postscript
    print -depsc figs2024\SiFSH1.eps;
end

%% Create right  panel of Figure 8: Under H1 rescaled 
plots.addBonfLine=true;
outFSH1=FSCorAna(RAWH1,'plots',plots,'init',round(n*0.25),'conflev',0.90,'resc',true);
title('')
% ylabel('Minimum MD')
% ylim([40 400])
if prin==1
    % print to postscript
    print -depsc figs2024\SiFSH1resc1.eps;
end

%% Figure 9: CA with two supplementary units
selsup={'TR', 'AE'};
selactive=setdiff(Ntable.Properties.RowNames,selsup);
Nsupr=Ntable(selsup,:);
Nactive=Ntable(selactive,:);
Sup=struct;
Sup.r=Nsupr;
outCorAna=CorAna(Nactive,'Sup',Sup);
title('')

if prin==1
    % print to postscript
    print -depsc figs2024\SiCAwithSU.eps;
end


%% ExportShifts
% clear
% conflevSim=0.99;
% load ExportShifts.mat
% Ntable=ExportShifts;
% N=Ntable;
% %% Create Figure 10
% conflev=0.95;
% % Example of findEmpiricalEnvelope a struct
% findEmp=struct;
% % Generate nsimul contingency tables
% findEmp.nsimul=1000;
% % Store the nsimul robust distance sorted (for each row)
% findEmp.StoreSim=true;
% % Under the null hypothesis of independence
% findEmp.underH0=true;
% % call to mcdCorAna
% RAWH0=mcdCorAna(N,'plots',0,'findEmpiricalEnvelope',findEmp,'bdp',0.5,'msg',0,'conflev',conflev);
% outFSH0=FSCorAna(RAWH0,'plots',0,'msg',0,'StoreSim',false,'conflev',conflevSim);
% 
% 
% % Under the alternative hypothesis of dependence
% findEmp.underH0=false;
% % call to mcdCorAna
% RAWH1=mcdCorAna(N,'plots',0,'findEmpiricalEnvelope',findEmp,'bdp',0.5,'msg',0,'conflev',conflev);
% outFSH1=FSCorAna(RAWH1,'plots',0,'msg',0,'StoreSim',false,'conflev',conflevSim);

%% Figure 10: Monitoring plots of inertia

lwdreal=1.5;

figure
subplot(2,2,1)
plot(outFSH0.ine(:,1),outFSH0.ine(:,2),'LineWidth',lwdreal)


hold on
% include specified tag in the current plot
gine=outFSH0.ine(:,[1 3:end]);

lwdenv=2;
% Superimpose 50% envelope
line(gine(:,1),gine(:,3),'LineWidth',lwdenv,'LineStyle','-.','Color','k','tag','env');
% Superimpose alpha1% and alpha2% envelope
line(gine(:,1),gine(:,2),'LineWidth',lwdenv,'LineStyle','-.','Color',[0.2 0.8 0.4],'tag','env');
line(gine(:,1),gine(:,4),'LineWidth',lwdenv,'LineStyle','-.','Color',[0.2 0.8 0.4],'tag','env');
ylabel('Inertia')
xlabel('Subset size m')

subplot(2,2,2)
plot(outFSH1.ine(:,1),outFSH1.ine(:,2),'LineWidth',lwdreal)

hold on
% include specified tag in the current plot
gine=outFSH1.ine(:,[1 3:end]);

lwdenv=2;
% Superimpose 50% envelope
line(gine(:,1),gine(:,3),'LineWidth',lwdenv,'LineStyle','-.','Color','k','tag','env');
% Superimpose alpha1% and alpha2% envelope
line(gine(:,1),gine(:,2),'LineWidth',lwdenv,'LineStyle','-.','Color',[0.2 0.8 0.4],'tag','env');
line(gine(:,1),gine(:,4),'LineWidth',lwdenv,'LineStyle','-.','Color',[0.2 0.8 0.4],'tag','env');
ylabel('Inertia')
xlabel('Subset size m')
prin=0;
if prin==1
    print -depsc figs\I6.eps;
end


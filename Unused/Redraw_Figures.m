close all
%% Read in SRII basis spectra from figures from Daniel 
clear all
fig1E = openfig('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\Figures\Fig1\SrII\fig_2_a_v2.fig');
fig1E = gcf;  axObjs = fig1E.Children; dataObjs = axObjs.Children;
SRII_Basis(2,:) = dataObjs(2).YData;        % 200 us
SRII_Basis(1,:) = dataObjs(1).YData;              % 20 us 
dataObjs = findobj(fig1E,'-property','XData');
qSRII = dataObjs(1).XData;

% Plot Figure 1E
% close all
figure('Position', [850 250 500 380])
plot(qSRII,SRII_Basis(1,:)+0.0005,'o','LineWidth',0.5,'MarkerSize',4,'color','#0072BD')
hold on
plot(qSRII,SRII_Basis(2,:),'o','LineWidth',0.5,'MarkerSize',4,'color','#4DBEEE')
axis([0.0 1.1 -6*10^-3 2*10^-3]);
xlabel('q (Å^{-1})')
ylabel('\DeltaS(q)')
box on; set(gca,'linewidth',3); set(gca,'FontSize',18);
xticks([0:0.2:1]); 
grid on
ax = gca; 
ax.GridLineWidth = 1;

%% Read in SRII basis spectra from figures from Daniel 
clear all
fig1F = openfig('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\Figures\Fig1\SrIIHtrII\fig_2_c_v2.fig');
fig1F = gcf;  axObjs = fig1F.Children; dataObjs = axObjs.Children;
SRIIHtrII_Basis(2,:) = dataObjs(2).YData;        % 200 us
SRIIHtrII_Basis(1,:) = dataObjs(1).YData;              % 20 us 
dataObjs = findobj(fig1F,'-property','XData');
qSRII = dataObjs(1).XData;

% Plot Figure 1F
figure('Position', [1050 250 500 380])
plot(qSRII,SRIIHtrII_Basis(1,:)-0.0005,'o','LineWidth',0.5,'MarkerSize',4,'color','#EDB120')
hold on
plot(qSRII,SRIIHtrII_Basis(2,:),'o','LineWidth',0.5,'MarkerSize',4,'color','#7E2F8E')
axis([0.0 1.1 -6*10^-3 2*10^-3]);
xlabel('q (Å^{-1})')
ylabel('\DeltaS(q)')
box on; set(gca,'linewidth',3); set(gca,'FontSize',18);
xticks([0:0.2:1]); 
grid on
ax = gca; 
ax.GridLineWidth = 1;

%% Read in time-resolved decomposition from Daniel's figures. 
clear all
fig1C = openfig('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\Figures\Fig1\SrII\fig_2_b_n_2_5.fig');
fig1C = gcf;  axObjs = fig1C.Children; dataObjs = axObjs.Children;
Component1 = dataObjs(1).YData;        
Component2 = dataObjs(2).YData;          
Component3 = dataObjs(3).YData; 
Component4 = dataObjs(4).YData; 

dataObjs = findobj(fig1C,'-property','XData');
Time = dataObjs(1).XData;
Time2 = dataObjs(3).XData;

% Plot Figure 1C
figure('Position', [1250 250 500 380])
plot(Time,Component1,'o','LineWidth',3.5,'MarkerSize',8,'color','#0072BD')
hold on
plot(Time,Component2,'o','LineWidth',3.5,'MarkerSize',8,'color','#4DBEEE')
plot(Time2,Component4,'LineWidth',3.5,'color','#4DBEEE')
plot(Time2,Component3,'LineWidth',3.5,'color','#0072BD')
axis([1 20000 0 1.1])
xticks([1 10 100 1000 10000])
xticklabels({'1','10','100','1000','10000'})
xlabel('Time (\mus)')
ylabel('Population')
box on; set(gca,'linewidth',3); set(gca,'FontSize',18);
set(gca,'xscale','log')

%% Read in time-resolved decomposition from Daniel's figures. 
clear all
fig1D = openfig('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\Figures\Fig1\SrIIHtrII\fig_2_d_n_2_5.fig');
fig1D = gcf;  axObjs = fig1D.Children; dataObjs = axObjs.Children;
Component1 = dataObjs(1).YData;        
Component2 = dataObjs(2).YData;          
Component3 = dataObjs(3).YData; 
Component4 = dataObjs(4).YData; 

dataObjs = findobj(fig1D,'-property','XData');
Time = dataObjs(1).XData;
Time2 = dataObjs(3).XData;

% Plot Figure 1D
figure('Position', [550 250 500 380])
plot(Time,Component1,'o','LineWidth',3.5,'MarkerSize',8,'color','#EDB120')
hold on
plot(Time,Component2,'o','LineWidth',3.5,'MarkerSize',8,'color','#7E2F8E')
plot(Time2,Component3,'LineWidth',3.5,'color','#EDB120')
plot(Time2,Component4,'LineWidth',3.5,'color','#7E2F8E')
axis([1 20000 0 1.1])
xticks([1 10 100 1000 10000])
xticklabels({'1','10','100','1000','10000'})
xlabel('Time (\mus)')
ylabel('Population')
box on; set(gca,'linewidth',3); set(gca,'FontSize',18);
set(gca,'xscale','log')


%% Read in time-resolved data SRII
clear all
fig1A = openfig('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\Figures\Fig1\SrII\fig_3a.fig');
fig1A = gcf;  axObjs = fig1A.Children; dataObjs = axObjs.Children;
dataObjs = findobj(fig1A,'-property','YData');
Data1 = dataObjs(1).YData;        
Data2 = dataObjs(2).YData;          
Data3 = dataObjs(4).YData; 
Data4 = dataObjs(5).YData; 
Data5 = dataObjs(7).YData; 
Data6 = dataObjs(8).YData; 
Data7 = dataObjs(10).YData;        
Data8 = dataObjs(11).YData;          
Data9 = dataObjs(13).YData; 
Data10 = dataObjs(14).YData; 
Data11 = dataObjs(16).YData; 
Data12 = dataObjs(17).YData; 
Data13 = dataObjs(19).YData; 
Data14 = dataObjs(20).YData;        
Data15 = dataObjs(22).YData;          
Data16 = dataObjs(23).YData; 
Data17 = dataObjs(25).YData; 
Data18 = dataObjs(26).YData; 

dataObjs = findobj(fig1A,'-property','XData');
qSRII = dataObjs(1).XData;

%close all
% Plot Figure 1B
figure('Position', [550 150 500 680])
plot(qSRII,Data2-0.0032,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
hold on
plot(qSRII,Data1-0.0032,'LineWidth',1.5,'color','red')
plot(qSRII,Data4-0.0028,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data3-0.0028,'LineWidth',1.5,'color','red')
plot(qSRII,Data6-0.0024,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data5-0.0024,'LineWidth',1.5,'color','red')
plot(qSRII,Data8-0.002,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data7-0.002,'LineWidth',1.5,'color','red')
plot(qSRII,Data10-0.0016,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data9-0.0016,'LineWidth',1.5,'color','red')
plot(qSRII,Data12-0.0012,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data11-0.0012,'LineWidth',1.5,'color','red')
plot(qSRII,Data14-0.0008,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data13-0.0008,'LineWidth',1.5,'color','red')
plot(qSRII,Data16-0.0004,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data15-0.0004,'LineWidth',1.5,'color','red')
plot(qSRII,Data18,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data17,'LineWidth',1.5,'color','red')

axis([0.0 1.1 -18*10^-3 2*10^-3]);
xlabel('q (Å^{-1})')
ylabel('\DeltaS(q)')
box on; set(gca,'linewidth',3); set(gca,'FontSize',18);
xticks([0:0.2:1]); 
grid on
ax = gca; 
ax.GridLineWidth = 1;

%% Read in time-resolved data SRII:HtrII
clear all
fig1B = openfig('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\Figures\Fig1\SrIIHtrII\fig_1_b_v2.fig');
fig1B = gcf;  axObjs = fig1B.Children; dataObjs = axObjs.Children;
dataObjs = findobj(fig1B,'-property','YData');
Data1 = dataObjs(1).YData;        
Data2 = dataObjs(2).YData;          
Data3 = dataObjs(4).YData; 
Data4 = dataObjs(5).YData; 
Data5 = dataObjs(7).YData; 
Data6 = dataObjs(8).YData; 
Data7 = dataObjs(10).YData;        
Data8 = dataObjs(11).YData;          
Data9 = dataObjs(13).YData; 
Data10 = dataObjs(14).YData; 
Data11 = dataObjs(16).YData; 
Data12 = dataObjs(17).YData; 
Data13 = dataObjs(19).YData; 
Data14 = dataObjs(20).YData;        
Data15 = dataObjs(22).YData;          
Data16 = dataObjs(23).YData; 
Data17 = dataObjs(25).YData; 
Data18 = dataObjs(26).YData; 

dataObjs = findobj(fig1B,'-property','XData');
qSRII = dataObjs(1).XData;

% Plot Figure 1B
figure('Position', [550 150 500 680])
plot(qSRII,Data2-0.0009,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
hold on
plot(qSRII,Data1-0.0009,'LineWidth',1.5,'color',[0.8350 0.1780 0.2840])
plot(qSRII,Data4-0.0008,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data3-0.0008,'LineWidth',1.5,'color','red')
plot(qSRII,Data6-0.0006,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data5-0.0006,'LineWidth',1.5,'color','red')
plot(qSRII,Data8-0.0005,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data7-0.0005,'LineWidth',1.5,'color','red')
plot(qSRII,Data10-0.0004,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data9-0.0004,'LineWidth',1.5,'color','red')
plot(qSRII,Data12-0.0003,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data11-0.0003,'LineWidth',1.5,'color','red')
plot(qSRII,Data14-0.0002,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data13-0.0002,'LineWidth',1.5,'color','red')
plot(qSRII,Data16-0.0001,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data15-0.0001,'LineWidth',1.5,'color','red')
plot(qSRII,Data18,'o','LineWidth',0.5,'MarkerSize',4,'color','black')
plot(qSRII,Data17,'LineWidth',1.5,'color','red')

axis([0.0 1.1 -18*10^-3 2*10^-3]);
xlabel('q (Å^{-1})')
ylabel('\DeltaS(q)')
box on; set(gca,'linewidth',3); set(gca,'FontSize',18);
xticks([0:0.2:1]); 
grid on
ax = gca; 
ax.GridLineWidth = 1;

%% Read in the comparision plots from Daniel's figures. 
clear all
fig3B = openfig('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\Figures\Fig3\SRII_BR_overlay.fig');
fig3B = gcf;  axObjs = fig3B.Children; dataObjs = axObjs.Children;
SRII_Basis(1,:) = dataObjs(2).YData;        % 200 us
bR_Basis(1,:) = dataObjs(1).YData;              % 20 us 
dataObjs = findobj(fig3B,'-property','XData');
qSRII = dataObjs(2).XData;
qbR = dataObjs(1).XData;

% Plot Figure 3B
figure('Position', [1050 250 500 380])
plot(qSRII,SRII_Basis(1,:),'o','LineWidth',0.5,'MarkerSize',4,'color','#0072BD')
hold on
plot(qbR,bR_Basis(1,:),'o','LineWidth',0.5,'MarkerSize',4,'color','k')
axis([0.0 1.1 -4*10^-3 2*10^-3]);
xlabel('q (Å^{-1})')
ylabel('\DeltaS(q)')
box on; set(gca,'linewidth',3); set(gca,'FontSize',18);
xticks([0:0.2:1]); 
grid on
ax = gca; 
ax.GridLineWidth = 1;

%% Read in the comparision plots from Daniel's figures. 
clear all
fig3A = openfig('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\Figures\Fig3\SRII_SRII_HtrII_overlay.fig');
fig3A = gcf;  axObjs = fig3A.Children; dataObjs = axObjs.Children;
SRII_Basis(1,:) = dataObjs(4).YData;        % 200 us
SRII_Basis(2,:) = dataObjs(3).YData;              % 20 us 
SRIIHtrII_Basis(2,:) = dataObjs(1).YData;        % 200 us
SRIIHtrII_Basis(1,:) = dataObjs(2).YData;        % 200 us

dataObjs = findobj(fig3A,'-property','XData');
qSRII = dataObjs(2).XData;

% Plot Figure 1F
figure('Position', [1050 250 500 380])
plot(qSRII,SRII_Basis(1,:),'o','LineWidth',0.5,'MarkerSize',4,'color','#4DBEEE')
hold on
plot(qSRII,SRII_Basis(2,:),'o','LineWidth',0.5,'MarkerSize',4,'color','#0072BD')
plot(qSRII,SRIIHtrII_Basis(1,:),'o','LineWidth',0.5,'MarkerSize',4,'color','#7E2F8E')
plot(qSRII,SRIIHtrII_Basis(2,:),'o','LineWidth',0.5,'MarkerSize',4,'color','#EDB120')

axis([0.0 1.1 -5*10^-3 2*10^-3]);
xlabel('q (Å^{-1})')
ylabel('\DeltaS(q)')
box on; set(gca,'linewidth',3); set(gca,'FontSize',18);
xticks([0:0.2:1]); 
grid on
ax = gca; 
ax.GridLineWidth = 1;



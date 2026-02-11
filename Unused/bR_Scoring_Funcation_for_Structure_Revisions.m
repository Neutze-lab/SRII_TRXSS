close all;
ReadData = 0;           % Slow part of the script is reading the theoretical data. Set to 1 to read in data. 
if ReadData ==  1; 
    clear all;  ReadData = 1; 
end
%% ==================================================================================
% Order the simulations as they have been done. 
stepnr = [201 126 102 127 103 128 104 129 117 130  ...
                105 131 106 132 107 133 108 134 118 135 ...
                109 136 110 137 111 138 112 139 119 140 ...
                113 141 114 142 115 143 116 144 120 145 ...
                121 146 122 147 123 148 124 149 125 150];             
                % This lists the simulations in the correct order for plotting. 
svdStart = 1501; svdEnd = 2001;                                            %                            
FitStart = 5; FitEnd = 410;                                                 % Low q and high q limits used to fit experimental data. 
Offset =  0.005*10^-3;
%% ===================================================================================
% Read in results from Daniel's Gromacs analysis 
if ReadData == 1; 
for kk = 1:size(stepnr,2)
        filename1 = ["C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\Daniel\MatLabScripts\grid-search-RICHARD\bR-fit-grid1-grid2\fit\bR-190BOG-101\input_fit\bR-190BOG-101-rest-mem-rest-prot\bR-190BOG-101-rest-mem-rest-prot"] ; 
        filename2 = ['C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\Daniel\MatLabScripts\grid-search-RICHARD\bR-fit-grid1-grid2\fit\bR-190BOG-' num2str(stepnr(kk)) '\input_fit\bR-190BOG-' num2str(stepnr(kk)) '-exc-mem-exc-prot\bR-190BOG-' num2str(stepnr(kk)) '-exc-mem-exc-prot'] ;
        filename3 = ['C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\Daniel\MatLabScripts\grid-search-RICHARD\bR-fit-grid1-grid2\fit\bR-190BOG-' num2str(stepnr(kk)) '\input_fit\bR-190BOG-' num2str(stepnr(kk)) '-rest-mem-exc-prot\bR-190BOG-' num2str(stepnr(kk)) '-rest-mem-exc-prot'] ;   
        filename4 = ['C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\Daniel\MatLabScripts\grid-search-RICHARD\bR-fit-grid1-grid2\fit\bR-190BOG-' num2str(stepnr(kk)) '\input_fit\bR-190BOG-' num2str(stepnr(kk)) '-exc-mem-rest-prot\bR-190BOG-' num2str(stepnr(kk)) '-exc-mem-rest-prot'] ;  
        filename = [filename1 filename2 filename3 filename4]; 
% Protein moving, membrane moving;   
        Load = [filename(1) filename(2)];   
        DataDI = LoadDiffData(Load);    % Function LoadDiffData defined at bottom.
%  Protein vacuum only term. 
        [U,S,V] = svd(squeeze(DataDI(6,:,svdStart:svdEnd)));
        dProtVac(kk,:) = mean(V(:,1))*U(:,1)*S(1,1);
% Load  data with protein changes and resting membrane held fixed.
        Load = [filename(1) filename(3)];
        DataDI = LoadDiffData(Load);   
        [U,S,V] = svd(squeeze(DataDI(2,:,svdStart:svdEnd)));
        dPmMfgvac(kk,:) = mean(V(:,1))*U(:,1)*S(1,1);
% Next load complementary data data with protein changes and excited membrane held fixed.
        Load = [filename(4) filename(2)];
        DataDI = LoadDiffData(Load); 
        [U,S,V] = svd(squeeze(DataDI(2,:,svdStart:svdEnd)));
        dPmMfevac(kk,:) = mean(V(:,1))*U(:,1)*S(1,1);
% Now look at the micelle alone. 
        Load = [filename(2) filename(1)];
        DataDI = LoadDiffData(Load); 
        [U,S,V] = svd(squeeze(DataDI(4,:,svdStart:svdEnd)));
        dPMicelle(kk,:) = mean(V(:,1))*U(:,1)*S(1,1);
end
end
        [U6,S6,V6]  = svd(transpose(dPMicelle));  
        dPMicelleSVD = mean(V6(:,1))*U6(:,1)*S6(1,1);
%% ===================================================================================
% Load additional data outside of the loop. 
% Load in Solvent Exclusion term from Solvent_Correction_Gromacs.m
        TEMP10 = load('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\Daniel\MatLabScripts\grid-search-RICHARD\SolventExclusionFactor.mat'); 
        SolventExclusionFactor = TEMP10.SolventExclusionFactor2;   %  TempExpt = transpose(TempExpt); 
% Read in experimental data and prepare vectors for fitting. 
      %  TEMP1 = load('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\Daniel\MatLabScripts\BR_Basis_Spectra\bRbasisRN.mat'); 
        TEMP1 = load('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\Daniel\MatLabScripts\grid-search-RICHARD\bRbasisRN.mat'); 
        TempExpt = TEMP1.Data2save; TempExpt = transpose(TempExpt); 
% Read in membrane correction from the BOG titration.         
        TEMP101 = load('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\Daniel\MatLabScripts\BOG-titration-RICHARD\MembraneComponent2.mat'); 
        SVDcomp2 = TEMP101.Data2save; 
% Read in membrane correction from the BOG titration.         
        TEMP102 = load('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\Daniel\MatLabScripts\grid-search-RICHARD\MicelleFluctuationsFinalQuarter.mat'); 
        SVDmicelle = TEMP102.Data2save;     
        SVDcomp2 = SVDmicelle*max(SVDcomp2(2,:))/max(SVDmicelle(2,:))   ; % Use this line if you want to fit using the micelle correction from last 500 pdb files. 
             
    dataSize = 475;            % DataSize of 480 corresponds to 1.0353 Å-1
                                         % If this is not constrained then get "NaN" when using interp1
    q20 = TempExpt(:,1); q20 = q20(1:dataSize);    % q20 has experimental data cut so interperlation does not fail. 
    ExptData20 = (TempExpt(:,2)+Offset)*10^8;             % Multiply by 10^8 so optimization routine works.  
    ExptData30 = (TempExpt(:,3)+Offset)*10^8;             % Multiply by 10^8 so optimization routine works. 
     
%% ===============================================================================================================
% The following loop does the fitting of the data & searchers for optimal parameters
    BestYet = 10^13;                        % Need a start value while fitting to minimizing below. 
% Following loop allows various parameters to be optimized stepwise. 
for ll = 1; %:21
        Bfactor = 100;       
        % Bfactor = 2.6 + ((ll-11)/40); 
        DampWidth = 0.18;
        % DampWidth = 0.18 + (ll-11)/100;     % Allows micell damping 
        qStretch = 0.99;   
        % qStretch = 1.0 + (ll-11)/200;   % Stretch the theoretical q to allow for errors in energy, distance etc.....
        WeightMicelle = 0.5;               % 0.5 gives a balance between dark & excite micell structures.  
        % Incorporate BOG correction into the analysis. 
        BOGcorrection =  115;  
        % BOGcorrection = 103 + (ll-11); %*4;   % Correction for deviations from 190 BOG simulation. 
for kk = 1:size(stepnr,2)
%% ==========================================================================
 % Prepare the matrices we use to score the fit.
     q = [0:0.01:2]*qStretch;          % q domain chosen to match Crysol calculations. 
     LowDamp = normpdf(q,0,DampWidth); LowDamp = LowDamp/max(LowDamp);  
     AllDamp = normpdf(q,0,Bfactor); AllDamp = AllDamp/max(AllDamp);  
     dProtMemVac(kk,:) = ((1-WeightMicelle)*dPmMfevac(kk,:) + WeightMicelle*dPmMfgvac(kk,:));               % This is what is used to fit our data. 
    dMemVac(kk,:) = (dProtMemVac(kk,:) + BOGcorrection*SVDcomp2(2,:) - dProtVac(kk,:)).*LowDamp;                        
    dProtMemVac(kk,:) =  AllDamp.*(dProtVac(kk,:) + dMemVac(kk,:)) ;      % Damp the micell correction.
%  Interpolate theoretical curves onto the data to allow fitting. 
    dProtMemVacInterp(:,kk) = interp1(q,dProtMemVac(kk,:),q20); 
    dProtVacInterp(:,kk) = interp1(q,dProtVac(kk,:),q20);     
% Shorten the experimental vectors to the chosen fitting lengths. 
    q2 = q20(FitStart:FitEnd); 
    ExptData2 = ExptData20(FitStart:FitEnd); 
    ExptData3 = ExptData30(FitStart:FitEnd); 
% Shorten the theoretical vectors accordingly    
    dProtMemVacInterp2(:,kk) = dProtMemVacInterp(FitStart:FitEnd,kk);
    dProtVacInterp2(:,kk) = dProtVacInterp(FitStart:FitEnd,kk);    
% Minimize against experimental data without low-q damping of membrane.  
    weight = q2.^2.*gaussmf(q2,[0.4 -0.3]);  weight = weight/max(weight);
    Contrast = interp1(q,SolventExclusionFactor,q2);
end 

%% Next step is to make the fit 
for jj = 1:size(stepnr,2)
for kk = 1:size(stepnr,2)
%% Version which allows scaling to be optimized for both pairs together. 
    fun = @(x)(sum(((x(1)*(dProtMemVacInterp2(:,jj)).*Contrast - ExptData2).*weight).^2 + ...
                                    ((x(1)*(dProtMemVacInterp2(:,kk)).*Contrast - ExptData3).*weight).^2));               
    x0 = [2];  x1 = fminsearch(fun,x0);
    scoreAll(jj,kk) = sum(((x1(1)*(dProtMemVacInterp2(:,jj)).*Contrast - ExptData2).*weight).^2 + ...
                                        ((x1(1)*(dProtMemVacInterp2(:,kk)).*Contrast - ExptData3).*weight).^2);
                                    % Note this score function is emphasizing the second component
    if  scoreAll(jj,kk)  < BestYet
            BestYet = scoreAll(jj,kk);        
            % [jj kk BestYet/10^10];
            Best2 = x1(1)*(dProtMemVacInterp2(:,jj)).*Contrast;
            Best3 = x1(1)*(dProtMemVacInterp2(:,kk)).*Contrast;
            BestProtOnly = x1(1)*(dProtVacInterp2(:,kk)).*Contrast;
            BestStretch = qStretch; 
            BestDampWidth = DampWidth; 
            BestWeightMicelle = WeightMicelle; 
            BestBfactor = Bfactor; 
            BestScale = x1(1); 
            BestBOGcorrection = BOGcorrection; 
            [jj kk]
    end           
%% ========================================================================
end 
end
end
%% Print the optimized values to screen.  
BestStretch
BestDampWidth
BestBfactor
BestBOGcorrection
BestYet/1e11
   
%% Plot the results of this analysis. 

% figure('Position', [1050 250 500 400])
figure('Position', [1050 250 500 380])
% plot(q2,q2.*ExptData2,'o','LineWidth',0.5,'MarkerSize',4,'color','k');
 plot(q2,q2.*ExptData3,'o','LineWidth',0.5,'MarkerSize',4,'color','k');
hold on; 
axis([0.0 1 -10*10^4 7*10^4]);
% plot(q2,q2.*Best2,'LineWidth',4);
plot(q2,q2.*Best3,'LineWidth',5,'color',[0.8500 0.3250 0.0980]);
plot(q2,q2.*BestProtOnly,'LineWidth',3,'Color',[0.9290 0.6940 0.1250]);
% title('Best Fit Component 1','FontSize',14)
xlabel('q (Å^{-1})')
ylabel('q \cdot \DeltaS(q)')
box on; set(gca,'linewidth',3); set(gca,'FontSize',18);
xticks([0:0.2:1]); 
grid on
ax = gca; 
ax.GridLineWidth = 1;

% plot(qSRII,SRII_Basis(1,:),'o','LineWidth',0.5,'MarkerSize',4,'color','#0072BD')
% hold on
% plot(qbR,bR_Basis(1,:),'o','LineWidth',0.5,'MarkerSize',4,'color','k')
% axis([0.0 1.1 -4*10^-3 2*10^-3]);
% xlabel('q (Å^{-1})')
% ylabel('\DeltaS(q)'


%% ========================================================================
% Functions used to load difference data from simulations. 
function DataTemp = LoadDiffData(yi)
run_name = yi; 
excited.protmem = load(char(run_name(2)+'_intensities_protmem.mat'));
excited.membrane = load(char(run_name(2)+'_intensities_membrane.mat'));
excited.protein = load(char(run_name(2)+'_intensities_protein.mat'));
ground.protmem = load(char(run_name(1)+'_intensities_protmem.mat'));
ground.membrane = load(char(run_name(1)+'_intensities_membrane.mat'));
ground.protein = load(char(run_name(1)+'_intensities_protein.mat'));
DI(1,:,:) = excited.protmem.intSol - ground.protmem.intSol;       
DI(2,:,:)  = excited.protmem.intVac - ground.protmem.intVac;       
DI(3,:,:) = excited.membrane.intSol - ground.membrane.intSol;   
DI(4,:,:)  = excited.membrane.intVac - ground.membrane.intVac;    
DI(5,:,:)  = excited.protein.intSol - ground.protein.intSol;      
DI(6,:,:)  = excited.protein.intVac - ground.protein.intVac;
DataTemp = DI; 
end
% ===================================
function scoreTwoDfit = TwoDfit(yi,yii,yiii,yiv)
 for i = 1: size(yi,1)
     for j = 1:size(yi,2)
%        Temp1(i,j) = yii(5)*(((i-1-abs(yii(1)))/3)^2/(yii(2)/3)^2 + ...
%                                        ((j-1-abs(yii(3)))/6)^2/(yii(4)/6)^2) + yii(6);             
       Temp1(i,j) = yii(3)*(((i-yiv(2)))^2/yii(1)^2 + ((j-yiv(1)))^2/yii(2)^2) + yii(4);    
     end 
 end 
 scoreTwoDfit = sum(sum((yi - Temp1.^yiii).^2));
end 

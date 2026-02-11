close all;
ReadData = 1;           % Slow part of the script is reading the theoretical data. Set to 1 to read in data. 
if ReadData == 1; 
clear all;  ReadData = 1; 
end
%% ==================================================================================
% Order the simulations so that the grid-plot makes sense. 
stepnr = [731 731 702 732 703 733 704 734 705 735 706 736 ...    
                707 737 708 738 709 739 710 740 711 741 712 742 ...
                713 743 714 744 715 745 716 746 717 747 718 748 ...
                 719 749 720 750 721 751 722 752 723 753 724 754 ...
                 725 755 726 756 727 757 728 758 729 759 730 760];  % Swap 728 for 757 for a bit. 
             
% This lists the simulations in the correct order for plotting.                                                                      
corrStart = 90; corrEnd = 425;       % Correlations on protein-only-vaccum curve with data q1(90) = 0.27 Å-1 < q < q1(425) 0 1.00 Å-1. 
FitStart = 30; FitEnd = 425;            % Low q and high q limits used to fit experimental data. 
ScaleCorr1 = 0.97;                 % Scales the q-axis so is in effect optimizing the detector to sample distance.    
q = [0:0.01:2]*17.7/18;         % Correct X-ray energy for the weighted average of undulator & not the undulator peak. 
KeepData = 99;                   % Number of curves of a possible 496 trajectories kept in this analysis. 
BestYet = 2500;                     % A value put in so that a scoring function would be lower than this. 
Offset =  -0.05*10^-3;       
BoltzmanWeight = 7;   % This tunes the choice of structures to be averaged away from the best structure. 

%% ===================================================================================
% Read in experimental data from Richard's analysis.       
        SRIIHtrIIold = load("Data\SRII_HtrII_BS_1.mat");
        Expt1(1,:) = transpose(SRIIHtrIIold.SRII_HtrII_BS_1) + Offset; 
        q1 = ScaleCorr1*transpose(SRIIHtrIIold.qIN); 
% Read in Solvent Exclusion term from Solvent_Correction_Gromacs.m
        TEMP10 = load('Data\SolventExclusionFactorSRIIHtrII.mat'); 
        SolventExclusionFactor = TEMP10.SolventExclusionFactor2;   %  TempExpt = transpose(TempExpt);  
% Read in Undulator spectrum     
       TEMP103 = load('Data\U17.mat');    
       U17 = TEMP103.Data2save;   
       
%% ===================================================================================
% Read in results from Daniel's Gromacs analysis 
if ReadData == 1; 
for kk = 1:size(stepnr,2)    
        filename1 = ["Data\SrIIHtrII-grid7\input_fit\SrIIHtrII-220BOG-701\input_fit\SrIIHtrII-220BOG-701-rest-mem-rest-prot\SrIIHtrII-220BOG-701-rest-mem-rest-prot"] ;
        filename2 = ['Data\SrIIHtrII-grid7\input_fit\SrIIHtrII-220BOG-' num2str(stepnr(kk)) '\input_fit\SrIIHtrII-220BOG-' num2str(stepnr(kk)) '-exc-mem-exc-prot\SrIIHtrII-220BOG-' num2str(stepnr(kk)) '-exc-mem-exc-prot'] ;
        filename3 = ['Data\SrIIHtrII-grid7\input_fit\SrIIHtrII-220BOG-' num2str(stepnr(kk)) '\input_fit\SrIIHtrII-220BOG-' num2str(stepnr(kk)) '-exc-mem-rest-prot\SrIIHtrII-220BOG-' num2str(stepnr(kk)) '-exc-mem-rest-prot'] ;
        filename4 = ['Data\SrIIHtrII-grid7\input_fit\SrIIHtrII-220BOG-' num2str(stepnr(kk)) '\input_fit\SrIIHtrII-220BOG-' num2str(stepnr(kk)) '-rest-mem-exc-prot\SrIIHtrII-220BOG-' num2str(stepnr(kk)) '-rest-mem-exc-prot'] ;            
        filename = [filename1 filename2 filename3 filename4];     
% Protein moving, membrane moving;   
        Load = [filename(1) filename(2)];
        DataDI = LoadDiffData(Load)*10^-8;    % Function LoadDiffData defined at bottom.
% Convolute these predictions with U17 undulator spectrum        
        for i = 1:size(DataDI,3)
            DataDI(6,:,i) = ConvU17(DataDI(6,:,i),q,U17);
        end              
%  Pick the conformations that we want to work with from the vacuum result.
        for i = 1:size(DataDI,3)
            test2(i,:) = [corr2(interp1(q,DataDI(6,:,i).*SolventExclusionFactor,q1(corrStart:corrEnd)),Expt1(1,corrStart:corrEnd)),i]; 
            test2(isnan(test2))=0;   % removes any instances of NaN and replaces with 0.
        end
        [Xsorted(kk,:,:),Ysorted(kk,:,:)] = sort(test2,1,'descend');
                        % This sorts out the most strongly correlated
                        % sub-sets of the data. Want to save this for
                        % drawing figures & getting mean structures. 
% Use this information to decide what sub-trajectories to keep
        for i=1:KeepData
            DSgood(i,:) = DataDI(6,:,Ysorted(kk,i,1));
        end   
           dProtVac(kk,:) = mean(DSgood);  
% Load the cross-term between micelle and protein. 
        Load = [filename(1) filename(4)];
        DataDI = LoadDiffData(Load)*10^-8;   
% Convolute these predictions with U17 undulator spectrum        
        for i = 1:size(DataDI,3)
            DataDI(2,:,i) = ConvU17(DataDI(2,:,i),q,U17);
        end      
% Keep the information chosen above.         
        for i=1:KeepData
            DSgood(i,:) = DataDI(2,:,Ysorted(kk,i,1));
        end   
            dPmMfgvac(kk,:) = mean(DSgood);          
% Load the complementary cross-term between micelle and protein.        
        Load = [filename(3) filename(2)];
        DataDI = LoadDiffData(Load)*10^-8; 
% Convolute these predictions with U17 undulator spectrum        
        for i = 1:size(DataDI,3)
            DataDI(2,:,i) = ConvU17(DataDI(2,:,i),q,U17);
        end      
% Keep the information chosen above.      
        for i=1:KeepData
            DSgood(i,:) = DataDI(2,:,Ysorted(kk,i,1));
        end    
            dPmMfevac(kk,:) = mean(DSgood);               
% Now look into the micelle correction term - not saved separately as in bR
        Load = [filename(2) filename(1)];
        DataDI = LoadDiffData(Load); 
        [U,S,V] = svd(squeeze(DataDI(4,:,:)));
        dPMicelle(kk,:) = mean(V(:,1))*U(:,1)*S(1,1);          
end 
        [U6,S6,V6]  = svd(transpose(dPMicelle));  
        FlucMicel = mean(V6(:,1))*U6(:,1)*S6(1,1);
        FlucMicel =    FlucMicel*1.0000e-10; 
end

%% ===============================================================================================================
% Loop does the fitting of the data & searchers for optimal parameters 
for ll = 1:21
        Bfactor = 100; % Bfactor = 3 + (ll-11)/5;   Almost never need to vary B-factor. Not imporant any more if large. 
         BestScale2 =  23.6; % 23.4
         % BestScale2 = 23.4 + (ll-11)/10;  
         DampWidth = 0.198;  % 0.278 At what q-domain are the micelle effects damped. 
         % DampWidth = 0.198  + (ll-11)/250;    
         dDampWidth = 1;   % 0.4 How sharp is the low-q damping of the micelle
         % dDampWidth = 0.96 + (ll-11)/100; 
         BOGcorr = -0.4; 
         % BOGcorr = - 0.4 + (ll-11)/20;    % Correction for deviations from 190 BOG simulation. 
for kk = 1:size(stepnr,2)
%% ==========================================================================
 % Prepare matrices used for fitting.
  %  LowDamp = normpdf(q,0,DampWidth); LowDamp = LowDamp/max(LowDamp);
                            % Alternative Gaussian model for low-q damping term. 
      LowDamp = (1-erf((q-DampWidth)/(DampWidth*dDampWidth)))/2;
                            % Use error function model for low-q damping term. 
     AllDamp = normpdf(q,0,Bfactor); AllDamp = AllDamp/max(AllDamp);  
                            % For all intents & purposes = 1 so could be removed. 
     dProtMemVacA(kk,:) = (dPmMfevac(kk,:) + dPmMfgvac(kk,:))/2;               
                             % Cancel the membrane-protein cross-term fluctuations as in Sarabi et al.  
     dMemVac(kk,:) = (dProtMemVacA(kk,:) + BOGcorr*FlucMicel(2,:) - dProtVac(kk,:)).*LowDamp;                        
                            % Subtract the protein component, add micell fluctuation term, and damp the lot at low-q.
     dProtMemVac(kk,:) =  AllDamp.*(dProtVac(kk,:) + dMemVac(kk,:));     
                             % Add the damped micell back to the protein-vacuum term to get the curve used for fittting.   
%% =====================================================================================
                             %  Interpolate theoretical curves onto the data to allow fitting. 
    dProtMemVacInterp(:,kk) = interp1(q,dProtMemVac(kk,:),q1(FitStart:FitEnd)); 
    Contrast = transpose(interp1(q,SolventExclusionFactor,q1(FitStart:FitEnd)));
                             % Shorten the experimental vectors to the chosen fitting lengths. 
    ExptData2 = transpose(Expt1(1,FitStart:FitEnd)); 
                              % Minimize against experimental data.  
    weight = q1(FitStart:FitEnd).^2.*gaussmf(q1(FitStart:FitEnd),[0.4 -0]);  
    weight = transpose(weight/max(weight));   % Normalizes the weighting function. 
    weight = weight./weight; % This line removes the weighting term to yield the R-factor. 
                              % Creates a weighting function for fitting. But perhaps not used?
                              % close all; plot(q1(FitStart:FitEnd),weight)
                              % score the correlation function. 
                score2(kk) =  corr2(dProtMemVacInterp(:,kk).*Contrast,ExptData2);     
% Find the best scaling values for the most correlated componet.             
                    fun = @(x)(sum((x(1)*dProtMemVacInterp(:,kk).*Contrast.*weight - ExptData2.*weight).^2));               
                    x0 = [20 0.0005];  x1 = fminsearch(fun,x0);     
% Write the correctly scaled function for this iteraction for plotting. 
            score300 = sum((ExptData2.*weight).^2); 
            Rfactor(kk) = 100*sqrt(sum((BestScale2*dProtMemVacInterp(:,kk).*Contrast.*weight - ExptData2.*weight).^2))/sqrt(score300);
% Check to see if it is an improvement, and if so then save values. 
            if  Rfactor(kk) < BestYet; 
                    BestYet = Rfactor(kk);  
                    BestScale2A = BestScale2; 
                    Best3(kk,:) = BestScale2*dProtMemVacInterp(:,kk).*Contrast;    
                    BestIndex = kk;       
                    BestDampWidth = DampWidth; 
                    BestDampSpread = dDampWidth;
                    BestBfactor = Bfactor; 
                    BestBOGcorrection = BOGcorr;
            end  
end 
end
% Print to screen the optimal values. 
fprintf('Best R-factor is %4.1f \n',min(Rfactor))
fprintf('Best correlation after all corrections is %4.1f \n',max(score2)*100)
% fprintf('Best B-factor is %4.2f \n',BestBfactor)
fprintf('Best damping midpoint is %4.3f \n',BestDampWidth)
fprintf('Best damping width is %4.3f \n',BestDampSpread)
fprintf('Best overall scaling factor is %4.2f \n',BestScale2A)
fprintf('Best BOG correction is %4.1f \n',BestBOGcorrection)

%% Select structures to analyze. 
ZplotRfactor = [Rfactor(1:12); Rfactor(13:24); Rfactor(25:36); Rfactor(37:48);  Rfactor(49:60)];      
       Rmin = min(min(ZplotRfactor)); 
ZplotWeighted = round(KeepData*exp(-(ZplotRfactor-Rmin)/BoltzmanWeight),0); 
ZplotW = [ZplotWeighted(1,1:12) ZplotWeighted(2,1:12) ZplotWeighted(3,1:12) ZplotWeighted(4,1:12) ZplotWeighted(5,1:12)]; 
    
clear TempTest; TempTest(1,1) = 0; 
for i = 1:size(stepnr,2)
    if ZplotW(1,i) > 0 
            for j = 1:ZplotW(1,i)
                TempTest(end+1,1) = [stepnr(i)];  % Ysorted(i,j,1)             
                TempTest(end,2) = [Ysorted(i,j,1)];       
            end 
    end     
end 
TempTest = TempTest(2:end,:); 
size(TempTest)
sum(ZplotW)
% save('Structures2Average_SRII.mat','TempTest');
dlmwrite('Structures2Average_SRIIHTRII_BS1.txt',TempTest, '\t');
Best4 = ZplotW*dProtMemVacInterp(:,:)';
Best4 = BestScale2*Best4.*Contrast'/sum(ZplotW);

%%
% figure('Position', [950 50 400 700])
figure('Position', [750 50 500 400])
% subplot(2,1,1)
 contour(ZplotRfactor,[0:10:200],'Linewidth',3)
xlabel('\delta (Helix E & F)')
ylabel('\gamma (Helix C, D & E)')
xticks([1:1:12]); xticklabels({'0','1/4','1/2','3/4','1','5/4','3/2','7/4','2','9/4','5/2','11/4'})
yticks([1:1:5]); yticklabels({'0','1/2','1','3/2','2'})
box on; set(gca,'linewidth',3); set(gca,'FontSize',18);
colorbar

% subplot(2,1,1)
%  heatmap(ZplotWeighted)
% subplot(2,1,2)
figure('Position', [750 250 500 400])
plot(q1(FitStart:FitEnd),q1(FitStart:FitEnd).*transpose(ExptData2),'o','LineWidth',0.5,'MarkerSize',4,'color','#7E2F8E');
hold on; 
plot(q1(FitStart:FitEnd),q1(FitStart:FitEnd).*(Best3(BestIndex,:)),'LineWidth',5,'color',[0.8500 0.3250 0.0980]);
plot(q1(FitStart:FitEnd),q1(FitStart:FitEnd).*(Best4(1,:)),'LineWidth',3,'color','k');   
axis([0.0 1.1 -5*10^-4 3*10^-4]);
xlabel('q (Å^{-1})')
ylabel('q \cdot \DeltaS(q)')
box on; set(gca,'linewidth',3); set(gca,'FontSize',18);
xticks([0:0.2:1]); grid on
ax = gca; 
ax.GridLineWidth = 1;

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
DI(1,:,:) = excited.protmem.intSol(:,1:496) - ground.protmem.intSol(:,1:496);       
DI(2,:,:)  = excited.protmem.intVac(:,1:496) - ground.protmem.intVac(:,1:496);       
DI(3,:,:) = excited.membrane.intSol(:,1:496) - ground.membrane.intSol(:,1:496);   
DI(4,:,:)  = excited.membrane.intVac(:,1:496) - ground.membrane.intVac(:,1:496);    
DI(5,:,:)  = excited.protein.intSol(:,1:496) - ground.protein.intSol(:,1:496);      
DI(6,:,:)  = excited.protein.intVac(:,1:496) - ground.protein.intVac(:,1:496);
DataTemp = DI; 
end
function U17conv = ConvU17(yi,yii,yiii) 
       U17 = yiii; 
       scaleU17 = sum(U17(:,2)); 
       U17CoM = sum(U17(:,1).*U17(:,2))/scaleU17;  
       Econv = 1+(U17(:,1)-U17CoM)/U17CoM; 
       q = yii; 
      for i = 1:round(size(U17,1)/5,0)-1
       TMP(i,:) = U17(5*i,2)*interp1(q*Econv(5*i,1),yi,q(1:135));
      end
      U17conv = [sum(TMP(:,:))/scaleU17  zeros(1,66)];      
%      close all
%      plot(q(1:135),U17conv*max(dProtVac(24,:))/max(U17conv),'color','k')
%      hold on
%      plot(q,dProtVac(24,:),'color','r')
end 

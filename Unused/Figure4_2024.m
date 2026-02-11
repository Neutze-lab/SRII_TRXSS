close all; 
clear all; 
binsize = 6; 
%% Read SRII basis spectrum 2 data in. 
SRIIBS2 = load("C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\BestFitSRIIBS2.mat"); 
FitSRIIBS2 = SRIIBS2.Data2save;

MeanDisplacementSRIIBS2 = mean(FitSRIIBS2,1);
ErrorSRIIBS2 = std(FitSRIIBS2);
tempsize = binsize*floor(size(MeanDisplacementSRIIBS2,2)/binsize);

MeanDisplacementSRIIBS2 = reshape(MeanDisplacementSRIIBS2(1:tempsize),binsize,[]); 
MeanDisplacementSRIIBS2 = mean(MeanDisplacementSRIIBS2,1); 
ErrorSRIIBS2 = reshape(ErrorSRIIBS2(1:tempsize),binsize,[]); 
ErrorSRIIBS2 = mean(ErrorSRIIBS2); 

%% Read SRII:HtrII basis spectrum 1 data in. 
SRIIHtrIIBS1 = load("C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\BestFitSRIIHtrIIBS1.mat");
FitSRIIHtrIIBS1 = SRIIHtrIIBS1.Data2save;

MeanDisplacementSRIIHtrIIBS1 = mean(FitSRIIHtrIIBS1,1);
ErrorSRIIHtrIIBS1 = std(FitSRIIHtrIIBS1);
MeanDisplacementHtrIIBS1 = MeanDisplacementSRIIHtrIIBS1(size(FitSRIIBS2,2):end);
ErrorHtrIIBS1 = ErrorSRIIHtrIIBS1(size(FitSRIIBS2,2):end);

MeanDisplacementSRIIHtrIIBS1 = reshape(MeanDisplacementSRIIHtrIIBS1(1:tempsize),binsize,[]); 
MeanDisplacementSRIIHtrIIBS1 = mean(MeanDisplacementSRIIHtrIIBS1,1); 
ErrorSRIIHtrIIBS1 = reshape(ErrorSRIIHtrIIBS1(1:tempsize),binsize,[]); 
ErrorSRIIHtrIIBS1 = mean(ErrorSRIIHtrIIBS1); 

tempsize2 = binsize*floor(size(MeanDisplacementHtrIIBS1,2)/binsize);
MeanDisplacementHtrIIBS1 = reshape(MeanDisplacementHtrIIBS1(1:tempsize2),binsize,[]); 
MeanDisplacementHtrIIBS1 = mean(MeanDisplacementHtrIIBS1,1); 
ErrorHtrIIBS1 = reshape(ErrorHtrIIBS1(1:tempsize2),binsize,[]); 
ErrorHtrIIBS1 = mean(ErrorHtrIIBS1); 


%% Read SRII:HtrII basis spectrum 2 data in. 
SRIIHtrIIBS2 = load("C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\BestFitSRIIHtrIIBS2.mat");
FitSRIIHtrIIBS2 = SRIIHtrIIBS2.Data2save;

MeanDisplacementSRIIHtrIIBS2 = mean(FitSRIIHtrIIBS2,1);
ErrorSRIIHtrIIBS2 = std(FitSRIIHtrIIBS2);
MeanDisplacementHtrIIBS2 = MeanDisplacementSRIIHtrIIBS2(size(FitSRIIBS2,2):end);
ErrorHtrIIBS2 = ErrorSRIIHtrIIBS2(size(FitSRIIBS2,2):end);

MeanDisplacementSRIIHtrIIBS2 = reshape(MeanDisplacementSRIIHtrIIBS2(1:tempsize),binsize,[]); 
MeanDisplacementSRIIHtrIIBS2 = mean(MeanDisplacementSRIIHtrIIBS2,1); 
ErrorSRIIHtrIIBS2 = reshape(ErrorSRIIHtrIIBS2(1:tempsize),binsize,[]); 
ErrorSRIIHtrIIBS2 = mean(ErrorSRIIHtrIIBS2); 

MeanDisplacementHtrIIBS2 = reshape(MeanDisplacementHtrIIBS2(1:tempsize2),binsize,[]); 
MeanDisplacementHtrIIBS2 = mean(MeanDisplacementHtrIIBS2,1); 
ErrorHtrIIBS2 = reshape(ErrorHtrIIBS2(1:tempsize2),binsize,[]); 
ErrorHtrIIBS2 = mean(ErrorHtrIIBS2); 


%% Figure for the paper 
figure('Position', [750 150 900 450])
box on; set(gca,'linewidth',3); set(gca,'FontSize',20);
plot(binsize*[1:size(MeanDisplacementSRIIBS2,2)],MeanDisplacementSRIIBS2,'o','linewidth',3,'color','#0072BD'); 
hold on
errorbar(binsize*[1:size(MeanDisplacementSRIIBS2,2)],MeanDisplacementSRIIBS2,ErrorSRIIBS2,'linewidth',1.2,'color','#0072BD'); 
plot(binsize*[1:size(MeanDisplacementSRIIBS2,2)],MeanDisplacementSRIIHtrIIBS1,'o','linewidth',3,'color','#7E2F8E'); 
errorbar(binsize*[1:size(MeanDisplacementSRIIBS2,2)],MeanDisplacementSRIIHtrIIBS1,ErrorSRIIHtrIIBS1,'linewidth',1.2,'color','#7E2F8E'); 
plot(binsize*[1:size(MeanDisplacementSRIIBS2,2)],MeanDisplacementSRIIHtrIIBS2,'o','linewidth',3,'color','#EDB120'); 
errorbar(binsize*[1:size(MeanDisplacementSRIIBS2,2)],MeanDisplacementSRIIHtrIIBS2,ErrorSRIIHtrIIBS2,'linewidth',1.2,'color','#EDB120'); 
plot(binsize*[230/binsize:(230+tempsize2-1)/binsize],MeanDisplacementHtrIIBS1,'o','linewidth',3,'color','#7E2F8E'); 
errorbar(binsize*[230/binsize:(230+tempsize2-1)/binsize],MeanDisplacementHtrIIBS1,ErrorHtrIIBS1,'linewidth',1.2,'color','#7E2F8E'); 
plot(binsize*[230/binsize:(230+tempsize2-1)/binsize],MeanDisplacementHtrIIBS2,'o','linewidth',3,'color','#EDB120'); 
errorbar(binsize*[230/binsize:(230+tempsize2-1)/binsize],MeanDisplacementHtrIIBS2,ErrorHtrIIBS2,'linewidth',1.2,'color','#EDB120'); 

axis([1 300 0 3.5]);
xlabel('Residue number')
ylabel('C\alpha rmsd (Å)')
box on; set(gca,'linewidth',3); set(gca,'FontSize',18);
LineThickness = 1.1; 
SRhelixColor = [0.6350 0.0780 0.1840 0.05]; 
rectangle('Position', [3 0 26-3 4], 'FaceColor', SRhelixColor,'linewidth',LineThickness)
rectangle('Position', [33 0 56-33 4], 'FaceColor',SRhelixColor,'linewidth',LineThickness)
rectangle('Position', [70 0 92-70 4], 'FaceColor',SRhelixColor,'linewidth',LineThickness)
rectangle('Position', [94 0 118-94 4], 'FaceColor',SRhelixColor,'linewidth',LineThickness)
rectangle('Position', [122 0 150-122 4], 'FaceColor',SRhelixColor,'linewidth',LineThickness)
rectangle('Position', [153 0 181-153 4], 'FaceColor',SRhelixColor,'linewidth',LineThickness)
rectangle('Position', [189 0 219-189 4], 'FaceColor',SRhelixColor,'linewidth',LineThickness)

rectangle('Position', [231 0 50-24 4], 'FaceColor',SRhelixColor,'linewidth',LineThickness)
rectangle('Position', [231+52-24 0 82-52 4], 'FaceColor',SRhelixColor,'linewidth',LineThickness)

% HELIX    1   1 GLY A    3  GLY A   26  1                                  24    
% HELIX    2   2 GLU A   33  LEU A   56  1                                  24    
% HELIX    3   3 ALA A   70  GLY A   92  1                                  23    
% HELIX    4   4 ASP A   94  VAL A  118  1                                  25    
% HELIX    5   5 GLU A  122  SER A  150  1                                  29    
% HELIX    6   7 SER A  153  GLY A  181  1                                  29    
% HELIX    7   8 THR A  189  LEU A  219  1                                  31    

% MD has residues 23 to 82
% HELIX   10 AB1 ALA B   24  THR B   50  1                                  27    
% HELIX   11 AB2 ASP B   52  LEU B   82  1                                  31    

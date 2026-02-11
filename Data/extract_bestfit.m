clear all;  
close all;
%% SRII BS2
% nelem = 460; % Number of protein structures in the selected folder. 
% addpath('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\SrII_BS2\SrII_BS2\SrII_BS2_dark')
% addpath('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\SrII_BS2\SrII_BS2\SrII_BS2_light')
% DarkDirectory = dir('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\SrII_BS2\SrII_BS2\SrII_BS2_dark\*.pdb');
% LightDirectory = dir('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\SrII_BS2\SrII_BS2\SrII_BS2_light\*.pdb');

%% SRII:HtrII BS1
nelem = 336;
addpath('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\SrIIHtrII_BS1\SrIIHtrII_BS1\SrIIHtrII_BS1_dark')
addpath('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\SrIIHtrII_BS1\SrIIHtrII_BS1\SrIIHtrII_BS1_light')
DarkDirectory = dir('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\SrIIHtrII_BS1\SrIIHtrII_BS1\SrIIHtrII_BS1_dark\*.pdb');
LightDirectory = dir('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\SrIIHtrII_BS1\SrIIHtrII_BS1\SrIIHtrII_BS1_light\*.pdb');

%% SRII:HtrII BS2
% nelem = 449;
% addpath('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\SrIIHtrII_BS2\SrIIHtrII_BS2\SrIIHtrII_BS2_dark')
% addpath('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\SrIIHtrII_BS2\SrIIHtrII_BS2\SrIIHtrII_BS2_light')
% DarkDirectory = dir('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\SrIIHtrII_BS2\SrIIHtrII_BS2\SrIIHtrII_BS2_dark\*.pdb');
% LightDirectory = dir('C:\Users\Richard\Desktop\Folders\Papers\TIME_RESOLVED\SRII_TRXSS\4Lucija\SrIIHtrII_BS2\SrIIHtrII_BS2\SrIIHtrII_BS2_light\*.pdb');

%% Now read in the pdb files from a specified directory & extract C-alpha values. 
    for i = 1:nelem
        infile = LightDirectory(i).name;
        fprintf('%s \n', infile)
% Read in all the pdb files in the selected light activated state directory
        pdb = pdbread(infile);
        coords1 = [[pdb.Model.Atom.X]' [pdb.Model.Atom.Y]' [pdb.Model.Atom.Z]'];  % note the ' is making a (complex conjugate) matrix transpose.  
        calphas1 = find(strcmp({pdb.Model.Atom.AtomName}','CA'));
        CalphaCoords1 = coords1(calphas1,:);
        CalphaTableLight(:,:,i) = CalphaCoords1(:,:); 
    end
 fprintf('\n Number of pdb files represented in the light-activated C-alpha matrix is %4.0f \n', size(CalphaTableLight,3))

     for i = 1:nelem
        infile = DarkDirectory(i).name;
        fprintf('%s \n', infile)
% Read in all the pdb files in the resting state simulation directory 
        pdb = pdbread(infile);
        coords2 = [[pdb.Model.Atom.X]' [pdb.Model.Atom.Y]' [pdb.Model.Atom.Z]'];  % note the ' is making a (complex conjugate) matrix transpose.  
        calphas2 = find(strcmp({pdb.Model.Atom.AtomName}','CA'));
        CalphaCoords2 = coords2(calphas2,:);
        CalphaTableDark(:,:,i) = CalphaCoords2(:,:); 
    end
 fprintf('\n Number of pdb files represented in the dark state C-alpha matrix is %4.0f \n', size(CalphaTableDark,3))   

%% Calculate distance between two C-alpha matrices
distance = @(xyz) sqrt(xyz(1)^2 + xyz(2)^2 + xyz(3)^2);
displacement = zeros(size(calphas1,1),1);
for j=1:nelem
    for i = 1:size(calphas1,1)
         displacement(i,j) = distance(CalphaTableLight(i,:,j)- CalphaTableDark(i,:,j));
    end
end

% Now correlate the pairs against the mean displacement to pick out the "best model" for drawing figures.
for k = 1:nelem
    temp2 = mean(displacement,2); 
    score2(k) = corr2(displacement(:,k),temp2); 
end
 [a b] = max(score2); 
 fprintf('\n Structures correlating best with the mean C-alpha displacement are pair number %4.0f \n', b)

%% ------------------------------------------------------------------------------------------------------
% Calculate Internal Distance Matrices for Light. 
for k = 1:nelem
    for i = 1:size(CalphaTableLight,1) 
        for j = 1:size(CalphaTableLight,1)
            InternalDistanceLight(i,j,k) = distance(CalphaTableLight(i,:,k) - CalphaTableLight(j,:,k));
        end
    end 
end
% Calculate Internal Distance Matrices for Dark. 
for k = 1:nelem
    for i = 1:size(CalphaTableDark,1) 
        for j = 1:size(CalphaTableDark,1)
            InternalDistanceDark(i,j,k) = distance(CalphaTableDark(i,:,k) - CalphaTableDark(j,:,k));
        end
    end 
end

% Subtract these matrices and project them
for k = 1:nelem
    DiffMatrices(:,:,k) = InternalDistanceLight(:,:,k) - InternalDistanceDark(:,:,k); 
    ProjectDiffMatrices(:,k) = mean(abs(DiffMatrices(:,:,k)),1);
end
MeanDiffMatrices = mean(DiffMatrices,3);

% Now correlate the pairs against the mean internal distance changes to pick out the "best model" for drawing figures.
for k = 1:nelem
    temp1 = squeeze(DiffMatrices(:,:,k));
    score(k) = corr2(temp1,MeanDiffMatrices); 
end
 [a b] = max(score); 
 fprintf('\n Structures correlating best with the mean Internal Distance Matrix are pair number %4.0f \n', b)

%% ======================================================================================================
%% Plot these displacements in an easily understood manner. 
figure('Position', [650 25 650 750])
subplot(2,1,1)
for j=1:nelem 
    plot([1:size(displacement,1)],displacement(:,j))
    hold on
end
% axis([1 219 0 1.5])
xlim([1 219])
title('C_\alpha displacements')
ylabel('Displacement of C_\alpha (Å)')
xlabel('Residue number')
box on; set(gca,'linewidth',2); set(gca,'FontSize',14);

subplot(2,1,2)
for j=1:nelem 
    plot([1:size(displacement,1)],ProjectDiffMatrices(:,j))
    hold on
end
% axis([1 219 0 1.5])
xlim([1 219])
title('Projection of Internal Distance Changes')
ylabel('Projected Distance Changes (Å)')
xlabel('Residue number')
box on; set(gca,'linewidth',2); set(gca,'FontSize',14);

%% SRII BS2 best fit - Write data to files. 
OutlierRemovedDisplacement = rmoutliers(displacement','median','ThresholdFactor',5);
% size(OutlierRemovedDisplacement)
% Data2save = mean(OutlierRemovedDisplacement,1);
Data2save = OutlierRemovedDisplacement; 
% save('BestFitSRIIBS2.mat','Data2save');
% save('BestFitSRIIHtrIIBS1.mat','Data2save');
% save('BestFitSRIIHtrIIBS2.mat','Data2save');

%% ============================================================================
close all; 
figure('Position', [650 25 650 750])
subplot(2,1,1)
for j=1:size(OutlierRemovedDisplacement,1) 
   % plot([1:size(displacement,1)],displacement(:,j))
     plot([1:size(displacement,1)],OutlierRemovedDisplacement(j,:))
    hold on
end
% xlim([1 219])
title('C_\alpha displacements')
ylabel('Displacement of C_\alpha (Å)')
xlabel('Residue number')
box on; set(gca,'linewidth',2); set(gca,'FontSize',14);

subplot(2,1,2)
for j=1:size(OutlierRemovedDisplacement,1) 
    plot([1:size(displacement,1)],OutlierRemovedDisplacement(j,:)-mean(OutlierRemovedDisplacement,1))
    hold on
end
%xlim([1 219])
title('C_\alpha displacements')
ylabel('Displacement of C_\alpha (Å)')
xlabel('Residue number')
box on; set(gca,'linewidth',2); set(gca,'FontSize',14);


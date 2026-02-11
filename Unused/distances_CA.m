clear all
close all
distance = @(xyz) sqrt(xyz(1)^2 + xyz(2)^2 + xyz(3)^2);
%% Read in continuous illumination file with two conformations
infile1 = "1aus.pdb";
pdb1 = pdbread(infile1);
fprintf("Loaded %s\n",infile1)
% CALCULATE DISTANCES
calphas1 = find(strcmp({pdb1.Model.Atom.AtomName}','CA'));
coords1 = [[pdb1.Model.Atom.X]' [pdb1.Model.Atom.Y]' [pdb1.Model.Atom.Z]'];
CAcoords1a = coords1(calphas1,:);
for i = 1:size(CAcoords1a,1) 
    for j = 1:size(CAcoords1a,1)
        d1(i,j) = distance(CAcoords1a(i,:) - CAcoords1a(j,:));
    end
end 

infile1 = "0054-01_refmac.pdb";
pdb1 = pdbread(infile1);
fprintf("Loaded %s\n",infile1)
% CALCULATE DISTANCES
calphas1 = find(strcmp({pdb1.Model.Atom.AtomName}','CA'));
coords1 = [[pdb1.Model.Atom.X]' [pdb1.Model.Atom.Y]' [pdb1.Model.Atom.Z]'];
CAcoords1b = coords1(calphas1,:);
for i = 1:size(CAcoords1b,1)
for j = 1:size(CAcoords1b,1)
    d2(i,j) = distance(CAcoords1b(i,:) - CAcoords1b(j,:));
end
end 

%% =====================================================================
close all
% Protein 1AUS.pdb - L subunit residues 333 to 337 missing and residues 1
% to 19. S subunit goes from 1 to 123.
zonA(1,:) = [1:439];   % L20 to L463 missing.... 
zonA(2,:) = [563:1001];  % M20 to M463
zonA(3,:) = [1125:1563];
zonA(4,:) = [1687:2125];

zonC(1,:) = [440:562]; 
zonC(2,:) = [1002:1124]; 
zonC(3,:) = [1564:1686]; 
zonC(4,:) = [2126:2248]; 

% Protein 0054-01_refmac.pdb 
zonB(1,:) = [1:439];   % L20 to L463
zonB(2,:) = [440:878];    % M20 to M463
zonB(3,:) = [879:1317];
zonB(4,:) = [1318:1756];

zonD(1,:) = [1757:1879]; 
zonD(2,:) = [1880:2002]; 
zonD(3,:) = [2003:2125]; 
zonD(4,:) = [2126:2248]; 

%% Compare with 1AUS
region1 = zonA(1,:); 
for i=1:4
region2 = zonB(i,:);
DiffDA = d1(squeeze(region1),squeeze(region1)) - d2(squeeze(region2),squeeze(region2)); 
temp10(i,:) = mean(mean(abs(DiffDA)));
ProjectDiffDA1(i,:) = mean(abs(DiffDA),1);
end 

region1 = zonC(1,:); 
for i=1:4
region2 = zonD(i,:);
DiffDA = d1(squeeze(region1),squeeze(region1)) - d2(squeeze(region2),squeeze(region2)); 
temp12(i,:) = mean(mean(abs(DiffDA)));
ProjectDiffDA3(i,:) = mean(abs(DiffDA),1);
end 

temp20 = [ProjectDiffDA1(1,:) ProjectDiffDA1(2,:) ProjectDiffDA1(3,:) ProjectDiffDA1(4,:) ProjectDiffDA3(1,:) ProjectDiffDA3(2,:) ProjectDiffDA3(3,:) ProjectDiffDA3(4,:)];
save('Bfactors.mat','temp20')

%% Compare internally
ct = 0;   % variable just to count for writing output. 
for i=1:4
    for j = i+1:4
        ct = ct + 1; 
        region1 = zonB(i,:);
        region2 = zonB(j,:);
        DiffDA = d2(squeeze(region1),squeeze(region1)) - d2(squeeze(region2),squeeze(region2)); 
        temp11(ct,:) = mean(mean(abs(DiffDA)));
        ProjectDiffDA2(ct,:) = mean(abs(DiffDA),1);
    end 
end

ct = 0;   % variable just to count for writing output. 
for i=1:4
    for j = i+1:4
        ct = ct + 1; 
        region1 = zonD(i,:);
        region2 = zonD(j,:);
        DiffDA = d2(squeeze(region1),squeeze(region1)) - d2(squeeze(region2),squeeze(region2)); 
        temp13(ct,:) = mean(mean(abs(DiffDA)));
        ProjectDiffDA4(ct,:) = mean(abs(DiffDA),1);
    end 
end

%% ==============================================================================
close all
LW = 1.0; 
residuesA = 20:332; 
residuesB = 338:438+20+5;
residuesC = 501:623

figure('Position', [40 140 700 550])
subplot(2,1,1)
    plot(residuesA,transpose(ProjectDiffDA1(1,1:313)),'color',[0, 0.4470, 0.7410],'linewidth',LW)
    hold on
    plot(residuesB,transpose(ProjectDiffDA1(1,314:439)),'color',[0, 0.4470, 0.7410],'linewidth',LW)
    plot(residuesC,transpose(ProjectDiffDA3(1,1:123)),'color',[0, 0.4470, 0.7410],'linewidth',LW)
    plot(residuesA,transpose(ProjectDiffDA1(2,1:313)),'color',[0.8500, 0.3250, 0.0980],'linewidth',LW)
    plot(residuesB,transpose(ProjectDiffDA1(2,314:439)),'color',[0.8500, 0.3250, 0.0980],'linewidth',LW)
    plot(residuesC,transpose(ProjectDiffDA3(2,1:123)),'color',[0.8500, 0.3250, 0.0980],'linewidth',LW)
    plot(residuesA,transpose(ProjectDiffDA1(3,1:313)),'color',[0.9290, 0.6940, 0.1250],'linewidth',LW)
    plot(residuesB,transpose(ProjectDiffDA1(3,314:439)),'color',[0.9290, 0.6940, 0.1250],'linewidth',LW)
    plot(residuesC,transpose(ProjectDiffDA3(3,1:123)),'color',[0.9290, 0.6940, 0.1250],'linewidth',LW)
    plot(residuesA,transpose(ProjectDiffDA1(4,1:313)),'color',[0.4940, 0.1840, 0.5560],'linewidth',LW)
    plot(residuesB,transpose(ProjectDiffDA1(4,314:439)),'color',[0.4940, 0.1840, 0.5560],'linewidth',LW)
    plot(residuesC,transpose(ProjectDiffDA3(4,1:123)),'color',[0.4940, 0.1840, 0.5560],'linewidth',LW)
    
axis([1 650 0 0.8])
ylabel('Internal distance change (Å)')
xlabel('Residue number')
box on; set(gca,'linewidth',2); set(gca,'FontSize',12);
title('Internal distance change on C\alpha  atoms','FontSize',14)

subplot(2,1,2)

    plot(residuesA,transpose(ProjectDiffDA2(1,1:313)),'color',[0, 0.4470, 0.7410],'linewidth',LW)
    hold on
    plot(residuesB,transpose(ProjectDiffDA2(1,314:439)),'color',[0, 0.4470, 0.7410],'linewidth',LW)
    plot(residuesC,transpose(ProjectDiffDA4(1,1:123)),'color',[0, 0.4470, 0.7410],'linewidth',LW)
    plot(residuesA,transpose(ProjectDiffDA2(2,1:313)),'color',[0.8500, 0.3250, 0.0980],'linewidth',LW)
    plot(residuesB,transpose(ProjectDiffDA2(2,314:439)),'color',[0.8500, 0.3250, 0.0980],'linewidth',LW)
    plot(residuesC,transpose(ProjectDiffDA4(2,1:123)),'color',[0.8500, 0.3250, 0.0980],'linewidth',LW)
    plot(residuesA,transpose(ProjectDiffDA2(3,1:313)),'color',[0.9290, 0.6940, 0.1250],'linewidth',LW)
    plot(residuesB,transpose(ProjectDiffDA2(3,314:439)),'color',[0.9290, 0.6940, 0.1250],'linewidth',LW)
    plot(residuesC,transpose(ProjectDiffDA4(3,1:123)),'color',[0.9290, 0.6940, 0.1250],'linewidth',LW)
    plot(residuesA,transpose(ProjectDiffDA2(4,1:313)),'color',[0.4940, 0.1840, 0.5560],'linewidth',LW)
    plot(residuesB,transpose(ProjectDiffDA2(4,314:439)),'color',[0.4940, 0.1840, 0.5560],'linewidth',LW)
    plot(residuesC,transpose(ProjectDiffDA4(4,1:123)),'color',[0.4940, 0.1840, 0.5560],'linewidth',LW)


    plot(residuesA,transpose(ProjectDiffDA2(5,1:313)),'color',[0.4660, 0.6740, 0.1880],'linewidth',LW)
    plot(residuesB,transpose(ProjectDiffDA2(5,314:439)),'color',[0.4660, 0.6740, 0.1880],'linewidth',LW)
    plot(residuesC,transpose(ProjectDiffDA4(5,1:123)),'color',[0.4660, 0.6740, 0.1880],'linewidth',LW)
    plot(residuesA,transpose(ProjectDiffDA2(6,1:313)),'color',[0.6350, 0.0780, 0.1840],'linewidth',LW)
    plot(residuesB,transpose(ProjectDiffDA2(6,314:439)),'color',[0.6350, 0.0780, 0.1840],'linewidth',LW)
    plot(residuesC,transpose(ProjectDiffDA4(6,1:123)),'color',[0.6350, 0.0780, 0.1840],'linewidth',LW)

axis([1 650 0 0.8])
ylabel('Internal distance change (Å)')
xlabel('Residue number')
box on; set(gca,'linewidth',2); set(gca,'FontSize',12);
title('Internal distance change on C\alpha  atoms','FontSize',14)
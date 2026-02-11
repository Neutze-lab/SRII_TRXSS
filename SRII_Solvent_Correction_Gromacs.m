% close all; clear all; 
filename(1) = ["Data\SrII-190BOG-protmem_intensities_protmem"]; 
filename(2) = ["Data\SrII-190BOG-watercut_intensities_protmem"]; 

Temp1 = load(filename(1)); q = Temp1.q;
DataProt = LoadTotalData(filename(1)); DataProt = mean(DataProt,3); 
DataWater = LoadTotalData(filename(2)); DataWater = mean(DataWater,3); 
% scalefactor = max(DataProt(3,:))/max(DataWater(2,:)); 
        % Here we divide the solvent exclusion Crysol predicts from the protein+micell by the solvent exclusion from MD simulations of water; 
        % A value of 1.2 means that Crysol predicts 20 % more water scattering than is happening. 
scalefactor = 0.7423;      % The water cage is somewhat too big. We need to down-scale the water scattering a bit.
        % Comes from ratio of 2 Å radius sphere (4331 water molecules) to the projection back to 0
        % Å radius (3215 water molecules) from Sarabi et al. 
DataWater = scalefactor*DataWater;

CosPhi = 0.5*(DataProt(2,:) + DataProt(3,:) - DataProt(1,:))./(sqrt(DataProt(2,:)).*sqrt(DataProt(3,:)));
        % Calculates cos-phi from Crysol since we cannot get it from the water-excluded simulation from Gromax without our own script. 
CosPhi2 = (1 - 0.053*(1 + erf((q-0.18)/0.043))/2); CosPhi2 = transpose(CosPhi2);
        % Create an approximation to Cos phi that is structure-less. Plotted in panel 1. 

Stotal2 = DataProt(2,:) + DataProt(3,:) -2*CosPhi.*(sqrt(DataProt(2,:)).*sqrt(DataProt(3,:)));
        % Should recover exactly Stotal from Crysol. 
Stotal2B = DataProt(2,:) + DataProt(3,:) -2*CosPhi2.*(sqrt(DataProt(2,:)).*sqrt(DataProt(3,:)));
        % Predicts the total intensity if the phase approximation CosPhi2 is used instead.  
Stotal3 = DataProt(2,:) + DataWater(2,:) -2*CosPhi2.*(sqrt(DataProt(2,:)).*sqrt(DataWater(2,:)));
        % Predicts the total intensity if the phase approximation CosPhi2 and the water-exclusion simulation are used. 

% Stotal3Filt = RemcoFilter(Stotal3./DataProt(2,:));
SolventExclusionFactor1B = 1 - sqrt(DataProt(3,:))./sqrt(DataProt(2,:)).*CosPhi; 
SolventExclusionFactor2B = 1 - sqrt(DataWater(2,:) )./sqrt(DataProt(2,:)).*CosPhi2; 
        % Approximations using expansion... which may not be valid. 
SolventExclusionFactor1 = Stotal2./DataProt(2,:); 
        % This recovers the solvent exclusion factor of Crysol. 
SolventExclusionFactor2 = Stotal3./DataProt(2,:); 
% SolventExclusionFactor2Filt = RemcoFilter(SolventExclusionFactor2); 

%% Write out the factor for solvent exclusion used in other script. 
save('SolventExclusionFactorSRII.mat','SolventExclusionFactor2');

%% Plot the figures. 
figure('Position', [950 20 450 750])
subplot(3,1,1)
plot(q(1:140),transpose(CosPhi(1:140)),'Linewidth',2)
axis([0 1.1 0.94 1.01]) 
hold on
plot(q(1:140),CosPhi2(1:140))
xlabel('q (Å^{-1})')
title('Cos Phi and approximation','FontSize',14)

subplot(3,1,2)
semilogy(q,DataProt(1,:),'Linewidth',3)
hold on
semilogy(q,Stotal2,'Linewidth',1.5)
semilogy(q,Stotal2B,'Linewidth',1)
semilogy(q,Stotal3,'Linewidth',1)
xlim([0 1.1])
xlabel('q (Å^{-1})')
ylabel('\DeltaS(q)')
title('Total scattering and approximations','FontSize',14)

subplot(3,1,3)
plot(q,1.8*SolventExclusionFactor2,'Linewidth',3)
hold on
% plot(q,SolventExclusionFactor2Filt,'Linewidth',2) 
plot(q,3.2*SolventExclusionFactor1,'Linewidth',1)
plot(q,1.5*SolventExclusionFactor1B,'Linewidth',1)
plot(q,SolventExclusionFactor2B,'Linewidth',2)
% plot(q,SolventExclusionFactor2*max(SolventExclusionFactor1(50:100))/max(SolventExclusionFactor2(5:100)) ,'Linewidth',2)
axis([0 1.1 -0.1 1])
title('Corrections to be used in other script','FontSize',14)

%% ===================================================================================
% Functions used to read or filter data. 

function DataGround = LoadTotalData(yi)
Temp1 = load(yi);
GI(1,:,:) = Temp1.intSol;      
GI(2,:,:) = Temp1.intVac;       
GI(3,:,:) = Temp1.solvScatt;   
DataGround = GI; 
end

function RWFilter = RemcoFilter(yi) 
    WO = 50; 
    signalwidth = 1;
    N = max(size(yi));    
    Yi2(1:2*N) = zeros(1,2*N);
    % Make the vector twice as long. 
    Yi2(1:N) = yi;
    % Calculate the fft and amplitude of the FFT vector.  
    Yft= fft(Yi2);
    YftA= sqrt((real(Yft)).^2 + (imag(Yft)).^2); 
   % Now create the filter function. 
   theta = [1:size(YftA,2)]; 
   Modulate = cos(0.5*(theta-N)*pi/WO);  
   Y1 = zeros(1,2*N);
   Y1(1:WO) = Modulate(N:N+WO-1); 
   Y1(2*N-WO+1:2*N) = Modulate(N-WO:N-1); 
   YftAm = Y1.*YftA; 
   % Now invert the fft. 
   RecovData1 = ifft(Yft); 
   RecovData2 = ifft(Y1.*(real(Yft) + sqrt(-1)*imag(Yft)));  
   RWFilter = real(RecovData2(1:N)); 
end

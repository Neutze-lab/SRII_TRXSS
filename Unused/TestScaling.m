clear all
temp1 = load("SRII_BS_2.mat");
temp1b = temp1.SRII_BS_2; 
temp2 = load("SRII_HtrII_BS_2.mat");
temp2b = temp2.SRII_HtrII_BS_2; 
q = temp1.qIN;

plot(q,temp1b)
hold on
plot(q,temp2b)

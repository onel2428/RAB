%This script allows generating a figure similar to [REF1, Fig.13]
%With the output data here, one can also generate a figure similar to [REF1, Fig.14]

%% System parameters
M = 4;        %Number of PB antennas
pT = 3;       %Number of transmit antennas at the PB
samples = 250;%Number of Monte Carlo samples (number of different deployment realizations)

%% Main section of the code script

NN = [1 2 4 8 16 32 64 128];
Full_CSI = zeros(1,length(NN));
RAB_NoPC = zeros(1,length(NN));
RAB_OptPC = zeros(1,length(NN));
RAB_AppPC = zeros(1,length(NN));

for z=5:length(NN)      %Results as a function of the number of EH devices
    N=NN(z);  
    
    full_csi = zeros(1,samples);
    rab_noPC = zeros(N,samples);
    rab_optPC = zeros(N,samples);
    rab_appPC = zeros(N,samples);
    for j=1:samples   
        [N j]  %checking progress of the running
        rand('seed',j)
        
        % Devices are randomly and uniformly distributed in a 10 m-radius
        % area around the PB
        radius = 10*sqrt(rand(N,1));
        theta = 2*pi*rand(N,1);
        
        beta = 1e-4*radius.^(-2);    %path loss model
        
        %LOS channel generation for each device
        h = zeros(N,M);
        h_rot = zeros(N,M,M);
        for i=1:N
            phi = -(0:(M-1))*pi*sin(theta(i));  %[REF1, eq.2]
            h(i,:) = exp(1i*pi/4)*exp(1i*phi);  %[REF1, eq.1 with varphi_0=pi/4]

            %For the rotary mechanism, the LOS channel changes for each beam
            %direction m
            for m=1:M
                phi_rot = -(0:(M-1))*pi*sin(theta(i)+pi*m/M)+mod(0:(M-1),2)*pi;
                h_rot(i,:,m) = exp(1i*pi/4)*exp(1i*phi_rot);                 
            end
        end
                        
        %Full-CSI (SDP-based) [REF2]
        [~, E] = CSIBeamf_SDP(M,N,h,beta);
        full_csi(j) = pT*E;  %minimum RF energy available available at the
                             %EH devices when using the full-CSI approach       
             
        % Optimal and suboptimal power control [REF1, S.IV-A]
        popt = power_control(radius,theta, M, pT);        

        %RAB performance for each EH device
        for i=1:N
            E=beta(i).*abs(sum(squeeze(h_rot(i,:,:))./sqrt(M))).^2;
            rab_noPC(i,j) = pT*mean(E);  % without power control
            
            rab_optPC(i,j) = sum(popt(1,:).*E*1e3); % without optimum power control  [REF1, eq.17]
            rab_appPC(i,j) = sum(popt(2,:).*E*1e3); % without sub-optimal power control [REF1, eq.19]
        end
    end
    
    % Results in dBm
    Full_CSI(z) = 10*log10(nanmean(full_csi)*1e3);  
    RAB_NoPC(z) = 10*log10(nanmean(min(rab_noPC,[],1))*1e3);  
    RAB_OptPC(z) = 10*log10(nanmean(min(rab_optPC,[],1))*1e3);  
    RAB_AppPC(z) = 10*log10(nanmean(min(rab_appPC,[],1))*1e3);  
end

%% Plot of the results
figure
set(gcf, 'Units', 'centimeters'); 
axesFontSize = 14;
legendFontSize = 12;
afFigurePosition = [2 7 14 10]; 
set(gcf, 'Position', afFigurePosition,'PaperSize',[14 10],'PaperPositionMode','auto'); % [left bottom width height], setting printing properties 

h1=semilogx(NN,Full_CSI,'-pr','LineWidth',1.5,'MarkerSize',7,'MarkerFaceColor','w');
hold on
h2=semilogx(NN,RAB_NoPC,'-x','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','Color',[255,140,0]/255);
hold on
h3=semilogx(NN,RAB_OptPC,'--.','LineWidth',1.5,'MarkerSize',15,'MarkerFaceColor','w','Color',[255,140,0]/255);
hold on
h4=semilogx(NN,RAB_AppPC,':.','LineWidth',1.5,'MarkerSize',15,'MarkerFaceColor','w','Color',[255,140,0]/255);

grid on
box on

set(gca,'FontSize',12)
xlabel('number of EH devices, $|\mathcal{S}|$','Interpreter', 'latex','fontsize',axesFontSize)
ylabel('average worst-case RF energy (dBm)','Interpreter', 'latex','fontsize',axesFontSize)


hl=legend([h1 h2 h3 h4],'full-CSI','$\mathrm{RAB}$ (without power control)',...
    '$\mathrm{RAB}$ using solution of LP in (17)','$\mathrm{RAB}$ using eq. (19)');
set(hl,'interpreter','latex','FontSize',legendFontSize);

xlim([1 128])
xticks(NN)
yticks(-30:4:-10)
ylim([-30 -10])


%% References:
%[REF1]   - O. L. A. López, H. Alves, S. Montejo-Sánchez, R. D. Souza and 
%           M. Latva-aho, "CSI-free Rotary Antenna Beamforming for Massive
%           RF Wireless Energy Transfer," in IEEE Internet of Things 
%           Journal, doi: 10.1109/JIOT.2021.3107222.
%[REF2]   - O. L. A. López, F. A. Monteiro, H. Alves, R. Zhang and M. Latva-Aho,
%           "A Low-Complexity Beamforming Design for Multiuser Wireless Energy 
%           Transfer," in IEEE Wireless Communications Letters, vol. 10,
%           no. 1, pp. 58-62, Jan. 2021, doi: 10.1109/LWC.2020.3020576.
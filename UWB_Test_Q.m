clc;
clear;
close all;

% Simulation Experiment Parameters
% fs : Sampling frequency of the pulses
% Ep : Used for adjusting the energy of the pulse
% tc : Width of the pulse
% Tf : Frame length reference (A Real Number)
% Nf : Pulse Repetition Info
% Td : Data pulse delay reference
% Nc : TH Code bin number reference. Usually taken as Nc=Tf
% AWGNSNR : SNR value for AWGN (in dB)
% Q : Interval

% Sample Experiment:
numberOfBits=1000;
fs=10e9;
Ep=1;
tc=1e-9;
Tf=25;
Nf=2;
Td=12;
Nc=Tf-Td; %Nc = Tf - delta
          %TH_Sequence = randi([1,13]);
AWGNSNR = 2;
Guard_Interval = 0*Tf;
Ts = Nf * Tf;

%Some Pre Calculations
ChipPoints = fs * tc;
DelayPoints = ChipPoints * Td;
SequenceLength = numberOfBits * Nf;
FramePoint = ChipPoints * Tf;
SymbolPoints = FramePoint * Nf;
SizeOfSignals = ChipPoints * (numberOfBits * Tf * Nf + (numberOfBits - 1) * Guard_Interval);


rng default



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pulse Choices

% 2 GHz Gaussian monopulse sampled at a rate of 100 GHz:
%{
fc = 2E9;  fs=100E9;              % center freq, sample freq
tc = gmonopuls('cutoff',fc);      % width of each pulse
t1  = -2*tc : 1/fs : 2*tc;
Y = gmonopuls(t1,fc);
%}
%fc = 1E9;              % center freq, sample freq
gmonopuls('cutoff',1/tc);      % width of each pulse
t_single_singal = -0.5e-9:1/fs:0.5e-9;
Impulse = gmonopuls(t_single_singal,1/tc);

% Gaussian Pulse 2nd Order Doublelets
%t_single_singal = -1e-9:1/fs:1e-9;
%Impulse = (1-(4*pi.*(t_single_singal.^2))/(tc/2.5)^2) .* exp(-2*pi.*(t_single_singal.^2)/(tc/2.5)^2) / sqrt(Ep);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get Power Spectral Density Of Impulse
N = length(Impulse);
xdft = fft(Impulse);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(Impulse):fs/2;


%Frequncy Threshold
%-185dB at 4GHz
%-260dB at 100Ghz
Cutoff_Frequency = 4e9;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Transmitted Pulsetrans

D1(1)=1;    % Set 1 to be initial Reference Locations
for i=2:numberOfBits
    D1(i)=D1(i-1)+(Ts*(i-1)-D1(i-1)-tc)+1;
    D1(i) = D1(i) +  Guard_Interval * (i-1);  %Add Guard Interval (2Frames)
end

for i=1:Nf
            myRand=randi([1,13]);
            TH_Sequence(i) = myRand;
end
if (Nf>0)
    k=1;
    for i=1:length(D1)
        for j=0:(Nf-1)
            newD1(k)=D1(i)+Tf*j-1+TH_Sequence(j+1)-1;
            %newD2(k)=D2(i)+Tf*j;
            k=k+1;
        end
    end
end

D1=newD1;
D2=D1+Td;
%D2=newD2; %Locations of the pulses with Nf

Dtr = D1' * tc;       % locations of the TR pulses in ns.
Ddata = D2' * tc;     % locations of the Data pulses in ns.

% Generate the TR Pulse Train
t  = 0 : 1/fs : numberOfBits*Ts*tc + (numberOfBits - 1) * Guard_Interval * tc; % signal evaluation time
SimulationTime=t;
PulseTran_Reference = pulstran(t,Dtr,Impulse,fs);

% Generate the Random Modulated Data and Data Pulses
DataSequence = randi([1 1000],1,numberOfBits);
DataSequence = mod(DataSequence,2);
Input_Sequence_UWB = zeros(1, SequenceLength);
%Repeated Data Impulses
for i = 1:numberOfBits
    for j = 1:Nf
        Input_Sequence_UWB(1,Nf*(i-1)+1:1:Nf*(i-1)+j) = DataSequence(1,i);
    end
end 

for i = 1:numberOfBits * Nf
    if(Input_Sequence_UWB(1,i) == 0)
        Input_Sequence_UWB(1,i) = -1;
    end
end

Ddata = [Ddata.' ; Input_Sequence_UWB]';

PulseTran_Data=pulstran(t,Ddata,Impulse,fs);
%---------------------
allPulses=PulseTran_Reference+PulseTran_Data;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Channel Model and Recevied Signals

% This model is based on LOS (0-4 m) channel measurements. (CM-1)
Lam = 0.02333;
lambda = 2.5;
Gam = 7.1;
gamma = 4.3; 
std_ln_1 = 3.4;
std_ln_2 = 3.4;
nlos = 0;
std_shdw = 3;
num_channels = 100;

[t_channel,h_channel,Tbar,Tba,RMS] = UWBCHANNEL(Lam,lambda,Gam,gamma,std_ln_1,std_ln_2,nlos,std_shdw,num_channels);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Iteration To Calculate BER-Q


Q_max = 20;
BER_Q_Conventional = zeros(1,Q_max);
BER_Q_Hybrid = zeros(1,Q_max);

for Q_Iteration  = 1:1:Q_max

    for channel_iteration = 1:1:num_channels

        
        
        
        
        %allpulses + Noise
        %Adding Noise
        h = 1;
        SNR_lin=10^(AWGNSNR/10);
        noisePower=sum(abs(h).^2)./SNR_lin;
        wNoise= sqrt(1*noisePower).*randn(1,length(allPulses));

        receive_signals_without_channel = allPulses + wNoise;
        
        
        
        %CM1 Received Signal
        CM1_channel_output = conv(receive_signals_without_channel,h_channel(1,num_channels)); 
        
        
        
        receive_signals_without_noise = CM1_channel_output(1,1:1:1+SizeOfSignals);
        
        Q_Interval = Q * ChipPoints/10*5;
        
        
        receive_signals = receive_signals_without_noise;
        
        %{
        %Recevied Signals + Noise
        %Adding Noise
        h = 1;
        SNR_lin=10^(AWGNSNR/10);
        noisePower=sum(abs(h).^2)./SNR_lin;
        wNoise= sqrt(1*noisePower).*randn(1,length(allPulses));

        receive_signals = receive_signals_without_noise + wNoise;
        %}



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Conventioanl TR Receiver

        %Conventional Lowpass FIR filter
        rng default
        Wn = (2/fs)*Cutoff_Frequency;
        Conventional_FIR_Filter = fir1(20,Wn,'low',kaiser(21,3));


        %Signal After Conventional FIR filter
        SACF = filter(Conventional_FIR_Filter,1,receive_signals);

        %{
        %Compare Received Signals before and after LPF
        figure
        plot(t,receive_signals,t,SACF);
        xlim([0 3e-7]);

        xlabel('Time (s)');
        ylabel('Amplitude');
        legend('Original Signal','Filtered Data');
        %}

        %Multiplication r(t-Td)*r(t) 
        %SACF = Signal After Conventional Filter
        SACF_Delayed = [zeros(1,DelayPoints),SACF];
        %Adding Delay and remove tail points
        SACF_Delayed = SACF_Delayed(1,1:SizeOfSignals+1);
        SACF_Multiplication = SACF.*SACF_Delayed;

        %Integration to form estimated data sequence per pulse
        %Integration Time: T_int
        %T_int = 1 chip duration
        %T_int = Q_Iteration * ChipPoints;
        SACF_Integrated = zeros(1, SequenceLength);

        for i = 1:SequenceLength
            temp1 = mod(i+Nf-1, Nf);
            temp2 = (i+Nf-1-temp1)/Nf;
            Start_Point = ((Guard_Interval+Ts)*(temp2-1) + Tf*temp1 + TH_Sequence(1,temp1+1)+Td)*ChipPoints;

            for j = 1:Q_Interval   
                if (Start_Point+j <= SizeOfSignals+1)
                    SACF_Integrated(1,i) = SACF_Integrated(1,i) + SACF_Multiplication(1,Start_Point+j);
                end 
            end
        end




        %Estimate data sequence per frame (From Nf*numberOfBits to numberOfBits)
        SACF_Estimated = zeros(1,numberOfBits);
        for i = 1:SequenceLength
            temp1 = mod(i+Nf-1, Nf);
            temp2 = (i+Nf-1-temp1)/Nf;
            if(SACF_Integrated(1,i) >=0)
                SACF_Estimated(temp2)=SACF_Estimated(temp2)+1;
            else
                SACF_Estimated(temp2)=SACF_Estimated(temp2)-1;
            end
        end
        for i = 1:numberOfBits
            if(SACF_Estimated(1,i) >=0)
                SACF_Estimated(1,i) = 1;
            else
                SACF_Estimated(1,i) = 0;
            end
        end

        %Calculate BER
        BER_Conventional = 0;
        for i = 1:numberOfBits
            if(SACF_Estimated(1,i) ~= DataSequence(1,i))
                BER_Conventional = BER_Conventional + 1;
            end
        end
        BER_Conventional = BER_Conventional/numberOfBits;
        if BER_Conventional > 0.5
            BER_Conventional = 0.5;
        end




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Hybrid TR Receiver

        %Matched FIR Filter
        Matched_Filter = PulseTran_Reference(1, 1+(TH_Sequence(1,1)-1)*ChipPoints: 1+(TH_Sequence(1,1))*ChipPoints);



        %SAMF = Signal After Matched Filter
        SAMF = filter(Matched_Filter,1,receive_signals);



        %Multiplication r(tao)*S(temp)(tao-t)
        SAMF_Delayed = [zeros(1,DelayPoints),SAMF];
        %SAMF = [SAMF,zeros(1,DelayPoints)];
        %Adding Delay and remove tail points
        SAMF_Delayed = SAMF_Delayed(1,1:1+SizeOfSignals);
        SAMF_Multiplication = SAMF.*SAMF_Delayed;



        SAMF_Integrated = zeros(1, numberOfBits);

        for i = 1:numberOfBits

                for k = 1:Nf
                Start_Point = ((Guard_Interval+Ts)*(i-1) + Tf*(k-1) + TH_Sequence(1,k)+Td)*ChipPoints;
                    for j = 1:Q_Interval   
                        if (Start_Point+j < SizeOfSignals+1)
                            SAMF_Integrated(1,i) = SAMF_Integrated(1,i) + SAMF_Multiplication(1,Start_Point+j);
                        end 
                    end
                end
        end


        for i = 1:numberOfBits
            if(SAMF_Integrated(1,i) >=0)
                SAMF_Integrated(1,i) = 1;
            else
                SAMF_Integrated(1,i) = 0;
            end
        end

        %Calculate BER
        BER_Hybrid = 0;
        for i = 1:numberOfBits
            if(SAMF_Integrated(1,i) ~= DataSequence(1,i))
                BER_Hybrid = BER_Hybrid + 1;
            end
        end
        BER_Hybrid = BER_Hybrid/numberOfBits;
        if BER_Hybrid > 0.5
            BER_Hybrid = 0.5;
        end



        %Record BER values in current loop


        BER_Q_Conventional(1,Q_Iteration) = BER_Q_Conventional(1,Q_Iteration) + BER_Conventional;
        BER_Q_Hybrid(1,Q_Iteration) = BER_Q_Hybrid(1,Q_Iteration) + BER_Hybrid;

    
    end
    
    
    
    
end

BER_Q_Conventional = BER_Q_Conventional./num_channels;
BER_Q_Hybrid = BER_Q_Hybrid./num_channels;


figure
semilogy(BER_Q_Conventional);


figure
title('Q')
xlabel('Q')
ylabel('Bit Error Probability')
semilogy(1:1:20,BER_Q_Conventional,'-o',1:1:20,BER_Q_Hybrid,'--*');

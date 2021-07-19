% QPSK receiver for over the air transmission..
%% Parameter initialization..
agc = comm.AGC;%automatic gain control
Downsampling = 2;%downsampling factor...
Upsampling = 4;
M = 4;%modulation order
Fs = 200e3;% 
Ts = 1/Fs;
FrameSize = 100;
MessageLength = 105;
BarkerLength = 13;
Filterorder = 4*10;
RCFiltSpan = 10;
Rolloff = 0.5;
RxBufferedFrames = 10;
% Square root raised cosine receive filter
hRxFilt = fdesign.decimator(Upsampling/Downsampling, ...
                'Square Root Raised Cosine', Upsampling, ...
                'N,Beta', Filterorder, Rolloff);
hDRxFilt = design(hRxFilt, 'SystemObject', true);
Filtercoef = hDRxFilt.Numerator; 
K = 1;
A = 1/sqrt(2);
% Look into model for details for details of PLL parameter choice. Refer equation 7.30 of "Digital Communications - A Discrete-Time Approach" by Michael Rice. 
PhaseErrorDetectorGain = 2*K*A^2+2*K*A^2; % K_p for Fine Frequency Compensation PLL, determined by 2KA^2 (for binary PAM), QPSK could be treated as two individual binary PAM
PhaseRecoveryGain = 1; % K_0 for Fine Frequency Compensation PLL
TimingErrorDetectorGain = 2.7*2*K*A^2+2.7*2*K*A^2; % K_p for Timing Recovery PLL, determined by 2KA^2*2.7 (for binary PAM), QPSK could be treated as two individual binary PAM, 2.7 is for raised cosine filter with roll-off factor 0.5
TimingRecoveryGain = -1; % K_0 for Timing Recovery PLL, fixed due to modulo-1 counter structure
CoarseCompFrequencyResolution = 50; % Frequency resolution for coarse frequency compensation
PhaseRecoveryLoopBandwidth = 0.01; % Normalized loop bandwidth for fine frequency compensation
PhaseRecoveryDampingFactor = 1; % Damping Factor for fine frequency compensation
TimingRecoveryLoopBandwidth = 0.01; % Normalized loop bandwidth for timing recovery
TimingRecoveryDampingFactor = 1; % Damping Factor for timing recovery
PostFilterOversampling = Upsampling/Downsampling;
fc = 450e6;
Gain = 70;
Decimationfactor = 100e6/Fs;
USRPFs = 1/Fs;
USRPFrameLength = Upsampling*FrameSize;
% FrameTime = USRPFrameLength/Fs;
% StopTime = 5;
%% Setting up the system objects...
Buffer = dsp.Buffer(FrameSize*2,FrameSize);
bbc = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1]; %bipolar barker code
ubc = ((bbc+1)/2)';%unipolar barker code..
temp = (repmat(ubc,1,2))';%repeating barker sequence for I-Q components...
Header = temp(:);

%Header = ubc; % for in-phase only constellation like dbpsk

qpsk = comm.QPSKModulator();
ModulatedHeader = step(qpsk,Header)';
qpskdemod = comm.QPSKDemodulator('PhaseOffset',pi/4);
order = 4;
% x = (0:order-1)';
% refConst = bpsk.step(x);
biterr = comm.ErrorRate;
Correlator = dsp.Crosscorrelator;
constdiagram = comm.ConstellationDiagram('SamplesPerSymbol',Upsampling/2, ...
    'XLimits',[-1 1],'YLimits',[-1 1],'SymbolsToDisplaySource',...
    'Property','SymbolsToDisplay',6000);
rxFilter = dsp.FIRDecimator(Downsampling,Filtercoef);

coarseFreqcomp = comm.CoarseFrequencyCompensator(...
    'Modulation','BPSK', ...
    'SampleRate',Fs * Upsampling/Downsampling, ...
    'FrequencyResolution',CoarseCompFrequencyResolution);

% coarseFreqcomp = QPSKCoarseFrequencyCompensator('ModulationOrder',order,...
% 	'CoarseCompFrequencyResolution', CoarseCompFrequencyResolution,...
% 	'SampleRate', Fs, 'DownsamplingFactor', Downsampling); %coarse Frequency Correction...

theta = PhaseRecoveryLoopBandwidth/(PhaseRecoveryDampingFactor + ...
                0.25/PhaseRecoveryDampingFactor)/PostFilterOversampling;
d = 1 + 2*PhaseRecoveryDampingFactor*theta+theta*theta;
K1 = (4*PhaseRecoveryDampingFactor*theta/d)/...
                (PhaseErrorDetectorGain*PhaseRecoveryGain);
K2 = (4*theta*theta/d)/...
                (PhaseErrorDetectorGain*PhaseRecoveryGain);
OldOutput = complex(0); % used to store past value
FineFreqCompensator = QPSKFineFrequencyCompensator(...
                'ProportionalGain', K1, ...
                'IntegratorGain', K2, ...
                'DigitalSynthesizerGain', -1*PhaseRecoveryGain);% Fine Frequency Correction..
% Refer C.57 to C.61 in Michael Rice's "Digital Communications 
% - A Discrete-Time Approach" for K1 and K2
theta = TimingRecoveryLoopBandwidth/...
                (TimingRecoveryDampingFactor + ...
                0.25/TimingRecoveryDampingFactor)/PostFilterOversampling;
d = 1 + 2*TimingRecoveryDampingFactor*theta + theta*theta;
K1 = (4*TimingRecoveryDampingFactor*theta/d)/...
                (TimingErrorDetectorGain*TimingRecoveryGain);
K2 = (4*theta*theta/d)/...
                (TimingErrorDetectorGain*TimingRecoveryGain);
TimingRec = QPSKTimingRecovery('ProportionalGain', K1,...
                'IntegratorGain', K2, ...
                'PostFilterOversampling', PostFilterOversampling, ...
                'BufferSize', FrameSize);% Timing Recovery...

%% Setting up the radio now..
% radio = comm.SDRuReceiver(...
%         'IPAddress', '192.168.90.2' , ...
%         'CenterFrequency',fc, ...
%         'Gain', Gain, ...
%         'DecimationFactor',100e6/Fs, ...
%         'SamplesPerFrame', USRPFrameLength, ...
%         'OutputDataType', 'double');

radio = sdrrx('Pluto', 'BasebandSampleRate', Fs,...
        'CenterFrequency',fc, ...
        'GainSource', 'Manual', ...
        'Gain', Gain, ...
        'SamplesPerFrame', USRPFrameLength, ...
        'OutputDataType', 'double');

currentTime = 0;

len = uint32(0);
BER = zeros(3,1);
corruptSignal = complex(zeros(4000,1));
temp = [];
[Count,Delay,Phase] = deal(0);
FrameIndex = 0;
SyncIndex = 0;
SyncFlag = true;
while true
    scatterplot(corruptSignal); pause(0.1);
    %scatterplot(fineCompSignal); pause(0.1);
	%Keep accesssing the SDRu system object output until it is valid...
	while len <= 0
		[corruptSignal, len] = step(radio);
	end
	if len > 0
        % Apply automatic gain control to the signal...
        AGCSignal = (1/sqrt(Upsampling))*step(agc,corruptSignal);
        % Pass the signal through square root raised cosine received
        % filter
        FiltSignal = step(rxFilter,AGCSignal);
        % Coarsely compensate for the frequency offset..
        coarseCompSignal = step(coarseFreqcomp,FiltSignal);
        coarseCompBuffer = complex(zeros(size(coarseCompSignal)));
        timingRecBuffer = complex(zeros(size(coarseCompSignal)));
        for i=1:length(coarseCompSignal)
            % Scalar processing for fine frequency compensation and timing
            % recovery
            fineCompSignal = step(FineFreqCompensator,...
                [OldOutput coarseCompSignal(i)]);
            coarseCompBuffer(i) = fineCompSignal;
            OldOutput = fineCompSignal;
            % Timing recovery of the received signal...
            [dataOut, isDataValid, timingRecBuffer(i)] = step(TimingRec, fineCompSignal);
            if isDataValid
                %Decoding the received data...
                rxData = step(Buffer,dataOut);
                % Get a frame of data aligned on the frame boundary
                Data = rxData(Delay+1:Delay+length(rxData)/2);
                % Phase Estimation..
                y = mean(conj(ModulatedHeader).*Data(1:BarkerLength));

                %compensating for the phase offset...
                if Data(1)~=0
                    ShiftedData = Data.*exp(-1j*Phase);
                else
                    ShiftedData = complex(zeros(size(Data)));
                end
                % Demodulate the phase recovered data
                demodOut = step(qpskdemod,ShiftedData);
                % Recover the message from the data..
                Received = demodOut(BarkerLength*log2(order)+1:FrameSize*log2(order));
                receivedmsg = Received(1:MessageLength);
                %Finding delay to achieve frame synchronization...
                z = abs(step(Correlator, ModulatedHeader, dataOut));
                [~,ind] = max(z);
                Delay = mod(length(dataOut)-ind,(length(dataOut)-1));
                %phase ambiguity correction...
                Phase = round(angle(y)*2/pi)/2*pi;
                % Print received frame and estimate the received frame index..
                msgascii = bin2text(receivedmsg');
                %if strcmp(msgascii,'Hello Worl')
                disp(msgascii);
                %end
            end
        %}
        end
    end
end








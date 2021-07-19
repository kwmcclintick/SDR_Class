% QPSK Transmitter with barker sequence..
%% Parameter Initialization...
Fs = 200e3;%sampling frequency
Upsampling = 4;%Upsampling factor..
Ts = 1/Fs;%sampling time..
FrameSize = 100;%Number of modulated symbols per frame
MessageLength = 105;
%FrameCount = 100;
RxBufferedFrame = 10;
RCFiltSpan = 10;
Rolloff = 0.5;
Filterorder = 10*4;
fc = 450e6;
Gain = 25;
USRPFramelength = Upsampling*FrameSize*RxBufferedFrame;
Frametime = USRPFramelength/Fs;
%Stoptime = 1000;

%% Bits Generation with Barker Sequence..
message = 'Hello World 000';
bbc = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];%bipolar barker code
ubc = ((bbc+1)/2)';%unipolar barker code..
temp = (repmat(ubc,1,2))';%repeating barker sequence for I-Q components...
Header = temp(:);
data = [Header;text2bin(message)';zeros(1,69)'];

%% Modulating the data
dbpsk = comm.DBPSKModulator();
hTxFilt = fdesign.interpolator(4,'Square Root Raised Cosine',4,...
    'N,Beta',40,0.5);
txFilterdes = design(hTxFilt,'SystemObject',true);
filtercoefficients = txFilterdes.Numerator/2;
txFilter = dsp.FIRInterpolator(Upsampling,filtercoefficients);
modulateddata = dbpsk(data);
transmittedsignal = txFilter(modulateddata);
% Visualization..
cdiag = comm.ConstellationDiagram('XLimits',[-1 1],'YLimits',[-1 1]);
cdiag(transmittedsignal);%QPSK constellation..

%% Transmitting packets over the air..
radio = sdrrx('Pluto', 'BasebandSampleRate', Fs,...
        'CenterFrequency',fc, ...
        'SamplesPerFrame', USRPFrameLength, ...
        'OutputDataType', 'double');
                         
currenttime = 0;
% Starting the transmission...
i=1;
while currenttime <Stoptime
   step(radio,transmittedsignal);
   disp(strcat('Frame# ',num2str(i),'done'));
   %currenttime = currenttime+Frametime;
   %i = i+1;
end







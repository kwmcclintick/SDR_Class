
Fs = 200e3; %sampling frequency
usr_carrier = 50e3; %another PLUTO user interfering
BW = 4e3; %user bandwidth
energy_conv = 0.5; %50% energy convergence
energy_threshold = 1; %amount of energy a carrier needs to be below to be desireble
numFrames = 10;
numSpurriousNoise = 5; %each time frame will have 3 time-varying ISM users

sym_T = 0.1; %half energy point is within 10% of the Fc

usr_freqs = usr_carrier-BW:0.1e3:usr_carrier+BW;
carriers = [902.15 902.35 902.55 902.75 902.95]*1e6; %carrier frequencies to be used

energy_per_frame = zeros(1, length(carriers));
energy_history = zeros(1, length(carriers));
for frame=1:numFrames
    t = 0:1/Fs:1; %time indecies over one second
    x = zeros(1, length(t));
    
    %randomly puts in some noise simulating other ISM users
    tv_noise_carriers = randsample(length(carriers),numSpurriousNoise)';
    noise_freqs = carriers(tv_noise_carriers(1))- ...
        BW:0.1e3:carriers(tv_noise_carriers(1))+BW;
    for j=2:length(tv_noise_carriers)
        noise_freqs = [noise_freqs carriers(tv_noise_carriers(j))- ...
            BW:0.1e3:carriers(tv_noise_carriers(j))+BW];
    end

    for i=1:length(usr_freqs)
        %add PLUTO usr noise, constant over frames
        x = x + cos(2*pi*t*usr_freqs(i));
    end
    for i=1:length(noise_freqs)
        %add ISM band noise, changes each frame
        x = x + randn(1,1)*cos(2*pi*t*noise_freqs(i));
    end
    
    x = x + randn(size(t)); %add AWGN noise
    nfft = 2^nextpow2(length(x)); %number of bins in the FFT
    Pxx = abs(fft(x,nfft)).^2/length(x)/Fs;
    Hpsd = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs); %PSD
    
    %find how much energy is in each carrier's BW
    energy_levels = zeros(1, length(carriers));
    for i=1:length(carriers)
        energy_levels(i) = sum(Hpsd.Data((i-1)*length(Hpsd.Data) ...
            /length(carriers) + 1:i*length(Hpsd.Data)/length(carriers)));
    end
    
    %update energy history to see if energy is converging over time
    if frame == 1
       energy_history = energy_levels;
    else
       energy_history = [energy_history; energy_levels]; 
    end
    
    for i=1:length(energy_per_frame)
        energy_per_frame(i) = mean(energy_history(i));
    end
    
    %PICK LOWEST E CARRIER
    sym = false;
    while sym == false
        
        if min is sym
            pick
            sym true
            index + 1=
        end
    end
    current_best_carrier = find(min(energy_levels));
    
    %PLOTTING
    low_energy_carriers = energy_levels < energy_threshold;
    subplot(2,1,1)
    plot(Hpsd); xlabel('Frequency (kHz)'); title('PSD');
    subplot(2,1,2)
    stem(low_energy_carriers); xlabel('Carrier'); ylabel('Low E?');
    pause(0.1);
end
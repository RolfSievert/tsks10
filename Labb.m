%% Load and plot signal
filename = 'signal-rolsi701.wav';
[y, Fs] = audioread(filename);

L = length(y);
T = 1/Fs;

% Plot in time domain
t = (0:L-1)*T;
plot(t, y);
xlabel('t (s)');
ylabel('y(t)');

%% Transform and plot
Y = fft(y);
% shift moves down the transform to be centered around origo?
Y_shift = fftshift(fft(y));
Y_abs = abs(Y_shift/L);

f1 = 45500;
f2 = 45501;

% Plot
% Shifted half length to the left to plot correct
f = Fs*(-L/2+1:L/2)/L;
plot(f, Y_abs);
xlabel('f (Hz)');
ylabel('abs(Y)/2');

%% Plot cut
% L/2 = 2600000
% frekvenser: 17, 36, 55, 74, 93, 112, 131, 150 kHz.
%carry_f = 55*10^3;
carry_f = 55*10^3;
% Hz per sample
Hz_per_sample = L/Fs;
% Bandwidth in Hz
B_f=9000;
% Cut start and end in frequencies
cut_start_f = carry_f - B_f;
cut_end_f = carry_f + B_f;

% Center of sampled list
center=L/2;

% Get Frequency per sample and multiply with carry to get index in list
cut_start = cut_start_f*Hz_per_sample;
cut_end = cut_end_f*Hz_per_sample;
cut_start_int = uint32(cut_start);
cut_end_int = uint32(cut_end);

% Bandpass
bp = zeros(1, L);
% fill the uncut with original signal frequencies
bp(center-cut_end_int:center-cut_start_int) = Y_abs(center-cut_end_int:center-cut_start_int);
bp(center+cut_start_int:center+cut_end_int) = Y_abs(center+cut_start_int:center+cut_end_int);

% Plot
%f = Fs*(1:L)/L;
plot(f, bp);
xlabel('f (Hz)');
ylabel('abs(Y)/2');

%% Apply filter to signal and plot
% filtered = bandpass(y, [140*10^3 160*10^3], Fs);
%cutoff = [cut_start/(L/2) cut_end/(L/2)];
carry_f = 55*10^3;
B_f = 9000;
cutoff = [carry_f - B_f/2 carry_f + B_f/2]/(Fs/2);
[b, a] = butter(3, cutoff, 'bandpass');
IQ = filter(b, a, y);

% Plot in time domain
t = (0:L-1)*T;
plot(t, IQ);
xlabel('t (s)');
ylabel('y(t)');

%% Remove eco
[corr, lags] = xcorr(IQ);
% Times T because xcorr returns samples
plot(lags*T, corr);
xlim([-1 1]);

% Eco is 0.38 seconds, 152000th sample
eco_t = 0.38;
% eco in samples
eco_s = eco_t*Fs;

for index = (1:L-1 - eco_s)
    IQ(index + eco_s) = IQ(index + eco_s) - 0.9*IQ(index);
end

%% IQ
% Bandwidth in samples
%B=B_f*Hz_per_sample;
fc = B_f/Fs;
%fc = carry_f*T;
[b, a] = butter(3, fc, 'low');
phase_shift = pi/4;
XI = filter(b, a, 2*IQ.*cos(2*pi*carry_f*t' + phase_shift));
XQ = -filter(b, a, 2*IQ.*sin(2*pi*carry_f*t' + phase_shift));


%% Sample down and play
down_scale = 45;
ds = downsample(XI, down_scale);
%ds = downsample(XQ, down_scale);
soundsc(ds, Fs/down_scale);
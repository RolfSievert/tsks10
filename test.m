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
%Y = fft(y);
% shift moves down the transform to be centered around origo?
Y = fft(y);
Y_abs = abs(Y/L);

% Plot
% Shifted half length to the left to plot correct
f = Fs*(-L/2+1:L/2)/L;
plot(f, Y_abs);
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
%%
cut = (B_f/2)/(Fs/2);
[b,a] = butter(3, cut, 'low');
phase = pi / 4;
I = filter(b,a, 2*IQ.*cos(phase + 2*pi*carry_f*t') );
Q = -filter(b,a, 2*IQ.*sin(phase + 2*pi*carry_f*t') );
%play(Q, Fs);
%% IQ
% Bandwidth in samples
%B=B_f*Hz_per_sample;
fc = (B_f)/(Fs);
%fc = carry_f*T;
[b, a] = butter(3, fc, 'low');
phase_shift = pi/4;
XI = filter(b, a, 2*IQ.*cos(2*pi*carry_f*t' + phase_shift));
XQ = -filter(b, a, 2*IQ.*sin(2*pi*fc*t' + phase_shift));
%%
downscale = 45;
soundsc(downsample(Q, downscale), Fs/downscale);
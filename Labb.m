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
Y_shift = fftshift(Y);
Y_abs = abs(Y_shift/L);

f1 = 45500;
f2 = 45501;

f = Fs*(-L/2+1:L/2)/L;
plot(f, Y_abs);
xlabel('f (Hz)');
ylabel('abs(Y)/2');

%% Apply filter to signal and plot
carry_f = 55*10^3;
B=10000;
cutoff = [carry_f - B/2 carry_f + B/2]/(Fs/2);
[b, a] = butter(3, cutoff, 'bandpass');
IQ = filter(b, a, y);

t = (0:L-1)*T;
plot(t, IQ);
xlabel('t (s)');
ylabel('y(t)');

%% Remove eco
[corr, lags] = xcorr(IQ);
plot(lags*T, corr); % Times T because xcorr returns samples
xlim([-1 1]);

eco_t = 0.38;
eco_s = eco_t*Fs;

for index = (1:L-1 - eco_s)
    IQ(index + eco_s) = IQ(index + eco_s) - 0.9*IQ(index);
end

%% IQ
% Bandwidth in samples
B_s = B/Fs;
[b, a] = butter(3, B_s, 'low');
phase_shift = pi/4;
I = filter(b, a, 2*IQ.*cos(2*pi*carry_f*t' + phase_shift));
Q = -filter(b, a, 2*IQ.*sin(2*pi*carry_f*t' + phase_shift));

%% Sample down and play
scale = 43;
Q_sound = downsample(Q, scale);
I_sound = downsample(I, scale);
%%
soundsc(Q_sound, Fs/scale);
%%
soundsc(I_sound, Fs/scale);
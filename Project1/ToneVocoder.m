% #ToneVocoder#
clear;clc;

% Import the sound_read

[sound_read, fs] = audioread('Sound.wav');
t = (0 : (length(sound_read) - 1)) / fs;


% Veriables Setting

SSN = 'N';  % Add SNR 'Y' or not 'N'
Mode = 1;  % Equal frequency: 0 | Average cochlea length: 1
N_sub = 2;  % Number of segments
N_BPF = 4;  % BPF order
Cf = 50;  % Cut-off frequence of LPF
N_LPF = 4;  % LPF order


% Tone Vocoder

% Add SNR
if SSN == 'Y'
    sound_read = getSSN(sound_read, fs);
end

% Passband Mode
% Mode 0 is Equal frequency
if Mode == 0

    Nmax = N_sub;
    fmin = 200;
    fmax = 7000;

    % Generate an array of buffer filter coefficients
    b = zeros(Nmax, 2 * N_BPF + 1);
    a = zeros(Nmax, 2 * N_BPF + 1);
    gap = (fmax - fmin) / Nmax;
    for j = 1 : Nmax
        [b(j, :), a(j, :)] = butter(N_BPF, [fmin + (j-1)*gap fmin + j*gap] / fs * 2);
    end

    % Generate sub-bands
    cuted = zeros(Nmax, length(sound_read));
    for j = 1 : Nmax
        cuted(j, :) = filter(b(j, :), a(j, :), sound_read);
    end

    % Generate the Envelope
    envelopes = Envelope(cuted(:, :), fs, N_LPF, Cf);

    % Carrier loaded
    loaded = zeros(Nmax, length(sound_read));
    gap = (fmax - fmin) / Nmax;
    for j = 1 : Nmax
        for k = 1 : length(sound_read)
            loaded(j, k) = envelopes(j, k) .* cos(2*pi * (fmin + (j-1/2)*gap) * (k-1) / fs);
        end
    end

    % Synthesis and normalization
    output = zeros(1, length(sound_read));
    for j = 1 : Nmax
        output = output + loaded(j, :);
    end

    % Display the output
    sound(output, fs);
    figure;
    plot(t, output);
    if SSN == 'Y'
        audiowrite(['./Output/EqualFrequency_WithSSN_',num2str(N_sub),'Bands_',num2str(Cf),'Hz.wav'], output, fs);
    end
    if SSN == 'N'
        audiowrite(['./Output/EqualFrequency_',num2str(N_sub),'Bands_',num2str(Cf),'Hz.wav'], output, fs);
    end


    % Mode 1 is Average cochlea length
elseif Mode == 1

    dmin = (50 * log(1827 / 827)) / (3 * log(10));
    dmax = (50 * log(35827 / 827)) / (3 * log(10));
    Nmax = N_sub;

    % Generate an array of buffer filter coefficients
    b = zeros(Nmax, 2 * N_BPF + 1);
    a = zeros(Nmax, 2 * N_BPF + 1);
    gap = (dmax - dmin) / Nmax;
    for j = 1 : Nmax
        [b(j, :), a(j, :)] = butter(N_BPF, [alter(dmin + (j-1)*gap) alter(dmin + j*gap)] / fs * 2);
    end

    % Generate sub-bands
    cuted = zeros(Nmax, length(sound_read));
    for j = 1 : Nmax
        cuted(j, :) = filter(b(j, :), a(j, :), sound_read);
    end

    % Generate the Envelope
    envelopes = Envelope(cuted(:, :), fs, N_LPF, Cf);

    % Carrier loaded
    loaded = zeros(Nmax, length(sound_read));
    gap = (dmax - dmin) / Nmax;
    for j = 1 : Nmax
        for k = 1 : length(sound_read)
            loaded(j, k) = envelopes(j, k) .* cos(2*pi * alter(dmin + (j-1/2)*gap) * (k-1) / fs);
        end
    end

    % Synthesis and normalization
    output = zeros(1, length(sound_read));
    for j = 1 : Nmax
        output = output + loaded(j, :);
    end

    % Display the outputs
%     sound(output, fs);
    figure;
    subplot(2, 1, 1), plot(t, output);

    ak = fftshift(abs(fft(output)));
    w = linspace(-7000, 7000, length(abs(fft(output))));
    subplot(2, 1, 2), plot(w, ak);

    % Save
    if SSN == 'Y'
        audiowrite(['./Output/Cochlea_WithSSN_',num2str(N_sub),'Bands_',num2str(Cf),'Hz.wav'], output, fs);
    end
    if SSN == 'N'
        audiowrite(['./Output/Cochlea_',num2str(N_sub),'Bands_',num2str(Cf),'Hz.wav'], output, fs);
    end

end


% Functions


% Generate SSN
function Out = getSSN(sound_read, fs)
sound_read = sound_read';
nfft = 512;
noverlap = nfft / 2;
Windows = hamming(nfft);
% Returns the power spectral density
[Pxx,w] = pwelch(repmat(sound_read, 1, 10), Windows, noverlap, nfft, fs);
% Generate filter
b = fir2(3000, w/(fs / 2), sqrt(Pxx / max(Pxx)));
% Generate white noise
N = length(sound_read);
noise = 1 - 2 * rand(1, N + length(b) - 1);
% Generate SSN
ssn = filter(b, 1, noise);
ssn = ssn(length(b) : end);
% Adjust the energy of SSN
ssn = ssn / norm(ssn) * norm(sound_read) * 10 ^ (1 / 4);
% Addition sound_read with ssn
Out = sound_read + ssn;
% Normalize out to sound_read
Out = Out / norm(Out) * norm(sound_read);
end

% Batch Extract Envelopes
function E = Envelope(Y, fs, N1, cf)
[b, a] = butter(N1, cf / (fs / 2), "low");
sizeY = size(Y);
m1 = sizeY(1);
m2 = sizeY(2);
Y = abs(Y);
E = zeros(m1, m2);
for n = 1 : m1
    E(n, :) = filter(b, a, Y(n, :));
end
end

% f-d Alter
function f = alter(d)
f=165.4*(10^(0.06*d)-1);
end
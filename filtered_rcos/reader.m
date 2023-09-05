filename = 'snr.bin';
fid = fopen(filename, 'rb');
data = fread(fid, 'float32');
fclose(fid);

fs = 32000;
N = length(data);
freq = linspace(0, fs/2, N/2+1);
fft_data = fft(data);
mag_data = abs(fft_data(1:N/2+1));
mag_data_dB = 20*log10(mag_data);

figure();
plot(freq, mag_data_dB);
grid on;











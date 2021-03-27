[data, fs] = audioread('crazilyLong.wav');
spectrum = 20 * log10( abs( fft(data(:, 1)) ));
NFFT = length(spectrum);
if mod(NFFT,2)==0
    Nout = (NFFT/2)+1;
else
    Nout = (NFFT+1)/2;
end
spectrum = spectrum(1:Nout);
freqs = ((0:Nout-1)'./NFFT).*fs;
Noct = 10.6;
noctaveSmoothStruct = james_initOctave(Noct, freqs);
ar3 = james_smoothOctave(noctaveSmoothStruct, spectrum');
plot(freqs, spectrum);
hold on;
plot(freqs, ar3);
hold off;
legend('Original', '1/N octave');
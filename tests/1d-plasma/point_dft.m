function [dft_real, dft_imag] = point_dft(dft_real, dft_imag, ex, freqs, t)

for idx = [1:length(freqs)]
  dft_real(idx) = dft_real(idx) + (ex*cos(2 * pi * freqs(idx) * t));
  dft_imag(idx) = dft_imag(idx) + (ex*sin(2 * pi * freqs(idx) * t));  
end


[y, Fs] = audioread('fuelfandango.wav');
T = 1/Fs;
Y = fft(y);
magnitude = abs(Y);
[maxMagnitude, maxIndex] = max(magnitude);
dominantFrequency = maxIndex * Fs / length(y); % Calculate dominant frequency
disp(['Dominant Frequency causing the problem: ' num2str(dominantFrequency) ' Hz']);

% Filter to remove the tone (Page 545 of TextBook).
omega_to_remove = (2*pi*dominantFrequency)*T;
pole1 = exp(1j * omega_to_remove);
pole2 = exp(-1j * omega_to_remove);
filter_zeros_coeff = [pole1, pole2];
filter_zeros = poly(filter_zeros_coeff);
kp = 0.96;
filter_pole_coeff = [kp*pole1, kp*pole2];
filter_pole = poly(filter_pole_coeff);

figure;
zplane(filter_zeros, filter_pole);
title('Pole-Zero Plot of Notch Filter');
xlabel('Real');
ylabel('Imaginary');

% Initialize filtered signal
filtered_y = [];

% Implement the filter using a for loop for both channels
for ch = 1:size(y, 2) % Iterate over channels
    x = y(:, ch); % Extract the channel data
    M = max(length(filter_pole), length(filter_zeros)); % Determine filter order
    filtered_x = [];
    for n = 1:length(x)
        temp_filtered = 0;
        for k = 1:M
            if n - k + 1 > 0
                temp_filtered = temp_filtered + filter_zeros(k) * x(n - k + 1);
            end
        end
        for k = 2:M
            if n - k + 1 > 0
                temp_filtered = temp_filtered - filter_pole(k) * filtered_x(n - k + 1);
            end
        end
        filtered_x(n) = temp_filtered;
    end
    filtered_y = [filtered_y, filtered_x]; % Concatenate filtered signal
end

% Reshape filtered signal to match original dimensions
filtered_y = reshape(filtered_y, size(y));

% Play the filtered audio
sound(filtered_y, Fs);

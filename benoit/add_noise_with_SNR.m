function [ noisy_data , noise_var ] = add_noise_with_SNR( y_stim , SNR_lin )

noisy_data = zeros(size(y_stim));
noise_var = zeros(1, size(y_stim, 1));

variance = var(y_stim, [], 2);
variance_m = mean( variance );

for k = 1 : size( y_stim , 1 )  
    noise_var(k) = 1/SNR_lin * (0.5 * variance_m + rand * variance_m);               %Here make the assumption that the variance of the noise on each electrode 
                                                                                    % vary of 50% around the mean value
    noisy_data(k,:) = sqrt( noise_var(k) ) * randn(1,size(y_stim,2));
end

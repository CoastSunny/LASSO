function [ Times_courses , ROIs_phase, data, noisy_data ] = create_timecourses( Source_distribution , EEG_fwd , noise_level, phase )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to define a topographic activity from a given source distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:           - Source_distribution: the amplitudes of the sources in the cortical tesselation.
%                   - EEG_fwd: the EEG lead field.
%                   - noise_level: the noise level defined by the ration between the signal variance and the noise 
%                   variance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:           - Times courses: the ditribution on the electrodes over time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Benoit Cottereau
% Date: 16/03/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = 0 : pi / 45 : 2 * pi;
data = zeros( size(EEG_fwd,2) , length(x) );

ROIs_amp = unique( Source_distribution );
ROIs_amp = ROIs_amp( 2 : end );
ROIs_phase = [];
for ndx_ROIs = 1 : length( ROIs_amp )
    activated_sources = find( Source_distribution == ROIs_amp( ndx_ROIs ) );
    y = cos( 4 * x - phase * pi / 180 ) + ...                   % New signal with a random phase delay between 0? and 360?
        cos( 8 * x - phase * pi / 180 ) + ...
        cos( 12 * x - phase * pi / 180 ) + ...
        cos( 16 * x - phase * pi / 180 );
    for ndx = 1 : length( activated_sources )
        data( activated_sources(ndx) , : ) = y * Source_distribution( activated_sources(ndx) );
    end
    ROIs_phase( ndx_ROIs ) = phase;                        % Get the phase of the cluster
    
end


y_stim = EEG_fwd * data;                                        % Creation of the associated new EEG dataset
Times_courses = y_stim;
noisy_data = [];
if noise_level > 0
    [noisy_data,noise_var] = add_noise_with_SNR( y_stim , noise_level ); 
    Times_courses = Times_courses + noisy_data;                            % Add some noise
end

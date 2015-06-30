function [signal, Source_distribution, data, noise] = GenerateData(ROIs, ndx_areas, VertConn, EEG_fwd, noise_level, phase)

ratio_va_outva = [ 3 0 ];           % Ratio between the number of activated clusters within the predefined visual areas and outside (i.e. [3 2] means
                                    % that 3 visual areas will have activated clusters and two other cluster will be created outside the visual ROIs).
cluster_size = 30;                 % Size (in per cent) of the activated clusters. For clusters outside the visual areas the size is obtained from                                   % the average of the ROI areas
% noise_level = 10;                   % the noise level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Visual_areas = zeros(1,length(VertConn));               % BRC: This part isn't so important 
for ndx = 1 : length(ROIs.ndx)              % Define the sources within the visual areas
    Visual_areas( ROIs.ndx{ndx} ) = ndx;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creation of a source vector
if ratio_va_outva(1)
    [ Source_distribution, strength ] = create_source_distribution( ROIs , ndx_areas, cluster_size , ratio_va_outva(1) , VertConn );
end
if ratio_va_outva(2)
    [ Source_distribution_out_va ] = create_source_distribution_out_va( Visual_areas , cluster_size , ratio_va_outva(2) , VertConn );
end
%Source_distribution = Source_distribution + Source_distribution_out_va;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creation of a source vector
[ signal, ~, data, noise ] = create_timecourses( Source_distribution , EEG_fwd , noise_level, phase );
%%%%%%%%%%%%%%%%%%
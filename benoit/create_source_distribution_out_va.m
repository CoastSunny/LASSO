function [ Sources_distribution ] = create_source_distribution_out_va( Unavailable_sources , cluster_size , cluster_nb , VertConn )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to define random activation in 'cluster_nb' areas of the visual system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:           - ROIs: cell containing the definition of the visual ROIs
%                   - cluster_size: size of the activetad area within each
%                   ROI (defined in percent).
%                   - cluster_nb: number of activated area.
%                   - VertConn: connectivity of the cortical tesselation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:           - Sources_distribution: the source ditributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Benoit Cottereau
% Date: 10/07/2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_size_ROI = 60;         % Define the mean number of sources within the visual areas (allows to have on average the same number of sources for the ROIs outside the visual areas)
Sources_distribution = zeros( 1 , length(VertConn) );

% Loop on all the activated clusters
for k  = 1 : cluster_nb

    ROI = setdiff( [1 : length(VertConn)] , Unavailable_sources);
    tmp = randperm(length(ROI));

    growing_seed = zeros(1,length(VertConn));
    growing_seed( ROI(tmp(1)) ) = 1;
    Sources_distribution_tmp = [];
    while (length( Sources_distribution_tmp ) < cluster_size / 100 * mean_size_ROI )
        growing_seed = dilatation( growing_seed , VertConn , 1 );
        Sources_distribution_tmp = intersect( ROI , find(growing_seed) );
    end
    if ( length( Sources_distribution_tmp ) > round( cluster_size / 100 * mean_size_ROI ) )
        Sources_distribution_tmp = Sources_distribution_tmp( 1 : round( cluster_size / 100 * mean_size_ROI ) );
    end
    
    % Random choice of the source amplitude within the ROI
    tmp = randperm(40);             
    ROI_amp = 1 + tmp(1) / 10;
    Sources_distribution( Sources_distribution_tmp ) = ROI_amp;
    Unavailable_sources = [ Unavailable_sources , find(Sources_distribution) ];
    
end
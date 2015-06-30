function [ Sources_distribution, strength ] = create_source_distribution( ROIs , ndx_areas, cluster_size , cluster_nb , VertConn )
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
tmp = randperm( length(ROIs.ndx) );
%ndx_areas = tmp( 1 : cluster_nb );
%ndx_areas = [ 3 5  ];

fixed_cluster_size = cluster_size;

Sources_distribution = zeros( 1 , length(VertConn) );
strength = zeros(1,numel(ndx_areas));
% Loop on all the activated clusters
for k  = 1 : length(ndx_areas)

    ROI = ROIs.ndx{ ndx_areas(k) };
%     if ((ROI == 15) | (ROI==16))
%         cluster_size = 10;
%     else
%         cluster_size = fixed_cluster_size;
%     end
    tmp = randperm(length(ROI));
    growing_seed = zeros(1,length(VertConn));
    growing_seed( ROI(tmp(1)) ) = 1;
    Sources_distribution_tmp = [];
    i = 1;
    while (length( Sources_distribution_tmp ) < cluster_size / 100 * length(ROI) )
        growing_seed = dilatation( growing_seed , VertConn , 1 );
        Sources_distribution_tmp = intersect( ROI , find(growing_seed) );
        i = i + 1;
        if i > 50
            i;
            break
        end
    end
    if ( length( Sources_distribution_tmp ) > round( cluster_size / 100 * length(ROI) ) )
        Sources_distribution_tmp = Sources_distribution_tmp( 1 : round( cluster_size / 100 * length(ROI) ) );
    end
    
    % Random choice of the source amplitude within the ROI
%     tmp = randperm(40);             
%     ROI_amp = 1 + tmp(1) / 10;

    tmp = randi(90);             
    ROI_amp = 1 + tmp / 10;
    strength(k) = ROI_amp;
%     if k > 2
%         ROI_amp = 2;
%     end
    Sources_distribution( Sources_distribution_tmp ) = ROI_amp;
    
end

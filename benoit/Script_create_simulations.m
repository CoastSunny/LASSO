% load Subject_48_initialization
% 
% rand('state',0);
% 
% idx = 1:14;
% ROIs.Tess_ndx = ROIs.Tess_ndx(idx);
% ROIs.corr = ROIs.corr(idx);
% ROIs.name = ROIs.name(idx);
% ROIs.ndx = ROIs.ndx(idx);
% 
idx = [ROIs.ndx{1:14}];
% 
% ndx_areas = [1 3];
% noise_level = 10;
% 
% [Y, signal] = GenerateData(ROIs, ndx_areas, VertConn, EEG_fwd, noise_level);
% Y = Y(:,1);

% Then, you can compute the minimum-norm (L2 solution) using: 

J_MN = minimum_norm( EEG_fwd(:,idx), Y,NUM_OF_LAMBDA );

% and the minimum-norm with FACE using the ROIs structure and: 
% 
% G = EEG_fwd;
% for ndx_area = 1 : length( ROIs.name )
%    G( : , ROIs.ndx{ ndx_area } ) =  G( : , ROIs.ndx{ ndx_area } ) * chol( ROIs.corr{ndx_area} )'; 
% end 
% [ J_FACE ] = minimum_norm( G , Topo_distribution ); 
% for ndx_area = 1 : length( ROIs.name )
%    J_FACE( ROIs.ndx{ ndx_area } , : ) =  chol( ROIs.corr{ndx_area} )' * J_FACE( ROIs.ndx{ ndx_area } , : ); 
% end
% Estimation of the results
[ AUC , AUC_close , AUC_far , n_MSE , n_DF , relative_energy ] = Inverse_perf_estimation( signal(idx) , J_MN , VertConn );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

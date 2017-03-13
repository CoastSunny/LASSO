function [ AUC , AUC_close , AUC_far , mse , n_DF , relative_energy ] = Inverse_perf_estimation( sourceDistribution , J , VertConn )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to define the AUC and relative amplitude error of an inverse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:           - Source_distribution: the real source distribution
%                   - J: the estimation of the source distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:           - AUC: the Area Under the ROC Curve
%                   - relative_amp_error: the error in the estimation of
%                   the relative amplitudes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Benoit Cottereau
% Date: 10/07/2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ac_sources = find( sourceDistribution );               % ndx to the activated sources

tmp = zeros( 1 , length(sourceDistribution) );
tmp( ac_sources ) = 1; 
AUC = rocArea( abs(J)' , tmp );

growing_seed = tmp;
while ( length( setdiff( find(growing_seed) ,ac_sources ) )  < length( ac_sources ) )
    growing_seed = dilatation( growing_seed , VertConn , 1 );
end
AUC_close = rocArea( abs( J( find(growing_seed) ) ) , tmp( find(growing_seed) ) );
[ a , far_sources_rnk ] = sort( abs( J( find( growing_seed == 0 ) ) ) );
ndx_sources = far_sources_rnk( end - length( find(tmp) ) : end );
AUC_far = rocArea( abs( J( [ find(tmp) , ndx_sources' ] ) ) , tmp( ( [ find(tmp) , ndx_sources' ] ) ) );
mse = norm(sourceDistribution-J, 'fro')^2 / numel(J);    % Mean square error
% DF = sum( ( Source_distribution(ac_sources) - J(ac_sources)' ).^2 ) / sum( (Source_distribution(ac_sources) ).^2 );  % Degree of Focalization
% 
n_Source_distribution = sourceDistribution / max(abs(sourceDistribution));
n_J = J / max( abs(J) );
%n_MSE = sum( (n_Source_distribution - n_J).^2 ) / sum( (n_Source_distribution).^2 );    % Normalized mean square error
n_DF = sum( (n_Source_distribution(ac_sources) - n_J(ac_sources)).^2 ) / sum( (n_Source_distribution(ac_sources)).^2 );    % Normalized mean square error

relative_energy = sum( abs( n_J(ac_sources) ) ) / sum( abs(n_J) );


    

function [hh,fg,cb] = plot_vector_field2( ph, varargin )
% *WAVE*
%
% PLOT VECTOR FIELD
%
% INPUT
% ph - complex-valued scalar field representing vector directions
%
% OUTPUT
% vector field plot
%

% checks
assert( ~isreal(ph), 'complex-valued input required, ph' )

% set defaults
if nargin > 1, plot_option = varargin{1}; else plot_option = 0; end

% init
[XX,YY] = meshgrid( 1:(size(ph,2)/2), 1:(size(ph,1)/2) );

Y = conv2(ph,[1,1;1,1],'valid');
Z = Y(1:2:end,1:2:end)/4

M = real( exp( 1i * angle(Z) ) ); 
N = imag( exp( 1i * angle(Z) ) );

% plotting
% fg = figure;
if ( plot_option == 0 )
    imagesc( angle(ph) ); axis image;%cb = colorbar;  %caxis( [-pi pi] ); 
    hold on;
    %set( get(cb,'ylabel'), 'string', 'Direction (rad)' )
    quiver( XX, YY, M, N, 0.5, 'r' );
    set( gca, 'fontname', 'arial', 'fontsize', 14, 'ydir', 'normal' ); hh = gca;
elseif ( plot_option == 1 )
%     imagesc( angle(ph) ); axis image; cb = colorbar; caxis( [-pi pi] ); 
%     hold on;
%     set( get(cb,'ylabel'), 'string', 'Direction (rad)' )
    quiver( 2*XX-0.5, 2*YY-0.5, M, N, 0.5, 'k', 'linewidth', 0.5);% 'MaxHeadSize', 0.2 );
%     set( gca, 'fontname', 'arial', 'fontsize', 14, 'ydir', 'normal' ); hh = gca;
    % delete( ih );  %axis off 
    %delete( cb );
end

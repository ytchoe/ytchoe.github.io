function [hh,fg,cb] = plot_vector_field_14May2021( ph, varargin )
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
[XX,YY] = meshgrid( 1:size(ph,2), 1:size(ph,1) );
M = real( exp( 1i * angle(ph) ) ); N = imag( exp( 1i * angle(ph) ) );

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
% 
%     % Stim Trial 4
%     [sx,sy] = meshgrid(10:18, 4:18);
% 	hlines_sens = streamline(stream2(XX,YY,M,N,sx, sy));
%     set(hlines_sens,'LineWidth',2,'Color','r')
%     hold on;

   % Baseline RL corner
%     [sx,sy] = meshgrid(18:32, 1:16);
        [sx,sy] = meshgrid(1:32, 1:32);

	hlines_sens = streamline(stream2(XX,YY,M,N,sx, sy));
    set(hlines_sens,'LineWidth',2,'Color','y')
    hold on;

    quiver( XX, YY, M, N, 0.75, 'k',  'linewidth', 1);% ,'MaxHeadSize', 2 );
    axis equal
%     set( gca, 'fontname', 'arial', 'fontsize', 14, 'ydir', 'normal' ); hh = gca;
    % delete( ih );  %axis off 
    %delete( cb );
end

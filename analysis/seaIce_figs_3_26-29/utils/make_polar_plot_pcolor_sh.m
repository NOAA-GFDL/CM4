function h=make_polar_plot_pcolor_sh(x,y,ifXY,data,radius,title_str,iscbar,cmin,cmax)

%x and y are square matrices with lat-lon coordinates
%data is a vector of length nD
%ifXY is the mask that defines the nD gridpoints
%radius is the radius of the polar stereographic plot

addpath('./m_map/')
load('cmap_jet3.mat','cmap')

	A=zeros(size(ifXY));
	A(ifXY)=data;
	A(~ifXY)=NaN;

x=[x;x(1,:)];
y=[y;y(1,:)];
A=[A;A(1,:)];


                m_proj('stereographic','lat',-90,'long',0,'radius',radius);
        	hold on
                h = m_pcolor( x, y, A );

                %shading interp
                shading flat
		colormap(cmap)
                caxis([cmin cmax])
                m_grid('xtick',[],'tickdir','out','ytick',[],'linest','-');
                m_coast('patch',rgb('Gainsboro'));
		%m_coast( 'line', 'linewidth', 1, 'color', 'k' );
                set(gcf,'color', 'w')
		title(title_str)
		if(iscbar)
		colorbar
		end
		drawnow





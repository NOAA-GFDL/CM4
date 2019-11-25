function h=make_polar_contour_plot(x,y,ifXY,data,radius,title_str,edge_color,line_width)

%x and y are square matrices with lat-lon coordinates
%data is a vector of length nD
%ifXY is the mask that defines the nD gridpoints
%radius is the radius of the polar stereographic plot

addpath('../m_map/')
addpath('./m_map/')
load('cmap_jet3.mat','cmap')

load('/work/Mitchell.Bushuk/raw_data/CM2.1R/ocean_grid.mat','x_full','y_full')

	A=zeros(size(ifXY));
	A(ifXY)=data;
	A(~ifXY)=NaN;

                m_proj('stereographic','lat',90,'long',0,'radius',radius);
        	hold on
		h = m_contour( x, y, A ,[0.15],edge_color,'LineWidth',line_width);
                %shading interp
                m_grid('xtick',[],'tickdir','out','ytick',[],'linest','-');
                m_coast('patch',rgb('Gainsboro'));
                set(gcf,'color', 'w')
		title(title_str)
		drawnow






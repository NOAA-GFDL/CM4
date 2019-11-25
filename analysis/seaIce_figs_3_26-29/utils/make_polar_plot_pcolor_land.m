function h=make_polar_plot_pcolor(x,y,ifXY,data,radius,title_str,iscbar,cmin,cmax)

%x and y are square matrices with lat-lon coordinates
%data is a vector of length nD
%ifXY is the mask that defines the nD gridpoints
%radius is the radius of the polar stereographic plot

addpath('../m_map/')
addpath('./m_map/')
%load('cmap_jet.mat','cmap')
%load('cmap_jet2.mat','cmap')
load('cmap_jet3.mat','cmap')
%load('cmap_binary.mat','cmap')

	A=zeros(size(ifXY));
	A(ifXY)=data;
	A(~ifXY)=NaN;
	%A(~ifXY)=-10;

A=[A;A(1,:)];
x=[x;x(1,:)];
y=[y;y(1,:)];

%A=[A,A(:,1)];
%x=[x,x(:,1)];
%y=[y,y(:,1)];

                m_proj('stereographic','lat',90,'long',0,'radius',radius);
        	hold on
                h = m_pcolor( x, y, A );
%                h = m_pcolor( x_ice, y_ice, A );

                
%		h = m_pcolor( x, y, A );

                %shading interp
                shading flat
		colormap(cmap)
                caxis([cmin cmax])
                m_grid('xtick',[],'tickdir','out','ytick',[],'linest','-');
                %m_coast('patch',rgb('Gainsboro'));
		m_coast( 'line', 'linewidth', 1, 'color', 'k' );
                set(gcf,'color', 'w')
		title(title_str)
		if(iscbar)
		 pos=get(gca,'position');

                        c=colorbar;

                 set(gca,'position',pos)

                 x_pos=get(c, 'Position');
                 y_pos=[x_pos(1)+0.1 x_pos(2) x_pos(3)/2 x_pos(4)];
                 %y_pos=[x_pos(1) x_pos(2)-0.08 6*x_pos(3) x_pos(4)];
                 set(c,'Position',y_pos)
		end
		drawnow


%F = scatteredInterpolant(x(ifXY),y(ifXY),data,'linear','nearest');
%F = scatteredInterpolant(x(ifXY),y(ifXY),data,'linear','linear');
%F = scatteredInterpolant(x(ifXY),y(ifXY),data,'nearest','nearest');
%xq=x_c(ifXY_c);
%yq=y_c(ifXY_c);
%vq=F(xq,yq);

%if(fillpole==1)
%
%	B=zeros(size(x_c));
%	%B(ifXY_c)=vq;
%	B(ifXY)=data;
%	B(~ifXY_c)=NaN;
%
%end
%
%		if(fillpole==1)
%	        h = m_pcolor( x_c, y_c, B );
%		end
%






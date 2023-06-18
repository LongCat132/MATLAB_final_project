function pixs=deg2pix(degree,inch,pwidth,vdist,ratio) 
% parameters: degree, inch (monitor size), pwidth (width in pixels), 
% vdist: viewsing distance
% ratio: ration = pheight/pwidth 高宽比
screenWidth = inch*2.54/sqrt(1+ratio^2);  
pix=screenWidth/pwidth; 
pixs = round(2*tan((degree/2)*pi/180) * vdist / pix); 
end
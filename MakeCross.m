function cross=MakeCross(wPtr,crosswidth,crosslength,pxlpdeg,RGB,bkRGB)

widthpxl=round(crosswidth*pxlpdeg);
lengthpxl=round(crosslength*pxlpdeg);
crossrect=zeros(lengthpxl,lengthpxl,3);
crossrect(:,:,1)=bkRGB(1);
crossrect(:,:,2)=bkRGB(2);
crossrect(:,:,3)=bkRGB(3);

vlinex1=round(lengthpxl/2-widthpxl/2);
vlinex2=round(lengthpxl/2+widthpxl/2);
vliney1=1;
vliney2=lengthpxl;
hlinex1=1;
hlinex2=lengthpxl;
hliney1=round(lengthpxl/2-widthpxl/2);
hliney2=round(lengthpxl/2+widthpxl/2);

crossrect(vliney1:vliney2,vlinex1:vlinex2,1)=RGB(1);
crossrect(vliney1:vliney2,vlinex1:vlinex2,2)=RGB(2);
crossrect(vliney1:vliney2,vlinex1:vlinex2,3)=RGB(3);
crossrect(hliney1:hliney2,hlinex1:hlinex2,1)=RGB(1);
crossrect(hliney1:hliney2,hlinex1:hlinex2,2)=RGB(2);
crossrect(hliney1:hliney2,hlinex1:hlinex2,3)=RGB(3);

cross=Screen('MakeTexture',wPtr,crossrect);





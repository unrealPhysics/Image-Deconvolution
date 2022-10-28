function cleanImage()

% SPLIT UNCLEAN DOT IMAGE INTO MANY SMALL PANELS OF RANDOM DOT IMAGES
% USE HIGH QUALITY IMAGE OF SAME SIZE AS SMALL PANELS (possibly multiple)
% TAKE MULTIPLE PHOTOS AT EVERY FOCUS LEVEL
    % USE SMALL CHANGES IN FOCUS LEVEL
% USE A GOOD CAMERA

cleanImage='page3.png';
blurImage='blurpage3.png';

cleanData=imread(cleanImage);
blurData=imread(blurImage);

% make greyscale and fix size
for colourIndex=1:3

% heck with stuff to make cutpoints same across images
[yloc1,yloc2,dy,xloc1,xloc2,dx]=generateCutPoints(cleanData);
[cg1,cg2,cg3,cg4,cg5,cg6,cg7,cg8]=colourCut(cleanData,yloc1,yloc2,dy, ...
    xloc1,xloc2,dx,colourIndex);
[bg1,bg2,bg3,bg4,bg5,bg6,bg7,bg8]=colourCut(blurData,yloc1,yloc2,dy, ...
    xloc1,xloc2,dx,colourIndex);

% take fourier transform

x1Spec=fft2(cg5);
x2Spec=fft2(cg6);
x3Spec=fft2(cg7);
x4Spec=fft2(cg8);

z1Spec=fft2(bg5);
z2Spec=fft2(bg6);
z3Spec=fft2(bg7);
z4Spec=fft2(bg8);

ZXconj=(z1Spec.*conj(x1Spec)+z2Spec.*conj(x2Spec) ...
    +z3Spec.*conj(x3Spec)+z4Spec.*conj(x4Spec))/4;
XXconj=(x1Spec.*conj(x1Spec)+x2Spec.*conj(x2Spec) ...
    +x3Spec.*conj(x3Spec)+x4Spec.*conj(x4Spec))/4;

[yLen,xLen]=size(ZXconj);

for x=1:xLen
    for y=1:yLen
        untransfer(y,x)=XXconj(y,x)/ZXconj(y,x);
    end
end

rg1=ifft2(fft2(bg1).*untransfer);
rg2=ifft2(fft2(bg2).*untransfer);
rg3=ifft2(fft2(bg3).*untransfer);
rg4=ifft2(fft2(bg4).*untransfer);

if colourIndex~=4
%     imageRecon(:,:,colourIndex)=rg1;
    imageRecon(:,:,colourIndex)=(rg1+rg2+rg3+rg4)/4;
%     imageClean(:,:,colourIndex)=cg1;
    imageClean(:,:,colourIndex)=(cg1+cg2+cg3+cg4)/4;
%     imageDirty(:,:,colourIndex)=bg1;
    imageDirty(:,:,colourIndex)=(bg1+bg2+bg3+bg4)/4;
    imageKazarel(:,:,colourIndex)=[cg1,cg2;cg3,cg4;cg5,cg6;cg7,cg8];
end
% rg2=ifft2(fft2(bg2).*untransfer);
% rg3=ifft2(fft2(bg3).*untransfer);
% rg4=ifft2(fft2(bg4).*untransfer);

% do each blur level on its own to see loss ratio
    % to calculate loss ratio take the average chromatic difference of all
    % pixels between the original and recovered ones
        % might not be a fantastic metric of comparison but it exists

% take the fourier transform (perhaps inverse, check with matlab functions)
% of each (line-by-line? globular? whole image?)

% calculate the modular transfer function (and thus blurring function)

% multiply by the transfer function (or inverse (either inverses)) take the
% inverse transform of this and thus obtain the image

% calculate loss ratio
end

figure(5)
subplot(2,2,1)
image(imageRecon)
subplot(2,2,2)
image(imageClean)
subplot(2,2,3)
image(imageDirty)

figure(1)
imagesc(imageKazarel)


end

% is functionally fixed at the specific style I chose to use
function [yloc1,yloc2,dy,xloc1,xloc2,dx]=generateCutPoints(data)

[yLen,xLen,cLen]=size(data);

greyData=sum(double(data),3); %flip maybe, idk

for y=2:(yLen-1)
    yDeviant(y) = abs(2*greyData(y,10) - greyData(y-1,10) - greyData(y+1,10) ...
        + 2*greyData(y,xLen-10) - greyData(y-1,xLen-10) - greyData(y+1,xLen-10));
end

[m,yloc1] = max(yDeviant(1:floor(length(yDeviant)/2)));
[m,yloc2] = max(yDeviant(floor(length(yDeviant)/2):end));
yloc1=yloc1+3;
yloc2=yloc2+floor(length(yDeviant)/2)-2;
yCanvas=yloc2-yloc1;
dy=floor(yCanvas/4); % hard coded and affixed

for x=2:(xLen-1)
    xDeviant(x) = abs(2*greyData(10,x) - greyData(10,x-1) - greyData(10,x+1) ...
        + 2*greyData(yLen-10,x) - greyData(yLen-10,x-1) - greyData(yLen-10,x+1));
end

[m,xloc1] = max(xDeviant(1:floor(length(xDeviant)/2)));
[m,xloc2] = max(xDeviant(floor(length(xDeviant)/2):end));
xloc1=xloc1+3;
xloc2=xloc2+floor(length(xDeviant)/2)-2;
xCanvas=xloc2-xloc1;
dx=floor(xCanvas/2); % hard coded and affixed

end

function [i1,i2,i3,i4,i5,i6,i7,i8]=colourCut(data,yloc1,yloc2,dy,xloc1,xloc2,dx,colourIndex)

greyData=double(data(:,:,colourIndex))/mean(data(:,:,colourIndex),'all');

i1=greyData(((yloc1):(yloc1+dy)),((xloc1):(xloc1+dx)));
i2=greyData(((yloc1):(yloc1+dy)),((xloc2-dx):(xloc2)));

i3=greyData(((yloc1+dy):(yloc1+2*dy)),((xloc1):(xloc1+dx)));
i4=greyData(((yloc1+dy):(yloc1+2*dy)),((xloc2-dx):(xloc2)));

i5=greyData(((yloc2-2*dy):(yloc2-dy)),((xloc1):(xloc1+dx)));
i6=greyData(((yloc2-2*dy):(yloc2-dy)),((xloc2-dx):(xloc2)));

i7=greyData(((yloc2-dy):(yloc2)),((xloc1):(xloc1+dx)));
i8=greyData(((yloc2-dy):(yloc2)),((xloc2-dx):(xloc2)));

end
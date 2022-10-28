function cleanImageNew()

% SPLIT UNCLEAN DOT IMAGE INTO MANY SMALL PANELS OF RANDOM DOT IMAGES
% USE HIGH QUALITY IMAGE OF SAME SIZE AS SMALL PANELS (possibly multiple)
% TAKE MULTIPLE PHOTOS AT EVERY FOCUS LEVEL
    % USE SMALL CHANGES IN FOCUS LEVEL
% USE A GOOD CAMERA

% clnImg=["./ph/DSC_0388.JPG", "./ph/DSC_0389.JPG"];
% blrImg=["./ph/resized408.JPG", "./ph/resized409.JPG"];

xSplits=1;
ySplits=2;
splitBoxes=xSplits*ySplits;

cleanData(:,:,:,1)=imread("./ph/DSC_0386","jpg");
cleanData(:,:,:,2)=imread("./ph/DSC_0387","jpg");
cleanData(:,:,:,3)=imread("./ph/DSC_0388","jpg");
cleanData(:,:,:,4)=imread("./ph/DSC_0389","jpg");
[~,~,~,cleanCopies]=size(cleanData);

% blurData(:,:,:,1)=imread("./ph/resized418","jpg");
% blurData(:,:,:,2)=imread("./ph/resized418b","jpg");
% blurData(:,:,:,3)=imread("./ph/resized418c","jpg");
% blurData(:,:,:,4)=imread("./ph/resized418d","jpg");

blurData(:,:,:,1)=imread("./ph/resized402","jpg");
blurData(:,:,:,2)=imread("./ph/resized403","jpg");
blurData(:,:,:,3)=imread("./ph/resized404","jpg");

% blurData(:,:,:,1)=imread("./ph/resized397","jpg");
% blurData(:,:,:,2)=imread("./ph/resized398","jpg");
% blurData(:,:,:,3)=imread("./ph/resized399","jpg");
[~,~,~,dirtyCopies]=size(blurData);

totalCellCount=splitBoxes*lcm(cleanCopies,dirtyCopies);

[dotRectangle,sceneRectangle]=generateCutPoints(cleanData(:,:,:,1));

% make monochrome and fix size
for colourIndex=1:3

% heck with stuff to make cutpoints same across images (done externally)
for k=1:cleanCopies
    [cd(:,:,(((k-1)*splitBoxes+1):(k*splitBoxes))), ...
     cs(:,:,(((k-1)*splitBoxes+1):(k*splitBoxes)))]=...
    colourCut(cleanData(:,:,:,k),dotRectangle,sceneRectangle,colourIndex, ...
        xSplits,ySplits);
end
for k=1:dirtyCopies
    [bd(:,:,(((k-1)*splitBoxes+1):(k*splitBoxes))), ...
     bs(:,:,(((k-1)*splitBoxes+1):(k*splitBoxes)))]=...
    colourCut(blurData(:,:,:,k),dotRectangle,sceneRectangle,colourIndex, ...
        xSplits,ySplits);
end

% take fourier transform

ZXconj=0;
XXconj=0;
% [yLen,xLen,~]=size(cd);
for k=1:(splitBoxes*cleanCopies)
    zSpec(:,:,k)=fft2(cd(:,:,k));
%     zSpec(:,:,k)=fft(fft(cd(:,:,k).*hanning(yLen)).'.*hanning(xLen)).';
end
for k=1:(splitBoxes*dirtyCopies)
    xSpec(:,:,k)=fft2(bd(:,:,k));
%     xSpec(:,:,k)=fft(fft(bd(:,:,k).*hanning(yLen)).'.*hanning(xLen)).';
end

for k=1:totalCellCount
    xk=k;
    zk=k;
    while xk>splitBoxes*dirtyCopies
        xk=xk-splitBoxes*dirtyCopies;
    end
    while zk>splitBoxes*cleanCopies
        zk=zk-splitBoxes*cleanCopies;
    end
    XXconj=XXconj+xSpec(:,:,xk).*conj(xSpec(:,:,xk));
    ZXconj=ZXconj+zSpec(:,:,zk).*conj(xSpec(:,:,xk));
end

untransfer=ZXconj./XXconj;

for k=1:splitBoxes
    pseudoSpec(:,:,k)=fft2(bs(:,:,k)).*abs(untransfer);
    pseudoSpek(:,:,k)=fft2(bs(:,:,k)).*untransfer;
    rs(:,:,k)=ifft2(pseudoSpec(:,:,k));
    rk(:,:,k)=ifft2(pseudoSpek(:,:,k));
%     ps(:,:,k)=ifft2(fft2(cs(:,:,k)));
end

if colourIndex~=4
    [yLen,xLen]=size(rs(:,:,1));
    imageRecon(:,:,colourIndex)=zeros(ySplits*yLen,xSplits*xLen);
    imageRekon(:,:,colourIndex)=zeros(ySplits*yLen,xSplits*xLen);

    k=0;
    for n=1:xSplits
        for m=1:ySplits
            k=k+1;
            imageRecon(((m-1)*yLen+1):(m*yLen), ...
                       ((n-1)*xLen+1):(n*xLen), ...
                       colourIndex)=rs(:,:,k);
            imageRekon(((m-1)*yLen+1):(m*yLen), ...
                       ((n-1)*xLen+1):(n*xLen), ...
                       colourIndex)=rk(:,:,k);
            imageClean(((m-1)*yLen+1):(m*yLen), ...
                       ((n-1)*xLen+1):(n*xLen), ...
                       colourIndex)=cs(:,:,k);
%             imagePseudo(((m-1)*yLen+1):(m*yLen), ...
%                         ((n-1)*xLen+1):(n*xLen), ...
%                         colourIndex)=ps(:,:,k);
            imageDirty(((m-1)*yLen+1):(m*yLen), ...
                       ((n-1)*xLen+1):(n*xLen), ...
                       colourIndex)=bs(:,:,k);
        end
    end
else
    [yLen,xLen]=size(rs(:,:,1));
    imageRetcon(:,:,3)=zeros(ySplits*yLen,xSplits*xLen);

    for falseColour=1:3
    k=0;
        for n=1:xSplits
            for m=1:ySplits
                k=k+1;
                imageRetcon(((m-1)*yLen+1):(m*yLen), ...
                            ((n-1)*xLen+1):(n*xLen), ...
                             falseColour)=rs(:,:,k);
            end
        end
    end
    imageRetcon=imageRetcon/3;
end

end

figure(2)
clf
image(imageRecon)
title("Reconstructed Image Ignoring Phase Differences")
figure(3)
clf
image(imageClean)
title("Original Clean Image")
figure(4)
clf
image(imageDirty)
title("Original Blurred Image")

% imageGreyset(:,:,1)=(imageRecon(:,:,1)+imageRecon(:,:,2)+imageRecon(:,:,3))/3;
% imageGreyset(:,:,2)=(imageRecon(:,:,1)+imageRecon(:,:,2)+imageRecon(:,:,3))/3;
% imageGreyset(:,:,3)=(imageRecon(:,:,1)+imageRecon(:,:,2)+imageRecon(:,:,3))/3;
% 
% figure(5)
% clf
% image(imageGreyset)

% figure(6)
% clf
% image(imageRetcon)

% figure(7)
% clf
% image(imageRetcon-imageGreyset)

figure(8)
clf
image(imageRekon)
title("Reconstructed Image Including Phase Differences")

% figure(9)
% clf
% image(imagePseudo)

chromaticFixNoPhase(:)=chromaticDifference(imageRecon,imageClean);
chromaticFixPhase(:)=chromaticDifference(imageRekon,imageClean);
chromaticLocks=["red  ", "green", "blue "];

for chromaticLoop=1:3
    fprintf("Difference in %s channel: %+.6f (no phase)\n", ...
        chromaticLocks(chromaticLoop),chromaticFixNoPhase(chromaticLoop))
    fprintf("Difference in %s channel: %+.6f (w/ phase)\n", ...
        chromaticLocks(chromaticLoop),chromaticFixPhase(chromaticLoop))
end

imwrite(imageClean,"./exports/clnScene3.jpg",'jpg');
imwrite(imageDirty,"./exports/blrScene3.jpg",'jpg');
imwrite(imageRekon,"./exports/otfScene3.jpg",'jpg');
imwrite(imageRecon,"./exports/mtfScene3.jpg",'jpg');

end

% is functionally fixed at the specific style I chose to use
function [dotRectangle,sceneRectangle] = generateCutPoints(dataClean)

cleanDouble=(double(dataClean))/255; %flip maybe, idk

confirmed = false;

while ~confirmed
    figure(1)
    imagesc(cleanDouble)
    title("Please select the top-left and bottom-right corners of the dot pattern")
    [x,y] = ginput(2);
    canvasWidth=round(x(2)-x(1));
    canvasHeight=round(y(2)-y(1));
    dotRectangle = [round(x(1)), round(y(1)), canvasWidth, canvasHeight];
    figure(1)
    hold on
    rectangle('Position',dotRectangle)
    answer = inputdlg("Is this rectangle correct? (y/n)", "Confirm rectangle",1,"y");
    if answer == "y"
        confirmed = true;
    end
end

confirmed = false;

while ~confirmed
    figure(1)
    imagesc(cleanDouble)
    title("Please select the top-left corner of the scene")
    [xS,yS] = ginput(1);
    sceneRectangle = [round(xS), round(yS), canvasWidth, canvasHeight];
    figure(1)
    hold on
    rectangle('Position',sceneRectangle)
    answer = inputdlg("Is this rectangle correct? (y/n)", "Confirm rectangle",1,"y");
    if answer == "y"
        confirmed = true;
    end
end

end

function [d,s]=colourCut(data,dotRectangle,scnRectangle,colourIndex,xSplits,ySplits)

if colourIndex==4
    greyData=sum(double(data),3);
else
    greyData=double(data(:,:,colourIndex))/255;
end

dx1=floor(dotRectangle(3)/xSplits);
dy1=floor(dotRectangle(4)/ySplits);
dx2=ceil(dotRectangle(3)/xSplits);
dy2=ceil(dotRectangle(4)/ySplits);
greyData(:,xSplits*dx1:xSplits*dx2,colourIndex)=0;
greyData(ySplits*dy1:ySplits*dy2,:,colourIndex)=0;
dx=ceil(dotRectangle(3)/xSplits);
dy=ceil(dotRectangle(4)/ySplits);

dotX=dotRectangle(1);
dotY=dotRectangle(2);
scnX=scnRectangle(1);
scnY=scnRectangle(2);

for n=0:xSplits
    dxBinds(n+1)=n*dx;
end

for m=0:ySplits
    dyBinds(m+1)=m*dy; %#ok<*AGROW> 
end

k=0;
for n=1:xSplits
    for m=1:ySplits
        k=k+1;
        d(:,:,k)=greyData((dotY+dyBinds(m)):(dotY+dyBinds(m+1))...
                         ,(dotX+dxBinds(n)):(dotX+dxBinds(n+1)));
        s(:,:,k)=greyData((scnY+dyBinds(m)):(scnY+dyBinds(m+1))...
                         ,(scnX+dxBinds(n)):(scnX+dxBinds(n+1)));
    end
end

%should actually do nothing
d=fliplr(d); % may cause problems, not easy to check w/out causing them

end

function averageChromaticDifference = chromaticDifference(image1,image2)

[xLim,yLim,~] = size(image1);

chromaticDifference(:,:,1)=image1(:,:,1)-image2(:,:,1);
chromaticDifference(:,:,2)=image1(:,:,2)-image2(:,:,2);
chromaticDifference(:,:,3)=image1(:,:,3)-image2(:,:,3);

% figure("Name","chromaticDifference")
% clf
% image(chromaticDifference)

averageChromaticDifference=sum(chromaticDifference,[1 2])/(xLim*yLim);

end

% NIST-developed software is provided by NIST as a public service. You may use,
% copy and distribute copies of the software in any medium, provided that you keep
% intact this entire notice. You may improve, modify and create derivative works of
%the software or any portion of the software, and you may copy and distribute such
% modifications or works. Modified works should carry a notice stating that you
% changed the software and should note the date and nature of any such change.
% Please explicitly acknowledge the National Institute of Standards and Technology
% as the source of the software.

% NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY OF
% ANY KIND, EXPRESS, IMPLIED, IN FACT OR ARISING BY OPERATION OF LAW, INCLUDING, WITHOUT
% LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
% NON-INFRINGEMENT AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE
% OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS
% WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE
% USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS,
% ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.

% You are solely responsible for determining the appropriateness of using and distributing
% the software and you assume all risks associated with its use, including but not limited
% to the risks and costs of program errors, compliance with applicable laws, damage to or
% loss of data, programs or equipment, and the unavailability or interruption of operation.
% This software is not intended to be used in any situation where a failure could cause risk
% of injury or damage to property. The software developed by NIST employees is not subject
% to copyright protection within the United States.

clear; clc
croppedpath = 'C:\Users\ona1\Desktop\Manual_Masked\';
maskedpath = 'C:\Users\ona1\Desktop\Manual_Masked\';

for n = 1:182
    filename = ['colony_',sprintf('%03d',n),['.tif']];

    I = imread([croppedpath, filename], 'tiff');
    I2 = imopen(I,strel('disk',15));
    mylevel = graythresh(I-I2);
    I3 = im2bw(I-I2, mylevel);
    bw = bwareaopen(I3, 2);
    cc = bwconncomp(bw, 8);
    colonyedges = regionprops(cc, 'all');


    I = imread([maskedpath, filename], 'tiff');
    cc2 = bwconncomp(I, 8);
    labeled = labelmatrix(cc2);
    colonymask = regionprops(cc2, 'all');

    myresult1(n) = structmax(colonyedges, 'Area');
    
    myresult2(n) = structmax(colonymask, 'Area');
end




for n = 1:182
    myfeatures(n).AreaRatio = myresult1(n).Area/myresult2(n).Area;
    myfeatures(n).CentroidsDistance = pdist([myresult1(n).Centroid;myresult2(n).Centroid]);
    myfeatures(n).EdgeFactor = myresult1(n).FilledArea/myresult2(n).FilledArea;
    myfeatures(n).HolesPerArea = abs(myresult1(n).EulerNumber)/myresult1(n).Area;
    myfeatures(n).EdgeQualityPerimeterFactor = myresult1(n).Perimeter/myresult2(n).Perimeter;
end

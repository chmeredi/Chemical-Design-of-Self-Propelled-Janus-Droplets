function output=tiff_stack_merger(fname)
%reads a sequence of .tif images in a subfolder with name fname. Then
%saves as a multipage .tif file in folder that contains subfolder.
%Names the multipage .tif fname.tif

if nargin<1
    [filen, foldername] = uigetfile('*','Select one of the images that you want to merge to a stack.');
    fname=[foldername filen];
end

%read in first image
A1=imread(fname);
imsize=size(A1);

%read in all images and save in a cell called A
all_images = dir(fullfile(foldername, '*.tif'));
nfiles = length(all_images);
A  = NaN(imsize(1),imsize(2),nfiles);
for i = 1 : nfiles
   A(:,:,i) = imread(foldername, all_images(i).name);
end

output=A;
%write all images in Matlab home folder
% imwrite(A{1},[foldername '.tif']);
% for i = 2 : nfiles
%     imwrite(A{i},[foldername '.tif'],'writemode','append');
% end

end
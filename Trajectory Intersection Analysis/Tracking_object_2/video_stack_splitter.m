function video_stack_splitter(fname)
%splits .avi or tiffstack into individual images in folder
%Make sure to have a folder ready to which to copy files
%The name of this folder needs to be the same as the original filename

if nargin==0
    [filen, foldername] = uigetfile('*','Select image series you wish to split in individual units.');
    fname=[foldername filen];
end

%Recognize extension to check whether to read tifstack or avi.
filen_type=fname(length(fname)-3:length(fname));
%name of folder to put single tifs
fsave = fname(1:length(fname)-4);
disp(['Saving individual files into folder ' fsave ])

%If dealing with tiffstack
if filen_type=='.tif' | filen_type=='tiff'
    
    %Get info on size
    info = imfinfo(fname);
    num_images = numel(info);
    
    %split and save individual images
    for i=1:num_images
        if ~mod(i/10,1)
            disp(['Only ' num2str(num_images-i) ' to go']);
        end
        A=imread(fname,i);
        filename=[fsave '\image_' num2str(i,'%04d') '.tif'];
        imwrite(A,filename);
    end
    
%If dealing with avi
elseif filen_type=='.avi'
    
    fname
    %get info on size
    v=VideoReader(fname);
    i=0;
    
    while hasFrame(v)
        i=i+1;
        if ~mod(i/10,1)
            disp(['Working on image ' num2str(i) '.']);
        end
        img = readFrame(v);
        filename=[fsave '\image_' num2str(i,'%04d') '.tif'];
        imwrite(img,filename);
    end
elseif filen_type=='.mp4'
    
    fname
    %get info on size
    v=VideoReader(fname);
    i=0;
    
    while hasFrame(v)
        i=i+1;
        if ~mod(i/10,1)
            disp(['Working on image ' num2str(i) '.']);
        end
        img = readFrame(v);
        filename=[fsave '\image_' num2str(i,'%04d') '.tif'];
        imwrite(img,filename);
    end    
else
    disp('The extension of your filename is not recognizes as either tif or avi. Cannot proceed.');
    disp('Try converting into tif or avi first or explicitly type the correct extension in the name.');
end

end %end function

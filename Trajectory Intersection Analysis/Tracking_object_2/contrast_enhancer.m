function contrast_enhancer()

[filename, foldername] = uigetfile('*.tif','Select tif stack');
A=[foldername filename];

A1=imread(A,1);

A1 = uint8((double(A1)-double(min(A1(:))))/(double(max(A1(:)))-double(min(A1(:))))*255);
figure
imagesc(A1)
colorbar
final='n';

while final~='y'
    A1t=A1;
    messagen='What do you want to use as a lower cutoff in pixel value?';
    lowl=get_numerical_input(messagen);
    A1t(find(A1t(:,:)<lowl))=lowl;
    
    figure
    imagesc(A1t)
    colorbar
    
    messagen='What do you want to use as an upper cutoff in pixel value?';
    highl=get_numerical_input(messagen);
    A1t(find(A1t(:,:)>highl))=highl;
       
    figure
    imagesc(A1t)
    colorbar
    
    message='Are you happy with the processing: (y or n) ?';
    acceptables=['y' 'n'];
    final=get_textual_input(message,acceptables);
end

close all

info = imfinfo(A);
num_images = numel(info);
fname=[foldername '2' filename];

imwrite(A1t,fname);
for n=2:num_images
    
    if ~mod(n/100,1)
        disp(['Working on frame: ' num2str(n) '.'])
    end
    
    A1t=imread(A,n);
    A1t = uint8((double(A1t)-double(min(A1t(:))))/(double(max(A1t(:)))-double(min(A1t(:))))*255);
    A1t(find(A1t(:,:)<lowl))=lowl;
    A1t(find(A1t(:,:)>highl))=highl;
    imwrite(A1t,fname,'WriteMode','append')
end


end
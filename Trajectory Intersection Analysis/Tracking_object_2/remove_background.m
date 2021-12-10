function remove_background(filename)
    %Takes as input the filename of a stack of tif images, averages those
    %images over time and removes the background
    %Outputs the same file with suffix -bs (background substracted)
    %
    
    %make sure input is provided
    if nargin<1
        [filen, foldername_temp] = uigetfile('*.tif','Select tif stack');
        filename=[foldername_temp filen];
    end
    
    %split up filename in parts
    startIndex=max(regexp(filename,'\'));
    foldername_temp=filename(1:startIndex);
    filen=filename(startIndex+1:end-4);
    
    %get number of images
    info = imfinfo(filename);
    num_images = numel(info);
    %num_images = 10;
    
    %get imagesize
    A1=imread(filename,1);
    sizeA1=size(A1);
    width=sizeA1(1);
    length=sizeA1(2);
    
    Asum=zeros(width,length);
    maxval=0;
    minval=255;
   
    disp('Calculating average image')
    %average images
    for p=1:num_images
        Ap=imread(filename,p);
        Ap=im2double(Ap);
        Asum=Asum+Ap;
        if ~mod(p/100,1)
            disp(['Working on frame: ' num2str(p) '.'])
        end
        %keep track of max and min value
        maxval_temp=max(max(Ap));
        minval_temp=min(min(Ap));
        if maxval_temp>maxval
            maxval=maxval_temp;
        end
        if minval_temp<minval
            minval=minval_temp;
        end
    end
    Amean=Asum./num_images;
    minmean=min(min(Amean));
    maxmean=max(max(Amean));
    minval=minval-minmean;
    maxval=maxval-maxmean;
    
    %substract from image1;
    A1=im2double(A1);
    A1fix=A1-Amean;
    %A1fix=A1fix-minval;
    A1fix=A1fix-min(min(A1));
    %A1fix= A1fix./maxval;
    A1fix=A1fix./(max(max(A1)));
    A1fix = A1fix*255;
    A1fix = uint8(A1fix);
    
    disp('Putting together tiffstack of substracted images')   
    %save image 1
    imwrite(A1fix,[foldername_temp filen '-bs.tif'])
    
    for p=2:num_images
        %substsract mean
        Ap=imread(filename,p);
        Ap=im2double(Ap);
        Apc = Ap-Amean;
    
        %enhance contrast
        %Apc = Apc-minval;
        Apc = Apc - min(min(Apc));
        %Apc = Apc./maxval;
        Apc = Apc./max(max(Apc));
        Apc = Apc*255;
        Apc = uint8(Apc);
        
        if ~mod(p/100,1)
            disp(['Working on frame: ' num2str(p) '.'])
        end
        
        try imwrite(Apc,[foldername_temp filen '-bs.tif'],'WriteMode','Append')
        catch
            try imwrite(Apc,[foldername_temp filen '-bs.tif'],'WriteMode','Append')
            catch
                disp(['Missing frame ' num2str(p) 'because of writing permission']);
            end
        end
        
    end
    
end
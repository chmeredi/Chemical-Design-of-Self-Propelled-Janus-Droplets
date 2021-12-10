function export_property_to_text(property, filename)
    %Acts on an individual object property of type particles, takes the indicated
    %property and saves them as a text file
        
    %Find out the number of columns
    count=length(property(1,:));
    
    name=[];
    name_t='%8.4f\t';
    name_end='\r\n';
    
    for i=1:count
        name=[name name_t];
    end
    
    name=[name name_end];
    
    %Save trace as textfile in designated folder
            
    fsave=filename;
    fileID = fopen(fsave,'w');
    fprintf(fileID, name , property');
    fclose(fileID);
        
end
function export_trace_to_text(objectlist, folder)
    %Acts on a 1xn cell of objects of type particles, takes the trace
    %property and saves them into individual text files
    %trace identities are organized as x y t ID
    %ID runs from first to last particle in one chain
    %textfile number indicates chain number.
    
    %For each item in list
    for p=1:length(objectlist)
        
        %Obtain trace
        object_p=objectlist{p};
        trace_p=object_p.tr;
        
        %Save trace as textfile in designated folder
        textfile=trace_p;
        fsave=[folder '\chain_' num2str(p,'%04d') '.txt' ];
        fileID = fopen(fsave,'w');
        fprintf(fileID , '%8.4f\t %8.4f\t %8.4f\t %8.4f\r\n' , textfile');
        fclose(fileID);
        
    end
end
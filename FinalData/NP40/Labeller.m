function [label_3110] = Labeller(filename);


[~,pedo_max,~,~,~] = pedo_extract(filename);
label_3110 = zeros(size(pedo_max));
%BigToe
for i=1:9
    for j=1:9
        if pedo_max(i,j)>0
            label_3110(i,j)=1;
        end
    end
end

%Medial
for i=11:14
    for j=1:11
        if pedo_max(i,j)>0
            label_3110(i,j)=2;
        end
    end
end
for i=15:21
    for j=1:10
        if pedo_max(i,j)>0
            label_3110(i,j)=2;
        end
    end
end
%Lateral
k=0;
for i=11:14
    for j=12:(13+k)
       if k<6
         k=k+1;
       end  
       if pedo_max(i,j)>0
          label_3110(i,j)=3;   
       end       
    end
end

%Lateral part 2
for i=15:23
    for j=11:(13+k)
       if k<6
         k=k+1;
       end  
       if pedo_max(i,j)>0
          label_3110(i,j)=3;   
       end       
    end
end

%Heel
for i=41:size(label_3110,1)
    for j=1:size(label_3110,2)
        if pedo_max(i,j)>0
            label_3110(i,j)=4;
        end
    end
end

disp(label_3110)
end



       
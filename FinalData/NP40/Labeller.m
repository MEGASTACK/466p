function [label_3108] = Labeller(filename);


[~,pedo_max,~,~,~] = pedo_extract(filename);
label_3108 = zeros(size(pedo_max));
%BigToe
for i=1:9
    for j=1:9
        if pedo_max(i,j)>0
            label_3108(i,j)=1;
        end
    end
end

%Medial
k=0
for i=11:22
    for j=1:12-k
        if pedo_max(i,j)>0
            label_3108(i,j)=2;
        if k<3
            k=k+1;
        end
    end
end

%Lateral
k=0;
r=0;
for i=13:22
     if i>14
        r=r+1;
     end
     if i>17
        r=r+1;
     end
     if i>21
        r=r+1;
     end
    for j=(13-r):(14+k)
        if j==-1 | j==0
            j=19;
        end
        if pedo_max(i,j)>0
            label_3108(i,j)=3;
            if k<5
                k=k+1;
            end           
        end
    end
end

%Lateral part 2
for i=23:26
    for j=11:size(label_3108,2)
        if pedo_max(i,j)>0
            label_3108(i,j)=3;
        end
    end
end

%Heel
for i=40:size(label_3108,1)
    for j=1:size(label_3108,2)
        if pedo_max(i,j)>0
            label_3108(i,j)=4;
        end
    end
end

disp(label_3108)
end



       
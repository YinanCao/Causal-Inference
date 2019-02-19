function out=errorbar_YC(x,y,se,w,lw,color)

if size(x,2)>1
    x=x';
end

if size(y,2)>1
    y=y';
end

if size(se,2)>1
    se=se';
end

z=[y-se,y+se];
for i=1:length(x)
    out=plot([x(i),x(i)],[z(i,1),z(i,2)],'-','color',color,'linewidth',lw);
    hold on;
end

m=[x-w/2, x+w/2];

for i=1:2
    y=z(:,i);
    
    for j=1:length(x)
        plot([m(j,1),m(j,2)],[y(j),y(j)],'color',color,'linewidth',lw)
        hold on;
    end
end


end
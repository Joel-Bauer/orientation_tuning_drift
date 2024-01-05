function out=circ_diff(x,dim)
% Uses the circ_dist.m function to calculate circualar distances between
% consecutive values in a vector

if dim==1
    out = nan(size(x,1)-1,size(x,2));
    dim_stay=2;
elseif dim==2
    out = nan(size(x,1),size(x,2)-1);
    dim_stay=1;
end

for i = 1:size(x,dim)-1
    for ii = 1:size(x,dim_stay)
        if dim==1
            out(i,ii)=circ_dist(x(i,ii),x(i+1,ii));
        elseif dim==2
            out(ii,i)=circ_dist(x(ii,i),x(ii,i+1));
        end
    end
end
end
function out = makeComb(input),
N = length(input);
searchsize = 1;
search = 1;
for mm=1:N,
    searchsize = searchsize*input(mm);
    search(mm+1) = search(mm)*input(mm);
end
for uu=1:searchsize, 
    for vv=1:N,
        comb(uu,vv) = mod(round(uu/search(vv)),input(vv))+1;
    end
end
out = comb;
end
function resp = Pertence_a(x,v)

pos = v==x;

if sum(pos)>0
    resp = true;
else
    
    resp = false;
end
function y=extracttask(stim,x)
    if x==0
        stim1=stim;
        stim1(stim1==x)=10;
        stim1(stim1<10)=0;
        stim1(stim1==10)=1;
        y=stim1;
    else
        stim1=stim;
        stim1(stim1>x)=0;
        stim1(stim1<x)=0;
        stim1(stim1==x)=1;
        y=stim1;
    end
    
    
end
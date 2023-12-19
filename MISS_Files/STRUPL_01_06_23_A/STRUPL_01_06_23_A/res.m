function [RS]=res(ny,ragnetto,icrk,ncrl, NCRNPL,theta,INFLEN,th,twonotch)
%
% computes vector of plastic resistances
%
    for i=1:ny
        RS(i,1)=0;
    end
    if ragnetto==1
        RS(1,1)=1;
        RS(2,1)=1;
        RS(3,1)=8;
        RS(4,1)=8;
     %   RS(5,1)=130;
        RS(5,1)=0;
        RS(6,1)=130;
        RS(7,1)=0;
    end %if ragnetto
    if icrk ==1
        sigmat=0.1 %(N/mm2)
        if twonotch==1
            sigmat=3.4
        end
% 
        index=0;
            for i=1:ncrl
                ncrnli=NCRNPL(i);
                    for j=1:ncrnli
                            res=sigmat*sin(theta)*th*INFLEN(i,j);
%
                            RS(index+1)=0.3*res/sin(theta);
                            if twonotch==1
                                RS(index+1)=sigmat*th*INFLEN(i,j);
                                RS(index+2)=0;
                            end
                            if twonotch==0
                        %    RS(index+2)=0;
                         %   RS(index+3)=res;
                            RS(index+2)=res;
                         %   RS(index+4)=0;      
                            RS(index+3)=res;
                         %   RS(index+5)=res;
                         %   RS(index+6)=0;
                             index=index+3;
                         %   index=index+6;
                            end %if twonotch
                         if twonotch==1
                             index=index+2;
                         end
                     end % for j
            end %for i
    end %if icrk
end

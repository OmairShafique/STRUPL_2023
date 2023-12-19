    function[ELSPERCRNODE]= elpercrnode(nodenumber,...
        COORD,ELNODES,ELSPERNODE,tetarad,maxelspernode,le)
%
% computes the element indeces belonging to one side of the crack at a
% given node
%
% First it computes the distance of the crak line with inclination theta,
% passing througth the given node, from the reference system origin;
%
% number of the considered node
%
    in= nodenumber;
%
% theta angle: angle, with respct horizontal x1 axis, of the normal versor
% to the crack direction (the versus is not important);
%
% le=dimensione tipica di un elemento
%
    si=sin(tetarad)
    co=cos(tetarad)
    xi1=COORD(in,1)
    xi2=COORD(in,2)
    erre=(xi1*co+xi2*si)
    inc=0;
    nodei=nodenumber;
     %
         for i=1:maxelspernode
            flag=0;
            iel=ELSPERNODE(nodei,i)
            for j=1:3
                    if flag==0
                        nofel=ELNODES(iel,j)
                            if nofel==nodei
                                flag2=1
                            end
                        if flag2==0
                             x1=COORD(nofel,1);
                             x2=COORD(nofel,2);
                             rnofel=((x1*co+x2*si))
                             dist=((rnofel-erre)/le)
                                flag4=0
                                if abs (dist)<0.0001
                                 flag4=1
                                end
                             if flag4==0
                                if dist>0
                                    inc=inc+1
                                    ELSPERCRNODE(inc)=iel
                                    flag=1;
                                end
                             end
                             
                         end
                      end
                    flag2=0
                end
            end
        end











  
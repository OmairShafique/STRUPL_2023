
function[ONESELPERCRN,maxinc]= oneselpercrn(ncrl,NCRNPL,NODESPERCRL,COORD,...
        ELNODES,ELSPERNODE,BETA,maxelpernode)
%
% computes the element numbers belonging to one side of the crack at a
% given node, matrix ONESELPERCRN(n.of crack nodes,max number od elements on
% one side of the cracks)
%
% First it computes the distance of the crak line with inclination theta,
% passing througth the given node, from the reference system origin;
%
% number of the considered node
%
    maxinc=0;
    ncolones=0;
    rowindex=0;
 for h=1:ncrl
        br=BETA(h);
        nofn=NCRNPL(h);
   for k=1:nofn % loop sui nodi della linea di frattura h considerata
            rowindex=rowindex+1;
            nodei= NODESPERCRL(h,k); % numero del nodo della frattura h considerato
            si=sin(br);
            co=cos(br);
            xi1=COORD(nodei,1);
            xi2=COORD(nodei,2);
            erre=(xi1*co+xi2*si);
            inc=0;
       for i=1:maxelpernode  %loop sugli elementi che convergono al nodo nodei
        iel=ELSPERNODE(nodei,i); %iel=numero di elemento che contiene il nodo nodei 
         if iel>0
            flag2=0;
           for j=1:3  %loop sui 3 nodi dell'elemento
             if flag2==0;
                 niel=ELNODES(iel,j); %node number of element iel
                 flag1=0;
                   for r=1:nofn %loop sui nodi della linea di frattura
                     if flag1==0
                      nnum=NODESPERCRL(h,r);
                        if nnum==niel; %il nodo considerato appartiene alla linea di frattura
                                        %si termina il loop sull'indice r e si aumenta l'indice j
                          flag1=1;
                         end %if nnum
                      end % if flag1=1
                    end %end for r
                        if flag1==0; % si è trovato un nodo non appartenete alla linea di fessura
                            x1=COORD(niel,1);
                            x2=COORD(niel,2);
                            rnofel=x1*co+x2*si;
                            dist=erre-rnofel;
                                if dist>0
                                    inc=inc+1;
                                    if inc>maxinc
                                        maxinc=inc;
                                    end
                                    flag2=1;
                                    ONESELPERCRN(rowindex,inc)=iel;
                                end %if dist
                        end % if flag1
               end %if flag2
             end % for j
          end % if iel
        end % for i
    end % for k
   end %for h
  end
    











  
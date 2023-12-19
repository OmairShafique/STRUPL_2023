%

ne=1
nn=3
ndf=2

COORD(1,1)=0
COORD(1,2)=1
COORD(2,1)=2
COORD(2,2)=1
COORD(3,1)=0
COORD(3,2)=0

ELNODES(1,1)=1 
ELNODES(1,2)=3
ELNODES(1,3)=2


ELTYPE(1)=3

gama=1800

VOL(1)=0.5*2*1*1

DIRG(1)=0
DIRG(2)=-1

[FS]=gravityloads( COORD,ELNODES,ELTYPE,nn,ndf,ne,gama,VOL,DIRG)

return


    
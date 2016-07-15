import sys
import numpy as np

def htmltab(plots, file=None, xtitle=None, ytitle=None, size=100, header=None) :

    p=np.array(plots)
    nx=p.shape[1]
    ny=p.shape[0]

    if file is not None : 
        f = open(file,'w') 
    else : 
        f=sys.stdout
    f.write('<HTML><BODY>\n')
    if header is not None :
        f.write(header+'<p>\n')
    f.write('<TABLE BORDER=2>\n')
    if xtitle is not None :
        f.write('<TR>\n')
        if ytitle is not None :
            f.write('<TD>\n')
        for ix in range(nx) :
            f.write('<TD>'+xtitle[ix]+'\n')
    for iy in range(ny) :
        f.write('<TR>\n')
        if ytitle is not None :
            f.write('<TD>'+ytitle[iy]+'\n')
        for ix in range(nx) :
            f.write('<TD>\n')
            f.write('<A HREF='+p[iy][ix]+'>'+
                    '<IMG SRC='+p[iy][ix]+' WIDTH='+str(size)+'%></A>\n')
            f.write('</TD>\n')
    f.write('</TABLE>\n')
    f.write('</BODY></HTML>\n')


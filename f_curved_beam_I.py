import numpy as np

def f_curved_beam_I(ri,ro,t,ti,to,bi,bo):
# Variables d'entrada (segons figura): ri, ro, t, ti, to, bi, bo
# Variables de sortida: radi centroide (rc), radi fibra neutre (rn), 
#                       àrea (A), i moment de segon ordre (I)    h = ro - ri
    h = ro-ri 
    A = ti*(bi-t) + to*(bo-t) + h*t

    rc= ri + (1/2*h**2*t + 1/2*ti**2*(bi-t) + to*(bo-t)*(h-to/2))/A
    rn= A/(bi*np.log((ri+t)/ri)+t*np.log((ro-to)/(ri+ti))+bo*np.log(ro/(ro-to)))

    I =(t*h**3/12       + (h*t)*(rc-ri-h/2)**2 + 
        ti**3*(bi-t)/12 + (ti*(bi-t))*(rc-ri-ti/2)**2 + 
        to**3*(bo-t)/12 + (to*(bo-t))*(rc-ri-(h-to/2))*2)
    return rc,rn,A,I
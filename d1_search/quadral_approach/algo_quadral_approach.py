
import numpy as np


def quadral_approach(func, single_interval, maxIter=None, eps_f=None, eps_x=None):
    '''
      func: function to search for
      single_interval: 'U' shape interval | type=tuple
      maxIter: maxinum iteration
      eps_f: tolerence of diff f
      eps_x: tolerence of diff x
    '''
    if not isinstance(single_interval, tuple):
       raise(Exception)

    # init
    x1, x2 = single_interval
    f1 = func(x1)
    f2 = func(x2) 

    x3 = (x1 + x2)/2.0
    f3 = func(x3)

    itera = 0 
    while 1:
   
       itera += 1
       if maxIter and itera > maxIter:
          break

       if eps_x and np.abs(x1-x2) <= eps_x:
          break

       if eps_f and np.abs(f1 - f2) <= eps_f:
          break
 
       x4 = quadral_minx(x1,f1,x2,f2,x3,f3)
       f4 = func(x4)

       if x3 > x4: # make sure x3 <= x4
          x3, x4 = x4, x3       # then swap x3 & x4, f3 & f4, use 'serial swap' of Python
          f3, f4 = f4, f3
       if f3 > f4:
          x1 = x3; f1 = f3
          x3 = x4; f3 = f4
       else:
          x2 = x4; f2 = f4
    
    return x4
              
def quadral_minx(x1,y1,x2,y2,x3,y3):
    '''
     use quadral function and x1,y1 x2,y2 x3,y3 pair to find min point of the simulated function 
    '''
    minx = 0.5 * \
           ((x2 * x2 - x3 * x3) * y1 + (x3 * x3 - x1 * x1) * y2 + (x1 * x1 - x2 * x2) * y3) \
         / ((x2 - x3) * y1 + (x3 - x1) * y2 + (x1 - x2) * y3)    
    return minx


if __name__ == "__main__":
   result = quadral_approach(np.sin, (4,5), maxIter=None, eps_f=None, eps_x=0.00001)
   print result

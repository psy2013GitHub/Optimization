
import numpy as np

def golden_seg(func, single_interval, maxIter=None, eps_f=None, eps_x=None):
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
    a, b = single_interval
    f_a = func(a)
    f_b = func(b)
    
    a1 = a + 0.382 * (b - a); f_a1 = func(a1)
    b1 = a + 0.618 * (b - a); f_b1 = func(b1)
    
    itera = 0   
    while 1:

       itera += 1
       if maxIter and itera > maxIter:
          break
   
       if eps_x and np.abs(a-b) <= eps_x:
          break
      
       if eps_f and np.abs(f_a - f_b) <= eps_f:
          break

       if f_a1 > f_b1: # del [a, a1]
          a = a1; f_a = f_a1
          a1 = b1; f_a1 = f_b1
          b1 = a + 0.618 * (b - a); f_b1 = func(b1)
       else:           # del [b1, b]
          b = b1; f_b = f_b1
          b1 = a1; f_b1 = f_a1
          a1 = a + 0.382 * (b - a); f_a1 = func(a1) 
    
    return (a+b)/2.0
    
if __name__ == "__main__":
   def func(x):
       result = np.power(x, 2) - 10 * x + 36
       return result
   print golden_seg(func, (-10,10), maxIter=None, eps_f=None, eps_x=0.000001)   

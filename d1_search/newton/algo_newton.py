
import numpy as np

def newton(func_diff1, func_diff2, start=0, func_raw=None, maxIter=None, eps_f=None, eps_x=None):
    '''
      func_diff1: diff func of func_raw
      func_diff2: diff func of func_diff1
      maxIter: maxinum iteration
      eps_f: tolerence of diff f
      eps_x: tolerence of diff x

    reason: 2nd-order Taylor approach
    pro: cmpred to other 1d search algo, newton didnot need to search single_interval
    con: need diff1 & diff2 function & start point
    '''

    x0 = start 
    itera = 0
    while 1:
   
       itera += 1
       if maxIter and itera > maxIter:
          break

       last_x0 = x0
       x0 -= func_diff1(x0) / float(func_diff2(x0))
               
       if eps_x and np.abs(last_x0-x0) <= eps_x:
          break

       if eps_f and np.abs(func(last_x0) - func(x0)) <= eps_f:
          break       

    return x0
              

if __name__ == "__main__":
   def func(x):
       return np.power(x,4) - 4 * np.power(x,3) - 6 * np.power(x,2) - 16 * x + 4
   def func_diff1(x):
       return 4 * np.power(x,3) - 12 * np.power(x,2) - 12 * x - 16
   def func_diff2(x):
       return 12 * np.power(x,2) - 24 * x - 12

   result = newton(func_diff1, func_diff2, start=3, func_raw=None, maxIter=None, eps_f=None, eps_x=0.0001)
   print result
   print func(result)

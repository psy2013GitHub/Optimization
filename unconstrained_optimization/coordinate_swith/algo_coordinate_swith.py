
import numpy as np


def algo_coodinate_switch(func, start=0, maxIter=None, eps_f=None, eps_x=None):
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
   def func(x1, x2):
       return 3 * np.power(x1, 2) + np.power(x2, 2) - 2 * x1 * x2 + 4 * x1 + 3 * x2
   

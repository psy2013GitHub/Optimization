
import numpy as np

import sys
sys.path.append("..")

def newton(diff_func_array, hessian_func_matrix, start=(0), method = "alpha_newton", func_raw=None, maxIter=None, eps_f=None, eps_x=None):
    '''
      diff_func_array: diff func array of func_raw
      hessian_func_matrix: hessian matrix of diff_func_matrix
      method: "alpha_newton" or "newton", see __init__.py
      maxIter: maxinum iteration
      eps_f: tolerence of diff f
      eps_x: tolerence of diff x
    '''
    
    nX = len(hessian_func_matrix)
   
    x0 = np.array(start, dtype=np.float) 
    itera = 0
    while 1:
   
       itera += 1
       if maxIter and itera > maxIter:
          break

       last_x0 = x0
       
       S = np.matrix([j(*x0) for i in hessian_func_matrix for j in i], dtype=np.float).reshape((nX,nX)).getI()  * \
                np.matrix([i(*x0) for i in diff_func_array], dtype=np.float).transpose() 

       if method=="newton":
          alpha = 1
       elif method=="alpha_newton":
          def func_alpha(func_raw, x0, S):
              # closure
              def func(alpha):
                  return func_raw(x0 - alpha * S)
              return func
          from d1_search.single_interval_search import algo_back_forth # fix me for import
          single_interval = algo_back_forth(func, start=0, h=0.1, alpha_forth=2, alpha_back=1, maxIter=10000, eps_f=None, eps_x=None) 
          from d1_search.golden_segment import algo_golden_segment
          alpha = golden_seg(func_alpha(func_raw, x0, S), single_interval, maxIter=None, eps_f=None, eps_x=None)
       else:
          print "only 'newton' or 'alpha_newton' supplied"  
          return None
       
       S = np.array(S).reshape(x0.shape) # convert to array
       x0 -= alpha * S
               
       if eps_x and np.all(np.abs(last_x0-x0) <= eps_x):
          break

       if eps_f and np.all(np.abs(func(last_x0) - func(x0)) <= eps_f):
          break       

    return x0
              

if __name__ == "__main__":

   # func
   def func(*nk_args):
       x1, x2, x3 = nk_args
       return np.power(x1 - x2 + x3, 2) + np.power(x2 - x1 + x3, 2) + np.power(x1 + x2 + x3, 2)
   # diff_func_array
   def diff_func_x1(*nk_args):
       x1, x2, x3 = nk_args
       return 2 * (3 * x1 - x2 + x3)
   def diff_func_x2(*nk_args):
       x1, x2, x3 = nk_args
       return 2 * (3 * x2 - x1 + x3)
   def diff_func_x3(*nk_args):
       x1, x2, x3 = nk_args
       return 2 * (3 * x3 + x1 + x2)
   diff_func_lst = [diff_func_x1, diff_func_x2, diff_func_x3]
   # hessian matrix
   def diff_func_x1_x1(*nk_args):
       x1, x2, x3 = nk_args
       return 6
   def diff_func_x1_x2(*nk_args):
       x1, x2, x3 = nk_args
       return -2
   def diff_func_x1_x3(*nk_args):
       x1, x2, x3 = nk_args
       return 2
   def diff_func_x2_x1(*nk_args):
       x1, x2, x3 = nk_args
       return -2
   def diff_func_x2_x2(*nk_args):
       x1, x2, x3 = nk_args
       return 6
   def diff_func_x2_x3(*nk_args):
       x1, x2, x3 = nk_args
       return 2
   def diff_func_x3_x1(*nk_args):
       x1, x2, x3 = nk_args
       return 2
   def diff_func_x3_x2(*nk_args):
       x1, x2, x3 = nk_args
       return 2
   def diff_func_x3_x3(*nk_args):
       x1, x2, x3 = nk_args
       return 6

   hessian_func_lst=[[diff_func_x1_x1,diff_func_x1_x2,diff_func_x1_x3],\
                     [diff_func_x2_x1,diff_func_x2_x2,diff_func_x2_x3],\
                     [diff_func_x3_x1,diff_func_x3_x2,diff_func_x3_x3]]



   result = newton(diff_func_lst, hessian_func_lst, start=(0,0,0), method = "newton", func_raw=None, maxIter=None, eps_f=None, eps_x=0.0001)
   print result

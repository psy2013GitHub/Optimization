
import numpy as np
import copy
from d1_search.single_interval_search import algo_back_forth # fix me for import
from d1_search.golden_segment import algo_golden_segment

def DFP(func, diff_func_array, start=(0), maxIter=None, eps_f=None, eps_x=None):
    '''
      diff_func_array: diff func array of func_raw
      hessian_func_matrix: hessian matrix of diff_func_matrix
      method: "Powell" or "modified_Powell", see __init__.py
      maxIter: maxinum iteration
      eps_f: tolerence of diff f
      eps_x: tolerence of diff x
    '''
    
    # init
    ndim = len(start)

    x0 = np.array(start, dtype=np.float).reshape((ndim,1)); last_x0 = x0
    M = np.eye(ndim, dtype=np.float)
    itera = 0 
    while 1:
       
       itera += 1
       if maxIter and itera > maxIter:
          break         

       # seek alpha which minimize f(x1) in the direction of S      
       #single_interval = algo_back_forth(func, start=0, h=0.1, alpha_forth=2, alpha_back=1, maxIter=10000, eps_f=None, eps_x=None)
       #alpha = golden_seg(func_alpha(func_raw, x0, S), single_interval, maxIter=None, eps_f=None, eps_x=0.0001)
       alpha = 0.3

       df_x0 =  np.array([i(*tuple(x0.reshape(ndim).tolist())) for i in diff_func_array], dtype = np.float).reshape((ndim,1)) 
       S = -1 * M.dot(df_x0)  
       x1 = x0 + alpha * S

       if itera != 1:
          
          df_x1 = np.array([i(*tuple(x1.reshape(ndim).tolist())) for i in diff_func_array], dtype = np.float).reshape((ndim,1))   
          delta_g = df_x1 - df_x0
          delta_x = x1 - x0          

          E = delta_x.dot(delta_x.transpose()) / float((delta_x.transpose().dot(delta_g))) - \
              M.dot(delta_g.dot(delta_g.transpose())).dot(M) / float(delta_g.transpose().dot(M).dot(delta_g))  

          M = M + E
        
       if eps_x and np.all(np.abs(x1 - x0)) <= eps_x:
          break
       if eps_f and np.abs(func(*tuple(x1.reshape(ndim).tolist())) - func(*tuple(x0.reshape(ndim).tolist()))) <= eps_f:
          break

       x0 = x1
    
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
 
   print DFP(func, diff_func_lst, start=(1,2,3), maxIter=None, eps_f=0.00001, eps_x=None)

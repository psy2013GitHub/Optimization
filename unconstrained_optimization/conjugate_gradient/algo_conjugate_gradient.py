
import numpy as np
from d1_search.single_interval_search import algo_back_forth # fix me for import
from d1_search.golden_segment import algo_golden_segment

def conjugate_gradient(diff_func_array, hessian_func_matrix, start=(0), method = "FR_conjugate_gradient", func_raw=None, maxIter=None, eps_f=None, eps_x=None):
    '''
      diff_func_array: diff func array of func_raw
      hessian_func_matrix: hessian matrix of diff_func_matrix
      method: "alpha_newton" or "newton", see __init__.py
      maxIter: maxinum iteration
      eps_f: tolerence of diff f
      eps_x: tolerence of diff x
    '''

    x0 = np.array(start, dtype=np.float); last_x0 = x0
    itera = 0 
    while 1:
       
       itera += 1
       if maxIter and itera > maxIter:
          break         

       # dierection S0, i.e. -diff_f(x0), the negative gradient of x0
       S0 = -1 * np.array([i(*x0) for i in diff_func_array], dtype=np.float)

       # seek alpha which minimize f(x1) in the direction of diff_f(x0)      
       single_interval = algo_back_forth(func, start=0, h=0.1, alpha_forth=2, alpha_back=1, maxIter=10000, eps_f=None, eps_x=None)
       alpha = golden_seg(func_alpha(func_raw, x0, S0), single_interval, maxIter=None, eps_f=None, eps_x=0.0001)

       # x1
       x1 = x0 + alpha * S0

       # direction S1
       if method == "conjugate_gradient":
          H_x1 = np.matrix([j(*x1) for i in hessian_func_matrix for j in i], dtype=np.float).reshape((nX,nX)) 
          A = np.matrix(S0).transpose() * H_x1
          S1 = # t(S0) * H_x1 * S1 = 0 => A * S1 = 0 => solve the linear equation  
       elif method == "FR_conjugate_gradient": # dnot need Hessian matrix
          gradient_x1 = 1 * np.array([i(*x1) for i in diff_func_array], dtype=np.float)
          S0_norm = np.sum(S0 * S0)
          if S0_norm < 0.00000001:
             return x0
          gradient_x1_norm = np.sum(gradient_x1 * gradient_x1)
          beta = gradient_x1_norm / S0_norm 
          S1 = -1 * gradient_x1 + beta * S0
       else:
          print "only 'conjugate_gradient' or 'FR_conjugate_gradient' supplied"
          return None  

       # seek alpha which minimize f(x1) in the direction of diff_f(x0)      
       single_interval = algo_back_forth(func, start=0, h=0.1, alpha_forth=2, alpha_back=1, maxIter=10000, eps_f=None, eps_x=None)
       alpha = golden_seg(func_alpha(func_raw, x1, S1), single_interval, maxIter=None, eps_f=None, eps_x=0.0001)

       # x_star
       x_star = x1 + alpha * S1

       if eps_x and np.all(np.abs(last_x0-x0) <= eps_x):
          break

       if eps_f and np.all(np.abs(func(last_x0) - func(x0)) <= eps_f):
          break
       
       last_x0 = x0
       x0 = x_star

    return x0

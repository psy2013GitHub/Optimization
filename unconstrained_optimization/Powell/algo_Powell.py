
import numpy as np
import copy
from d1_search.single_interval_search import algo_back_forth # fix me for import
from d1_search.golden_segment import algo_golden_segment

def Powell(func, start=(0), method = "modified_Powell", maxIter=None, eps_f=None, eps_x=None):
    '''
      diff_func_array: diff func array of func_raw
      hessian_func_matrix: hessian matrix of diff_func_matrix
      method: "Powell" or "modified_Powell", see __init__.py
      maxIter: maxinum iteration
      eps_f: tolerence of diff f
      eps_x: tolerence of diff x
    '''
    
    # init
    if method != "Powell" and method != "modified_Powell":
       print "only 'Powell' or 'modified_Powell' supplied"
       return 
       
    ndim = len(start)
    Direction_mat = np.eye(ndim, dtype=np.float) # linear independent vectors  

    x0 = np.array(start, dtype=np.float); last_x0 = x0
    itera = 0 
    while 1:
       
       itera += 1
       if maxIter and itera > maxIter:
          break         

       # got all n directions
       for loop in xrange(ndim):
           x = x0; max_f_descend_idx = None; delta = None
           for d in xrange(ndim):
           
               direction = Direction[:,d]
       
               # seek alpha which minimize f(x1) in the direction of S      
               single_interval = algo_back_forth(func, start=0, h=0.1, alpha_forth=2, alpha_back=1, maxIter=10000, eps_f=None, eps_x=None)
               alpha = golden_seg(func_alpha(func_raw, x0, direction), single_interval, maxIter=None, eps_f=None, eps_x=0.0001)

               # update x
               raw_f_x = func(*tuple(x.tolist()))
               x += alpha * direction 
               new_f_x = func(*tuple(x.tolist()))
               f_descend = new_f_x - raw_f_x
               if not max_f_descend_idx or delta < f_descend: # used for "modified_Powell"
                  max_f_descend_idx = d
                  delta = f_descend

            if eps_x and np.sqrt(np.sum(np.power(x - x0,2))) <= eps_x:
               return x
            if eps_f and np.abs(func(*tuple(x.tolist())) - func(*tuple(x0.tolist()))) <= eps_f:
               return x
              
           # got S
           S = x - x0

           # del first direction & append S according to the method used
           if method == "Powell":
              Direction_mat[:,:-2] = Direction_mat[:,1:-1]
              Direction_mat[:,-1] = S
              
              # seek alpha which minimize f(x1) in the direction of S      
              single_interval = algo_back_forth(func, start=0, h=0.1, alpha_forth=2, alpha_back=1, maxIter=10000, eps_f=None, eps_x=None)
              alpha = golden_seg(func_alpha(func_raw, x0, S), single_interval, maxIter=None, eps_f=None, eps_x=0.0001)
              
              x0 = x0 + alpha * S

           else:
              f1 = func(*tuple(x0.tolist()))
              f2 = func(*tuple(x.tolist()))
              x_refection = 2 * x - x0; f3 = func(*tuple(x_refection.tolist())) # reflection is the point of 'x0' with regards to 'x'
               
              if f3 > f1 and (f1 - 2 * f2 + f3) * np.power((f1 - f2 - delta), 2) >= delta * np.power((f1 - f3), 2) / 2.0:
                 x = x_reflection if f3 < f2 else x
                 continue
              else: # del max f descend direction & append S according to the method used 
                 Direction_mat[:,max_f_descend_idx:-2] = Direction_mat[:,max_f_descend_idx+1:-1]
                 Direction_mat[:,-1] = S
              
                 # seek alpha which minimize f(x1) in the direction of S      
                 single_interval = algo_back_forth(func, start=0, h=0.1, alpha_forth=2, alpha_back=1, maxIter=10000, eps_f=None, eps_x=None)
                 alpha = golden_seg(func_alpha(func_raw, x0, S), single_interval, maxIter=None, eps_f=None, eps_x=0.0001)
              
                 x0 = x0 + alpha * S


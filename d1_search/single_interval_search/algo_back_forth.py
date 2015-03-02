

import numpy as np


def algo_back_forth(func, start=0, h=0.1, alpha_forth=2, alpha_back=1, maxIter=None, eps_f=None, eps_x=None): 
    '''
     main algorithm for interval search 

     start: start point for interval search
     h: step
     alpha_forth: alpha for forth, need to >= 1
     alpha_back:  alpha for back,  need to <= 1
    '''
    
    result_lst = []
    itera = 0
    l_end = start; r_start = start


    # determine if the interval which includes "start" point is "single interval"
    tpl = search_min_interval(start-h, start, h, alpha_forth, alpha_back, direction="back")
    if tpl:
       result_lst.append(tpl)     
    else:
       tpl = search_min_interval(start-h, start, h, alpha_forth, alpha_back, direction="back")
       if tpl:
          result_lst.append(tpl)     

    # from start to -inf & +inf
    while 1:

      itera += 1
      if maxIter and itera > maxIter:
         break               

      l_start = l_end - h
      tpl = search_min_interval(l_start, l_end, h, alpha_forth, alpha_back, direction="back")
      if tpl:
         result_lst.append(tpl)
         l_end = l_start - alpha_back * h
      else:
         l_end = l_start
         
      
      r_end = r_start + h
      tpl= search_min_interval(r_start, r_end, h, alpha_forth, alpha_back, direction="forth") 
      if tpl:
         result_lst.append(tpl) 
         r_start = r_start + alpha_forth * h
      else:
         r_start = r_end 

    return result_lst

def search_min_interval(start, end, h, alpha_forth, alpha_back, direction=None):
    if func(start) > func(end):
       if direction == "forth":
          middle = end
          end = start + alpha_forth * h   # step forth
          if func(middle) < func(end):
             return (start, middle, end)
    else:  
       if direction == "back":
          middle = start
          start = start - alpha_back * h  # step back
          if func(middle) < func(start):
             return (start, middle, end)
    return None  

if __name__ == "__main__":
  def func(x):
      return np.power(x, 4) - np.power(x, 2) - 2 * x + 5
  result_lst = algo_back_forth(func, start=0, h=0.1, alpha_forth=2, alpha_back=1, maxIter=1000, eps_f=None, eps_x=None) 
  print result_lst
  print func(result_lst[0][0])
  print func(result_lst[0][1])
  print func(result_lst[0][2])
  # x = np.arange(-0.5,2,0.1)
  # y = func(x)
  # import matplotlib.pyplot as plt
  # plt.plot(x,y,'o')        
  # plt.show()

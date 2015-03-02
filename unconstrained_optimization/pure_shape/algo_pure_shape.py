
import numpy as np

def pure_shape(func, ndim=None, X=None, reflection_alpha=1.2, extend_gamma=2, squeeze_beta=0.3, shorten_edge_sita=0.5, maxIter=None, eps_f=None, eps_x=None):
    '''
      Usage: 
            1, pure_shape(func, ndim=3, reflection_alpha=1.2, extend_gamma=2, squeeze_beta=0.3, maxIter=*, eps_f=*, eps_x=*)
            2, pure_shape(func, X=[[x1,y1],[x2,y2],...], reflection_alpha=1.2, extend_gamma=2, squeeze_beta=0.3, maxIter=*, eps_f=*, eps_x=*)
    '''
    # init X
    if reflection_alpha < 1:
       print "reflection_alpha need > 1"
       return
    if extend_gamma <= reflection_alpha:
       print "extend_gamma need > reflection_alpha"
       return
    if squeeze_beta >= 1:
       print "squeeze_beta need < 1"
       return
    if shorten_edge_sita >= 1:
       print "shorten_edge_sita need < 1"
       return
 
    if not len(X):
       if not ndim:
          print "please suplly ndim"
          return
       X = np.zeros((ndim, ndim+1), dtype=np.float)
       X[:,0:ndim] = np.eye(ndim, dtype=np.float); X[:,ndim] = X[:,ndim] + 1
       Y = np.zeros((ndim,), dtype=np.float)
       for d in xrange(ndim+1): # ndim + 1 columns
           Y[d] = func(*tuple(X[:,d].tolist()))
    else:
       ndim = len(X[0])
       Y = np.array([func(*x) for x in X], dtype=np.float)
       X = np.array(X, dtype=np.float).transpose()
       if X.shape[1] != ndim + 1:
          print "please supply ndim + 1 vectors"
          return
    
    # iterate
    itera = 0
    while 1:
       
       itera += 1
       if maxIter and itera > maxIter:
          return X[:,0]
       
       # sort
       idx = np.argsort(Y); X = X[:, idx]; Y = Y[idx]
       X_l = X[:, 0]; f_X_l = Y[0] # min point
       X_g = X[:,-2]; f_X_g = Y[-2] # next_max point
       X_h = X[:,-1]; f_X_h = Y[-1] # max point  

       # weight center i.e. mass point except X_h
       X_f = np.sum(X[:,:-1], axis=1) / float(ndim)
       f_X_f = func(*tuple(X_f.tolist()))

       # break while 1
       if itera!=1:
          if eps_f and np.sqrt(np.sum(np.power(Y - f_X_f, 2)) / float(ndim + 1)) <= eps_f:
             return X[:,0]
          if eps_x and np.abs(np.sum(np.power(X - X[:,0].reshape((ndim,1)), 2))) <= eps_x:
             return X[:,0]
       
       # reflection
       X_r = X_f + reflection_alpha * (X_f - X_h); f_X_r = func(*tuple(X_r.tolist()))
        
       # go futher along reflection or squeeze ?
       if f_X_r < f_X_l: # go further
          X_e = X_f + extend_gamma * (X_f - X_h); f_X_e = func(*tuple(X_e.tolist()))
          if f_X_e < f_X_r:
             X[:, -1] = X_e; Y[-1] = f_X_e # replace X_h with X_e
             continue 
          else: # stay still
             X[:, -1] = X_r; Y[-1] = f_X_r # replace X_h with X_r
             continue
       else: # squeeze
          if f_X_r < f_X_g: # f_X_l <= f_X_r < f_X_g, replace X_h with X_r
             X[:, -1] = X_r; Y[-1] = f_X_r 
             continue
          else: # f_X_r >= f_X_g
             if f_X_r < f_X_h: # f_X_g <= f_X_r < f_X_h, squeeze
                X_s = X_f + squeeze_beta * (X_f - X_h); f_X_s = func(*tuple(X_s.tolist()))
             else: # f_X_r >= f_X_h, squeeze
                X_s = X_f - squeeze_beta * (X_f - X_h); f_X_s = func(*tuple(X_s.tolist()))
             if f_X_s < np.min((f_X_r, f_X_h)):
                X[:, -1] = X_s; Y[-1] = f_X_s # replace X_h with X_s
                continue
             else: # shorten edge
                X[:, 1:] = X[:, 0].reshape((ndim,1)) + shorten_edge_sita * (X[:, 1:] - X[:, 0].reshape((ndim,1))) # remeber: X[:,0] is always a row vector, so reshape
                continue   
                         


if __name__ == "__main__":
   def func(*nk_args):
       x1, x2 = nk_args
       return (x1 - 5) * (x1 - 5) + (x2 - 6) * (x2 - 6)                   
   result = pure_shape(func, X=[[8,9],[10,11],[8,11]], maxIter=None, eps_f=0.001, eps_x=None)
   print result 
   print func(*result)
              

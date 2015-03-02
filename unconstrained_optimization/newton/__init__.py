
# extend 1d_newton to nd_newton, then got "newton" method

# diff between "newton" & "alpha_newton" is the choise of alpha
# "newton" set alpha=1, thus suffer from the hazard of maximize f at some iteration
# "alpha_newton" set alpha to guarentee minimize f at each iteration


# pros: 
#      fast if diff matrix & inverse Hessian matrix easy to calc

# cons:
#      diff matrix & inverse Hessian matrix

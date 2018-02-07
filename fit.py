import numpy as np
from argparse import ArgumentParser
from lmfit import Parameters, Model
import matplotlib.pyplot as plt

def sigmoidal(x,rhol,rhov,d,l): 
  return (rhol + rhov)/2 - (rhol - rhov)/2 * np.tanh(2*(x-l)/d)

def generate_random_data(f,N,s,bounds,params):

  x = np.linspace(*bounds,num=N)
  y = f(x,**params) + s * (np.random.random(N) - 0.5)
 
  return x,y

def main():
  "Test lmfit package"
  mod = Model(sigmoidal)
  params = Parameters()
  params.add('rhol', value=0.5)
  params.add('rhov', value=0.1)
  params.add('d',    value=0.5)
  params.add('l',    value=0.5)

  x,y = generate_random_data(sigmoidal,300,0.3,(-2,2),{'rhol':0.85,'rhov':0,'d':0.5,'l':0})

  init = mod.eval(params, x=x)
  result = mod.fit(y,params,x=x)
  dely = result.eval_uncertainty(sigma=3)
  ci = result.conf_interval()

  fig, fit = plt.subplots()
  fit.plot(x,y,'o')
  fit.plot(x,result.best_fit)
  fit.fill_between(x, result.best_fit-dely, result.best_fit+dely, color='#ABABAB')
  fit.errorbar(result.best_values['l'], result.best_values['rhol']/2, xerr=ci['l'][6][1])
  plt.show()

  print(result.fit_report())
  print(result.ci_report())

if __name__ == '__main__':
  main()

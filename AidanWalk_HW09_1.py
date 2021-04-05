#Aidan Walk
#Astr 260-001
#HW09 Problem 1
#31 March 2021, 17:00

import numpy as np
import matplotlib.pyplot as plt

def analytic_solution(t, K=None, P0=None, r=None):
    """K = maximum carrying capacity of ecosystem
       P0 = initial population at time t=0
       r = population growth rate
       t = time """
    numerator = K*P0*np.exp(r*t)
    denominator = K + P0*(np.exp(r*t) - 1)
    return numerator / denominator

def dP_dt(P=None):
    '''P = Population size
       K = maximum carrying capacity of ecosystem
       r = population growth rate'''
    return r*P*(1 - P/K)

def forward_euler(timestep = None,
                  max_time = None,
                  initial_val = None,
                  deriv = None):
    '''Forward Euler's Method for Population Growth
       returns times, values
       timestep = 25
       max_time = 500
       initial_val = 1 billion
       deriv = dP/dt'''
    y = initial_val #initialize y 
    times = np.arange(0, max_time, timestep) #array of time steps 
    
    ypoints = []
    for t in times:
        ypoints.append(y)
        #compute new/next value of y
        y = y + timestep*deriv(P=y)
    
    return times, ypoints

def rk4(timestep = None,
        max_time = None,
        initial_val = None,
        deriv = None):
    '''Runge-Kutta 4 method for population growth
       returns times, vals
       timestep = 10
       max_time = 500
       initial_val = 1 billion
       deriv = dP/dt'''
    y = initial_val #initialize y 
    times = np.arange(0, max_time, timestep) #array from 0 to total time in steps of timestep
    h = timestep #for better readability
    ypoints=[]
    #compute y-val for each time in times
    for t in times:
        ypoints.append(y)
        #compute k1->4
        k1 = h*deriv(P= y)
        k2 = h*deriv(P= y + 1/2*k1)            
        k3 = h*deriv(P= y + 1/2*k2)
        k4 = h*deriv(P= y+k3)
        #update y position
        y = y + 1/6 * (k1 + 2*k2 + 2*k3 + k4)
        
    return times, ypoints

def createVisualization():
    '''Plot of Population Growth
       plots analytic solution, Eulers method, RK4 method
       Makes plot look pretty
       saves figure to local directory as AidanWalk_HW09_1_Plot.png'''
    #plot analytic solution
    plt.plot(years_since_start+1800, analytic_sol,
             label='Analytic Solution', color='#fe4010', linewidth=3)
    #plot Eulers method
    plt.plot(eulersTimes+1800, eulersVals,
             label="Eulers Method", color='#8c8c8c', linewidth=2)
    #plot RK4
    plt.plot(RK4Times+1800, RK4Vals,
             color='black', label='Runge-Kutta 4', linestyle=':')
             
    plt.title('Population Growth') #assigns plot title
    plt.xlabel('Years') #label x-axis
    #plt.xticks(yearRange)
    plt.ylabel('Population Size *10^9') #label y-axis
    legend = plt.legend(loc='lower right') #Display legend
    legend.get_frame().set_facecolor('#f8f9f9') #color legend background
    #save plot
    plotFName='AidanWalk_HW09_1_Plot'
    plt.savefig(plotFName +'.png')
    print('Plot saved to:', plotFName)
    
    return None

if __name__ == "__main__":
    K = 10 #10 billion
    P0 = 1 #billion, population at year 1800
    r = 0.014 #1.4% growth rate
    start_year = 1800 #Year of t0
    max_year = 2300 #year of tfinal
    max_time = max_year - start_year #total time interval in years
    
    #calculate analytic solution
    years_since_start = np.arange(0, max_time)
    analytic_sol = analytic_solution(years_since_start, 
                                     K=K, P0=P0, r=r)
    #calculate Euler's Method
    eulersTimes, eulersVals = forward_euler(initial_val=P0,
                                            timestep=25,
                                            max_time=max_time, deriv=dP_dt)
    #calculate Runge-Kutta 4 method
    RK4Times, RK4Vals = rk4(timestep=25, 
                            max_time=max_time,
                            initial_val = P0,
                            deriv = dP_dt)

    #make plot look pretty and save plot
    createVisualization()
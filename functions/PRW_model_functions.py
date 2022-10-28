import numpy as np
import pandas as pd

#calculate smoothed velocity
#input is array of x or y positions over time of one cell (output from PRW_polaritybias() function)
#returns array of smoothed velocity in either the x or y direction
def calc_vel(arr):
    series = pd.Series(arr)
    series_smooth = series.rolling(3,center=True).mean()
    series_smooth_dx = series_smooth.diff()
    vel = series_smooth_dx.rolling(2).mean().shift(-1)
    vel[np.isnan(vel)] = 0
    return vel

#run simulation for one cell over time
#input is standard devaition for theta angle, stepsize (v), timestep (dt), and total time of simulation
#return time series of x and y positions and list of times that correspond to each position
def PRW(std_dev, v, dt, time):
    #initialize time
    t = 0

    #initialize x and y positions
    pos_x = 0
    pos_y = 0 

    #create lists to store cell x and y positions and add initial positions
    all_pos_x = [pos_x]
    all_pos_y = [pos_y]

    #initialize theta (cell direction)
    prev_theta = np.random.uniform(0,2*np.pi)

    #create list to store theta values and add initial value
    all_theta = [prev_theta]

    #get new x and y position and add to x and y postion lists
    pos_x += v * np.sin(prev_theta)
    pos_y += v * np.cos(prev_theta)
    all_pos_x.append(pos_x)
    all_pos_y.append(pos_y)

    #progress time
    t += dt

    #loop through time and get x and y positions and theta for each time step
    while t<time:
        d_theta = np.random.normal(0,std_dev)
        theta = prev_theta + d_theta
        pos_x += v * np.sin(theta)
        pos_y += v * np.cos(theta)
        t += dt
        all_pos_x.append(pos_x)
        all_pos_y.append(pos_y)
        all_theta.append(theta)
        prev_theta = theta
        
    return all_pos_x, all_pos_y, theta

#function to run PRW simulation for many cells
#input is the number of cells to simulate (Nwalkers), time step (dt), total time of simulation (t),
#and standard deviation for theta angle
#returns a list of dataframes, where each cell has its own dataframe that contains x and y positions,
#x and y components of velocity, and theta
def run_PRW_sim(Nwalkers, dt, time, std_dev):
    #initialize list to store dataframes for each cell
    data_sim = []
    #loop through each cell and run simulation
    for i in range(0,Nwalkers):
        v = np.random.exponential(3)
        x,y,theta = PRW(std_dev, v, dt, time)
        x_vel = calc_vel(x)
        y_vel = calc_vel(y)
        vel = np.sqrt(x_vel**2 + y_vel**2)
        r = np.sqrt(np.array(x)**2 + np.array(y)**2)
        stepsize = np.sqrt((np.array(x[0:-1])-np.array(x[1:]))**2 +(np.array(y[0:-1])-np.array(y[1:]))**2)
        net_disp = abs(r[-1] - r[0])
        #tot_path_len = np.sum(r)
        onewalker = {'x': x, 'y': y, 'vx': x_vel, 'vy': y_vel, 'v': vel, 'theta': theta, 'DoverT': net_disp/np.sum(stepsize)}
        onewalker_df = pd.DataFrame(data=onewalker)
        data_sim.append(onewalker_df)
    return data_sim
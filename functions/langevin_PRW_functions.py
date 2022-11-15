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
#input is cell speed (S), cell persisence (P), timestep (dt), and total time of simulation
#return time series of x and y positions and list of times that correspond to each position
def PRW_langevin(S, P, dt, time):
    #simplify coefficients
    alpha = 1 - (dt/P)
    F = np.sqrt(((S**2)*(dt**3))/P)
    #initialize time
    t = 0

    #initialize x and y positions
    pos_x = 0
    pos_y = 0 

    #create lists to store cell x and y positions and add initial positions
    all_pos_x = [pos_x]
    all_pos_y = [pos_y]

    #initialize previous step directions/sizes
    dx_prev = np.random.uniform(-1,1)
    dy_prev = np.random.uniform(-1,1)

    #create lists to store cell dx and dy steps and add initial steps
    all_dx = [dx_prev]
    all_dy = [dy_prev]

    #create list to store theta values and add initial value
    all_theta = [np.arctan2(dy_prev,dx_prev)]


    #loop through time and get x and y positions and theta for each time step
    while t<time:
        dx = alpha * dx_prev + F * np.random.normal()
        dy = alpha * dy_prev + F * np.random.normal()
        pos_x += dx
        pos_y += dy
        theta = np.arctan2(pos_y,pos_x)
        dx_prev = dx
        dy_prev = dy
        t += dt
        all_pos_x.append(pos_x)
        all_pos_y.append(pos_y)
        all_dx.append(dx)
        all_dy.append(dy)
        all_theta.append(theta)
        
    return all_pos_x, all_pos_y, all_dx, all_dy, theta

#function to run Langevin PRW simulation for many cells
#input is the number of cells to simulate (Nwalkers), time step (dt), total time of simulation (t),
#and standard deviation for theta angle
#returns a list of dataframes, where each cell has its own dataframe that contains x and y positions,
#x and y components of velocity, and theta
def run_PRW_langevin_sim(Nwalkers, dt, time, S, P):
    #initialize list to store dataframes for each cell
    data_sim = []
    #loop through each cell and run simulation
    for i in range(0,Nwalkers):
        x,y,dx,dy,theta = PRW_langevin(S, P, dt, time)
        x_vel = calc_vel(x)
        y_vel = calc_vel(y)
        vel = np.sqrt(x_vel**2 + y_vel**2)
        r = np.sqrt(np.array(x)**2 + np.array(y)**2)
        stepsize = np.sqrt((np.array(x[0:-1])-np.array(x[1:]))**2 +(np.array(y[0:-1])-np.array(y[1:]))**2)
        net_disp = abs(r[-1] - r[0])
        #tot_path_len = np.sum(r)
        onewalker = {'x': x, 'y': y, 'dx': dx, 'dy':dy, 'vx': x_vel, 'vy': y_vel, 'v': vel, 'theta': theta, 'DoverT': net_disp/np.sum(stepsize)}
        onewalker_df = pd.DataFrame(data=onewalker)
        data_sim.append(onewalker_df)
    return data_sim
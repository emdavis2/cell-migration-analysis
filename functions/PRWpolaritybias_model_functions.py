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

#used to progress internal polarity bias angle, omega, at each time step
#input is omega from previous time step and the standard deviation of the normal distriubution from which we add noise to omega
#output is new omega
def update_omega(omega_prev, std_dev):
    return omega_prev + np.random.normal(0, std_dev)

#run simulation for one cell over time
#input is initial omega value, theta value, time step size, total time length of simulation, cell stepsize (v), and 
#corresponding standard deviations for omega and theta
#returns cell x and y positions, internal polarity bias x and y positions, omega and theta angles as time series data
def PRW_polaritybias(w_prev, theta_prev, dt, time, v, std_dev_w, std_dev_theta):
    #initialize time
    t = 0
    #initialize x and y position of cell
    x_pos = 0
    y_pos = 0
    
    #initialize internal polarity x and y position
    pol_x = np.cos(w_prev)
    pol_y = np.sin(w_prev)

    #create lists to store internal polarity bias angle (omega) and cell direction of movement (theta) values over time
    #initialize omega and theta
    omega = [w_prev]
    theta = [theta_prev]

    #create list to store x and y positions of cell over time
    #add initial x and y positions to lists
    all_x_pos = [x_pos]
    all_y_pos = [y_pos]

    #create list to store internal polarity bias positions over time
    #add initial polarity bias x and y positions to lists
    all_pol_x = [pol_x]
    all_pol_y = [pol_y]

    #progress time 
    t += dt

    #begin loop to update angles and positions over time for each time step 
    while t < time:
        #update internal polarity bias angle, omega, with standard deviation corresponding to omega
        w_next = update_omega(w_prev, std_dev_w)
        #update overall cell position angle such that the new angle is chosen from a normal distribution
        #that is centered on the internal polarity bias angle with standard deviation corresponding to theta
        theta_next = np.random.normal(w_next, std_dev_theta)

        #update x and y positions according to new theta values
        #(stepsize (v) is constant here and is an argument to this function)
        x_pos += v*np.cos(theta_next)
        y_pos += v*np.sin(theta_next)

        #update x and y internal polairty positions
        pol_x += v*np.cos(w_next)
        pol_y += v*np.sin(w_next)

        #add updated values to corresponding lists
        omega.append(w_next)
        theta.append(theta_next)

        all_x_pos.append(x_pos)
        all_y_pos.append(y_pos)

        all_pol_x.append(pol_x)
        all_pol_y.append(pol_y)

        #set up for next time step
        w_prev = w_next
        theta_prev = theta_next

        #progress time
        t += dt

    return all_x_pos, all_y_pos, all_pol_x, all_pol_y, omega, theta

#function to run PRW_polaritybias simulation for many cells
#input is the number of cells to simulate (Nwalkers), time step (dt), total time of simulation (t),
#and standard deviation omega and standard deviation theta
#returns a list of dataframes, where each cell has its own dataframe that contains x and y positions, x and y internal polarity
#positions, x and y components of velocity, and omega and theta
def run_PRWpolaritybias_sim(Nwalkers, dt, time, std_dev_w, std_dev_theta):
    #initialize list to store dataframes for each cell
    data_sim = []
    #loop through each cell and run simulation
    for i in range(0,Nwalkers):
        omega_prev = np.random.uniform(0,2*np.pi)
        theta_prev = np.random.uniform(0,2*np.pi)
        v = np.random.exponential(3)
        x,y,x_pol,y_pol,omega,theta = PRW_polaritybias(omega_prev, theta_prev, dt, time, v, std_dev_w, std_dev_theta)
        x_vel = calc_vel(x)
        y_vel = calc_vel(y)
        vel = np.sqrt(x_vel**2 + y_vel**2)
        stepsize = np.sqrt((np.array(x[0:-1])-np.array(x[1:]))**2 +(np.array(y[0:-1])-np.array(y[1:]))**2)
        r = np.sqrt(np.array(x)**2 + np.array(y)**2)
        net_disp = abs(r[-1] - r[0])
        #tot_path_len = np.sum(r)
        onewalker = {'x': x, 'y': y, 'x_pol':x_pol, 'y_pol': y_pol, 'vx': x_vel, 'vy': y_vel, 'v': vel, 'omega': omega, 'theta': theta, 'DoverT': net_disp/np.sum(stepsize)}
        onewalker_df = pd.DataFrame(data=onewalker)
        data_sim.append(onewalker_df)
    return data_sim
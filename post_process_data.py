# -*- coding: utf-8 -*-
"""
Spyder Editor

"""

# -------------------------------------------------------
from sys import stderr
from numpy import pi
import numpy as np
import scipy.io
import matplotlib.pyplot as plt


def derivative(f,a,method='central',h=0.01, order=1):
    '''Compute the difference formula for f'(a) with step size h.

    Parameters
    ----------
    f : function
        Vectorized function of one variable
    a : number
        Compute derivative at x = a
    method : string
        Difference formula: 'forward', 'backward' or 'central'
    h : number
        Step size in difference formula

    Returns
    -------
    float
        Difference formula:
            central: f(a+h) - f(a-h))/2h
            forward: f(a+h) - f(a))/h
            backward: f(a) - f(a-h))/h            
    '''
    fa = f[a]
    
    for o in range(1,order+1):
    
        fa__h = 0;
        if a - o < 0: # for initialization purposes
            fa__h += f[a]
        else:
            fa__h += f[a - o]
        
        fa_h = 0;
        if a + o >= len(f): # for initialization purposes
            fa_h += f[a]
        else:
            fa_h += f[a + o]
        
    
    if a - 2 < 0: # for initialization purposes
        fa__2h = f[a]
    else:
        fa__2h = f[a - 2]
    
    if a + 2 >= len(f): # for end purposes
        fa_2h = f[a]
    else:
        fa_2h = f[a + 2]
    
    if method == 'central':
        return (fa_h - fa__h)/(2*h*o), (fa_h - 2*fa + fa__h)/(h**2)
    elif method == 'forward':
        return (fa_h - fa)/(h*o), (fa_2h - 2*fa_h + fa)/(h**2)
    elif method == 'backward':
        return (fa - fa__h)/(h*o), (fa - 2*fa_h + fa__2h)/(h**2)
    else:
        raise ValueError("Method must be 'central', 'forward' or 'backward'.")

def Differentiate_raw_data(filename, t_p, x_p, y_p, z_p, display_out):
    resultsfile=open(filename,'w') 
    resultsfile.close()
    
    # differentiate parsed data
    i = 0; vx = []; vy = []; vz = []; ax = []; ay = []; az = [];
    for t, x, y, z in zip(t_p, x_p, y_p, z_p):
        
        #if i == 0:
        #    pos_prev = np.array([x_p[i], y_p[i], z_p[i]]); t_prev = 0.0;
        #    vel_prev = np.array( [0,0,0] ); 
        #else:
        #    pos_prev = np.array([x_p[i-1], y_p[i-1], z_p[i-1]])
        #    vel_prev = velocity[i-1]
        #    t_prev = t_p[i-1];
        #    
        #pos = np.array([x, y, z]);
        #
        #v = (pos - pos_prev) / (0.01)
        #
        #vx += [v[0]]; vy += [v[1]]; vz += [v[2]];
        #velocity += [v]
        #
        #if i == 0:
        #    a = [0.0] * 3;
        #else:
        #    a = (v - vel_prev) / (t - t_prev)
        #
        #acceleration += [a[0]]
        
        [vx_d,ax_d] = derivative(x_p,i,method='central',h=0.01, order=20)
        [vy_d,ay_d] = derivative(y_p,i,method='central',h=0.01, order=20)
        [vz_d,az_d] = derivative(z_p,i,method='central',h=0.01, order=20)
        
        vx += [vx_d]; vy += [vy_d]; vz += [vz_d];
        ax += [ax_d]; ay += [ay_d]; az += [az_d];
        
        # SAVE TO FILE
      
        data = [t, x, y, z, vx_d, vy_d, vz_d, ax_d, ay_d, az_d]
        data = ['{0:.8e}'.format(flt) for flt in data]
        results = " ".join(map(str,data))
        
        resultsfile=open(filename,'a')
        resultsfile.write(' ' + results + '\n')
        i += 1;
    
    # Write data to file
    resultsfile=open(filename,'r')
    lines = resultsfile.readlines()
    resultsfile.close()
    
    # REMOVE LAST EMPTY LINE VERY IMP for C++
    w = open(filename,'w')
    lines[-1] = lines[-1][0:-1]
    w.writelines([item for item in lines])
    w.close()
    # !!
    
    diff_data = [vx, vy, vz]
    
    if display_out:
    
        # Plot down sampled Raw data
        my_dpi = 100
        plt.figure(figsize=(1200/my_dpi, 500/my_dpi), dpi=my_dpi);
        plt.title("x y z comparision");
        plt.subplot(1, 3, 1);
        plt.plot(t_p, vx, "r", label = "Groundtruth");
        
        plt.title("vx");
        plt.legend();
        
        plt.subplot(1, 3, 2);
        plt.plot(t_p, vy, "r", label = "Groundtruth");
        
        plt.title("vy");
        plt.legend();
        
        plt.subplot(1, 3, 3);
        plt.plot(t_p, vz, "r", label = "Groundtruth");
        
        plt.title("vz");
        plt.legend();
    
    return diff_data 
    

def load_raw_data(filename):
    
    mat = scipy.io.loadmat(filename) # get optitrack data
    data = mat['optitack_imu_data']
    data = data[:,0::]
    
    [t_raw_a, t_a, x_a, y_a, z_a, roll_a, pitch_a, yaw_a, ax_a, ay_a, az_a, gx_a, gy_a, gz_a] = data
    
    ax_out = []; ay_out = []; az_out = [];
    mx_out = []; my_out = []; mz_out = [];
    
    g = 9.81;
    
    for t_raw, t, x, y, z, roll, pitch, yaw, ax, ay, az, gx, gy, gz in zip(t_raw_a, t_a, x_a, y_a, z_a, roll_a, pitch_a, yaw_a, ax_a, ay_a, az_a, gx_a, gy_a, gz_a):
    
        ax_out += [ax/g]; ay_out += [ay/g]; az_out += [az/g];
        mx_out += [0.0]; my_out += [0.0]; mz_out += [0.0];

    ## FLUSH TO STERR
    #stderr.write("""\r t = %s s |  (ax,ay,az) = (%f,%f,%f) m/s^2
    #                \r t = %s s |  (gx,gy,gz) = (%f,%f,%f) rad/s
    #                \r t = %s s |  (mx,my,mz) = (%f,%f,%f) m/s^2
    #                \r t = %s s |  (roll,pitch,yaw) = (%f,%f,%f) deg"""
    #             % (t, ax, ay, az, t, gx, gy, gz, t, mx, my, mz, t, roll * (180.0 / pi), pitch * (180.0 / pi), yaw * (180.0 / pi)) )
    #stderr.flush()
    
    data_out = [t_raw_a, t_a, x_a, y_a, z_a, roll_a, pitch_a, yaw_a, ax_out, ay_out, az_out, gx_a, gy_a, gz_a, mx_out, my_out, mz_out]
    
    return data_out

def parse_data(Raw_position_data, display_result):
    # loop to parse optitrack data
    
    [t_raw_a, x_a, y_a, z_a, roll_a, pitch_a, yaw_a] = Raw_position_data
    
    x_prev = 100.0; # dummy index to initialize first value
    
    j = 0; x_p = []; y_p = []; z_p = []; t_p = [];
    for t, x, y, z in zip( t_raw_a, x_a, y_a, z_a ):
        
        if x != x_prev:
            x_p += [x]
            y_p += [y]
            z_p += [z]
            t_p += [t]
            
        x_prev = x
    
    if display_result:
        # Plot parsed Raw data
        my_dpi = 100
        plt.figure(figsize=(1200/my_dpi, 500/my_dpi), dpi=my_dpi);
        plt.title("x y z comparision");
        plt.subplot(1, 3, 1);
        plt.plot(t_p, x_p, "k", label = "Mahony");
        #plt.plot("Groundtruth", t_a, x_a, "r");
        
        plt.title("x");
        plt.legend();
        
        plt.subplot(1, 3, 2);
        plt.plot(t_p, y_p, "k", label = "Mahony");
        #plt.plot("Groundtruth", t_a, x_a, "r");
        
        plt.title("y");
        plt.legend();
        
        plt.subplot(1, 3, 3);
        plt.plot(t_p, z_p, "k", label = "Mahony");
        #plt.plot("Groundtruth", t_a, x_a, "r");
        
        plt.title("z");
        plt.legend();
    
    Parsed_data = [t_p, x_p, y_p, z_p]
    return Parsed_data

def unparse_data(Raw_position_data, Parsed_data, display_result):
    # loop to parse optitrack data
    
    [t_raw_a, x_a, y_a, z_a, roll_a, pitch_a, yaw_a] = Raw_position_data
    [vx_p, vy_p, vz_p] = Parsed_data
    x_prev = 100.0; # dummy index to initialize first value
    
    j = 0; vx_a = []; vy_a = []; vz_a = [];
    for t, x, y, z in zip( t_raw_a, x_a, y_a, z_a ):
        
        if x != x_prev:
            vx_hold = vx_p[j]
            vy_hold = vy_p[j]
            vz_hold = vz_p[j]
            
            j += 1
            
        vx_a += [vx_hold]
        vy_a += [vy_hold]
        vz_a += [vz_hold]
            
        x_prev = x
    
    if display_result:
        # Plot parsed Raw data
        my_dpi = 100
        plt.figure(figsize=(1200/my_dpi, 500/my_dpi), dpi=my_dpi);
        plt.title("x y z comparision");
        plt.subplot(1, 3, 1);
        plt.plot(t_raw_a, vx_a, "k", label = "Mahony");
        #plt.plot("Groundtruth", t_a, x_a, "r");
        
        plt.title("vx");
        plt.legend();
        
        plt.subplot(1, 3, 2);
        plt.plot(t_raw_a, vy_a, "k", label = "Mahony");
        #plt.plot("Groundtruth", t_a, x_a, "r");
        
        plt.title("vy");
        plt.legend();
        
        plt.subplot(1, 3, 3);
        plt.plot(t_raw_a, vz_a, "k", label = "Mahony");
        #plt.plot("Groundtruth", t_a, x_a, "r");
        
        plt.title("vz");
        plt.legend();
    
    unparsed_data = [vx_a, vy_a, vz_a]
    return unparsed_data

def delay_signal(delay,Raw_position_data):
    # delay the downsampled signal
    [t_raw_a, x_a, y_a, z_a, roll_a, pitch_a, yaw_a] = Raw_position_data
    
    i = 0; 
    x_delay = []; y_delay = []; z_delay = []; roll_delay = []; pitch_delay = []; yaw_delay = [];
    for n in range(len(t_raw_a)):
        
        if (i-delay) < 0:
            roll_delay += [roll_a[0]]; pitch_delay += [pitch_a[0]]; yaw_delay += [yaw_a[0]]
            x_delay += [x_a[0]]; y_delay += [y_a[0]]; z_delay += [z_a[0]]; 
        else:
            roll_delay += [roll_a[i-delay]]; pitch_delay += [pitch_a[i-delay]]; yaw_delay += [yaw_a[i-delay]]
            x_delay += [x_a[i-delay]]; y_delay += [y_a[i-delay]]; z_delay += [z_a[i-delay]]; 
        
        i += 1;
    
    delayed_data = [roll_delay,pitch_delay,yaw_delay,x_delay,y_delay,z_delay]
        
    return delayed_data

def downsample_data(downsample_frequency,delayed_data):
        
    [roll_delay,pitch_delay,yaw_delay,x_delay,y_delay,z_delay] = delayed_data
    
    
    i = 0;
    x_ds = []; y_ds = []; z_ds = []; roll_ds = []; pitch_ds = []; yaw_ds =[];
    for x, y, z, roll, pitch, yaw in zip(x_delay, y_delay, z_delay, roll_delay, pitch_delay, yaw_delay):
        
        if (i) % downsample_frequency == 0: # downsample the motion capture data
            [roll_hold , pitch_hold, yaw_hold, x_hold, y_hold, z_hold] = [roll, pitch, yaw, x, y, z];
        
        roll_ds += [roll_hold]; pitch_ds += [pitch_hold]; yaw_ds += [yaw_hold]
        x_ds += [x_hold]; y_ds += [y_hold]; z_ds += [z_hold];
        i += 1;
        
    ds_data = [roll_ds,pitch_ds,yaw_ds,x_ds,y_ds,z_ds]
        
    return ds_data

def export_modified_data(filename,IMU_data,Raw_position_data,Raw_velocity_data,delayed_data,ds_data,display_out):
    
    [t_a, ax_a, ay_a, az_a, gx_a, gy_a, gz_a, mx_a, my_a, mz_a] = IMU_data
    [t_raw_a, x_a, y_a, z_a, roll_a, pitch_a, yaw_a] = Raw_position_data
    [vx_a, vy_a, vz_a] = Raw_velocity_data
    [roll_delay,pitch_delay,yaw_delay,x_delay,y_delay,z_delay] = delayed_data
    [roll_ds,pitch_ds,yaw_ds,x_ds,y_ds,z_ds] = ds_data
    
    data_write_out = [ax_a, ay_a, az_a, gx_a, gy_a, gz_a, mx_a, my_a, mz_a, roll_ds, pitch_ds, yaw_ds, x_ds, y_ds, z_ds, roll_a, pitch_a, yaw_a, x_a, y_a, z_a, vx_a, vy_a, vz_a, t_a]
    data_write_out = list(map(list, zip(*data_write_out))) # transpose a list of lists
    
    with open(filename, "w", newline=None) as f:
        np.savetxt(filename, data_write_out, delimiter=" ", fmt = "%.8e", newline='\n')
    
    # read data from file
    resultsfile=open(filename,'r')
    lines = resultsfile.readlines()
    
    # REMOVE LAST EMPTY LINE VERY IMP for C++
    resultsfile.close()
    w = open(filename,'w')
    
    lines[-1] = lines[-1][0:-1]
    w.writelines([item for item in lines])
    
    w.close()
    # !!!!
    
    if display_out:
    
        # Plot down sampled Raw data
        my_dpi = 100
        plt.figure(figsize=(1200/my_dpi, 500/my_dpi), dpi=my_dpi);
        plt.title("x y z comparision");
        plt.subplot(1, 3, 1);
        plt.plot(t_a, x_ds, "k", label = "downsample");
        plt.plot(t_a, x_a, "r", label = "Groundtruth");
        
        plt.title("x");
        plt.legend();
        
        plt.subplot(1, 3, 2);
        plt.plot(t_a, y_ds, "k", label = "Mahony");
        plt.plot(t_a, y_a, "r", label = "Groundtruth");
        
        plt.title("y");
        plt.legend();
        
        plt.subplot(1, 3, 3);
        plt.plot(t_a, z_ds, "k", label = "Mahony");
        plt.plot(t_a, z_a, "r", label = "Groundtruth");
        
        plt.title("z");
        plt.legend();

def main():
    
    plt.close('all') # close all existing figures
    
    filename = 'optitrack_Data.mat'
    [t_raw_a, t_a, x_a, y_a, z_a, roll_a, pitch_a, yaw_a, ax_a, ay_a, az_a, gx_a, gy_a, gz_a, mx_a, my_a, mz_a] = load_raw_data(filename)
    
    IMU_data = [t_a, ax_a, ay_a, az_a, gx_a, gy_a, gz_a, mx_a, my_a, mz_a]
    Raw_position_data = [t_raw_a, x_a, y_a, z_a, roll_a, pitch_a, yaw_a]
    
    # parse data for differentiation
    Parsed_data = parse_data(Raw_position_data, False)
    [t_p, x_p, y_p, z_p] = Parsed_data
    
    # differentiate data
    filename = "NAV3_optitrack_data.bin"
    diff_data = Differentiate_raw_data(filename, t_p, x_p, y_p, z_p, False)
    
    Raw_velocity_data = unparse_data(Raw_position_data, diff_data, True)
    
    # delay raw signal
    delay = 20; #samples
    delayed_data = delay_signal(delay,Raw_position_data)
    
    # Downsample signal
    # 1000 = 1hz
    # 200 = 5hz
    # 40 = 25 Hz
    downsample_frequency_units = 40;
    ds_data = downsample_data(downsample_frequency_units,delayed_data)
    
    # export modified data to file
    filename = "NAV3_data.bin"
    export_modified_data(filename,IMU_data,Raw_position_data,Raw_velocity_data,delayed_data,ds_data,True)
    filename = "NAV3_data_ds_%i_delay_%i.bin" %(downsample_frequency_units,delay)
    export_modified_data(filename,IMU_data,Raw_position_data,Raw_velocity_data,delayed_data,ds_data,False)
    
if __name__ == "__main__":
    main()

#-------------------------------------------------------
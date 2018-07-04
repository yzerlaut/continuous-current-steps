import sys, os
sys.path.append('../../')
from data_analysis.IO import hdf5 # download data_analysis module at https://bitbucket.org/yzerlaut/data_analysis
import numpy as np
from graphs.my_graph import *  # download graphs module at https://bitbucket.org/yzerlaut/graphs

##################################################################
######## FETCHING DATA ###########################################
##################################################################

def sort_keys(data):
    keys = list(data.keys())
    for key in keys:
        if len(key.split('Iout'))>1:
            data['Iout'] = key
        if len(key.split('Vm'))>1:
            data['Vm'] = key
        if len(key.split('LFP'))>1:
            data['LFP'] = key

def get_Iout_Flag_for_Inactive():
    Iout_Flag_for_Inactive = 0
    with open(os.path.realpath(__file__).replace('.py', '.h')) as f:
        content = f.readlines()
        for line in content:
            sp = line.split('const double Iout_Flag_for_Inactive = ')
            if len(sp)>1:
                Iout_Flag_for_Inactive = float(sp[1].replace(';', ''))
    if Iout_Flag_for_Inactive==0:
        print('---------------------------------------------')
        print('Failed to find the baseline level in the C++ file')
        print('---------------------------------------------')
    return Iout_Flag_for_Inactive
            
def get_onset_stimulation(data):
    Iout0 = data['Iout_Flag_for_Inactive']
    data['protocol_start'] = np.arange(len(data[data['Vm']]))[data[data['Iout']]!=(Iout0*1e-12)][0]*data['dt']

def from_filename_to_time_stamp(filename, extension='.RTXI.h5'):
    new_filename = filename.split(os.path.sep)[-1] # to be sure to have the last extension
    time_stamp = 0
    for val, factor in zip(new_filename.split(extension)[0].split(':'), [3600., 60., 1.]):
        time_stamp += factor*float(val)
    return time_stamp

def find_stimulation_file_and_add_parameters(data):
    """

    """
    file_directory = os.path.dirname(data['filename'])
    time_stamp_filename = from_filename_to_time_stamp(data['filename'], extension='.RTXI.h5')
    file_list= os.listdir(file_directory)
    json_list, time_stamps = [], []
    for ff in file_list:
        if len(ff.split('.csteps.JSON'))>1:
            json_list.append(ff)
            time_stamps.append(from_filename_to_time_stamp(ff, extension='.csteps.JSON'))
    # the right file is the one right before the recording
    ii = np.argmin(np.abs(np.array(time_stamps)-data['protocol_start']-time_stamp_filename))
    param_filename = json_list[ii]
    with open(file_directory+os.path.sep+param_filename, 'r') as fn:
        exec('dd='+fn.read()) # the
    # we loop through the dd dictionary and set its key and values in the data dictionary
    for key in locals()['dd']:
        data[key] = locals()['dd'][key]
    data['start_vector'] = np.array(data['start_vector'])+data['protocol_start']
    data['stop_vector'] = np.array(data['stop_vector'])+data['protocol_start']
    data['amplitude_vector'] = np.array(data['amplitude_vector'])
        
def load_data(filename):
    
    # load HDF5 (see data_analysis.IO module)
    data = hdf5.load_continous_RTXI_recording(filename)
    data['filename'] = filename
    
    # protocol-dependent processing of the input
    sort_keys(data)
    data['Iout_Flag_for_Inactive'] = get_Iout_Flag_for_Inactive()
    get_onset_stimulation(data)
    find_stimulation_file_and_add_parameters(data)
    data['istim_protocol_start'] = np.arange(len(data['start_vector']))[data['start_vector']>data['protocol_start']][0]
    print(data['protocol_start'], data['start_vector'])
    # other...
    data['t'] = np.arange(len(data[data['Vm']]))*data['dt']
    
    return data
    
##################################################################
######## ANALYSIS STARTS HERE ####################################
##################################################################

def get_response_patterns(data,
                          muV_levels=np.linspace(-70e-3, -50e-3, 3),
                          pre_window = 100e-3,
                          spike_threshold = -50e-3,
                          muV_discret=3):

    # levels of amplitudes
    amp_levels = np.unique(data['amplitude_vector'])
    
    # levels of muV
    if muV_levels is None:
        muV_levels_min = np.histogram(data[data['Vm']], bins=20)[1][1]
        muV_levels_max = np.histogram(data[data['Vm']], bins=20)[1][-2]
        muV_levels = np.linspace(muV_levels_min, muV_levels_max, muV_discret+1)
        
    # time window
    t_window = np.arange(int((data['stop_vector'][0]-data['start_vector'][0]+2*pre_window)/data['dt']))*data['dt']-pre_window
    
    RESPONSES = [[ [] for i in range(len(amp_levels))] for j in range(len(muV_levels)-1)]
    # i.e. will corespond to RESPONSES[muV_level][amp]
        
    i=data['istim_protocol_start']
    while (data['stop_vector'][i]<data['t'][-1]) and (data[data['Iout']][int(data['start_vector'][i]/data['dt'])]!=(data['Iout_Flag_for_Inactive']*1e-12)):

        cond_pre = (data['t']>data['start_vector'][i]-pre_window) &\
                   (data['t']<data['start_vector'][i])
        muV_before = np.mean(data[data['Vm']][cond_pre])
        # find the pre-depolarization level of this episode
        j = np.argwhere((muV_levels[:-1]<muV_before) & (muV_levels[1:]>=muV_before))
        # find the current amplitude of this episode
        k = np.argwhere(amp_levels==data['amplitude_vector'][i])
        # print(k, j, muV_before)
        if (len(k)>0) and (len(j)>0):
            cond = (data['t']>(data['start_vector'][i]-pre_window)) &\
                   (data['t']<=(data['stop_vector'][i]+pre_window))
            vec = 0*t_window+data[data['Vm']][cond][-1]
            vec[:np.min([len(vec), len(data[data['Vm']][cond])])] = data[data['Vm']][cond][:np.min([len(vec), len(data[data['Vm']][cond])])]
            RESPONSES[j[0][0]][k[0][0]].append(vec)
            RESPONSES[j[0][0]][k[0][0]][-1][RESPONSES[j[0][0]][k[0][0]][-1]>spike_threshold] = spike_threshold
        else:
            print('muV_level not classified: ', 1e3*muV_before)
        i+=1

    return t_window, RESPONSES, amp_levels, muV_levels

    
##################################################################
######## PLOTTING STARTS HERE ####################################
##################################################################

def make_raw_data_figure(data,
                         args,
                         figsize=(.8,.13),
                         Vm_color=Blue,
                         Iinj_color=Orange,
                         Vpeak = -10e-3,
                         Tbar = 1, Ibar = 100,
                         Vm_enhancement_factor=3.,
                         spike_ms=4.):
    """

    """
    fig, ax = figure(figsize=figsize, left=.1, bottom=.1)
    # time conditions
    cond = (data['t']>args.tzoom[0]) & (data['t']<args.tzoom[1])
    # from Iinj to Vm
    ax.plot(data['t'][cond], data[data['Iout']][cond], color=Iinj_color, lw=1)
    Imin, Imax = np.min(data[data['Iout']][cond]), np.max(data[data['Iout']][cond])
    Vmin, Vmax = np.min(data[data['Vm']][cond]), np.max(data[data['Vm']][cond])
    ispikes = np.argwhere((data[data['Vm']][cond][1:]>Vpeak) & (data[data['Vm']][cond][:-1]<=Vpeak)).flatten()
    for ii in ispikes:
        ax.plot([data['t'][ii]], [(Vpeak-Vmin)*Vm_enhancement_factor/(Vmax-Vmin)*(Imax-Imin)+Imax], '*',  color=Vm_color, ms=spike_ms)
    data[data['Vm']][data[data['Vm']]>Vpeak] = Vpeak
    ax.plot(data['t'][cond], Imax+(data[data['Vm']][cond]-Vmin)*Vm_enhancement_factor/(Vmax-Vmin)*(Imax-Imin), color=Vm_color, lw=1)

    if args.debug:
        for k in range(data['istim_protocol_start'], 25):
            ax.plot([data['start_vector'][k]], [1e-12*data['amplitude_vector'][k]], 'k*')
    set_plot(ax, [], xlim=[data['t'][cond][0], data['t'][cond][-1]])
    ax.annotate('$I_{inj}$', (0, Imin), color=Iinj_color)
    ax.annotate('$V_{m}$', (0, Imax), color=Vm_color)
    ax.plot([1,1+Tbar], [Imin, Imin], 'k-', lw=2)
    ax.plot([1,1], [Imin, Imin+Ibar*1e-12], 'k-', lw=2)
    ax.annotate(str(int(Tbar))+'s', (1.2, Imin), color='k')
    ax.annotate(str(Ibar)+'pA', (1., Imin+Ibar*1e-12), color=Iinj_color, rotation=90)
    ax.annotate(str(round(1e3*Ibar*1e-12/Vm_enhancement_factor*(Vmax-Vmin)/(Imax-Imin),1))+'mV',\
                (0., Imin+Ibar*1e-12), color=Vm_color, rotation=90)
    return fig, ax


def make_trial_average_figure(data, args,
                              ylim=[-75,-45]):
    """

    """

    t_window, RESPONSES, amp_levels, muV_levels = get_response_patterns(data,
                                                   muV_levels=1e-3*np.array(args.muV_levels))

    fig, AX = figure(figsize=(.2*len(amp_levels),.2),
                     axes=(1, len(amp_levels)),
                     left=0.1/len(amp_levels), top=.8)
    
    for a in range(len(amp_levels)):
        for m in range(len(muV_levels)-1):
            RESP = []
            for resp in RESPONSES[m][a]:
                RESP.append(resp)
            AX[0][a].plot(1e3*t_window, 1e3*np.mean(RESP, axis=0),
                          '-', lw=0.5, color=viridis(m/(len(muV_levels)-1)))
        AX[0][a].plot([0,0], ylim, 'w.')
        AX[0][a].set_title('$I_{inj}$='+str(amp_levels[a])+'pA')
        set_plot(AX[0][a], xlabel='time from stim. (ms)')
            
    return fig, AX


if __name__ == '__main__':

    import matplotlib.pylab as plt


    import argparse
    # First a nice documentation 
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--filename", '-f', help="filename",type=str,
                        default='/home/yann/DATA/2018_07_04/15:57:47.RTXI.h5')
    parser.add_argument("--mean",help="", type=float, default=5.)
    parser.add_argument("--tzoom",help="", type=float, nargs=2, default=[0,20])
    parser.add_argument("--muV_levels", help="", type=float, nargs='*', default=[-70, -60, -50])
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument("-d", "--debug", help="debug option", action="store_true")
    parser.add_argument("-s", "--save", help="save the figures", action="store_true")
    parser.add_argument("-ta", "--trial_average_analysis", action="store_true")
    args = parser.parse_args()


    if not os.path.isfile(args.filename):
        print('---------------------------------------------------------')
        print('you should provide a hdf5 file as a "--filename" argument')
        print('---------------------------------------------------------')
    else:
        data = load_data(args.filename)
        print(data.keys())
        if args.trial_average_analysis:
            make_trial_average_figure(data, args)
            # get_response_patterns(data)
        else:
            make_raw_data_figure(data, args)
        show()

        


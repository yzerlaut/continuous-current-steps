import sys, os
sys.path.append('../../')
from data_analysis.IO import rtxi # download data_analysis module at https://bitbucket.org/yzerlaut/data_analysis
import numpy as np
from graphs.my_graph import *  # download graphs module at https://bitbucket.org/yzerlaut/graphs
from data_analysis.manipulation.files import * # download data_analysis module at https://bitbucket.org/yzerlaut/
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
    data = rtxi.load_continous_RTXI_recording(filename)
    data['filename'] = filename
    
    # protocol-dependent processing of the input
    sort_keys(data)
    data['Iout_Flag_for_Inactive'] = get_Iout_Flag_for_Inactive()
    get_onset_stimulation(data)
    find_stimulation_file_and_add_parameters(data)
    data['istim_protocol_start'] = np.arange(len(data['start_vector']))[data['start_vector']>data['protocol_start']][0]
    
    # other...
    data['t'] = np.arange(len(data[data['Vm']]))*data['dt']
    # stim duration is constant in this protocol
    data['stim_duration'] = data['stop_vector'][0]-data['start_vector'][0]
    
    return data
    
##################################################################
######## ANALYSIS STARTS HERE ####################################
##################################################################

def quick_conductance_estimate(data, args):
    
    # levels of amplitudes
    amp_levels = np.unique(data['amplitude_vector'])
    data['amp_levels'] = amp_levels
        
    # time window
    pre_window = data['stim_duration']
    t_window = np.arange(int((data['stop_vector'][0]-data['start_vector'][0]+2*pre_window)/data['dt']))*data['dt']-pre_window
    discard_time = 1e-3*args.discard_window_for_Vm_rise

    PRE_MUV, DEFLECT_MUV, I_AMP = [], [], []

    i=data['istim_protocol_start']
    while (data['stop_vector'][i]<data['t'][-1]) and (data[data['Iout']][int(data['start_vector'][i]/data['dt'])]!=(data['Iout_Flag_for_Inactive']*1e-12)):

        # only if above a given 
        if np.abs(data['amplitude_vector'][i])>=args.min_I_for_quick_conductance_estimate:
            
            # pre-stimulus condition
            cond_pre = (data['t']>data['start_vector'][i]-data['stim_duration']+discard_time) &\
                       (data['t']<data['start_vector'][i])
            muV_before = np.mean(data[data['Vm']][cond_pre])
            # post-stim cond
            cond_post = (data['t']>data['start_vector'][i]+discard_time) &\
                        (data['t']<data['start_vector'][i]+data['stim_duration'])
            muV_after = np.mean(data[data['Vm']][cond_post])
            PRE_MUV.append(muV_before)
            DEFLECT_MUV.append(muV_after-muV_before)
            I_AMP.append(data['amplitude_vector'][i])
        i+=1

    PRE_MUV, DEFLECT_MUV, I_AMP = np.array(PRE_MUV), np.array(DEFLECT_MUV), np.array(I_AMP)
    imuV_sorted = np.argsort(PRE_MUV)[:args.N_for_quick_conductance_estimate]
    conductance_estimate = round(np.mean(I_AMP[imuV_sorted])/1e3/np.mean(DEFLECT_MUV[imuV_sorted]),2)
    print('---------------------------')
    print('conductance estimate: ', conductance_estimate)
    f = open(os.path.join(os.getenv("HOME"), 'DATA', 'cell_conductance.txt'), 'w')
    f.write(str(conductance_estimate))
    f.close()
    print('writen in ~/DATA/cell_conductance.txt')
    print('---------------------------')

    
def get_response_patterns(data, args,
                          spike_threshold = -50e-3,
                          muV_discret=3):
    
    # levels of amplitudes
    amp_levels = np.unique(data['amplitude_vector'])
    
    # levels of muV
    muV_levels=1e-3*np.array(args.muV_levels)    
    data['amp_levels'] = amp_levels
    data['muV_levels'] = muV_levels
        
    # time window
    pre_window = data['stim_duration'] # np.max([1e-3*args.pre_window, data['stim_duration']]) # at least stim duration, switched to s !!
    t_window = np.arange(int((data['stop_vector'][0]-data['start_vector'][0]+2*pre_window)/data['dt']))*data['dt']-pre_window
    
    RESPONSES = [[ [] for i in range(len(amp_levels))] for j in range(len(muV_levels)-1)]
    PRE_MUV = [[ [] for i in range(len(amp_levels))] for j in range(len(muV_levels)-1)]
    POST_MUV = [[ [] for i in range(len(amp_levels))] for j in range(len(muV_levels)-1)]
    DEFLECT_MUV = [[ [] for i in range(len(amp_levels))] for j in range(len(muV_levels)-1)]
    CONDUCTANCES = [[ [] for i in range(len(amp_levels))] for j in range(len(muV_levels)-1)]
    
    discard_time = 1e-3*args.discard_window_for_Vm_rise

    # i.e. will corespond to RESPONSES[muV_level][amp]
        
    i=data['istim_protocol_start']
    while (data['stop_vector'][i]<data['t'][-1]) and (data[data['Iout']][int(data['start_vector'][i]/data['dt'])]!=(data['Iout_Flag_for_Inactive']*1e-12)):

        # pre-stimulus condition
        cond_pre = (data['t']>data['start_vector'][i]-data['stim_duration']+discard_time) &\
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
            PRE_MUV[j[0][0]][k[0][0]].append(muV_before)
            
            # post-stimulus condition (rise time discarded)
            cond_post = (data['t']>data['start_vector'][i]+discard_time) &\
                       (data['t']<data['start_vector'][i]+data['stim_duration'])
            POST_MUV[j[0][0]][k[0][0]].append(np.mean(data[data['Vm']][cond_post]))
            DEFLECT_MUV[j[0][0]][k[0][0]].append(POST_MUV[j[0][0]][k[0][0]][-1]-PRE_MUV[j[0][0]][k[0][0]][-1])
        else:
            print('muV_level not classified: ', 1e3*muV_before)
        i+=1

    data['t_window'] = t_window
    data['RESPONSES'] = RESPONSES
    data['PRE_MUV'] = PRE_MUV
    data['POST_MUV'] = POST_MUV
    data['DEFLECT_MUV'] = DEFLECT_MUV
    data['CONDUCTANCES'] = CONDUCTANCES


##################################################################
######## PLOTTING STARTS HERE ####################################
##################################################################

def make_raw_data_figure(data,
                         args,
                         figsize=(.8,.13),
                         Vm_color=Blue,
                         Iinj_color=Orange,
                         Vpeak = -10e-3,
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
        ax.plot([data['t'][cond][ii]], [(Vpeak-Vmin)*Vm_enhancement_factor/(Vmax-Vmin)*(Imax-Imin)+Imax], '*',  color=Vm_color, ms=spike_ms)
    data[data['Vm']][data[data['Vm']]>Vpeak] = Vpeak
    ax.plot(data['t'][cond], Imax+(data[data['Vm']][cond]-Vmin)*Vm_enhancement_factor/(Vmax-Vmin)*(Imax-Imin), color=Vm_color, lw=1)

    if args.debug:
        for k in range(data['istim_protocol_start'], 25):
            ax.plot([data['start_vector'][k]], [1e-12*data['amplitude_vector'][k]], 'k*')
            
    set_plot(ax, [], xlim=[data['t'][cond][0], data['t'][cond][-1]])
    ax.annotate('$I_{inj}$', (args.tzoom[0], Imin), color=Iinj_color)
    ax.annotate('$V_{m}$', (args.tzoom[0], Imax), color=Vm_color)
    ax.plot([args.tzoom[0]+1,args.tzoom[0]+1+args.Tbar], [Imin, Imin], 'k-', lw=2)
    ax.plot([args.tzoom[0]+1,args.tzoom[0]+1], [Imin, Imin+args.Ibar*1e-12], 'k-', lw=2)
    ax.annotate(str(int(args.Tbar))+'s', (args.tzoom[0]+1.2, Imin), color='k')
    ax.annotate(str(args.Ibar)+'pA', (args.tzoom[0]+1, Imin+args.Ibar*1e-12), color=Iinj_color, rotation=90)
    ax.annotate(str(round(1e3*args.Ibar*1e-12/Vm_enhancement_factor*(Vmax-Vmin)/(Imax-Imin),1))+'mV',\
                (args.tzoom[0], Imin+args.Ibar*1e-12), color=Vm_color, rotation=90)
    return fig, ax


def make_trial_average_figure(data, args):
    """

    """

    get_response_patterns(data, args)

    if len(data['amp_levels'])==1:
        fig, ax = figure(figsize=(.25,.2), left=0.9, top=.8)
        AX = [[ax]]
    else:
        fig, AX = figure(figsize=(.2*len(data['amp_levels']),.2),
                         axes=(1, len(data['amp_levels'])),
                         wspace=0.2,
                         left=0.25, top=.8)
    
    for a in range(len(data['amp_levels'])):
        for m in range(len(data['muV_levels'])-1):
            RESP = []
            for resp in data['RESPONSES'][m][a]:
                RESP.append(resp)
            if len(RESP)>0:
                AX[0][a].plot(1e3*data['t_window'], 1e3*np.mean(RESP, axis=0),
                          '-', lw=0.5, color=viridis(m/(len(data['muV_levels'])-1)))
            else:
                print('interval :', data['muV_levels'][m], data['muV_levels'][m+1],\
                      'does not contain data')
        AX[0][a].plot([0,0], args.Vm_lim, 'w.')
        AX[0][a].set_title('$I_{inj}$='+str(data['amp_levels'][a])+'pA')
        if (a==0):
            set_plot(AX[0][a], xlabel='time from stim. (ms)', ylabel='Vm (mV)')
        else:
            set_plot(AX[0][a], xlabel='time from stim. (ms)', yticks_labels=[])
            
    return fig, AX

def make_conductance_fig(data, args):
    """

    """
    
    get_response_patterns(data, args)

    if len(data['amp_levels'])==1:
        fig, ax = figure(figsize=(.25,.2), left=0.9, top=.8)
        AX = [[ax]]
    else:
        fig, AX = figure(figsize=(.2*len(data['amp_levels']),.2),
                         axes=(1, len(data['amp_levels'])),
                         left=0.1/len(data['amp_levels']), top=.8)
    
    for a in range(len(data['amp_levels'])):
        for m in range(len(data['muV_levels'])-1):
            RESP = []
            for resp in data['RESPONSES'][m][a]:
                RESP.append(resp)
            AX[0][a].plot(1e3*data['t_window'], 1e3*np.mean(RESP, axis=0),
                          '-', lw=0.5, color=viridis(m/(len(data['muV_levels'])-1)))
        AX[0][a].plot([0,0], args.Vm_lim, 'w.')
        AX[0][a].set_title('$I_{inj}$='+str(data['amp_levels'][a])+'pA')
        if (a==0):
            set_plot(AX[0][a], xlabel='time from stim. (ms)', ylabel='Vm (mV)')
        else:
            set_plot(AX[0][a], xlabel='time from stim. (ms)')
            
    return fig, AX

def make_depol_fig(data, args):
    """

    """
    
    get_response_patterns(data, args)
    
    if len(data['amp_levels'])==1:
        fig, ax = figure(figsize=(.25,.2), left=0.9, top=.8)
        AX = [[ax]]
    else:
        fig, AX = figure(figsize=(.2*len(data['amp_levels']),.2),
                         axes=(1, len(data['amp_levels'])),
                         wspace=0.2,
                         left=0.25, top=.8)

    Iinj = [[] for m in range(len(data['muV_levels'])-1)]
    Depol = [[] for m in range(len(data['muV_levels'])-1)]

    for a in range(len(data['amp_levels'])):
        for m in range(len(data['muV_levels'])-1):
            for pre_muV, depol in zip(data['PRE_MUV'][m][a], data['DEFLECT_MUV'][m][a]):
                AX[0][a].plot([1e3*pre_muV], [1e3*depol], 'ko', ms=args.ms)
                Iinj[m].append(data['amp_levels'][a])
                Depol[m].append(1e3*depol)
                
        AX[0][a].plot(args.Vm_lim, args.Depol_lim, 'w.')
        AX[0][a].set_title('$I_{inj}$='+str(data['amp_levels'][a])+'pA')
        if (a==0):
            set_plot(AX[0][a], xlabel='pre-stim $\mu_V$ (mV)', ylabel='stim-evoked \n depol. (mV)')
        else:
            set_plot(AX[0][a], xlabel='pre-stim $\mu_V$ (mV)')

    fig2, ax2 = figure(figsize=(.25,.2), left=0.9, top=.8)
    for m in range(len(data['muV_levels'])-1):
        x, y, sy = [], [], []
        for a in np.unique(Iinj[m]):
            x.append(a)
            cond = np.array(Iinj[m])==a
            y.append(np.mean(np.array(Depol[m])[cond]))
            sy.append(np.std(np.array(Depol[m])[cond]))
        ax2.plot(x, y, lw=3, color=viridis(m/(len(data['muV_levels'])-1)))
        ax2.fill_between(x, np.array(y)-np.array(sy), np.array(y)+np.array(sy), lw=0, alpha=.2,
                         color=viridis(m/(len(data['muV_levels'])-1)))
        
    return fig, fig2



if __name__ == '__main__':

    import matplotlib.pylab as plt


    import argparse
    # First a nice documentation 
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--filename", '-f', help="filename",type=str, default='')
    parser.add_argument("--mean",help="", type=float, default=5.)

    # analysis setting
    parser.add_argument("-nqc", "--N_for_quick_conductance_estimate",
                        help="",
                        type=int, default=20)
    parser.add_argument('-mIqc', "--min_I_for_quick_conductance_estimate",
                        help="",
                        type=float, default=60)
    parser.add_argument("--pre_window",
                        help="window length to evaluate pre-stimulus values (in ms)",
                        type=float, default=200)
    parser.add_argument("--discard_window_for_Vm_rise",
                        help="window length to discard Vm rise (in ms)",
                        type=float, default=20)
    
    
    # graphs settings
    parser.add_argument("--tzoom",help="", type=float, nargs=2, default=[0,20])
    parser.add_argument("--Vm_lim",help="", type=float, nargs=2, default=[-75,-49])
    parser.add_argument("--Depol_lim",help="", type=float, nargs=2, default=[-13,13])
    parser.add_argument("--Ibar",help="", type=float, default=50)
    parser.add_argument("--Tbar",help="", type=float, default=1)
    parser.add_argument("--muV_levels", help="", type=float, nargs='*', default=[-70, -60, -50])
    parser.add_argument("--ms",help="marker size", type=float, default=4)
    
    # protocol types
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument("-d", "--debug", help="debug option", action="store_true")
    parser.add_argument("-s", "--save", help="save the figures", action="store_true")
    parser.add_argument("-rd", "--raw_data", action="store_true")
    parser.add_argument("-ta", "--trial_average_analysis", action="store_true")
    parser.add_argument("-ca", "--conductance_analysis", action="store_true")
    parser.add_argument("-da", "--depol_analysis", action="store_true")
    args = parser.parse_args()


    if not os.path.isfile(args.filename):
        print('------------------------------------------------------------')
        print('you should provide a hdf5 file as a "--filename (-f)" argument')
        print('------------------------------------------------------------')
        print('as you didnt, you should pick up a file from:')
        last_dir = get_directories_ordered_by_creation(os.path.join(os.getenv("HOME"), 'DATA'))[-1]
        args.filename = choose_a_file_based_on_keyboard_input(last_dir, extension='RTXI.h5', Nmax=5)

    print('[...] loading data')
    data = load_data(args.filename)
    print('[...] analyzing')
    if args.trial_average_analysis:
        make_trial_average_figure(data, args);show()
    elif args.conductance_analysis:
        make_conductance_fig(data, args);show()
    elif args.depol_analysis:
        make_depol_fig(data, args);show()
    elif args.raw_data:
        make_raw_data_figure(data, args);show()
    else:
        quick_conductance_estimate(data, args)

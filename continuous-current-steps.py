import sys, os
sys.path.append('../../')
from data_analysis.IO import rtxi # download data_analysis module at https://bitbucket.org/yzerlaut/data_analysis
import numpy as np
from graphs.my_graph import *  # download graphs module at https://bitbucket.org/yzerlaut/graphs
from data_analysis.manipulation.files import * # download data_analysis module at https://bitbucket.org/yzerlaut/
sys.path.append('../synaptic-barrages/')
from synaptic_barrages import compute_network_states_and_responses

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

        
def load_data(filename):
    
    # load HDF5 (see data_analysis.IO module)
    data = rtxi.load_continous_RTXI_recording(filename, with_metadata=True)
    data['filename'] = filename
    
    # protocol-dependent processing of the input
    sort_keys(data)

    data['start_vector'] = np.array(data['start_vector'])
    data['stop_vector'] = np.array(data['stop_vector'])
    data['amplitude_vector'] = np.array(data['amplitude_vector'])
    
    # other...
    data['t'] = np.arange(len(data[data['Vm']]))*data['dt']
    # stim duration is constant in this protocol
    for key in data['params']:
        if len(key.split('Increments'))>1:
            Nfreq = data['params'][key]
        if len(key.split('Duration'))>1:
            data['stim_duration'] = 1e-3*data['params'][key]
    
    return data
    
##################################################################
######## ANALYSIS STARTS HERE ####################################
##################################################################

def quick_conductance_estimate(data, args):

    compute_network_states_and_responses(data, args,
                                         keys=['amplitude_vector'])


    ilow_LFP = np.argsort(-data['LFP_levels'])[:args.N_for_quick_conductance_estimate]
    Depol = np.mean(data['Depol_levels_post'][ilow_LFP]-data['Depol_levels_pre'][ilow_LFP])
    conductance_estimate = 1e9*1e-12*data['AMPLITUDE_VECTOR'][0]/Depol # in nS
    
    print('---------------------------')
    print('conductance estimate: ', conductance_estimate)
    f = open(os.path.join(os.getenv("HOME"), 'DATA', 'cell_conductance.txt'), 'w')
    f.write(str(conductance_estimate))
    f.close()
    print('writen in ~/DATA/cell_conductance.txt')
    print('---------------------------')

##################################################################
######## PLOTTING STARTS HERE ####################################
##################################################################

def make_raw_data_figure(data,
                         args,
                         figsize=(.8,.16),
                         Vm_color=Blue,
                         Iinj_color=Orange,
                         LFP_color=Grey,
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

    # # Vm plot with spikes
    Vmin, Vmax = np.min(data[data['Vm']][cond]), np.max(data[data['Vm']][cond])
    ispikes = np.argwhere((data[data['Vm']][cond][1:]>Vpeak) & (data[data['Vm']][cond][:-1]<=Vpeak)).flatten()
    for ii in ispikes:
        ax.plot([data['t'][cond][ii]], [(Vpeak-Vmin)*Vm_enhancement_factor/(Vmax-Vmin)*(Imax-Imin)+2*(Imax-Imin)+Imin], '*',  color=Vm_color, ms=spike_ms)
    data[data['Vm']][data[data['Vm']]>Vpeak] = Vpeak
    ax.plot(data['t'][cond], (data[data['Vm']][cond]-Vmin)*Vm_enhancement_factor/(Vmax-Vmin)*(Imax-Imin)+2*(Imax-Imin)+Imin, color=Vm_color, lw=1)

    # # LFP plot
    LFPmin, LFPmax = np.min(data[data['LFP']][cond]), np.max(data[data['LFP']][cond])
    ax.plot(data['t'][cond], (data[data['LFP']][cond]-LFPmin)/(LFPmax-LFPmin)*(Imax-Imin)+(1.6+Vm_enhancement_factor)*(Imax-Imin)+Imin, color=LFP_color, lw=1)

    condD = np.array(data['start_vector'])<args.tzoom[1]
    for ts, te, fe in zip(data['start_vector'][condD], data['stop_vector'][condD], data['amplitude_vector'][condD]):
        ax.fill_between([ts, te], Imin*np.ones(2), np.ones(2)*2*(Imax-Imin)+Imin-.5*(Imax-Imin), color=Pink, alpha=0.3, lw=0)
        if args.debug:
            ax.annotate(r'$I_{stim}$='+str(int(fe))+'pA',  (te, Imax+Imin))
                
    set_plot(ax, [], xlim=[data['t'][cond][0], data['t'][cond][-1]])
    ax.annotate('$I_{inj}$', (args.tzoom[0], Imin), color=Iinj_color)
    ax.annotate('$V_{m}$', (args.tzoom[0], 2*(Imax-Imin)+Imin), color=Vm_color)
    ax.annotate('LFP', (args.tzoom[0], 5*(Imax-Imin)+Imin), color=LFP_color)
    ax.plot([args.tzoom[0]+.1*np.diff(args.tzoom)[0],args.tzoom[0]+.1*np.diff(args.tzoom)[0]+args.Tbar], [Imin, Imin], 'k-', lw=2)
    ax.plot([args.tzoom[0]+.1*np.diff(args.tzoom)[0],args.tzoom[0]+.1*np.diff(args.tzoom)[0]], [Imin, Imin+args.Ibar*1e-12], 'k-', lw=2)
    if args.Tbar<1:
        ax.annotate(str(int(1e3*args.Tbar))+'ms', (args.tzoom[0]+.1*np.diff(args.tzoom)[0], Imin), color='k')
    else:
        ax.annotate(str(np.round(args.Tbar,1))+'s', (args.tzoom[0]+.1*np.diff(args.tzoom)[0], Imin), color='k')
    ax.annotate(str(int(args.Ibar))+'pA', (args.tzoom[0]+.1*np.diff(args.tzoom)[0], Imin+args.Ibar*1e-12), color=Iinj_color, rotation=90)
    ax.annotate(str(np.round(1e3*args.Ibar*1e-12/Vm_enhancement_factor*(Vmax-Vmin)/(Imax-Imin),1))+'mV',\
                (args.tzoom[0], Imin+Imax+args.Ibar*1e-12), color=Vm_color, rotation=90)
    ax.annotate(str(np.round(1e6*args.Ibar*1e-12*(LFPmax-LFPmin)/(Imax-Imin),1))+'uV',\
                (args.tzoom[0], Imin+3*Imax+args.Ibar*1e-12), color=LFP_color, rotation=90)
    return fig, ax


def make_trial_average_figure(data, args):
    """

    """
    compute_network_states_and_responses(data, args, keys=['amplitude_vector'])

    data['amplitude_levels'] = np.unique(data['AMPLITUDE_VECTOR'])
    
    fig, AX = figure(figsize=(.2*len(data['amplitude_levels']),.4),
                     axes=(2, len(data['amplitude_levels'])),
                     wspace=0.2,
                     left=0.25, top=.8)
    
    number_of_common_trials = 1000
    for a, f in enumerate(data['amplitude_levels']):
        # loop over frequency levels
        cond = (data['AMPLITUDE_VECTOR']==f)
        for i in range(args.N_state_discretization):
            true_cond = data['cond_state_'+str(i+1)] & cond
            AX[1][a].plot(1e3*data['t_window'],
               1e3*data['Vm_Responses'][true_cond,:].mean(axis=0),
                          '-', color=COLORS[i], lw=2)
            AX[1][a].fill_between(1e3*data['t_window'],
                                  1e3*data['Vm_Responses'][true_cond,:].mean(axis=0)+1e3*data['Vm_Responses'][true_cond,:].std(axis=0),
                                  1e3*data['Vm_Responses'][true_cond,:].mean(axis=0)-1e3*data['Vm_Responses'][true_cond,:].std(axis=0),
                                  lw=0., color=COLORS[i], alpha=.3)
            # for the raster plot, we want a vcommon trial number
            number_of_common_trials = np.min([number_of_common_trials,\
                                              len(data['AMPLITUDE_VECTOR'][true_cond])])

    for a, f in enumerate(data['amplitude_levels']):
        # loop over frequency levels
        cond = (data['AMPLITUDE_VECTOR']==f)
        for i in range(args.N_state_discretization):
            true_cond = data['cond_state_'+str(i+1)] & cond
            for k, s in enumerate(np.arange(len(true_cond))[true_cond][:number_of_common_trials]):
                spk_train = data['Spike_Responses'][s]
                AX[0][a].plot(1e3*spk_train, 0*spk_train+k+i*(number_of_common_trials+2), 'o',
                              color=COLORS[i], ms=args.ms)

            AX[0][a].fill_between([1e3*data['t_window'][0],1e3*data['t_window'][-1]],
                                  i*(number_of_common_trials+2)*np.ones(2)-1, 
                                  i*(number_of_common_trials+2)*np.ones(2)+number_of_common_trials, 
                                  color=COLORS[i], alpha=.3, lw=0)
                
        AX[0][a].set_title(r'$I_{stim}$='+str(int(f))+'pA')
        AX[1][a].plot([0,0], args.Vm_lim, 'w.', ms=1e-8, alpha=0)
        if (a==0):
            AX[0][a].plot(1e3*data['t_window'][0]*np.ones(2),
                      args.N_state_discretization*(number_of_common_trials+2)-np.arange(2)*number_of_common_trials-2,
                      'k-', lw=1)
            AX[0][a].annotate(str(number_of_common_trials)+'trials', (1e3*data['t_window'][0],
                                                                      args.N_state_discretization*(number_of_common_trials+2)))
            set_plot(AX[1][a], xlabel='time from stim. (ms)', ylabel='Vm (mV)', ylim=args.Vm_lim)
            set_plot(AX[0][a], ['bottom'], ylabel='Spikes', ylim =[-3, AX[0][a].get_ylim()[1]+3])
        else:
            set_plot(AX[0][a], ['bottom'])
            set_plot(AX[1][a], xlabel='time from stim. (ms)', yticks_labels=[], ylim=args.Vm_lim)
            
    return fig, AX

def make_conductance_fig(data, args):
    """

    """
    
    compute_network_states_and_responses(data, args, keys=['amplitude_vector'])

    fig, ax = figure()
    for i in range(args.N_state_discretization):
        true_cond = data['cond_state_'+str(i+1)]
        Depol = np.mean(data['Depol_levels_post'][true_cond]-data['Depol_levels_pre'][true_cond])
        muG = 1e9*np.mean(1e-12*data['AMPLITUDE_VECTOR'][0]/Depol)
        ax.bar([i], [muG], color=COLORS[i])

    set_plot(ax, ylabel='input \n resistance (nS)',
             xticks=range(args.N_state_discretization),
             xticks_labels=['BA', 'IA', 'SA'])
    return fig, ax

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
    parser.add_argument("--Vspike", help="", type=float, default=-30)
    parser.add_argument("--N_state_discretization", help="", type=int, default=3)
    
    
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
    COLORS = [get_linear_colormap(Orange, Blue)(i/(args.N_state_discretization-1)) for i in range(args.N_state_discretization)]

    print('[...] loading data')
    data = load_data(args.filename)
    print('[...] analyzing')
    if args.trial_average_analysis:
        fig, _ = make_trial_average_figure(data, args)
    elif args.conductance_analysis:
        fig, _ = make_conductance_fig(data, args)
    elif args.depol_analysis:
        fig, _ = make_depol_fig(data, args)
    elif args.raw_data:
        fig, _ = make_raw_data_figure(data, args)
    else:
        quick_conductance_estimate(data, args)

    if args.save:
        fig.savefig(desktop+'fig.svg')
    else:
        show()
        

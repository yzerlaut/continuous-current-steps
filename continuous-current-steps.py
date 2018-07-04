import sys
sys.path.append('../../')
from data_analysis.IO import hdf5 # dowload data_analysis module at https://bitbucket.org/yzerlaut/data_analysis
import numpy as np

if __name__ == '__main__':
    if len(sys.argv)==1:
        print('---------------------------------------------')
        print('you should provide a hdf5 file as an argument')
        print('---------------------------------------------')
    else:
        filename = sys.argv[-1]
        data = hdf5.load_continous_RTXI_recording(filename)
        print(data)
    



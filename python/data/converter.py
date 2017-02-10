"""
Converts a .mat file to the appropriate JSON format
"""
import scipy.io as sio
import os
import json

def convert_data_set(lfile):
    """
    Converts a .mat file to the appropriate JSON format
    """
    directory = os.path.dirname(os.path.abspath(__file__))
    try:
        data = sio.loadmat(lfile)
    except:
        raise Exception(("Could not load file, check check that %s exists in"+
                         "the directory %s.") % (lfile,directory))


    ouput = {
        'info':'%s - %s' %(data['__header__'],data['data'][0][0][10][0]),
        'dimensions':{
            'numberOfAnchors': data['data'][0][0][0].shape[0],
            'numberOfMeasurements':data['data'][0][0][0].shape[1]
        },
        'ranges':data['data'][0][0][0].tolist(),
        'acceleration': None,
        'rates': None,
        'thrust': None,
        'torques': None,
        'times': None
    }
    
    try:
        # Saves data in the JSON format in the filename in which it was loaded
        path = os.path.join(directory, lfile[:-4] + '.json')
        with open(path, 'w') as outfile:
            json.dump(ouput, outfile , separators=(',', ':'), sort_keys=True, indent=4)
    except:
        raise Exception('Could not save the file')
    

if __name__ == "__main__":
    loadFileName = 'data4.mat'
    convertDataSet(loadFileName)
    
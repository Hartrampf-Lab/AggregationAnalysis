###########################################################################
#
#
#   AFPS Data Parsing Program - afps_dataiter.py
#   Pietro Luigi Willi, Bálint Tamas, Nina Hartrampf
#   TITLE OF PAPER
#   UZH, Institute of Chemistry.
#   23/06/2023, Zurich.
#
#
###########################################################################

import csv
import os
from collections import namedtuple
import matplotlib.pyplot as plt
import matplotlib as mpl
from itertools import chain
import time
import numpy as np
import peakutils
import pandas as pd
import numpy as np
from datetime import datetime
from scipy.optimize import curve_fit
import math 
import warnings
warnings.simplefilter("ignore")


mpl.rcParams["xtick.top"] = False
mpl.rcParams["ytick.right"] = False
mpl.rcParams['lines.linewidth'] = 1


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print ('%r  %2.2f ms' % \
                  (method.__name__, (te-ts) * 1000))
        return result

    return timed

class Synthesis:
    """
    This data analysis tool was developed to analyse the data produced by the 
    Automated Fast-flow Peptide Synthesiser (AFPS). This script is an edited 
    version of afps.py (https://github.com/amijalis/afps-integration/tree/0.1)
    written by Alexander J. Mijalis, published as part of the paper "A fully 
    automated flow-based approach for accelerated peptide synthesis"
    (doi: 10.1038/nchembio.2318).

    The class has the ability of opening the .pep files produced by the AFPS 
    by instantiating a Synthesis object from the serial number of the synthesis
    using the class method from_serial. This allows it to gain access to the 
    raw data of the AFPS, which is contianed in a file which is named used the 
    following automated naming convention:
    YYYYMMDD_HHMMSS_RAW_AFPS80_SN{6 digit serial number}

    for example:
    20230524_160826_RAW_AFPS80_SN001030

    This file is composed of various columns containing the data from all the
    temperature and pressure sensors, the position of the valves, pump commands
    and aruments, current step number and most importantly, the absorption of
    the UV-Vis module. The data is sampled at 50 Hz. 
    
    The synthesiser is automating SPPS, a series of recursive coupling and
    deprotection reactions used to synthesise peptides. By plotting the raw UV 
    trace, two distinct type of peaks should become aparent: a large and vastly
    oversaturated peak belonging to the coupling step followed by a smaller peak
    belonging to deprotection. The deprotection peka is the peak of interest for 
    determining the success of the coupling and the influence and onset of 
    aggregation.

    The deprotection steps can be isolated by querying the corresponding step 
    number. The first deprotection happens at step 36 and then every 22 steps.
    The detect_depro_steps method takes care of this querying and returns the 
    steps numbers corresponding to deprotection steps. Some parameters are then 
    computed from the isolated peak. These are the maximum, the full width at 
    half maximum, the area under the peak and the angle. 
    
    The angle is a novity introduced by our paper and approach. It differs from
    the other parameters due to the complexity required to compute it. 
    Its complexity is derived from the fitting of a gaussian function onto the 
    peak and from the need to increase the robustness of the method. The whole 
    angle computation is split into 2 main methods, process_peak and get_angle. 
    The process_peak method is given the raw deprotection peak as input. The peak
    is first split into front and tail by calling the static method split_peak. 
    A threshold is set so that any absorbance values larger than the threshold are
    set to the threshold + 0.001. This ensures that all the oversaturated parts of
    the peak have the same values. The middle of the peak is found by finding
    the median time of the maximum absorbance. The thresholding makes it easier to
    find the middle of the peak in case of oversaturation. The peak is then split
    into front and tail along the maxima. Using the static method mirror_peak, 
    the two half peaks are then individually mirrored and concatenated with the 
    original halfpeak, resulting in two full gaussian peaks. Finally the static 
    method trim_peak, removes the oversaturated part of the peak by simply 
    deleting the points with an absorbance higher than the threshold by their 
    index. In addition it also trims the peaks at their minima, making the peak 
    more gaussian-like. These three static methods are all called by the 
    process_peak method in the order described to yield the x and y values of the
    front and tail peaks. 
    The get_angle method then calls the process_peak method and fits a parametrized
    gaussian function with the following formula:
    f(x)=a*e^((-(x-b)^2)/2c^2 )
    The curve fitting yields the parameters a, b and c for the front and the tail 
    peak. The parameters are then used to find the maximum gradient by of the 
    function f(x). The maximum gradient corresponds to the gradient at the point of
    inflection. The gradient of the function f(x) is defined by its first order 
    derivative, f'(x):
    f'(x)=((b-x)/c^2)*a*e^((-(x-b)^2)/2c^2 ) 
    In an isosceles triangle the tangent of an angle of a corner is defined as the
    opposite edge over the adjacent edge. Therefore the arctangent of one over the 
    gradient at the inflection will yield the angle of interest. After performing 
    the angle computation on the front and tail peak, the two resultant angles are 
    added together to yield the peak angle. 
    The get_integrals method iterates through the steps that are defined as 
    deprotections and computes the various summary statistics. The angle is 
    normalized using the resin mass. The mass normalization is performed to remove 
    the peak angle's dependence on resin mass. The standardization normalizes the 
    peak angle to the standard resin mass, 150 mg, using the following formula:
    θ_st=180-150/m(180-θ)

    --------------------------------------
    Attributes:
    pepfile : str
        name and path of the .pep file. this is automatically assigned by calling 
        the Syntheisis class with the from_serial class method with the serial 
        number of interst as the argument (See Usage). 

    --------------------------------------
    Usage:
    - Define a Synthesis instance by calling the from_serial classmethod on the 
      Synthesis class, with the serial number of the synthesis you want to analyse 
      as the argument.
    >>> syn = Synthesis.from_serial(serial_n)
    
    - Plot the raw uv trace by calling plot_uv on the Synthesis instance you 
      have created:
    >>> syn.plot_uv()

    - Detect the deprotection steps and iterate through them, calculating all the 
      peak parameters, like angle area height and width, by calling get_integrals
      on the Synthesis instance. This will return a list of tuples where every tuple
      corresponds to the data extracted from the deprotection of one amino acid.
      The length of the list will be equivalent to the number of amino acids in the 
      sequence:
    >>> integrals = syn.get_integrals() 

    
    """
    def __init__(self,pepfile):
        self.stdblue = "#265952"
        self.stdred = "#DD7965"
        self.stdyellow = "#E1C009"
        self.pepfile = pepfile
        self.serial_no, self.sequence, self.datafile  = self.get_afps_sequence(pepfile)
        self.afps_data = self.parse_afps_datafile()
        

    @timeit
    def parse_afps_datafile(self):
        "Parse the afps datafile, returning a list of dictionaries containing each timepont's data"
        with open(self.datafile) as csvfile:
            csvreader = csv.reader(csvfile)
            next(csvreader,None)
            sequence = ''.join(next(csvreader)).split(": ")[1]

            for i in range(3):
                if i == 0:
                    self.resinmass = next(csvreader,None)
                    print(self.resinmass)
                if i == 1:
                    self.notenum = next(csvreader,None)
                    print(self.notenum)
                if i == 2: 
                    self.stepnum = next(csvreader,None)
                    print(self.stepnum)

            headers = next(csvreader)
            # Put the data in a dictionary. date and time values should be treated as strings,
            # step_no as an int, and the rest as floats. 
            afps_data = [{key:value if key in ['DATE_YYYYMMD', 'TIME_HHMMSSS.SSS'] \
                        else int(value.split('.')[0]) if key == 'STEP_NO'\
                        else float(value) for (key,value) in zip(headers,row)} \
                        for row in csvreader]

            for row in afps_data:
                hrs = int(row['TIME_HHMMSSS.SSS'][0:2])
                mins = int(row['TIME_HHMMSSS.SSS'][2:4])
                secs = int(row['TIME_HHMMSSS.SSS'][4:6])
                ms = int(row['TIME_HHMMSSS.SSS'][7:10])
                float_ms = float(hrs*60*60*1000 + mins*60*1000 + secs*1000 + ms)
                row['TIME_MS'] = float_ms
            return afps_data
    
    @staticmethod
    def get_afps_sequence(pepfile):
        "Returns (sequence, serial_no, datafile) for a .pep file input."

        if not pepfile.endswith('/'): pepfile += '/'
        AFPS_File = namedtuple('AFPS_File', 'serial_no sequence afps_datafile_path')
        for afps_file in os.listdir(pepfile):
            if afps_file.split('.')[-1] == 'afps':
                serial_no = int(afps_file.split('.')[0].split('_')[-1].strip('SN'))
                afps_datafile_path = pepfile + afps_file
                
                with open(afps_datafile_path) as csvfile:
                    csvreader = csv.reader(csvfile)
                    next(csvreader,None)
                    sequence = ''.join(next(csvreader)).split(": ")[1]
                
                return AFPS_File(serial_no = serial_no, \
                                sequence = sequence,   \
                                afps_datafile_path = afps_datafile_path)

    @staticmethod
    def get_pep_files(path):
        "Yields .pep files in path." 
        if not path.endswith('/'): path = path + '/'
        for afps_dir in os.listdir(path):
            if afps_dir.split('.')[-1] == 'pep':
                yield(path + afps_dir)

    @classmethod
    def from_serial(cls, serial_no, path='data/'):
        "Create a Synthesis object from serial number."
        if not path.endswith('/'): path += '/'
        for item in os.listdir(path):
            extension = item.split('.')[-1]
            name = item.split('.')[0]
            if extension == 'pep':
                if int(name.split('_')[-1].strip('SN')) == serial_no:
                    return cls(path+item)
        print('Serial Number', serial_no, 'not found')
        return None

    def plot_uv(self):
        "Plots the raw UV data using matplotlib"
        try:            
            time, uv, step = zip(*((row['TIME_MS'], row['UV-VIS'], row['STEP_NO']) for row in self.afps_data))
        except KeyError:
            time, uv, step = zip(*((row['TIME_MS'], row['UV_VIS'], row['STEP_NO']) for row in self.afps_data))
        fig, line = plt.subplots()
        line.plot(time, uv, 'black')
        plt.xlabel('Time')
        plt.ylabel('Intesity')
        plt.show()
        
    def plot_press(self):
        "Plots the raw pressure data using matplotlib"
        try:            
            time, pres1, pres2, step = zip(*((row['TIME_MS'], row['PRESSURE_1'], row['PRESSURE_2'], row['STEP_NO']) for row in self.afps_data))#select and assign the columns of interest
        except KeyError:
            time, pres1, pres2, step = zip(*((row['TIME_MS'], row['PRESSURE_1'], row['PRESSURE_2'], row['STEP_NO']) for row in self.afps_data))#selects the pressure columns from the raw afps data and returns it as a tuple
        fig, line = plt.subplots()
        line.plot(time, pres1, label="Pressure 1")
        line.plot(time, pres2, label="Pressure 2")
        plt.xlabel('Time')
        plt.ylabel('Pressure')
        plt.legend()
        plt.show()
        
    def plot_temp(self):
        "Plots the raw temperature data using matplotlib"
        try:            
            time, T_1, T_2, T_3, T_4, T_5, T_6, step = zip(*((row['TIME_MS'], row["TC_1"], row["TC_2"], row["TC_3"], row["TC_4"], row["TC_5"], row["TC_6"], row['STEP_NO']) for row in self.afps_data))#select and assign the columns of interest from the raw afps
        except KeyError:
            time, T_1, T_2, T_3, T_4, T_5, T_6, step = zip(*((row['TIME_MS'], row["TC_1"], row["TC_2"], row["TC_3"], row["TC_4"], row["TC_5"], row["TC_6"], row['STEP_NO']) for row in self.afps_data))#select and assign the columns of interest
        #here all the temperatures were extracted
        fig, line = plt.subplots()
        line.plot(time, T_1, label="Temperature 1")
        line.plot(time, T_2, label="Temperature 2")
        line.plot(time, T_3, label="Temperature 3")
        line.plot(time, T_4, label="Temperature 4")
        line.plot(time, T_5, label="Temperature 5")
        line.plot(time, T_6, label="Temperature 6")
        plt.xlabel('Time')
        plt.ylabel('Temperature')
        plt.ylim(0,150)
        plt.legend()
        plt.show()
    
    def plot_depro_temp(self): #lets check how stable the T is between deprotection steps
        try:            
            time, T_1, T_2, T_3, T_4, T_5, T_6, step = zip(*((row['TIME_MS'], row["TC_1"], row["TC_2"], row["TC_3"], row["TC_4"], row["TC_5"], row["TC_6"], row['STEP_NO']) for row in self.all_depro_data))#select and assign the columns of interest
        except KeyError:
            time, T_1, T_2, T_3, T_4, T_5, T_6, step = zip(*((row['TIME_MS'], row["TC_1"], row["TC_2"], row["TC_3"], row["TC_4"], row["TC_5"], row["TC_6"], row['STEP_NO']) for row in self.all_depro_data))#select and assign the columns of interest
        #what we finally want is the temperature inside the reactor.
        fig, line = plt.subplots()
        line.plot(time, T_1, label="Temperature 1")
        line.plot(time, T_2, label="Temperature 2")
        line.plot(time, T_3, label="Temperature 3")
        line.plot(time, T_4, label="Temperature 4")
        line.plot(time, T_5, label="Temperature 5")
        line.plot(time, T_6, label="Temperature 6")
        plt.xlabel('Time')
        plt.ylabel('Temperature')
        plt.ylim(15,100)
        plt.legend()
        plt.show()
    
    def plot_depro_press(self):
        "Plots the raw pressure data using matplotlib"
        try:            
            time, pres1, pres2, step = zip(*((row['TIME_MS'], row['PRESSURE_1'],row['PRESSURE_2'], row['STEP_NO']) for row in self.all_depro_data))#select and assign the columns of interest
        except KeyError:
            time, pres1, pres2, step = zip(*((row['TIME_MS'], row['PRESSURE_1'],row['PRESSURE_2'], row['STEP_NO']) for row in self.all_depro_data))
        fig, line = plt.subplots()
        line.plot(time, pres1, label="Pressure 1")
        line.plot(time, pres2, label="Pressure 2")
        plt.xlabel('Time')
        plt.ylabel('Pressure')
        plt.legend()
        plt.show()
        
    def plot_depro_temp(self):
        '''
        Look at VALVE_8_POS in order to determine which preheating loop was used.
        Then check the temperature of the preheating loop that was effectively used.
        
        HX ID,POS (8),TEMP [DEGC],LENGTH [FT]
        1    ,1      ,30         ,1
        2    ,2      ,60         ,2.5
        3    ,3      ,90         ,5
        4    ,4      ,90         ,10
        
        TC_1: reactor nozzle
        TC_2: reactor
        TC_3: 1ft loop (at 30° C)
        TC_4: 2.5ft loop (60° C)
        TC_5: 5ft loop (90° C)
        TC_6: 10ft loop (90° C)
        '''
        pos_TC = {0:"TC_6", 1:"TC_3", 2:"TC_4", 3:"TC_5", 4:"TC_6"} #maps the position of valve 8 to the name of the column contianing the temperature of the respective heating loop. 
        try:            
            time, temp, step = zip(*((row['TIME_MS'], row[pos_TC[row["VALVE_8_POS"]]], row['STEP_NO']) for row in self.all_depro_data))#row["VALVE_8_POS"] selects the position of valve 8, this one is then used as the index on the pos_TC dictionary which outputs a string which corresponds to the column contining the temperature of the loop which was effectively used 
        except KeyError:
            time, temp, step = zip(*((row['TIME_MS'],row[pos_TC[row["VALVE_8_POS"]]], row['STEP_NO']) for row in self.all_depro_data))
        print(len(temp))
        fig, line = plt.subplots()
        line.plot(time, temp, label="Activation Loop Temperature")
        plt.xlabel('Time')
        plt.ylabel('Temperature')
        plt.legend()
        plt.show()
        
        
    def integral_temp(self):
        pos_TC = {0:"TC_6", 1:"TC_3", 2:"TC_4", 3:"TC_5", 4:"TC_6"} #maps the position of valve 8 to the name of the column contianing the temperature of the respective heating loop. 
        #what we finally want is the temperature of the preheating loop which was effectively used.
        #for this we need to look at the temperature of the preheating loop selected by valve 8
        try:            
            time, temp, step = zip(*((row['TIME_MS'], row[pos_TC[row["VALVE_8_POS"]]], row['STEP_NO']) for row in self.all_depro_data)) #row["VALVE_8_POS"] selects the position of valve 8, this one is then used as the index on the pos_TC dictionary which outputs a string which corresponds to the column contining the temperature of the loop which was effectively used 
        except KeyError: #if it fails, just do it again XD.
            time, temp, step = zip(*((row['TIME_MS'], row[pos_TC[row["VALVE_8_POS"]]], row['STEP_NO']) for row in self.all_depro_data))
        temp_df = np.array([time, temp, step]).T
        return temp_df
    
    def plot_uv_title(self, title):
        "Plots the raw UV data using matplotlib"
        try:            
            time, uv, step = zip(*((row['TIME_MS'], row['UV-VIS'], row['STEP_NO']) for row in self.afps_data))
        except KeyError:
            time, uv, step = zip(*((row['TIME_MS'], row['UV_VIS'], row['STEP_NO']) for row in self.afps_data))
        fig, line = plt.subplots()
        fig.suptitle(title)
        line.plot(time, uv)
        plt.xlabel('Time')
        plt.ylabel('Intensity')
        plt.show()

    def plot_raw(self, field):
        time, field, step =  zip(*((row['TIME_MS'], row[field], row['STEP_NO']) for row in self.afps_data))
        fig, line = plt.subplots()
        line.plot(time, field)
        plt.show()

    def plot_multiple(self, field1, field2, title):
        time, field1, field2, = zip(*((row['TIME_MS'], row[field1], row[field2]) for row in self.afps_data))
        fig, ax1 = plt.subplots()

        color = 'k'
        ax1.set_xlabel('time (s)')
        ax1.set_ylabel('field1', color=color)
        ax1.plot(time, field1, color=color)
        ax1.tick_params(axis='y', labelcolor=color)

        ax2 = ax1.twinx() #instantiate a second axes that shares the same x-axis
        color = 'tab:red'
        ax2.set_ylabel('field2', color=color) #we already handled the x-label with ax1
        ax2.plot(time, field2, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        fig.tight_layout() #otherwise the right y-label is slightly clipped
        plt.title(title)
        plt.legend([field1, field2])
        plt.show()

    def get_data(self, *fields):
        # Returns the specified columns
        return tuple(zip(*([row[field] for field in fields] for row in self.afps_data)))

    def correct_uv_baseline(self, order=3):
        "Returns a baseline corrected UV trace. Requires peakutils and numpy."
        try:
            time, step_no, numpy_uv = np.array(self.get_data('TIME_MS', 'STEP_NO', 'UV_VIS'))
        except KeyError:
            time, step_no, numpy_uv = np.array(self.get_data('TIME_MS', 'STEP_NO', 'UV-VIS'))
        baseline = peakutils.baseline(numpy_uv, order)
        numpy_uv_corrected = numpy_uv - baseline
        return list(zip(step_no, time, numpy_uv_corrected))
    
    def detect_depro_steps(self):
        # Parse out UV peaks by step number. 42, 64, 88....
        first_depro_step = 35
        depro_step_lag = 22
        self.all_depro_data = [row for row in self.afps_data if (row['STEP_NO']-first_depro_step) % depro_step_lag == 0 and row['STEP_NO'] >= first_depro_step] #define 22, 35 as a variable

        # Create a set of all unique step numbers captured here (one for each deprotection peak)
        depro_steps = {row['STEP_NO'] for row in self.all_depro_data}
        self.depro_steps = depro_steps

        # This puts all of the data belonging to each depro peak
        # in a dictionary indexed by the step number it belongs to.
        
        print(len(depro_steps), 'deprotection peaks found')# How many steps did we get?
        print(len(self.sequence), 'residues in sequence')# How many residues in our sequence?
        return depro_steps
    
    def clean_up(self):
        """this function is used to collect the sequence length, number of deprotections, synthesis file name, n. of steps and the resin mass. These metadata are used by the full_data_converter program to make decision whether to include or exclude a particular sequence into the converted dataset."""
        return len(self.detect_depro_steps()), len(self.sequence), self.notenum, self.stepnum, self.resinmass 
    
    def integrate(self, peak):
        # Performs the riemann sum underneath the given curve with format [*[time,data]].
        # Sort the data by time
        peak.sort(key=lambda x: x[0])
        peak = np.array(peak)
        #Vectorization speeds up computation
        dt = peak[1:,0] - peak[:-1,0]
        avg_uv = (peak[1:,1] + peak[:-1,1]) / 2 #midpoint rule
        area = dt * avg_uv
        integral = np.sum(area)
        return integral

    def get_max(self, peak):
        time,uv = zip(*peak)
        maximum = max(uv)
        half_max = maximum / 2
        right = 0
        left = 0
        for time, uv in peak:
            if uv - half_max > 0:
                left = time

        for time, uv in reversed(peak):
            if uv - half_max > 0:
                right = time
        # time is in ms; return in units of seconds
        width = abs(right-left) / 1000
        if width == 0:
            print('Warning: Peak width 0')
        return maximum, width
    
    def get_arclength(self, peak):
        """
        This method calculates the distance between each subsequent pair of points
        using the pythagorean equation: c**2 = a**2 + b**2
        """
        # Sort the data by time
        peak.sort(key=lambda x: x[0])
        peak = np.array(peak)
        dt = peak[1:,0] - peak[:-1,0]
        dy = (peak[1:,1] - peak[:-1,1])*10**6
        length =  np.sqrt(dt**2+dy**2)
        arclength = np.sum(length)
        return arclength
    
    @staticmethod
    def split_peak(df, threshold):
        """
        The split peak method is used by the process_peak method to split the peak into
        a front and tail part.
        It cuts the peak in half at the index of the median of the maxima (index = time).
        Why not just cut the peak in half at the maxima?
        There are often several maxima with the same value,especially when the signal
        is oversaturated.
        it is thus necessary to take the median of the indices of all the maxima.
        To make the method robust to outliers, UV-absorption values above the 
        threshold value are set to the threshold value + 0.001. This makes it 
        easier to find the middle of the peak in case of oversaturation.
        This method also removes all the data points that corresponds to oversaturated
        maxima in order to make the fitting in the get_angle method more accurate.
        """
        tzero = datetime.now() #only important for timing the whole process
        df_np = np.array(df)
        #df_np = self.normalizer(df_)
        df_np[np.argwhere(df_np[:, 1] >= threshold).flatten(),1] = threshold + 0.001 #anything larger or equal to 0.9 is regarder as overasaturated and set equal to 0.901 to make it easier to remove these points later, but also makes it easier to find the real median
        maxUV_i = (np.argwhere(df_np[:, 1] == np.amax(df_np[:,1]))).flatten() #find all the indices corresponding to the uv maxima.
        maxTime = df_np[:,0][maxUV_i] #index the time with the indices of the maxima in order to find the time of the maxima. 
        maxTime_median = np.median(maxTime) #take the median of the time of the maxima
        maxTime_median_i = int(np.median(maxUV_i)) #get the mean of the indices of the maxima

        peak_front = df_np[:maxTime_median_i+1] #splits the peak by the median. +1 so that the median is included in both the tail and the front
        peak_tail = df_np[maxTime_median_i:] #splits the peak by the median.
        front_corr = np.expand_dims((peak_front[: ,0]-maxTime_median), axis=1) 
        tail_corr = np.expand_dims((peak_tail[:, 0]-maxTime_median), axis=1)#corrects the time so that the median time = 0. expand dims adds a new axis to the corrected time array so that an array of shape (n,) becomes of shape (n,1) or (1,n). This necessary for the concatenation that follows.
        peak_front = np.concatenate((peak_front, front_corr),axis=1)
        peak_tail = np.concatenate((peak_tail, tail_corr),axis=1)
        return peak_front, peak_tail
    
    @staticmethod
    def mirror_peak(halfpeak):
        """
        The mirror_peak methdod takes the half-peak split by the split_peak method
        and mirrors it into a full peak. The mirroring is done by flipping the 
        half-peak and concatenating it with the original in such a way that the 
        maxima are merged side to side.  
        """
        peak = halfpeak.copy()
        peak[:, 2] = -peak[:, 2] #necessary for mirroring of peak. we are taking the negative of the time (corrected and uncorrected)
        peak[:, 0] = -peak[:, 0]
        peak = np.flipud(peak) #flipping of the negated time component
        if halfpeak[0,1] > halfpeak[-1,1]:
            mirrored_peak = np.concatenate([peak, halfpeak],axis=0)
        else:
            mirrored_peak = np.concatenate([halfpeak, peak],axis=0)
        x_peak = mirrored_peak[:, 2]
        y_peak = mirrored_peak[:, 1]
        return x_peak, y_peak
    
    @staticmethod
    def trim_peak(x_peak, y_peak, threshold):
        """
        The trim_peak method, is used by the process_peak method to remove the part of the peak
        that is oversaturated. In the split_peak method, the parts of the peak that are larger 
        than the threshold are set to 0.001 + threshold. In this way the index of the oversaturated 
        parts of the peak are very easy to find using argwhere. Using the index, these parts are then 
        removed from the x and y component. This is necessary to improve the quality of the fit.
        """
        peak_saturation = np.argwhere(y_peak > threshold) #remember that all the values larger than 0.9 were set to 0.901. now we find the index of all these values...
        x_peak = np.delete(x_peak, peak_saturation)
        y_peak = np.delete(y_peak, peak_saturation) #...and delete the datapoints using the index
        y_tail_min = min(y_peak)
        y_peak -= y_tail_min #make baseline (min of tail is usually lower than front)
        x_peak = x_peak[int(np.argwhere(y_peak == min(y_peak))[0]):int(np.argwhere(y_peak == min(y_peak))[1])+1] #select all values between the two minima. makes a nicer gaussian shape.
        y_peak = y_peak[int(np.argwhere(y_peak == min(y_peak))[0]):int(np.argwhere(y_peak == min(y_peak))[1])+1]
        return x_peak, y_peak

    def process_peak(self, df, threshold):
        """
        The process_peak method calls the split_peak, mirror_peak and trim_peak methods. 
        The second part of the code in this method takes care of the fact that the bottom of the 
        tailing part of the peak deviates from the gaussian character required to achieve an accurate
        fit at the middel of the peak, which is where the point of inflection occurs.
        """
        peak_front, peak_tail = self.split_peak(df=df, threshold=threshold)
        x_front, y_front = self.mirror_peak(peak_front)
        x_tail, y_tail = (self.mirror_peak(peak_tail))
        x_front, y_front = self.trim_peak(x_front, y_front, threshold)
        x_tail, y_tail = self.trim_peak(x_tail, y_tail, threshold)
        bottom_trimm = 0.03
        yzero = np.argmin(np.abs(y_tail-(max(y_tail)*bottom_trimm))) #heuristics to make the shape more gaussian-like. the slope of the peak front is often too large at the front (where the gradient should be 0). this part is hereby removed. the y-value which was equal to 0.03 of the maximum becomes 0 and thanks to the absolute it becomes also the minimum.
        y_tail = y_tail[yzero:-yzero] #negative indexing because peak is symmetric
        x_tail = x_tail[yzero:-yzero]
        return x_front, y_front, x_tail, y_tail
    
    @staticmethod
    def gauss_function(x, a, b, c):
        """
        This is the gaussian function fitted by the scipy curve_fit module in the
        get_angle method.
        """
        x = np.array(x)
        return a*np.exp(-(x-b)**2/(2*c**2))
    
    @staticmethod
    def gauss_1st_deriv(x, a, b, c):
        """
        This is the first derivative of the fitted gaussian function gauss_function.
        It is used in the computation of the angle. The angle is the arctangent of 
        1 over the gradient at the inflection point. The gradient at the inflection
        point can be found by simply taking the max of the first derivative.
        """
        return (a*np.exp(-(x-b)**2/(2*c**2))*((b-x)/c**2))
    
    @staticmethod
    def compute_angle(grad):
        angle = (180/math.pi) * np.arctan(1/grad)
        return angle

    def get_angle(self, df, params_front, params_tail, threshold=0.9):
        """
        This is is the method that is called to compute the angles of the deprotection
        peaks. It computes the angles of the peak by first calling the split peak 
        method, which splits the deprotection peak into two parts. A gaussian function
        is then fitted onto each part separately. Then the angle is calculated by 
        taking the arctangent of 1 over the gradient of the inflection point. The 
        angles of the two parts are then added together. 
        """
        x_front, y_front, x_tail, y_tail = self.process_peak(df,threshold)
        y_tail *= 1000 
        y_front *= 1000 # the curve optimization algorithm performs best when the x and y values are approximately of the same magnitude.

        params_front = curve_fit(f=self.gauss_function, xdata=x_front, ydata=y_front, p0=params_front)[0] #p0 is the initial guess, params_front are the parameters of the function that will be optimized by the algorithm 
        
        params_tail = curve_fit(f=self.gauss_function, xdata=x_tail, ydata=y_tail, p0=params_tail)[0]
        degree_front = self.compute_angle(np.max(self.gauss_1st_deriv(np.linspace(min(x_front), max(x_front),100000), *params_front)))
        degree_tail = self.compute_angle(np.max(self.gauss_1st_deriv(np.linspace(min(x_tail), max(x_tail),100000), *params_tail))) #angle calculation requires finding the gradient of the point of inversion. the point of inversion will have the highest gradient and therefore be the maximum of the 1st derivative. tan-1(1/m) where m is the gradient of the point of inversion equals to the angle.
        degree_peak = degree_front + degree_tail
        asymetry_factor = degree_front / degree_tail
        return degree_peak, degree_front, degree_tail, asymetry_factor, params_front, params_tail
    
    def get_integrals(self):
        """
        This method iterates through the steps only selecting the steps which contain 
        the deprotection peak. This is ensured by the detect_depro_steps() method. An
        error handling loop (while loop) ensures that any anomalous synthesis is still
        usable in terms of data by considering previous steps. Various parameters, such 
        as the area (integral), maximum, FWHM (width) and angle, are then calculated.
        """
        steps = self.detect_depro_steps()
        corrected_uv_data = self.correct_uv_baseline() 
        integrals = []
        #peak_list = []
        peak = []
        temp = self.integral_temp()
        p_front = np.array([1, 0, 1000])
        p_tail = np.array([1, 0, 1000])
        resinmass = float(self.resinmass[0].split("RESIN AMOUNT [mg]:")[1]) #mg
        standard_mass = 150 #mg
        extension = 0.3

        for index,step in enumerate(sorted(steps)):
            t_i = np.argwhere(temp[:, -1] == step) #find the index of a particular step number (temp[:,2] corresponds to the step number)
            t_loop = np.array(temp[t_i, 1]) #use this index to find the temperature of the loop for that particular step (temp[:,1] is the temperature).
            t_reac = np.array(temp[t_i, 2]) #use this index to find the temperature of the reactor for that particular step (temp[:,1] is the temperature).
            t_loop_mean = np.mean(t_loop) #take the mean of all the temperature readings within the step.
            t_reac_mean = np.mean(t_reac)
            
            try:
                peak = [ [i[1], i[2]] for i in corrected_uv_data if i[0] == step ] #i[1] is time and i[2] is uv data, while i[0] is the step number.
                maximum,width = self.get_max(peak)
                if index == 0:
                    p_front[0] = p_tail[0] = maximum * 1000 #improving estimate of the initial guess. this corresponds to the amplitude. 
                angle, angle_front, angle_tail, asymm, params_front, params_tail = self.get_angle(peak,p_front, p_tail,threshold=0.95)
            except Exception:
                peak = [ [i[1], i[2]] for i in corrected_uv_data if i[0] == step-3] #i[1] is time and i[2] is uv data, while i[0] is the step number.
                peak.extend([ [i[1], i[2]] for i in corrected_uv_data if i[0] == step-2])
                peak.extend([ [i[1], i[2]] for i in corrected_uv_data if i[0] == step-1])
                peak.extend([ [i[1], i[2]] for i in corrected_uv_data if i[0] == step])
                maximum,width = self.get_max(peak)
                if index == 0:
                    p_front[0] = p_tail[0] = maximum * 1000 #improving estimate of the initial guess. this corresponds to the amplitude. 
                angle, angle_front, angle_tail, asymm, params_front, params_tail = self.get_angle(peak,p_front, p_tail,threshold=0.95)
            angle_mass_norm = 180 - (standard_mass/resinmass) * (180-angle) 
            integral = self.integrate(peak)
    
            integrals.append((index, self.sequence[-1-index], integral, maximum, width, angle, angle_front, angle_tail,asymm, t_loop_mean, t_reac_mean, resinmass, angle_mass_norm))
        return integrals
    
    def find_params(self):
        """
        THis method finds the file containing additional synthesis parameters such 
        as flow rate, coupling agent, strokes, etc... and converts it to a csv file.
        the file essentially contains the syntheisi parameters to be used for every
        amino acid.
        """
        digits = 6
        path = f"Data/AFPS80_SN{'0'*(digits-len(str(self.serial_no)))}{str(self.serial_no)}.pep"
        for file in os.listdir(path):
            if 'AMINO_PARAMS' in file:
                params = pd.read_csv(os.path.join(path, file))
                return params
    
    def add_params(self, df):
        """
        This method calls the list of parameters for all the amino acids using
        the amino acid in the sequence by iterating through the sequence. In 
        this way it is able to add the synthesis parameters ot the synthesis
        dataframe.
        """
        params = self.find_params()
        params["AA"] = params[' ID AA']
        for i, aa in enumerate(list(df["AA Name"])):
            df.loc[i,['coupling_agent', 'coupling_strokes', 'deprotection_strokes', 'flow_rate']] = params.query(f"AA == '{aa}'").loc[:,['ACT COUPLING', 'NSTROKES DIEA COUPLING', 'NSTROKES DEPRO', 'Q PRIME AND COUPLING']].values[0]
        return df  

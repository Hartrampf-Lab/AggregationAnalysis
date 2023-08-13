from collections import namedtuple
import matplotlib.pyplot as plt
from itertools import chain
import numpy as np
import pandas as pd
import numpy as np
import warnings
warnings.simplefilter("ignore")


class SynthesisPlotter():
    def __init__(self, integrals):
        self.stdblue = "#265952"
        self.stdred = "#DD7965"
        self.stdyellow = "#E1C009"
        self.stdbrown = "#82695C"
        self.stdgreen = "olivedrab"
        self.integrals = integrals

    def plot_integrals(self,save=None):
        fig = plt.figure(figsize=(10, 6))
        num, residue,area, height, width,angle, angle_front, angle_tail, asymm, t_loop_mean, t_reac_mean, resinmass, angle_mass_norm = zip(*self.integrals)
        height_norm = [i/height[0] for i in height] #Normalize to the first deprotection signal
        width_norm = [i/width[0] for i in width]
        area_norm = [i/area[0] for i in area]
        aggregation = np.array(width_norm) - np.array(height_norm)
        x = np.linspace(0, len(aggregation)-1, len(angle),)
        lim_u_blue = max(aggregation)
        lim_l_blue = min(aggregation)
        lim_u_orange = max(angle)
        lim_l_orange = min(angle)

        ax = fig.add_subplot(111)
        ax1 = ax.twinx()
        #ax2 = ax.twinx()
        #ax2.spines.right.set_position(("axes",1.002))
        ax.plot(x, aggregation, label="Aggregation Factor", c=self.stdblue)
        ax1.plot(x, angle, label="Peak Angle", c=self.stdred)
        ax.set_ylabel("Aggregation Factor (w-h)")
        ax1.set_ylabel("Peak Angle")
        #ax2.set_ylabel("Fake Angle")
        x_ticks = ['g' if amino == 'e' else amino for amino in residue]
        ax.set_xticks(num, x_ticks)
        ax.set_xlabel("Residue")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax1.get_legend_handles_labels()
        ax.legend(h1+h2, l1+l2, loc="best", frameon=False)
        plt.tight_layout()
        if save != None:
            plt.savefig(save+"_integralplot.pdf") #best if saved as pdf 
        plt.show()
        
        print("\n\nTemperature profile")
        plt.figure(figsize=(10, 5))
        plt.plot(t_loop_mean)
        plt.ylabel("Activation Loop (in usage) Temperature (° C)")
        plt.xlabel("Residue")
        plt.xticks(num, x_ticks)
        if save != None:
            plt.savefig(save+"_tempplot.pdf") #best if saved as pdf 
        plt.show()
    
    def plot_angle(self, save=None):
        fig = plt.figure(figsize=(10, 6))
        num, residue,area, height, width, angle, angle_front, angle_tail, asymm, t_loop_mean, t_reac_mean, resinmass, angle_mass_norm = zip(*self.integrals)
        ax = fig.add_subplot(111)
        x = np.linspace(0, len(angle)-1, len(angle))
        lim_u_blue = max(chain(angle, angle_mass_norm))
        lim_l_blue = min(chain(angle, angle_mass_norm))
        lim_u_orange = max(chain(angle_front, angle_tail))
        lim_l_orange = min(chain(angle_front ,angle_tail))

        ax1 = ax.twinx()
        ax.plot(x, angle, label="Peak Angle", c=self.stdblue)
        ax.plot(x, angle_mass_norm, label="Mass Normalized Peak Angle", c=self.stdred)
        ax1.plot(x, angle_front, label="Front Angle", c=self.stdgreen)
        ax1.plot(x, angle_tail, label="Tail Angle", c=self.stdyellow)
        ax.set_ylabel("Angle")
        ax1.set_ylabel("Partial Angle")
        x_ticks = ['g' if amino == 'e' else amino for amino in residue]
        ax.set_xticks(num,x_ticks)
        ax.set_xlabel("Residue")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax1.get_legend_handles_labels()
        ax.legend(h1+h2, l1+l2, loc="best", frameon=False)
        plt.tight_layout()
        if save != None:
            plt.savefig(save+"_angleplot.pdf") #best if saved as pdf 
        plt.show()
    
    def plot_4_integrals(self, save=None):
        fig = plt.figure(figsize=(10, 6))
        num, residue,area, height, width, angle, angle_front, angle_tail, asymm, t_loop_mean, t_reac_mean, resinmass, angle_mass_norm = zip(*self.integrals)
        height_norm = [1*i/height[0] for i in height] #Normalize to the first deprotection signal
        width_norm = [1*i/width[0] for i in width]
        area_norm = [1*i/area[0] for i in area]
        
        aggregation = np.array(width_norm) - np.array(height_norm)
        
        ax = fig.add_subplot(111)

        blue_var = aggregation
        orange_var = angle
        yellow_var = height_norm
        purple_var = width_norm
        gray_var = area_norm

        x = np.linspace(0, len(blue_var)-1, len(blue_var),)
        y1 = blue_var
        y2 = orange_var
        y3 = yellow_var
        y4 = purple_var 
        y5 = gray_var
        
        lim_u_blue = max(np.max(y1), np.max(y3), np.max(y4), np.max(y5))
        lim_l_blue = min(np.min(y1), np.min(y3), np.min(y4), np.min(y5))
        lim_u_orange = np.max(orange_var)
        lim_l_orange = np.min(orange_var)
        #lim_u_yellow = 1
        #lim_l_yellow = 0
        ax1 = ax.twinx()
        #ax2 = ax.twinx()
        #ax2.spines.right.set_position(("axes",1.002))
        ax.plot(x, y1, label="Aggregation Factor", c=self.stdblue)
        ax1.plot(x, y2, label="Peak Angle", c=self.stdred)
        ax.plot(x, y3, label="Height (norm)", c=self.stdyellow)
        ax.plot(x, y4, label="Width (norm)", c="purple")
        ax.plot(x, y5, label="Area (norm)", c=self.stdgreen)
        ax.set_ylabel("Param. normalized to 1st peak")
        ax1.set_ylabel("Angle")
        x_ticks = ['g' if amino == 'e' else amino for amino in residue]
        ax.set_xticks(num,x_ticks)
        ax.set_xlabel("Residue")
        #ax.set_ylim(bottom=lim_l_blue,top=lim_u_blue)
        #ax1.set_ylim(bottom=lim_l_orange,top=lim_u_orange)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax1.get_legend_handles_labels()
        ax.legend(h1+h2, l1+l2, loc="best", frameon=False)
        plt.tight_layout()
        if save != None:
            plt.savefig(save+"_4integralplot.pdf") #best if saved as pdf 
        plt.show()
        
        
        print("\n\nTemperature profile")
        plt.figure(figsize=(5, 5))
        plt.plot(t_loop_mean)
        plt.ylabel("Activation Loop (in usage) Temperature (° C)")
        plt.xlabel("Peak index")
        if save != None:
            plt.savefig(save+".pdf") #best if saved as pdf 
        plt.show()
    
    @staticmethod
    def plot_compare(combined, param, serial, sn="all", save=None):
        
        if sn == "all":
            sn = combined["serial_n"].unique()

        plt.figure()
        for n  in sn:
            plt.plot(combined.query(f"serial_n == {n}")["AA number"], combined.query(f"serial_n == {n}")[param],label=f"{n}: {serial[n]}")
        plt.ylabel(param)
        plt.legend(frameon=False, loc="best")
        if save != None:
            plt.savefig(save+f"_{param}_compare.pdf") #best if saved as pdf 
        plt.show()

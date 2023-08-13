###########################################################################
#
#
#   AFPS Aggregation Finding Program - aggregationfinder.py
#   Pietro Luigi Willi, Bálint Tamas, Nina Hartrampf
#   TITLE OF PAPER
#   UZH, Institute of Chemistry.
#   23/06/2023, Zurich.
#
#
###########################################################################


import scipy
from  sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import pandas as pd
from tqdm.notebook import tqdm
import numpy as np

class InvalidSerialNumber(Exception):
        "Raised when the serial number called by sn doesn't exist in serial_n or is not an integer"
        pass

class AggregationFinder:
    """
    This data analysis tool was developed for finding and identifying
    points of aggregation based on the peak angle, i.e. the angle of
    the deprotection peak. The tool provides two different readouts
    for aggregation, i.e. two different methods for finding
    aggregation. One method is based on summing the rows of the
    slope matrix. Here the point of aggregation is definend as the
    argmax (index of max) of the cumulative gradient sum while the
    maximum gradient is the extent of aggregation. It has been
    observed that a gradient value higher than 10 is usually
    consistent with aggregation. The other method is based on fitting
    a sigmoid on the angle trace:

    (a/(1+exp(d*x-b)))+c

    The fitting provides the 4 parameters of the sigmoid, which we
    can use to determine where the aggregation will be happening
    (b/d), how much (a) and how fast (d) it is happening. In order
    for the two methods not to be confused by sharp peaks introduced
    by temperature anomalies, any peaks resulting from deprotections
    performed with activation temperatures 20% higher or lower than
    the local median (median of tempertures within the synthesis)
    will be trimmed. The trimming happens by setting the peak angle
    to the average between the next and previous peak angle.

    --------------------------------------
    Attributes:
    df : pd.DataFrame
        the dataframe to operate on.
        Must include columns param, "serial_n", "AA Name", "AA number", "t_loop"
    param : str
        the parameter onto which the aggregation computations will performed.
        This must be the name of one of the columns in df.


    --------------------------------------
    Usage:
    - Define an AggregationFinder instance. This must include a df and a param
      as shown below. In this example the peak angle will be used as the
      target of the aggreggation computations. It is strongly advised to use
      the Peak Angle. However a different parameter such as "first_diff"
      (aggregation factor) or "Area" could also be used:
    >>> param = "Peak Angle"
    >>> aggr_finder = AggregationFinder(df,param)

    - Use the find_aggregation method to compute the sigmoid fitting and the
      cumulative slope. This method will return the DataFrame which you
      provided as attribute concatenated with the gradient and sigmoid prediction,
      aswell as a new DataFrame containing the parameters of the sigmoid and the
      max&argmax of the cumulative slope. You can recycle the df if you want
      to change param, but do not recycle the DataFrame containing the predicted
      paramameters as it will be overwritten. On the find_aggregation method you
      can specify the trimming interval using the interval parameter. In other
      words, you can specify how much percentage difference from the median the
      deprotection temperture (t_loop) needs to be for the peak to be trimmed.
      This is set to 0.2 (20%) as default.
    >>> df, aggr_angle_df = aggr_finder.find_aggregation(interval=0.2)

    - Finally you can plot the aggregation prediction of the 2 methods using the
      plot_aggregation method. here you need to specify the serial number as
      parameter sn. This is necessary to distinguish between the various sequences
      in the dataset. Even if you may only have one sequence, you still need to
      specify it. the method will plot two subplots. In the upper you will see the
      trimmed param in red and the fitted sigmoid in blue. On this subplot you will also
      see a brown vertical line corresponding to the prediction of aggregation of the
      cumulative slope method and a yellow vertical line corresponding to the prediction of
      aggregation of the sigmoid method. In the bottom plot you will see the 
      cumulative slope in blue and the same brown line. The parameters returned
      by the two aggregation finidng methods are printed below the plot but can also
      be found on the aforementioned DataFrame. You can also save this figure by
      specifying the name of the file as a string followed by the .pdf or .jpg suffix:
    >>> aggr_finder.plot_aggregation(sn=sn, save="aggr_plot.pdf")
    """

    def __init__(self, df, param):
        self.df = df
        self.param = param
        self.stdblue = "#265952"
        self.stdred = "#DD7965"
        self.stdyellow = "#E1C009"
        self.stdbrown = "#82695C"

    def accumulate_slope(self, sn, suppl=" trimmed"):
        """
        Returns the row-wise sum of the slope matrix as a numpy array.
        The slope matrix is the slope of every point with respect to 
        every other point.
        """

        seq_i = list(np.argwhere(self.df["serial_n"].to_numpy()==sn).flatten())
        seq = self.df.loc[seq_i, str(self.param)+suppl]
        x = np.linspace(0, len(seq)-1, len(seq))
        y = seq.to_numpy()
        grad = []
        for i in range(len(x)):
            theta = sum(((y-y[i])/(x-x[i]+10**-22)))#sums the rows of the slope matrix, in the denominator we add a very small number to avoid divison by 0.
            grad.append(theta)
        grad = np.array(grad) / len(grad)
        return grad

    def sigmoid(self, x, a, b, c, d):
        return (a/(1+np.exp(d*x-b))) + c

    def fit_sigmoid(self, seq, func, suppl=" trimmed"):
        """
        Fits the sigmoid onto the peak angle and returns the
        y values of the sigmoid function, its parameters and
        the R2 score. The categorization of aggregation works
        best if the sigmoid function is fitted onto the trimmed
        parameter, thus the suppl parameter defaults to trimmed.
        """

        y = np.array(seq[str(self.param)+suppl])
        #the parameters below represent inital guesse which are used to reduce the runtime of the optimization algorithm
        a_i = min(y)-max(y)
        d_i = y[np.argmax(y)] - y[np.argmax(y)-1]
        b_i = np.argmax(y) * d_i
        c_i = max(y)
        p0 = np.array([a_i, b_i, c_i, d_i]) # inital guess
        x = np.linspace(0, len(y)-1, len(y))

        try:
            params = scipy.optimize.curve_fit(func, x, y, p0=p0, maxfev=1000000)[0] #fits the parameters of the sigmoid using the initial guesses
            self.params = params
            #maxfev determines the maximum calls of the function (max optim. iterations)
            #other parameters of interest to input as kwargs:
            #epsfcn: this is the step length of the forward difference approximation
            pred = func(x, params[0], params[1], params[2], params[3])
            r2 = r2_score(y, pred,)
        except RuntimeError:
            print("Max optimizer calls reached: Runtime Error")
            params = [None, None, None, None]
            self.params = params
            pred = x.fill(None)
        return params, pred, r2
    
    def add_sequence(self):
        """
        Makes two new rows. One containing the whole sequence
        of the peptide being synthesised on every row and one
        containig the sequence of the growing peptide chain.
        """

        seq = []
        seq_cumul = []
        sn_list = self.df.serial_n.unique()
        for sn in tqdm(sn_list):
            seq_n = self.df.query(f"serial_n == {sn}")
            seq_iter = "".join(list(seq_n["AA Name"]))
            for i in seq_n.index:
                seq.append(seq_iter)
                seq_cumul.append("".join(list(seq_n.loc[:i,"AA Name"])))
        self.df["sequence"] = seq
        self.df["chain"] = seq_cumul
        return self.df

    def trim_peak(self, interval):
        """
        In this function the sharp peaks caused by large temperture differences will be
        flattened in order to make the detection and categorization of aggregation
        easier. The flattening will be done by finding all the temperature anomalies.
        temperature anomalies are deprotection temperture of the activation loop  which
        are outside of the 20% interval of the median (20% higher or lower than the median)
        and their angle value will be replaced by the mean between the angle of the
        previous and following residue.

        interval: specifies how far from the average a peak temperature needs to be
        in order to be trimmed. for example if the interval is 0.2, then peaks with
        temperature readings 20% higher or lower than the average are trimmed.
        """

        self.df[f"{self.param} trimmed"] = self.df[self.param]
        for sn in tqdm(pd.unique(self.df["serial_n"])):
            local_df = self.df.query(f"serial_n == {sn}")
            local_t_median = local_df["t_loop"].median()
            anomaly_i_list = local_df.query(f"t_loop >= {local_t_median}*{1+interval} or t_loop <= {local_t_median}*{1-interval}").index
            for i in anomaly_i_list:
                if self.df.loc[i, "AA number"] == 0:
                    self.df.loc[i, f"{self.param} trimmed"] = self.df.loc[i+1, self.param] #if the C/H residue is the first in the sequence, sets its angle equal to the next residue
                elif self.df.loc[i, "AA number"] == len(local_df) - 1:
                    self.df.loc[i, f"{self.param} trimmed"] = self.df.loc[i-1, self.param] #if the C/H residue is the first in the sequence, sets its angle equal to the first residue#if the C/H residue is the last in the sequence, set its angle equal to the previous residue
                else:
                    self.df.loc[i, f"{self.param} trimmed"] = (self.df.loc[i-1, self.param] + self.df.loc[i+1, self.param])/2 #takes the mean between next and previous residue param.
        return self.df

    def find_aggregation(self,interval=0.2):
        """
        This function applies both the cumulative slope method and the sigmoid method
        to detect the start of aggregation and its extent.
        """

        self.df = self.trim_peak(interval=interval)
        self.df = self.add_sequence()
        self.aggregation = pd.DataFrame(columns=("serial_n", 
                                                "sequence", 
                                                "aggregation_index", 
                                                "normalized_gradient", 
                                                "sigmoid_a", 
                                                "sigmoid_b", 
                                                "sigmoid_c", 
                                                "sigmoid_d", 
                                                "r2", 
                                                "max", 
                                                "argmax",
                                                "max_temp",
                                                "min_temp",
                                                "med_temp", 
                                                "max_coupling", 
                                                "min_coupling", 
                                                "med_coupling", 
                                                "max_depro", 
                                                "min_depro", 
                                                "med_depro", 
                                                "max_flowrate", 
                                                "min_flowrate", 
                                                "med_flowrate")) #makes a new dataset with aggregation parameters based on self.param for all the sequences in self.df.
        gradient = []
        sigline = []
        
        for sn in tqdm(pd.unique(self.df["serial_n"])):
            local_df = self.df.query(f"serial_n == {sn}")#isolates individual sequences by querying their serial number
            grad = self.accumulate_slope(sn)
            g_max = max(grad)
            g_argmax = np.argmax(grad) #finds the index at which the gradient is largest: aggregation
            if len(local_df) > 3:
                par_sig,sig,r2 = self.fit_sigmoid(local_df,self.sigmoid) #fits a sigmoid on every sequence in the dataset as long as it has more than 3 residues
                self.aggregation.loc[len(self.aggregation)] = ([sn,
                                                                local_df["sequence"].unique()[0],
                                                                g_argmax,
                                                                g_max,
                                                                par_sig[0],
                                                                par_sig[1],
                                                                par_sig[2],
                                                                par_sig[3],
                                                                r2,
                                                                max(local_df[self.param]),
                                                                np.argmax(local_df[self.param]),
                                                                max(local_df["t_loop"]),
                                                                min(local_df["t_loop"]),
                                                                np.median(local_df["t_loop"]),
                                                                max(local_df["coupling_strokes"]),
                                                                min(local_df["coupling_strokes"]),
                                                                np.median(local_df["coupling_strokes"]),
                                                                max(local_df["deprotection_strokes"]),
                                                                min(local_df["deprotection_strokes"]),
                                                                np.median(local_df["deprotection_strokes"]),
                                                                max(local_df["flow_rate"]),
                                                                min(local_df["flow_rate"]),
                                                                np.median(local_df["flow_rate"])])
            else:
                sig = np.zeros(len(local_df))
                self.aggregation.loc[len(self.aggregation)] = ([sn,
                                                                local_df["sequence"].unique()[0],
                                                                g_argmax,
                                                                g_max,
                                                                None,
                                                                None,
                                                                None,
                                                                None,
                                                                None,
                                                                max(local_df[self.param]),
                                                                np.argmax(local_df[self.param]),
                                                                max(local_df["t_loop"]),
                                                                min(local_df["t_loop"]),
                                                                np.median(local_df["t_loop"]),
                                                                max(local_df["coupling_strokes"]),
                                                                min(local_df["coupling_strokes"]),
                                                                np.median(local_df["coupling_strokes"]),
                                                                max(local_df["deprotection_strokes"]),
                                                                min(local_df["deprotection_strokes"]),
                                                                np.median(local_df["deprotection_strokes"]),
                                                                max(local_df["flow_rate"]),
                                                                min(local_df["flow_rate"]),
                                                                np.median(local_df["flow_rate"])])
            gradient.append(grad)
            sigline.append(sig)
        gradient = np.concatenate(gradient)
        sigline = np.concatenate(sigline)
        self.df[f"{self.param} Slope"] = gradient
        self.df[f"{self.param} Sigmoid"] = sigline
        return self.df, self.aggregation

    def plot_aggregation(self, sn, save=None):
        if sn not in self.df["serial_n"].unique() or type(sn) != int:
            print("The serial number passed is not valid or doesn't exist")
            raise InvalidSerialNumber
        local_df = self.df.query(f"serial_n == {sn}")
        if len(self.aggregation) > 1:
            local_params = self.aggregation.query(f"serial_n == {sn}")
        else:
            local_params = self.aggregation
        local_df["index"] = list(range(len(local_df)))
        local_df = local_df.set_index("index")
        local_params["sigmoid index"] = local_params["sigmoid_b"] / local_params["sigmoid_d"]
        x = np.linspace(0, len(local_df)-1, len(local_df))
        x_s = np.linspace(0, len(local_df)-1, 1000)
        sigmoid_pred = self.sigmoid(x_s, 
                                    local_params["sigmoid_a"].values, 
                                    local_params["sigmoid_b"].values, 
                                    local_params["sigmoid_c"].values, 
                                    local_params["sigmoid_d"].values)
        fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(7,5))
        ax1.set_title(f"Sequence {sn}")
        
        ax1.plot(x, local_df[self.param], label=self.param, c="grey", alpha=0.5)
        ax1.plot(x, local_df[f"{self.param} trimmed"], label=f"{self.param} trimmed", c="black")
        ax1.plot(x_s, sigmoid_pred, label="Sigmoid Fit", c=self.stdred)
        ax1.set_ylabel(f"{self.param}")
        ax1.set_xlabel("Residue Index")
        #ax1.axvline(x=float(local_params["aggregation_index"]), color=self.stdbrown, linestyle="--")
        if float(local_params["sigmoid index"]) > 0 and float(local_params["sigmoid index"]) < len(local_df):
            ax1.axvline(x=float(local_params["sigmoid index"]), color=self.stdyellow, linestyle=":", label="Sigmoid aggregation")
        ax1.set_ylim(min(local_df[self.param])-0.5, max(local_df[self.param])+0.5)
        ax1.legend(frameon=False)
        ax3 = ax2.twinx()
        ax2.plot(x, local_df[f"{self.param}"], label=f"{self.param}", c="grey", alpha=0.5)
        ax2.plot(x, local_df[f"{self.param} trimmed"], label=f"{self.param} trimmed", c="black")
        ax3.plot(x, local_df[f"{self.param} Slope"], label=f"Cumulative slope of sn: {sn}", c=self.stdred)
        ax3.set_ylabel(f"Cumulative {self.param} Slope")
        ax2.set_ylabel(f"{self.param}")
        ax2.set_xlabel("Residue Index")
        ax2.axvline(x=float(local_params["aggregation_index"]), color=self.stdbrown, linestyle="--", label="Aggregation index")
        ax3.set_ylim(min(local_df[f"{self.param} Slope"])-0.5, max(local_df[f"{self.param} Slope"])+0.5)
        ax2.legend(frameon=False)
        plt.tight_layout()
        

        if float(local_params["sigmoid index"]) > 0 and float(local_params["sigmoid index"]) < len(local_df):
            over = ''
        else:
            over = '(out of range)'

        print(f"\nCumulative Slope Aggregation Prediction:\nIndex of aggregation = {float(local_params['aggregation_index'])}\nMagnitude of aggregation = {float(local_params['normalized_gradient']):0.2f}\n\nSigmoid Aggregation Prediction:\nx-value of inflection point = {float(local_params['sigmoid index']):0.2f}\t   {over}\nAmplitude of sigmoid = {float(local_params['sigmoid_a']):0.2f}\t  (magnitude of aggregation)\nGradient at inflection point = {float(local_params['sigmoid_d']):0.2f}\nR2 = {float(local_params['r2']):0.4f}")

        if save != None:
            plt.savefig(save) #best if saved as pdf
        plt.show()
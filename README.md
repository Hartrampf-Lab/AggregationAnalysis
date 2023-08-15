# Aggregation Analysis
Aggregation Analysis is a tool for investigating and detecting aggregation for the Automated Fast-Flow Peptide Synthesiser. 
It contains the code to implement the methods and concepts introduced in our paper:
<br>
<br>
**A robust analytical method to investigate sequence dependence in flow-based peptide synthesis** 
<br>
(doi: tba) 
<br>
## Conda environment
A conda environment YML file is included in order to successfully run the code.
Create the dataanlysis environment by running the following command: 
<br>
```bash
conda env create -f environment.yml
```
In order to be able to use the environement kernel on jupyter notebook, run the following command:
```bash
conda activate dataanalysis
python -m ipykernel install --user --name dataanalysis
```

## Description of the data processing and analysis method.
The description is taken from the supporting information of the paper "A robust analytical method to investigate sequence dependence in flow-based peptide synthesis", which can be found on doi.org....

The repository is composed of four python scripts: The main script
(**afps_dataiter.py**) iterates through the files containing the
synthesis data (.pep files) and extracts the synthesis data. The output
of the main script is plotted using the **plotting.py** script. Using
the areas computed by the **afps_dataiter.py** script, the
**permutator.py** script combinatorically iterates through all the
deletion permutations and calculates their likelihood based on the
areas. The **aggregationfinder.py** script uses the output of the main
script to find the position and magnitude of aggregation based on the
change in peak angle. The three scripts all contain a class called by an
IPython notebook, **AFPSDataExport.ipynb**.

### **Peak isolation and computation of area, height, and width**

The afps_dataiter contains the synthesis class, which consists of
various methods enabling it to open the raw synthesis (.pep) files,
isolate the deprotection peaks and then calculate a set of summary
statistics. First, the baseline is corrected by subtracting a fitted
third-order polynomial. The deprotection peaks are then isolated by
querying the steps that correspond to deprotection. New deprotection
starts every 22 steps and the first deprotection step occurs at step 36.
The summary statistics that are then calculated from the deprotection
peaks are area, maximum height, full width at half maximum and angle.
The area is calculated by taking the Riemann sum under the peak using
the following formula:

$$a = \sum_{}^{}\frac{0.5(A_{n + 1} + A_{n})}{t_{n + 1} - t_{n}}\ $$

Where *a* stands for the area under the peak, *A* for absorption and *t*
for time (ms). The height of the peak is simply computed by just taking
the maximum value of the peak. The width is found by calculating the
difference in time between the two half maxima at either side of the
peak.

### **Angle computation**

The angle computation differs in complexity from the previous summary
statistics in 2.1. Its complexity is derived from the fitting of a
Gaussian function onto the peak and from the need to increase the
robustness of the method. The whole angle computation is split into two
main methods, process_peak and get_angle. The process_peak method is
given the raw deprotection peak as input. The peak is first split into
front and tail at its maximum by calling the static method split_peak. A
threshold is set so that any absorbance values larger than the threshold
(0.90 AU) are set to the threshold + 0.001. This ensures that all the
oversaturated parts of the peak have the same values. Then the middle of
the peak is found by finding the median time of the maximum absorbance.
The thresholding makes it easier to find the middle of the peak in case
of oversaturation. The peak is then split into front and tail along the
median of the maxima. Using the static method mirror_peak, the two half
peaks are then individually mirrored and concatenated with the original
half peak, resulting in two full Gaussian-like peaks. Finally, the
static method trim_peak, removes the oversaturated part of the peak by
simply deleting the points with an absorbance higher than the threshold
by their index. In addition, it also trims the peaks at their minima,
making the peak more Gaussian-like. These three static methods are all
called by the process_peak method in the order described to yield the
coordinates of the front and tail peaks.

The get_angle method then calls the process_peak method and fits a
parametrized Gaussian function with the following formula:

$$A = ae^{\frac{- {(t - b)}^{2}}{{2c}^{2}}}$$

Where *t* stands for time and *A* for absorption. The curve fitting
yields the parameters *a, b* and *c* for the front and the tail peak.
The parameters are then used to find the maximum gradient of the
function *G*. The maximum gradient corresponds to the gradient at the
point of inflection. The gradient of the function *G* is defined by its
first order derivative, $\frac{\text{dA}}{\text{dt}}$:

$$\frac{\text{dA}}{\text{dt}} = \frac{- (t - b)}{c^{2}}ae^{\frac{- {(t - b)}^{2}}{{2c}^{2}}}$$

The arctangent of one over the gradient at the inflection will yield the
angle at the top of the peak. After performing the angle computation
both on the front and tail peak, the two resultant angles are summed
together to yield the peak angle using the following formulas:

$$\alpha_{\text{front}} = \arctan\left( \frac{1}{max\left( \frac{\text{dA}}{\text{dt}}_{\text{front}} \right)} \right)$$

$$\alpha_{\text{tail}} = \arctan\left( \frac{1}{max\left( \frac{\text{dA}}{\text{dt}}_{\text{tail}} \right)} \right)$$

$$\alpha = \ \alpha_{\text{front}} + \alpha_{\text{tail}}$$

The get_integrals method iterates through the steps that are defined as
deprotections and computes the various summary statistics.

### **Angle normalization using the mass**

The angle is standardized using the resin mass. The mass normalization
is performed to remove the peak angle’s dependence on resin mass (but
not resin loading). The standardization normalizes the peak angle to a
defined standard resin mass of *m<sub>st</sub>* = 150 mg using the
following formula:

$$\alpha_{\text{st}} = 180{^\circ} - \frac{m_{\text{st}}}{m}(180{^\circ} - \alpha)$$

Where *m* is mass of the resin used during the synthesis in mg,
*m<sub>st</sub>* the standard resin mass in mg, *α* the measured angle
and *α<sub>st</sub>* is the standardized angle.

### **Trimming of outlying temperatures**

The trim_peak method performs a temperature correction. To reduce
outliers, points belonging to deprotections that have been performed
outside of a +/-20% temperature range compared to the median temperature
of the synthesis are set to the average of the previous and next points
carried out within the temperature range. If a temperature anomaly
occurs as the first or last deprotection this point is set equal to only
one neighbor.

### **Aggregation Detection**

After the listed descriptors of all deprotection peaks have been
calculated and added to a pandas dataframe, an object of the
AggregationFinder class contained by the aggregationfinder.py script is
instantiated using the dataframe created by the Synthesis class and the
name of the parameter onto which the method of the AggregationFinder
class is to be applied. The methods were designed to be used with the
peak angle, but they were also designed to work with other parameters.
This class executes the automated detection of the position and
magnitude of the aggregation, based on the detection of a sudden
permanent increase of the peak angle. Two different methods of
aggregation detection are applied. The cumulative slope method computes
the average slope of the summary statistic trace between a point and all
the other point for every point, using the following formula:

$$C_{k} = \ \frac{\sum_{\begin{matrix}
i = 1 \\
i \neq k \\
\end{matrix}}^{n}\frac{y_{i} - y_{k}}{x_{i} - x_{k}}}{n}$$

Where *y<sub>k</sub>* is the angle belonging to the peak of choice,
*y<sub>i</sub>* the angle belonging to every other peak, *x<sub>k</sub>*
is the position of the peptide in a reversed CàN peptide, and
*x<sub>i</sub>* is the position of all the other amino acids. n is the
length of the synthesized peptide, to normalize amongst different
peptide lengths. The maximum of the cumulative slope corresponds to the
position and magnitude of aggregation for a particular synthesis.

The other method fits a parametrized sigmoid curve onto the summary
statistic trace using the following formula and uses the point of
inflection of the sigmoid as the position of aggregation.

$$\sigma(x) = \frac{a}{1 + e^{dx - b}} + c$$

The $\frac{a}{c}$ indicates the magnitude of aggregation and
$\frac{b}{d}$ the position of it.

The cumulative slope is calculated by the accumulate_slope method,
yielding the gradient array, while the sigmoid is fit by the fit_sigmoid
method.

The find_aggregation method then iterates through a dataframe containing
the synthesis summary statistics for one or multiple syntheses. It
differentiates between different syntheses by querying them by their
serial numbers sequentially. It first calls the trim_peak method, then
the add_sequence method. The add_sequence method computes the full
sequence and the growing peptide chain sequence of the synthesis. By
iterating through the unique serial numbers of the syntheses, the
add_sigmoid method and the accumulate_slope method are clled on each
synthesis individually.

The afps_dataiter script, more specifically the parts performing the
parsing of the .pep files were adapted from the
literature.<sup>SI1,2</sup>

###  **Loading determination**

Using the areas of the deprotection peaks, the loading of the solid
support can be determined. Loading, *L*, is defined as mmol of binding
points, *n*, per mass unit of resin, *m*. In most techniques to
determine the loading, the number of binding points are approximated
with the amount of Fmoc removed after the first coupling. Using this
assumption, the determination is done using the following formula:

$$L = \frac{n}{m}\ $$

$$n = \ \dot{V}\int_{t_{\text{start}}}^{t_{\text{end}}}c\text{dt}$$

$$A = l\varepsilon c$$

Combining these formulas results in:

$$L = \ \frac{\dot{V}}{\text{mlε}}\int_{t_{\text{start}}}^{t_{\text{end}}}A\text{dt}$$

In this final equation $\dot{V}$ stands for flowrate, *l* for part
length in the flow cell of the UV-detector, ε is the extinction
coefficient of the deprotection solution and the integral of the
absorption (*A*) over time is the area of the peak. Using this simple
formula, we can automatically determine the loading of the resins with
ease, eliminating the cumbersome extra step its measurement required.

### **Deletion Computation using peak integrals**

The deletion probability calculation assumes that the first amino acid
was coupled with 100% coupling efficiency. Under this assumption we can
calculate the coupling efficiency of all the other amino acids by
dividing their deprotection peak area by the area of the first
deprotection. The different deletion permutations are enumerated by
using the combinatorial property of binary numbers. The number of
different deletion combinations is equal to (n - 1)<sup>2</sup> – 1,
where n is the length of the sequence. By counting down from the maximum
number of combinations in binary, all the various deletion permutations
are enumerated. In the binary number the 0 digits correspond to a
deletion at that particular position. This same property that is used to
enumerate the permutations is also used to calculate the probability.
The probability of a particular deletion permutation can be calculated
using the following formula:

$$P = \prod_{i = 0}^{n}\left( b_{n}A_{n} + b_{n}^{- 1}A_{n}^{- 1} \right)$$

Where *b* is the binary array that represents a particular deletion
permutation, *b<sub>n</sub>* is the digit at position *n* of the binary
array, *A<sub>n</sub>* is the area of the n-th deprotection normalized
by the area of the first deprotection, $A_{n}^{- 1} = 1 - A_{n}$, and
$b_{n}^{- 1}$ is the ones’ complement of the binary array with the first
digit always set to 1 as it corresponds to the first deprotection. This
is all implemented in the DeletionPermutator class of the permutator.py
script, which is provided with the pandas dataframe containing the
processed data from the synthesis class and the serial number of the
synthesis of interest. The class has three methods: the
convert_to_binary_with_leading_zeros static method converts numbers into
binary with leading zeros and one in the first position, the
permute_sequence method iterates enumerates the deletion permutations by
calling the convert_to_binary_with_leading_zeros and the compute_mass
calculates the mass for every deletion permutation.

### **Comparison of the efficiency of data analytical methods for oversaturated signals** 

To investigate the capability of the data analytical methods in
capturing aggregation in signals that have oversaturated peaks we
generated the oversaturation *in silico*. Taking the standard synthesis
(150 mg, normal loading, barstar\[75–90\] synthesis) we gradually cut
down the signal. Once the value of the cutting down reached the tallest
peak, that is where oversaturation was considered 0%. For every value of
oversaturation, we analyzed the traces using both the peak angle and
aggregation factor and compared the similarity of the plots to the
original analysis using R<sup>2</sup>. An example for this can be seen
in **Figure S1** at approximately 50% oversaturation.

![**Figure S1A)** An example for the artificially generated
oversaturated plot at the UV absorption value of 0.4, which corresponds
to an approximate 50% oversaturation, compared to the original plot](https://github.com/Hartrampf-Lab/AggregationAnalysis/edit/main/plots/artificial_oversaturation.png?raw=true)
![**Figure S1B)** The analysis of the oversaturated plot using both peak angle and
aggregation factor. It is clearly visible that while peak angle
maintains its characteristics and covers the aggregation point
accurately aggregation factor loses accuracy. By calculating and
plotting the R2 we can compare the capability of the two methods in
handling oversaturation
](https://github.com/Hartrampf-Lab/AggregationAnalysis/edit/main/plots/oversaturation_PA_AF_trace.png?raw=true)


## Requirements
* <a href='https://www.python.org/downloads/release/python-3110/'>Python 3.11</a>
* <a href='https://www.rdkit.org/'>RDKit</a>


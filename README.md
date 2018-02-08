# PisCES
## Introduction
This is code of paper "Global Spectral clustering in dynamic networks" (Fuchen Liu, David Choi, Lu Xie and Kathryn Roeder, 2018, http://www.pnas.org/content/early/2018/01/10/1718449115).  Codes are included for the main algorithm--PisCES, and for our visualization tool--sankey plots.

In this paper, we build a global spectral clustering method to perform community detection in dynamic networks (or network series).We implement degreecorrected spectral clustering, with a smoothing term to promote similarity across time periods, and iterate until a fixed point is achieved. Specifically, this global spectral clustering approach combines the current network with the leading eigenvector of both the previous and future results. The combination is formed as an optimization problem that can be solved globally under moderate levels of smoothing when the number of communities is known. We also utilized data-driven method to choose appropriate levels of both smoothing and model order, as well as to balance regularization
with “letting the data speak". After finding the community assignments, we use sankey plots to show the flows (or changes) of communities over time.

<p align="center">
 <img src="graphs/T1toT6.png" width="370"/> 
</p>

## Usage
After extracting a list of adjacency (or correlation) matrices from your networks, save them as a 3-dimension array A with dimensions as NNT, where N is the number of nodes and T is the number of times. Then run PisCES to get a list of memberships. More detailed usage and options can be found in code folder: PisCES_code.
```
memberships = PisCES(A)
```
After getting the memberships, we also provide the code for sankey plot to visualize the results. A detailed script can be found in code folder: sankey_plot.

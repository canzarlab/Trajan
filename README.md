# Trajan #

### Alignment of complex trajectories from single-cell RNAseq experiments. ###

![Trajan Overview](img/trajan_overview.png)

## Installing Trajan ##

All dependencies are bundled with Trajan. In particular, Trajan implements its own non-linear solver and does 
not rely on any external (I)LP solver. To ```convert``` different input formats accepted by Trajan you need 
to install the Boost graph library. To build Trajan, simply run

```
make
```

## Running Trajan ##

### To begin ###

First you will need to obtain a pair of trajectories using any with any of the available trajectory-building techniques (Saelens et al. 2018), over 50 trajectory inference methods have been developed since 2014. Trajan does not make any assumptions on the type of methods used to build your trajectory, one should take of converting the output in into a suitable input format for Trajan. In our examplary workflow we will describe how to transform the output obtained from Monocle (Trapnell, C. et al. 2014) into a suitable input for Trajan.

In general one should provide for each trajectory, the following files:
1. __Edge Set__: t.tree
2. __Map__: t.map

additionally a Distance Matrix between all pairs of edges:
  - __Cost Matrix__: cost_matrix_0_1.csv

A typical input Folder should look like this: 
  - t_0.map
  - t_1.map
  - t_0.tree
  - t_1.tree
  - cost_matrix_0_1.csv  


### Usage ###

Once compiled. Trajan can be run easily with the following command:

```
trajan <tree_1> <map_1> <tree_1> <map_1> <cost_matrix> <align> <weightfunc> <k> <vareps> <coneps> <solver>
```

For our typical input:

```
trajan input/t_0.tree input/t_0.map input/t_1.tree input/t_1.map input/cost_matrix_0_1.csv results/output_0_1.csv 2 e 0 0.02 0.01 2 
```

#### Input/Output formats

Input: 

`<tree>`
  : input tree file in the following format (use convert.cpp to convert from .dot):        
  [child node] [parent node] default [newline] ...

Output:

`<align>`
  : output the final alignment in file in the following format:
  [node in first graph] [node in second graph] [fractional solution] [newline] ...


## Example Workflow using Monocle generated Trajectories ## 

TODO: 
<!---
Differentiation vs Reprogramming (Cacchiarelli, D. & Trapnell, C. et al. 2018)
-->

## References

<div id="refs" class="references">


<div id="ref-trapnell:2014">
  
The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells
Cole Trapnell  et al. Nature Biotechnology volume 32, pages 381–386 (2014)
<https://doi.org/10.1038/nbt.2859>

</div>


<div id="ref-trapnell:2018">

Aligning Single-Cell Developmental and Reprogramming Trajectories Identifies Molecular Determinants of Myogenic Reprogramming Outcome.
Cacchiarelli, D. & Trapnell, C. et al. Cell Syst. 7(3), 258-268 (2018).
<https://doi.org/10.1016/j.cels.2018.07.006>.

</div>

<div id="ref-alpert:2018">

Alignment of single-cell trajectories to compare cellular expression dynamics
Ayelet Alpert, Lindsay S Moore, Tania Dubovik & Shai S Shen-Orr.
Nature Methods volume 15, pages 267–270 (2018)
<https://doi.org/10.1038/nmeth.4628>.

</div>

<div id="ref-Saelens:2018">

A comparison of single-cell trajectory inference methods: towards more accurate and robust tools
Wouter Saelens, Robrecht Cannoodt, Helena Todorov, Yvan Saeys.
bioRxiv (2018)
<https://doi.org/10.1101/276907>.

</div>


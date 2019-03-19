# HINTRA
An algorithm for collaborative intra-tumor heterogeneity detection


## Abstract
Despite the remarkable advances in sequencing and computational techniques, noise in the data and complexity of the underlying biological mechanisms render deconvolution of the phylogenetic relationships between cancer mutations difficult. Besides that, the majority of the existing datasets consist of bulk sequencing data of single tumor sample of an individual. Accurate inference of the phylogenetic order of mutations is particularly challenging in these cases and the existing methods are faced with several theoretical limitations. To overcome these limitations, new methods are required for integrating and harnessing the full potential of the existing data.

We introduce a method called Hintra for intra-tumor heterogeneity detection. Hintra integrates sequencing data for a cohort of tumors and infers tumor phylogeny for each individual based on the evolutionary information shared between different tumors. Through an iterative process, Hintra learns the repeating evolutionary patterns and uses this information for resolving the phylogenetic ambiguities of individual tumors. The results of synthetic experiments show an improved performance compared to two state-of-the-art methods. The experimental results with a recent Breast Cancer dataset are consistent with the existing knowledge and provide potentially interesting findings.


## Inputs
Assume that the name of the dataset is `example`. Then, the input files include:

* `example.Rcounts`: This file contains reference read counts for tumor samples. Each line of this text file corresponds to a sample and is a sequence of integer values separated by tabs. Each value corresponds to the reference read count of a gene. The order of genes is consistent for all samples (lines). This file should not contain any zeros. For the variants that does not exist in a sample, the reference read count can be set to an arbitrary number (e.g. 100). An example file for 5 samples and 4 genes is provided below:

```
70  100 30  150
100 100 51 90
95 82 75  59
0 39 100 100
55 201 100 100
```

* `example.Vcounts`: This file contains variant read counts for tumor samples. Each line of this text file corresponds to a sample and is a sequence of integer values separated by tabs. Each value corresponds to the variant read count of a gene. The order of genes is consistent for all samples (lines). Below is an example corresponding to the above `example.Rcounts` file:

```
10  0 25  50
0 0 49 60
25 41 11  6
0 35 0 0
55 35 13 0
```

* `example.Tids`: This file contains tumor IDs for samples. Each line of this text file contains an integer positive number corresponding to the ID of a tumor from which a sample is taken. The order of IDs should follow the same order as samples appear in the `example.Rcounts` and `example.Vcounts` files. The IDs should be starting from 1 and should be consecutive. In other words, the largest ID number is equal to the number of tumors. The samples should be ordered in the `example.Rcounts` and `example.Vcounts` files such that samples from each tumor are beside each other and they are sorted by the tumor ID number. Below is an example corresponding to the above `example.Rcounts` file:

```
1
1
2
2
3
```

* The folder `PhylogenySet`: This folder contains the set of all valid phylogenetic structures (i.e. structures without branching root node) for different numbers of mutations. It contains files named as `ParentVectors_K-?.txt`, where `?` indicates the number of mutations. Currently, this folder contains files for between 1 and 6 mutations. Larger topologies can be generated using the provided R code `TopologyEnumerator.R`. Each file contains one or more lines with each line representing a phylogenetic tree in the parent vector format (i.e. a sequence of numbers each indicating the parent of the corresponding node with 0 being the root node).

* `example.margs` (if applicable): This file contains pre-processed marginal likelihood data for tumors (not samples). This will be used if available. Otherwise, Hintra will produce this file as described later. Each line of this tab-separated file is dedicated to one tumor and is a sequence of floating values. Each value is a marginal log-likelihood of a tree topology. The order of values follows the order of topologies provided in the folder `PhylogenySet`. The order of lines follows the order of tumor IDs (1, 2,...).

* `example.maxs` (if applicable): Similar to `example.margs`, but contains maximal log-likelihoods.


## Outputs
For the given inputs for a dataset named `example`, Hintra produces the following outputs:

* `example.margs` (if applicable): If not already existing, or if the command line argument `-u` is used as described in section **Usage**, this file will be computed. See the **Inputs** section above for the description of contents.

* `example.maxs` (if applicable): If not already existing, or if the command line argument `-u` is used as described in section **Usage**, this file will be computed. See the **Inputs** section above for the description of contents.

*Note: If either of the files `example.margs` or `example.maxs` does not exist, both will be regenarated by Hintra*

* `example.trees`: This tab-separated file contains the detected phylogenetic trees. Each line is dedicated to one tumor in the same order as provided in the input. For each tumor, the file shows the corresponding phylogenetic tree in the parent vector format. The values of vector elements indicate the index of the parent gene (following the input order). The germline root node is referred to by `0` and non-mutated genes are indicated by `-1`. For example, a line with values `3  -1  0  3` is a tree with germline root getting a mutation at gene 3 and then the resulting cell getting simultaneous mutations at genes 1 and 4.

* `example.facts`: This tab-separated file contains the learned parameter *Beta*. A possible example of this file is shown below. The file contains multiple columns. The first three columns are *'Ancestors'*, *'Support'*, and *'Support (Discont'd)'*. Column *'Ancestors'* indicate the *ancestory sets* found based on the inputs. Every ancestry set starts with 0 (the root) and the consequent mutations/genes follow by order of occurrence and sparated by semi-colon (`;`). The second column *'Support'* is the support or evidence for the corresponding ancestry set. The third column *'Support (Discont'd)'* is the support or evidence for the corresponding ancestry set without having any consequent mutations (i.e. when the last mutation of the ancstry set is a leaf). The following columns each indicate one of the input genes. The values below these columns show the evidence for the corresponding ancestry set resulting in/being followed by a mutation in those genes. For more details on how these evidences are computed see the definition of parameter *Beta* in the manuscript provided at the **Reference** section.

```
Ancestors  Support   Support (Discont'd)  2     1     3     4
0;3;1;     50        1.3                  20.8  0     0     27.9
0;3;4;     12.4      3.5                  3.5   5.4   0     0
0;3;       59        0                    1.4   48.7  0     8.9
0;         122       0                    50.1  5.7   59    5.2
...
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*Note: In this example, the spaces are adjusted for convenience. In real outputs the values are separated by a single tab. The numbers are also rounded to one decimal in this example. In real outputs the values are double-floating point values.*

* `example.post`: The structure of this file is very similar to `example.facts` except that it does not contain the *Support (Discont'd)* column and the evidences are computed based on the maximum *a posteriori* estimation trees provided in `example.trees` instead of being computed using a Bayesian approach (integrating over all possible topologies). In other words, this file contains summary statistics based on the inferred trees in file `example.trees`.


## Usage

We provide both source code (`Hintra.cpp`) and executable files. The executable `Hintra-Lin` file is compiled using

`clang++ -std=c++11 -fopenmp -lm Hintra.cpp -o Hintra-Lin`

on Linux. The `Hintra-Win.exe` file is the Windows executable produced with Visual Studio version 15. If any of these executables does not work on your system, you should compile the source code on your system. On Linux, if you do not have access to *clang*, you can replace `clang++` in the above command with `g++`, assuming that you have *g++* installed. *Please note that this code uses OpenMP API and this library should be linked (e.g. with `-fopenmp`) when compiling.*

The command for using Hintra is as follows:

`./Hintra-<X> [-u] <path_to_input> <no_samples> <no_genes> <delta> <EM_iterations> <no_cores>`

* `<X>` can be either `Lin` or `Win.exe` depending on operating system.

* `-u` is an optional argument and if used, it forces the program to recompute the `.margs` and `.maxs` files. This is useful when changing the `<delta>` parameter or when the content of the read count files changes without their names being changed.

* `<path_to_input>` is the complete path to the input dataset. For example if the `example.Rcounts` and `example.Vcounts` files are in path `C:/Project/`, then `<path_to_input>` should be set to `C:/Project/example`. Please note that no extenion is used for this argument and the `.Rcounts` and `.Vcounts` files (as well as `.margs` and `.maxs` files if they are provided) should have the same name (e.g. `example`) and be in the same folder indicated by `<path_to_input>`.

* `<no_samples>` is the number of samples which is equal to the number of lines in the `.Rcounts` or `.Vcounts` files.

* `<no_genes>` is the number of genes which is equal to the number of elements in each line of the `.Rcounts` or `.Vcounts` files.

* `<delta>` is the discretization parameter (see the reference).

* `<EM_iterations>` is the maximum number of iterations of the Expectation Maximization algorithm.

* `<no_cores>` is the number of CPU cores used for parallelization.

**Example:**

`./Hintra-Win.exe -u C:/Project/example 5 4 0.1 50 4`

`./Hintra-Lin -u C:/Project/example 5 4 0.1 50 4`

**Note: the `PhylogenySet` folder should be unzipped to the same folder as the executable.**


## Reference
{to be added}

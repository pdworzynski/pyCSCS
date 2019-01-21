pyCSCS

### Installation
```
conda install -c askerdb pycscs
```
### Examples
Download the small tests
```
wget https://raw.githubusercontent.com/askerdb/pyCSCS/master/tests/data/small_GNPS_buckettable.tsv
wget https://raw.githubusercontent.com/askerdb/pyCSCS/master/tests/data/small_GNPS_edges.tsv
```

Import the library and run the wrapper function that loads GNPS formatted formatted data
```python
>>> import pyCSCS
>>> pyCSCS.cscs_from_files("small_GNPS_buckettable.tsv", "small_GNPS_edges.tsv")
          Sample1   Sample2   Sample3   Sample4   Sample5   Sample6
Sample1  0.000000  0.350203  1.000000  1.000000  1.000000  1.000000
Sample2  0.350203  0.000000  0.440548  1.000000  1.000000  1.000000
Sample3  1.000000  0.440548  0.000000  1.000000  1.000000  1.000000
Sample4  1.000000  1.000000  1.000000  0.000000  0.601044  1.000000
Sample5  1.000000  1.000000  1.000000  0.601044  0.000000  0.182195
Sample6  1.000000  1.000000  1.000000  1.000000  0.182195  0.000000
```

### Other
To Download datasets from GNPS see the illustration in the [q2-cscs user guide](https://github.com/madeleineernst/q2-cscs#2-compute-the-chemical-structural-and-compositional-dissimilarity-for-a-real-world-dataset)

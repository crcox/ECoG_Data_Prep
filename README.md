Preparing the ECoG dataset for analysis
=======================================

Overview of the data
--------------------
The ECoG dataset has data for 12 subjects. Subjects `1--10` have
coordinates associated with their electrode labels. Subjects `11` and
`12` do not, and so solutions obtained for these subjects cannot (at
this time) be visualized. So, let's stick to working with the first 10
subjects for now.

Each subject has at least one `.mat` file that contains the timeseries
data (local field potentials) collected at `2000 Hz`. Subjects `01, 05,
06, 07, 08,` and `10` have a second `.mat` file, which is "referenced"
to a particular electrode. These referenced datasets may be preferable,
because the reference electrode may have been placed to have a good read
on autonomic responses of no interest. 

NOTE: We should figure out if the other subjects also had a reference
electrode placed, and if so if we can use it to clean up those subject's
data in a similar way.

Below is a list of all the data files for the first `10` subjects.

```bash
Pt01/namingERP_Pt01.mat
Pt01/namingERP_Pt01_refD14.mat
Pt02/namingERP_Pt02.mat
Pt03/namingERP_Pt03.mat
Pt04/namingERP_PtMA_REF4.mat
Pt05/namingERP_Pt05.mat
Pt05/namingERP_Pt05_ref02.mat
Pt06/namingERP_Pt06.mat
Pt06/namingERP_Pt06_ref01.mat
Pt07/namingERPdata_Pt07.mat
Pt07/namingERPdata_Pt07_ref02.mat
Pt08/namingERPdata_Pt08.mat
Pt08/namingERPdata_Pt08_ref03.mat
Pt09/namingERPdata_Pt09.mat
Pt10/namingERPdata_Pt10.mat
Pt10/namingERPdata_Pt10_ref02.mat
```

Source data structure
---------------------
The source data contained in the `.mat` files above appear to be
sparsely populated, but contain important information. At the top level,
the `DATA` field contains a `time x electrode` array.

```matlab
% top level

    SUBJECT: ''
       NAVE: []
       DATA: [7130000x20 double]
        DIM: [1x2 struct]
       DATE: 0
    HISTORY: {}
       MISC: []
```

The `DIM` field is a `2` element structure. The first element contains
information about the time domain, and the second element contains
information about the electrode domain. `interval` encodes the
time, in seconds, between samples (i.e., time points; rows in `DATA`).
`5e-4` is half a millisecond. `1/5e-4 = 2000 Hz`. `label` contains the
electode labels, where are critical when relating to the coordinate
space.

```matlab
.DIM(1)

        type: 'interval'
        name: 'time'
    interval: 5.0000e-04
       scale: [1x7130000 double]
       label: []
```
```matlab
.DIM(2)

        type: 'nominal'
        name: 'name'
    interval: []
       scale: []
       label: [64x4 char]
```

Supplemental data
-----------------
Supplemental data that the script depends on can be found in the
`stimuli` and `coords` directories within this repository. The
information represented in those files clear upon inspection.


Running the script
==================
The main script is `setup_data.m`. You will need to update the
`DATA_DIR` and `META_DIR` paths. `DATA_DIR` should be the path to the
directory containing the subject folders `Pt01`, `Pt02`, etc. `META_DIR`
is the path to a directory that includes the supplemental data. You may
want to update the script so that it refers to the two directories in
this repository.

When developing the script, I ran the code on a single subject. You will
need to step through the code section by section to ensure that it runs
with other subjects.

There are a handful of constants that you need to set at the top of the
script:

- AverageOverSessions
- BoxCarSize
- WindowStartInMilliseconds
- WindowSizeInMilliseconds

`AverageOverSessions` controls whether you want to average repetitions
of the same item (presented in different sessions) together into a
single time series for each item.

`BoxCarSize` give you the option to average consecutive time points
together in a very simple way. Given the vector `x = [1,2,3,4,5,6]` and
`BoxCarSize = 3`, the results of `boxcarmean(x, BoxCarSize)` would be
`[2,5]`.

The final two constants allow you to define a window, relative to each
stimulus onset, of time points that you want to retain as part of your
dataset. If you set the start of the window to `0` and the size to
`inf`, all time points will be included, up to the next stimulus onset.

`WindowStartInMillisecond` defines when the window should start, given a
number of milliseconds post stimulus onset.

`WindowSizeInMilliseconds` defines the duration of the time window.

Script output
-------------
Because the decisions to average, boxcar average, or define a window
change the dimensions of the data, each configuration of these values
can be thought of as creating a new dataset. These choices will be
codified in a directory tree, such as the one displayed below. Here you
can see I chose to average over sessions (`avg/`), I boxcar averaged
with 20 tick bins, and I defined a window that began at 0ms and was
1000ms long.

```bash
avg
└── BoxCar
    └── 020
        └── WindowStart
            └── 0000
                └── WindowSize
                    └── 1000
                        ├── metadata.mat
                        └── s01.mat
```

TO DO Items
===========
See Issues page.

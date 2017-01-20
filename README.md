# LocoPositioningSystemSlam

Matlap code for the Loco Positioning System SLAM (LPS) or Local Positioning System.
The aim of the project is to develop tools for simultaneous estimation of anchor positions
and measurement movements using a set of distance measurements.

By measuring distances 

This is a challenging parameter estimation problem that often has both noise, missing data and outliers. 
For the project we also develop tools that can be of use for other TOA type
problems, e.g. minimal solvers, ransac for initialization, non-linear L2 and L1 optimization, Wiberg algorithms,
motion priors, etc.

The code is built on several papers, please cite

Batstone, K. J., Oskarsson, M., & Åström, K. (2016). Robust Time-of-Arrival Self Calibration with Missing Data and Outliers. I 2016 24th European Signal Processing Conference (EUSIPCO). (s. 2370-2374). Institute of Electrical and Electronics Engineers Inc.. DOI: 10.1109/EUSIPCO.2016.7760673

@inbook{6505724206a7482c8bec0010f5a11197,
title = "Robust Time-of-Arrival Self Calibration with Missing Data and Outliers",
author = "Batstone, {Kenneth John} and Magnus Oskarsson and Karl Åström",
year = "2016",
month = "9",
doi = "10.1109/EUSIPCO.2016.7760673",
pages = "2370--2374",
booktitle = "2016 24th European Signal Processing Conference (EUSIPCO)",
publisher = "Institute of Electrical and Electronics Engineers Inc.",
address = "United States",
}

Other relevant papers are:

@inbook{9daca234c0a54724a83ad5ba89b10704,
  title     = "Robust time-of-arrival self calibration and indoor localization using Wi-Fi round-trip time measurements",
  author    = "Kenneth Batstone and Magnus Oskarsson and Kalle {\AA}str{\"o}m",
  year      = "2016",
  month     = "7",
  doi       = "10.1109/ICCW.2016.7503759",
  pages     = "26--31",
  booktitle = "Communications Workshops (ICC), 2016 IEEE International Conference on",
  publisher = "Institute of Electrical and Electronics Engineers Inc.",
  address   = "United States",
}

Minimal problems when there are pairs or triplets of measurement units

@article{dae5bcf45c304e0197b051f4cc5fb186,
  title     = "TOA-Based Self-Calibration of Dual-Microphone Array",
  keywords  = "Time-of-arrival (TOA), dual-microphone array, self-calibration, minimal, solver",
  author    = "Zhayida Simayijiang and Simon Burgess and Yubin Kuang and Karl {\AA}str{\"o}m",
  year      = "2015",
  doi       = "10.1109/JSTSP.2015.2417117",
  volume    = "9",
  pages     = "791--801",
  journal   = "IEEE Journal on Selected Topics in Signal Processing",
  issn      = "1941-0484",
  publisher = "IEEE--Institute of Electrical and Electronics Engineers Inc.",
  number    = "5",
}

Minimal and overdetermined algorithms for the problem when say anchors lie in a plane, while measurement units move in 3D:

@article{d90f6a6a564a4a8b84f3f7615bd4957e,
  title     = "TOA sensor network self-calibration for receiver and transmitter spaces with difference in dimension",
  keywords  = "TOA, Array calibration, Minimal problem, Ad hoc microphone arrays",
  author    = "Simon Burgess and Yubin Kuang and Karl {\AA}str{\"o}m",
  year      = "2015",
  doi       = "10.1016/j.sigpro.2014.05.034",
  volume    = "107",
  pages     = "33--42",
  journal   = "Signal Processing",
  issn      = "0165-1684",
  publisher = "Elsevier",
  number    = "Online 11 June 2014",
}

Minimal problems when there is missing (or extra measurements):

@inbook{b96bb484340b4fb1b113b1096d820dbc,
  title     = "Prime Rigid Graphs and Multidimensional Scaling with Missing Data",
  author    = "Magnus Oskarsson and Karl {\AA}str{\"o}m and Anna Torstensson",
  year      = "2014",
  doi       = "10.1109/ICPR.2014.139",
  publisher = "IEEE--Institute of Electrical and Electronics Engineers Inc.",
  pages     = "750--755",
  booktitle = "Pattern Recognition (ICPR), 2014 22nd International Conference on",
}

First paper on solving the minimal 4-6 and 5-5 problems:

@inbook{e90f122e669e44c88da841fae8645aae,
  title     = "A Complete Characterization and Solution to the Microphone Position Self-Calibration Problem",
  keywords  = "Sensor Network Calibration, TOA, Localization",
  author    = "Yubin Kuang and Simon Burgess and Anna Torstensson and Karl {\AA}str{\"o}m",
  note      = "The paper is to appear in the IEEE conference proceedings {"}Acoustics, Speech and Signal Processing (ICASSP), 2013 IEEE International Conference on{"}",
  year      = "2013",
  publisher = "IEEE--Institute of Electrical and Electronics Engineers Inc.",
  pages     = "3875--3879",
  booktitle = "A",
}

Paper on the same problem in the far-field setting:

@inbook{eaa3fa5a8944475684e9ad8dd941e1bf,
  title     = "Understanding TOA and TDOA Network Calibration using Far Field Approximation as Initial Estimate",
  keywords  = "TOA, TDOA, anchor free, localization, sensor networks",
  author    = "Yubin Kuang and Erik Ask and Simon Burgess and Karl {\AA}str{\"o}m",
  year      = "2012",
  isbn      = "978-9-898425-99-7",
  pages     = "590--596",
  booktitle = "ICPRAM 2012 - Proceedings of the 1st International Conference on Pattern Recognition Applications and Methods, Volume 2",
  publisher = "SciTePress",
}


The problem has some similarites to the low rank approximation problem:

@inbook{7f903c057cc04b659922a65dc6f60e68,
  title     = "Trust No One: Low Rank Matrix Factorization Using Hierarchical RANSAC",
  author    = "Magnus Oskarsson and Kenneth Batstone and Kalle {\AA}str{\"o}m",
  year      = "2016",
  month     = "6",
  pages     = "5820--5829",
  booktitle = "2016 IEEE Conference on Computer Vision and Pattern Recognition (CVPR), proceedings of",
  publisher = "Computer Vision Foundation",
}

@inbook{a190a386ae3d4078be7b86c26f57472e,
  title     = "On the Minimal Problems of Low-Rank Matrix Factorization",
  keywords  = "Computer vision, low rank matrix factorization, minimal problems, robust methods",
  author    = "Fangyuan Jiang and Magnus Oskarsson and Karl {\AA}str{\"o}m",
  year      = "2015",
  doi       = "10.1109/CVPR.2015.7298870",
  isbn      = "978-1-4673-6963-3",
  pages     = "2549--2557",
  editor    = "Kristen Grauman and Erik Learned-Miller and Antonio Torralba and Andrew Zisserman",
  booktitle = "Computer Vision and Pattern Recognition (CVPR), 2015 IEEE Conference on",
  publisher = "IEEE--Institute of Electrical and Electronics Engineers Inc.",
}



## Code

Run main.m 

## Data

Example data for the system can be downloaded from the
LocalPositioningSystemDatabase (lpsdb)
at [http://vision.maths.lth.se/lpsdb/](http://vision.maths.lth.se/lpsdb/)

Briefly the datasets are:

1.  Simulated - several simulations (3DR, 3DS, Discrete)
2.  bitcraze - experiments made at xxx (3DR, 3DS, Cont, Echo, Single, Distinct)

Here
* 3DR - Dimension of affine hull of receivers is 3. (Anchors not in a plane)
* 2DR - Dimension of affine hull of receivers is 2. (Anchors in a plane)
* 1DR - Dimension of affine hull of receivers is 1. (Anchors on a line)
* (0DR - Dimension of affine hull of receivers is 0. (Anchors all in the same position)
* 3DS- Dimension of affine hull of sound events is 3. (Measuring unit not in a plane)
* 2DS- Dimension of affine hull of sound events is 2. (Measuring unit in a plane)
* 1DS - Dimension of affine hull of sound events is 1. (Measuring unit on a line)
* (0DS Dimension of affine hull of sound events is 0. (Measuring unit all in the same position)
* Discrete - Measurement made at several discrete points, no motion prior on the measuring unit
* Cont - Continuous measurements - the measurement unit is moving smoothly.

### Experiments in lunch room  September 18, 2014

The microphone setup is illustrated in the following figure. 

![Lunch Room](/tex/images/IMG_2283.JPG "Lunch Room")

We made 3 experiments (bassh1, bassh2, bassh3). There are several still 
images of the microphone setup and one stationary few film recording of the
experiments.


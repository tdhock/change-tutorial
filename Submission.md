---
title: "Introduction to optimal changepoint detection algorithms"
author: |
   | Toby Dylan Hocking^1^ and Rebecca Killick^2^
   |
   | 1. McGill University, Montreal, Canada
   | 2. Lancaster University, UK
institute: 
   - $^1$McGill University, Montreal, Canada
   - $^2$Lancaster University, UK
output: html_document
---

**Keywords**: time series, changepoint, segmentation, optimization, machine learning.

## Goals/Aims/Learning objectives

Following the course participants will be able to:

* recognize datasets that potentially contain changepoints.
* identify appropriate changepoint methods dependent on the type of
  change suspected.
* create labels that indicate presence or absence of changepoints, for
  supervised analysis.
* perform changepoint analyses using unsupervised and supervised methods.
* check assumptions made within a changepoint analysis.
* summarize and evaluate results of a changepoint analysis.
* compare the accuracy of different labeled changepoint detection methods
  using cross-validation.
 
## Justification

Changepoint detection is a part of time series data analysis which is
important in fields such as finance, genomics and environment. There
are many R packages that implement changepoint detection algorithms,
but each package has a different interface. We will explain the
differences between these packages, and show examples using several
real data sets. We will emphasize coding exercises, so that
participants can get familiar with using these packages. This tutorial
will be ideal for useRs with time series data that do not yet know how
to perform changepoint analysis in R.

## Brief description of Tutorial
 
Changepoint analysis is used to model time series data that has abrupt
changes in the statistical distribution. The points in time when the
statistical properties change are referred to as changepoints, and
there are many algorithms available for computing the optimal
changepoints for a given data set. This course introduces participants
to changepoint analysis (also known as time series segmentation or
structural change detection) and available optimal algorithms. It is highly interactive and uses
packages available on CRAN and GitHub.

Techniques covered in this course include: likelihood and
nonparametric methods for detecting changes in mean and
variance. We will also cover supervised changepoint detection methods
for labeled time series data sets. Each topic is explained
theoretically and there will be break outs where participants will use
the techniques on real data sets from a variety of application areas
including finance, genomics and the environment.
 
## Detailed outline of tutorial content

Note that we will include interactive exercises using
http://rcloud.social/, each of about 2-5 minutes. These exercises
will permit students to interactively experiment and learn about
changepoint detection from the R command line.

### Rebecca, Unsupervised changepoint detection, 90 minutes

#### What is changepoint analysis and the different types of changepoints? 10 minutes

#### Methods for detecting changepoints. 60 minutes

* Single and multiple changepoint detection algorithms
* Likelihood based approaches
* Exercise: Use the cpt.mean function to see if there is evidence for a change in
mean in the Nile river data.  If you identify a change, where is it and what are the pre and post change means?
* Exercise: Use the cpt.var function to see if there is evidence for changes in
variance in the FTSE100 data from the changepoint package.  If you identify changes, where are they and what are the variances in each segment?
* Exercise: Use the cpt.meanvar function to identify regions with different C+G
content in the HC1 data within the changepoint package.
* Non parametric approaches
* Exercise: Look at the HeartRate data from the changepoint.np package. Use one of the non-parametric functions to see if there is evidence for changes in heart rate.
* Choosing the number of changes
* Exercise: Look at the FTSE100 data again and use the CROPS technique to
determine an appropriate number of changes.
* Demonstration using changepoint and changepoint.np packages

#### Checking assumptions and summarizing results of a changepoint analysis. 20 minutes

* In changepoint detection you cannot check assumptions such as Normality prior to analysis as the changes influence any diagnostics you may perform.
* Demonstation assumption checking using changepoint package and previous class exercises.
* Exercise: Check the assumptions you have made on the simulated, Nile, FTSE100 and HeartRate data using either the segment or residual check.

### Toby, Supervised changepoint detection, 70 minutes

#### What is the difference between unsupervised and supervised changepoint detection? 10 minutes

* In supervised changepoint detection, there are labels which
  indicate presence and absence of changepoints in particular data
  subsets. These labels can be used for choosing the best model and
  parameters.
* For a given set of 6 profiles (CRAN package neuroblastoma), plot
  noisy data, then superimpose labels, without showing predicted
  changepoints.
* Exercise: plot data and labels for a different set of profiles.
* Labels can be created using prior knowledge or visual inspection.
* Exercise: create a set of labels via visual inspection for one
  un-labeled segmentation problem in the neuroblastoma data set. Make
  sure that there is at least one positive and one negative label.

#### Computing the number of incorrect labels, 20 minutes

* For a given labeled segmentation problem, compute optimal Gaussian
  changepoint models for 1 to 10 segments (CRAN package
  Segmentor3IsBack).
* Compute number of incorrect labels for each model. 
* Choose the number of segments by minimizing the number of incorrect
  labels.
* Compare supervised versus unsupervised changepoint detection: many
  versus one data set, quantitative versus qualitative
  evaluation.
* Exercise: perform the same analysis on the segmentation problem that
  you labeled in the last section. Which models are optimal? (in terms
  of number of incorrect labels)

#### Supervised penalty learning, 20 minutes

* Compute the target interval of penalty values that select changepoint
  models with minimal incorrect labels. 
* Compute a feature vector for each segmentation problem, and a
  feature matrix for each labeled set of related segmentation
  problems.
* Learn an affine function f(feature vector)=log(penalty).
* Exercise: to learn the coefficients of the BIC penalty, what feature
  vector should be used?
* Un-regularized interval regression (survival package). Learns
  weights for a given set of features, but may overfit if
  non-relevant features are used.
* Elastic net regularized interval regression (anujkhare/iregnet
  package on GitHub). Simultaneously learns weights and performs
  feature selection. Avoids overfitting by setting some feature
  weights to zero.

#### Cross-validation experiments, 20 minutes

* K-fold cross-validation can be used to compare prediction accuracy
  of supervised and unsupervised changepoint detection.
* Compute test error and ROC curves for unsupervised and supervised
  penalty functions: BIC, 1 feature un-regularized, multi-feature
  un-regularized, multi-feature regularized. Which penalty function is
  most accurate?
* Exercise: perform cross-validation to compare Gaussian and Logistic
  models for log(penalty) values. Which distribution results in more
  accurate penalty functions?
 
## Pre-requisite background knowledge and packages

* Basic knowledge of R; reading in data, working with vectors and functions. 
* Basic knowledge of likelihood and hypothesis testing / model choice would be useful.
 
We will provide a script that will automatically install all packages
required for this tutorial:

* changepoint: parametric changepoint models.
* changepoint.np: non-parametric changepoint models.
* neuroblastoma: labeled data for supervised changepoint detection.
* Segmentor3IsBack: parametric changepoint models.
* survival: supervised penalty learning via un-regularized interval regression.
* iregnet: supervised penalty learning via elastic net regularized
  interval regression.

## Potential attendees

Time series data are rather common in many fields (finance, genomics, 
environment), so changepoint detection should be a rather popular
topic. We expect an audience of about 50 people, so please reserve a
large classroom or small lecture hall.

## Instructor Biographies

Toby Dylan Hocking (McGill University, Montreal, Canada;
toby.hocking@r-project.org, https://github.com/tdhock) is a
post-doctoral researcher, working on new machine learning models for
genomic data. He has implemented several R packages for changepoint
detection (neuroblastoma, bams, PeakSegDP, PeakSegJoint, coseg). He
has also implemented several R graphics packages, including
directlabels (won best student poster at useR2011) and animint
(presented in a JSM2015 invited session and a useR2016 tutorial).

Rebecca Killick (Lancaster University, UK; r.killick@lancs.ac.uk,
http://www.lancs.ac.uk/~killick) is a Lecturer in Statistics at the
University of Lancaster. Her research is in developing methodology for
the analysis of nonstationary time series to address real world
problems.  She has taught a range of courses over the last 10 years
from first year undergraduate introductory courses to PhD level
courses both theoretical and practical. She has created and
contributed to several R packages on changepoint detection including
changepoint, changepoint.np, EnvCpt and delivered a workshop on
changepoint detection at eRum2016. Her code has been adapted for
delivery in commercial software and is available as part of the NAG
libraries.


## Tutorial title: An introduction to changepoint models using R

## Instructors

Rebecca Killick (Lancaster University, UK)

Toby Dylan Hocking (McGill University, Montreal, Canada) toby.hocking@r-project.org

## Short Instructor Biography

Rebecca Killick is a Lecturer in Statistics at the University of
Lancaster. Her research is in developing methodology for the analysis
of nonstationary time series to address real world problems.  She has
taught a range of courses over the last 10 years from first year
undergraduate introductory courses to PhD level courses including
theoretical and practical courses. **TODO** discuss R packages you have developed.

Toby Dylan Hocking is a post-doctoral researcher, working on new
machine learning models for genomic data. He has implemented several R
packages for change-point detection (neuroblastoma, bams, PeakSegDP,
PeakSegJoint, coseg). He has also implemented several R graphics
packages, including directlabels (won best student poster at useR2011)
and animint (presented in a JSM2015 invited session and a useR2016
tutorial).

## Brief description of Tutorial
 
**Description of proposed workshop**

More data is being collected in the world today than ever before
resulting in larger and longer data sets. With the length of data sets
increasing, the statistical properties of the data are likely to
change over time. Traditional statistical methods often make the
assumption that the statistical properties do not change over time and
thus can lead to incorrect conclusions when changes occur.  A common
method for relaxing this assumption is to segment the data into
smaller periods within which the statistical properties do not
change. One advantage of this technique is that the traditional
(stationary) methods can be utilized on the individual segments. The
points in time when the statistical properties change are referred to
as changepoints.

This course introduces participants to the analysis of changepoint
models (also known as time series segmentation or structural
changes). The course is aimed at those with an interest in discovering
methods for models that include changepoints. It is interactive and
uses packages available on CRAN and GitHub.

Techniques covered in this course include: likelihood and
nonparametric methods for detecting changes in mean, regression and
variance. Each topic is explained theoretically and there will be
break outs where participants will use the techniques on real data
sets from a variety of application areas including finance, genomics
and environment.
 
## Goals

Following the course participants will be able to:
* recognize datasets that potentially contain changepoints
* identify appropriate changepoint methods dependent on the type of change suspected
* perform changepoint analyses using a variety of techniques in R
* summarize and evaluate results of a changepoint analysis
* check assumptions made within a changepoint analysis
 
## Detailed Outline

### Rebecca, 80 minutes

* What is changepoint analysis?
* Types of changepoints
* Methods for detecting changepoints (bulk of course)
- Nonparametric
- Likelihood
* Checking assumptions
* Summarizing results of a changepoint analysis

### Toby, Supervised change-point detection, 80 minutes

#### What is the difference between unsupervised and supervised change-point detection? 10 minutes

* In supervised change-point detection, there are labels which
  indicate presence and absence of change-points in particular data
  subsets. These labels can be used for choosing the best model and
  parameters.
* Plot noisy data, then superimpose labels, without showing
  predicted change-points.
* Demonstration via neuroblastoma data set in CRAN package
  neuroblastoma. 

#### Computing the number of incorrect labels, 25 minutes

* For a given labeled segmentation problem, compute change-point models
  with different penalty parameters.
* Compute number of incorrect labels for each model. 
* Choose the penalty parameter by minimizing the number of incorrect
  labels.
* Compare supervised versus unsupervised change-point detection: many
  versus one data set, quantitative versus qualitative
  evaluation.
* Demonstration via CRAN package cghseg.

#### Supervised penalty learning, 25 minutes

* Compute the target interval of penalty values that select change-point
  models with minimal incorrect labels. 
* Compute a feature vector for each segmentation problem, and a
  feature matrix for each labeled set of related segmentation
  problems. 
* Learn an affine function f(features)=penalty.
* Demonstration via iregnet package on GitHub.

#### Cross-validation experiments, 20 minutes

* K-fold cross-validation can be used to compare prediction accuracy
  of supervised and unsupervised change-point detection.
* Compute test error and ROC curves for BIC (unsupervised) and learned
  (supervised) penalty functions.
 
## Justification

## Background knowledge

* Basic knowledge of R; reading in data, working with vectors and functions. 
* Basic knowledge of likelihood and hypothesis testing / model choice would be useful.
 
**Required packages**

If you want to do exercises using a local copy of R, make sure to
install these packages before the tutorial:

* changepoint (and dependencies)
* changepoint.np
* TODO Toby's packages + R installation script.

TODO http://rcloud.social/ 

## Expected number of attendees

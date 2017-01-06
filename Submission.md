## Tutorial title: Change-point detection

## Instructors

Toby Dylan Hocking (McGill University, Montreal, Canada) toby.hocking@r-project.org

Rebecca Killick (Lancaster University, UK) r.killick@lancs.ac.uk

## Short Instructor Biography

Toby Dylan Hocking is a post-doctoral researcher, working on new
machine learning models for genomic data. He has implemented several R
packages for change-point detection (neuroblastoma, bams, PeakSegDP,
PeakSegJoint, coseg). He has also implemented several R graphics
packages, including directlabels (won best student poster at useR2011)
and animint (presented in a JSM2015 invited session and a useR2016
tutorial).

Rebecca Killick is a Lecturer in Statistics at the University of
Lancaster. Her research is in developing methodology for the analysis
of nonstationary time series to address real world problems.  She has
taught a range of courses over the last 10 years from first year
undergraduate introductory courses to PhD level courses both
theoretical and practical. She has created and contributed to several R packages on changepoint detection including changepoint, changepoint.np, EnvCpt and delivered a workshop on changepoint detection at eRum2016. Her code has been adapted for delivery in commercial software and is available as part of the NAG libraries.

## Brief description of Tutorial
 
TODO: should we shorten this first paragraph to make this description more "brief"?
RK: agreed, do you want to take a stab at this as i've written this paragraph so many times already :-)

More data is being collected in the world today than ever before
resulting in larger and longer data sets. With the length of data sets
increasing, the statistical properties of the data are likely to
change over time. Traditional statistical methods often make the
assumption that the statistical properties do not change over time and
thus can lead to incorrect inferences when changes are present. A common
method for relaxing this assumption is to segment the data into
smaller periods within which the statistical properties do not
change. One advantage of this technique is that the traditional
(stationary) methods can be utilized on the individual segments. The
points in time when the statistical properties change are referred to
as changepoints.

This course introduces participants to the analysis of changepoint
models (also known as time series segmentation or structural
changes). It is highly interactive and uses packages available on CRAN and
GitHub.

Techniques covered in this course include: likelihood and
nonparametric methods for detecting changes in mean and
variance. We will also cover supervised changepoint detection methods
for labeled time series data sets. Each topic is explained
theoretically and there will be break outs where participants will use
the techniques on real data sets from a variety of application areas
including finance, genomics and the environment.
 
## Goals

TODO: update the following with Toby's goals

Following the course participants will be able to:
* recognize datasets that potentially contain changepoints
* identify appropriate changepoint methods dependent on the type of change suspected
* perform changepoint analyses using a variety of techniques in R
* check assumptions made within a changepoint analysis
* summarize and evaluate results of a changepoint analysis
 
## Detailed Outline

RK: Any chance of doing a 70/90 split given i'm covering the "what is a changepoint" material too?  I did the last session in 3 hours so trying to cram the same content in less than half the time isn't going to work.

### Rebecca, 80 minutes

#### What is changepoint analysis and the different types of changepoints? 10 minutes

#### Methods for detecting changepoints. 55 minutes

* Single and multiple changepoint detection algorithms
* Likelihood based approaches
* Non parametric approaches
* Choosing the number of changes
* Demonstration using changepoint and changepoint.np packages

#### Checking assumptions and summarizing results of a changepoint analysis. 15 minutes

* In changepoint detection you cannot check assumptions such as Normality prior to analysis as the changes influence any diagnostics you may perform.
* Demonstation assumption checking using changepoint package and previous class exercises.



### Toby, Supervised changepoint detection, 80 minutes

#### What is the difference between unsupervised and supervised changepoint detection? 10 minutes

* In supervised changepoint detection, there are labels which
  indicate presence and absence of changepoints in particular data
  subsets. These labels can be used for choosing the best model and
  parameters.
* Plot noisy data, then superimpose labels, without showing
  predicted changepoints.
* Demonstration via neuroblastoma data set in CRAN package
  neuroblastoma. 

#### Computing the number of incorrect labels, 25 minutes

* For a given labeled segmentation problem, compute changepoint models
  with different penalty parameters.
* Compute number of incorrect labels for each model. 
* Choose the penalty parameter by minimizing the number of incorrect
  labels.
* Compare supervised versus unsupervised changepoint detection: many
  versus one data set, quantitative versus qualitative
  evaluation.
* Demonstration via CRAN package Segmentor3IsBack.

#### Supervised penalty learning, 25 minutes

* Compute the target interval of penalty values that select changepoint
  models with minimal incorrect labels. 
* Compute a feature vector for each segmentation problem, and a
  feature matrix for each labeled set of related segmentation
  problems. 
* Learn an affine function f(features)=penalty.
* Demonstration via iregnet package on GitHub.

#### Cross-validation experiments, 20 minutes

* K-fold cross-validation can be used to compare prediction accuracy
  of supervised and unsupervised changepoint detection.
* Compute test error and ROC curves for BIC (unsupervised) and learned
  (supervised) penalty functions.
 
## Justification

changepoint detection is a part of time series data analysis which is
important in fields such as finance, genomics and environment. There
are many R packages that implement changepoint detection algorithms,
but each package has a different interface. We will explain the
differences between these packages, and show examples using several
real data sets. We will emphasize coding exercises, so that
participants can get familiar with using these packages. This tutorial
will be ideal for useRs with time series data that do not yet know how
to perform changepoint analysis in R.

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

Time series data are rather common, so changepoint detection should be
a rather popular topic. We expect an audience of about 50 people,
so please reserve a large classroom or small lecture hall.

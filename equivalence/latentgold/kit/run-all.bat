@echo off
rem Estimate every model of the multilpa comparison kit, in order.
rem
rem Set LG to your Latent GOLD 6.1 executable. The batch flags below (/b for
rem batch, /o for the output listing) are the documented command-line form;
rem check them against the manual for your build before the first run.
rem Step-one models must run before the Step-3 models that read their output.

set LG="C:\Program Files\LatentGOLD6.1\lg61.exe"
cd /d "%~dp0"

echo Estimating c01_continuous_varying.lgs
%LG% c01_continuous_varying.lgs /b /o c01_continuous_varying.lst
echo Estimating c02_continuous_equal.lgs
%LG% c02_continuous_equal.lgs /b /o c02_continuous_equal.lst
echo Estimating c03_full_covariance.lgs
%LG% c03_full_covariance.lgs /b /o c03_full_covariance.lst
echo Estimating c04_categorical.lgs
%LG% c04_categorical.lgs /b /o c04_categorical.lst
echo Estimating c05_mixed.lgs
%LG% c05_mixed.lgs /b /o c05_mixed.lst
echo Estimating c06_missing_fiml.lgs
%LG% c06_missing_fiml.lgs /b /o c06_missing_fiml.lst
echo Estimating c07_covariates.lgs
%LG% c07_covariates.lgs /b /o c07_covariates.lst
echo Estimating c08_bivariate_residuals.lgs
%LG% c08_bivariate_residuals.lgs /b /o c08_bivariate_residuals.lst
echo Estimating c12_course_engagement.lgs
%LG% c12_course_engagement.lgs /b /o c12_course_engagement.lst
echo Estimating c09_step1.lgs
%LG% c09_step1.lgs /b /o c09_step1.lst
echo Estimating c09_distal_bch.lgs
%LG% c09_distal_bch.lgs /b /o c09_distal_bch.lst
echo Estimating c09_distal_modal.lgs
%LG% c09_distal_modal.lgs /b /o c09_distal_modal.lst
echo Estimating c09_distal_proportional.lgs
%LG% c09_distal_proportional.lgs /b /o c09_distal_proportional.lst
echo Estimating c09_covariate_ml.lgs
%LG% c09_covariate_ml.lgs /b /o c09_covariate_ml.lst
echo Estimating c10_lta.lgs
%LG% c10_lta.lgs /b /o c10_lta.lst
echo Estimating c11_lta_mixture.lgs
%LG% c11_lta_mixture.lgs /b /o c11_lta_mixture.lst

echo Done. Copy every .lst and *_posteriors.txt file back to equivalence\latentgold\returned\

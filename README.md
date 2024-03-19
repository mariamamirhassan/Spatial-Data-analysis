# Spatial-Data-analysis ( Predicting California Temperature)
## Introduction
The Environmental Statistics Project focuses on spatial interpolations to estimate temperature values in California. Spatial interpolation techniques allow us to predict temperature values at unsampled locations based on observed data from sampled locations. By analyzing temperature data from 456 locations across California, this project aims to uncover spatial patterns and trends that can inform decision-making processes in various sectors.

## Main Objectives
The primary objectives of this project are:

Apply Inverse Distance Weighting (IDW) interpolation to estimate temperature values and determine the optimal power parameter.
Utilize Trend Surface Models (TSM) to model the spatial variation of average temperature and identify the most suitable polynomial degree.
Employ kriging to estimate temperature values across the study area based on the chosen theoretical variogram model.
Data Description & Exploration
The temperature data used in this project covers the fourth quarter of the year (October, November, December) and is derived from 456 locations across California. Through data exploration, various visualizations and statistical measures were utilized to understand the distribution, skewness, spatial autocorrelation, and relationships within the dataset.

## Spatial Interpolations
Spatial interpolation techniques including Inverse Distance Weighting, Trend Surface Models, and kriging were employed to estimate temperature values across California.

## Inverse Distance Weighting part
Inverse Distance Weighting (IDW) interpolation was utilized to estimate temperature values at unsampled locations. The optimal power parameter of 3.991 was determined through grid search, indicating moderate spatial dependence in the dataset.

## Trend Surface Model
Trend Surface Models were employed to capture the spatial variation of average temperature. The 5th degree polynomial was identified as the best fitting model based on the lowest Akaike Information Criterion (AIC) value.

## Kriging
Kriging, a spatial interpolation technique, was utilized to estimate temperature values across the study area. The Mattern model was selected based on its ability to capture spatial correlation patterns observed in the temperature data.

## Conclusion
This project has provided valuable insights into the spatial patterns and variations of temperature across California through the application of spatial interpolation techniques. The results obtained can inform decision-making processes and contribute to further studies related to climate, environmental management, and spatial analysis in the region.





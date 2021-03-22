# Forecasting univariate stock indices with parametric models and neural networks: A comparison
This is the GitHub project for my master thesis. In here you find the relevant code and the writing itself. 

The aim of this thesis is to survey the accuracy and practicability of different univariate forecasting methods for different stock indices around the world. The methods used in this study are split into two groups based on their inner workings. The first group encom- passes models with a well-established statistical foundation i.e. Exponential state-space and ARIMA(p,d,q) models, while the latter focuses on artificial neural networks (ANNs) which can be considered nonparametric as they do not need statistical assumptions before estimation but rather derive a structure directly from the underlying data.

The models are configured to provide a seven-day-ahead forecast which will be compared using scale invariant error metrics. The models with statistical foundation will also be compared among themselves using maximum entropy estimators or their derivatives.

## Usage

All the relevant parts for this thesis can be accessed via the `/R/main.R` script file. The comments inside will guide you through my work. 

It is encouraged to use _R-Studio_, however all other editors work just fine. Please make sure to download the project in its entirety as the relative pathing won't work otherwise and thus the helper functions cannot be loaded.
 

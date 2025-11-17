
# ---------------------------------------------------------------
#           INTEGRATED AUTOCORRELATION TIME
# ---------------------------------------------------------------

import numpy as np

def autocorr_func(data, size, num_lags):
      """
      Compute the autocorrelation function and 
      the integrated autocorrelation time among data 

      Parameters:
            data (a Numpy array): the original dataset.
            size (int): the length of the dataset.
            num_lags (int): the number of autocorrelation lags.

      Returns:
            corr (a Numpy array): the autocorrelation function
            tau_int (a Numpy array): the integrated autocorrelation time
      """

      # Compute the mean value of the dataset
      mean = np.mean(data)

      # Compute the sample variance (with Bessel's correction 
      # (N - 1) in denominator -> ddof = 1)
      sigma_squared = np.var(data, ddof=1)

      # Create the autocorrelation function and integrated autocorrelation time 
      corr = np.zeros(num_lags)
      tau_int = np.zeros(num_lags)

      for k in range(num_lags):
            for i in range(size - k):
                  corr[k] += (data[i] - mean)*(data[i + k] - mean)
      
            corr[k] /= (float(size - k) * sigma_squared)

      for i in range(num_lags):
            for j in range (i):
                  tau_int[i] += corr[j]

      return corr, tau_int
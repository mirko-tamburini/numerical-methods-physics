
# ----------------------------------------------------
#           AUTOCORRELATION FUNCTION
#                    AND
#        INTEGRATED AUTOCORRELATION TIME
# ----------------------------------------------------

import numpy as np

def autocorr_func(data, lags):
      """
      Compute the autocorrelation function of a dataset 

      Parameters:
            data (a Numpy array): the original dataset.
            lags (int): the number of autocorrelation lags.

      Returns:
            corr (a Numpy array): the autocorrelation function
      """

      data = np.asarray(data)
      n = len(data)

      mu = np.mean(data)
      data_centered = data - mu

      sigma2 = np.var(data, ddof=1)
      corr = np.zeros(lags)

      for lag in range(lags):
            if lag == 0:
                  corr[lag] = 1.0
                  continue

            c = np.dot(data_centered[:n-lag], data_centered[lag:])
            corr[lag] = c / ((n - lag) * sigma2)

      return corr

def int_autocorr_time(corr):
      """
      Compute the integrated autocorrelation time of a dataset 

      Parameters:
            corr (a Numpy array): the autocorrelation function

      Returns:
            t_int (a Numpy array): the integrated autocorrelation time
      """
      t_int = np.zeros_like(corr)
      t_int[0] = 0.5
      t_int[1:] = 0.5 + np.cumsum(corr[1:])
      return t_int

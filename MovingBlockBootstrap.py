
# ------------------------------------------------------------------
#     MOVING BLOCK BOOTSTRAP (MBB) = technique of binned resampling
# ------------------------------------------------------------------

import numpy as np

def moving_block_bootstrap(data, func, num_bin, num_resamples = 100):
      """
      Compute the standard error using the Moving Block Bootstrap (MBB).
    
      Parameters:
            data (array-like): The original dataset.
            func (function): Function to calculate an observable.
            num_bin (int): Block size. (Typically few times the integrated autocorrelation time)
            num_resamples (int): Number of bootstrap resamples.
      
      Returns:
            float: Estimated standard error.
      """

      # Array of resampled observable
      resampled_obs = np.array([binned_resampling(data, func, num_bin) for _ in range(num_resamples)])

      # Compute standard deviation of resampled observable 
      return np.std(resampled_obs, ddof=0)


def binned_resampling(data, func, num_bin):
      """
      Generate one binned resample using the Moving Block Bootstrap.
      
      Parameters:
            data (array-like): The original dataset.
            func (function): Function to calculate an observable.
            num_bin (int): Block size.
      
      Returns:
            float: Mean of the resampled dataset.
      """

      N = len(data)
      rng = np.random.default_rng()  # Random generator

      if num_bin > N:
            raise ValueError("Block size (num_bin) cannot be larger than the dataset size.")
      
      # Preallocate space for the resampled data
      resampled_data = np.empty(N)

      # Generate random starting indices for the blocks
      block_starts = rng.integers(0, N - num_bin + 1, size = N // num_bin)

      # Copy blocks efficiently in the resampled array
      for i, start in enumerate(block_starts):
            resampled_data[i * num_bin : (i + 1) * num_bin] = data[start : start + num_bin]
      
      # Handle remainder efficiently (if any)
      remainder = N % num_bin
      if remainder:
            start = rng.integers(0, N - remainder + 1)
            resampled_data[-remainder:] = data[start : start + remainder]

      return np.func(resampled_data)


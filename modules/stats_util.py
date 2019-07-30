import numpy as np

def is_outlier(points, thresh=3.5):
	""" 
	Returns a boolean array with True if points are outliers and False
	otherwise.     

	Parameters:
	-----------
		points : An numobservations by numdimensions array of observations
		thresh : The modified z-score to use as a threshold. Observations with
			a modified z-score (based on the median absolute deviation) greater
			than this value will be classified as outliers.

	Returns:
	--------
		mask : A numobservations-length boolean array.

	References:
	----------
		Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
		Handle Outliers", The ASQC Basic References in Quality Control:
		Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
	"""

	if len(points.shape) == 1:         
		points = points[:,None]
	#print('points:')
	#print(points)

	median = np.median(points, axis=0)     
	#print('median:')
	#print(median)
	diff = np.sum((points - median)**2, axis=-1)     
	#print('diff:')
	#print(diff)
	diff = np.sqrt(diff)
	#print('diff:')
	#print(diff)

	med_abs_deviation = np.median(diff)
	if len(points) == 1 and med_abs_deviation == 0.0:
		#print([False])
		return np.array([False])

	modified_z_score = 0.6745 * diff / med_abs_deviation

	return modified_z_score > thresh

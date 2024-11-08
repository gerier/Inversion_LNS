import numpy as np
import numpy.linalg as npl
import sys

def get_norm(signal, order, choice_norm):
	# from signal extract the norm
	
	# if not choice_norm['per_rec']:
	# 	if choice_norm['ord'] == np.inf:
	# 		signal_allrec = [signal[r][t] for r in range(len(signal)) for t in range(len(signal[r]))] 
	# 		norm_signal = npl.norm(signal_allrec,  ord=np.inf)
	# 	else:
	# 		print("ERROR: cannot compute the norm 1 of sum_r^receivers d_t^r")
	# 		sys.exit()
	# else:
	# 	norm_signal = [ npl.norm(signal[r], ord=choice_norm['ord']) for r in range(len(signal))] 
    
	if choice_norm ==1:
		signal_allrec = [signal[t][r] for t in range(len(signal)) for r in range(len(signal[t]))] 
		norm_signal = npl.norm(signal_allrec,  ord=np.inf)
	if choice_norm ==2:
		norm_signal = [ npl.norm([signal[t][r] for t in range(len(signal))], ord=order) for r in range(len(signal[0]))] 
	if choice_norm ==3:
		norm_signal = [ [ abs(signal[t][r]) for r in range(len(signal[t]))] for t in range(len(signal))] 
	return norm_signal
    
    

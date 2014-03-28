__kernel void normalisation( __global double *iA, __global double *test, __global int *current, __global int *nbr, __global int *nbrItems){
	int i = get_global_id(0);
	int j = get_global_id(1);
	int p = get_global_size(0);
	int r = get_global_size(1);
	int Q = *nbrItems;
	if(i != *current) {
		double lValue = iA[i * p + *current];
		double lPrec = iA[(i-1) * p + j];
		iA[i * p + j] -= lValue * lPrec;
	}
	
	
}
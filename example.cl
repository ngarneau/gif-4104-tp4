__kernel void findBiggest( __global double *iA, __global int *index, __global int *col, __global int *numCol, __global float *maxLocal){
    int x = get_global_id(0);
    int local_x = get_local_id(0);
    if(x%*numCol == *col){
	    if(*maxLocal < iA[x]){
	    	*index=x/(*numCol);
	    	*maxLocal = iA[x];
		}
	}
}
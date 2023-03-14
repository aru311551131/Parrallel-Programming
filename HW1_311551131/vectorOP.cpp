#include "PPintrin.h"

// implementation of absSerial(), but it is vectorized using PP intrinsics
void absVector(float *values, float *output, int N)
{
  __pp_vec_float x;
  __pp_vec_float result;
  __pp_vec_float zero = _pp_vset_float(0.f);
  __pp_mask maskAll, maskIsNegative, maskIsNotNegative;

  //  Note: Take a careful look at this loop indexing.  This example
  //  code is not guaranteed to work when (N % VECTOR_WIDTH) != 0.
  //  Why is that the case?
  for (int i = 0; i < N; i += VECTOR_WIDTH)
  {

    // All ones
    maskAll = _pp_init_ones();

    // All zeros
    maskIsNegative = _pp_init_ones(0);

    // Load vector of values from contiguous memory addresses
    _pp_vload_float(x, values + i, maskAll); // x = values[i];

    // Set mask according to predicate
    _pp_vlt_float(maskIsNegative, x, zero, maskAll); // if (x < 0) {

    // Execute instruction using mask ("if" clause)
    _pp_vsub_float(result, zero, x, maskIsNegative); //   output[i] = -x;

    // Inverse maskIsNegative to generate "else" mask
    maskIsNotNegative = _pp_mask_not(maskIsNegative); // } else {

    // Execute instruction ("else" clause)
    _pp_vload_float(result, values + i, maskIsNotNegative); //   output[i] = x; }

    // Write results back to memory
    _pp_vstore_float(output + i, result, maskAll);
  }
}

void clampedExpVector(float *values, int *exponents, float *output, int N)
{
  //
  // PP STUDENTS TODO: Implement your vectorized version of
  // clampedExpSerial() here.
  //
  // Your solution should work for any value of
  // N and VECTOR_WIDTH, not just when VECTOR_WIDTH divides N
  //
  __pp_vec_float val;
  __pp_vec_int exp;
  __pp_vec_float result;
  __pp_vec_float nine;
  __pp_vec_float zero = _pp_vset_float(0.f);
  __pp_vec_float oneee = _pp_vset_float(1.f);
  __pp_vec_int zeroo = _pp_vset_int(0);
  __pp_vec_int one = _pp_vset_int(1);
  __pp_mask maskAll, maskAll0, masktmp, masktmp2, masktmp3;

  maskAll = _pp_init_ones();
  maskAll0 = _pp_init_ones(0);
  _pp_vset_float(nine, 9.999999, maskAll);
  _pp_vset_float(result, 1.0f, maskAll);
    
  for (int i = 0; i < N; i += VECTOR_WIDTH)
  {
    _pp_vset_float(result, 1.0f, maskAll);
    int max_exp = 0;
    for (int j = 0; j < VECTOR_WIDTH && j + i < N; j++) {
      max_exp = max(max_exp, exponents[j+i]);
    }
    masktmp = _pp_init_ones();
    if( N - i >= VECTOR_WIDTH) {
      
      masktmp2 = _pp_mask_and(masktmp, masktmp);
      _pp_vload_int(exp, exponents + i, masktmp); 
      _pp_vload_float(val, values + i, masktmp); 
    }
      
    else {
      __pp_vec_float nval = _pp_vset_float(0.f);
      __pp_vec_int nexp = _pp_vset_int(1);
      for (int j = 0; j < N - i; j++){
	      nval.value[j] = values[i+j];
	      nexp.value[j] = exponents[i+j];
	    }
      _pp_vmove_int(exp, nexp, masktmp); 
      _pp_vmove_float(val, nval, masktmp); 

    }
      

    
    // check exp = 0
    _pp_veq_int(masktmp3, exp, zeroo, masktmp);


    while( _pp_cntbits(masktmp) ) {
      _pp_vmult_float(result, result, val, masktmp);
      // exp - 1 && change mask
      _pp_vsub_int(exp, exp, one, masktmp);
      _pp_vgt_int(masktmp, exp, zeroo, masktmp);
    }
    // set exp = 0
    _pp_vset_float(result, 1.0f, masktmp3);

    //store
    _pp_vstore_float(output + i, result, masktmp2);

    // > 9.9999 ??
    _pp_vgt_float(masktmp, result, nine, masktmp2);
    _pp_vstore_float(output + i, nine, masktmp);
    _pp_vstore_float(output + i, oneee, masktmp3);

    

  }
}

// returns the sum of all elements in values
// You can assume N is a multiple of VECTOR_WIDTH
// You can assume VECTOR_WIDTH is a power of 2
float arraySumVector(float *values, int N)
{

  //
  // PP STUDENTS TODO: Implement your vectorized version of arraySumSerial here
  //
  __pp_vec_float val;
  __pp_vec_float exp;
  __pp_vec_float result = _pp_vset_float(0.f);
  __pp_vec_float zero = _pp_vset_float(0.f);
  __pp_vec_float one = _pp_vset_float(1.0);
  __pp_mask maskAll, maskAll0, masktmp;
  maskAll = _pp_init_ones();
  for (int i = 0; i < N; i += VECTOR_WIDTH)
  {

    _pp_vload_float(val, values + i, maskAll); 
    _pp_vadd_float(result, result, val, maskAll);
  }
  _pp_hadd_float(result, result);
  int a = VECTOR_WIDTH/2;
  while (a > 1) {
    _pp_interleave_float(result, result);
    _pp_hadd_float(result, result);
    a = a/2;
  }
  return result.value[0];
}
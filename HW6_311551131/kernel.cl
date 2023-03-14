__kernel void convolution(const __global float *img, __global float *out, const int width, 
                          const int height, __global float *filter, const int filterWidth)
{
   int index = get_global_id(0);
   int col = index % width;
   int row = index / width;
   int half_filter = filterWidth / 2;
   float sum = 0;

   for(int k = -half_filter; k <= half_filter; k++) 
   {
        if(row + k < 0 || row + k >= height)
            continue;
        for(int l = -half_filter; l <= half_filter; l++)
        {
            
            if (col + l >= 0 && col + l < width)
            {
                sum += img[(row + k) * width + col + l] *
                       filter[(k + half_filter) * filterWidth + l + half_filter];
            }
        }
   }
   
   out[row * width + col] = sum;
}

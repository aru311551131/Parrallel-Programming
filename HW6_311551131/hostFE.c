#include <stdio.h>
#include <stdlib.h>
#include "hostFE.h"
#include "helper.h"

void hostFE(int filterWidth, float *filter, int imageHeight, int imageWidth,
            float *inputImage, float *outputImage, cl_device_id *device,
            cl_context *context, cl_program *program)
{
    cl_int status;
    int filterSize = filterWidth * filterWidth;

    // initialize device
    cl_command_queue myqueue = clCreateCommandQueue(*context, *device, 0,
                                                    &status);
    // create buffer
    cl_mem img_input = clCreateBuffer(*context, CL_MEM_USE_HOST_PTR, 
                                    imageHeight * imageWidth * sizeof(float),
                                    inputImage, &status);

    cl_mem img_output = clCreateBuffer(*context, CL_MEM_WRITE_ONLY, 
                                       imageHeight * imageWidth * sizeof(float),
                                       NULL, &status);

    cl_mem filter_mem = clCreateBuffer(*context, CL_MEM_USE_HOST_PTR, 
                                       filterWidth * filterWidth * sizeof(float),
                                       filter, &status);
    
    // Select Kernel
    cl_kernel kernel = clCreateKernel(*program, "convolution", &status);
    CHECK(status, "clCreateKernel");

    // Set Arguments, Enqueue Kernel
    clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&img_input);
    clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&img_output);
    clSetKernelArg(kernel, 2, sizeof(cl_int), &imageWidth);
    clSetKernelArg(kernel, 3, sizeof(cl_int), &imageHeight);
    clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&filter_mem);
    clSetKernelArg(kernel, 5, sizeof(cl_int), &filterWidth);

    size_t global_size = imageHeight * imageWidth;
    size_t local_size = 64;
    status = clEnqueueNDRangeKernel(myqueue, kernel, 1, 0, &global_size,
                                    &local_size, 0, NULL, NULL);
    CHECK(status, "clEnqueueReadBuffer");

    //clFinish(myqueue);
    clEnqueueReadBuffer(myqueue, img_output, CL_TRUE, 0, 
                        imageHeight * imageWidth * sizeof(float), outputImage,
                        0, NULL, NULL);
    CHECK(status, "clEnqueueReadBuffer");
}
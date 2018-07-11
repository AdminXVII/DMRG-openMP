#include <stdlib.h>
#include "dmrg_vbatch.h"

#ifdef USE_MAGMA
#include "cuda.h"
#include "cuda_runtime.h"
#endif

void *dmrg_malloc(const size_t alloc_size )
{
void *ptr = NULL;
#ifdef USE_MAGMA
   {
   const unsigned int flags = cudaMemAttachGlobal;
   cudaError_t ierr = cudaMallocManaged(  &ptr, alloc_size, flags );
   int isok = (ierr == cudaSuccess );
   if (!isok) {
      fprintf(stderr,"dmrg_malloc: CUDA ERROR %s\n", cudaGetErrorString(ierr) );


      if (ierr == cudaErrorMemoryAllocation) {
         fprintf(stderr,"dmrg_malloc:cudaErrorMemoryAllocation, alloc_size=%ld\n",alloc_size);
         }
      else if (ierr == cudaErrorNotSupported) {
         fprintf(stderr,"dmrg_malloc:cudaErrorNotSupported, alloc_size=%ld\n",alloc_size);
         }
      else if (ierr == cudaErrorInvalidValue) {
         fprintf(stderr,"dmrg_malloc:cudaErrorInvalidValue, alloc_size=%ld\n",alloc_size);
         };
     };

   assert( ierr == cudaSuccess );
   }
#else
   ptr = (void *) malloc( alloc_size );
#endif
  return (ptr);
}

void dmrg_free( void *ptr )
{
#ifdef USE_MAGMA
  {
  

  cudaError_t ierr = cudaFree( ptr );
  int isok = (ierr == cudaSuccess);
  if (!isok) {
   fprintf(stderr,"dmrg_free: CUDA ERROR %s\n", cudaGetErrorString(ierr) );

   if (ierr == cudaErrorInvalidDevicePointer) {
    fprintf(stderr,"dmrg_free: cudaErrorInvalidDevicePointer \n");
    }
   else if (ierr == cudaErrorIllegalAddress) {
    fprintf(stderr,"dmrg_free: cudaErrorIllegalAddress \n");
    }
   else if (ierr == cudaErrorInitializationError) {
    fprintf(stderr,"dmrg_free: cudaErrorInitializationError \n");
    };
  };
    
  assert( ierr == cudaSuccess );
   

  };
#else
  free( ptr );
#endif
}

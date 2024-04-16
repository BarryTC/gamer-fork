#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUDA_CheckError.h"
#endif

#ifdef GRAVITY



/********************************************************
1. Point-mass external potential
   --> It can be regarded as a template for implementing
       other external potential

2. This file is shared by both CPU and GPU

   GPU_Poisson/CUPOT_ExtPot_Global_Friction_Evolution.cu -> CPU_Poisson/CPU_ExtPot_Global_Friction_Evolution.cpp

3. Three steps are required to implement external potential

   I.   Set auxiliary arrays
        --> SetExtPotAuxArray_Global_Friction_Evolution()

   II.  Specify external potential
        --> ExtPot_Global_Friction_Evolution()

   III. Set initialization functions
        --> SetGPUExtPot_Global_Friction_Evolution()
            SetCPUExtPot_Global_Friction_Evolution()
            Init_ExtPot_Global_Friction_Evolution()

4. The external potential major routine, ExtPot_Global_Friction_Evolution(),
   must be thread-safe and not use any global variable

5. Reference: https://github.com/gamer-project/gamer/wiki/Gravity#external-accelerationpotential
********************************************************/



// =================================
// I. Set auxiliary arrays
// =================================

#ifndef __CUDACC__

extern double Hosthalo_M_vir;
extern double Hosthalo_R_scale;
extern double Hosthalo_fc;


//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtPotAuxArray_Global_Friction_Evolution
// Description :  Set the auxiliary arrays ExtPot_AuxArray_Flt/Int[] used by ExtPot_Global_Friction_Evolution()
//
// Note        :  1. Invoked by Init_ExtPot_Global_Friction_Evolution()
//                2. AuxArray_Flt/Int[] have the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
void SetExtPotAuxArray_Global_Friction_Evolution( double AuxArray_Flt[], int AuxArray_Int[] )
{

// parameter unit setting

   AuxArray_Flt[0] = 0.5*amr->BoxSize[0];                         // x coordinate of the external potential center
   AuxArray_Flt[1] = 0.5*amr->BoxSize[1];                         // y ...
   AuxArray_Flt[2] = 0.5*amr->BoxSize[2];                         // z ...
   AuxArray_Flt[3] = NEWTON_G;                                    // gravitational_constant*Msun
   AuxArray_Flt[4] = Hosthalo_M_vir;                              // host halo virial radius
   AuxArray_Flt[5] = Hosthalo_R_scale;                            // host halo scale radius
   AuxArray_Flt[6] = Hosthalo_fc;                                 // host halo fc function

} // FUNCTION : SetExtPotAuxArray_Global_Friction_Evolution
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external potential
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtPot_Global_Friction_Evolution
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary arrays UserArray_Flt/Int[] are set by SetExtPotAuxArray_Global_Friction_Evolution(), where
//                      UserArray_Flt[0] = x coordinate of the external potential center
//                      UserArray_Flt[1] = y ...
//                      UserArray_Flt[2] = z ..
//                      UserArray_Flt[3] = gravitational_constant*point_source_mass
//                3. Currently it does not support the soften length
//                4. GenePtr has the size of EXT_POT_NGENE_MAX defined in Macro.h (default = 6)
//
// Parameter   :  x/y/z             : Target spatial coordinates
//                Time              : Target physical time
//                UserArray_Flt/Int : User-provided floating-point/integer auxiliary arrays
//                Usage             : Different usages of external potential when computing total potential on level Lv
//                                    --> EXT_POT_USAGE_ADD     : add external potential on Lv
//                                        EXT_POT_USAGE_SUB     : subtract external potential for preparing self-gravity potential on Lv-1
//                                        EXT_POT_USAGE_SUB_TINT: like SUB but for temporal interpolation
//                                    --> This parameter is useless in most cases
//                PotTable          : 3D potential table used by EXT_POT_TABLE
//                GenePtr           : Array of pointers for general potential tables
//
// Return      :  External potential at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real ExtPot_Global_Friction_Evolution( const double x, const double y, const double z, const double Time,
                              const double UserArray_Flt[], const int UserArray_Int[],
                              const ExtPotUsage_t Usage, const real PotTable[], void **GenePtr )
{
   const real   Hosthalo_M_vir   = (real)UserArray_Flt[4];

   if ( Hosthalo_M_vir > 0.0 )
   {
      const double Cen[3]           = { UserArray_Flt[0], UserArray_Flt[1], UserArray_Flt[2] };
      const real   GMsun            = (real)UserArray_Flt[3];
      const real   Hosthalo_R_scale = (real)UserArray_Flt[5];
      const real   Hosthalo_fc      = (real)UserArray_Flt[6];
      const real   dx               = (real)(x - Cen[0]);
      const real   dy               = (real)(y - Cen[1]);
      const real   dz               = (real)(z - Cen[2]);

      // potential components
      float Hosthalo_pot = (real)(-(Hosthalo_M_vir/(SQRT(dx*dx + dy*dy + dz*dz)*Hosthalo_fc))*log(1.0+SQRT(dx*dx + dy*dy + dz*dz)/Hosthalo_R_scale));

      return GMsun*Hosthalo_pot;
   }
   else
   {
      return 0.0;
   }

} // FUNCTION : ExtPot_Global_Friction_Evolution



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtPot_t ExtPot_Ptr = ExtPot_Global_Friction_Evolution;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtPot_Global_Friction_Evolution
// Description :  Return the function pointers of the CPU/GPU external potential routines
//
// Note        :  1. Invoked by Init_ExtPot_Global_Friction_Evolution()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      SetExtPot_Global_Friction_Evolution( ExtPot_t &CPUExtPot_Ptr, ExtPot_t &GPUExtPot_Ptr )
//
// Parameter   :  CPU/GPUExtPot_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtPot_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtPot_Global_Friction_Evolution( ExtPot_t &GPUExtPot_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtPot_Ptr, ExtPot_Ptr, sizeof(ExtPot_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtPot_Global_Friction_Evolution( ExtPot_t &CPUExtPot_Ptr )
{
   CPUExtPot_Ptr = ExtPot_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtPotAuxArray_Global_Friction_Evolution( double [], int [] );
void SetCPUExtPot_Global_Friction_Evolution( ExtPot_t & );
#ifdef GPU
void SetGPUExtPot_Global_Friction_Evolution( ExtPot_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtPot_Global_Friction_Evolution
// Description :  Initialize external potential
//
// Note        :  1. Set auxiliary arrays by invoking SetExtPotAuxArray_*()
//                   --> They will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU external potential major routines by invoking SetCPU/GPUExtPot_*()
//                3. Invoked by Init_ExtAccPot()
//                   --> Enable it by linking to the function pointer "Init_ExtPot_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Init_ExtPot_Global_Friction_Evolution()
{

   SetExtPotAuxArray_Global_Friction_Evolution( ExtPot_AuxArray_Flt, ExtPot_AuxArray_Int );
   SetCPUExtPot_Global_Friction_Evolution( CPUExtPot_Ptr );
#  ifdef GPU
   SetGPUExtPot_Global_Friction_Evolution( GPUExtPot_Ptr );
#  endif

} // FUNCTION : Init_ExtPot_Global_Friction_Evolution

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
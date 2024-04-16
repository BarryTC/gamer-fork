#include "GAMER.h"
#include "TestProb.h"

// problem-specific global variables and function prototypes
// =======================================================================================
static int      Idx_ParTag    = Idx_Undefined;           // (0: unbound; > 0: bound)
static int      Idx_ParEnergy = Idx_Undefined;           // individual particle specific total energy
static int      Idx_ParAngMom = Idx_Undefined;           // individual particle specific angular momentum
static int      Idx_ParEnergyBac = Idx_Undefined;        // individual particle specific total energy (+background)
static int      Idx_ParID     = Idx_Undefined;           // individual particle unique ID
static double   M_gas_disk;                              // total mass of gas disk
static double   r_disk_scale;                            // scale radius of gas disk
static double   r_disk_cut;                              // hard truncation radius of gas disk
static double   h_disk_scale;                            // scale height of gas disk
static double   h_disk_cut;                              // hard truncation height of gas disk
static double   T_disk;                                  // initial temperature of gas disk
static double   T_ambient;                               // ambient temperature of gas background
static double   Dens_central;                            // gas disk central density
static double   Dens_ambient;                            // gas background density
static char     RotVel_Table_File[MAX_STRING];           // filename of the DM+stars rotation curve table
static double  *Table_R=NULL;                            // DM+stars rotation curve table [radius]
static double  *Table_V=NULL;                            // DM+stars rotation curve table [velocity]
static int      RotVel_Table_NBin;                       // number of bins in RotVel_Table_Data[]
static bool     RotVel_Table_UnitConv;                   // convert from [kpc, km/sec] to code units; OFF(0)/ON(1)
static bool     OPT_Auto_AMR;                            // whether to enable "OPT__FLAG_USER"; OFF(0)/ON(1)

#ifdef PARTICLE
void AddNewParticleAttribute_Disk_Gas_Stream();
#endif
# ifdef GRAVITY
bool Flag_Disk_Gas_Stream( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );
#endif
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


// errors
#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   if ( amr->Par->Interp != PAR_INTERP_TSC )
      Aux_Error( ERROR_INFO, "please set PAR_INTERP = 3 (TSC) !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      for (int f=0; f<6; f++)
      if ( OPT__BC_FLU[f] == BC_FLU_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for fluid is not recommended for this test !!\n" );

#     ifdef GRAVITY
      if ( OPT__BC_POT == BC_POT_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for gravity is not recommended for this test !!\n" );
#     endif

      //if ( !OPT__RECORD_USER )
      //   Aux_Error( ERROR_INFO, "please set OPT__RECORD_USER = 1 (ON) !!\n" );
   } // if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//                3. Must call EoS_Init() before calling any other EoS routine
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",           &VARIABLE,                   DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "M_gas_disk",                &M_gas_disk,                 2.3e10,        0.0,              NoMax_double      );
   ReadPara->Add( "r_disk_scale",              &r_disk_scale,               6.0e3,         0.0,              NoMax_double      );
   ReadPara->Add( "r_disk_cut",                &r_disk_cut,                 1.5e4,         0.0,              NoMax_double      );
   ReadPara->Add( "h_disk_scale",              &h_disk_scale,               4.0e2,         0.0,              NoMax_double      );
   ReadPara->Add( "h_disk_cut",                &h_disk_cut,                 3.0e3,         0.0,              NoMax_double      );
   ReadPara->Add( "T_disk",                    &T_disk,                     1.0e4,         1.0,              NoMax_double      );
   ReadPara->Add( "T_ambient",                 &T_ambient,                  1.0e6,         1.0,              NoMax_double      );
   ReadPara->Add( "RotVel_Table_File",         RotVel_Table_File,           NoDef_str,     Useless_str,      Useless_str       );
   ReadPara->Add( "RotVel_Table_UnitConv",     &RotVel_Table_UnitConv,      true,          Useless_bool,     Useless_bool      );
   ReadPara->Add( "OPT_Auto_AMR",              &OPT_Auto_AMR,               true,          Useless_bool,     Useless_bool      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check and reset the runtime parameters


// (2) set the problem-specific derived parameters
Dens_central = M_gas_disk / (4.0*M_PI*SQR(r_disk_scale)*h_disk_scale*(1.0 - exp(-r_disk_cut/r_disk_scale)*(1.0 + r_disk_cut/r_disk_scale))*(1.0-exp(-h_disk_cut/h_disk_scale)));
Dens_ambient = Dens_central*exp(-r_disk_cut/r_disk_scale)*exp(-h_disk_cut/h_disk_scale)*1.0e-4;


// (3) load DM+stars rotation curve for gas disk construction
   if ( OPT__INIT != INIT_BY_RESTART )
   {
      const bool RowMajor_No  = false;       // load data into the column-major order
      const bool AllocMem_Yes = true;        // allocate memory for Merger_Prof1/2
      const int  NCol         = 2;           // total number of columns to load
      const int  Col[NCol]    = {0, 1};      // target columns: (radius, velocity)
      double    *RotVel_Table_Data = NULL;   // DM+stars rotation curve table [radius/velocity]

      RotVel_Table_NBin = Aux_LoadTable( RotVel_Table_Data, RotVel_Table_File, NCol, Col, RowMajor_No, AllocMem_Yes );
//    convert input table of rotation curve to code units
      Table_R = RotVel_Table_Data + 0*RotVel_Table_NBin;
      Table_V = RotVel_Table_Data + 1*RotVel_Table_NBin;
      if ( RotVel_Table_UnitConv )
      {
         for (int b=0; b<RotVel_Table_NBin; b++)
         {
            Table_R[b] *= Const_kpc / UNIT_L;
            Table_V[b] *= (Const_km/Const_s) / UNIT_V;
         }
      }

   } // if ( OPT__INIT != INIT_BY_RESTART )


// (4) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 10.0*Const_Gyr/UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }


// (5) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                                            = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  total mass of gas disk                                     = %13.7e\n", M_gas_disk );
      Aux_Message( stdout, "  scale radius of gas disk                                   = %13.7e\n", r_disk_scale );
      Aux_Message( stdout, "  hard truncation radius of gas disk                         = %13.7e\n", r_disk_cut );
      Aux_Message( stdout, "  scale height of gas disk                                   = %13.7e\n", h_disk_scale );
      Aux_Message( stdout, "  hard truncation height of gas disk                         = %13.7e\n", h_disk_cut );
      Aux_Message( stdout, "  initial temperature of gas disk                            = %13.7e\n", T_disk );
      Aux_Message( stdout, "  ambient temperature of gas background                      = %13.7e\n", T_ambient );
      Aux_Message( stdout, "  DM+stars rotation curve Table_File                         = %s\n",     RotVel_Table_File );
      Aux_Message( stdout, "  convert from [kpc, km/sec] to code units; OFF(0)/ON(1)     = %d\n",     RotVel_Table_UnitConv );
      Aux_Message( stdout, "  set radius-dependent AMR on-the-fly; OFF(0)/ON(1)          = %d\n",     OPT_Auto_AMR );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter




//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  *Modeled after Omry's RAMSES source code (..\patch\init\agama\condinit.f90)
//                1. This function sets up the initial gas disk on grid
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{
      const double dcell    = amr->dh[lv];                     // cell size at level "lv"
      const double x_dis    = x - 0.5*amr->BoxSize[0];         // x coordinate wrt simulation box center
      const double y_dis    = y - 0.5*amr->BoxSize[1];         // y coordinate wrt simulation box center
      const double z_dis    = z - 0.5*amr->BoxSize[2];         // z coordinate wrt simulation box center
      double       r_dis    = sqrt(SQR(x_dis) + SQR(y_dis));   // r coordinate wrt simulation box center

//    set up the gas disk on grid
      if ( abs(z_dis) < h_disk_cut  &&  r_dis < r_disk_cut ){
//       part of gas disk component
         const real V_cir   = Mis_InterpolateFromTable( RotVel_Table_NBin, Table_R, Table_V, r_dis );
         const real weight  = (FMIN(r_dis+dcell/2.0,r_disk_cut) - (r_dis-dcell/2.0))/dcell;      // smoothing factor for cells overlap with truncation radii
         if (weight != 1.0)  r_dis = r_dis + (weight-1.0)*dcell/2.0;                             // modify "r_dis" the a cell extends across truncation radii

         fluid[DENS]        = FMAX(weight*(Dens_central*exp(-r_dis/r_disk_scale)*exp(-abs(z_dis)/h_disk_scale)), Dens_ambient);
         const real Pres    = EoS_DensTemp2Pres_CPUPtr( fluid[DENS], T_disk, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
         const real Eint    = EoS_DensPres2Eint_CPUPtr( fluid[DENS], Pres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
         const real Cs      = sqrt(GAMMA*Pres/fluid[DENS]);                                      // gas sound speed of this cell
         const real Vel_rot = sqrt(FMAX(SQR(V_cir) - SQR(Cs)*(r_dis/r_disk_scale), 0.0));        // gas total (rotational) velocity of this cell
         fluid[MOMX]        = fluid[DENS]*(-y_dis/r_dis)*Vel_rot;
         fluid[MOMY]        = fluid[DENS]*(x_dis/r_dis)*Vel_rot;
         fluid[MOMZ]        = 0.0;
         fluid[ENGY]        = Hydro_ConEint2Etot( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], Eint, 0.0 );
      }
      else{
//       part of ambient gas background
         fluid[DENS]        = Dens_ambient;
         fluid[MOMX]        = (real) 0.0;
         fluid[MOMY]        = (real) 0.0;
         fluid[MOMZ]        = (real) 0.0;
         const real Pres    = EoS_DensTemp2Pres_CPUPtr( fluid[DENS], T_ambient, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
         const real Eint    = EoS_DensPres2Eint_CPUPtr( fluid[DENS], Pres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
         fluid[ENGY]        = Hydro_ConEint2Etot( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], Eint, 0.0 );
      }
} // FUNCTION : SetGridIC
#  endif // #if ( MODEL == HYDRO )




//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewParticleAttribute_Disk_Gas_Stream
// Description :  Add test-problem-specific particle attributes
//
// Note        :  1. "Idx_ParTag" has "0" for unbound particles and positive integers for bound particles
//                2. "Idx_ParEnergy" stores the total energy of individual particles
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewParticleAttribute_Disk_Gas_Stream()
{
   if ( Idx_ParTag          == Idx_Undefined )
      Idx_ParTag            = AddParticleAttribute( "ParTag" );
   if ( Idx_ParEnergy       == Idx_Undefined )
      Idx_ParEnergy         = AddParticleAttribute( "ParEnergy" );
   if ( Idx_ParAngMom       == Idx_Undefined )
      Idx_ParAngMom         = AddParticleAttribute( "ParAngMom" );
   if ( Idx_ParEnergyBac    == Idx_Undefined )
      Idx_ParEnergyBac      = AddParticleAttribute( "ParEnergyBac" );
   if ( Idx_ParID           == Idx_Undefined )
      Idx_ParID             = AddParticleAttribute( "ParID" );
} // FUNCTION : AddNewParticleAttribute_Disk_Gas_Stream




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Disk_Gas_Stream
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Disk_Gas_Stream()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();
   Init_Function_User_Ptr  = SetGridIC;
   //if (OPT_Auto_AMR)    Flag_User_Ptr = Flag_Disk_Gas_Stream;

#  ifdef PARTICLE
   Par_Init_Attribute_User_Ptr = AddNewParticleAttribute_Disk_Gas_Stream;
#  endif
#  endif // #if ( MODEL == HYDRO )

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Disk_Gas_Stream

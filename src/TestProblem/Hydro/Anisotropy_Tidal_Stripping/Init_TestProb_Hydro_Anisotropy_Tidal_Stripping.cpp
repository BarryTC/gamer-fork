#include "GAMER.h"
#include "TestProb.h"

// problem-specific global variables and function prototypes
// =======================================================================================
static int      Idx_ParTag    = Idx_Undefined;           // (0: unbound; > 0: bound)
static int      Idx_ParEnergy = Idx_Undefined;           // individual particle specific total energy
static int      Idx_ParAngMom = Idx_Undefined;           // individual particle specific angular momentum
static int      Idx_ParEnergyBac = Idx_Undefined;        // individual particle specific total energy (+background)
static int      Idx_ParID     = Idx_Undefined;           // individual particle unique ID
const  int      NParBoundMin = 10;                       // the minimum number of bound particles considered when determining CM position and velocity
       int      NParBound_Sum;                           // total number of bound particles across all ranks
       int      Ter_NPar;                                // applicable when "OPT_Auto_Termination == 1"
static bool     Sorting_Step  = true;                    // whether to execute particle sorting in this step
static bool     OPT_beta_compute;                        // whether to compute velocity anisotropy on-the-fly
static bool     OPT_Auto_Termination;                    // whether to automatically terminate simulation when "N_{par, bound} < 1"
static bool     OPT_Auto_AMR;                            // whether to enable "OPT__FLAG_USER"
static double   Particle_Mass;                           // simulated (unmodified) individual particle mass
static double   Subhalo_Ini_R_vir;                       // subhalo initial virial radius
static double   Subhalo_Ini_vel_x;                       // subhalo initial x-component velocity
static double   Subhalo_Ini_vel_y;                       // subhalo initial y-component velocity
static double   Subhalo_Ini_vel_z;                       // subhalo initial z-component velocity
static double   Subhalo_Ini_vel_cir;                     // circular velocity at the subhalo virial radius "= SQRT(G*M_vir/Subhalo_Ini_R_vir)"
       double   Subhalo_CM[3] = { __DBL_MAX__, __DBL_MAX__, __DBL_MAX__ };  // subhalo center of mass coordinates
static double   Subhalo_vel[3];                          // subhalo center of mass velocity
static double   Hosthalo_R_vir;                          // host halo virial radius [code units]
static double   Hosthalo_c;                              // host halo concentration parameter
static double   Hosthalo_CM[3];                          // host halo center of mass coordinates
       double   Hosthalo_M_vir;                          // host halo virial mass [code units]
       double   Hosthalo_R_scale;
       double   Hosthalo_fc;

#ifdef PARTICLE
void Init_ExtPot_Anisotropy_Tidal_Stripping();
void AddNewParticleAttribute_Anisotropy_Tidal_Stripping();
void Compute_Den_Max_On_Grid( double DenMaxOnGridPos[] );
void Output_UserWorkBeforeOutput_Anisotropy_Tidal_Stripping();
void Aux_Record_Anisotropy_Tidal_Stripping();
#endif
# ifdef GRAVITY
bool Flag_Anisotropy_Tidal_Stripping( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );
void Recompute_Potential_Anisotropy_Tidal_Stripping( const int OptExtPot_Input );
void ParEnergy_Calculation_Anisotropy_Tidal_Stripping( const bool OPT_Self_Total_Pot );
void ParEnergy_Sorting_Anisotropy_Tidal_Stripping( const bool FirstTime_MassRest, long& NParBound_Sum_5percent,
                        double& CM_Pos_Relative_Error, double& CM_Vel_Relative_Error );
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

      if ( !OPT__RECORD_USER )
         Aux_Error( ERROR_INFO, "please set OPT__RECORD_USER = 1 (ON) !!\n" );
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
   ReadPara->Add( "Particle_Mass",             &Particle_Mass,              1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Subhalo_Ini_R_vir",         &Subhalo_Ini_R_vir,          1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Subhalo_Ini_vel_x",         &Subhalo_Ini_vel_x,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Subhalo_Ini_vel_y",         &Subhalo_Ini_vel_y,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Subhalo_Ini_vel_z",         &Subhalo_Ini_vel_z,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Subhalo_Ini_vel_cir",       &Subhalo_Ini_vel_cir,        0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Hosthalo_M_vir",            &Hosthalo_M_vir,             1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Hosthalo_R_vir",            &Hosthalo_R_vir,             1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Hosthalo_c",                &Hosthalo_c,                 5.0,           0.0,              NoMax_double      );
   ReadPara->Add( "OPT_beta_compute",          &OPT_beta_compute,           true,          Useless_bool,     Useless_bool      );
   ReadPara->Add( "OPT_Auto_Termination",      &OPT_Auto_Termination,       true,          Useless_bool,     Useless_bool      );
   ReadPara->Add( "Ter_NPar",                  &Ter_NPar,                   1,             0,                NoMax_int         );
   ReadPara->Add( "OPT_Auto_AMR",              &OPT_Auto_AMR,               true,          Useless_bool,     Useless_bool      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   Subhalo_vel[0] = Subhalo_Ini_vel_x;
   Subhalo_vel[1] = Subhalo_Ini_vel_y;
   Subhalo_vel[2] = Subhalo_Ini_vel_z;
   for (int i=0; i<3; i++)  Hosthalo_CM[i] = 0.5*amr->BoxSize[i];

// (1-3) check and reset the runtime parameters


// (2) set the problem-specific derived parameters
   Hosthalo_R_scale = (real)(Hosthalo_R_vir/Hosthalo_c);
   Hosthalo_fc      = (real)(log(1.0+Hosthalo_c)-1.0/(1.0+(1.0/Hosthalo_c)));


// (3) reset other general-purpose parameters
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


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                                            = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  individual subhalo particle mass                           = %13.7e\n", Particle_Mass );
      Aux_Message( stdout, "  subhalo initial virial radius                              = %13.7e\n", Subhalo_Ini_R_vir );
      Aux_Message( stdout, "  subhalo initial x velocity                                 = %13.7e\n", Subhalo_Ini_vel_x );
      Aux_Message( stdout, "  subhalo initial y velocity                                 = %13.7e\n", Subhalo_Ini_vel_y );
      Aux_Message( stdout, "  subhalo initial z velocity                                 = %13.7e\n", Subhalo_Ini_vel_z );
      Aux_Message( stdout, "  subhalo initial circular velocity at R_vir                 = %13.7e\n", Subhalo_Ini_vel_cir );
      Aux_Message( stdout, "  host halo virial mass                                      = %13.7e\n", Hosthalo_M_vir );
      Aux_Message( stdout, "  host halo virial radius                                    = %13.7e\n", Hosthalo_R_vir );
      Aux_Message( stdout, "  host halo concentration                                    = %13.7e\n", Hosthalo_c );
      Aux_Message( stdout, "  output velocity anisotropy on-the-fly; OFF(0)/ON(1)        = %d\n",     OPT_beta_compute );
      Aux_Message( stdout, "  terminate simulation when N_{par} < Ter_NPar; OFF(0)/ON(1) = %d\n",     OPT_Auto_Termination );
      Aux_Message( stdout, "  applicable when OPT_Auto_Termination == 1                  = %d\n",     Ter_NPar );
      Aux_Message( stdout, "  set radius-dependent AMR on-the-fly; OFF(0)/ON(1)          = %d\n",     OPT_Auto_AMR );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
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

      double Dens, MomX, MomY, MomZ, Eint, Etot;
//    gas density and energy cannot be zero so set them to an extremely small value
      Dens = 10.0*TINY_NUMBER;
      MomX = 0.0;
      MomY = 0.0;
      MomZ = 0.0;
      Eint = 10.0*TINY_NUMBER;
      Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );

      fluid[DENS] = Dens;
      fluid[MOMX] = MomX;
      fluid[MOMY] = MomY;
      fluid[MOMZ] = MomZ;
      fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC
#  endif // #if ( MODEL == HYDRO )




//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewParticleAttribute_Anisotropy_Tidal_Stripping
// Description :  Add test-problem-specific particle attributes
//
// Note        :  1. "Idx_ParTag" has "0" for unbound particles and positive integers for bound particles
//                2. "Idx_ParEnergy" stores the total energy of individual particles
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewParticleAttribute_Anisotropy_Tidal_Stripping()
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
} // FUNCTION : AddNewParticleAttribute_Anisotropy_Tidal_Stripping




//-------------------------------------------------------------------------------------------------------
// Function    :  Recompute_Potential_Anisotropy_Tidal_Stripping
// Description :  Recompute gravitational potential on all grids
//
// Note        :  See "src/Fluid/Flu_CorrAfterAllSync.cpp" and the options for "OptExtPot_t" in "include/Typedef.h"
//
// Parameter   :  OptExtPot_Input : External potential setting
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Recompute_Potential_Anisotropy_Tidal_Stripping( const int OptExtPot_Input )
{
# ifdef GRAVITY
   // enable or disable the external potential
   OPT__EXT_POT = OptExtPot_Input;

   // re-evaluate on-grid gravitational potential
   for (int lv=0; lv<=MAX_LEVEL; lv++)
   {
      if ( NPatchTotal[lv] == 0 )   break;

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "      recalculate potential at Lv %2d        ... ", lv );

      if ( lv > 0 )
      Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, NULL_INT, DATA_GENERAL, _DENS, _NONE, Rho_ParaBuf, USELB_YES );

      Gra_AdvanceDt( lv, Time[lv], NULL_REAL, NULL_REAL, NULL_INT, amr->PotSg[lv], true, false, false, false, false );

      if ( lv > 0 )
      Buf_GetBufferData( lv, NULL_INT, NULL_INT, amr->PotSg[lv], POT_FOR_POISSON, _POTE, _NONE, Pot_ParaBuf, USELB_YES );

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

   }
#  endif // #ifdef GRAVITY
} // FUNCTION : Recompute_Potential_Anisotropy_Tidal_Stripping




//-------------------------------------------------------------------------------------------------------
// Function    :  Compute_Den_Max_On_Grid
// Description :  Compute the grid coordinates of on-grid maximum density
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Compute_Den_Max_On_Grid( double DenMaxOnGridPos[] )
{
   const int    CountMPI          = 4;
   double dens, max_dens_loc=-__DBL_MAX__, max_dens_pos_loc[3];
   double send[CountMPI], (*recv)[CountMPI] = new double [MPI_NRank][CountMPI];

   const bool   IntPhase_No       = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const real   MinTemp_No        = -1.0;
   const real   MinEntr_No        = -1.0;
   const bool   DE_Consistency_No = false;
#  ifdef PARTICLE
   const bool   TimingSendPar_No  = false;
   const bool   PredictParPos_No  = false;
   const bool   JustCountNPar_No  = false;
#  ifdef LOAD_BALANCE
   const bool   SibBufPatch       = true;
   const bool   FaSibBufPatch     = true;
#  else
   const bool   SibBufPatch       = NULL_BOOL;
   const bool   FaSibBufPatch     = NULL_BOOL;
#  endif
#  endif // #ifdef PARTICLE

   for (int lv=0; lv<NLEVEL; lv++)
   {
//    initialize the particle density array (rho_ext) and collect particles to the target level
#     ifdef PARTICLE
      Par_CollectParticle2OneLevel( lv, _PAR_MASS|_PAR_POSX|_PAR_POSY|_PAR_POSZ|_PAR_TYPE, PredictParPos_No, NULL_REAL,
                                    SibBufPatch, FaSibBufPatch, JustCountNPar_No, TimingSendPar_No );

      Prepare_PatchData_InitParticleDensityArray( lv, Time[lv] );
#     endif

//    get the total density on grids
      real (*TotalDens)[PS1][PS1][PS1] = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
      int   *PID0List                  = new int  [ amr->NPatchComma[lv][1]/8 ];

      for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

      Prepare_PatchData( lv, Time[lv], TotalDens[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List, _PAR_DENS, _NONE,
                         OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                         MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

      delete [] PID0List;

//    free memory for collecting particles from other ranks and levels, and free density arrays with ghost zones (rho_ext)
#     ifdef PARTICLE
      Par_CollectParticle2OneLevel( lv, _PAR_MASS|_PAR_POSX|_PAR_POSY|_PAR_POSZ|_PAR_TYPE, PredictParPos_No, NULL_REAL,
                                    SibBufPatch, FaSibBufPatch, JustCountNPar_No, TimingSendPar_No );

      Prepare_PatchData_FreeParticleDensityArray( lv );
#     endif

//    determine the coordinates of on-grid maximum density in each process
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       skip non-leaf patches
         if ( amr->patch[0][lv][PID]->son != -1 )  continue;

         for (int k=0; k<PS1; k++)  {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*amr->dh[lv];
         for (int j=0; j<PS1; j++)  {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*amr->dh[lv];
         for (int i=0; i<PS1; i++)  {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*amr->dh[lv];

            dens = TotalDens[PID][k][j][i];

            if ( dens > max_dens_loc )
            {
               max_dens_loc        = dens;
               max_dens_pos_loc[0] = x;
               max_dens_pos_loc[1] = y;
               max_dens_pos_loc[2] = z;
            }
         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

      delete [] TotalDens;
   } // for (int lv=0; lv<NLEVEL; lv++)


// gather data to the root rank
   send[0] = max_dens_loc;
   send[1] = max_dens_pos_loc[0];
   send[2] = max_dens_pos_loc[1];
   send[3] = max_dens_pos_loc[2];

   MPI_Gather( send, CountMPI, MPI_DOUBLE, recv[0], CountMPI, MPI_DOUBLE, 0, MPI_COMM_WORLD );

// determine the maximum density coordinates
   double max_dens      = -__DBL_MAX__;
   int    max_dens_rank = -1;

   if ( MPI_Rank == 0 )
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         if ( recv[r][0] > max_dens )
         {
            max_dens      = recv[r][0];
            max_dens_rank = r;
         }
      }

      if ( max_dens_rank < 0  ||  max_dens_rank >= MPI_NRank )
         Aux_Error( ERROR_INFO, "incorrect max_dens_rank (%d) !!\n", max_dens_rank );
   } // if ( MPI_Rank == 0 )

// return the coordinates of on-grid maximum density
   if ( MPI_Rank == 0 )
      for (int d=0; d<3; d++)    DenMaxOnGridPos[d] = recv[max_dens_rank][1+d];
} // FUNCTION : Compute_Den_Max_On_Grid




//-------------------------------------------------------------------------------------------------------
// Function    :  ParEnergy_Calculation_Anisotropy_Tidal_Stripping
// Description :  Compute individual particle total energy
//
// Note        :  (1) self_potential + kinetic or (2) boosted + kinetic
//                [Part II] loop thorough all particles to (1) assign unbound particles negative mass and (2) determine the total energy of individual particles
//                The code block is modeled after src/Particle/Par_Aux_GetConservedQuantity.cpp
//
// Parameter   :  None
//
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void ParEnergy_Calculation_Anisotropy_Tidal_Stripping( const bool OPT_Self_Total_Pot )
{
   const real *Pos[3]               = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   const real *Vel[3]               = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };
   const real *PMass                = amr->Par->Mass;
   const real *PType                = amr->Par->Type;
   const long  ParN_ThisRank        = amr->Par->NPar_AcPlusInac;

// 1. kinematic energy
   for (long p=0; p<ParN_ThisRank; p++)
   {
//    skip inactive (with negative mass), and tracer particles
      if ( PMass[p] > (real)0.0  &&  PType[p] != PTYPE_TRACER )
      {
//       compute and record the particle kinetic energy
         if (OPT_Self_Total_Pot)
         {
            amr->Par->Attribute[Idx_ParEnergy][p] = 0.5*( SQR(Vel[0][p]-Subhalo_vel[0]) + SQR(Vel[1][p]-Subhalo_vel[1]) + SQR(Vel[2][p]-Subhalo_vel[2]) );
//          mark unbound particles as (practically) massless
            if ( amr->Par->Attribute[Idx_ParTag][p] == (real)0.0 )    amr->Par->Mass[p] = 10.0*TINY_NUMBER;
         }
         else  amr->Par->Attribute[Idx_ParEnergyBac][p] = 0.5*( SQR(Vel[0][p]) + SQR(Vel[1][p]) + SQR(Vel[2][p]) );
      }
   }

// 2. potential energy
   if (OPT_Self_Total_Pot)  Recompute_Potential_Anisotropy_Tidal_Stripping( EXT_POT_NONE ); // temporarily exclude the background static gravitational potential and recompute potential on all grids

#  ifdef GRAVITY
   const ParInterp_t IntScheme  = amr->Par->Interp;
   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinDens_No        = -1.0;
   const real MinPres_No        = -1.0;
   const real MinTemp_No        = -1.0;
   const real MinEntr_No        = -1.0;

   const int  PotGhost          = amr->Par->GhostSize;
   const int  PotSize           = PS1 + 2*PotGhost;

   double Ep_Coeff = 1.0;
   double PrepPotTime, dh, _dh;
   int    PotSg;

// check
// Par->ImproveAcc only supports particle interpolation schemes with ParGhost == 1 (CIC & TSC)
   if ( amr->Par->ImproveAcc  &&  PotGhost != GRA_GHOST_SIZE-1 )
      Aux_Error( ERROR_INFO, "PotGhost (%d) != GRA_GHOST_SIZE-1 (%d) for amr->Par->ImproveAcc !!\n",
               PotGhost, GRA_GHOST_SIZE-1 );

// Par->ImproveAcc must work with STORE_POT_GHOST
#  ifndef STORE_POT_GHOST
   if ( amr->Par->ImproveAcc )
      Aux_Error( ERROR_INFO, "amr->Par->ImproveAcc must work with STORE_POT_GHOST !!\n" );
#  endif

// loop over particles in all leaf patches to compute potential energy
   else
   for (int lv=0; lv<NLEVEL; lv++)
   {
      dh          = amr->dh[lv];
      _dh         = 1.0/dh;
      PrepPotTime = Time[lv];

//    determine PotSg for STORE_POT_GHOST
#     ifdef STORE_POT_GHOST
      if ( amr->Par->ImproveAcc )
      {
         if      (  Mis_CompareRealValue( PrepPotTime, amr->PotSgTime[lv][   amr->PotSg[lv] ], NULL, false )  )
            PotSg =   amr->PotSg[lv];

         else if (  Mis_CompareRealValue( PrepPotTime, amr->PotSgTime[lv][ 1-amr->PotSg[lv] ], NULL, false )  )
            PotSg = 1-amr->PotSg[lv];

         else
            Aux_Error( ERROR_INFO, "Cannot determine PotSg (lv %d, PrepTime %20.14e, SgTime[0] %20.14e, SgTime[1] %20.14e !!\n",
                     lv, PrepPotTime, amr->PotSgTime[lv][0], amr->PotSgTime[lv][1] );
      }
#     endif


//    OpenMP parallel region
#     pragma omp parallel
      {

//    per-thread variables
      bool   GotYou;
      long   ParID;

      real *Pot = new real [ 8*CUBE(PotSize) ];    // 8: number of patches per patch group
      typedef real (*vla)[PotSize][PotSize][PotSize];
      vla Pot3D = ( vla )Pot;

#     pragma omp for schedule( runtime )
      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
      {
//       2-1. find the patch groups with particles
//       --> use patch group as the calculation unit since Prepare_PatchData() only work with patch group
//       --> some patches may not have particles ...
         GotYou = false;

         for (int PID=PID0; PID<PID0+8; PID++)
         {
            if ( amr->patch[0][lv][PID]->NPar > 0 )
            {
               GotYou = true;
               break;
            }
         }


//       nothing to do if there are no particles in the target patch group
         if ( !GotYou )    continue;


//       2-2. prepare the potential data for the patch group with particles (need NSIDE_26 for ParGhost>0)
#        ifdef STORE_POT_GHOST
         if ( amr->Par->ImproveAcc )
         {
            const int didx = 1;  // assuming GRA_GHOST_SIZE - Pot_Size = 1

            for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
            {
               if ( amr->patch[0][lv][PID]->NPar == 0 )  continue;   // skip patches with no particles

               for (int k=didx; k<GRA_NXT-didx; k++)
               for (int j=didx; j<GRA_NXT-didx; j++)
               for (int i=didx; i<GRA_NXT-didx; i++)
                  Pot3D[P][k-didx][j-didx][i-didx] = amr->patch[PotSg][lv][PID]->pot_ext[k][j][i];
            }
         }

         else
#        endif
            Prepare_PatchData( lv, PrepPotTime, Pot, NULL, PotGhost, 1, &PID0, _POTE, _NONE,
                              OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCH, (PotGhost==0)?NSIDE_00:NSIDE_26, IntPhase_No,
                              OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );


//       2-3. calculate potential energy
         for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
         {
            if ( amr->patch[0][lv][PID]->NPar == 0 )  continue;   // skip patches with no particles

            switch ( IntScheme )
            {
//             2-3-3. TSC
               case ( PAR_INTERP_TSC ):
               {
                  int    idxLCR[3][3];    // array index of the left/central/right cells (idxLCR[0/1/2][d])
                  double dr       [3];    // distance to the left edge of the central cell
                  double Frac  [3][3];    // weighting of the left/central/right cells (Frac[0/1/2][d])

                  for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
                  {
                     ParID = amr->patch[0][lv][PID]->ParList[p];

//                   skip tracer and inactive (with negative mass) particles
                     if ( PType[ParID] == PTYPE_TRACER  ||  PMass[ParID] < (real)0.0 )
                        continue;

                     for (int d=0; d<3; d++)
                     {
//                      calculate the array index of the left, central, and right cells
                        dr       [d]  = ( Pos[d][ParID] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh + PotGhost;
                        idxLCR[1][d]  = int( dr[d] );
                        idxLCR[0][d]  = idxLCR[1][d] - 1;
                        idxLCR[2][d]  = idxLCR[1][d] + 1;
                        dr       [d] -= (double)idxLCR[1][d];

//                      prevent from round-off errors (especially for NGP and TSC)
                        if ( idxLCR[0][d] < 0 )
                        {
#                          ifdef DEBUG_PARTICLE
                           if (  ! Mis_CompareRealValue( Pos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                           Aux_Error( ERROR_INFO, "index outside the pot array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                                    d, Pos[d][ParID], amr->patch[0][lv][PID]->EdgeL[d], idxLCR[0][d], idxLCR[2][d] );
#                          endif

                           idxLCR[0][d] = 0;
                           idxLCR[1][d] = 1;
                           idxLCR[2][d] = 2;
                        }

                        else if ( idxLCR[2][d] >= PotSize )
                        {
#                          ifdef DEBUG_PARTICLE
                           if (  ! Mis_CompareRealValue( Pos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                           Aux_Error( ERROR_INFO, "index outside the pot array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                                    d, Pos[d][ParID], amr->patch[0][lv][PID]->EdgeR[d], idxLCR[0][d], idxLCR[2][d] );
#                          endif

                           idxLCR[0][d] = PotSize - 3;
                           idxLCR[1][d] = PotSize - 2;
                           idxLCR[2][d] = PotSize - 1;
                        }

//                      get the weighting of the nearby 27 cells
                        Frac[0][d] = 0.5*SQR( 1.0 - dr[d] );
                        Frac[1][d] = 0.5*( 1.0 + 2.0*dr[d] - 2.0*SQR(dr[d]) );
                        Frac[2][d] = 0.5*SQR( dr[d] );
                     } // for (int d=0; d<3; d++)

//                   get potential energy
                     for (int k=0; k<3; k++)
                     for (int j=0; j<3; j++)
                     for (int i=0; i<3; i++)
                     if (OPT_Self_Total_Pot)  amr->Par->Attribute[Idx_ParEnergy][ParID] += Ep_Coeff*Pot3D[P][ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
                                                                  *Frac[i][0]*Frac[j][1]*Frac[k][2];
                     else  amr->Par->Attribute[Idx_ParEnergyBac][ParID] += Ep_Coeff*Pot3D[P][ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
                                                                  *Frac[i][0]*Frac[j][1]*Frac[k][2];

                  } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               } // PAR_INTERP_TSC
               break;

               default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
            } // switch ( IntScheme )
         } // for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)

//    2-4. free memory
      delete [] Pot;

      } // end of OpenMP parallel region
   } // for (int lv=0; lv<NLEVEL; lv++)
#  endif // #ifdef GRAVITY
//Aux_Message( stdout, "      Potential Computation Completes at MPI_Rank = %10d\n", MPI_Rank);
} // FUNCTION : ParEnergy_Calculation_Anisotropy_Tidal_Stripping




//-------------------------------------------------------------------------------------------------------
// Function    :  ParEnergy_Sorting_Anisotropy_Tidal_Stripping
// Description :  Compute individual particle total energy
//
// Note        :  (1) self_potential + kinetic or (2) boosted + kinetic
//                [Part II] loop thorough all particles to (1) assign unbound particles negative mass and (2) determine the total energy of individual particles
//                The code block is modeled after src/Particle/Par_Aux_GetConservedQuantity.cpp
//
// Parameter   :
//
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void ParEnergy_Sorting_Anisotropy_Tidal_Stripping( const int NParBoundMin, long& NParBound_Sum_5percent,
                        double& CM_Pos_Relative_Error, double& CM_Vel_Relative_Error )
{
// [Part III] relabel bound/unbound particles, and perform non-parallel sorting based on the particle total specific energy
   const real *Pos[3]           = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   const real *Vel[3]           = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };
   const real *PMass            = amr->Par->Mass;
   const real *PType            = amr->Par->Type;
   const long  ParN_ThisRank    = amr->Par->NPar_AcPlusInac;
   const long  ParNAc_ThisRank  = amr->Par->NPar_Active;        // number of active (non-negative mass and non-tracer) particles
   long  NParBound_ThisRank     = 0;                            // number of bound particles in this rank
   long  NParUnbound_ThisRank   = 0;                            // number of unbound particles in this rank
   long *ParBoundp_ThisRank     = new long [ ParNAc_ThisRank ]; // record bound particle indices "p" in this rank
   long *ParUnboundp_ThisRank   = new long [ ParNAc_ThisRank ]; // record unbound particle indices "p" in this rank
   long  Par_Idx_Now;
   int   NParUnbound_Sum        = 0;                            // total number of unbound particles across all ranks
   double CM_Pos[3]             = { 0.0, 0.0, 0.0 };
   double CM_Vel[3]             = { 0.0, 0.0, 0.0 };

   for (long p=0; p<ParN_ThisRank; p++)
   {
//    skip inactive (with negative mass), and tracer particles
      if ( PMass[p] > (real)0.0  &&  PType[p] != PTYPE_TRACER )
      {
//       unbound particles: set "Idx_ParTag = 0.0"
         if ( amr->Par->Attribute[Idx_ParEnergy][p] >= (real)0.0 )
         {
            ParUnboundp_ThisRank[NParUnbound_ThisRank] = p;
            amr->Par->Attribute[Idx_ParTag][p]         = 0.0;

            ++NParUnbound_ThisRank;
         }
//       bound particles: collect entry indices and restore particle mass (if previously marked as massless)
         else
         {
            ParBoundp_ThisRank[NParBound_ThisRank] = p;
            if ( PMass[p] > (real)0.0  &&  PMass[p] < 0.5*Particle_Mass )    amr->Par->Mass[p] = Particle_Mass;

            ++NParBound_ThisRank;
         }
      }
   }
   //Aux_Message( stdout, "      NParUnbound_ThisRank = %10d, NParBound_ThisRank = %10d\n", NParUnbound_ThisRank, NParBound_ThisRank);


// gather bound particle attributes to the root rank
   double *ParBoundAtt_ThisRank[7];                                  // (in order) PosX, PosY, PosZ, VelX, VelY, VelZ, ParEnergy
   long *ParBoundpIdx_ThisRank   = new long [ NParBound_ThisRank ];  // index "p"
   long *ParBoundEngIdx_ThisRank = new long [ NParBound_ThisRank ];  // energy-sorted indices
   for (int v=0; v<7; v++)  ParBoundAtt_ThisRank[v] = new double [ NParBound_ThisRank ];

   for (long p=0; p<NParBound_ThisRank; p++)
   {
      Par_Idx_Now                = ParBoundp_ThisRank[p];

      ParBoundAtt_ThisRank[0][p] = Pos[0][Par_Idx_Now];  // PosX
      ParBoundAtt_ThisRank[1][p] = Pos[1][Par_Idx_Now];  // PosY
      ParBoundAtt_ThisRank[2][p] = Pos[2][Par_Idx_Now];  // PosZ
      ParBoundAtt_ThisRank[3][p] = Vel[0][Par_Idx_Now];  // VelX
      ParBoundAtt_ThisRank[4][p] = Vel[1][Par_Idx_Now];  // VelY
      ParBoundAtt_ThisRank[5][p] = Vel[2][Par_Idx_Now];  // VelZ
      ParBoundAtt_ThisRank[6][p] = amr->Par->Attribute[Idx_ParEnergy][Par_Idx_Now]; // ParEnergy
      ParBoundpIdx_ThisRank[p]   = Par_Idx_Now;          // index "p"
   }

// free memory
   delete [] ParBoundp_ThisRank;

// get the number of bound particles in each rank to the rook rank
// see "src/Particle/Par_Aux_Record_ParticleCount.cpp" and "../Par_ScatterParticleData.cpp"
   long NParBound_AllRank[MPI_NRank];
   int Count[MPI_NRank], Disp[MPI_NRank];
   MPI_Barrier( MPI_COMM_WORLD );
   MPI_Gather( &NParBound_ThisRank, 1, MPI_LONG, NParBound_AllRank, 1, MPI_LONG, 0, MPI_COMM_WORLD );

// collect the bound particle attributes to the root rank
// see "src/SelfGravity/Poi_GetAverageDensity.cpp"
   if ( MPI_Rank == 0 )
   {
      Disp[0] = 0;
      for (int r=0; r<MPI_NRank; r++)
      {
         NParBound_Sum += NParBound_AllRank[r];
         Count[r]       = NParBound_AllRank[r];
      }
      for (int r=1; r<MPI_NRank; r++)  Disp[r] = Disp[r-1] + Count[r-1];
   } // if ( MPI_Rank == 0 )

   NParBound_Sum_5percent = MAX(NParBoundMin, round(0.05*NParBound_Sum));       // 5% of the total bound particle counts
   double *ParBoundAtt_AllRank[7];                                              // (in order) PosX, PosY, PosZ, VelX, VelY, VelZ, ParEnergy
   for (int v=0; v<7; v++)  ParBoundAtt_AllRank[v] = new double [ NParBound_Sum ];
   long *ParBoundpIdx_AllRank = new long [ NParBound_Sum ];                     // ParTag
   MPI_Bcast( &NParBound_Sum, 1, MPI_INT, 0, MPI_COMM_WORLD );
   MPI_Bcast( &NParBound_Sum_5percent, 1, MPI_LONG, 0, MPI_COMM_WORLD );
   //Aux_Message( stdout, "      NParBound_Sum = %10d\n", NParBound_Sum);

// gather bound particle information to the root rank
   for (int v=0; v<7; v++)  MPI_Gatherv( ParBoundAtt_ThisRank[v], NParBound_ThisRank, MPI_DOUBLE, ParBoundAtt_AllRank[v], Count, Disp, MPI_DOUBLE, 0, MPI_COMM_WORLD );

// gather additional least unbound particle information if "NParBoundMin > NParBound_Sum"
   if (NParBoundMin > NParBound_Sum)
   {
      Aux_Message( stdout, "      NParBoundMin = %10d, NParBound_Sum = %10d, NParBound_Sum_5percent = %10d\n", NParBoundMin, NParBound_Sum, NParBound_Sum_5percent );
//    gather unbound particle attributes to the root rank
      double *ParUnboundAtt_ThisRank[7];                                      // (in order) PosX, PosY, PosZ, VelX, VelY, VelZ, ParEnergy
      for (int v=0; v<7; v++)  ParUnboundAtt_ThisRank[v] = new double [ NParUnbound_ThisRank ];
      for (long p=0; p<NParUnbound_ThisRank; p++)
      {
         Par_Idx_Now                  = ParUnboundp_ThisRank[p];

         ParUnboundAtt_ThisRank[0][p] = Pos[0][Par_Idx_Now];  // PosX
         ParUnboundAtt_ThisRank[1][p] = Pos[1][Par_Idx_Now];  // PosY
         ParUnboundAtt_ThisRank[2][p] = Pos[2][Par_Idx_Now];  // PosZ
         ParUnboundAtt_ThisRank[3][p] = Vel[0][Par_Idx_Now];  // VelX
         ParUnboundAtt_ThisRank[4][p] = Vel[1][Par_Idx_Now];  // VelY
         ParUnboundAtt_ThisRank[5][p] = Vel[2][Par_Idx_Now];  // VelZ
         ParUnboundAtt_ThisRank[6][p] = amr->Par->Attribute[Idx_ParEnergy][Par_Idx_Now]; // ParEnergy
      }

//    get the number of unbound particles in each rank to the rook rank
      long NParUnbound_AllRank[MPI_NRank];
      int Count_UB[MPI_NRank], Disp_UB[MPI_NRank];
      MPI_Gather( &NParUnbound_ThisRank, 1, MPI_LONG, NParUnbound_AllRank, 1, MPI_LONG, 0, MPI_COMM_WORLD );

//    collect the unbound particle attributes to the root rank
      if ( MPI_Rank == 0 )
      {
         Disp_UB[0] = 0;
         for (int r=0; r<MPI_NRank; r++)
         {
            NParUnbound_Sum += NParUnbound_AllRank[r];
            Count_UB[r]      = NParUnbound_AllRank[r];
         }
         for (int r=1; r<MPI_NRank; r++)  Disp_UB[r] = Disp_UB[r-1] + Count_UB[r-1];
      } // if ( MPI_Rank == 0 )

      double *ParUnboundAtt_AllRank[7];                   // (in order) PosX, PosY, PosZ, VelX, VelY, VelZ, ParEnergy
      for (int v=0; v<7; v++)  ParUnboundAtt_AllRank[v] = new double [ NParUnbound_Sum ];
//    gather bound particle information to the root rank
      for (int v=0; v<7; v++)  MPI_Gatherv( ParUnboundAtt_ThisRank[v], NParUnbound_ThisRank, MPI_DOUBLE, ParUnboundAtt_AllRank[v], Count_UB, Disp_UB, MPI_DOUBLE, 0, MPI_COMM_WORLD );

//    sort particle total energy entries;
      int *SortByEnergy_IdxTable_UB = new int [ NParUnbound_Sum ];
      Mis_Heapsort( NParUnbound_Sum, ParUnboundAtt_AllRank[6], SortByEnergy_IdxTable_UB );

//    compute the CM position and velocity from the 5% most bound particles
      if ( MPI_Rank == 0 )
      {
         for (int p=0; p<(NParBoundMin-NParBound_Sum); p++)
         {
            Par_Idx_Now = SortByEnergy_IdxTable_UB[p];

            CM_Pos[0]  += ParUnboundAtt_AllRank[0][Par_Idx_Now];  // PosX
            CM_Pos[1]  += ParUnboundAtt_AllRank[1][Par_Idx_Now];  // PosY
            CM_Pos[2]  += ParUnboundAtt_AllRank[2][Par_Idx_Now];  // PosZ
            CM_Vel[0]  += ParUnboundAtt_AllRank[3][Par_Idx_Now];  // VelX
            CM_Vel[1]  += ParUnboundAtt_AllRank[4][Par_Idx_Now];  // VelY
            CM_Vel[2]  += ParUnboundAtt_AllRank[5][Par_Idx_Now];  // VelZ
         }
      }

//    free memory
      delete [] SortByEnergy_IdxTable_UB;
      for (int v=0; v<7; v++)  delete [] ParUnboundAtt_ThisRank[v];
      for (int v=0; v<7; v++)  delete [] ParUnboundAtt_AllRank[v];
   } // if (NParBoundMin > NParBound_Sum)


// free memory
   delete [] ParUnboundp_ThisRank;


   if ( MPI_Rank == 0 )
   {
//    sort particle total energy entries; see "src/Miscellaneous/Mis_Heapsort.cpp" and "/tool/analysis/gamer_compare_data/GAMER_Functions/SortParticle.cpp"
      int *SortByEnergy_IdxTable = new int [ NParBound_Sum ];
      Mis_Heapsort( NParBound_Sum, ParBoundAtt_AllRank[6], SortByEnergy_IdxTable );

      if (NParBoundMin > NParBound_Sum)
      {
         for (int p=0; p<NParBound_Sum; p++)
         {
            Par_Idx_Now = SortByEnergy_IdxTable[p];

            CM_Pos[0]  += ParBoundAtt_AllRank[0][Par_Idx_Now];  // PosX
            CM_Pos[1]  += ParBoundAtt_AllRank[1][Par_Idx_Now];  // PosY
            CM_Pos[2]  += ParBoundAtt_AllRank[2][Par_Idx_Now];  // PosZ
            CM_Vel[0]  += ParBoundAtt_AllRank[3][Par_Idx_Now];  // VelX
            CM_Vel[1]  += ParBoundAtt_AllRank[4][Par_Idx_Now];  // VelY
            CM_Vel[2]  += ParBoundAtt_AllRank[5][Par_Idx_Now];  // VelZ
         }
      }
      else
      {
         for (int p=0; p<NParBound_Sum_5percent; p++)
         {
            Par_Idx_Now = SortByEnergy_IdxTable[p];

            CM_Pos[0]  += ParBoundAtt_AllRank[0][Par_Idx_Now];  // PosX
            CM_Pos[1]  += ParBoundAtt_AllRank[1][Par_Idx_Now];  // PosY
            CM_Pos[2]  += ParBoundAtt_AllRank[2][Par_Idx_Now];  // PosZ
            CM_Vel[0]  += ParBoundAtt_AllRank[3][Par_Idx_Now];  // VelX
            CM_Vel[1]  += ParBoundAtt_AllRank[4][Par_Idx_Now];  // VelY
            CM_Vel[2]  += ParBoundAtt_AllRank[5][Par_Idx_Now];  // VelZ
         }
      }

//    update the subhalo CM position and velocity
      for (int v=0; v<3; v++)
      {
         CM_Pos[v]  /= (double)NParBound_Sum_5percent;
         CM_Vel[v]  /= (double)NParBound_Sum_5percent;
      }

      CM_Pos_Relative_Error = SQRT(SQR(CM_Pos[0]-Subhalo_CM[0])+SQR(CM_Pos[1]-Subhalo_CM[1])+SQR(CM_Pos[2]-Subhalo_CM[2]))/Subhalo_Ini_R_vir;
      CM_Vel_Relative_Error = SQRT(SQR(CM_Vel[0]-Subhalo_vel[0])+SQR(CM_Vel[1]-Subhalo_vel[1])+SQR(CM_Vel[2]-Subhalo_vel[2]))/Subhalo_Ini_vel_cir;

      for (int v=0; v<3; v++)
      {
         Subhalo_CM[v]  = CM_Pos[v];
         Subhalo_vel[v] = CM_Vel[v];
      }

//    record the ParTag list
      for (int p=0; p<NParBound_Sum; p++)  ParBoundpIdx_AllRank[SortByEnergy_IdxTable[p]] = p;

//    free memory
      delete [] SortByEnergy_IdxTable;
   } // if ( MPI_Rank == 0 )

   if (NParBound_Sum > 0)
   {
//    send energy-sorted indices from the root rank to all ranks
      MPI_Scatterv( ParBoundpIdx_AllRank, Count, Disp, MPI_LONG, ParBoundEngIdx_ThisRank, NParBound_ThisRank, MPI_LONG, 0, MPI_COMM_WORLD );

//    map energy-sorted indicies back to "ParTag" in each rank
      for (long p=0; p<NParBound_ThisRank; p++)
      {
         Par_Idx_Now                                  = ParBoundpIdx_ThisRank[p];  // index "p"

         amr->Par->Attribute[Idx_ParTag][Par_Idx_Now] = ParBoundEngIdx_ThisRank[p] + 1;
      }
   }

// broadcast
   MPI_Bcast( &CM_Pos_Relative_Error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Bcast( &CM_Vel_Relative_Error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Bcast( Subhalo_CM, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Bcast( Subhalo_vel, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

// free memory
   for (int v=0; v<7; v++)  delete [] ParBoundAtt_AllRank[v];
   for (int v=0; v<7; v++)  delete [] ParBoundAtt_ThisRank[v];
   delete [] ParBoundpIdx_AllRank;
   delete [] ParBoundpIdx_ThisRank;
   delete [] ParBoundEngIdx_ThisRank;

} // FUNCTION : ParEnergy_Sorting_Anisotropy_Tidal_Stripping




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_UserWorkBeforeOutput_Anisotropy_Tidal_Stripping
// Description :  Iteratively determine the simulated particle CM and rank particles based on the
//                respective total energy based on the method of van den Bosch and Go, MNRAS (2018)
//
//
// Note        :  1. The function iteratively determine the subhalo CM by selecting a fixed percentage
//                   (5%) of the most energetically bound particles
//                2. Unbound particles are labeled as "ParTag = 0." Bound particles all have positive
//                   integer values for "ParTag" starting from "1," order ranked such that the more
//                   bound particles have smaller "ParTag" values
//                3. In each iteration, the particle self-potential (in evaluating the particle total
//                   energy) includes contributions from bound particles ONLY
//
// Procedure   :  1. Unbound particles are temporarily marked as massless and the background potential is
//                   also temporarily disabled, after the subhalo self-potential is recomputed accounting
//                   only for bound members of the subhalo
//                2. Compute the specific energy (total energy per unit mass) for individual particles,
//                   and the values are stored in the column "ParEnergy"
//                3. Unbound particles are labeled as "ParTag = 0." Bound particles, if previously marked
//                   as massless, have their mass restored. All bound particles are then ranked based on
//                   the respective specific energy (non-parallel sorting) such that "ParTag = 1 (2) ..."
//                   is assigned to the (second) most bound particle and so on so forth
//                4. Recompute the "Subhalo_CM[3]" and "Subhalo_vel[3]" from the 5% most bound particles,
//                   and then repeat again from step 1. The iteration stops when the relative changes in
//                   CM position and velocity are smaller than 1e-4
//                5. Before exiting the function, we re-enable the background potential, restore particle
//                   mass, re-compute the total gravitational potential
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Output_UserWorkBeforeOutput_Anisotropy_Tidal_Stripping()
{
   const char  filename_center[]      = "Record__Center";
   const char  filename_beta[]        = "Record__Beta";
   const real *Pos[3]                 = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   const real *Vel[3]                 = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };
   const real *PMass                  = amr->Par->Mass;
   const real *PType                  = amr->Par->Type;
   const long  ParN_ThisRank          = amr->Par->NPar_AcPlusInac;
   const long  ParNAc_ThisRank        = amr->Par->NPar_Active;   // number of active (non-negative mass and non-tracer) particles
   const int   NParCMEst              = 100;                     // the number of most bound particles considered when estimating CM position and velocity
   const int   NIterMax               = 50;
   static long NParBound_Sum_5percent = 0;                       // 5% of the total bound particle counts
   int    NIter                       = 0;
   double CM_Pos_Relative_Error       = 1.0;
   double CM_Vel_Relative_Error       = 1.0;
   bool   FirstTime_MassRest          = true;
   static bool FirstTime              = true;
   long Par_Idx_Now;

// ============================================================================================================
// [Part I] Pre-iteration setups

// adopt the on-grid maximum density position as the initial guess of the subhalo CM
   if ( Step == (real)0 )  Compute_Den_Max_On_Grid( Subhalo_CM );
// adopt the CM position and velocity inferred from the most bound particles in the previous output
   else if (NParBound_Sum_5percent > NParBoundMin)
   {
      Aux_Message( stdout, "      NParBound_Sum_5percent = %10d, Step = %10d\n", NParBound_Sum_5percent, Step);
      //Aux_Message( stdout, "      Before: Subhalo_vel[0] = %14.7e, Subhalo_vel[1] = %14.7e, Subhalo_vel[2] = %14.7e\n", Subhalo_vel[0], Subhalo_vel[1], Subhalo_vel[2]);
      double CM_Pos_Est_ThisRank[3] = { 0.0, 0.0, 0.0 };
      double CM_Pos_Est_AllRank[3]  = { 0.0, 0.0, 0.0 };
      double CM_Vel_Est_ThisRank[3] = { 0.0, 0.0, 0.0 };
      double CM_Vel_Est_AllRank[3]  = { 0.0, 0.0, 0.0 };
      long ParVel_Sum_ThisRank      = 0;
      long ParVel_Sum_AllRank       = 0;

      for (long p=0; p<ParN_ThisRank; p++)
      {
//       skip inactive (with negative mass), and tracer particles
         if ( (amr->Par->Attribute[Idx_ParTag][p] > (real)0.0) && (amr->Par->Attribute[Idx_ParTag][p] <= NParBound_Sum_5percent)
            && (PMass[p] > (real)0.0) &&  (PType[p] != PTYPE_TRACER) )
         {
               CM_Pos_Est_ThisRank[0]  += Pos[0][p];  // PosX
               CM_Pos_Est_ThisRank[1]  += Pos[1][p];  // PosY
               CM_Pos_Est_ThisRank[2]  += Pos[2][p];  // PosZ
               CM_Vel_Est_ThisRank[0]  += Vel[0][p];  // VelX
               CM_Vel_Est_ThisRank[1]  += Vel[1][p];  // VelY
               CM_Vel_Est_ThisRank[2]  += Vel[2][p];  // VelZ
               ++ParVel_Sum_ThisRank;
         }
      }
      MPI_Barrier( MPI_COMM_WORLD );
      MPI_Reduce( CM_Pos_Est_ThisRank, CM_Pos_Est_AllRank, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( CM_Vel_Est_ThisRank, CM_Vel_Est_AllRank, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &ParVel_Sum_ThisRank, &ParVel_Sum_AllRank, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );

      if ( MPI_Rank == 0 )
      {
         //Aux_Message( stdout, "      NParBound_Sum_5percent = %10d, ParVel_Sum_AllRank = %10d\n", NParBound_Sum_5percent, ParVel_Sum_AllRank);
         for (int v=0; v<3; v++)
         {
            Subhalo_CM[v]  = CM_Vel_Est_AllRank[v]/(double)ParVel_Sum_AllRank;
            Subhalo_vel[v] = CM_Vel_Est_AllRank[v]/(double)ParVel_Sum_AllRank;
         }
         Aux_Message( stdout, "      After: Subhalo_CM[0] = %14.7e,  Subhalo_CM[1] = %14.7e,  Subhalo_CM[2] = %14.7e\n", Subhalo_CM[0], Subhalo_CM[1], Subhalo_CM[2]);
         Aux_Message( stdout, "      After: Subhalo_vel[0] = %14.7e, Subhalo_vel[1] = %14.7e, Subhalo_vel[2] = %14.7e\n", Subhalo_vel[0], Subhalo_vel[1], Subhalo_vel[2]);
      }

      MPI_Bcast( Subhalo_CM, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Bcast( Subhalo_vel, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   }
// adopt the CM position and velocity inferred from the "NParCMEst" most bound particles
   else
   {
//    gather all particle attributes to the root rank
      long ParMAc_Counter       = 0;
      double *ParAllAtt_ThisRank[7];                            // (in order) PosX, PosY, PosZ, VelX, VelY, VelZ, ParEnergy
      for (int v=0; v<7; v++)  ParAllAtt_ThisRank[v] = new double [ ParNAc_ThisRank ];

      for (long p=0; p<ParN_ThisRank; p++)
      {
//       skip inactive (with negative mass), and tracer particles
         if ( (PMass[p] > (real)0.0) &&  (PType[p] != PTYPE_TRACER) )
         {
            ParAllAtt_ThisRank[0][ParMAc_Counter] = Pos[0][p];  // PosX
            ParAllAtt_ThisRank[1][ParMAc_Counter] = Pos[1][p];  // PosY
            ParAllAtt_ThisRank[2][ParMAc_Counter] = Pos[2][p];  // PosZ
            ParAllAtt_ThisRank[3][ParMAc_Counter] = Vel[0][p];  // VelX
            ParAllAtt_ThisRank[4][ParMAc_Counter] = Vel[1][p];  // VelY
            ParAllAtt_ThisRank[5][ParMAc_Counter] = Vel[2][p];  // VelZ
            ParAllAtt_ThisRank[6][ParMAc_Counter] = amr->Par->Attribute[Idx_ParEnergy][p]; // ParEnergy

            ++ParMAc_Counter;
         }
      }

      if (ParMAc_Counter != ParNAc_ThisRank)  Aux_Message( stdout, "      At MPI_Rank = %10d: ParMAc_Counter (%10d) != ParNAc_ThisRank(%10d)\n", MPI_Rank, ParMAc_Counter, ParNAc_ThisRank);

//    get the number of active particles in each rank to the rook rank
//    see "src/Particle/Par_Aux_Record_ParticleCount.cpp" and "../Par_ScatterParticleData.cpp"
      long NParAll_AllRank[MPI_NRank];
      int NParAll_Sum = 0;
      int Count_All[MPI_NRank], Disp_All[MPI_NRank];
      MPI_Gather( &ParNAc_ThisRank, 1, MPI_LONG, NParAll_AllRank, 1, MPI_LONG, 0, MPI_COMM_WORLD );

//    collect the bound particle attributes to the root rank
//    see "src/SelfGravity/Poi_GetAverageDensity.cpp"
      if ( MPI_Rank == 0 )
      {
         Disp_All[0] = 0;
         for (int r=0; r<MPI_NRank; r++)
         {
            NParAll_Sum += NParAll_AllRank[r];
            Count_All[r] = NParAll_AllRank[r];
         }
         for (int r=1; r<MPI_NRank; r++)  Disp_All[r] = Disp_All[r-1] + Count_All[r-1];
      } // if ( MPI_Rank == 0 )

      double *ParAllAtt_AllRank[7];                                              // (in order) PosX, PosY, PosZ, VelX, VelY, VelZ, ParEnergy
      for (int v=0; v<7; v++)  ParAllAtt_AllRank[v] = new double [ NParAll_Sum ];
      MPI_Bcast( &NParAll_Sum, 1, MPI_INT, 0, MPI_COMM_WORLD );

//    gather all particle information to the root rank
      for (int v=0; v<7; v++)  MPI_Gatherv( ParAllAtt_ThisRank[v], ParNAc_ThisRank, MPI_DOUBLE, ParAllAtt_AllRank[v], Count_All, Disp_All, MPI_DOUBLE, 0, MPI_COMM_WORLD );

      if ( MPI_Rank == 0 )
      {
//       sort particle total energy entries
         double CM_Pos_Est[3] = { 0.0, 0.0, 0.0 };
         double CM_Vel_Est[3] = { 0.0, 0.0, 0.0 };
         int *SortByEnergy_IdxTable_All = new int [ NParAll_Sum ];
         Mis_Heapsort( NParAll_Sum, ParAllAtt_AllRank[6], SortByEnergy_IdxTable_All );

         for (int p=0; p<NParCMEst; p++)
         {
            Par_Idx_Now     = SortByEnergy_IdxTable_All[p];

            CM_Pos_Est[0]  += ParAllAtt_AllRank[0][Par_Idx_Now];  // PosX
            CM_Pos_Est[1]  += ParAllAtt_AllRank[1][Par_Idx_Now];  // PosY
            CM_Pos_Est[2]  += ParAllAtt_AllRank[2][Par_Idx_Now];  // PosZ
            CM_Vel_Est[0]  += ParAllAtt_AllRank[3][Par_Idx_Now];  // VelX
            CM_Vel_Est[1]  += ParAllAtt_AllRank[4][Par_Idx_Now];  // VelY
            CM_Vel_Est[2]  += ParAllAtt_AllRank[5][Par_Idx_Now];  // VelZ
            //Aux_Message( stdout, "      p = %10d, Par_Idx_Now = %10d, Par_Energy = %14.7e\n", p, Par_Idx_Now, ParAllAtt_AllRank[6][Par_Idx_Now]);
         }

         for (int v=0; v<3; v++)
         {
            Subhalo_CM[v]  = CM_Pos_Est[v]/(double)NParCMEst;
            Subhalo_vel[v] = CM_Vel_Est[v]/(double)NParCMEst;
         }
         //Aux_Message( stdout, "      (MPI_Rank = %10d) After: Subhalo_CM[0] = %14.7e,  Subhalo_CM[1] = %14.7e,  Subhalo_CM[2] = %14.7e\n", MPI_Rank, Subhalo_CM[0], Subhalo_CM[1], Subhalo_CM[2]);
         //Aux_Message( stdout, "      (MPI_Rank = %10d) After: Subhalo_vel[0] = %14.7e, Subhalo_vel[1] = %14.7e, Subhalo_vel[2] = %14.7e\n", MPI_Rank, Subhalo_vel[0], Subhalo_vel[1], Subhalo_vel[2]);

//       free memory
         delete [] SortByEnergy_IdxTable_All;
      }
      MPI_Bcast( Subhalo_CM, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Bcast( Subhalo_vel, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

//    free memory
      for (int v=0; v<7; v++)  delete [] ParAllAtt_ThisRank[v];
      for (int v=0; v<7; v++)  delete [] ParAllAtt_AllRank[v];
   }

   //Aux_Message( stdout, "      Start of an iteration at MPI_Rank = %10d, NPar_Tot = %10d, %10d\n", MPI_Rank, ParN_ThisRank, amr->Par->NPar_Active);

// iteratively update CM estimates untill convergence
   if ( MPI_Rank == 0 )  Aux_Message( stdout, "      Self-potential Analysis Starts: \n");
   while ( true )
   {
//    increment the iteration counter
      ++NIter;
      //Aux_Message( stdout, "      Start of an iteration at MPI_Rank = %10d, NPar_Tot = %10d, %10d\n", MPI_Rank, ParN_ThisRank, amr->Par->NPar_Active);

//    ============================================================================================================
//    [Part II] update individual particle energy entries
      ParEnergy_Calculation_Anisotropy_Tidal_Stripping( true );

//    ============================================================================================================
//    [Part III] relabel bound/unbound particles, and perform non-parallel sorting based on the particle total specific energy
      NParBound_Sum = 0;                                // initialize/reset the value
      ParEnergy_Sorting_Anisotropy_Tidal_Stripping( NParBoundMin, NParBound_Sum_5percent, CM_Pos_Relative_Error, CM_Vel_Relative_Error );
      MPI_Barrier( MPI_COMM_WORLD );

      if ( MPI_Rank == 0 )
      {
         Aux_Message( stdout, "      NIter =  %10d,", NIter );
         Aux_Message( stdout, "      CM_Pos_Relative_Error = %14.7e, CM_Vel_Relative_Error = %14.7e,", CM_Pos_Relative_Error, CM_Vel_Relative_Error );
         Aux_Message( stdout, "      Bound ParN = %10d\n", NParBound_Sum);

         if ( (MAX(CM_Pos_Relative_Error, CM_Vel_Relative_Error) <= 0.0001  ||  NIter >= NIterMax) && !((NParBound_Sum < NParBoundMin) && FirstTime_MassRest) )
         {
//          [Part IV] output CM information before exiting the loop
            if ( FirstTime )
            {
               if ( Aux_CheckFileExist(filename_center) )
                  Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_center );
               else
               {
                  FILE *file_center = fopen( filename_center, "w" );
                  fprintf( file_center, "#%19s  %10s  %14s  %10s  %14s  %14s  %14s  %14s  %14s  %14s  %14s %14s\n",
                           "Time", "Step", "NPar_Bound", "NIter", "CM_x", "CM_y", "CM_z", "CM_Velx", "CM_Vely", "CM_Velz", "Pos_RelErr", "Vel_RelErr");
                  fclose( file_center );
               }

               if ( OPT_beta_compute )
               {
                  if ( Aux_CheckFileExist(filename_beta) )
                     Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_beta );
                  else
                  {
                     FILE *file_beta = fopen( filename_beta, "w" );
                     fprintf( file_beta, "#%19s  %10s  %14s  %14s  %14s  %14s\n",
                              "Time", "Step", "NPar_Bound", "NPar_Bound_Half", "beta_Total", "beta_Half");
                     fclose( file_beta );
                  }
               }

               FirstTime = false;
            }

            FILE *file_center = fopen( filename_center, "a" );
            fprintf( file_center, "%20.14e  %10ld  %14d  %10d  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e %14.7e\n",
                     Time[0], Step, NParBound_Sum, NIter, Subhalo_CM[0], Subhalo_CM[1], Subhalo_CM[2], Subhalo_vel[0], Subhalo_vel[1], Subhalo_vel[2], CM_Pos_Relative_Error, CM_Vel_Relative_Error );
            fclose( file_center );

//          exit the loop
            //Aux_Message( stdout, "      Exit the loop at Time = %20.14e\n", Time[0] );
         } // if ( MAX(CM_Pos_Relative_Error, CM_Vel_Relative_Error) <= 0.0001  ||  NIter >= NIterMax )
      } // if ( MPI_Rank == 0 )

//    exit the loop when applicable
      if ( MAX(CM_Pos_Relative_Error, CM_Vel_Relative_Error) <= 0.0001  ||  NIter >= NIterMax )
      {
//       restore all particle masses and repeat the calculations as necessary
         if ((NParBound_Sum < NParBoundMin) && FirstTime_MassRest)
         {
            for (long p=0; p<ParN_ThisRank; p++)
            {
               if ( (PMass[p] > (real)0.0) &&  (PType[p] != PTYPE_TRACER) )
                  {
                     amr->Par->Attribute[Idx_ParTag][p] = 1.0;
                     amr->Par->Mass[p]                  = Particle_Mass;
                  }
            }

            FirstTime_MassRest = false;
         }
         else  break;
      }
   } // while ( true )


// [Part V] re-enable the external potential and recompute potential on all grids
// Restore the particle mass of all unbound particles.
   for (long p=0; p<ParN_ThisRank; p++)
   {
      if ( (PMass[p] > (real)0.0 && PMass[p] < 0.5*Particle_Mass) &&  PType[p] != PTYPE_TRACER )
         {
            amr->Par->Mass[p] = Particle_Mass;
         }
   }

// Recompute the gravitational potential on all grids
   Recompute_Potential_Anisotropy_Tidal_Stripping( EXT_POT_FUNC );


// [Part VI] update the specific angular momentum AND total energy (+background) of all active particles
   for (long p=0; p<ParN_ThisRank; p++)
   {
//    skip inactive (with negative mass), and tracer particles
      if ( PMass[p] > (real)0.0  &&  PType[p] != PTYPE_TRACER )
      {
//       compute and record the specific angular momentum for all particles
         amr->Par->Attribute[Idx_ParAngMom][p] = SQRT(SQR((Pos[1][p]-Subhalo_CM[1])*(Vel[2][p]-Subhalo_vel[2])
                                                 - (Pos[2][p]-Subhalo_CM[2])*(Vel[1][p]-Subhalo_vel[1])) +
                                                 SQR((Pos[2][p]-Subhalo_CM[2])*(Vel[0][p]-Subhalo_vel[0])
                                                 - (Pos[0][p]-Subhalo_CM[0])*(Vel[2][p]-Subhalo_vel[2])) +
                                                 SQR((Pos[0][p]-Subhalo_CM[0])*(Vel[1][p]-Subhalo_vel[1])
                                                 - (Pos[1][p]-Subhalo_CM[1])*(Vel[0][p]-Subhalo_vel[0])));
      }
   }


   if ( OPT_beta_compute )
   {
//    [Part VII] compute and record velocity anisotropy of the bound remnant
      double rx_p               = 0.0;
      double ry_p               = 0.0;
      double rz_p               = 0.0;
      double r_norm_p           = 0.0;
      double v_radial_sq_p      = 0.0;
      double v_tangential_sq_p  = 0.0;
      int NParBound_Sum_Half    = NParBound_Sum/2;
      long double v_radial_sq_half       = 0.0;
      long double v_tangential_sq_half   = 0.0;
      long double v_radial_sq_total      = 0.0;
      long double v_tangential_sq_total  = 0.0;
      long double v_radial_sq_half_AllRank      = 0.0;
      long double v_tangential_sq_half_AllRank  = 0.0;
      long double v_radial_sq_total_AllRank     = 0.0;
      long double v_tangential_sq_total_AllRank = 0.0;


      for (long p=0; p<ParN_ThisRank; p++)
      {
//       skip inactive (with negative mass), and tracer particles
         if ( PMass[p] > (real)0.0  &&  PType[p] != PTYPE_TRACER )
         {
//          record radial and tangential components of the velocity squared for ALL bound particles
            rx_p              = Pos[0][p]-Subhalo_CM[0];
            ry_p              = Pos[1][p]-Subhalo_CM[1];
            rz_p              = Pos[2][p]-Subhalo_CM[2];
            r_norm_p          = SQRT(SQR(rx_p) + SQR(ry_p) + SQR(rz_p));
            rx_p             /= r_norm_p;
            ry_p             /= r_norm_p;
            rz_p             /= r_norm_p;
            v_radial_sq_p     = SQR(rx_p*(Vel[0][p]-Subhalo_vel[0]) + ry_p*(Vel[1][p]-Subhalo_vel[1]) + rz_p*(Vel[2][p]-Subhalo_vel[2]));
            v_tangential_sq_p = (SQR(Vel[0][p]-Subhalo_vel[0]) + SQR(Vel[1][p]-Subhalo_vel[1]) + SQR(Vel[2][p]-Subhalo_vel[2])) - v_radial_sq_p;

//          record radial and tangential components of the velocity squared of ALL bound particles
            if ( amr->Par->Attribute[Idx_ParTag][p] > (real)0.0 )
            {
               v_radial_sq_total       += v_radial_sq_p;
               v_tangential_sq_total   += v_tangential_sq_p;

//             record radial and tangential components of the velocity squared for 50% most bound particles
               if ( amr->Par->Attribute[Idx_ParTag][p] <= NParBound_Sum_Half )
               {
                  v_radial_sq_half     += v_radial_sq_p;
                  v_tangential_sq_half += v_tangential_sq_p;
               }
            }
         }
      }

      MPI_Barrier( MPI_COMM_WORLD );
      MPI_Reduce( &v_radial_sq_half, &v_radial_sq_half_AllRank, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &v_tangential_sq_half, &v_tangential_sq_half_AllRank, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &v_radial_sq_total, &v_radial_sq_total_AllRank, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &v_tangential_sq_total, &v_tangential_sq_total_AllRank, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

      if ( MPI_Rank == 0 )
      {
         long double beta_half              = 1.0-(v_tangential_sq_half_AllRank/(2.0*v_radial_sq_half_AllRank));
         long double beta_total             = 1.0-(v_tangential_sq_total_AllRank/(2.0*v_radial_sq_total_AllRank));

         FILE *file_beta = fopen( filename_beta, "a" );
         fprintf( file_beta, "%20.14e  %10ld  %14d  %14d  %14.7e  %14.7e\n",
                  Time[0], Step, NParBound_Sum, NParBound_Sum_Half, (double)beta_total, (double)beta_half );
         fclose( file_beta );
      } // if ( MPI_Rank == 0 )
   } // if ( OPT_beta_compute )


   ParEnergy_Calculation_Anisotropy_Tidal_Stripping( false );

   Sorting_Step = false;
} // FUNCTION : Output_UserWorkBeforeOutput_Anisotropy_Tidal_Stripping




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_Anisotropy_Tidal_Stripping
// Description :  Call Output_UserWorkBeforeOutput_Anisotropy_Tidal_Stripping() at every global time-step
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_Anisotropy_Tidal_Stripping()
{
   if ( Sorting_Step )   Output_UserWorkBeforeOutput_Anisotropy_Tidal_Stripping();
   Sorting_Step = true;

// if the automatic simulation termination check is enabled
   if ( OPT_Auto_Termination )
   {
//    terminate the system if the subhalo is disrupted
      if ( NParBound_Sum < Ter_NPar )    End_GAMER();
   }

} // FUNCTION : Aux_Record_Anisotropy_Tidal_Stripping




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Anisotropy_Tidal_Stripping
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Anisotropy_Tidal_Stripping()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();
   Init_ExtPot_Ptr         = Init_ExtPot_Anisotropy_Tidal_Stripping;
   Init_Function_User_Ptr  = SetGridIC;
   if (OPT_Auto_AMR)    Flag_User_Ptr = Flag_Anisotropy_Tidal_Stripping;

#  ifdef PARTICLE
   Par_Init_Attribute_User_Ptr = AddNewParticleAttribute_Anisotropy_Tidal_Stripping;
#  endif
#  endif // #if ( MODEL == HYDRO )

   Output_UserWorkBeforeOutput_Ptr = Output_UserWorkBeforeOutput_Anisotropy_Tidal_Stripping;
   Aux_Record_User_Ptr             = Aux_Record_Anisotropy_Tidal_Stripping;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Anisotropy_Tidal_Stripping
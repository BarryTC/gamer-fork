#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined GRAVITY )

extern double Central_CM[3];


//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_BH_Buoyancy
// Description :  Flag cells for refinement for the BH Buoyancy test problem
//
// Note        :  1. Linked to the function pointer "Flag_User_Ptr" by Init_TestProb_Hydro_BH_Buoyancy()
//                2. Please turn on the runtime option "OPT__FLAG_USER"
//                3. A cell will be flagged for refinement if any of the refinement criteria is satisfied and
//                   this function decides whether or not a cell on level lv will be refined to level lv+1
//                4. The function is modeled after "Flag_AGORA.cpp" and "Flag_Bondi.cpp"
//
// Parameter   :  i,j,k     : Indices of the targeted element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv        : Refinement level of the targeted patch
//                PID       : ID of the targeted patch
//                Threshold : User-provided threshold for the flag operation, which is loaded from the
//                            file "Input__Flag_User"
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_BH_Buoyancy( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

//  affect AMR only for "lv > 4"
    if ( lv > 4 )
    {
        const double dh     = amr->dh[lv];
        const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                                amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                                amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

        bool Flag = false;

//      flag cells ture (for refinement) if lying within the target radius at each level
        const double dr[3]   = { Pos[0]-Central_CM[0], Pos[1]-Central_CM[1], Pos[2]-Central_CM[2] };
        const double Radius  = SQRT( SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]) );
        Flag = Radius < Threshold[0];
//      Aux_Message( stdout, "      lv =  %10d; Threshold[0] = %14.7e; Threshold[lv] = %14.7e; Pos[0] =  %14.7e; Radius =  %14.7e\n", lv, Threshold[0], Threshold[lv], Pos[0], Radius );
        return Flag;
    }
//  leave patches with "lv < 4" unaffected
    else return false;

} // FUNCTION : Flag_BH_Buoyancy



#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )
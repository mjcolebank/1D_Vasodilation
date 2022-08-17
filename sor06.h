/***************************************************************************/
/*                                                                         */
/*  Program: sor06.h                                                       */
/*  Version: 2.0                                                           */
/*  By: Mette Olufsen, Math-Tech                                           */
/*  Date: 14. Jan. 1997                                                    */
/*                                                                         */
/*  This header file defines the global parameters                         */
/*                                                                         */
/***************************************************************************/

// $Id: sor06.h,v 1.8 2010-10-20 15:38:28 mette Exp $

#ifndef _SOR06_H
#define _SOR06_H

#include <cmath>
#include "tools.h"

int    nbrves;               // Number of vessels in the tree.

int    tmstps = 8192,                 // The number of timesteps per period.
       plts   = 512;                 // Number of plots per period.

int max_D = 4; // MJC: For passing bifurcation pointer
const int WALL_MODEL = 2; // Define which wall model to use.
//const char *CO_filename = "qC6_8192.dat";     // Input flow file at the heart for the control mouse.
const char *CO_filename = "qH1_8192.dat";     // Input flow file at the heart for the hypoxic mouse.
// const char *CO_filename = "qH1_16384.dat";     // Input flow file at the heart for the hypoxic mouse.


double conv   = 1332.220,              // Conversion from mmHg to SI-units.
       rho    = 1.055,                // Density of blood [g/cm^3].
       mu     = 0.049,                // Viscosity of blood [g/cm/s].
       mu_pl  = mu,                   // Viscosity of blood [g/cm/s].
       nu     = mu/rho,

       Tper   = 0.11,                 // The period of one heart beat [s].
                                      // for the friction.
       Lr     = 0.1,                  // Characteristic radius of the
                                      // vessels in the tree [cm].
       Lr2    = sq(Lr),               // The squared radius [cm2].
       Lr3    = cu(Lr),               // The radius to the third power [cm^3].
       g      = 981.0,                // The gravitational force [cm/s^2].
       q      = 10.0*Lr2,             // The characteristic flow [cm^3/s].
       Fr2    = sq(q)/g/pow(Lr,5),    // The squared Froudes number.
       Re     = q*rho/mu/Lr,          // Reynolds number.
       bound_thick = sqrt(nu*Tper/(2.0*M_PI))/Lr, // Boundary layer thickness (MJC)
       Period = Tper*q/Lr3,           // The dimension-less period.
       k      = Period/tmstps,        // Length of a timestep.
       Deltat = Period/plts,          // Interval between each point plottet.
       p0     = 0.0*rho/g/Lr*conv,   // Ensures a certain diastolic pressure.
    // MJC: for stenosis, not used
    inert_perm = 0.10,
    ALF = 1.0e-4,

       *fjac[24],   // Work space used by bound_bif.
       *fjac_ST[4], fvec_ST, // MJC: for Lax-Wendroff in bound_right
       xr, f, df;                     // Work space used by bound_right.

#endif

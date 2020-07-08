

/*
 * Include Files
 *
 */
#if defined(MATLAB_MEX_FILE)
#include "tmwtypes.h"
#include "simstruc_types.h"
#else
#include "rtwtypes.h"
#endif

/* %%%-SFUNWIZ_wrapper_includes_Changes_BEGIN --- EDIT HERE TO _END */
#include <math.h>
/* %%%-SFUNWIZ_wrapper_includes_Changes_END --- EDIT HERE TO _BEGIN */
#define u_width 1
#define y_width 1
/*
 * Create external references here.  
 *
 */
#include "MyCode\inc\EllipseFitting.h"
#include "MyCode\scr\EllipseFitting.c"
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
/* extern double func(double a); */
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Output functions
 *
 */
void SigIdentifyFunc_Outputs_wrapper(const real_T *RunCmd,
			const real_T *SignalA,
			const real_T *SignalB,
			real_T *Sin_drift,
			real_T *Cos_drift,
			real_T *Gain,
			real_T *Theta)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
/* This sample sets the output equal to the input
      y0[0] = u0[0]; 
 For complex signals use: y0[0].re = u0[0].re; 
      y0[0].im = u0[0].im;
      y1[0].re = u1[0].re;
      y1[0].im = u1[0].im;
*/
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
    gEllipseSC_AB.Cmd.bit.StartIteration = *RunCmd;
    SC_AngleChoose(&gEllipseSC_AB,*SignalA,*SignalB,ANGLE_LENGTH_60DG);
    EllipseFittingCalcu(&gEllipseSC_AB,*SignalA,*SignalB);
    *Sin_drift = gEllipseSC_AB.Sin_drift;
    *Cos_drift = gEllipseSC_AB.Cos_drift;
    *Gain = gEllipseSC_AB.Gain;
    *Theta = gEllipseSC_AB.Theta;
}

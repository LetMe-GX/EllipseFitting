/****************************************************************************
 * Funciton: Description: EllipseFitting For SC Encoder
 * Data: 20191029
 * History: 
****************************************************************************/
#ifndef ELLIPSE_FITTING_H
#define ELLIPSE_FITTING_H

#ifdef __cplusplus
extern "C" {
#endif

/* Includes -------------------------------------------------------------- */
#include "DataTypeDef.h"
    
/* Private typedef ------------------------------------------------------- */
typedef struct ELLIPSE_FITTING_COMMAND_STRUCT_DEF{
    u16 CalcuOver:1;          //1:end
    u16 Err:1;                //1:error
    u16 State:2;
    u16 PreSampleflag:1;
    u16 StartIteration:1;     //0:stop 1:run
    u16 reserved:10;
}ELLIPSE_FITTING_COMMAND_STRUCT;

typedef union ELLIPSE_FITTING_COMMAND_UNION_DEF{
    u16 all;
    ELLIPSE_FITTING_COMMAND_STRUCT bit;
}ELLIPSE_FITTING_COMMAND_UNION;
        
typedef struct SC_SIG_FITTING_STRUCT_DEF{
    ELLIPSE_FITTING_COMMAND_UNION Cmd;
    
    s64 Pre_Total[5];
    s16 Pre_Cnt;
    s32 Pre_Data[5];
    
    u16 Sec;
    u16 UseFlag;
    s16 SinShadow;
    s16 CosShadow;
    
    f32 A;
    f32 B;
    f32 C;
    f32 D;
    
    f32 P[16];
    f32 K[4];
    
    f32 x1;              //sin*cos
    f32 x2;              //cos*cos
    f32 x3;              //sin
    f32 x4;              //cos
    f32 y;               //-sin*sin
    
    s16 Sin_drift;       //unit:Q1  1->1
    s16 Cos_drift;       //unit:Q1  1->1
    s16 Gain;            //unit:Q12 1->1024
    s16 Theta;           //unit:Q16 2pi->65536
}SC_SIG_FITTING_STRUCT;

/* Private define -------------------------------------------------------- */
#define ANGLE_LENGTH_60DG    20000

/* Private function prototypes ------------------------------------------- */
extern SC_SIG_FITTING_STRUCT gEllipseSC_AB;

extern void EllipseFittingCalcu(SC_SIG_FITTING_STRUCT *p,s16 InputSin,s16 InputCos);
extern void SC_AngleChoose(SC_SIG_FITTING_STRUCT *p,s16 InputSin,s16 InputCos,s16 mThreshold);

#ifdef __cplusplus
}
#endif /* extern "C" */

#endif  //end of definition

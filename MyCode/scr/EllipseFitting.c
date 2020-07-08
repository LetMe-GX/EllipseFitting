/****************************************************************************
 * Funciton: Description: EllipseFitting For SC Encoder 
 *                       £¨Recursive Least Square Algorithm£©
 * Data: 20191029
 * History: 
****************************************************************************/

/* Includes -------------------------------------------------------------- */
#include "MyCode\inc\EllipseFitting.h"

/* Private variables ----------------------------------------------------- */
SC_SIG_FITTING_STRUCT gEllipseSC_AB;

/* Private functions ----------------------------------------------------- */

/****************************************************************************
 * Function Name: InitEllipseFitFunc
 * Description  : Parameters Initialization
 * Input        : p
 * Output       : None
 * Return       : None
****************************************************************************/
void InitEllipseFitFunc(SC_SIG_FITTING_STRUCT *p)
{
    u16 i;
    
    p->Cmd.all = 0;
    
    p->A = 0.0;
    p->B = 0.0;
    p->C = 0.0;
    p->D = 0.0;
    
    for(i=0;i<4;i++)
    {
        p->K[i] = 0;
    }
    
    for(i=0;i<16;i++)
    {
        if((i==0)||(i==5)||(i==10)||(i==15))
        {
            p->P[i] = 1;
        }
        else
        {
            p->P[i] = 0;
        }
    }
    
    p->Pre_Cnt = 0;
    for(i=0;i<5;i++)
    {
        p->Pre_Data[i] = 0;
        p->Pre_Total[i] = 0;
    }
    
    p->SinShadow = 0;
    p->CosShadow = 0;
    p->Sec = 0;
    p->UseFlag = 0;
}

/****************************************************************************
 * Function Name: SC_AngleChoose
 * Description  : Sampling Points Selection
 * Input        : p,InputSin,InputCos,mThreshold
 * Output       : None
 * Return       : None
****************************************************************************/
void SC_AngleChoose(SC_SIG_FITTING_STRUCT *p,s16 InputSin,s16 InputCos,s16 mThreshold)
{
    s16 mData1,mData2;
    
    if(p->UseFlag == 1)
    {
        return;
    }
    
    mData1 = abs(InputSin - p->SinShadow);
    mData2 = abs(InputCos - p->CosShadow);
    
    if(((mData1 > mThreshold)||(mData2 > mThreshold))
      ||(p->Cmd.bit.State != 2))
    {
        p->SinShadow = InputSin;
        p->CosShadow = InputCos;
        p->UseFlag = 1;
    }
}

/****************************************************************************
 * Function Name: PreSampleForSC
 * Description  : Presampling(For Reducing Order)
 * Input        : p,InputSin,InputCos
 * Output       : None
 * Return       : None
****************************************************************************/
void PreSampleForSC(SC_SIG_FITTING_STRUCT *p,s16 InputSin,s16 InputCos)
{
    s32 mDataA,mDataB;
    s16 i;
    
    mDataA = (s32)InputSin;
    mDataB = (s32)InputCos;
    
    if(p->Pre_Cnt >= 64)
    {
        for(i=0;i<5;i++)
        {
            p->Pre_Data[i] = p->Pre_Total[i] >> 6;
            p->Pre_Total[i] = 0;
            p->Cmd.bit.PreSampleflag = 1;
        }
        p->Pre_Cnt = 0;
    }
    else
    {
        p->Pre_Total[0] = p->Pre_Total[0] - mDataA * mDataA;
        p->Pre_Total[1] = p->Pre_Total[1] - mDataA * mDataB;
        p->Pre_Total[2] = p->Pre_Total[2] - mDataB * mDataB;
        p->Pre_Total[3] = p->Pre_Total[3] - mDataA;
        p->Pre_Total[4] = p->Pre_Total[4] - mDataB;
        p->Pre_Cnt++;
    }
}

/****************************************************************************
 * Function Name: Ellipse_K_Calcu
 * Description  : Calculate K Value 
 * Input        : p,InputSin,InputCos
 * Output       : None
 * Return       : None
****************************************************************************/
void Ellipse_K_Calcu(SC_SIG_FITTING_STRUCT *p,s16 InputSin,s16 InputCos)
{
    f32 C_Sin,D_Cos;
    f32 m_P_Matrix[4];
    f32 m_Data;
    u16 i;
    
    C_Sin = (f32)InputSin;
    D_Cos = (f32)InputCos;
    
    p->x1 = C_Sin * D_Cos + (f32)p->Pre_Data[1];
    p->x2 = D_Cos * D_Cos + (f32)p->Pre_Data[2];
    p->x3 = C_Sin + (f32)p->Pre_Data[3];
    p->x4 = D_Cos + (f32)p->Pre_Data[4];
    p->y = -C_Sin * C_Sin - (f32)p->Pre_Data[0];
    
    for(i=0;i<4;i++)
    {
        m_P_Matrix[i] = p->P[i*4] * p->x1;
        m_P_Matrix[i] = m_P_Matrix[i] + p->P[i*4+1] * p->x2;
        m_P_Matrix[i] = m_P_Matrix[i] + p->P[i*4+2] * p->x3;
        m_P_Matrix[i] = m_P_Matrix[i] + p->P[i*4+3] * p->x4;
    }
    
    m_Data = 1 + m_P_Matrix[0] * p->x1;
    m_Data = m_Data + m_P_Matrix[1] * p->x2;
    m_Data = m_Data + m_P_Matrix[2] * p->x3;
    m_Data = m_Data + m_P_Matrix[3] * p->x4;
    
    if((m_Data < 0.00001)&&(m_Data > -0.00001))
    {
        for(i=0;i<4;i++)
        {
            p->K[i] = 0.0;
        }
    }
    else
    {
        for(i=0;i<4;i++)
        {
            p->K[i] = m_P_Matrix[i] / m_Data;
        }
    }
}

/****************************************************************************
 * Function Name: Ellipse_Output_Calcu
 * Description  : Calculate A,B,C,D
 * Input        : p
 * Output       : None
 * Return       : None
****************************************************************************/
void Ellipse_Output_Calcu(SC_SIG_FITTING_STRUCT *p)
{
    f32 mDeltaData;
    f32 mPData[4];
    u16 i;
    
    mDeltaData = p->A * p->x1;
    mDeltaData = mDeltaData + p->B * p->x2;
    mDeltaData = mDeltaData + p->C * p->x3;
    mDeltaData = mDeltaData + p->D * p->x4;
    mDeltaData = p->y - mDeltaData;
    
    p->A = p->A + p->K[0] * mDeltaData;
    p->B = p->B + p->K[1] * mDeltaData;
    p->C = p->C + p->K[2] * mDeltaData;
    p->D = p->D + p->K[3] * mDeltaData;
    
    for(i=0;i<4;i++)
    {
        mPData[i] = p->P[i] * p->x1;
        mPData[i] = mPData[i] + p->P[i+4] * p->x2;
        mPData[i] = mPData[i] + p->P[i+8] * p->x3;
        mPData[i] = mPData[i] + p->P[i+12] * p->x4;
    }
    
    for(i=0;i<4;i++)
    {
        p->P[i*4+0] = p->P[i*4+0] - p->K[i] * mPData[0];
        p->P[i*4+1] = p->P[i*4+1] - p->K[i] * mPData[1];
        p->P[i*4+2] = p->P[i*4+2] - p->K[i] * mPData[2];
        p->P[i*4+3] = p->P[i*4+3] - p->K[i] * mPData[3];
    }
}

/****************************************************************************
 * Function Name: SCEncoderCalcu
 * Description  : Calculate drift,gain,phase 
 * Input        : p 
 * Output       : None
 * Return       : None
****************************************************************************/
void SCEncoderCalcu(SC_SIG_FITTING_STRUCT *p)
{
    f32 m_Data;
    u32 m_DataInt;
    
    m_Data = p->A * p->A - 4.0 * p->B;
    if(abs(m_Data) < 0.00001)
    {
        p->Sin_drift = 0;
        p->Cos_drift = 0;
        p->Cmd.bit.Err = 1;
    }
    else
    {
        p->Sin_drift = (s16)((2.0 * p->B * p->C - p->A * p->D) / m_Data);
        p->Cos_drift = (s16)((2.0 * p->D - p->A * p->C) / m_Data);
    }
    
    m_DataInt = abs(p->B * 1048576.0);
    p->Gain = sqrt(m_DataInt);    // 1/1024  //DSP-qsqrt  Matlab-sqrt
    
    if(p->Gain == 0)
    {
        p->Theta = 0;
        p->Cmd.bit.Err = 1;
    }
    else
    {
        m_Data = (-512.0) * p->A;
        m_Data = m_Data / (f32)p->Gain;
        
        if(abs(m_Data) < 0.20943951)
        {
            p->Theta = (s16)(m_Data * 10430.37835);   // rad->Q16(360)
        }
        else
        {
            p->Theta = 0;
        }
    }
}

/****************************************************************************
 * Function Name: EllipseFittingCalcu
 * Description  : Scheduling
 * Input        : p,InputSin,InputCos
 * Output       : None
 * Return       : None
****************************************************************************/
void EllipseFittingCalcu(SC_SIG_FITTING_STRUCT *p,s16 InputSin,s16 InputCos)
{
    switch(p->Cmd.bit.State)
    {
        case 0:
            InitEllipseFitFunc(p);
            p->Cmd.bit.State = 1;
            break;
        case 1:
            PreSampleForSC(p,InputSin,InputCos);
            p->UseFlag = 0;
            if(p->Cmd.bit.StartIteration != 0)
            {
                p->Cmd.bit.State = 2;
            }
            break;
        case 2:
            if(p->UseFlag == 1)
            {
                Ellipse_K_Calcu(p,InputSin,InputCos);
                Ellipse_Output_Calcu(p);
                p->UseFlag = 0;
            }
            if(p->Cmd.bit.StartIteration == 0)
            {
                p->Cmd.bit.State = 3;
            }
            break;
        case 3:
            SCEncoderCalcu(p);
            p->Cmd.bit.CalcuOver = 1;
            p->Cmd.bit.State = 0;
            break;
    }
}

/***************************** END OF FILE *********************************/

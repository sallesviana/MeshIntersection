#include "math.h"

#include "WolframRTL.h"

static WolframCompileLibrary_Functions funStructCompile;

static void * E0 = 0;


static mint I0_0;

static mbool initialize = 1;

#include "list.h"

DLLEXPORT int Initialize_list(WolframLibraryData libData)
{
if( initialize)
{
funStructCompile = libData->compileLibraryFunctions;
I0_0 = (mint) 3;
{
E0 = funStructCompile->getExpressionFunctionPointer(libData, "Hold[Function[List[t00x, t00y, t00z, t01x, t01y, t01z, t02x, t02y, t02z, r00x, r00y, r00z, r01x, r01y, r01z], list]]");
}
if( E0 == 0)
{
return LIBRARY_FUNCTION_ERROR;
}
initialize = 0;
}
return 0;
}

DLLEXPORT void Uninitialize_list(WolframLibraryData libData)
{
if( !initialize)
{
initialize = 1;
}
}

DLLEXPORT int list(WolframLibraryData libData, mreal A1, mreal A2, mreal A3, mreal A4, mreal A5, mreal A6, mreal A7, mreal A8, mreal A9, mreal A10, mreal A11, mreal A12, mreal A13, mreal A14, mreal A15, MTensor *Res)
{
mreal R0_0;
mreal R0_1;
mreal R0_2;
mreal R0_3;
mreal R0_4;
mreal R0_5;
mreal R0_6;
mreal R0_7;
mreal R0_8;
mreal R0_9;
mreal R0_10;
mreal R0_11;
mreal R0_12;
mreal R0_13;
mreal R0_14;
mreal R0_15;
mreal R0_16;
mreal R0_17;
mreal R0_18;
mreal R0_19;
mreal R0_20;
mreal R0_21;
mreal R0_22;
mreal R0_23;
mreal R0_24;
mreal R0_25;
mreal R0_26;
mreal R0_27;
mreal R0_28;
mreal R0_29;
mreal R0_30;
mreal R0_31;
mreal R0_32;
mreal R0_33;
mreal R0_34;
mreal R0_35;
mreal R0_36;
mreal R0_37;
mreal R0_38;
mreal R0_39;
mreal R0_40;
mreal R0_41;
mreal R0_42;
mreal R0_43;
mreal R0_44;
mreal R0_45;
mreal R0_46;
mreal R0_47;
mreal R0_48;
mreal R0_49;
mreal R0_50;
mreal R0_51;
mreal R0_52;
mreal R0_53;
mreal R0_54;
mreal R0_55;
MTensor* T0_0;
MTensor* T0_1;
MTensorInitializationData Tinit;
int err = 0;
Tinit = funStructCompile->GetInitializedMTensors(libData, 2);
T0_0 = MTensorInitializationData_getTensor(Tinit, 0);
T0_1 = MTensorInitializationData_getTensor(Tinit, 1);
R0_0 = A1;
R0_1 = A2;
R0_2 = A3;
R0_3 = A4;
R0_4 = A5;
R0_5 = A6;
R0_6 = A7;
R0_7 = A8;
R0_8 = A9;
R0_9 = A10;
R0_10 = A11;
R0_11 = A12;
R0_12 = A13;
R0_13 = A14;
R0_14 = A15;
{
int S0[15];
void * S1[15];
S0[0] = 3;
S1[0] = (void*) (&R0_0);
S0[1] = 3;
S1[1] = (void*) (&R0_1);
S0[2] = 3;
S1[2] = (void*) (&R0_2);
S0[3] = 3;
S1[3] = (void*) (&R0_3);
S0[4] = 3;
S1[4] = (void*) (&R0_4);
S0[5] = 3;
S1[5] = (void*) (&R0_5);
S0[6] = 3;
S1[6] = (void*) (&R0_6);
S0[7] = 3;
S1[7] = (void*) (&R0_7);
S0[8] = 3;
S1[8] = (void*) (&R0_8);
S0[9] = 3;
S1[9] = (void*) (&R0_9);
S0[10] = 3;
S1[10] = (void*) (&R0_10);
S0[11] = 3;
S1[11] = (void*) (&R0_11);
S0[12] = 3;
S1[12] = (void*) (&R0_12);
S0[13] = 3;
S1[13] = (void*) (&R0_13);
S0[14] = 3;
S1[14] = (void*) (&R0_14);
err = funStructCompile->evaluateFunctionExpression(libData, E0, 0, 0, 15, S0, S1, 3, 2, (void*) T0_1);
if( err)
{
goto error_label;
}
}
err = funStructCompile->MTensor_getMTensorInitialized(T0_0, *T0_1, &I0_0, 1);
if( err)
{
goto error_label;
}
funStructCompile->MTensor_copy(Res, *T0_0);
error_label:
funStructCompile->ReleaseInitializedMTensors(Tinit);
funStructCompile->WolframLibraryData_cleanUp(libData, 1);
return err;
}


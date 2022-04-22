uk_0_eval(x)= @evalpoly x 1.0 

uk_1_eval(x)= @evalpoly x 0.0  0.125  0.0  -0.20833333333333334 

uk_2_eval(x)= @evalpoly x 0.0  0.0  0.0703125  0.0  -0.4010416666666667  0.0  0.3342013888888889 

uk_3_eval(x)= @evalpoly x 0.0  0.0  0.0  0.0732421875  0.0  -0.8912109375  0.0  1.8464626736111112  0.0  -1.0258125964506173 

uk_4_eval(x)= @evalpoly x 0.0  0.0  0.0  0.0  0.112152099609375  0.0  -2.3640869140625  0.0  8.78912353515625  0.0  -11.207002616222995  0.0  4.669584423426247 

uk_5_eval(x)= @evalpoly x 0.0  0.0  0.0  0.0  0.0  0.22710800170898438  0.0  -7.368794359479631  0.0  42.53499874538846  0.0  -91.81824154324003  0.0  84.63621767460074  0.0  -28.212072558200244 

uk_6_eval(x)= @evalpoly x 0.0  0.0  0.0  0.0  0.0  0.0  0.5725014209747314  0.0  -26.491430486951554  0.0  218.1905117442116  0.0  -699.5796273761327  0.0  1059.9904525279999  0.0  -765.2524681411816  0.0  212.5701300392171 

uk_7_eval(x)= @evalpoly x 0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.7277275025844574  0.0  -108.09091978839464  0.0  1200.9029132163525  0.0  -5305.646978613405  0.0  11655.393336864536  0.0  -13586.550006434136  0.0  8061.722181737308  0.0  -1919.4576623184068 

uk_8_eval(x)= @evalpoly x 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  6.074042001273483  0.0  -493.915304773088  0.0  7109.514302489364  0.0  -41192.65496889756  0.0  122200.46498301747  0.0  -203400.17728041555  0.0  192547.0012325315  0.0  -96980.5983886375  0.0  20204.29133096615 

uk_9_eval(x)= @evalpoly x 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  24.380529699556064  0.0  -2499.830481811209  0.0  45218.76898136274  0.0  -331645.1724845636  0.0  1.2683652733216248e6  0.0  -2.813563226586534e6  0.0  3.763271297656404e6  0.0  -2.998015918538106e6  0.0  1.311763614662977e6  0.0  -242919.18790055133 

const uk_10_coef=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 110.01714026924674, 0.0, -13886.089753717039, 0.0, 308186.40461266245, 0.0, -2.785618128086455e6, 0.0, 1.328876716642182e7, 0.0, -3.756717666076335e7, 0.0, 6.634451227472903e7, 0.0, -7.410514821153264e7, 0.0, 5.095260249266463e7, 0.0, -1.970681911843222e7, 0.0, 3.2844698530720375e6, ]
uk_10_eval(x) = evalpoly(x, uk_10_coef)

const uk_11_coef=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 551.3358961220206, 0.0, -84005.43360302408, 0.0, 2.24376817792245e6, 0.0, -2.4474062725738734e7, 0.0, 1.420629077975331e8, 0.0, -4.958897842750303e8, 0.0, 1.1068428168230145e9, 0.0, -1.621080552108337e9, 0.0, 1.5535968995705795e9, 0.0, -9.39462359681578e8, 0.0, 3.255730741857656e8, 0.0, -4.932925366450995e7, ]
uk_11_eval(x) = evalpoly(x, uk_11_coef)

const uk_12_coef=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3038.0905109223845, 0.0, -549842.3275722886, 0.0, 1.739510755397817e7, 0.0, -2.2510566188941535e8, 0.0, 1.5592798648792577e9, 0.0, -6.563293792619284e9, 0.0, 1.79542137311556e10, 0.0, -3.302659974980072e10, 0.0, 4.128018557975397e10, 0.0, -3.463204338815877e10, 0.0, 1.868820750929582e10, 0.0, -5.866481492051846e9, 0.0, 8.14789096118312e8, ]
uk_12_eval(x) = evalpoly(x, uk_12_coef)

const uk_13_coef=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18257.75547429318, 0.0, -3.871833442572612e6, 0.0, 1.4315787671888906e8, 0.0, -2.167164983223796e9, 0.0, 1.763473060683497e10, 0.0, -8.786707217802325e10, 0.0, 2.879006499061506e11, 0.0, -6.453648692453765e11, 0.0, 1.008158106865382e12, 0.0, -1.098375156081223e12, 0.0, 8.19218669548577e11, 0.0, -3.990961752244664e11, 0.0, 1.1449823773202577e11, 0.0, -1.4679261247695614e10, ]
uk_13_eval(x) = evalpoly(x, uk_13_coef)

const uk_14_coef=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 118838.42625678328, 0.0, -2.918838812222081e7, 0.0, 1.247009293512711e9, 0.0, -2.1822927757529232e10, 0.0, 2.0591450323241003e11, 0.0, -1.1965528801961816e12, 0.0, 4.612725780849132e12, 0.0, -1.2320491305598287e13, 0.0, 2.334836404458184e13, 0.0, -3.1667088584785152e13, 0.0, 3.056512551993531e13, 0.0, -2.051689941093443e13, 0.0, 9.109341185239896e12, 0.0, -2.406297900028503e12, 0.0, 2.8646403571767896e11, ]
uk_14_eval(x) = evalpoly(x, uk_14_coef)

const uk_15_coef=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 832859.3040162894, 0.0, -2.3455796352225152e8, 0.0, 1.1465754899448242e10, 0.0, -2.2961937296824658e11, 0.0, 2.4850009280340854e12, 0.0, -1.663482472489248e13, 0.0, 7.437312290867914e13, 0.0, -2.3260483118893994e14, 0.0, 5.230548825784446e14, 0.0, -8.574610329828949e14, 0.0, 1.0269551960827622e15, 0.0, -8.894969398810261e14, 0.0, 5.427396649876595e14, 0.0, -2.2134963870252512e14, 0.0, 5.417751075510603e13, 0.0, -6.019723417234003e12, ]
uk_15_eval(x) = evalpoly(x, uk_15_coef)

const uk_16_coef=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.252951493434797e6, 0.0, -2.0016469281917765e9, 0.0, 1.1099740513917906e11, 0.0, -2.5215584749128555e12, 0.0, 3.1007436472896465e13, 0.0, -2.3665253045164925e14, 0.0, 1.2126758042503475e15, 0.0, -4.3793258383640155e15, 0.0, 1.1486706978449752e16, 0.0, -2.226822513391114e16, 0.0, 3.213827526858623e16, 0.0, -3.4447226006485136e16, 0.0, 2.705471130619707e16, 0.0, -1.5129826322457674e16, 0.0, 5.705782159023669e15, 0.0, -1.301012723549699e15, 0.0, 1.3552215870309362e14, ]
uk_16_eval(x) = evalpoly(x, uk_16_coef)

const uk_17_coef=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0069589531988926e7, 0.0, -1.807822038465807e10, 0.0, 1.1287091454108745e12, 0.0, -2.886383763141477e13, 0.0, 4.000444570430363e14, 0.0, -3.450385511846272e15, 0.0, 2.0064271476309532e16, 0.0, -8.270945651585064e16, 0.0, 2.4960365126160426e17, 0.0, -5.62631788074636e17, 0.0, 9.575335098169137e17, 0.0, -1.233611693196069e18, 0.0, 1.1961991142756303e18, 0.0, -8.592577980317544e17, 0.0, 4.4347954614171885e17, 0.0, -1.5552983504313898e17, 0.0, 3.3192764720355212e16, 0.0, -3.2541926196426675e15, ]
uk_17_eval(x) = evalpoly(x, uk_17_coef)

const uk_18_coef=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.2593921650476694e8, 0.0, -1.7228323871735056e11, 0.0, 1.2030115826419195e13, 0.0, -3.4396530474307606e14, 0.0, 5.33510697870884e15, 0.0, -5.160509319348522e16, 0.0, 3.376676249790609e17, 0.0, -1.5736434765189596e18, 0.0, 5.402894876715981e18, 0.0, -1.3970803516443374e19, 0.0, 2.7572829816505184e19, 0.0, -4.1788614446568374e19, 0.0, 4.859942729324835e19, 0.0, -4.301555703831442e19, 0.0, 2.8465212251676553e19, 0.0, -1.3639420410571586e19, 0.0, 4.4702009640123085e18, 0.0, -8.966114215270461e17, 0.0, 8.301957606731907e16, ]
uk_18_eval(x) = evalpoly(x, uk_18_coef)

const uk_19_coef=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.836255180230434e9, 0.0, -1.7277040123530002e12, 0.0, 1.3412416915180642e14, 0.0, -4.2619355104269e15, 0.0, 7.351663610930973e16, 0.0, -7.921651119323832e17, 0.0, 5.789887667664653e18, 0.0, -3.0255665989903716e19, 0.0, 1.1707490535797255e20, 0.0, -3.434621399768417e20, 0.0, 7.756704953461136e20, 0.0, -1.3602037772849937e21, 0.0, 1.8571089321463448e21, 0.0, -1.9677247077053117e21, 0.0, 1.601689857369359e21, 0.0, -9.824438427689853e20, 0.0, 4.39279220088871e20, 0.0, -1.3512175034359957e20, 0.0, 2.556380296052923e19, 0.0, -2.242438856186774e18, ]
uk_19_eval(x) = evalpoly(x, uk_19_coef)

const uk_20_coef=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.646840080706557e10, 0.0, -1.8187262038511043e13, 0.0, 1.5613123930484675e15, 0.0, -5.484033603883292e16, 0.0, 1.0461721131134348e18, 0.0, -1.2483700995047234e19, 0.0, 1.0126774169536592e20, 0.0, -5.891794135069496e20, 0.0, 2.548961114664971e21, 0.0, -8.40591581710835e21, 0.0, 2.1487414815055883e22, 0.0, -4.302534303482378e22, 0.0, 6.7836616429518815e22, 0.0, -8.423222750084318e22, 0.0, 8.194331005435126e22, 0.0, -6.173206302884411e22, 0.0, 3.5284358439034075e22, 0.0, -1.478774352843361e22, 0.0, 4.285296082829493e21, 0.0, -7.671943936729004e20, 0.0, 6.393286613940834e19, ]
uk_20_eval(x) = evalpoly(x, uk_20_coef)

const UK_POLYS = [uk_0_eval, uk_1_eval, uk_2_eval, uk_3_eval, uk_4_eval,
    uk_5_eval, uk_6_eval, uk_7_eval, uk_8_eval, uk_9_eval, uk_10_eval,
    uk_11_eval, uk_12_eval, uk_13_eval, uk_14_eval, uk_15_eval, uk_16_eval,
    uk_17_eval, uk_18_eval, uk_19_eval, uk_20_eval]

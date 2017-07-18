import FWCore.ParameterSet.Config as cms

# coefficient evaluated from SinglePhoton pT=50 GeV considering all the layers in EE
dEdX_weights = cms.vdouble(0.0,
                           0.0158115,
                           0.0150911,
                           0.0232459,
                           0.00583132,
                           0.00740258,
                           0.00878334,
                           0.00543954,
                           0.00740718,
                           0.0114131,
                           0.00791042,
                           0.0101713,
                           0.00819247,
                           0.0121763,
                           0.00709762,
                           0.0122922,
                           0.00735609,
                           0.0136739,
                           0.00694588,
                           0.0139695,
                           0.00701347,
                           0.016785,
                           0.00965742,
                           0.0110348,
                           0.00789397,
                           0.0132299,
                           0.0134486,
                           0.043594,
                           0.0477878,
                           0.,
                           0.,
                           0.,
                           0.,
                           0.,
                           0.,
                           0.,
                           0.,
                           0.,
                           0.,
                           0.,
                           0.
                           )

# coefficient evaluated from SinglePhoton pT=50 GeV considering only the trigger layer (i.e. only odd layer numbers in EE)
dEdX_weights_trgLayer = cms.vdouble(0.0,
                                    0.0183664,
                                    0.,
                                    0.0353655,
                                    0.,
                                    0.0131052,
                                    0.,
                                    0.0119768,
                                    0.,
                                    0.0193766,
                                    0.,
                                    0.0170963,
                                    0.,
                                    0.0191943,
                                    0.,
                                    0.0189112,
                                    0.,
                                    0.0202912,
                                    0.,
                                    0.0211648,
                                    0.,
                                    0.0247953,
                                    0.,
                                    0.020009,
                                    0.,
                                    0.0196586,
                                    0.,
                                    0.0447745,
                                    0.,
                                    0.0603181,
                                    0.,
                                    0.,
                                    0.,
                                    0.,
                                    0.,
                                    0.,
                                    0.,
                                    0.,
                                    0.,
                                    0.,
                                    0.,
                                    0.
                                    )

# coefficient evaluated from SinglePion pT=50 GeV considering all the layers
dEdx_weights_hadron = cms.vdouble(0.0,
                                  0.030774,
                                  0.0222183,
                                  0.0682369,
                                  0.0604082,
                                  0.0332034,
                                  0.0110801,
                                  0.0278163,
                                  0.0338501,
                                  -0.00667608,
                                  0.00536669,
                                  0.00591378,
                                  0.00566522,
                                  0.00928049,
                                  0.0140489,
                                  0.0193408,
                                  0.0161543,
                                  0.0151341,
                                  0.0233381,
                                  -0.00248578,
                                  0.0270719,
                                  0.0242825,
                                  0.0254389,
                                  0.0113491,
                                  0.0205752,
                                  0.0154181,
                                  0.0312033,
                                  0.0230816,
                                  0.0421248,
                                  0.0613817,
                                  0.113306,
                                  0.0915309,
                                  0.0844658,
                                  0.10493,
                                  0.11685,
                                  0.0845081,
                                  0.143039,
                                  0.0883324,
                                  0.109259,
                                  0.0849568,
                                  0.233339
                                  )


# coefficient evaluated from SinglePion pT=50 GeV considering only the trigger layer (i.e. only odd layer numbers in EE + all layers in the FH)
dEdX_weights_hadrons_trgLayer = cms.vdouble(0.0,
                                   0.0398721,
                                   0.,
                                   0.096048,
                                   0.,
                                   0.0800138,
                                   0.,
                                   0.0459438,
                                   0.,
                                   0.0125788,
                                   0.,
                                   0.00736786,
                                   0.,
                                   0.0142515,
                                   0.,
                                   0.0368964,
                                   0.,
                                   0.0381923,
                                   0.,
                                   0.0189196,
                                   0.,
                                   0.0446335,
                                   0.,
                                   0.0388062,
                                   0.,
                                   0.042066,
                                   0.,
                                   0.0304025,
                                   0.0456096,
                                   0.0678159,
                                   0.108399,
                                   0.105739,
                                   0.0651398,
                                   0.127217,
                                   0.0833317,
                                   0.0976129,
                                   0.155841,
                                   0.0565189,
                                   0.119126,
                                   0.239373
                                   )



# coefficient evaluated from SinglePion pT=20 GeV considering all the layers
#dEdx_weights_hadron = cms.vdouble(0.0,
#                                  0.0444183,
#                                  0.0235909,
#                                  0.0537158,
#                                  0.0176413,
#                                  0.0287962,
#                                  0.0335657,
#                                  0.0067575,
#                                  0.0213418,
#                                  -0.00445101,
#                                  0.0240096,
#                                  0.00856574,
#                                  0.0101149,
#                                  0.00727469,
#                                  0.0255024,
#                                  0.0180101,
#                                  0.0225153,
#                                  0.00711793,
#                                  0.0250128,
#                                  0.00467858,
#                                  0.0314869,
#                                  0.0212221,
#                                  0.0185467,
#                                  0.0247948,
#                                  0.0215014,
#                                  0.0225237,
#                                  0.0224248,
#                                  0.0360396,
#                                  0.05312,
#                                  0.0867038,
#                                  0.0907029,
#                                  0.0850328,
#                                  0.0856558,
#                                  0.097269,
#                                  0.0680684,
#                                  0.0976288,
#                                  0.0932654,
#                                  0.0911714,
#                                  0.0787013,
#                                  0.0968706,
#                                  0.15873                               
#)

# coefficient evaluated from SinglePion pT=50 GeV considering all the layers 
# relaxing also the TE (clusterization) threshold down to 1 transverse-mip
dEdx_weights_hadron_se5te1 = cms.vdouble(0.0,
                                         0.027708,
                                         0.0260267,
                                         0.0544282,
                                         0.0559544,
                                         0.0439143,
                                         0.0210649,
                                         0.0229192,
                                         0.00942202,
                                         0.00337046,
                                         0.000310888,
                                         0.000701685,
                                         0.00985089,
                                         0.014931,
                                         0.0147767,
                                         0.017622,
                                         0.0151429,
                                         0.0169176,
                                         0.0205255,
                                         -0.0020214,
                                         0.0209482,
                                         0.0233454,
                                         0.0277382,
                                         0.0141124,
                                         0.0185298,
                                         0.012053,
                                         0.0253528,
                                         0.0223149,
                                         0.0449065,
                                         0.0823653,
                                         0.0868177,
                                         0.057514,
                                         0.0778446,
                                         0.0832075,
                                         0.0640331,
                                         0.073178,
                                         0.0897546,
                                         0.0790907,
                                         0.0909747,
                                         0.0246653,
                                         0.256982
                                         )




# coefficient evaluated from Tau pT=50 GeV considering all the layers 
# dedicated to the decay mode: 1-prong
dEdx_weights_tau_1p = cms.vdouble(0.0,
                                         0.0263336,
                                         0.0207405,
                                         0.0699111,
                                         0.026906,
                                         0.037967,
                                         0.0196117,
                                         0.0102783,
                                         0.0132995,
                                         0.00953104,
                                         7.04488e-05,
                                         0.0037915,
                                         0.0219398,
                                         -0.00587652,
                                         0.0175544,
                                         0.0255551,
                                         0.0176376,
                                         0.0110859,
                                         0.0238254,
                                         0.0056354,
                                         0.0198025,
                                         0.0134421,
                                         0.0272093,
                                         0.0152909,
                                         0.033458,
                                         0.0128845,
                                         0.0302987,
                                         0.0213093,
                                         0.0451539,
                                         0.078241,
                                         0.0850126,
                                         0.0747919,
                                         0.0886989,
                                         0.0595063,
                                         0.0924794,
                                         0.0773141,
                                         0.0763301,
                                         0.117302,
                                         0.0581958,
                                         0.113706,
                                         0.18393              
                                         )

# coefficient evaluated from Tau pT=50 GeV considering all the layers 
# dedicated to the decay mode: 1-prong+piZeros
dEdx_weights_tau_1ppz = cms.vdouble(0.0,
                                         0.0251373,
                                         0.0281176,
                                         0.0376747,
                                         0.0254283,
                                         0.0200089,
                                         0.0058648,
                                         0.0167053,
                                         0.00199537,
                                         0.0203277,
                                         0.00290897,
                                         0.0212336,
                                         0.00315886,
                                         0.0218388,
                                         -0.0023309,
                                         0.0186982,
                                         0.00326528,
                                         0.0187691,
                                         0.000959069,
                                         0.0184387,
                                         0.00977819,
                                         0.0241289,
                                         0.0130964,
                                         0.0230427,
                                         0.0307277,
                                         0.0244709,
                                         0.0323741,
                                         0.0344123,
                                         0.0500692,
                                         0.0951707,
                                         0.0962284,
                                         0.0983047,
                                         0.0856944,
                                         0.102888,
                                         0.0811987,
                                         0.0993765,
                                         0.0998725,
                                         0.106711,
                                         0.0946023,
                                         0.127984,
                                         0.179042
                                         )


# coefficient evaluated from Tau pT=50 GeV considering all the layers 
# dedicated to the decay mode: 3-prongs
dEdx_weights_tau_3p = cms.vdouble(0.0,
                                         0.0680721,
                                         0.075155,
                                         0.0614739,
                                         0.0449299,
                                         0.0451136,
                                         0.0309278,
                                         0.00298436,
                                         0.0185999,
                                         0.00380071,
                                         0.0172034,
                                         0.0218986,
                                         0.00803686,
                                         0.0217498,
                                         0.0198373,
                                         0.0286634,
                                         0.0187655,
                                         0.0132061,
                                         0.0291525,
                                         0.0305885,
                                         0.0216399,
                                         0.0301103,
                                         0.0211445,
                                         0.042797,
                                         0.0278589,
                                         0.0181009,
                                         0.0328307,
                                         0.0408944,
                                         0.0613162,
                                         0.151007,
                                         0.114219,
                                         0.102304,
                                         0.13269,
                                         0.100985,
                                         0.143111,
                                         0.0953097,
                                         0.154416,
                                         0.115139,
                                         0.11433,
                                         0.127114,
                                         0.242156
                                         )







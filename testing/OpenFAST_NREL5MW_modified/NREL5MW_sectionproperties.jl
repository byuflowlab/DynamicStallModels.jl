using GXBeam



xaf = [1.00000000, 0.99619582, 0.98515158, 0.96764209, 0.94421447, 0.91510964, 0.88074158, 0.84177999, 0.79894110, 0.75297076, 0.70461763, 0.65461515, 0.60366461, 0.55242353, 0.50149950, 0.45144530, 0.40276150, 0.35589801, 0.31131449, 0.26917194, 0.22927064, 0.19167283, 0.15672257, 0.12469599, 0.09585870, 0.07046974, 0.04874337, 0.03081405, 0.01681379, 0.00687971, 0.00143518, 0.00053606, 0.00006572, 0.00001249, 0.00023032, 0.00079945, 0.00170287, 0.00354717, 0.00592084, 0.01810144, 0.03471169, 0.05589286, 0.08132751, 0.11073805, 0.14391397, 0.18067874, 0.22089879, 0.26433734, 0.31062190, 0.35933893, 0.40999990, 0.46204424, 0.51483073, 0.56767889, 0.61998250, 0.67114514, 0.72054815, 0.76758733, 0.81168064, 0.85227225, 0.88883823, 0.92088961, 0.94797259, 0.96977487, 0.98607009, 0.99640466, 1.00000000]
yaf = [0.00000000, 0.00017047, 0.00100213, 0.00285474, 0.00556001, 0.00906779, 0.01357364, 0.01916802, 0.02580144, 0.03334313, 0.04158593, 0.05026338, 0.05906756, 0.06766426, 0.07571157, 0.08287416, 0.08882939, 0.09329359, 0.09592864, 0.09626763, 0.09424396, 0.09023579, 0.08451656, 0.07727756, 0.06875796, 0.05918984, 0.04880096, 0.03786904, 0.02676332, 0.01592385, 0.00647946, 0.00370956, 0.00112514, -0.00046881, -0.00191488, -0.00329201, -0.00470585, -0.00688469, -0.00912202, -0.01720842, -0.02488211, -0.03226730, -0.03908459, -0.04503763, -0.04986836, -0.05338180, -0.05551392, -0.05636585, -0.05605816, -0.05472399, -0.05254383, -0.04969990, -0.04637175, -0.04264894, -0.03859653, -0.03433153, -0.02996944, -0.02560890, -0.02134397, -0.01726049, -0.01343567, -0.00993849, -0.00679919, -0.00402321, -0.00180118, -0.00044469, 0.00000000]

# Next, we define the chord length, the twist angle, and the pitch axis.  The pitch axis is the x-coordinate, 
# measured back from the leading edge, normalized by chord, that defines where the airfoil should be twisted about
# as seen in the figure below.  
# So a pitch axis of 0.25 means that the airfoil rotates about the quarter chord.  Positive twist is in the 
# direction of increasing angle of attack. In our case the twist is zero so the pitch axis is irrelevant though 
# we use the value noted in the paper.


chord = 1.9
twist = 0.0*pi/180
paxis = 0.4750 / chord


uni = Material(37.00e9, 9.00e9, 9.00e9, 4.00e9, 4.00e9, 4.00e9, 0.28, 0.28, 0.28, 1.86e3)
double = Material(10.30e9, 10.30e9, 10.30e9, 8.00e9, 8.00e9, 8.00e9, 0.30, 0.30, 0.30, 1.83e3)
gelcoat = Material(1e1, 1e1, 1e1, 1.0, 1.0, 1.0, 0.30, 0.30, 0.30, 1.83e3)
nexus = Material(10.30e9, 10.30e9, 10.30e9, 8.00e9, 8.00e9, 8.00e9, 0.30, 0.30, 0.30, 1.664e3)
balsa = Material(0.01e9, 0.01e9, 0.01e9, 2e5, 2e5, 2e5, 0.30, 0.30, 0.30, 0.128e3)
mat = [uni, double, gelcoat, nexus, balsa]


# The heart of the parmeterization is defining the ply sequence.  A single ply (or multiple plys of the same 
# material and orientation) is given by a Layer.  A Layer is defined by a material, a thickness, and a fiber 
# orientation. Thicknesses are absolute (not normalized), and orientation angles should be 
# in radians.  For example:

t = 0.001
theta = 20*pi/180
layer = Layer(balsa, t, theta)
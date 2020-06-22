#NDIM
2 
#PRIM
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
#SUPER
8 0 
0 8
#ORB
s0  0.0d0  0.0d0  0.0d0 #0
#HAMILT          tup  tdn   U
0 0 1.0 0.0 0.0  1.0  1.0  0.0 
0 0 0.0 1.0 0.0  1.0  1.0  0.0
0 0 0.0 0.0 0.0  0.0  0.0  4.0 
#SYMM
d  0.0d0 0.0d0 0.0d0 1.0d0 0.0d0 0.d0
d  0.0d0 0.0d0 0.0d0 0.0d0 1.0d0 0.d0
c4 0.0d0 0.0d0 0.0d0 0.0d0 0.0d0 1.d0
#PHASE
 1 1
-1 1
s0 0.0 0.0 0.0  1.0
s0 0.0 1.0 0.0 -1.0
# The first two entries in each line specify the orbital type where the up and
# down electrons are created. In this case we only have one orbital type and so
# both entry are the same and equal to 0 (the counter adopts a C-like convention 
# staring from 0). The next three entries (dx,dy,dz) are the displacement
# indicating the actual distance  between the orbitals on which the pairs is
# created.
#
# So the second line, for example, means that an the up and down
# electron will be created displaced by "1 0 0" both on orbital type 0. 
#
# What follows the "#" sign on each line is not read by the code but it is useful 
# for specifying the following PAIR field.
#
#BONDS   
0 0  0.0  0.0  0.0 # 1
0 0  1.0  0.0  0.0 # 2 -2
0 0  1.0  1.0  0.0 # 3 -3
0 0  0.0  1.0  0.0 # 4 -4
0 0 -1.0  1.0  0.0 # 5 -5
#
# This specifies how the bonds in the previous section have to be combined when
# computing pairing amplitude. The very first line lists the bonds the we want to use.
# In this case all of those specified above. Note that each of the entry in BONDS
# has two integers associated with it, apart from the case where the bond describes 
# creation and annihiliation *exactly* on the same site (as for the first entry above).
# In all the other case the positive integer (say 2,3,4 and 5) indicates 
# that we are creating an up electron on orbital type 0 and a down electron on
# an orbital (still of type 0) which is displaced by (dx, dy, dz) from the first. 
# The negative integers (-2,-3,-4 and -5) indicate that we are creating a down electron 
# on the given orbital type and an up electron on a site displaced by (-dx,-dy,-dz).
# Perhaps this is better explained if we just comment the next entries. Each line
# represents the actual linear combination of pair states that we want to compute the
# amplitude for. So, the first colum contains the label of the pairing function
# and its followed by the coefficient with which bonds have to be combined. Take, 
# for example, the d-wave. We create a pair states combining (1 0 0) (-1 0 0) with "+" 
# sign and (0 1 0) (0 -1 0) with "-" sign.
#PAIR     
           1     2    -2    3    -3     4     -4     5    -5
s-wave    1.0   0.0   0.0  0.0   0.0   0.0    0.0   0.0   0.0
s*-wave   0.0   1.0   1.0  0.0   0.0   1.0    1.0   0.0   0.0
s**-wave  0.0   0.0   0.0  1.0   1.0   0.0    0.0   1.0   1.0
d-wave    0.0   1.0   1.0  0.0   0.0  -1.0   -1.0   0.0   0.0
d*-wave   0.0   0.0   0.0  1.0   1.0   0.0    0.0  -1.0  -1.0
#END


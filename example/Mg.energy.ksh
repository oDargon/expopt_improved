#! /bin/ksh
################################################################################
#                                                                              #
# This is a sample shell script for running a MOLCAS job.                      #
#                                                                              #
#------------------------------------------------------------------------------#
#                                                                              #
# Author:  Per-Olof Widmark                                                    #
#          STS, IBM Sweden                                                     #
#                                                                              #
# Written: March 1992                                                          #
#                                                                              #
#------------------------------------------------------------------------------#
#                                                                              #
# History: none                                                                #
#                                                                              #
################################################################################
#set -v
#------------------------------------------------------------------------------#
# Job definitions                                                              #
#------------------------------------------------------------------------------#
#export MOLCAS_PROPERTIES=LONG
Home=`pwd`
WorkDir=$Home/tmp
Project='Mg'
export Project WorkDir
print 'Start of job:' $Project
print 'Current directory:' $Home
print 'Scratch directory:' $WorkDir
#------------------------------------------------------------------------------#
# Prologue                                                                     #
#------------------------------------------------------------------------------#
trap 'exit' ERR
rm -fr $WorkDir
mkdir $WorkDir
cd $WorkDir
#------------------------------------------------------------------------------#
# Input                                                                        #
#------------------------------------------------------------------------------#
Input="$Project.input"
cat <<EOF > $Input

>>> EXPORT DOCI=0

 &SEWARD  &END
Expert
Title
 X atom
Symmetry
 x y z
*RX2C
*R08O08
*NODKH
#DKH#
Finite
*OneOnly
*----------------------------------------------------------------
Basis set
Mg.... / inline 
#include <basis.17s12p.x2c.MB>
Mg    0.00000   0.00000   0.00000
End of basis
*----------------------------------------------------------------
End of input

*#include <basis>
*#include <basis.17s12p.nr.prim>
*#include <basis.17s12p.dkh22.prim>
*#include <basis.17s12p.dkh44.prim>
*#include <basis.17s12p.dkh66.prim>
*#include <basis.17s12p.dkh88.prim>
*#include <basis.17s12p.x2c.prim>
*#include <basis.17s12p.nr.MB>
*#include <basis.17s12p.dkh22.MB>
*#include <basis.17s12p.dkh44.MB>
*#include <basis.17s12p.dkh66.MB>
*#include <basis.17s12p.dkh88.MB>
*#include <basis.17s12p.x2c.MB>


 &GUESSORB &END
sThr
 0.0d0
tThr
 1.0d12
*PrMO's
* 4 5.0
*PrPopulation
End of Input


 &SCF &END
Occupied
  3  1  1  0  1  0  0  0
PrOrbitals
 2 1.0 3
End of input

>>COPY  Mg.ScfOrb Mg.1S.MO


>LINK -FORCE  Mg.1S.MO    NAT001
>LINK -FORCE  Mg.OneInt   ONE001
>LINK -FORCE  Mg.RunFile  RUN001

*>LINK -FORCE  Mg.1S.NO    NAT002
>LINK -FORCE  Mg.OneInt   ONE002
>LINK -FORCE  Mg.RunFile  RUN002

 &GENANO &END
Title
 Mg atom
sets
 1
Center
Mg
Weights
 1.000
Lift degeneracy
end of input

Orbital weights
 1.0 1.0


EOF
#------------------------------------------------------------------------------#
# Start executing molcas job                                                   #
#------------------------------------------------------------------------------#
cp $Home/BASIS.DAT basis
cp $Home/basis.* .
cp $Home/*.input .
#---
#cat $Input | include > tmp.input
#molcas tmp.input
#---
for DKH in NODKH  RX2C
do
   cat $Input | sed "s/#DKH#/$DKH/" | include > tmp.input
   molcas tmp.input
done
#------------------------------------------------------------------------------#
# Epilogue                                                                     #
#------------------------------------------------------------------------------#
cd -
# rm -r $WorkDir
echo "Calculation is finished!"
#------------------------------------------------------------------------------#
# Finished, now cleanup.                                                       #
#------------------------------------------------------------------------------#
exit

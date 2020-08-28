#! /bin/ksh
X=Mg
./$X.energy.ksh > $X.energy.log 2> $X.energy.err
cp $X.energy.log $X.energy.bak
#---
rm -f energy.tmp
grep -i energy $X.energy.log | sed 's/:/ /g' > tmp.energy.log
awk '/^\ *Total SCF energy/{print $4}' tmp.energy.log >> energy.tmp
awk '/^\ *Average CI energy/{print $4}' tmp.energy.log >> energy.tmp
# awk '/^\ *CI ENERGY/{print $3}' tmp.energy.log >> energy.tmp
# awk '/^\ *ACPF ENERGY/{print $3}' tmp.energy.log >> energy.tmp
awk '/^\ *SDCI ENERGY/{print $3}' tmp.energy.log >> energy.tmp
# awk '/^\ *CI CORRELATION ENERGY/{print $4}' tmp.energy.log >> energy.tmp
#---
awk -f average.awk energy.tmp > ENERGY.DAT
#---
cat energy.tmp >> energy.save
cat ENERGY.DAT >> energy.save
grep real $X.energy.err >> time.log
#---
exit 0

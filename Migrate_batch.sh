#!/bin/sh
for i in A1 A2 ML MO PN Y2:do
    cd /Volumes/home/Users/stsmall/Desktop/Wuchereria\ bancrofti/Wb_CO1/Analysis/Migrate/IInds/$i
    /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/migrate-3.3.2/migrate-n paramfile.within -nomenu
done

for i in A1_A2 A1_ML A1_MO A1_Y1 A1_Y2 A2_ML A2_MO A2_Y1 A2_Y2 ML_MO PN_A1 PN_A2 PN_ML PN_MO PN_Y1 PN_Y2 Y1_ML Y1_MO Y1_Y2 Y2_ML Y2_MO:do
   cd /Volumes/home/Users/stsmall/Desktop/Wuchereria\ bancrofti/Wb_CO1/Analysis/Migrate/Village/$i
    /Volumes/home/Users/stsmall/Desktop/PopGen_programs/Software/migrate-3.3.2/migrate-n paramfile.pairs -nomenu
done

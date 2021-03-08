
for ModeIx in  1 2 3 4 5 6 7 8 9 10
do
	../iMOD_v1.04_Linux/bin/imove SpikeTrimer.NoGlyc_model.pdb SpikeTrimer.NoGlyc_ic.evec Spike.NoGlyc_Mode${ModeIx}.pdb ${ModeIx}
done

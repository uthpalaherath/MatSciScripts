#!/usr/bin/env python

from abipy.waves.wfkfile import WfkFile


def readWFK(infile=None):
    """readWFK.

    This method reads the WFK.nc file from Abinit and
    stores the Wavefunction

    Parameters
    ----------
    infile :
        infile
    """

    wfk = abipy.waves.WfkFile("foo_WFK.nc")

    # Get a wavefunction.
    wave = wfk.get_wave(spin=0, kpoint=[0, 0, 0], band=0)


if __name__ == "__main__":
    readWFK("si_scf_WFK.nc")

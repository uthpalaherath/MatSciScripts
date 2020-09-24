#!/usr/bin/env python

""" Wannier Manifold Electron Occupation Calculater.

This is an extension to the original code developed by Xingyu Liao in Dr. Park's
group at UIC.
It calculates the total number of electrons in the Wannier manifold.
This is a requirement to know to perform DMFT calculations.

This requires that the DMFTwDFT/bin directory is added to $PYTHONPATH.

- Uthpala Herath

"""


from scipy import *
from os import path
from argparse import RawTextHelpFormatter
import argparse
import sys, subprocess, os
import re

import Struct
import VASP
from INPUT import *


class ElectronOccupation:
    """This class has the methods to calculate the electrons in the
    correlated subspace."""

    def __init__(self, args):

        self.dft = args.dft
        self.np = args.np


        # import the VASP class. This can be used for other DFT codes as well.
        self.create_DFTmu()
        self.DFT = VASP.VASP_class(dft=self.dft)

        # Initial guess for Fermi energy
        self.DFT.EFERMI = 7.0

        # Wannier parameters
        self.kmeshtol = 1e-06
        self.dos_kmesh = 64

        # create bands_plot directory
        if not os.path.exists("bands_plot"):
            os.makedirs("bands_plot")

        # Starting the calculation
        self.gen_win()
        self.read_num_wann()
        self.dft_run()
        self.update_win()
        self.run_wan90()
        self.postw90_run()
        self.calculator()

    def create_DFTmu(self):
        """
        This function creates a DFT_mu.out with an initial
        guess for the chemical potential. It will be updated once
        the DFT calculation is finished.
        """
        mu = 7.0

        if os.path.exists("DFT_mu.out"):
            os.remove("DFT_mu.out")
        f = open("DFT_mu.out", "w")
        f.write("%f" % mu)
        f.close()

    def gen_win(self):
        """
        This method generates wannier90.win for the initial DFT run.
        """
        # generating wannier90.win
        TB = Struct.TBstructure("POSCAR", p["atomnames"], p["orbs"])
        TB.Compute_cor_idx(p["cor_at"], p["cor_orb"])

        # Read number of bands from DFT input file
        try:
            if self.dft == "vasp":
                fi = open("INCAR", "r")
                data = fi.read()
                fi.close()
                self.DFT.NBANDS = int(
                    re.findall(r"\n\s*NBANDS\s*=\s*([\d\s]*)", data)[0]
                )
                print ("Number of bands read from INCAR = %d " % self.DFT.NBANDS)

        except:
            self.DFT.NBANDS = 100
            print ("WARNING: Number of bands not set in DFT input file!")

        # Setting num_bands in .win file.
        # If set to False num_bands is set to number of DFT bands.
        if list(p.keys()).count("num_bands_win"):
            if p["num_bands_win"]:
                self.wanbands = p["num_bands_win"]
            else:
                self.wanbands = self.DFT.NBANDS
        else:
            self.wanbands = self.DFT.NBANDS

        self.DFT.Create_win(
            TB,
            p["atomnames"],
            p["orbs"],
            p["L_rot"],
            self.wanbands,
            # Initially DFT.EFERMI is taken from DFT_mu.out but will
            # be updated later once the DFT calculation is complete.
            self.DFT.EFERMI + p["ewin"][0],
            self.DFT.EFERMI + p["ewin"][1],
            self.kmeshtol,
        )

        # If exclude_bands are to be included in the .win file.
        # Appending to current .win file.
        if list(p.keys()).count("exclude_bands"):
            if p["exclude_bands"]:
                f = open("wannier90.win", "a")
                f.write("\nexclude_bands :\t")
                f.write(", ".join(str(x) for x in p["exclude_bands"]))
                f.write("\n")
                f.close()

            else:
                pass
        else:
            pass

    def read_num_wann(self):
        """This reads the number of wannier bands from the generatied wannier.win
        file."""

        fi = open("wannier90.win", "r")
        data = fi.readlines()
        fi.close()

        for line in data:
            if re.match("num_wann", line):
                self.num_wann = line.split()[-1]

    def update_win(self):
        """
        This updates the wannier90.win file with the number of bands and fermi energy
        from the initial DFT calculation.
	"""

        # Updating wannier90.win with the number of DFT bands
        self.DFT.Read_NBANDS()
        self.DFT.Read_EFERMI()
        self.DFT.Update_win(
            self.DFT.NBANDS,
            self.DFT.EFERMI + p["ewin"][0],
            self.DFT.EFERMI + p["ewin"][1],
        )

        # Update wannier90.win file with the dos block.
        f = open("wannier90.win", "a")
        f.write("\ndos=true")
        f.write("\ndos_kmesh=%s" % self.dos_kmesh)
        f.write("\ndos_project=%s" % self.num_wann)
        f.close()

        print ("wannier90.win updated.")

        # find line number of dis_project in wannier90.win.
        fi = open("wannier90.win", "r")
        data = fi.readlines()
        fi.close()

        for ln, line in enumerate(data):
            if re.match("dos_project", line):
                self.dos_project_line = ln + 1

    def dft_run(self):
        """ DFT runner.
        This method performs the initial DFT calculation.
        """

        # VASP
        if self.dft == "vasp":
            # initial VASP run
            print ("Running VASP...")
            self.vasp_exec = "vasp_std"
            cmd = "mpirun -np "+str(self.np) + " " + self.vasp_exec  # + " > dft.out 2> dft.error"
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            if err:
                print ("DFT calculation failed! Check dft.error for details.\n")
                errdir = "dft.error"
                f = open(errdir, "wb")
                f.write(err)
                f.close()
                sys.exit()
            else:
                print ("DFT calculation complete.\n")
                outdir = "dft.out"
                f = open(outdir, "wb")
                f.write(out)
                f.close()

    def run_wan90(self, filename="wannier90"):
        """
        Running wannier90.x to generate .chk file.
        """

        print ("Running wannier90.x ...")
        cmd = "mpirun -np" + " " + str(self.np) + " " + "wannier90.x" + " " + filename
        out, err = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()
        if err:
            print ("wannier90 calculation failed!")
            print (err.decode("utf-8"))
            sys.exit()
        else:
            print ("wannier90 calculation complete.")
            print (out.decode("utf-8"))

    def postw90_run(self):
        NUM_OF_ORBT = int(self.num_wann)  ## NUMBER OF WANNIER BAND
        LINE_dos = self.dos_project_line  ### LINENUMBER WHICH DEFINE THE dos_project
        for i in range(NUM_OF_ORBT):
            print i + 1, "th iteration is running"
            fi = open("wannier90.win", "r")
            WIN = array(fi.readlines())
            fi.close()
            WIN[LINE_dos - 1] = "dos_project=" + str(i + 1) + "\n"
            fi = open("wannier90.win", "w")
            for j, winline in enumerate(WIN):
                fi.write(winline)
            fi.close()
            print "postw90 is running"
            cmd = "mpirun -np " + str(self.np) + " postw90.x wannier90"
            out, err = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            ).communicate()
            os.system("cp wannier90-dos.dat bands_plot/dos" + str(i + 1))
            print "-------------------------------------------------"

    def Load_DOS(self, filename, Fermi):
        fi = open(filename, "r")
        dos = fi.readlines()
        xy = zeros((len(dos), 2))
        for i, den in enumerate(dos):
            xy[i] = array(array(den.split()).astype(float))
            xy[i][0] -= Fermi
        fi.close()
        return xy

    def Integration(self, xy):
        total = 0.0
        for i in range(len(xy) - 1):
            total = (
                total + abs(xy[i + 1][0] - xy[i][0]) * (xy[i][1] + xy[i + 1][1]) / 2.0
            )
        return total

    def calculator(self):
        occ = 0.0
        for i in range(int(self.num_wann)):
            if not path.exists("dos" + str(i + 1)):
                print ("File does not exist!")
                exit()
            sepE = 0
            speI = 0
            dos = self.Load_DOS("dos" + str(i + 1), self.DFT.EFERMI)
            for j, ele in enumerate(dos):
                if ele[0] > sepE:
                    speI = j
                    break
            print i + 1, ":", self.Integration(dos[:speI, :])

            occ = self.Integration(dos[:speI, :]) + occ
        print ("occupancy is: ", occ)


if __name__ == "__main__":
    args = sys.argv[1:]
    if args:

        des = (
            "This script counts the total number of electrons in the Wannier sub space."
        )
        parser = argparse.ArgumentParser(
            description=des, formatter_class=RawTextHelpFormatter
        )

        parser.add_argument(
            "-dft",
            default="vasp",
            type=str,
            help="Choice of DFT code.",
            choices=["vasp", "siesta", "qe"],
        )

        parser.add_argument("-np", default=1, type=int, help="Number of processors.")

        args = parser.parse_args()
        ElectronOccupation(args)

    else:
        print ("Usage: electron_count.py -h")

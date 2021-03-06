#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 16:54:54 2019

@author: similarities
"""

import matplotlib.pyplot as plt
import numpy as np
import random
import math

marker_style = [',', '+', '.', 'v', '*', "P", "H", 'x', "D", ">"]


class Open_and_Plot_Picture:
    def __init__(self, filename, xmin, xmax, lambdaL, ROI_y, filedescription):

        self.filename = filename
        # px size full picture * usually 0 - 2048
        self.ymin = 0

        self.ymax = 2048
        # defined Roi in x to avoid boarder effects
        self.xmin = xmin

        self.xmax = xmax
        # integration ROI y for each HHG line
        self.ROI_y = ROI_y

        self.w0 = 6 * 1E-9  # radius of focus FWHM

        self.picture = np.empty([])

        self.integrated = np.empty([])

        self.x_backsubstracted = np.empty([2048, 2048])

        self.x_axis_in_nm = np.empty([2048, 1])

        self.x_axis_nm_to_px = np.empty([50, 1])

        self.lambdaL = lambdaL

        self.lineout = np.zeros([2048, 2048])

        self.lineout_x = np.empty([2048, 1])

        self.FWHM_for_N = np.zeros([2048, 3])

        self.C = 17.5 / 2048

        self.count_N = 30

        self.N_list = np.zeros([self.count_N, 1])

        self.N_diffraction_limit = np.zeros([self.count_N, 1])

        self.filedescription = filedescription

    def open_file(self):

        self.picture = plt.imread(self.filename)

        return self.picture

    def background(self):

        back_mean = np.mean(self.picture[:, 1780:1948], axis=1)

        N = len(self.picture) - 1

        for i in range(1, N):
            self.x_backsubstracted[::, i] = self.picture[::, i] - back_mean[i]

            i = i + 1

        plt.figure(3)

        plt.imshow(self.x_backsubstracted)

    def grating_function(self):
        # create x axis in separate list.
        N = 2048

        for i in range(0, N):
            self.x_axis_in_nm[i] = 1.24679344e-06 * i ** 2 - 1.65566701e-02 * i + 5.22598053e+01

            # i = i+1

    def N_borders_in_px(self):

        # maximum and minimum harmonic orders on the picture

        self.N_min = self.lambdaL / self.x_axis_in_nm[0] + 1

        self.N_max = self.lambdaL / self.x_axis_in_nm[-1] - 6

        self.N_count = int(self.N_max - self.N_min)

        print(int(self.N_min), int(self.N_max), "Nmin, Nmax on this shot")

        # borders to px

        for i in range(self.N_min, self.N_max):
            x = self.lambdaL / i
            # print(x, "N", i)
            self.x_axis_nm_to_px[i] = np.rint(4.71439193e-01 * x ** 2 - 1.06651902e+02 * x + 4.29603367e+03)

    def N_HHG_ROI_horizontal(self):

        A = int(self.N_min)

        N = int(self.N_max)

        # visualisation of taken lineouts on picture
        # integration over a certain amount of lines self.roi_y
        # creates lineout (sum of lineouts per N)
        for i in range(A, N):

            x = 0

            a = int(self.x_axis_nm_to_px[i])

            self.FWHM_for_N[a + x, 0] = i

            self.FWHM_for_N[a + x, 1] = a + x

            plt.figure(3)

            plt.hlines(a + self.ROI_y, xmin=0, xmax=2048, color="m", linewidth=0.5)

            plt.hlines(self.x_axis_nm_to_px[i], xmin=0, xmax=2048, color="w", linewidth=0.5)

            for x in range(0, self.ROI_y):
                self.lineout[a, ::] = self.lineout[a, ::] + self.x_backsubstracted[a + x, ::]

        # deletes 0 entrys in this array
        self.FWHM_for_N = np.delete(self.FWHM_for_N, np.where(~ self.FWHM_for_N.any(axis=1))[0], axis=0)
        # print(self.FWHM_for_N, "List_of_line_outs_N")

        plt.savefig(self.filedescription + "_divergence" + ".png", bbox_inches="tight", dpi=1000)

        # print(self.FWHM_for_N, "x-achse- oders: N, px,0")
        return self.lineout, self.FWHM_for_N

    def find_FWHM(self):

        # creates x-axis array

        self.lineout_x = np.arange(self.xmin, self.xmax)

        i = 0

        # determines maximum of each lineout
        # determines FWHM via stepfunction for Maximum/2
        # writes FWHM int the created array to each corresponding HHG number

        for i in range(0, (len(self.FWHM_for_N))):
            j = int(self.FWHM_for_N[i, 1])
            # print (len(self.FWHM_for_N), "laenge i", j, "j")

            plt.figure(1)

            plt.plot(self.lineout_x, self.lineout[j, self.xmin:self.xmax])

            maximum = np.amax(self.lineout[j, self.xmin:self.xmax])
            # offset of background....
            minimum = 1000

            half_max = (maximum - minimum) / 2

            # gives for one peak stepfunction
            # width of step function is FWHM
            d = np.sign(half_max - self.lineout[j, self.xmin:self.xmax]) - 1

            plt.figure(2)

            plt.plot(self.lineout_x, d)

            FWHM = np.amax(np.nonzero(d)) - np.amin(np.nonzero(d))

            self.FWHM_for_N[i, 2] = self.C * (np.amax(np.nonzero(d)) - np.amin(np.nonzero(d)))

        print("output_list: N, px, FWHM", self.FWHM_for_N)

        return self.FWHM_for_N

    def plot_diffraction_limit(self):

        C1 = 1 / (math.pi * self.w0)

        denting_constant = 0

        for x in range(0, self.count_N - 1):
            self.N_list[x] = int(self.N_min) + x

            self.N_diffraction_limit[x] = C1 * (
            ((4 * math.pi * denting_constant) ** 2 + (self.lambdaL * 1E-9 / self.N_list[x]) ** 2)) ** 0.5

        plt.figure(5)

        plt.scatter(self.N_list, self.N_diffraction_limit, marker="<", color="c", label="theory D0 = 0 nm", alpha=0.5)

        plt.xlim(15, 28)

        plt.ylim(1, 11)

        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        plt.savefig("20190123_divergence_mrad" + ".png", bbox_inches="tight", dpi=1000)

    def plot_result(self):

        marker2 = random.choice(marker_style)

        marker_style.remove(marker2)

        plt.figure(5)

        plt.ylabel("FWHM [mrad]")

        plt.xlabel("harmonic order [N]")

        plt.scatter(self.FWHM_for_N[::, 0], self.FWHM_for_N[::, 2], label=self.filedescription, linewidth=0.2,
                    marker=marker2)

        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    def save_data(self):
        np.savetxt(self.filedescription + "_divergence_mrad" + ".txt", self.FWHM_for_N, delimiter=' ', fmt='%1.4e')


# class gives the following params: ("filepath/filename',xmin [px], xmax [px],fundamental wavelength[um],ROI_y, filename )


Picture1 = Open_and_Plot_Picture('rotated/spectro1__Wed Jan 23 2019_16.50.14_62.tif', 150, 1500, 807.5, 40,
                                 "20190123_62 z2400 GVD900")
Picture1.open_file()
Picture1.background()
Picture1.grating_function()
Picture1.N_borders_in_px()
Picture1.N_HHG_ROI_horizontal()
Picture1.find_FWHM()
Picture1.plot_result()

Picture1 = Open_and_Plot_Picture('rotated/spectro1__Wed Jan 23 2019_14.51.31_17.tif', 150, 1500, 813.5, 40,
                                 "20190123_17 z2000 GVD900")
Picture1.open_file()
Picture1.background()
Picture1.grating_function()
Picture1.N_borders_in_px()
Picture1.N_HHG_ROI_horizontal()
Picture1.find_FWHM()
Picture1.plot_result()

Picture1 = Open_and_Plot_Picture('rotated/spectro1__Wed Jan 23 2019_16.15.56_46.tif', 150, 1500, 807.5, 50,
                                 "20190123_46 z3600 GVD900")
Picture1.open_file()
Picture1.background()
Picture1.grating_function()
Picture1.N_borders_in_px()
Picture1.N_HHG_ROI_horizontal()
Picture1.find_FWHM()
Picture1.plot_result()

Picture1 = Open_and_Plot_Picture('rotated/spectro1__Wed Jan 23 2019_16.13.00_44.tif', 150, 1500, 807.5, 50,
                                 "20190123_44 z3200 GVD900")
Picture1.open_file()
Picture1.background()
Picture1.grating_function()
Picture1.N_borders_in_px()
Picture1.N_HHG_ROI_horizontal()
Picture1.find_FWHM()
Picture1.plot_result()

Picture1 = Open_and_Plot_Picture('rotated/spectro1__Wed Jan 23 2019_16.51.09_63.tif', 100, 1400, 808.5, 40,
                                 "20190123_63 z2000 GVD900")
Picture1.open_file()
Picture1.background()
Picture1.grating_function()
Picture1.N_borders_in_px()
Picture1.N_HHG_ROI_horizontal()
Picture1.find_FWHM()
Picture1.plot_result()

Picture1 = Open_and_Plot_Picture('rotated/spectro1__Wed Jan 23 2019_14.46.17_14.tif', 200, 1400, 812.5, 50,
                                 "20190123_14 z1600 GVD900")
Picture1.open_file()
Picture1.background()
Picture1.grating_function()
Picture1.N_borders_in_px()
Picture1.N_HHG_ROI_horizontal()
Picture1.find_FWHM()
Picture1.plot_result()
Picture1.plot_diffraction_limit()

Picture1.save_data()

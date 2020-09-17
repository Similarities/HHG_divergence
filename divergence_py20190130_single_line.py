
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 16:54:54 2019

@author: similarities
"""

import matplotlib.pyplot as plt
import numpy as np
import random

marker_style = [',', '+', '.', 'v', '*', "P", "H", 'x', "D", ">"]


class FwhmImageProcessing:

    def __init__(self, filename, x_min, x_max, lambda_fundamental, px_range, harmonic_number, file_discription):
        self.filename = filename
        # px size full picture * usually 0 - 2048
        self.y_min = 0
        self.y_max = 2048
        # defined Roi in x to avoid boarder effects
        self.x_min = x_min
        self.x_max = x_max
        # integration ROI y for each HHG line
        self.pixel_range = px_range
        self.picture = np.empty([])
        self.integrated = np.empty([])
        self.x_backsubstracted = np.empty([2048, 2048])
        self.lambda_fundamental = lambda_fundamental
        self.line_out = np.zeros([self.y_max, 1])
        self.line_out_x = np.arange(self.x_min, self.x_max)
        self.fwhm_for_harmonic = np.zeros([2048, 3])
        self.calibration_to_msr = 17.5 / 2048
        self.full_divergence = 17.5
        self.harmonic_counter = int
        self.integrated_signal = np.zeros([20, 4])  # '?'
        self.normalization_factor_mrad = np.zeros([20, 1])
        self.file_description = file_discription
        self.harmonic_selected = harmonic_number
        self.result_array = np.zeros([1, 4])
        self.border_up = self.select_harmonic_in_px()


    def open_file(self):
        self.picture = plt.imread(self.filename)
        return self.picture

    def background(self):
        back_mean = np.mean(self.picture[:, 1780:1948], axis=1)
        for x in range(0, self.y_max):
            self.x_backsubstracted[::, x] = self.picture[::, x] - back_mean[x]
        plt.figure(1)
        plt.ylim(100, 1000)
        plt.imshow(self.x_backsubstracted)

    def grating_function(self):
        # for testing, inverse function of nm-to px
        x_axis_in_nm = np.empty([2048, 1])
        for x in range(0, self.y_max):
            x_axis_in_nm[x] = 1.27877896e-06 * x ** 2 - 1.37081526e-02 * x + 3.46785380e+01
        return x_axis_in_nm

    def select_harmonic_in_px(self):
        # visualisation of taken lineouts on picture
        # borders in nm
        harmonic_in_nm = self.lambda_fundamental / self.harmonic_selected
        # print(self.lambda_fundamental, self.harmonic_selected, "fundamental wavelength, selected harmonic number N")
        a = int(7.79104482e-01 * harmonic_in_nm ** 2 - 1.24499534e+02 * harmonic_in_nm + 3.38549944e+03)
        self.border_up = int(a - self.pixel_range / 2)
        self.plot_roi_on_image(0, 2048)
        print(self.border_up, self.border_up+self.pixel_range, "ROI in px")
        self.initialize_result_array(a)
        return self.border_up

    def plot_roi_on_image(self, xmin, xmax):
        plt.figure(1)
        plt.hlines(self.border_up, xmin=xmin, xmax=xmax, color="m", linewidth=0.5)
        plt.hlines(self.border_up + self.pixel_range, xmin=xmin, xmax=xmax, color="w", linewidth=1.)
        plt.show()

    def initialize_result_array(self, a):
        print('xxxxxxxxxxxxx', self.result_array)
        self.result_array[0, 3] = self.lambda_fundamental / self.harmonic_selected
        self.result_array[0, 0] = self.harmonic_selected
        self.result_array[0, 1] = a
        return self.result_array

    def sum_over_pixel_range_y(self):
        self.line_out = self.x_backsubstracted[self.border_up: self.border_up + self.pixel_range, ::]
        self.line_out = np.sum(self.line_out, axis=0)
        self.line_out = self.line_out[self.x_min:self.x_max]
        return self.line_out

    def correction_background(self, value):
        self.line_out[::] = self.line_out[::] - value
        return self.line_out

    def integrated_signal_in_lineout(self):
        integrated = np.sum(self.line_out, axis=0)
        print(integrated, 'counts in ROI')
        return integrated

    def plot_x_y(self, x, y, name, plot_number, axis_x_name, axis_y_name):
        plt.figure(plot_number)
        plt.plot(x, y, label=name)
        plt.xlabel(str(axis_x_name))
        plt.ylabel(str(axis_y_name))

    def calibrate_px_to_msr(self):
        self.line_out_x[::] = self.line_out_x[::] * self.calibration_to_msr
        return self.line_out_x

    def step_function_for_fwhm(self):
        self.sum_over_pixel_range_y()
        self.plot_x_y(self.line_out_x, self.line_out, 'line-out un-corrected', 2, 'px', 'counts')
        maximum = np.amax(self.line_out[::])
        minimum = np.amin(self.line_out[::])
        half_max = (maximum - minimum) / 2
        self.correction_background(minimum)
        self.plot_x_y(self.line_out_x, self.line_out, 'line-out corrected', 2, 'px', 'counts')
        # width of step function is FWHM
        d = np.sign(half_max - self.line_out[::]) - 1
        self.calibrate_px_to_msr()
        self.plot_x_y(self.line_out_x, d, 'stepfunction', 3, 'mrad', 'value')
        self.result_array[0, 2] = (np.amax(np.nonzero(d)) - np.amin(np.nonzero(d)))
        print(self.result_array)
        return self.result_array


# insert the following ('filepath/picture_name.tif', border in picture L(int), border in picture R (int), fundamental frequency (float), ROI_y(px), harmonic number (int), "picture name for plot")
Picture1 = FwhmImageProcessing('rotated\spectro1__Wed Jan 30 2019_11.13.52_14.tif', 100, 1200, 800., 50,
                               25, "20190123_xx")
Picture1.open_file()
Picture1.background()

Picture1.select_harmonic_in_px()
Picture1.step_function_for_fwhm()
Picture1.integrated_signal_in_lineout()
plt.show()


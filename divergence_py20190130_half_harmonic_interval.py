# !/usr/bin/env python2
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
        self.y_min = 0
        self.y_max = 2048
        self.x_min = x_min
        self.x_max = x_max
        self.pixel_range = px_range
        self.picture = np.empty([])
        # self.integrated = np.empty([])
        self.x_backsubstracted = np.empty([2048, 2048])
        self.lambda_fundamental = lambda_fundamental
        self.line_out = np.zeros([self.y_max, 1])
        self.line_out_x = np.arange(self.x_min, self.x_max)
        self.fwhm_for_harmonic = np.zeros([2048, 3])
        self.calibration_to_msr = 17.5 / 2048
        self.full_divergence = 17.5
        self.normalization_factor_mrad = np.zeros([20, 1])
        self.file_description = self.filename
        self.harmonic_selected = harmonic_number
        self.result_array = np.zeros([1, 5])
        self.border_up, self.border_down = self.energy_range()

    def open_file(self):
        self.picture = plt.imread(self.filename)
        return self.picture

    def background(self):
        back_mean = np.mean(self.picture[:, 1780:1948], axis=1)
        for x in range(0, self.y_max):
            self.x_backsubstracted[::, x] = self.picture[::, x] - back_mean[x]
        plt.figure(1)
        plt.imshow(self.x_backsubstracted)

    def grating_function(self):
        # for testing, inverse function of nm-to px
        x_axis_in_nm = np.empty([2048, 1])
        for x in range(0, self.y_max):
            x_axis_in_nm[x] = 1.27877896e-06 * x ** 2 - 1.37081526e-02 * x + 3.46785380e+01
        return x_axis_in_nm

    def nm_in_px(self, energy_nm):
        a = int(7.79104482e-01 * energy_nm ** 2 - 1.24499534e+02 * energy_nm + 3.38549944e+03)
        return a

    def energy_range(self):
        previous_harmonic = self.lambda_fundamental / (self.harmonic_selected - 0.5)
        next_harmonic = self.lambda_fundamental / (self.harmonic_selected + 0.5)
        self.border_up = np.int(self.nm_in_px(previous_harmonic))
        self.border_down = np.int(self.nm_in_px(next_harmonic))
        print(self.border_up, self.border_down, "ROI in px")
        self.pixel_range = np.int(self.border_down - self.border_up)
        print(self.pixel_range, 'ROI in pixel range')
        self.plot_roi_on_image(0, 2048)

        return self.border_up, self.border_down

    def plot_roi_on_image(self, xmin, xmax):
        plt.figure(1)
        plt.hlines(self.border_up, xmin=xmin, xmax=xmax, color="y", linewidth=1)
        plt.hlines(self.border_down, xmin=xmin, xmax=xmax, color="w", linewidth=1.)
        harmonic_in_nm = self.lambda_fundamental / self.harmonic_selected
        plt.hlines(np.int(self.nm_in_px(harmonic_in_nm)), xmin=xmin, xmax=xmax,
                  color='r', linewidth=0.5)

    def initialize_result_array(self):
        self.result_array[0, 3] = self.lambda_fundamental / self.harmonic_selected
        self.result_array[0, 0] = self.harmonic_selected
        return self.result_array

    def sum_over_pixel_range_y(self):
        self.line_out = self.x_backsubstracted[self.border_up: self.border_down, ::]
        self.line_out = np.sum(self.line_out, axis=0)
        self.line_out = self.line_out[self.x_min:self.x_max]
        return self.line_out

    def correction_background(self, value):
        self.line_out[::] = self.line_out[::] - value
        return self.line_out

    def integrated_signal_in_lineout(self):
        integrated = np.sum(self.line_out[::])
        self.result_array[0, 2] = integrated
        return integrated, self.result_array

    def plot_x_y(self, x, y, name, plot_number, axis_x_name, axis_y_name):
        plt.figure(plot_number)
        plt.plot(x, y, label=name)
        plt.xlabel(str(axis_x_name))
        plt.ylabel(str(axis_y_name))

    def calibrate_px_to_msr(self, array_x):
        array_x[::] = array_x[::] * self.calibration_to_msr
        return array_x

    def prepare_for_stepfunction(self):
        self.sum_over_pixel_range_y()
        self.plot_x_y(self.line_out_x, self.line_out, 'line-out un-corrected', 2, 'px', 'counts')
        maximum = np.amax(self.line_out[::])
        minimum = np.amin(self.line_out[::])
        half_max = (maximum - minimum) / 2
        self.correction_background(minimum)
        self.plot_x_y(self.line_out_x, self.line_out, 'line-out corrected', 2, 'px', 'counts')
        return half_max

    def step_function_for_fwhm(self):
        half_max = self.prepare_for_stepfunction()
        # width of step function is FWHM
        d = np.sign(half_max - self.line_out[::]) - 1
        self.line_out_x = self.calibrate_px_to_msr(self.line_out_x)
        self.plot_x_y(self.line_out_x, d, 'stepfunction', 3, 'mrad', 'value')
        self.result_array[0, 1] = self.calibration_to_msr * (np.amax(np.nonzero(d)) - np.amin(np.nonzero(d)))
        print(self.calibration_to_msr * (np.amax(np.nonzero(d)) - np.amin(np.nonzero(d))), 'FWHM')
        return self.result_array

    def px_in_nm(self, px_number):
        return 1.24679344e-06 * px_number ** 2 - 1.65566701e-02 * px_number + 5.22598053e+01

    def delta_energy(self):
        lower = int(self.border_up)
        upper = int(self.border_down)
        delta = self.px_in_nm(lower) - self.px_in_nm(upper)
        # Delta Energy / harmonic all in nm
        delta_vs_energy = delta / (self.lambda_fundamental / self.harmonic_selected)
        self.result_array[0, 4] = delta_vs_energy
        return delta_vs_energy, delta, self.result_array

    def prepare_header(self):
        self.initialize_result_array()
        self.integrated_signal_in_lineout()
        self.delta_energy()
        # insert header line and change index
        header_names = (['harmonic number', 'mrad', 'integrated counts in delta E', 'harmonic in nm', 'delta E/E'])
        parameter_info = (
            ['fundamental_nm:', str(self.lambda_fundamental), 'pixel_range:', str(self.pixel_range), 'xxxx'])
        return np.vstack((header_names, self.result_array, parameter_info))

    def save_data(self):
        result = self.prepare_header()
        print(self.result_array)
        print('saved data')
        np.savetxt(self.file_description[33:42] + self.file_description[-6:-4] + 'divergence' + ".txt", result,
                   delimiter=' ',
                   header='string', comments='', fmt='%s')


# insert the following ('filepath/picture_name.tif', border in picture L(int), border in picture R (int), fundamental frequency (float), ROI_y(px), harmonic number (int), "picture name for plot")
Picture1 = FwhmImageProcessing('rotated_20190130_1\spectro1__Wed Jan 30 2019_10.27.43_1.tif', 100, 1200, 808., 50,
                               25, "20190123_xx")
Picture1.open_file()
Picture1.background()
Picture1.step_function_for_fwhm()
Picture1.save_data()
plt.show()

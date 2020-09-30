# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 16:54:54 2019
@author: similarities
"""

import matplotlib.pyplot as plt
import numpy as np




class FwhmImageProcessing:

    def __init__(self, filename, x_min, x_max, lambda_fundamental, px_range, harmonic_number, file_discription):
        self.filename = filename
        self.filedescription = self.filename[31:42] + '_' + self.filename[-6:-4]
        self.y_min = 0
        self.y_max = 2048
        self.x_min = 150
        self.x_max = 900
        self.pixel_range = px_range
        self.picture = np.empty([])
        self.harmonic_selected = harmonic_number
        # self.integrated = np.empty([])
        self.x_backsubstracted = np.empty([2048, 2048])
        self.lambda_fundamental = lambda_fundamental
        self.line_out = np.zeros([self.y_max, 1])
        self.line_out_x = np.arange(self.x_min, self.x_max)
        #self.fwhm_for_harmonic = np.zeros([2048, 3])
        self.calibration_to_msr = 17.5 / 2048
        self.full_divergence = 17.5
        self.normalization_factor_mrad = np.zeros([20, 1])
        self.border_up = self.select_harmonic_in_px()
        self.maximum_harmonic = 28
        self.result_array = np.zeros([self.maximum_harmonic, 5])


    def open_file(self):
        self.picture = plt.imread(self.filename)
        return self.picture

    def background(self):
        back_mean = np.mean(self.picture[:, 1000:1200], axis=1)
        for x in range(0, self.y_max):
            self.x_backsubstracted[::, x] = self.picture[::, x] - back_mean[x]
        plt.figure(1)
        #plt.ylim(100, 1000)
        plt.imshow(self.x_backsubstracted)

    def select_harmonic_in_px(self):
        harmonic_in_nm = self.lambda_fundamental / self.harmonic_selected
        self.border_up = int(self.nm_in_px(harmonic_in_nm) - self.pixel_range / 2)
        self.plot_roi_on_image(0, 2048)
        print(self.border_up, self.border_up + self.pixel_range, "ROI in px")
        return self.border_up

    def nm_in_px(self, px_in):
        return int(7.79104482e-01 * px_in ** 2 - 1.24499534e+02 * px_in + 3.38549944e+03)

    def plot_roi_on_image(self, xmin, xmax):
        plt.figure(1)
        plt.hlines(self.border_up, xmin=xmin, xmax=xmax, color="m", linewidth=0.5)
        plt.hlines(self.border_up + self.pixel_range, xmin=xmin, xmax=xmax, color="w", linewidth=1.)


    def sum_over_pixel_range_y(self):
        self.line_out = self.x_backsubstracted[self.border_up: self.border_up + self.pixel_range, ::]
        self.line_out = np.sum(self.line_out, axis=0)
        self.line_out = self.line_out[self.x_min:self.x_max]
        return self.line_out

    def correction_background(self, value):
        self.line_out[::] = self.line_out[::] - value
        return self.line_out

    def integrated_signal_in_lineout(self):
        integrated = np.sum(self.line_out[::])
        return integrated

    def plot_x_y(self, x, y, name, plot_number, axis_x_name, axis_y_name):
        plt.figure(plot_number)
        plt.plot(x, y, label=name)
        plt.xlabel(str(axis_x_name))
        plt.ylabel(str(axis_y_name))
        plt.legend()

    def calibrate_px_to_msr(self, array_x):
        array_x[::] = array_x[::] * self.calibration_to_msr
        return array_x

    def prepare_for_stepfunction(self):
        self.sum_over_pixel_range_y()
        maximum = np.amax(self.line_out[::])
        minimum = np.amin(self.line_out[::])
        half_max = (maximum - minimum) / 2
        self.correction_background(minimum)
        #self.plot_x_y(self.line_out_x, self.line_out, 'linout_corrected', 2, 'px', 'counts')
        self.plot_x_y(self.line_out_x, self.line_out, str(self.harmonic_selected), 2, 'px', 'counts')
        return half_max

    def step_function_for_fwhm(self):
        half_max = self.prepare_for_stepfunction()
        # width of step function is FWHM
        d = np.sign(half_max - self.line_out[::]) - 1
        self.line_out_x = self.calibrate_px_to_msr(self.line_out_x)
        self.plot_x_y(self.line_out_x, d, 'stepfunction', 3, 'mrad', 'value')
        self.line_out_x = np.arange(self.x_min, self.x_max)
        result_FWHM = self.calibration_to_msr * (np.amax(np.nonzero(d)) - np.amin(np.nonzero(d)))
        return result_FWHM

    def px_in_nm(self, px_number):
        return 1.24679344e-06 * px_number ** 2 - 1.65566701e-02 * px_number + 5.22598053e+01

    def delta_energy(self):
        lower = int(self.border_up + self.pixel_range / 2)
        upper = int(self.border_up + self.pixel_range)
        delta = self.px_in_nm(lower) - self.px_in_nm(upper)
        energy_nm = (self.lambda_fundamental/self.harmonic_selected)
        delta_vs_energy = delta / energy_nm
        return energy_nm, delta_vs_energy

    def batch_over_N(self):
        for x in range(self.harmonic_selected, self.maximum_harmonic):
            self.result_array[x, 0] = x
            self.harmonic_selected = x
            self.select_harmonic_in_px()
            self.result_array[x, 1] = self.step_function_for_fwhm()
            self.result_array[x, 2] = np.sum(self.line_out[::])
            self.result_array[x, 4], self.result_array[x, 3] = self.delta_energy()
            # clean for empty entries
        self.result_array = np.delete(self.result_array, np.where(~self.result_array.any(axis=1))[0],
                                             axis=0)
        print(self.result_array)
        self.plot_scatter(self.result_array[::, 0], self.result_array[::, 1], self.filedescription,
                              'harmonic number N', 'divergence in mrad', 5)
        self.save_data()
        return self.result_array

    def plot_scatter(self, x,y, name, axis_name_x, axis_name_y, plot_number):
        plt.figure(plot_number)
        plt.scatter(x,y, label = name)
        plt.xlabel(axis_name_x)
        plt.ylabel(axis_name_y)

    def prepare_header(self):
        self.integrated_signal_in_lineout()
        self.delta_energy()
        # insert header line and change index
        header_names = (['harmonic number', 'mrad', 'integrated counts in delta E', 'harmonic in nm', 'delta E/E'])
        parameter_info = (
            ['fundamental_nm:', str(self.lambda_fundamental), 'pixel_range:', str(self.pixel_range), 'xxxx'])
        return np.vstack((header_names, self.result_array, parameter_info))

    def save_data(self):
        result = self.prepare_header()
        #print(self.result_array)
        print('saved data')
        np.savetxt(self.filedescription[31:42] + '_' + self.filedescription[-6:-4] + ".txt", self.result_array, delimiter=' ',
                   header='string', comments='',
                   fmt='%s')


# insert the following ('filepath/picture_name.tif', border in picture L(int), border in picture R (int), fundamental frequency (float), ROI_y(px), harmonic number (int), "picture name for plot")
Picture1 = FwhmImageProcessing('rotated_20190130\spectro1__Wed Jan 30 2019_11.49.36_23.tif', 100, 1200, 800., 10,
                               25, "20190123_xx")
Picture1.open_file()
Picture1.background()
Picture1.select_harmonic_in_px()
Picture1.batch_over_N()
Picture1.save_data()
plt.show()
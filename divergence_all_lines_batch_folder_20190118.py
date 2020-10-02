# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 16:54:54 2019
@author: similarities
"""

import matplotlib.pyplot as plt
import numpy as np
import os


class FwhmImageProcessing:

    def __init__(self, filename, lambda_fundamental, maximum_harmonic, harmonic_number):
        self.filename = filename
        self.filedescription = self.filename[31:42] + '_' + self.filename[-6:-4]
        self.y_min = 0
        self.y_max = 2048
        self.x_min = 150
        self.x_max = 1500
        self.picture = np.empty([])
        self.harmonic_selected = harmonic_number
        self.x_backsubstracted = np.empty([2048, 2048])
        self.lambda_fundamental = lambda_fundamental
        self.line_out = np.zeros([self.y_max, 1])
        self.line_out_x = np.arange(self.x_min, self.x_max)
        self.calibration_to_msr = 17.5 / 2048
        self.full_divergence = 17.5
        self.normalization_factor_mrad = np.zeros([20, 1])
        self.border_up, self.border_down = self.energy_range()
        self.maximum_harmonic = maximum_harmonic
        self.result_array = np.zeros([self.maximum_harmonic, 5])

    def open_file(self):
        self.picture = plt.imread(self.filename)
        return self.picture

    def background_y(self):
        back_mean = np.mean(self.picture[:, 1600:1700], axis=1)
        for x in range(0, self.y_max):
            self.x_backsubstracted[::, x] = self.picture[::, x] - back_mean[x]
        self.post_back_ground()
        plt.figure(1)
        # plt.ylim(100, 1000)
        plt.imshow(self.x_backsubstracted)
        plt.vlines(self.x_min, 0, 2048)
        plt.vlines(self.x_max, 0, 2048)
        return self.x_backsubstracted

    def background_x(self):
        back_mean = np.mean(self.x_backsubstracted[1880:1948, :], axis=0)
        for x in range(0, 2048):
            self.x_backsubstracted[x, ::] = self.picture[x, ::] - back_mean[x]
        return self.x_backsubstracted

    def post_back_ground(self):
        back_post = np.mean(self.x_backsubstracted[:, self.x_max + 100: self.x_max + 250], axis=1)
        for counter in range(0, 2048):
            self.x_backsubstracted[::, counter] = self.x_backsubstracted[::, counter] - back_post[counter]
        return self.x_backsubstracted

    def energy_range(self):
        # print(self.harmonic_selected, ':')
        previous_harmonic = self.lambda_fundamental / (self.harmonic_selected - 0.3)
        next_harmonic = self.lambda_fundamental / (self.harmonic_selected + 0.3)
        self.border_up = np.int(self.nm_in_px(previous_harmonic))
        self.border_down = np.int(self.nm_in_px(next_harmonic))
        # print(self.border_up, self.border_down, "ROI in px")
        self.pixel_range = np.int(self.border_down - self.border_up)
        # print(self.pixel_range, 'ROI in pixel range')
        self.plot_roi_on_image()
        return self.border_up, self.border_down

    def nm_in_px(self, wavelength_in):
        return int(3.71756520e-01 * wavelength_in ** 2 - 9.89219185e+01 * wavelength_in + 4.50769013e+03)

    def plot_roi_on_image(self):
        plt.figure(1)
        plt.hlines(self.border_up, xmin=0, xmax=2048, color="w", linewidth=0.1)
        plt.hlines(self.border_down, xmin=0, xmax=2048, color="g", linewidth=0.1)

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

        if minimum < 0:
            self.correction_background(minimum)
            maximum = np.amax(self.line_out[::])
            minimum = np.amin(self.line_out[::])

        half_max = (maximum - minimum) / 2
        # self.plot_x_y(self.line_out_x, self.line_out, 'linout_corrected', 2, 'px', 'counts')
        self.plot_x_y(self.line_out_x, self.line_out, str(self.harmonic_selected), 2, 'px', 'counts')
        return half_max

    def step_function_for_fwhm(self):
        half_max = self.prepare_for_stepfunction()
        # width of step function is FWHM
        d = np.sign(half_max - self.line_out[::]) - 1
        self.line_out_x = self.calibrate_px_to_msr(self.line_out_x)
        self.plot_x_y(self.line_out_x, d, 'stepfunction', 3, 'mrad', 'value')
        self.line_out_x = np.arange(self.x_min, self.x_max)
        result_FWHM = 1. * self.calibration_to_msr * (np.amax(np.nonzero(d)) - np.amin(np.nonzero(d)))
        return result_FWHM

    def px_in_nm(self, px_number):
        return 1.22447518e-06 * px_number ** 2 - 1.73729829e-02 * px_number + 5.82820234e+01

    def delta_energy(self):
        delta = self.px_in_nm(self.border_up) - self.px_in_nm(self.border_down)
        energy_nm = (self.lambda_fundamental / self.harmonic_selected)
        delta_vs_energy = delta / energy_nm
        return energy_nm, delta_vs_energy

    def batch_over_N(self):
        for x in range(self.harmonic_selected, self.maximum_harmonic):
            self.result_array[x, 0] = x
            self.harmonic_selected = x
            self.energy_range()
            self.result_array[x, 1] = self.step_function_for_fwhm()
            self.result_array[x, 2] = np.sum(self.line_out[::])
            self.result_array[x, 4], self.result_array[x, 3] = self.delta_energy()
            # clean for empty entries
        self.result_array = np.delete(self.result_array, np.where(~self.result_array.any(axis=1))[0],
                                      axis=0)
        self.plot_scatter(self.result_array[::, 0], self.result_array[::, 1], self.filedescription,
                          'harmonic number N', 'divergence in mrad', 5)
        self.save_data()
        return self.result_array

    def plot_scatter(self, x, y, name, axis_name_x, axis_name_y, plot_number):
        plt.figure(plot_number)
        plt.scatter(x, y, label=name)
        plt.xlabel(axis_name_x)
        plt.ylabel(axis_name_y)
        # plt.legend()

    def prepare_header(self):
        self.integrated_signal_in_lineout()
        self.delta_energy()
        # insert header line and change index
        header_names = (['harmonic_number', 'mrad', 'integrated_counts_in_delta_E', 'harmonic_in_nm', 'delta_E/E'])
        parameter_info = (
            ['fundamental_nm:', str(self.lambda_fundamental), 'pixel_range:', str(self.border_down - self.border_up),
             'x_range: ' + str(self.x_min) + ':' + str(self.x_max)])
        return np.vstack((header_names, self.result_array, parameter_info))

    def save_data(self):
        result = self.prepare_header()
        print('....saving:', self.filedescription)
        plt.figure(1)
        plt.savefig(self.filedescription + "_raw_roi_" + ".png", bbox_inches="tight", dpi=1000)
        plt.figure(2)
        plt.savefig(self.filedescription + "_integrated_lineout" + ".png", bbox_inches="tight", dpi=1000)
        plt.figure(5)
        plt.savefig(self.filedescription + "_div_mrad_FWHM" + ".png", bbox_inches="tight", dpi=1000)
        print('saved data')
        np.savetxt(self.filedescription + ".txt", result, delimiter=' ',
                   header='string', comments='',
                   fmt='%s')


def get_file_list(path_picture):
    tif_files = []
    counter = 0
    for file in os.listdir(path_picture):
        print(file)
        try:
            if file.endswith(".tif"):
                tif_files.append(str(file))
                counter = counter + 1
            else:
                print("only other files found")
        except Exception as e:
            raise e
            print("no files found here")
    return tif_files


def process_files(my_files, path):
    for x in range(0, len(my_files)):
        file = path + '/' + my_files[x]
        Processing_Picture = FwhmImageProcessing(file, 800, 30, 17)
        Processing_Picture.open_file()
        Processing_Picture.background_y()
        Processing_Picture.batch_over_N()
        Processing_Picture.save_data()
        plt.close(1)
        plt.close(2)
        plt.close(5)


my_files = get_file_list('rotated_20190118')
process_files(my_files, 'rotated_20190118')

# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 13:33:12 2019

@author: similarity
"""

import numpy as np
import matplotlib.pyplot as plt
import math

class full_HHG_beamsize:

    def __init__(self, Nmax, beam_diameter, focal_lenght, distance_source_grating, capturing_r_V, capturing_r_H, D0, w0):

        self.Nmax = Nmax

        self.beam_diameter = beam_diameter

        self.focal_length = focal_lenght

        self.distance_source_grating = distance_source_grating

        self.capturing_r_V = capturing_r_V

        self.capturing_r_H = capturing_r_H

        self.D0 = D0

        self.w0 = w0 *1E-6

        self.lambdaL = 800 *1E-9

        self.beam_radius_at_grating = (self.beam_diameter / self.focal_length) *self.distance_source_grating
        #print(self.beam_radius_at_grating, "radius at grating")

        self.switch = 0

        self.a = [i for i in range (1, Nmax +1 )]

        self.b = [i for i in range (1, Nmax +1 )]

        self.c = [i for i in range (1, Nmax +1 )]






    def divergence_HHG_with_denting(self):


        for x in range(0, self.Nmax):

            C1 = 1/(math.pi *self.w0)

            self.b[x] = C1 *( ( (4 *math.pi *self.denting_constant()) **2 + (self.lambdaL / self.a[x] ) **2) ) **0.5

            self.b[x] = self.b[x] * 1000


        print( "denting for the given parameters in [m]:", self.denting_constant() )

        print("w0: [um] ", self.w0)

        print("denting D0: [nm] ", self.D0)


        #self.plot_results(self.a, self.b, "harmonic order N", "HHG divergence in [mrad]", str(self.D0)+"nm ")

        return self.b



    def denting_constant(self):

        xi = self.D0*1E-9

        #print("maximum denting", xi)

        return xi


    def theoretical_HHG_beam_radius_at_distance(self):


        self.divergence_HHG_with_denting()

        for i in range(1, self.Nmax):

            self.b[i] = 0.5 * self.distance_source_grating * self.b[i] / 1000

            #print( self.b[i], "in mm radius", i, "harmonic N", self.beam_radius_at_grating/i, "compare simple equation")

        return self.b




    def test_switch(self, var_detector):

        for i in range(0, self.Nmax):

            a = self.b[i]
            #print(var_detector, str(var_detector), a, "theoretical")


            if a > var_detector:

                print("capturing angle smaller than HHG beam divergence, for N :", self.a[i])

                print("theoretical: ", a,  " captured: " , var_detector )

                switch_ = self.a[i]


            else:

                switch_ = self.a[i]

                return switch_












        #print(switch_, "switch")
        return switch_



        #print(self.switch, "from here the capturing angle is bigger as HHG radius")

        #return self.switch







    def calculate_relative_area_detected(self):

        switch_V_ = self.test_switch( self.capturing_r_V)

        switch_H_ = self.test_switch(self.capturing_r_H)

        print( "div smaller detector V for N:", switch_V_, self.capturing_r_V)

        print( "div smaller detector H for N:", switch_H_, self.capturing_r_H)






        for x in range(switch_H_, switch_V_):

            captured_area = 3.142 * (self.capturing_r_V * self.b[x])

            print(captured_area, "captured Area for N :", x)

            theoretical_Area = 3.142 * ( self.b[x] ) ** 2

            print(theoretical_Area, "theorectial area at grating for N: ", x+1 )

            relative_factor = captured_area / theoretical_Area

            print('relative captured Area for N: ', relative_factor)

           # print(x, "index")

            self.c[x] = relative_factor



        for x in range(0, switch_H_ ):

            captured_area = 3.142 * self.capturing_r_V * self.capturing_r_H

            theoretical_Area = 3.142 * ( self.b[x] ) ** 2


            relative_factor = captured_area / theoretical_Area

            self.c[x] = relative_factor







        if switch_V_ == switch_H_:

            for x in range(switch_H_-1, self.Nmax):

                captured_area = 3.142 * (self.capturing_r_V * self.capturing_r_H)

                theoretical_Area = 3.142 * ( self.b[x] ) ** 2

                print(theoretical_Area, self.a[x], captured_area)

                relative_factor = captured_area / theoretical_Area

                self.c[x] = relative_factor

                print("set self.c", self.c[x])



        elif switch_V_ < self.Nmax:

            for x in range(switch_V_, self.Nmax):

                relative_factor = 1

                self.c[x] = relative_factor
                


        self.c = self.reciprocal_(self.c)

        self.plot_results(self.a[10 ::], self.c[10 ::], "N", "theoretical/ captured Area", "Denting "+str(self.D0)+' [nm]')

        return self.a, self.c




    def reciprocal_(self, array1D_):


        for x in range(0, len(array1D_)):

            #print(len(array1D))



            if array1D_[x] == 0:
                #print(x)

                BaseException



            else:


                array1D_[x] = 1/ array1D_[x]



        return array1D_




    def plot_results(self, x_, y_,  xlabel_, ylabel_, legend_):


        plt.scatter(x_, y_, label = legend_)

        plt.xlabel(xlabel_)

        plt.ylabel(ylabel_)

        plt.legend()




#full_HHG_beamsize(self, Nmax, beam_diameter, focal_lenght, distance_source_grating, capturing_r_V, capturing_r_H, Denting[nm], beam waist [um])
# all in [mm]

Full_beam0 = full_HHG_beamsize(40, 60, 1500, 1416, (3.66 *0.5), (25 *0.5), 0, 6)


Full_beam0.divergence_HHG_with_denting()
Full_beam0.theoretical_HHG_beam_radius_at_distance()
Full_beam0.calculate_relative_area_detected()


Full_beam5 = full_HHG_beamsize(40, 60, 1500, 1416, (3.66 *0.5), (25 *0.5), 5, 6)
Full_beam5.divergence_HHG_with_denting()
Full_beam5.theoretical_HHG_beam_radius_at_distance()
Full_beam5.calculate_relative_area_detected()

Full_beam20 = full_HHG_beamsize(40, 60, 1500, 1416, (3.66 *0.5), (25 *0.5), 20, 6)
Full_beam20.divergence_HHG_with_denting()
Full_beam20.theoretical_HHG_beam_radius_at_distance()
Full_beam20.calculate_relative_area_detected()


Full_beam40 = full_HHG_beamsize(40, 60, 1500, 1416, (3.66 *0.5), (25 *0.5), 40, 6)
Full_beam40.divergence_HHG_with_denting()
Full_beam40.theoretical_HHG_beam_radius_at_distance()
Full_beam40.calculate_relative_area_detected()


plt.savefig("denting_div_loss_per_N"+".png",  bbox_inches="tight", dpi = 1000)
plt.show()







    


    

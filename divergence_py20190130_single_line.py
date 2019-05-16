#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 16:54:54 2019

@author: similarities
"""

import matplotlib.pyplot as plt
import numpy as np
import random


marker_style = [',', '+', '.', 'v', '*', "P", "H", 'x', "D",">"]



class Open_and_Plot_Picture:
    
    
        def __init__(self, filename,xmin, xmax,lambdaL,ROI_y,N_select,filedescription):
            self.filename = filename
            # px size full picture * usually 0 - 2048
            self.ymin = 0
            self.ymax = 2048
            # defined Roi in x to avoid boarder effects
            self.xmin = xmin
            self.xmax = xmax
            # integration ROI y for each HHG line
            self.ROI_y = ROI_y
            
            self.picture = np.empty([])
            self.integrated= np.empty([])
            self.x_backsubstracted=np.empty([2048, 2048])
            self.x_axis_in_nm =np.empty([2048,1])
            self.x_axis_nm_to_px = np.empty([50,1])
            self.lambdaL = lambdaL
            self.lineout = np.zeros([2048,2048])
            self.lineout_x = np.empty([2048,1])
            self.FWHM_for_N = np.zeros([2048,3])
            self.C = 17.5/2048
            self.full_divergence = 17.5
            self.N_count = int
            self.N_integrated=np.zeros([20,4])
            
            self.normalization_factor_mrad = np.zeros([20,1])
            
            self.filedescription = filedescription
            self.N_select = N_select
            self.FWHM_for_single_N = np.zeros([1,4])



        def open_file(self):
            
            self.picture = plt.imread(self.filename)
            
            return self.picture
            
      
        
        def background(self):
            
            back_mean=np.mean(self.picture[:, 1780:1948], axis = 1)
            
            
            i=1
            
            N=len(self.picture)-1
            
            while i<= N:
                
                self.x_backsubstracted[::,i] = self.picture[::,i]- back_mean[i]

                i = i+1
                
                
            plt.figure(3)
            plt.ylim(100,1000)
            
            plt.imshow(self.x_backsubstracted)



        
        def grating_function(self):
            #create x axis in separate list.
            N = 2048
            i = 0

            while i <= N-1:
                
                self.x_axis_in_nm[i] = 1.27877896e-06*i**2 -1.37081526e-02*i +  3.46785380e+01

                i = i+1
   
    
    
   
                
                

                

        def N_selected_ROI(self):
        
            
            
            
            # visualisation of taken lineouts on picture


            # borders in nm
            self.FWHM_for_single_N[0,3] = self.lambdaL/ self.N_select
            print(self.FWHM_for_single_N)

            aa = self.lambdaL/ self.N_select
            
            print(self.lambdaL, self.N_select, "fundamental wavelength, selected harmonic number N")
                
            a = int(7.79104482e-01*aa**2 -1.24499534e+02*aa + 3.38549944e+03) 
            a2 = 1.22447518e-06*a**2 -1.73729829e-02*a+  5.82820234e+01
            print(a, "selected N in px")
            print(a2, "selected N in nm backwards")
            border_up = int(a - self.ROI_y/2)
            border_down = int(a + self.ROI_y/2)
            
            print(border_up, border_down, "ROI in px")
            
            
            #self.FWHM_for_single_N[0] = self.N_select
                
            #3self.FWHM_for_single_N[0,1] = a
            
            self.FWHM_for_N[0,0] = self.N_select
            self.FWHM_for_N[0,1] = a
            self.FWHM_for_N= np.delete(self.FWHM_for_N,np.where(~self.FWHM_for_N.any(axis=1))[0], axis=0)
                
            plt.figure(3)
                
                
            plt.hlines(border_up, xmin = 0, xmax = 2048, color ="m", linewidth =0.5)
                
            plt.hlines(border_down, xmin = 0, xmax = 2048, color ="w", linewidth =1.)
        


            for x in range(0,int(self.ROI_y)):
            
                self.lineout[border_up,::] =self.lineout[border_up, ::] + self.x_backsubstracted[border_up+x,::]
                

            return self.lineout, self.FWHM_for_N

    
    
          
      
   

         
            
        def find_FWHM(self):
            
            #creates x-axis array
            
            self.lineout_x = np.arange(self.xmin,self.xmax)
            
            i = 0
            
            # determines maximum of each lineout
            # determines FWHM via stepfunction for Maximum/2
            # writes FWHM int the created array to each corresponding HHG number
            
            for i in range(0,(len(self.FWHM_for_N))):

                j = int(self.FWHM_for_N[i,1]-self.ROI_y/2)
           #print (len(self.FWHM_for_N), "laenge i", j, "j")


                plt.figure(1)
                
                plt.plot(self.lineout_x*self.C, self.lineout[j,self.xmin:self.xmax])
                plt.xlabel("divergence [mrad]")
                plt.show()

                
                maximum = np.amax(self.lineout[j,self.xmin:self.xmax])
                # offset of background....
                print(self.lineout[j,::])
                minimum = self.lineout[j,1100]
                print(minimum, "minimum")

                half_max = (maximum-minimum)/2
                
                    # gives for one peak stepfunction 
                    # width of step function is FWHM
                d = np.sign(half_max - self.lineout[j,self.xmin:self.xmax] )-1
                              
                plt.figure(2)
                
                plt.plot(self.lineout_x, d)
                plt.show()
                
        
                #FWHM = np.amax(np.nonzero(d))-np.amin(np.nonzero(d))

                self.FWHM_for_N[i,2] = self.C*(np.amax(np.nonzero(d))-np.amin(np.nonzero(d)))

            
            print("output_list: N, px, FWHM", self.FWHM_for_N)    
            
            return self.FWHM_for_N
  


        


            
     
            






# insert the following ('filepath/picture_name.tif', border in picture L(int), border in picture R (int), fundamental frequency (float), ROI_y(px), harmonic number (int), "picture name for plot")
#@ mingyuan: "picture name for plot" is unneccessary, try different ROI_y ranges, set boundaryL close to the onset of the detection range, put 
# put the correct fundamental frequency in
Picture1=Open_and_Plot_Picture('second_run/spectro1_secondrun__Wed Jan 30 2019_16.23.08_30.tif',100, 1200, 800., 50, 25,"20190123_12 z1600 GVD900 PPin")
Picture1.open_file()
Picture1.background()

Picture1.N_selected_ROI()
Picture1.find_FWHM()










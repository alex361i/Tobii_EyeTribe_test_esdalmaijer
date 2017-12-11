# -*- coding: utf-8 -*-
import time
import datetime
import tobii_research as tr
import os
import json
from pygaze.display import Display
from pygaze.screen import Screen
from pygaze.eyetracker import EyeTracker
import pygaze.libtime as timer


#disp = Display()
#scr = Screen()
#scr.clear()
#disp.close()

#scr.draw_text("Preparing experiment...", 
#fontsize=20)
#disp.fill(scr)
#disp.show()

#address = "tobii-ttp://is404-100106241134"
#this adress is for the glasses"
address = "tobii-ttp://is404-100106241134"
eyetracker = tr.EyeTracker(address)
print("Model: " + eyetracker.model)

with open("testing here we should have the tobi license file", "rb") as f:
    license = f.read()

failed_licenses_applied_as_list_of_keys = eyetracker.apply_licenses([tr.LicenseKey(license)])
failed_licenses_applied_as_list_of_bytes = eyetracker.apply_licenses([license])
failed_licenses_applied_as_key = eyetracker.apply_licenses(tr.LicenseKey(license))
failed_licenses_applied_as_bytes = eyetracker.apply_licenses(license)

#if len(failed_licenses_applied_as_list_of_keys) == 0:
    #print "Successfully applied license from list of keys."
#else:
    #print "Failed to apply license from list of keys. Validation result: {0}.".\
                                         #format(failed_licenses_applied_as_list_of_keys[0].validation_result)
global_gaze_data = None
fil=open("dumb_data.txt","w")

def gaze_data_callback(gaze_data):
    global global_gaze_data
    global_gaze_data = gaze_data
    #print gaze_data['device_time_stamp']
    
def gaze_data(eyetracker):
    global global_gaze_data
    
    print "Subscribing to gaze data for eye tracker with serial number{0}.".format(eyetracker.serial_number)
    eyetracker.subscribe_to(tr.EYETRACKER_GAZE_DATA, gaze_data_callback,as_dictionary=True)
    
    time.sleep(2)  
    
    eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA, gaze_data_callback)
    print "Unsubscribed from gaze data."
    
    print "Last received gaze package:"
    print global_gaze_data
    print "here it ends"+"\n"
    fil.write(json.dumps(global_gaze_data))
gaze_data(eyetracker)

#t0 = timer.get_time()
#while timer.get_time() - t0 < 5000:
    #gaze_data(eyetracker)
    #scr.clear()
    #scr.draw_fixation(fixtype='dot', pos=gazepos)
    #disp.fill(scr)
    #disp.show()

#eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA, gaze_data_callback)


#gaze_data(eyetracker)   



def notification_callback(notification, data):
    print "Notification {0} received at time stamp {1}.".format(notification,
                        data.system_time_stamp)

def notifications(eyetracker):
    all_notifications =\
        (tr.EYETRACKER_NOTIFICATION_CONNECTION_LOST,
         tr.EYETRACKER_NOTIFICATION_CONNECTION_RESTORED,
         tr.EYETRACKER_NOTIFICATION_CALIBRATION_MODE_ENTERED,
         tr.EYETRACKER_NOTIFICATION_CALIBRATION_MODE_LEFT,
         tr.EYETRACKER_NOTIFICATION_TRACK_BOX_CHANGED,
         tr.EYETRACKER_NOTIFICATION_DISPLAY_AREA_CHANGED,
         tr.EYETRACKER_NOTIFICATION_GAZE_OUTPUT_FREQUENCY_CHANGED)
        
    for notification in all_notifications:
        eyetracker.subscribe_to(notification,
                                lambda x, notification=notification:notification_callback(notification, x))
        print "Subscribed to {0} for eye tracker with serial number {1}.".format(notification,
                             eyetracker.serial_number)
    # Trigger some notifications
    calibration = tr.ScreenBasedCalibration(eyetracker)
    calibration.enter_calibration_mode()
    
    calibration.leave_calibration_mode()
    
    for notification in all_notifications:
        eyetracker.unsubscribe_from(notification)
        print "Unsubscribed from {0}.".format(notification)
        

#print notifications(eyetracker)               


def time_synchronization_data_callback(time_synchronization_data):
    print time_synchronization_data

def time_synchronization_data(eyetracker):
    print "Subscribing to time synchronization data for eye tracker with serial number {0}.".\
                                                                                       format(eyetracker.serial_number)
    eyetracker.subscribe_to(tr.EYETRACKER_TIME_SYNCHRONIZATION_DATA,time_synchronization_data_callback, as_dictionary=True)
    # Wait while some time synchronization data is collected.
    time.sleep(2)
    eyetracker.unsubscribe_from(tr.EYETRACKER_TIME_SYNCHRONIZATION_DATA,time_synchronization_data_callback)
    print "Unsubscribed from time synchronization data."

#time_synchronization_data(eyetracker)                                                                                  

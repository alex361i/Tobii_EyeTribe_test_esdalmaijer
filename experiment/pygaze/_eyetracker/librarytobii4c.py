# -*- coding: utf-8 -*-
#
# This file is part of PyGaze - the open-source toolbox for eye tracking
#
# PyGaze is a Python module for easily creating gaze contingent experiments
# or other software (as well as non-gaze contingent experiments/software)
# Copyright (C) 2012-2013 Edwin S. Dalmaijer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

import psychopy.visual
import psychopy.event
import re

# TobiiTracker
import copy
import math
import numpy 
import pandas 
import Queue
import itertools

from pygaze import settings
from pygaze.libtime import clock
import pygaze
from pygaze.screen import Screen
from pygaze.keyboard import Keyboard
from pygaze.sound import Sound

from pygaze._eyetracker.baseeyetracker import BaseEyeTracker
# we try importing the copy_docstr function, but as we do not really need it
# for a proper functioning of the code, we simply ignore it when it fails to
# be imported correctly
try:
    from pygaze._misc.misc import copy_docstr
except:
    pass


## letobii
#from letobii import TobiiController
#

# TobiiController
from tobii.eye_tracking_io.basic import EyetrackerException

import os
import datetime

import tobii.eye_tracking_io.mainloop
import tobii.eye_tracking_io.browsing
import tobii.eye_tracking_io.eyetracker
#import tobii.eye_tracking_io.time.clock
#import tobii.eye_tracking_io.time.sync

from tobii.eye_tracking_io.types import Point2D, Blob



import tobii_research as tr
# try importing PIL
try:
    from PIL import Image
    from PIL import ImageDraw
except:
    try:
        import Image
        import ImageDraw
    except:
        print("Failed to import PIL.")


# # # # #
# functions

def deg2pix(cmdist, angle, pixpercm):

    """Returns the value in pixels for given values (internal use)

    arguments
    cmdist    -- distance to display in centimeters
    angle        -- size of stimulus in visual angle
    pixpercm    -- amount of pixels per centimeter for display

    returns
    pixelsize    -- stimulus size in pixels (calculation based on size in
               visual angle on display with given properties)
    """

    cmsize = math.tan(math.radians(angle)) * float(cmdist)
    return cmsize * pixpercm


# # # # #
# classes

class TobiiTracker(BaseEyeTracker):

    """A class for Tobii EyeTracker objects"""

    def __init__(self, display, logfile=settings.LOGFILE,
        eventdetection=settings.EVENTDETECTION, saccade_velocity_threshold=35,
        saccade_acceleration_threshold=9500, blink_threshold=settings.BLINKTHRESH, **args):

        """Initializes a TobiiTracker instance

        arguments
        display    --    a pygaze.display.Display instance

        keyword arguments
        None
        """

        # try to copy docstrings (but ignore it if it fails, as we do
        # not need it for actual functioning of the code)
        try:
            copy_docstr(BaseEyeTracker, TobiiTracker)
        except:
            # we're not even going to show a warning, since the copied
            # docstring is useful for code editors; these load the docs
            # in a non-verbose manner, so warning messages would be lost
            pass

        # object properties
        self.disp = display
        self.screen = Screen()
        self.dispsize = settings.DISPSIZE # display size in pixels
        self.screensize = settings.SCREENSIZE # display size in cm
        self.screendist = settings.SCREENDIST # distance between participant and screen in cm
        self.pixpercm = (self.dispsize[0]/float(self.screensize[0]) + self.dispsize[1]/float(self.screensize[1])) / 2.0
        self.kb = Keyboard(keylist=['space', 'escape', 'q'], timeout=1)
        self.errorbeep = Sound(osc='saw',freq=100, length=100)

        # output file properties
        self.outputfile = logfile
        self.description = "experiment" # TODO: EXPERIMENT NAME
        self.participant = "participant" # TODO: PP NAME

        # eye tracker properties
        self.connected = False
        self.recording = False
        self.eye_used = 0 # 0=left, 1=right, 2=binocular
        self.left_eye = 0
        self.right_eye = 1
        self.binocular = 2
        self.errdist = 2 # degrees; maximal error for drift correction
        self.pxerrdist = deg2pix(self.screendist, self.errdist, self.pixpercm)
        self.maxtries = 100 # number of samples obtained before giving up (for obtaining accuracy and tracker distance information, as well as starting or stopping recording)
        self.prevsample = (-1,-1)

        # validation properties
        self.nvalsamples = 2000 # samples for one validation point

        # event detection properties
        self.fixtresh = 1.5 # degrees; maximal distance from fixation start (if gaze wanders beyond this, fixation has stopped)
        self.fixtimetresh = 100 # milliseconds; amount of time gaze has to linger within self.fixtresh to be marked as a fixation
        self.spdtresh = saccade_velocity_threshold # degrees per second; saccade velocity threshold
        self.accthresh = saccade_acceleration_threshold # degrees per second**2; saccade acceleration threshold
        self.blinkthresh = blink_threshold # milliseconds; blink detection threshold used in PyGaze method
        self.eventdetection = eventdetection
        self.set_detection_type(self.eventdetection)
        self.weightdist = 10 # weighted distance, used for determining whether a movement is due to measurement error (1 is ok, higher is more conservative and will result in only larger saccades to be detected)

        # initialize controller
        self.controller = TobiiController(display)
        self.controller.setDataFile("%s_TOBII_output.tsv" % logfile)
        #self.controller.waitForFindEyeTracker()
        self.controller.activate(self.controller.eyetracker)

        # initiation report
        self.controller.datafile.write("pygaze initiation report start\n")
        self.controller.datafile.write("display resolution: %sx%s\n" % (self.dispsize[0],self.dispsize[1]))
        self.controller.datafile.write("display size in cm: %sx%s\n" % (self.screensize[0],self.screensize[1]))
        self.controller.datafile.write("fixation threshold: %s degrees\n" % self.fixtresh)
        self.controller.datafile.write("speed threshold: %s degrees/second\n" % self.spdtresh)
        self.controller.datafile.write("acceleration threshold: %s degrees/second**2\n" % self.accthresh)
        self.controller.datafile.write("pygaze initiation report end\n")


    def calibrate(self, calibrate=True, validate=True):

        """Calibrates the eye tracking system

        arguments
        None

        keyword arguments
        calibrate    --    Boolean indicating if calibration should be
                    performed (default = True)
        validate    --    Boolean indicating if validation should be performed
                    (default = True)

        returns
        success    --    returns True if calibration succeeded, or False if
                    not; in addition a calibration log is added to the
                    log file and some properties are updated (i.e. the
                    thresholds for detection algorithms)
        """

        # calibration and validation points
        lb = int(0.1 * self.dispsize[0]) # left bound
        xc = int(0.5 * self.dispsize[0]) # horizontal center
        rb = int(0.9 * self.dispsize[0]) # right bound
        ub = int(0.1 * self.dispsize[1]) # upper bound
        yc = int(0.5 * self.dispsize[1]) # vertical center
        bb = int(0.9 * self.dispsize[1]) # bottom bound

        calpos = [(lb,ub), (rb,ub) , (xc,yc), (lb,bb), (rb,bb)]


        # # # # #
        # calibration

        while True:

            ret = self.controller.doCalibration(calpos)

            if ret == 'accept':
                calibrated =  True
                break

            elif ret == 'abort':
                calibrated =  False
                return calibrated


        # # # # #
        # validation

        # show menu
        self.screen.draw_text(text="Press space to start validation")
        self.disp.fill(self.screen)
        self.disp.show()
        self.screen.clear()

        # wait for spacepress
        self.kb.get_key(keylist=['space'],timeout=None)

        # start recording
        self.start_recording()

        # arrays for data storage
        lxacc = numpy.zeros(len(calpos))
        lyacc = numpy.zeros(len(calpos))
        rxacc = numpy.zeros(len(calpos))
        ryacc = numpy.zeros(len(calpos))

        # loop through all calibration positions
        for pos in calpos:
            # show validation point
            self.screen.draw_fixation(fixtype='dot', pos=pos)
            self.disp.fill(self.screen)
            self.disp.show()
            self.screen.clear()
            # allow user some time to gaze at dot
            clock.pause(1000)
            # collect samples
            lxsamples = numpy.zeros(self.nvalsamples)
            lysamples = numpy.zeros(self.nvalsamples)
            rxsamples = numpy.zeros(self.nvalsamples)
            rysamples = numpy.zeros(self.nvalsamples)
            prevsample = ()
            for i in range(0,self.nvalsamples):
                # add new sample to list
                newsample = self.controller.getCurrentGazePosition()
                if newsample != prevsample and newsample != (None,None,None,None):
                    lx, ly, rx, ry = newsample
                    lxsamples[i] = lx
                    lysamples[i] = ly
                    rxsamples[i] = rx
                    rysamples[i] = ry
                    prevsample = newsample[:]
            # calculate mean deviation
            lxdev = numpy.mean(abs(lxsamples-pos[0]))
            lydev = numpy.mean(abs(lysamples-pos[1]))
            rxdev = numpy.mean(abs(rxsamples-pos[0]))
            rydev = numpy.mean(abs(rysamples-pos[1]))
            i = calpos.index(pos)
            lxacc[i] = lxdev
            lyacc[i] = lydev
            rxacc[i] = rxdev
            ryacc[i] = rydev
            # wait for a bit to slow down validation process a bit
            clock.pause(2000)
            
        # AFTERMATH
        	# store some variables
        #pixpercm = (self.dispsize[0]/float(self.screensize[0]) + self.dispsize[1]/float(self.screensize[1])) / 2
        #screendist = settings.SCREENDIST
        # calculate mean accuracy
        self.pxaccuracy = [(numpy.mean(lxacc), numpy.mean(lyacc)), (numpy.mean(rxacc), numpy.mean(ryacc))]
        #self.pxaccuracy = ((deg2pix(screendist, self.accuracy[0][0], pixpercm),deg2pix(screendist, self.accuracy[0][1], pixpercm)), (deg2pix(screendist, self.accuracy[1][0], pixpercm),deg2pix(screendist, self.accuracy[1][1], pixpercm)))

        # # # # #
        # RMS noise

        # present instructions
        self.screen.draw_text(text="Noise calibration: please look at the dot\n\n(press space to start)", pos=(self.dispsize[0]/2, int(self.dispsize[1]*0.2)), center=True)
        self.screen.draw_fixation(fixtype='dot')
        self.disp.fill(self.screen)
        self.disp.show()
        self.screen.clear()

        # wait for spacepress
        self.kb.get_key(keylist=['space'], timeout=None)

        # show fixation
        self.screen.draw_fixation(fixtype='dot')
        self.disp.fill(self.screen)
        self.disp.show()
        self.screen.clear()

        # wait for a bit, to allow participant to fixate
        clock.pause(500)

        # get samples
        sl = [self.sample()] # samplelist, prefilled with 1 sample to prevent sl[-1] from producing an error; first sample will be ignored for RMS calculation
        t0 = clock.get_time() # starting time
        while clock.get_time() - t0 < 1000:
            s = self.sample() # sample
            if s != sl[-1] and s != (-1,-1) and s != (0,0):
                sl.append(s)

        # calculate RMS noise
        Xvar = []
        Yvar = []
        for i in range(2,len(sl)):
            Xvar.append((sl[i][0]-sl[i-1][0])**2)
            Yvar.append((sl[i][1]-sl[i-1][1])**2)
        XRMS = (sum(Xvar) / len(Xvar))**0.5
        YRMS = (sum(Yvar) / len(Yvar))**0.5
        self.pxdsttresh = (XRMS, YRMS)


        # # # # #
        # sample rate

        # calculate intersample times
        t0 = self.controller.gazeData[0]['system_time_stamp']
        ist = numpy.zeros(len(self.controller.gazeData)-1)
        for i in range(0,len(self.controller.gazeData)-1):
            ist[i] = (self.controller.gazeData[i+1]['system_time_stamp'] - self.controller.gazeData[i]['system_time_stamp']) / 1000.0

        # mean intersample time
        self.sampletime = numpy.mean(ist)
        self.samplerate = int(1000.0 / self.sampletime)

        # stop recording WITHOUT saving gaze data
        #self.controller.stopTracking()
        #self.controller.eyetracker.events.EYETRACKER_GAZE_DATA  -= self.controller.on_gazedata
        #self.eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA, self.on_gazedata)
        #self.controller.eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA,self.controller.on_gazedata_all)
        #self.controller.eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA,self.controller.on_gazedata)
        self.controller.gazeData = []
        self.controller.eventData = []
        self.recording = False


        # # # # #
        # calibration report

        # recalculate thresholds (degrees to pixels)
        self.pxfixtresh = deg2pix(self.screendist, self.fixtresh, self.pixpercm)
        self.pxspdtresh = deg2pix(self.screendist, self.spdtresh/1000.0, self.pixpercm) # in pixels per millisecons
        self.pxacctresh = deg2pix(self.screendist, self.accthresh/1000.0, self.pixpercm) # in pixels per millisecond**2

        # write report to log
        self.controller.datafile.write("pygaze calibration report start\n")
        self.controller.datafile.write("samplerate: %s Hz\n" % self.samplerate)
        self.controller.datafile.write("sampletime: %s ms\n" % self.sampletime)
        self.controller.datafile.write("accuracy (in pixels): LX=%s, LY=%s, RX=%s, RY=%s\n" % (self.pxaccuracy[0][0],self.pxaccuracy[0][1],self.pxaccuracy[1][0],self.pxaccuracy[1][1]))
        #self.controller.datafile.write("accuracy (in pixels): LX=%s, LY=%s, RX=%s, RY=%s\n" % ((deg2pix(screendist, self.pxaccuracy[0][0], pixpercm),deg2pix(screendist, self.accuracy[0][1], pixpercm)), (deg2pix(screendist, self.accuracy[1][0], pixpercm),deg2pix(screendist, self.accuracy[1][1], pixpercm)))
        self.controller.datafile.write("precision (RMS noise in pixels): X=%s, Y=%s\n" % (self.pxdsttresh[0],self.pxdsttresh[1]))
        self.controller.datafile.write("distance between participant and display: %s cm\n" % self.screendist)
        self.controller.datafile.write("fixation threshold: %s pixels\n" % self.pxfixtresh)
        self.controller.datafile.write("speed threshold: %s pixels/ms\n" % self.pxspdtresh)
        self.controller.datafile.write("acceleration threshold: %s pixels/ms**2\n" % self.pxacctresh)
        self.controller.datafile.write("pygaze calibration report end\n")

        return True


    def close(self):

        """Neatly close connection to tracker

        arguments
        None

        returns
        None        --    saves data and sets self.connected to False
        """

        # stop tracking
        if self.recording:
            self.stop_recording()

        # save data
        self.controller.closeDataFile()

        # close connection
        self.controller.destroy()
        self.connected = False


    def connected(self):

        """Checks if the tracker is connected

        arguments
        None

        returns
        connected    --    True if connection is established, False if not
        """

        return self.connected


    def drift_correction(self, pos=None, fix_triggered=False):

        """Performs a drift check

        arguments
        None

        keyword arguments
        pos            -- (x, y) position of the fixation dot or None for
                       a central fixation (default = None)
        fix_triggered    -- Boolean indicating if drift check should be
                       performed based on gaze position (fix_triggered
                       = True) or on spacepress (fix_triggered =
                       False) (default = False)

        returns
        checked        -- Boolaan indicating if drift check is ok (True)
                       or not (False); or calls self.calibrate if 'q'
                       or 'escape' is pressed
        """

        if fix_triggered:
            return self.fix_triggered_drift_correction(pos)

        if pos == None:
            pos = self.dispsize[0] / 2, self.dispsize[1] / 2

        # start recording if recording has not yet started
        if not self.recording:
            self.start_recording()
            stoprec = True
        else:
            stoprec = False

        result = False
        pressed = False
        while not pressed:
            pressed, presstime = self.kb.get_key()
            if pressed:
                if pressed == 'escape' or pressed == 'q':
                    print("libtobii.TobiiTracker.drift_correction: 'q' or 'escape' pressed")
                    return self.calibrate(calibrate=True, validate=True)
                gazepos = self.sample()
                if ((gazepos[0]-pos[0])**2  + (gazepos[1]-pos[1])**2)**0.5 < self.pxerrdist:
                    result = True
                else:
                    self.errorbeep.play()

        # stop recording WITHOUT saving gaze data
        if stoprec:
            #self.controller.stopTracking()
            #self.controller.eyetracker.events.EYETRACKER_GAZE_DATA -= self.controller.on_gazedata
            #self.controller.eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA,self.controller.on_gazedata_all)
            #self.controller.eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA,self.controller.on_gazedata)
            self.controller.gazeData = []
            self.controller.eventData = []
            self.recording = False

        return result


    def fix_triggered_drift_correction(self, pos=None, min_samples=10, max_dev=60, reset_threshold=30):

        """Performs a fixation triggered drift correction by collecting
        a number of samples and calculating the average distance from the
        fixation position

        arguments
        None

        keyword arguments
        pos            -- (x, y) position of the fixation dot or None for
                       a central fixation (default = None)
        min_samples        -- minimal amount of samples after which an
                       average deviation is calculated (default = 10)
        max_dev        -- maximal deviation from fixation in pixels
                       (default = 60)
        reset_threshold    -- if the horizontal or vertical distance in
                       pixels between two consecutive samples is
                       larger than this threshold, the sample
                       collection is reset (default = 30)

        returns
        checked        -- Boolaan indicating if drift check is ok (True)
                       or not (False); or calls self.calibrate if 'q'
                       or 'escape' is pressed
        """

        if pos == None:
            pos = self.dispsize[0] / 2, self.dispsize[1] / 2

        # loop until we have sufficient samples
        lx = []
        ly = []
        while len(lx) < min_samples:

            # pressing escape enters the calibration screen
            if self.kb.get_key()[0] in ['escape','q']:
                print("libtobii.TobiiTracker.fix_triggered_drift_correction: 'q' or 'escape' pressed")
                return self.calibrate(calibrate=True, validate=True)

            # collect a sample
            x, y = self.sample()

            if len(lx) == 0 or x != lx[-1] or y != ly[-1]:

                # if present sample deviates too much from previous sample, reset counting
                if len(lx) > 0 and (abs(x - lx[-1]) > reset_threshold or abs(y - ly[-1]) > reset_threshold):
                    lx = []
                    ly = []

                # collect samples
                else:
                    lx.append(x)
                    ly.append(y)

            if len(lx) == min_samples:

                avg_x = sum(lx) / len(lx)
                avg_y = sum(ly) / len(ly)
                d = ((avg_x - pos[0]) ** 2 + (avg_y - pos[1]) ** 2)**0.5

                if d < max_dev:
                    return True
                else:
                    lx = []
                    ly = []


    def get_eyetracker_clock_async(self):

        """Retrieve difference between tracker time and experiment time

        arguments
        None

        keyword arguments
        None

        returns
        timediff    --    tracker time minus experiment time
        """

        return self.controller.get_device_time() - clock.get_time()#self.controller.syncmanager.convert_from_local_to_remote(self.controller.clock.get_time()) - clock.get_time()


    def log(self, msg):

        """Writes a message to the log file

        arguments
        ms        -- a string to include in the log file

        returns
        Nothing    -- uses native log function of iViewX to include a line
                   in the log file
        """

        self.controller.recordEvent(msg)


    def log_var(self, var, val):

        """Writes a variable to the log file

        arguments
        var        -- variable name
        val        -- variable value

        returns
        Nothing    -- uses native log function of iViewX to include a line
                   in the log file in a "var NAME VALUE" layout
        """

        msg = "var %s %s" % (var, val)
        self.log(msg)


    def prepare_backdrop(self):

        """Not supported for TobiiTracker (yet)"""

        print("function not supported yet")


    def prepare_drift_correction(self, pos):

        """Not supported for TobiiTracker (yet)"""

        print("function not supported yet")


    def pupil_size(self):

        """Returns newest available pupil size

        arguments
        None

        returns
        pupilsize    --    a float or -1 on an error
        """

        # get new pupil size
        pupilsize = self.controller.getCurrentPupilSize()

        # if pupil size is invalid, return missing value
        if None in pupilsize:
            return -1

        # return pupil size
        if self.eye_used == self.left_eye or self.eye_used == self.binocular:
            return pupilsize[0]
        else:
            return pupilsize[1]


    def sample(self):

        """Returns newest available gaze position

        arguments
        None

        returns
        sample    -- an (x,y) tuple or a (-1,-1) on an error
        """

        # get new sample
        gazepos = self.controller.getCurrentGazePosition()

        # if sample is invalid, return missing value
        if None in gazepos:
            return (-1,-1)

        # return gaze x,y
        if self.eye_used == self.left_eye or self.eye_used == self.binocular:
            return gazepos[0], gazepos[1]
        else:
            return gazepos[2], gazepos[3]


    def send_command(self, cmd):

        """Not supported for TobiiTracker (yet)"""

        print("function not supported yet")


    def set_backdrop(self):

        """Not supported for TobiiTracker (yet)"""

        print("function not supported yet")


    def set_eye_used(self):

        """Logs the eye_used variable, based on which eye was specified
        (if both eyes are being tracked, the left eye is used)

        arguments
        None

        returns
        Nothing    -- logs which eye is used by calling self.log_var, e.g.
                   self.log_var("eye_used", "right")
        """

        if self.eye_used == self.right_eye:
            self.log_var("eye_used", "right")

        else:
            self.log_var("eye_used", "left")


    def start_recording(self):

        """Starts recording eye position

        arguments
        None

        returns
        None        -- sets self.recording to True when recording is
                   successfully started
        """

        if not self.recording:
            try:
                self.controller.startTracking()
                self.recording = True
            except:
                self.recording = False
                raise Exception("Error in libtobii.TobiiTracker.start_recording: failed to start recording")

        else:
            print("WARNING! libtobii.TobiiTracker.start_recording: already recording!")


    def status_msg(self, msg):

        """Not supported for TobiiTracker (yet)"""

        print("function not supported yet")


    def stop_recording(self):

        """Stop recording eye position

        arguments
        None

        returns
        Nothing    -- sets self.recording to False when recording is
                   successfully started
        """

        if self.recording:
            try:
                self.controller.stopTracking()
                self.recording = False
            except:
                self.recording = True
                raise Exception("Error in libtobii.TobiiTracker.stop_recording: failed to stop recording")

        else:
            print("WARNING! libtobii.TobiiTracker.stop_recording: recording has not started yet!")


    def set_detection_type(self, eventdetection):

        """Set the event detection type to either PyGaze algorithms, or
        native algorithms as provided by the manufacturer (only if
        available: detection type will default to PyGaze if no native
        functions are available)

        arguments
        eventdetection    --    a string indicating which detection type
                        should be employed: either 'pygaze' for
                        PyGaze event detection algorithms or
                        'native' for manufacturers algorithms (only
                        if available; will default to 'pygaze' if no
                        native event detection is available)
        returns        --    detection type for saccades, fixations and
                        blinks in a tuple, e.g.
                        ('pygaze','native','native') when 'native'
                        was passed, but native detection was not
                        available for saccade detection
        """

        # warn if detection is set to native
        if eventdetection == 'native':
            print("WARNING! 'native' event detection has been selected, \
                but Tobii does not provide detection algorithms; PyGaze \
                algorithm will be used instead")

        # set event detection methods to PyGaze
        self.eventdetection = 'pygaze'

        return (self.eventdetection,self.eventdetection,self.eventdetection)


    def wait_for_event(self, event):

        """Waits for event

        arguments
        event        -- an integer event code, one of the following:
                    3 = STARTBLINK
                    4 = ENDBLINK
                    5 = STARTSACC
                    6 = ENDSACC
                    7 = STARTFIX
                    8 = ENDFIX

        returns
        outcome    -- a self.wait_for_* method is called, depending on the
                   specified event; the return values of corresponding
                   method are returned
        """

        if event == 5:
            outcome = self.wait_for_saccade_start()
        elif event == 6:
            outcome = self.wait_for_saccade_end()
        elif event == 7:
            outcome = self.wait_for_fixation_start()
        elif event == 8:
            outcome = self.wait_for_fixation_end()
        elif event == 3:
            outcome = self.wait_for_blink_start()
        elif event == 4:
            outcome = self.wait_for_blink_end()
        else:
            raise Exception("Error in libtobii.TobiiTracker.wait_for_event: eventcode %s is not supported" % event)

        return outcome


    def wait_for_blink_end(self):

        """Waits for a blink end and returns the blink ending time

        arguments
        None

        returns
        timestamp        --    blink ending time in milliseconds, as
                        measured from experiment begin time
        """


        # # # # #
        # Tobii method

        if self.eventdetection == 'native':

            # print warning, since Tobii does not have a blink detection
            # built into their API

            print("WARNING! 'native' event detection has been selected, \
                but Tobii does not offer blink detection; PyGaze algorithm \
                will be used")

        # # # # #
        # PyGaze method

        blinking = True

        # loop while there is a blink
        while blinking:
            # get newest sample
            gazepos = self.sample()
            # check if it's valid
            if self.is_valid_sample(gazepos):
                # if it is a valid sample, blinking has stopped
                blinking = False

        # return timestamp of blink end
        return clock.get_time()


    def wait_for_blink_start(self):

        """Waits for a blink start and returns the blink starting time

        arguments
        None

        returns
        timestamp        --    blink starting time in milliseconds, as
                        measured from experiment begin time
        """

        # # # # #
        # Tobii method

        if self.eventdetection == 'native':

            # print warning, since Tobii does not have a blink detection
            # built into their API

            print("WARNING! 'native' event detection has been selected, \
                but Tobii does not offer blink detection; PyGaze algorithm \
                will be used")

        # # # # #
        # PyGaze method

        blinking = False

        # loop until there is a blink
        while not blinking:
            # get newest sample
            gazepos = self.sample()
            # check if it's a valid sample
            if not self.is_valid_sample(gazepos):
                # get timestamp for possible blink start
                t0 = clock.get_time()
                # loop until a blink is determined, or a valid sample occurs
                while not self.is_valid_sample(self.sample()):
                    # check if time has surpassed BLINKTHRESH
                    if clock.get_time()-t0 >= self.blinkthresh:
                        # return timestamp of blink start
                        return t0


    def wait_for_fixation_end(self):

        """Returns time and gaze position when a fixation has ended;
        function assumes that a 'fixation' has ended when a deviation of
        more than self.pxfixtresh from the initial fixation position has
        been detected (self.pxfixtresh is created in self.calibration,
        based on self.fixtresh, a property defined in self.__init__)

        arguments
        None

        returns
        time, gazepos    -- time is the starting time in milliseconds (from
                       expstart), gazepos is a (x,y) gaze position
                       tuple of the position from which the fixation
                       was initiated
        """

        # # # # #
        # Tobii method

        if self.eventdetection == 'native':

            # print warning, since Tobii does not have a fixation detection
            # built into their API

            print("WARNING! 'native' event detection has been selected, \
                but Tobii does not offer fixation detection; PyGaze algorithm \
                will be used")

        # # # # #
        # PyGaze method

        # function assumes that a 'fixation' has ended when a deviation of more than fixtresh
        # from the initial 'fixation' position has been detected

        # get starting time and position
        stime, spos = self.wait_for_fixation_start()

        # loop until fixation has ended
        while True:
            # get new sample
            npos = self.sample() # get newest sample
            # check if sample is valid
            if self.is_valid_sample(npos):
                # check if sample deviates to much from starting position
                if (npos[0]-spos[0])**2 + (npos[1]-spos[1])**2 > self.pxfixtresh**2: # Pythagoras
                    # break loop if deviation is too high
                    break

        return clock.get_time(), spos


    def wait_for_fixation_start(self):

        """Returns starting time and position when a fixation is started;
        function assumes a 'fixation' has started when gaze position
        remains reasonably stable (i.e. when most deviant samples are
        within self.pxfixtresh) for five samples in a row (self.pxfixtresh
        is created in self.calibration, based on self.fixtresh, a property
        defined in self.__init__)

        arguments
        None

        returns
        time, gazepos    -- time is the starting time in milliseconds (from
                       expstart), gazepos is a (x,y) gaze position
                       tuple of the position from which the fixation
                       was initiated
        """

        # # # # #
        # Tobii method

        if self.eventdetection == 'native':

            # print warning, since Tobii does not have a fixation start
            # detection built into their API (only ending)

            print("WARNING! 'native' event detection has been selected, \
                but Tobii does not offer fixation detection; PyGaze \
                algorithm will be used")


        # # # # #
        # PyGaze method

        # function assumes a 'fixation' has started when gaze position
        # remains reasonably stable for self.fixtimetresh

        # get starting position
        spos = self.sample()
        while not self.is_valid_sample(spos):
            spos = self.sample()

        # get starting time
        t0 = clock.get_time()

        # wait for reasonably stable position
        moving = True
        while moving:
            # get new sample
            npos = self.sample()
            # check if sample is valid
            if self.is_valid_sample(npos):
                # check if new sample is too far from starting position
                if (npos[0]-spos[0])**2 + (npos[1]-spos[1])**2 > self.pxfixtresh**2: # Pythagoras
                    # if not, reset starting position and time
                    spos = copy.copy(npos)
                    t0 = clock.get_time()
                # if new sample is close to starting sample
                else:
                    # get timestamp
                    t1 = clock.get_time()
                    # check if fixation time threshold has been surpassed
                    if t1 - t0 >= self.fixtimetresh:
                        # return time and starting position
                        return t0, spos


    def wait_for_saccade_end(self):

        """Returns ending time, starting and end position when a saccade is
        ended; based on Dalmaijer et al. (2013) online saccade detection
        algorithm

        arguments
        None

        returns
        endtime, startpos, endpos    -- endtime in milliseconds (from
                               expbegintime); startpos and endpos
                               are (x,y) gaze position tuples
        """

        # # # # #
        # Tobii method

        if self.eventdetection == 'native':

            # print warning, since Tobii does not have a blink detection
            # built into their API

            print("WARNING! 'native' event detection has been selected, \
                but Tobii does not offer saccade detection; PyGaze \
                algorithm will be used")

        # # # # #
        # PyGaze method

        # get starting position (no blinks)
        t0, spos = self.wait_for_saccade_start()
        # get valid sample
        prevpos = self.sample()
        while not self.is_valid_sample(prevpos):
            prevpos = self.sample()
        # get starting time, intersample distance, and velocity
        t1 = clock.get_time()
        s = ((prevpos[0]-spos[0])**2 + (prevpos[1]-spos[1])**2)**0.5 # = intersample distance = speed in px/sample
        v0 = s / (t1-t0)

        # run until velocity and acceleration go below threshold
        saccadic = True
        while saccadic:
            # get new sample
            newpos = self.sample()
            t1 = clock.get_time()
            if self.is_valid_sample(newpos) and newpos != prevpos:
                # calculate distance
                s = ((newpos[0]-prevpos[0])**2 + (newpos[1]-prevpos[1])**2)**0.5 # = speed in pixels/sample
                # calculate velocity
                v1 = s / (t1-t0)
                # calculate acceleration
                a = (v1-v0) / (t1-t0) # acceleration in pixels/sample**2 (actually is v1-v0 / t1-t0; but t1-t0 = 1 sample)
                # check if velocity and acceleration are below threshold
                if v1 < self.pxspdtresh and (a > -1*self.pxacctresh and a < 0):
                    saccadic = False
                    epos = newpos[:]
                    etime = clock.get_time()
                # update previous values
                t0 = copy.copy(t1)
                v0 = copy.copy(v1)
            # udate previous sample
            prevpos = newpos[:]

        return etime, spos, epos


    def wait_for_saccade_start(self):

        """Returns starting time and starting position when a saccade is
        started; based on Dalmaijer et al. (2013) online saccade detection
        algorithm

        arguments
        None

        returns
        endtime, startpos    -- endtime in milliseconds (from expbegintime);
                       startpos is an (x,y) gaze position tuple
        """

        # # # # #
        # Tobii method

        if self.eventdetection == 'native':

            # print warning, since Tobii does not have a blink detection
            # built into their API

            print("WARNING! 'native' event detection has been selected, \
                but Tobii does not offer saccade detection; PyGaze \
                algorithm will be used")

        # # # # #
        # PyGaze method

        # get starting position (no blinks)
        newpos = self.sample()
        while not self.is_valid_sample(newpos):
            newpos = self.sample()
        # get starting time, position, intersampledistance, and velocity
        t0 = clock.get_time()
        prevpos = newpos[:]
        s = 0
        v0 = 0

        # get samples
        saccadic = False
        while not saccadic:
            # get new sample
            newpos = self.sample()
            t1 = clock.get_time()
            if self.is_valid_sample(newpos) and newpos != prevpos:
                # check if distance is larger than precision error
                sx = newpos[0]-prevpos[0]; sy = newpos[1]-prevpos[1]
                if (sx/self.pxdsttresh[0])**2 + (sy/self.pxdsttresh[1])**2 > self.weightdist: # weigthed distance: (sx/tx)**2 + (sy/ty)**2 > 1 means movement larger than RMS noise
                    # calculate distance
                    s = ((sx)**2 + (sy)**2)**0.5 # intersampledistance = speed in pixels/ms
                    # calculate velocity
                    v1 = s / (t1-t0)
                    # calculate acceleration
                    a = (v1-v0) / (t1-t0) # acceleration in pixels/ms**2
                    # check if either velocity or acceleration are above threshold values
                    if v1 > self.pxspdtresh or a > self.pxacctresh:
                        saccadic = True
                        spos = prevpos[:]
                        stime = clock.get_time()
                    # update previous values
                    t0 = copy.copy(t1)
                    v0 = copy.copy(v1)

                # udate previous sample
                prevpos = newpos[:]

        return stime, spos


    def is_valid_sample(self, gazepos):

        """Checks if the sample provided is valid, based on Tobii specific
        criteria (for internal use)

        arguments
        gazepos        --    a (x,y) gaze position tuple, as returned by
                        self.sample()

        returns
        valid            --    a Boolean: True on a valid sample, False on
                        an invalid sample
        """

        # return False if a sample is invalid
        if gazepos == (-1,-1):
            return False
        # sometimes, on Tobii devices, invalid samples are the negative
        # display resolution value
        elif gazepos[0] == -1*self.dispsize[0] or gazepos[1] == -1*self.dispsize[1]:
            return False

        # in any other case, the sample is valid
        return True



# The code below is based on Hiroyuki Sogo's (Ehime University)
# TobiiController for PsychoPy. The code has been modified by Edwin Dalmaijer
# to work with PyGaze display routines; he added documentation as well.
# Furthermore the code has been modified By Alex Petrisor Giurca (Denmark Technical University) to interact with a Tobii4C and the latest tobii_pro_sdk
# more information: http://www.s12600.net/psy/etc/python.html#TobiiController
# original code: http://www.s12600.net/psy/etc/TobiiController/TobiiControllerP.py


class TobiiController:

    """Class to handle communication to Tobii eye trackers, as well as some
    display operations"""

    def __init__(self, disp):

        """Initializes TobiiController instance

        arguments
        disp        --    a pygaze.display.Display instance

        keyword arguments
        None
        """

        # visuals and interaction
        self.disp = disp
        self.screen = Screen()
        self.kb = Keyboard(keylist=None, timeout=None)

        # eye tracking
        self.address="tobii-ttp://is404-100106241134"
        self.eyetracker = tr.EyeTracker(self.address)
        #self.eyetrackers = tr.find_all_eyetrackers()
        self.gazeData = []
        self.eventData = []
        self.gazeDataAll = []
        #self.queue = Queue.Queue()
        #global activeTracking 
        self.stopRecording = False
        self.recordingState = False
        
        
        
        #callibration
        self.calibration = tr.ScreenBasedCalibration(self.eyetracker)
        # initialize communications
        
       
        
        with open("Here it should be the licence file name", "rb") as f:
            license = f.read()
        failed_licenses_applied_as_list_of_keys= self.eyetracker.apply_licenses([tr.LicenseKey(license)])
        if len(failed_licenses_applied_as_list_of_keys) == 0:
            print "Successfully applied license from list of keys."
        else:
            print "something went wrong"
        #self.eyetracker.activate(self)
            

    def doCalibration(self,calibrationPoints):
        
        """Performs a calibration; displaying points and the calibration
        menu and keyboard input are handled by PyGaze routines, calibration
        is handled by Tobii routines
        
        arguments
        calibrationPoints    --    a list of (x,y) typles, specifying the
                        coordinates for the calibration points
                        (coordinates should be in PyGaze notation,
                        where (0,0) is the topleft and coordinates
                        are specified in pixels, e.g. (1024,768))
        
        keyword arguments
        None
        
        returns
        None or retval    --    returns None if no tracker is connected
                        returns retval when a tracker is connected;
                        retval can be one of three string values:
                            'accept'
                            'retry'
                            'abort'
        """
        
        # immediately return when no eyetracker is connected
        if self.eyetracker is None:
            return
        
        # set some properties
        self.points = calibrationPoints
        self.point_index = -1
        
        # visuals
        img = Image.new('RGB',self.disp.dispsize)
        draw = ImageDraw.Draw(img)
        
        self.calin = {'colour':(0,0,0), 'pos':(int(self.disp.dispsize[0]/2),int(self.disp.dispsize[1]/2)), 'r':2}
        self.calout = {'colour':(128,255,128), 'pos':(int(self.disp.dispsize[0]/2),int(self.disp.dispsize[1]/2)), 'r':64}
        self.calresult = {'img':img}
        self.calresultmsg = {'text':"",'pos':(int(self.disp.dispsize[0]/2),int(self.disp.dispsize[1]/4))}
        
        # start calibration
        self.initcalibration_completed = False
        print "StartCalibration"
        self.calibration.enter_calibration_mode()
        self.initcalibration_completed = True
        while not self.initcalibration_completed:
            pass
        
        # draw central target
        self.screen.clear()
        self.screen.draw_circle(colour=self.calout['colour'], pos=self.calout['pos'], r=self.calout['r'], fill=False)
        self.screen.draw_circle(colour=self.calin['colour'], pos=self.calin['pos'], r=self.calin['r'], fill=True)
        self.disp.fill(self.screen)
        self.disp.show()
        # wait for start command
        self.kb.get_key(keylist=['space'],timeout=None)
        
        # run through all points
        for self.point_index in range(len(self.points)):
            # create tobii.eye_tracking_io.types 2D point
            px, py = self.points[self.point_index]
            p = Point2D()
            p.x, p.y = float(px)/self.disp.dispsize[0], float(py)/self.disp.dispsize[1]
            # recalculate to psycho coordinates
            self.calin['pos'] = (int(px),int(py))
            self.calout['pos'] = (int(px),int(py))

            # show target while decreasing its size for 1.5 seconds
            t0 = clock.get_time()
            currentTime = (clock.get_time() - t0) / 1000.0
            while currentTime < 1.5:
                # reduce size of the outer ring, as time passes
                self.calout['r'] = int(40*(1.5-(currentTime))+4)
                # check for input (should this even be here?)
                self.kb.get_key(keylist=None, timeout=1)
                # draw calibration point
                self.screen.clear()
                self.screen.draw_circle(colour=self.calout['colour'], pos=self.calout['pos'], r=self.calout['r'], fill=False)
                self.screen.draw_circle(colour=self.calin['colour'], pos=self.calin['pos'], r=self.calin['r'], fill=True)
                self.disp.fill(self.screen)
                self.disp.show()
                # get time
                currentTime = (clock.get_time() - t0) / 1000.0
            
            # wait for point calibration to succeed
            self.add_point_completed = False
            print p.x
            self.calibration.collect_data(p.x,p.y)
            print "Collecting data at {0}.".format(self.point_index)
            self.add_point_completed = True
            while not self.add_point_completed:
                # TODO: why would you continuously show the same stuff and poll the keyboard without using the input?
#                psychopy.event.getKeys()
#                self.calout.draw()
#                self.calin.draw()
#                win.flip()
                pass
         
        # wait for calibration to be complete
        #self.computeCalibration_completed = False
        #self.computeCalibration_succeeded = False
        self.calibration_result = self.calibration.compute_and_apply()
        print "Compute and apply returned {0} and collected at {1} points.".\
                                          format(self.calibration_result.status, len(self.calibration_result.calibration_points))
        self.computeCalibration_completed = True
        while not self.computeCalibration_completed:
            pass
        self.calibration.leave_calibration_mode()

        # reset display (same seems to be done below: what's the use?)
        self.disp.show()
        
        # get calibration info
        #self.getcalibration_completed = False
        #self.getcalibration_completed = True
        self.calib = tr.CalibrationPoint
        #while not self.getcalibration_completed:
            #pass
        
        # fill screen with half-gray
        self.screen.clear(colour=(128,128,128))
        
        # show calibration info
        if not tr.CALIBRATION_STATUS_SUCCESS:
            # computeCalibration failed.
            self.calresultmsg['text'] = 'Not enough data was collected (Retry:r/Abort:ESC)'
            
        elif self.calib == None:
            # no calibration data
            self.calresultmsg['text'] = 'No calibration data (Retry:r/Abort:ESC)'
            
        else:
            # show the calibration accuracy
            print "here we should have the plot"
            # alternative approach (Dalmaijer): use PyGaze drawing operations on self.screen, then present self.screen
            self.screen.draw_text(text=self.calresultmsg['text'],pos=self.calresultmsg['pos'])
            self.disp.fill(self.screen)
            self.disp.show()
        # original approach (Sogo): draw an image, then show that image via PscyhoPy
        self.calresult['img'] = img
        
        # show menu
        self.screen.draw_text(text="Press A to accept the calibration")
        self.disp.fill(self.screen)
        self.disp.show()
        
        
        
        # wait for keyboard input
        key, presstime = self.kb.get_key(keylist=['a','r','escape'], timeout=None)
        if key == 'a':
            retval = 'accept'
        elif key == 'r':
            retval = 'retry'
        elif key == 'escape':
            retval = 'abort'
        
        return retval
    

    def destroy(self):

        """Removes eye tracker and stops all operations

        arguments
        None

        keyword arguments
        None

        returns
        None        --    sets self.eyetracker and self.browser to None;
                    stops browser and the
                    tobii.eye_tracking_io.mainloop.MainloopThread
        """
        
        
        self.eyetracker = None
    


    ############################################################################
    # activation methods
    ############################################################################

    def activate(self,eyetracker):

        """Connects to specified eye tracker

        arguments
        eyetracker    --    key for the self.eyetracker dict under which the
                    eye tracker to which you want to connect is found

        keyword arguments
        None

        returns
        None        --    calls TobiiController.on_eyetracker_created, then
                    sets self.syncmanager
        """

        #eyetracker_info = self.eyetracker
        
        print "Connecting to:", self.eyetracker
        #self.eyetracker.subscribe_to(tr.EYETRACKER_GAZE_DATA, self.on_gazedata ,as_dictionary=True)
        self.recordingState = True
        #for having all the data in the right order, events and coordinates
        self.eyetracker.subscribe_to(tr.EYETRACKER_GAZE_DATA,self.on_gazedata_all,as_dictionary=True)
        #for calibration
        self.eyetracker.subscribe_to(tr.EYETRACKER_GAZE_DATA,self.on_gazedata,as_dictionary=True)


        while self.eyetracker==None:
            print "something went wrong activate"
            pass   
    

    ############################################################################
    # calibration methods
    ############################################################################
    
    
    def startTracking(self):

        """Starts the collection of gaze data

        arguments
        None

        keyword arguments
        None

        returns
        None        --    resets both self.gazeData and self.eventData, then
                    sets TobiiTracker.on_gazedata as an event callback
                    for self.eyetracker.events.OnGazeDataReceived and
                    calls self.eyetracker.StartTracking()
        """
        #global activeTracking
        self.gazeDataAll = []
        self.gazeData = []
        self.eventData = []
        
        #self.eyetracker.subscribe_to(tr.EYETRACKER_GAZE_DATA,self.on_gazedata,as_dictionary=True)
        if self.recordingState:
            self.stopRecording = True
            #self.recordingState = True
        else:
            self.stopRecording = False
        
    

    def stopTracking(self):

        """Starts the collection of gaze data

        arguments
        None

        keyword arguments
        None

        returns
        None        --    calls self.eyetracker.StopTracking(), then unsets
                    TobiiTracker.on_gazedata as an event callback for
                    self.eyetracker.events.OnGazeDataReceived, and
                    calls TobiiTracker.flushData before resetting both
                    self.gazeData and self.eventData
        """
        #global activeTracking
        #self.eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA, self.on_gazedata)
        
        #self.eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA, self.on_gazedata_all)
        #self.eyetracker.unsubscribe_from(tr.EYETRACKER_TIME_SYNCHRONIZATION_DATA,self.on_timedata)
        #self.eyetracker.unsubscribe_to(tr.EYETRACKER_GAZE_DATA, self.on_gazedata)
        #self.eyetracker.EYETRACKER_GAZE_DATA -= self.on_gazedata
        #print self.eyetracker['left_eye']
        #self.activeRecording = False
        
        self.sortData()
        self.gazeData = []
        self.eventData = []
        self.gazeDataAll = []
        if self.stopRecording:
            self.recordingState =False

    def on_gazedata(self,gaze):

        """Adds new data point to the data collection (self.gazeData)

        arguments
        error        --    some Tobii error message, isn't used in function
        gaze        --    Tobii gaze data struct

        keyword arguments
        None

        returns
        None        --    appends gaze to self.gazeData list
        """

        self.gazeData.append(gaze)
        
    def on_gazedata_all(self,gaze):

        """Adds new data point to the data collection (self.gazeData)

        arguments
        error        --    some Tobii error message, isn't used in function
        gaze        --    Tobii gaze data struct

        keyword arguments
        None

        returns
        None        --    appends gaze to self.gazeData list
        """
        #global activeTracking
        #print self.activeRecording
        
        self.gazeDataAll.append(gaze)
        
        
  

    def getPupilSize(self,gaze):

        """Extracts the pupil size of both eyes from the Tobii gaze struct.

        arguments
        gaze        --    Tobii gaze struct

        keyword aguments
        None

        returns
        pupilsiez    --    a (L,R) tuple for the pupil size in mm of both eyes
        """

        return (gaze['left_pupil_diameter'],
                gaze['right_pupil_diameter'])

    def getCurrentPupilSize(self):

        """Provides the newest pupil size

        arguments
        None

        keyword argument
        None

        returns
        pupilsize    --    a (L,R) tuple for the pupil size in mm of both eyes
                    or (None,None) if no new size is available
        """

        if len(self.gazeData) == 0:
            return (None,None)
        else:
            return self.getPupilSize(self.gazeData[-1])

    def getGazePosition(self,gaze):

        """Extracts the gaze positions of both eyes from the Tobii gaze
        struct and recalculates them to PyGaze coordinates

        arguments
        gaze        --    Tobii gaze struct

        keyword arguments
        None

        returns
        gazepos    --    a (Lx,Ly,Rx,Ry) tuple for the gaze positions
                    of both eyes
        """
        #print gaze['left_gaze_point_on_display_area']
        return (gaze['left_gaze_point_on_display_area'][0]*self.disp.dispsize[0], #if gaze['left_gaze_point_on_display_area'][0]!='nan' else -1,0,
                gaze['left_gaze_point_on_display_area'][1]*self.disp.dispsize[1], # if gaze['left_gaze_point_on_display_area'][1]!='nan'else -1,0,
                gaze['right_gaze_point_on_display_area'][0]*self.disp.dispsize[0], # if gaze['left_gaze_point_on_display_area'][0]!='nan' else -1,0,
                gaze['right_gaze_point_on_display_area'][1]*self.disp.dispsize[1])
            #(int(gaze['left_gaze_point_on_display_area']*self.disp.dispsize[0]))

    def getCurrentGazePosition(self):

        """Provides the newest gaze position sample

        arguments
        None

        keyword arguments
        None

        returns
        gazepos    --    a (Lx,Ly,Rx,Ry) tuple for the gaze positions
                    of both eyes or (None,None,None,None) if no new
                    sample is available
        """
        if len(self.gazeData)==0:
            return (None,None,None,None)
        else:
            return self.getGazePosition(self.gazeData[-1])

    def setDataFile(self,filename):

        """Opens a new textfile and writes date, time and screen resolution
        to it

        arguments
        filename    --    a string containing the filename, including an
                    extension

        keyword arguments
        None

        returns
        None        --    sets self.datafile to an open textfile
        """

        print 'set datafile ' + filename
        self.datafile = open(filename,'w')
        self.datafile.write('Recording date:\t'+datetime.datetime.now().strftime('%Y/%m/%d')+'\n')
        self.datafile.write('Recording time:\t'+datetime.datetime.now().strftime('%H:%M:%S')+'\n')
        self.datafile.write('Recording resolution\t%d x %d\n\n' % tuple(self.disp.dispsize))


    def closeDataFile(self):

        """Closes the datafile after writing the last data to it

        arguments
        None

        keyword arguments
        None

        returns
        None        --    calls TobiiController.flushData, then closes
                    self.datafile and sets self.datafile to None
        """

        print 'datafile closed'
        if self.datafile != None:
            self.sortData()
            self.eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA, self.on_gazedata_all)
            self.eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA, self.on_gazedata)

            self.datafile.close()

        self.datafile = None


    def recordEvent(self,event):

        """Adds an event to the event data

        arguments
        event        --    a string containing an event description

        keyword arguments
        None

        returns
        None        --    appends a (timestamp,event) tuple to
                    self.eventData
        """
        #t = self.gazeData[0]['device_time_stamp']
        
        t = tr.get_system_time_stamp()
        t = t/1000.0
        t= round(t,1)
       
        #t = self.syncmanager.convert_from_local_to_remote(self.clock.get_time())
        
        self.gazeDataAll.append((t,event))
        
       

    def flushData(self):
        if self.datafile == None:
            print 'data file is not set.'
            return
        
        if len(self.gazeData)==0:
            print("WARNING! libtobii.TobiiController.flushData: no data to write to file.")
            return
        
        self.datafile.write('\t'.join(['TimeStamp',
                                       'GazePointXLeft',
                                       'GazePointYLeft',
                                       'ValidityLeft',
                                       'GazePointXRight',
                                       'GazePointYRight',
                                       'ValidityRight',
                                       'GazePointX',
                                       'GazePointY',
                                       'Event'])+'\n')
        timeStampStart = self.gazeData[0]['system_time_stamp']
        for g in self.gazeData:
            
            self.datafile.write('%.1f\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t%d'%(
                                (g['system_time_stamp']-timeStampStart)/1000.0,
                                g['left_gaze_point_on_display_area'][0]*self.disp.dispsize[0] if g['left_pupil_validity']!=4 else -1.0,
                                 g['left_gaze_point_on_display_area'][1]*self.disp.dispsize[1] if g['left_pupil_validity']!=4 else -1.0,
                                  g['left_pupil_validity'],
                                g['right_gaze_point_on_display_area'][0]*self.disp.dispsize[0] if g['right_pupil_validity']!=4 else -1.0,
                                 g['right_gaze_point_on_display_area'][1]*self.disp.dispsize[1] if g['right_pupil_validity']!=4 else -1.0,
                                  g['right_pupil_validity']))
            
            # if no correct sample is available, data is missing
            if g['left_pupil_validity'] == 4 and g['right_pupil_validity'] == 4: #not detected
                ave = (-1.0,-1.0)
            # if the right sample is unavailable, use left sample
            elif g['left_pupil_validity'] == 4:
                ave = (g['right_gaze_point_on_display_area'][0],g['right_gaze_point_on_display_area'][1])
                # if the left sample is unavailable, use right sample
            elif g['right_pupil_validity'] == 4:
                ave = (g['left_gaze_point_on_display_area'][0],g['left_gaze_point_on_display_area'][1])
            # if we have both samples, use both samples
            else:
                ave = ((g['left_gaze_point_on_display_area'][0] + g['right_gaze_point_on_display_area'][0]) / 2.0,
					   (g['left_gaze_point_on_display_area'][1] + g['right_gaze_point_on_display_area'][1]) / 2.0)
            print "this is the system_time_stamp"
            print g['system_time_stamp']
            
            # write gaze position to the datafile, based on the selected sample(s)
            self.datafile.write('\t%.4f\t%.4f\t'%ave)
            self.datafile.write('\n')
            
            #ave = ((g['left_gaze_point_on_display_area'] + g['left_gaze_point_on_display_area']) / 2.0,
                   #g['right_point_on_display_area'] + g['right_gaze_point_on_display_area'])
                
            #self.datafile.write('\t%.4f' % ave)
            #      self.datafile.write('\n')
        
        formatstr = '%.1f'+'\t'*9+'%s\n'
        for e in self.eventData:
            self.datafile.write(formatstr % ((e[0]-timeStampStart)/1000.0,e[1]))
        
        

        self.datafile.flush() # internal buffer to RAM
        os.fsync(self.datafile.fileno()) # RAM file cache to disk
        
    def sortData(self):
        if self.datafile == None:
            print ('data file is not set.')
            return
        
        if len(self.gazeDataAll)==0:
            print("WARNING! libtobii.TobiiController.flushData: no data to write to file.")
            return
        
        self.datafile.write('\t'.join(['TimeStamp',
                                       'GazePointXLeft',
                                       'GazePointYLeft',
                                       'ValidityLeft',
                                       'GazePointXRight',
                                       'GazePointYRight',
                                       'ValidityRight',
                                       'GazePointX',
                                       'GazePointY',
                                       'AvgPupilDiamX',
                                       'AvgPupilDiamY',
                                       'Event'])+'\n')
        
        valid = ["validation","recording","image","imgname","trial","pupdata","baseline","fixation","target"]
        for g in self.gazeDataAll:
            #g = str(g).lower()
            if not any(x in str(g).lower() for x in valid):
                self.datafile.write('%.1f\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t%d'%(
                                (g['system_time_stamp'])/1000.0,
                                g['left_gaze_point_on_display_area'][0]*self.disp.dispsize[0] if g['left_pupil_validity']!=0 else -1.0,
                                 g['left_gaze_point_on_display_area'][1]*self.disp.dispsize[1] if g['left_pupil_validity']!=0 else -1.0,
                                  g['left_pupil_validity'],
                                g['right_gaze_point_on_display_area'][0]*self.disp.dispsize[0] if g['right_pupil_validity']!=0 else -1.0,
                                 g['right_gaze_point_on_display_area'][1]*self.disp.dispsize[1] if g['right_pupil_validity']!=0 else -1.0,
                                  g['right_pupil_validity']))
                    # if no correct sample is available, data is missing
                if g['left_pupil_validity'] == 0 and g['right_pupil_validity'] == 0: #not detected
                    ave = (-1.0,-1.0)
                # if the right sample is unavailable, use left sample
                elif g['left_pupil_validity'] == 0:
                    ave = (g['right_gaze_point_on_display_area'][0],g['right_gaze_point_on_display_area'][1])
                    # if the left sample is unavailable, use right sample
                elif g['right_pupil_validity'] == 0:
                    ave = (g['left_gaze_point_on_display_area'][0],g['left_gaze_point_on_display_area'][1])
                # if we have both samples, use both samples
                else:
                    ave = ((g['left_gaze_point_on_display_area'][0] + g['right_gaze_point_on_display_area'][0]) / 2.0,
    					   (g['left_gaze_point_on_display_area'][1] + g['right_gaze_point_on_display_area'][1]) / 2.0)
                
                
                # write gaze position to the datafile, based on the selected sample(s)
                self.datafile.write('\t%.4f\t%.4f\t'%ave)
                
                     # if no correct sample is available, data is missing
                if g['left_pupil_validity'] == 0 and g['right_pupil_validity'] == 0: #not detected
                    avePup = (-1.0,-1.0)
                # if the right sample is unavailable, use left sample
                elif g['left_pupil_validity'] == 0:
                    avePup = (g['right_pupil_diameter'],g['right_pupil_diameter'])
                    # if the left sample is unavailable, use right sample
                elif g['right_pupil_validity'] == 0:
                    avePup = (g['left_pupil_diameter'],g['left_pupil_diameter'])
                # if we have both samples, use both samples
                else:
                    avePup = ((g['left_pupil_diameter'] + g['right_pupil_diameter']) / 2.0,
    					   (g['left_pupil_diameter'] + g['right_pupil_diameter']) / 2.0)
                # write gaze position to the datafile, based on the selected sample(s)
                self.datafile.write('\t%.4f\t%.4f\t'%avePup)
                #print (str(g))
                
                
            else:
                e1 = str(g)
                
                e1 = re.sub('[^a-zA-Z0-9-_\=*.]', ' ', e1)
                e = e1.split()
                
                
                
                self.datafile.write("MSG")
                #print(str(e[0]))
                for each in e:
                    self.datafile.write("\t" + str(each))
            self.datafile.write("\n")
            #print(str(g))
            
        self.datafile.flush()
        os.fsync(self.datafile.fileno())

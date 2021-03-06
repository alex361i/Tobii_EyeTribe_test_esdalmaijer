import os

# main
DIR = os.path.abspath(os.path.dirname(__file__))
IMGDIR = os.path.join(DIR, 'imgs')
DATADIR = os.path.join(DIR, 'data')
DATADIRTX = os.path.join(DIR, 'data\experiment4ctx300')
DATADIRTXFINAL = os.path.join(DIR, 'data\experimentTobiiTx300Tobii4CEyeTribe')
#LOGFILENAME = raw_input("Participant: ")
#LOGFILE = os.path.join(DATADIR, LOGFILENAME)

# display
DISPTYPE = 'pygame'
#DISPSIZE = (1366,768)
DISPSIZE = (1920,1080)


# validation points
PXY = [0.1, 0.5, 0.9]
CALIBPOINTS = []
for x in PXY:
	for y in PXY:
		CALIBPOINTS.append((int(x*DISPSIZE[0]),(int(y*DISPSIZE[1]))))

# images
IMAGES = []
for imgname in os.listdir(IMGDIR):
	IMAGES.append(os.path.join(IMGDIR, imgname))

# light-dark
PUPTRIALS = 50
SACTRIALS = 50

# timing
ITI = 1000
POINTTIME = 2000
IMGTIME = 10000
BASELINETIME = 200
PUPTRIALTIME = 2500

# trackers
TRACKERTYPE = 'tobiitx'
DUMMYMODE = False

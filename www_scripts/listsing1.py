
import os

urlbase = 'https://s3.amazonaws.com/nist-ics-wireless/measurements/'

myDirPaths = [
	"active_meas",
	"cloud_meas",
	"oats",
	"passive_meas/NIST MachineShop", 
	"passive_meas/TransmissionAssembly/TX1 position M18", 
	"passive_meas/TransmissionAssembly/TX2 position P10_Q10"]
print(myDirPaths)

# This section lists groups into a single file
ndirs = len(myDirPaths)
myOutFile = "downloads.htm"
f = open(myOutFile, "w")
myGroupHeaders = [
	"Active Mobile", 
	"Active Stationary", 
	"Active OATS", 
	"Passive Machine Shop", 
	"Passive Automotive 1", 
	"Passive Automotive 2" ]
anchors = ["active","cloud","oats","passive_shop","passive_m18","passive_P10_Q10" ]
f.write('<!DOCTYPE htm PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">')
f.write('<htm>\n<head>\n')
  
f.write('</head>\n<meta http-equiv="Content-Type" content="text/htm; charset=utf-8">\n')
f.write('<body>\n<ul>\n')
f.write('<style>table, th, td { border: 1px solid black; border-collapse: collapse; padding: 8px;}</style>\n')
f.write('<h1>Industrial Wireless Downloads</h1>')
introtext = 'This page contains links to all of the files currently available for download.  \
 To download, click the link or right click and select Save-As.  Our downloadable \
 files are stored on an external server which may cause some download difficulties. \
 Utilities such as DownloadThemAll! for FireFox are available to make downloading entire sections easier. \
 <p>We have tried to make accessing the data as easy as possible given our resources. Please \
 report broken links, missing data, or other requests to our project leader \
 listed on our <a href="https://www.nist.gov/programs-projects/wireless-systems-industrial-environments">project home page</a>.</p>'
f.write(introtext)


f.write('<h3><u>Table of Contents</u></h3>\n')
for ii in range(len(anchors)):
	f.write('<li><a href=#%s>%s</a><br></li>\n' % (anchors[ii], myGroupHeaders[ii]))
f.write('<br><hr>\n')
'''
f.write('<a name="%s"></a>' % 'reduced_tap_delay_profiles')
f.write('<u><h2>Reduced Tap Channel Impulse Response Files</h2></u>')
f.write('<p> Tap reduction algorithm is described in: </p><blockquote>C. Mehlfuhrer and M. Rupp, "Approximation and resampling of tapped delay \
		 line channel models with guaranteed channel properties," in 2008 IEEE \
		 International Conference on Acoustics, Speech and Signal Processing, 2008, \
		 no. 2, pp. 2869-2872.</blockquote>')'''
f.write('<u><h2>Measurement Files</h2></u>')
f.write('<h3>META Data</h3>')
f.writelines(['<a href="%s%s">%s</a>...&nbsp' % (urlbase, "SmartManufacturing_MetaData.xlsx", "SmartManufacturing_MetaData.xlsx")])
f.write('This file provides all of the situational data about instrumentation settings and antenna position.\n')
for ii in range(ndirs):
	print(ii)
	myDirPath = myDirPaths[ii]
	print("Working on %s ==> %s" % (myDirPath, myOutFile))
	myFiles = os.listdir(myDirPath)
	#print(myFiles)
	f.write('<a name="%s"></a>' % anchors[ii])
	f.write('<h3> %s </h3>' % myGroupHeaders[ii])
	f.writelines(['<a href="%s%s/%s">%s</a>...&nbsp\n' % (urlbase, myDirPath, fx, fx) for fx in myFiles if fx.endswith('.mat') or fx.endswith('.csv') or fx.endswith('.txt')])

f.write('<p>End of Downloads</p>\n')
f.write('<p><br></p>\n')
f.write('<p><br></p>\n')
f.write('</ul>\n</body>\n</htm>\n')
f.close()
